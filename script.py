import mdtraj as md
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.preprocessing import MinMaxScaler

def load_structure(file_path):
    return md.load(file_path)

def identify_cleavage_site(structure, site_sequence):
    protein_seq = structure.top.to_fasta()[0]
    return protein_seq.find(site_sequence)

def find_accessible_region(structure, cleavage_site, distance_threshold):
    atom_indices = [atom.index for atom in structure.top.atoms if atom.element.symbol != 'H']
    query_indices = np.array([cleavage_site])
    neighbors = md.compute_neighbors(structure, distance_threshold, query_indices, atom_indices)
    return np.unique(np.concatenate(neighbors))

def design_candidate_peptides(structure, accessible_region, peptide_length_range):
    protein_seq = structure.top.to_fasta()[0]
    candidate_peptides = []
    for start_idx in accessible_region:
        for peptide_length in peptide_length_range:
            end_idx = start_idx + peptide_length
            if end_idx < len(protein_seq):
                candidate_peptides.append(protein_seq[start_idx:end_idx])
    return candidate_peptides

def evaluate_immunogenicity(peptides):
    # Define the hydrophobicity scale
    hydrophobicity_scale = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
        'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
        'W': -0.9, 'Y': -1.3
    }
    
    scores = []
    for peptide in peptides:
        seq = Seq(peptide)
        pa = ProteinAnalysis(str(seq))
        
        # Calculate instability index
        instability_index = pa.instability_index()

        # Calculate hydrophobicity using the hydrophobicity scale
        hydrophobicity = np.mean([hydrophobicity_scale.get(aa, 0) for aa in peptide])
        
        # Calculate charge
        charge = pa.charge_at_pH(7.4)
        
        # Calculate propensity to form secondary structures
        helix, turn, sheet = pa.secondary_structure_fraction()
        secondary_structure_score = helix + sheet

        # Combine the scores using your preferred weighting
        combined_score = instability_index + hydrophobicity + charge + secondary_structure_score
        scores.append(combined_score)

    scaler = MinMaxScaler()
    return scaler.fit_transform(np.array(scores).reshape(-1, 1)).flatten()

def optimize_peptide_sequence(peptides, immunogenicity_scores, top_n=1):
    top_peptide_indices = np.argsort(immunogenicity_scores)[-top_n:]
    top_peptides = [peptides[i] for i in top_peptide_indices]
    return top_peptides


if __name__ == "__main__":
    file_path = "c5.pdb"
    site_sequence = "VNND"
    distance_threshold = 10
    peptide_length_range = range(10, 21)  # Generate peptide lengths from 10 to 20 inclusive
    
    structure = load_structure(file_path)
    protein_seq = structure.top.to_fasta()[0]

    print(f"Protein sequence in FASTA format:\n>{structure.topology.chain(0)}\n{protein_seq}")
    
    #cleavage_site = identify_cleavage_site(structure, site_sequence)

    cleavage_site = 73  # Directly set the cleavage site index

    # Set the number of top peptides you'd like to select
    top_n = 5

    if cleavage_site == -1:
        print("Cleavage site sequence not found in the protein. Please ensure the site_sequence variable is correct.")
    else:
        accessible_region = find_accessible_region(structure, cleavage_site, distance_threshold)
        candidate_peptides = design_candidate_peptides(structure, accessible_region, peptide_length_range)
        immunogenicity_scores = evaluate_immunogenicity(candidate_peptides)
        optimized_peptides = optimize_peptide_sequence(candidate_peptides, immunogenicity_scores, top_n)
        
        print(f"Top {top_n} optimized peptides:")
        for i, peptide in enumerate(optimized_peptides, start=1):
            print(f"{i}. {peptide}")