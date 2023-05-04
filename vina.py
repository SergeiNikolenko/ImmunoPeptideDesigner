import os
import glob
import subprocess
import shutil


def get_next_folder_name(base_folder, folder_prefix):
    folder_number = 1
    while os.path.exists(os.path.join(base_folder, f"{folder_prefix}{folder_number}")):
        folder_number += 1
    return f"{folder_prefix}{folder_number}"

def read_params_from_file(filename):
    params = {}
    with open(filename, 'r') as file:
        for line in file:
            if '=' in line:
                key, value = line.strip().split('=')
                params[key] = value
    return params

# Автоматически задайте имя папки с лигандами и рецепторами
ligand_folder = 'ligands'
receptor_folder = 'receptors'

# Проверьте наличие папки ligands
if not os.path.exists(ligand_folder):
    print(f"Directory '{ligand_folder}' not found. Please check the current working directory and folder name.")
    exit()


params_filename = 'params.txt'
params = read_params_from_file(params_filename)

# Добавьте определение переменных здесь
tries = int(params.get('tries', 10))
xtn = int(params.get('exhaustiveness', 8))
nt = int(params.get('cpu', 16))
num_mod = int(params.get('num_mod', 20))
energy_range = int(params.get('energy_range', 20))
os.chdir(ligand_folder)
vinastart = open('vinastart.sh', 'w')
vinastart.write('#!/bin/bash' + '\n')
mols = glob.glob('./*.pdbqt')
count = len(mols)
h = 1

# Загрузите имена всех файлов PDBQT рецепторов из папки рецепторов
os.chdir('../' + receptor_folder)
receptor_files = glob.glob('./*.pdbqt')
receptor_names = [os.path.splitext(os.path.basename(r))[0] for r in receptor_files]

# Задайте имя общего лог-файла
combined_log_file = 'combined_log.txt'

os.chdir('../' + ligand_folder)

for receptor in receptor_names:
    i = 0
    for i in range(0, count):
        t = 1
        while t <= tries:
            f = open(mols[i])
            NAME = str(mols[i])
            name = (NAME[2:-6])
            config = open(receptor + "_" + name + '_' + str(t) + '.cfg', 'w')
            config.write('receptor=../' + receptor_folder + '/' + receptor + '.pdbqt' + '\n' + 'ligand=' + name + '.pdbqt' + '\n' + ' ' + '\n' + 'center_x=' + params['cx'] + '\n' + 'center_y=' + params['cy'] + '\n' + 'center_z=' + params['cz'] + '\n' + ' ' + '\n' + 'size_x=' + params['sx'] + '\n' + 'size_y=' + params['sy'] + '\n' + 'size_z=' + params['sz'] + '\n' + 'exhaustiveness=' + str(xtn) + '\n' + 'num_modes=' + str(num_mod) + '\n' + 'energy_range=' + str(energy_range) + '\n' + 'cpu=' + str(nt) + '\n' + 'out=' + receptor + '_' + name + '_' + str(t) + '_out.pdbqt')
            #print('Config file for ' + mols[i] + ' ' + 'is written!')
            vinastart.write('vina' + ' ' + '--config' + ' ' + receptor + '_' + name + '_' + str(t) + '.cfg' + ' >> ' + combined_log_file + '\n')
            #print('Run command for ' + mols[i] + ' ' + 'is written to vinastart!')
            t += 1
vinastart.close()
print('Config files generation is complete!')
print('The vinastart is ready to run!')

# Добавление команды удаления файлов конфигурации в файл vinastart.sh
with open("vinastart.sh", "a") as vinastart:
    vinastart.write("echo 'Deleting configuration files...'\n")
    for receptor in receptor_names:
        i = 0
        for i in range(0, count):
            t = 1
            while t <= tries:
                NAME = str(mols[i])
                name = (NAME[2:-6])
                config_file = receptor + "_" + name + '_' + str(t) + '.cfg'
                vinastart.write(f"rm {config_file}\n")
                t += 1
    vinastart.write("echo 'Configuration files deleted.'\n")
# Создание папки resultN
results_folder = get_next_folder_name(os.path.dirname(os.path.abspath(__file__)), "result")
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), results_folder))

# Добавление перемещения файлов с результатами и общего лога в папку resultN после выполнения vinastart.sh
with open("vinastart.sh", "a") as vinastart:
    vinastart.write("echo 'Moving result files and combined log to results folder...'\n")
    for receptor in receptor_names:
        i = 0
        for i in range(0, count):
            t = 1
            while t <= tries:
                NAME = str(mols[i])
                name = (NAME[2:-6])
                output_file = receptor + "_" + name + '_' + str(t) + '_out.pdbqt'
                vinastart.write(f"mv {output_file} ../{results_folder}/\n")
                t += 1
    vinastart.write(f"mv {combined_log_file} ../{results_folder}/\n")
    vinastart.write("echo 'Result files and combined log moved to results folder.'\n")


os.chmod('vinastart.sh', 0o755)
subprocess.run('./vinastart.sh', shell=True, check=True)
