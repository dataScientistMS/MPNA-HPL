#!/bin/bash

# Configuration du compte et du binaire à tester
ACCOUNT="r250142"
BINARY="./xhpl_custom"

# Valeurs par défaut (notre base de comparaison)
B_N=115000; B_NB=224; B_NODES=2

# Cette fonction fabrique un dossier et un script de soumission Slurm pour chaque test
generate_test() {
    local series=$1; local folder=$2; local N=$3; local NB=$4; local NODES=$5
    local CORES=$((NODES * 64)) # On utilise les 64 coeurs par noeud de Romeo
    
    mkdir -p "CUSTOM_$series/$folder"
    
    # Création du fichier .slurm à la volée
    cat > "CUSTOM_$series/$folder/run.slurm" << EOS
#!/bin/bash
#SBATCH --account=$ACCOUNT
#SBATCH --constraint=x64cpu
#SBATCH --partition=short
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$CORES
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --job-name=C_$folder
#SBATCH --output=job.out

# Chargement des bibliothèques nécessaires sur le nœud de calcul
romeo_load_x64cpu_env
module load aocl-4.2.0 openmpi/aocc-4.2.0/4.1.7

# On s'assure que le multithreading ne parasite pas nos tests MPI
export OMP_NUM_THREADS=1

# Exécution du binaire avec les paramètres N et NB choisis
srun ../../$BINARY $N $NB
EOS
}

echo "Préparation de l'étude comparative..."

# Série S1 : On regarde comment le code réagit quand la matrice grossit
for val in 60000 80000 100000 115000 130000; do
    generate_test "S1_N" "N_$val" $val $B_NB $B_NODES
done

# Série S2 : On cherche la taille de bloc optimale pour le cache
for val in 128 192 224 256 384; do
    generate_test "S2_NB" "NB_$val" $B_N $val $B_NODES
done

# Série S5 : Test de montée en charge (scalabilité)
# Note : N augmente avec le nombre de nœuds pour garder la charge par nœud équilibrée
generate_test "S5_Scale" "1Node"  80000  $B_NB 1
generate_test "S5_Scale" "2Nodes" 115000 $B_NB 2
generate_test "S5_Scale" "4Nodes" 160000 $B_NB 4
generate_test "S5_Scale" "6Nodes" 200000 $B_NB 6

echo "Tous les dossiers sont prêts."
echo "Astuce pour tout lancer : for d in CUSTOM_S*/*/; do cd \$d && sbatch run.slurm && cd ../..; done"