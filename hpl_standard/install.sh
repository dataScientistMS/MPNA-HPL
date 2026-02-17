#!/bin/bash

# Le script définit d'abord les variables de chemin et les versions de référence
ACCOUNT="r250142"
HPL_VERSION="2.3"
INSTALL_DIR="$HOME/hpl-${HPL_VERSION}"
CUSTOM_SOURCE="hpl_custom.c"
CUSTOM_EXE="xhpl_custom"

# Il charge l'environnement logiciel nécessaire à l'exploitation du cluster Romeo
echo "Chargement des modules système..."
romeo_load_x64cpu_env
module load aocl-4.2.0 openmpi/aocc-4.2.0/4.1.7

# Cette fonction gère l'installation de la version officielle de Netlib
install_standard_hpl() {
    cd $HOME
    if [ ! -d "$INSTALL_DIR" ]; then
        echo "Récupération du code source Netlib..."
        wget http://www.netlib.org/benchmark/hpl/hpl-${HPL_VERSION}.tar.gz
        tar xf hpl-${HPL_VERSION}.tar.gz
    fi
    
    cd $INSTALL_DIR
    # Le script génère un fichier Make.romeo adapté aux processeurs Zen 3 du cluster
    cat > Make.romeo << 'EOF_MAKE'
ARCH         = romeo
TOPdir       = $(HOME)/hpl-2.3
INCdir       = $(TOPdir)/include
BINdir       = $(TOPdir)/bin/$(ARCH)
LIBdir       = $(TOPdir)/lib/$(ARCH)
HPLlib       = $(LIBdir)/libhpl.a
MPdir        = $(MPI_ROOT)
MPinc        = -I$(MPdir)/include
MPlib        = -L$(MPdir)/lib -lmpi
LAdir        = $(AOCL_ROOT)
LAinc        = -I$(LAdir)/include
LAlib        = -L$(LAdir)/lib -lblis
CC           = mpicc
CCFLAGS      = $(HPL_DEFS) -O3 -march=znver3 -fomit-frame-pointer
LINKER       = $(CC)
LINKFLAGS    = $(CCFLAGS)
EOF_MAKE

    echo "Compilation du binaire standard..."
    make arch=romeo
}

# Il compile ici l'implémentation personnalisée en la liant à la bibliothèque BLIS
compile_custom_hpl() {
    if [ -f "$CUSTOM_SOURCE" ]; then
        echo "Compilation du binaire personnalisé (1D)..."
        mpicc -O3 -march=znver3 "$CUSTOM_SOURCE" -o "$CUSTOM_EXE" -lblis -lm
    else
        echo "Erreur : le fichier source est introuvable."
    fi
}

# Cette fonction automatise la création du plan d'expérience (25 scénarios)
generate_test_plan() {
    # Définition des paramètres de la baseline (2 nœuds de calcul)
    B_N=115000; B_NB=224; B_P=8; B_Q=16; B_NODES=2; B_PFACT=2

    # Le script prépare un dossier et un fichier de soumission Slurm pour chaque test
    make_slurm() {
        local folder=$1; local N=$2; local NB=$3; local P=$4; local Q=$5; local NODES=$6
        mkdir -p "$folder"
        
        # Il écrit le fichier HPL.dat indispensable au fonctionnement du benchmark
        cat > "$folder/HPL.dat" << EOD
HPLinpack benchmark
HPL.out
6
1
$N
1
$NB
0
1
$P
$Q
16.0
1
$B_PFACT
1
4
1
2
1
1
1
1
1
1
2
64
0
0
1
8
EOD

        # Le fichier Slurm est configuré pour un accès exclusif aux nœuds de calcul
        cat > "$folder/run_hpl.slurm" << EOS
#!/bin/bash
#SBATCH --account=$ACCOUNT
#SBATCH --constraint=x64cpu
#SBATCH --partition=short
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$((P*Q))
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=job.out

romeo_load_x64cpu_env
module load aocl-4.2.0 openmpi/aocc-4.2.0/4.1.7

# Le script bride le multithreading pour garantir la fiabilité des mesures MPI
export OMP_NUM_THREADS=1

srun ../xhpl
EOS
    }

    # Il boucle sur les différentes valeurs de N pour construire la série S1
    echo "Génération de l'arborescence des tests..."
    for val in 60000 80000 100000 115000 130000; do
        make_slurm "S1_N_$val" $val $B_NB $B_P $B_Q $B_NODES
    done
}

# Cette commande parcourt les journaux de sortie pour extraire les performances en GFlops
extract_results() {
    echo "Récupération des mesures de performance :"
    grep -r "WR" S*/job.out | awk -F'[: ]+' '{print "Cas de test :", $1, "->", $8, "GFlops"}'
}

# Le point d'entrée oriente l'exécution selon l'argument fourni en ligne de commande
case "$1" in
    install) install_standard_hpl ;;
    custom)  compile_custom_hpl ;;
    plan)    generate_test_plan ;;
    results) extract_results ;;
    *) echo "Usage : $0 {install|custom|plan|results}" ;;
esac