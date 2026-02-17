#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <cblas.h>
#include <string.h>

#define IDX(row, col, lda) ((col) * (lda) + (row))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* * Calcule combien de colonnes chaque processus doit stocker localement.
 * On utilise une distribution cyclique par blocs (1D).
 */
int get_local_num_columns(int N, int NB, int size, int rank) {
    int num_blocks_global = (N + NB - 1) / NB;
    int num_blocks_local = num_blocks_global / size;
    int remainder = num_blocks_global % size;
    if (rank < remainder) num_blocks_local++;
    
    int total_cols = num_blocks_local * NB;
    int last_block_owner = (num_blocks_global - 1) % size;
    
    if (rank == last_block_owner) {
        int last_block_size = N - (num_blocks_global - 1) * NB;
        total_cols -= (NB - last_block_size);
    }
    return total_cols;
}

/* Traduit une colonne globale en son emplacement local (qui possède quoi et où) */
void global_to_local(int global_col, int NB, int size, int *owner, int *local_col) {
    int global_block = global_col / NB;
    int offset_in_block = global_col % NB;
    *owner = global_block % size;
    int local_block = global_block / size;
    *local_col = local_block * NB + offset_in_block;
}

/* Retrouve l'indice global d'une colonne stockée localement */
int local_to_global(int local_col, int NB, int size, int rank) {
    int local_block = local_col / NB;
    int offset_in_block = local_col % NB;
    int global_block = local_block * size + rank;
    return global_block * NB + offset_in_block;
}

/* * --- BLOC DE VÉRIFICATION DE LA FACTORISATION ---
 * * void verify_LU(int N, int local_cols, int NB, int size, int rank, double *A_factorized, double *A_original) {
 * // génère un vecteur X identique sur tous les processus pour servir de base au test
 * double *X = (double*)malloc(N * sizeof(double));
 * for(int i=0; i<N; i++) X[i] = sin(i); 
 *
 * // Étape 1 : calcule le vecteur de référence Y_ref = A * X
 * double *Y_ref_local = (double*)calloc(N, sizeof(double));
 * double *Y_ref = (double*)malloc(N * sizeof(double));
 * for (int j_loc = 0; j_loc < local_cols; j_loc++) {
 * int j_glob = local_to_global(j_loc, NB, size, rank);
 * for (int i = 0; i < N; i++) Y_ref_local[i] += A_original[IDX(i, j_loc, N)] * X[j_glob];
 * }
 * MPI_Allreduce(Y_ref_local, Y_ref, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 *
 * // Étape 2 : calcule Y_calc = L * (U * x) en utilisant la matrice factorisée
 * // Il s'occupe d'abord de la partie triangulaire supérieure (U)
 * double *T_local = (double*)calloc(N, sizeof(double));
 * double *T = (double*)malloc(N * sizeof(double));
 * for (int j_loc = 0; j_loc < local_cols; j_loc++) {
 * int j_glob = local_to_global(j_loc, NB, size, rank);
 * for (int i = 0; i <= j_glob; i++) T_local[i] += A_factorized[IDX(i, j_loc, N)] * X[j_glob];
 * }
 * MPI_Allreduce(T_local, T, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 *
 * // applique ensuite la partie triangulaire inférieure (L) avec une diagonale unité
 * double *Y_calc_local = (double*)calloc(N, sizeof(double));
 * double *Y_calc = (double*)malloc(N * sizeof(double));
 * for (int j_loc = 0; j_loc < local_cols; j_loc++) {
 * int j_glob = local_to_global(j_loc, NB, size, rank);
 * for (int i = j_glob + 1; i < N; i++) Y_calc_local[i] += A_factorized[IDX(i, j_loc, N)] * T[j_glob];
 * }
 * MPI_Allreduce(Y_calc_local, Y_calc, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 * for(int i=0; i<N; i++) Y_calc[i] += T[i];
 *
 * // Étape 3 : compare les deux résultats pour obtenir la norme de l'erreur
 * double error = 0.0, norm_ref = 0.0;
 * for(int i=0; i<N; i++) {
 * double diff = Y_ref[i] - Y_calc[i];
 * error += diff * diff;
 * norm_ref += Y_ref[i] * Y_ref[i];
 * }
 * if(rank == 0) printf("Vérification : Erreur relative ||Ax - LUx|| / ||Ax|| = %e\n", sqrt(error)/sqrt(norm_ref));
 *
 * free(X); free(Y_ref); free(Y_ref_local); free(T); free(T_local); free(Y_calc); free(Y_calc_local);
 * }
 */

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) printf("Usage: %s <N> <NB>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]);
    int NB = atoi(argv[2]);

    int local_cols = get_local_num_columns(N, NB, size, rank);
    int local_rows = N;
    long local_size = (long)local_cols * (long)local_rows;

    if (rank == 0) printf("Lancement HPL Custom | N=%d | NB=%d | Procs=%d\n", N, NB, size);

    double *A = (double*)malloc(local_size * sizeof(double));
    double *panel_buffer = (double*)malloc(local_rows * NB * sizeof(double));
    /* double *A_copy = (double*)malloc(local_size * sizeof(double)); //  réserve de la RAM pour sauvegarder A */

    if (!A || !panel_buffer) {
        printf("Erreur d'allocation sur le rang %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* * Remplissage de la matrice avec des valeurs aléatoires.
     * On booste la diagonale pour être sûr que la factorisation soit stable numériquement.
     */
    srand(rank + 1); 
    for (int j = 0; j < local_cols; ++j) {
        int global_col = local_to_global(j, NB, size, rank);
        for (int i = 0; i < local_rows; ++i) {
            double val = (double)rand() / RAND_MAX;
            if (i == global_col) val += N * 1.0; 
            A[IDX(i, j, local_rows)] = val;
            /* A_copy[IDX(i, j, local_rows)] = val; // duplique les données avant qu'elles ne soient écrasées */
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    /* --- DEBUT DE LA FACTORISATION LU --- */
    for (int k = 0; k < N; k += NB) {
        int jb = MIN(NB, N - k);
        int global_block_idx = k / NB;
        int owner = global_block_idx % size;

        /* Etape 1 : Seul le processus qui possède le bloc courant travaille sur le "panel" */
        if (rank == owner) {
            int local_col_start, dummy;
            global_to_local(k, NB, size, &dummy, &local_col_start);
            
            /* On extrait les colonnes locales vers un buffer temporaire pour faciliter les calculs */
            for (int j = 0; j < jb; ++j) {
                 for(int i=0; i<N; ++i) panel_buffer[IDX(i, j, N)] = A[IDX(i, local_col_start+j, local_rows)];
            }

            /* * Algorithme LU classique sur le bloc vertical (Panel).
             * On utilise BLAS (dscal et dger) pour aller plus vite. 
             */
            for (int j = 0; j < jb; ++j) {
                double pivot = panel_buffer[IDX(k+j, j, N)];
                int len = N - (k + j + 1);
                if(len > 0) cblas_dscal(len, 1.0/pivot, &panel_buffer[IDX(k+j+1, j, N)], 1);
                if (j + 1 < jb) {
                   cblas_dger(CblasColMajor, N - (k + j + 1), jb - (j + 1), -1.0,
                              &panel_buffer[IDX(k+j+1, j, N)], 1,
                              &panel_buffer[IDX(k+j, j+1, N)], N,
                              &panel_buffer[IDX(k+j+1, j+1, N)], N);
                }
            }

            /* On réinjecte le panel factorisé dans la matrice principale */
            for (int j = 0; j < jb; ++j) {
                for(int i=0; i<N; ++i) A[IDX(i, local_col_start+j, local_rows)] = panel_buffer[IDX(i, j, N)];
            }
        }

        /* Etape 2 : Le proprio envoie son panel factorisé à tout le monde pour la mise à jour */
        MPI_Bcast(panel_buffer, N * jb, MPI_DOUBLE, owner, MPI_COMM_WORLD);
        
        /* Etape 3 : Mise à jour du reste de la matrice (Trailing Matrix Update) */
        int blocks_passed_global = global_block_idx + 1;
        int my_blocks_passed = blocks_passed_global / size;
        if (rank < (blocks_passed_global % size)) my_blocks_passed++;
        
        int local_update_start_col = my_blocks_passed * NB;
        int cols_to_update = local_cols - local_update_start_col;

        if (cols_to_update > 0) {
            /* On prépare le bloc U (partie haute) avant le produit de matrice GEMM */
            double *U_block = (double*)malloc(jb * cols_to_update * sizeof(double));
            for(int c=0; c<cols_to_update; ++c) {
                for(int r=0; r<jb; ++r) {
                    U_block[IDX(r, c, jb)] = A[IDX(k+r, local_update_start_col+c, local_rows)];
                }
            }
            
            double *L_ptr = &panel_buffer[IDX(k+jb, 0, N)];
            double *A_ptr = &A[IDX(k+jb, local_update_start_col, local_rows)];
            
            /* * C'est ici que le gros du calcul se passe.
             * On met à jour toutes les colonnes à droite du panel actuel.
             */
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        N - (k + jb), cols_to_update, jb, -1.0,
                        L_ptr, N, U_block, jb, 1.0, A_ptr, local_rows);
            free(U_block);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double time_sec = end_time - start_time;
        /* Formule standard pour compter les opérations flottantes d'une LU */
        double gflops = (2.0 * pow(N,3) / 3.0) / time_sec / 1e9;
        printf("--- FIN DU CALCUL ---\n");
        printf("GFlops : %.4f | Temps : %.4fs\n", gflops, time_sec);
    }
    /* verify_LU(N, local_cols, NB, size, rank, A, A_copy); // lance la procédure de test final */

    free(A);
    free(panel_buffer);
    /* free(A_copy); */
    MPI_Finalize();
    return 0;
}