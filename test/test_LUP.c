#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(N, N,-1,1);

    LUP* lup = solve_lup(a);

    print_matrix(a);
    printf("----\n");
    print_matrix(lup->L);
    printf("----\n");
    print_matrix(lup->U);
    printf("----\n");
    print_matrix(lup->P);
    printf("----\n");
    Matrix* LU = multiply_matrices(lup->L,lup->U);
    Matrix* PtLU = multiply_matrices(transpose(lup->P),LU);
    print_matrix(subtract_matrices(PtLU,a));
    printf("error: %lf\n",f_norm(subtract_matrices(PtLU,a)));


    // Free the matrix
    free_matrix(a);
    free_LUP(lup);
    free_matrix(LU);
    free_matrix(PtLU);

    return 0;
}