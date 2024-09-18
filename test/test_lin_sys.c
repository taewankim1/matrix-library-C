#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 3
int main() {
    Matrix* A = create_random_matrix(N, N,-1,1);
    Matrix* A_diag = create_identity_matrix(N);
    multiply_with_scalar_r(A_diag,3.0);
    Matrix* b = create_random_matrix(N,1,-1,1);

    Matrix* x = solve_ls(A,b);
    Matrix* x_diag = solve_ls(A_diag,b);

    printf("---A---\n");
    print_matrix(A);
    printf("---b---\n");
    print_matrix(b);
    printf("---x---\n");
    print_matrix(x);
    printf("---A_diag---\n");
    print_matrix(A_diag);
    printf("---x_diag---\n");
    print_matrix(x_diag);

    free_matrix(A);
    free_matrix(b);
    free_matrix(A_diag);
    free_matrix(x);
    free_matrix(x_diag);

    return 0;
}