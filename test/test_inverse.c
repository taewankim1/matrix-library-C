#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    Matrix* A = create_random_matrix(N, N,-1,1);
    Matrix* Imat = create_identity_matrix(N);
    multiply_with_scalar_r(Imat,1e-6);
    A = add_matrices(A,Imat);
    Matrix* A_inv = inverse(A);
    Matrix* tmp1 = multiply_matrices(A,A_inv);
    Matrix* tmp2 = multiply_matrices(A_inv,A);

    printf("---A---\n");
    print_matrix(A);
    printf("---A_inv---\n");
    print_matrix(A_inv);
    printf("---A * A_inv---\n");
    print_matrix(tmp1);
    printf("---A_inv * A---\n");
    print_matrix(tmp2);

    free_matrix(A);
    free_matrix(Imat);
    free_matrix(A_inv);
    free_matrix(tmp1);
    free_matrix(tmp2);

    return 0;
}