#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(N, N,-1,1);

    QR* qr = solve_qr(a);

    print_matrix(a);
    printf("----\n");
    print_matrix(qr->Q);
    printf("----\n");
    print_matrix(qr->R);
    printf("----\n");
    Matrix* QmulR = multiply_matrices(qr->Q,qr->R);
    print_matrix(QmulR);
    printf("----\n");
    Matrix* error = subtract_matrices(a,QmulR);
    print_matrix(error);
    printf("error: %lf\n",f_norm(error));


    // Free the matrix
    free_matrix(a);
    free_matrix(QmulR);
    free_matrix(error);
    free_QR(qr);

    return 0;
}