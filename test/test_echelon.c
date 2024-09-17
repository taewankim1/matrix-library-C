#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a1 = create_random_matrix(3, 6,-1,1);
    Matrix* a2 = create_random_matrix(6, 3,-1,1);
    Matrix* a3 = create_random_matrix(3, 3,-1,1);
    Matrix* a4 = create_identity_matrix(3);
    set_value(a4,0,0,0);
    Matrix* a5 = create_matrix(3, 3);
    // multiply_with_scalar_r(a5,0.0);
    // Matrix* a = create_identity_matrix(5);

    Matrix* b1 = row_echelon_form(a1);
    Matrix* b2 = row_echelon_form(a2);
    Matrix* b3 = row_echelon_form(a3);
    Matrix* b4 = row_echelon_form(a4);
    Matrix* b5 = row_echelon_form(a5);

    print_matrix(a1);
    printf("----\n");
    print_matrix(b1);
    printf("----\n");
    print_matrix(a2);
    printf("----\n");
    print_matrix(b2);
    printf("----\n");
    print_matrix(a3);
    printf("----\n");
    print_matrix(b3);
    printf("----\n");
    print_matrix(a4);
    printf("----\n");
    print_matrix(b4);
    printf("----\n");
    print_matrix(a5);
    printf("----\n");
    print_matrix(b5);
    printf("----\n");

    // Free the matrix
    free_matrix(a1);
    free_matrix(b1);
    free_matrix(a2);
    free_matrix(b2);
    free_matrix(a3);
    free_matrix(b3);
    free_matrix(a4);
    free_matrix(b4);
    free_matrix(a5);
    free_matrix(b5);

    return 0;
}