#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(3, 5,-1,1);
    Matrix* b = transpose(a);
    Matrix* c = transpose(b);

    print_matrix(a);
    printf("--\n");
    print_matrix(b);
    printf("--\n");
    print_matrix(subtract_matrices(a,c));

    // Free the matrix
    free_matrix(a);
    free_matrix(b);
    free_matrix(c);

    return 0;
}