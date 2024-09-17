#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 5,-1,1);
    // Matrix* a = create_identity_matrix(5);

    Matrix* a_ = add_two_rows(a,2,3,2);

    print_matrix(a);
    printf("----\n");
    print_matrix(a_);

    // Free the matrix
    free_matrix(a);
    free_matrix(a_);

    return 0;
}