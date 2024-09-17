#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 5,-1,1);

    Matrix* a_ = multiply_row_with_scalar(a,2,3);
    Matrix* a__ = multiply_col_with_scalar(a_,2,3);
    Matrix* a___ = multiply_with_scalar(a__,3.0);

    print_matrix(a);
    printf("----\n");
    print_matrix(a_);
    printf("----\n");
    print_matrix(a__);
    printf("----\n");
    print_matrix(a___);

    // Free the matrix
    free_matrix(a);
    free_matrix(a_);
    free_matrix(a__);
    free_matrix(a___);

    return 0;
}