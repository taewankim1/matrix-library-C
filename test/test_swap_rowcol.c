#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 5,-1,1);
    Matrix* a_ = swap_col(a,0,4);
    Matrix* a__ = swap_row(a_,0,3);

    print_matrix(a);
    printf("----\n");
    print_matrix(a_);
    printf("----\n");
    print_matrix(a__);
    // printf("----\n");
    // print_matrix(a___);

    // Free the matrix
    free_matrix(a);
    free_matrix(a_);
    free_matrix(a__);
    // free_matrix(a___);

    return 0;
}