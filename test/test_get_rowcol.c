#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 5,-1,1);

    Matrix* a_col = get_col(a,2);
    Matrix* a_row = get_row(a,3);

    print_matrix(a);
    print_matrix(a_col);
    print_matrix(a_row);

    // Free the matrix
    free_matrix(a);
    free_matrix(a_col);
    free_matrix(a_row);

    return 0;
}