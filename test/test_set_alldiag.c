#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 5,-1,1);
    set_all(a,3.2);
    print_matrix(a);

    Matrix* b = create_random_matrix(5, 5,-1,1);
    set_all(b,5.0);
    print_matrix(b);

    // Free the matrix
    free_matrix(a);
    free_matrix(b);

    return 0;
}