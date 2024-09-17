#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 108
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(N, N,-1,1);
    Matrix* b = create_identity_matrix(N);
    Matrix* c = multiply_matrices(a,b);

    // Check if they are equal
    bool check_by_f_norm = f_norm(subtract_matrices(a,c)) <= 1.e-9 * f_norm(a);
    bool check_by_c_equality = check_equality(a,c,1e-9);
    printf("Is a equal to c? %s\n",check_by_f_norm ? "true" : "false");
    printf("Is a equal to c? %s\n",check_by_c_equality ? "true" : "false");

    // Free the matrix
    free_matrix(a);
    free_matrix(b);
    free_matrix(c);

    return 0;
}