#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    // Matrix* a = create_identity_matrix(N);
    char* path = "./matrix_tmp1.data";
    FILE* f = fopen(path, "rb"); // "rb" <- reading it in binary mode
    Matrix* a = create_matrix_fromfilef(f);

    path = "./matrix_tmp2.data";
    f = fopen(path, "rb"); // "rb" <- reading it in binary mode
    Matrix* b = create_matrix_fromfilef(f);

    print_matrix(a);
    print_matrix(b);

    printf("%s\n",check_equality(a,b,1e-9) ? "true" : "false");

    free_matrix(a);
    free_matrix(b);

    return 0;
}