#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 5
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(4, 3,-1,1);
    Matrix* b = create_random_matrix(5, 3,-1,1);
    Matrix* c = create_random_matrix(3, 3,-1,1);

    Matrix* marr[3] = {a,b,c};

    Matrix* vcat_mat = vcat_matrices(3,marr);


    print_matrix(a);
    printf("---\n");
    print_matrix(b);
    printf("---\n");
    print_matrix(c);
    printf("---\n");
    print_matrix(vcat_mat);
    printf("---\n");

    Matrix* d = create_random_matrix(3, 4,-1,1);
    Matrix* e = create_random_matrix(3, 5,-1,1);
    Matrix* f = create_random_matrix(3, 3,-1,1);

    Matrix* marrr[3] = {d,e,f};

    Matrix* hcat_mat = hcat_matrices(3,marrr);

    print_matrix(d);
    printf("---\n");
    print_matrix(e);
    printf("---\n");
    print_matrix(f);
    printf("---\n");
    print_matrix(hcat_mat);

    // Free the matrix
    free_matrix(a);
    free_matrix(b);
    free_matrix(c);
    free_matrix(d);
    free_matrix(e);
    free_matrix(f);
    free_matrix(vcat_mat);
    free_matrix(hcat_mat);

    return 0;
}