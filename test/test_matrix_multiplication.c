#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include "../include/matrix.h"

#define N 1024
int main() {
    // Create a 3x3 matrix
    Matrix* a = create_random_matrix(N, N,-1,1);
    Matrix* b = create_random_matrix(N, N,-1,1);
    // Matrix* a = create_random_matrix(3, 5,-5,5);
    // Matrix* b = create_random_matrix(5, 7,-5,5);
    // Matrix* a = create_random_matrix(7, 15,-5,5);
    // Matrix* b = create_random_matrix(15, 13,-5,5);

    // Values for FLOP
    long long Gflop = (long long) N * (long long) N * 2 * (long long) N / 1e9;

    // Mat_0
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    Matrix* c0 = multiply_matrices0(a,b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time_sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Mat_0 : elapsed_time %.9f GFLOPS\n",Gflop / elapsed_time_sec);

    // Mat_1
    clock_gettime(CLOCK_MONOTONIC, &start);
    Matrix* c1 = multiply_matrices1(a,b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Mat_1 : elapsed_time %.9f GFLOPS\n",Gflop / elapsed_time_sec);

    // Mat_2
    clock_gettime(CLOCK_MONOTONIC, &start);
    Matrix* c2 = multiply_matrices2(a,b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Mat_2 : elapsed_time %.9f GFLOPS\n",Gflop / elapsed_time_sec);

    // Mat_3
    clock_gettime(CLOCK_MONOTONIC, &start);
    Matrix* c3 = multiply_matrices3(a,b);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time_sec = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Mat_3 : elapsed_time %.9f GFLOPS\n",Gflop / elapsed_time_sec);

    // Check if they are equal
    bool check_c0_c1 = f_norm(subtract_matrices(c0,c1)) <= 1.e-9 * f_norm(c0);
    bool check_c0_c2 = f_norm(subtract_matrices(c0,c2)) <= 1.e-9 * f_norm(c0);
    bool check_c0_c3 = f_norm(subtract_matrices(c0,c3)) <= 1.e-9 * f_norm(c0);
    printf("Is c1 equal to c0? %s\n",check_c0_c1 ? "true" : "false");
    printf("Is c2 equal to c0? %s\n",check_c0_c2 ? "true" : "false");
    printf("Is c3 equal to c0? %s\n",check_c0_c3 ? "true" : "false");
    // print_matrix(c0);
    // print_matrix(c2);
    // print_matrix(c3);
    // printf("%f, %f\n",f_norm(subtract_matrices(c0,c1)), 1.e-3 * f_norm(c0));
    // printf("\n");
    // print_matrix(a);
    // printf("\n");
    // print_matrix(b);
    // printf("\n");
    // print_matrix(c);
    // printf("\n");

    // Free the matrix
    free_matrix(a);

    return 0;
}