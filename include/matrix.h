#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

// typedef int32_t b32;

typedef struct {
    u32 rows;
    u32 cols;
    u8 is_square;
    double* data;
} Matrix;

Matrix* create_matrix(const u32 rows, const u32 cols);
void free_matrix(Matrix* mat);

Matrix* create_random_matrix(unsigned int num_rows, unsigned int num_cols, double min, double max);
Matrix* create_identity_matrix(unsigned int size);
Matrix* create_matrix_fromfilef(FILE* f);
Matrix* copy_matrix(const Matrix* a);

bool check_equality(Matrix* m1, Matrix* m2, double tolerance);

void print_matrix(const Matrix* mat);

Matrix* add_matrices(const Matrix* a, const Matrix* b);
Matrix* subtract_matrices(const Matrix* a, const Matrix* b);

double get_value(const Matrix* mat, const unsigned int row, const unsigned int col);
void set_value(Matrix* mat, const unsigned int row, const unsigned int col, double value);

Matrix* multiply_matrices(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices0(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices1(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices2(const Matrix* a, const Matrix* b);
Matrix* multiply_matrices3(const Matrix* a, const Matrix* b);

double f_norm(const Matrix* mat); // frobenius norm

Matrix* get_col(Matrix* mat, u32 col);
Matrix* get_row(Matrix* mat, u32 row);

void set_all(Matrix* mat, double value);
void set_diag(Matrix* mat, double value);

bool multiply_row_with_scalar_r(Matrix* a, u32 row, double num);
Matrix* multiply_row_with_scalar(const Matrix* a, u32 row, double num);
bool multiply_col_with_scalar_r(Matrix* a, u32 col, double num);
Matrix* multiply_col_with_scalar(const Matrix* a, u32 col, double num);
bool multiply_with_scalar_r(Matrix* a, double num);
Matrix* multiply_with_scalar(const Matrix* a, double num);

bool add_two_rows_r(Matrix* mat, u32 where, u32 row, double multiplier);
Matrix* add_two_rows(Matrix* mat, u32 where, u32 row, double multiplier);

Matrix* remove_col(Matrix* mat, u32 col);
Matrix* remove_row(Matrix* mat, u32 row);

bool swap_col_r(Matrix* mat, u32 col1, u32 col2);
Matrix* swap_col(Matrix* mat, u32 col1, u32 col2);
bool swap_row_r(Matrix* mat, u32 row1, u32 row2);
Matrix* swap_row(Matrix* mat, u32 row1, u32 row2);

Matrix* vcat_matrices(u32 mnum, Matrix** marr);
Matrix* hcat_matrices(u32 mnum, Matrix** marr);

int find_pivotidx(Matrix* mat, u32 col, u32 row);
Matrix* row_echelon_form(Matrix* mat);
int find_pivotmaxidx(Matrix* mat, u32 col, u32 row);
Matrix* reduced_row_echelon_form(Matrix* mat);

#endif