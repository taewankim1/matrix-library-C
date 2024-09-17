#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "matrix.h"
#include "mat_utils.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

Matrix* create_matrix(const u32 rows, const u32 cols){
    if (rows <= 0){
        fprintf(stderr, "Invalid rows\n");
        return NULL;
    }
    if (cols <= 0){
        fprintf(stderr, "Invalid cols\n");
        return NULL;
    }
    Matrix* mat = (Matrix*) malloc(sizeof(Matrix));
    if (mat == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix structure.\n");
        return NULL;
    }
    mat->rows = rows;
    mat->cols = cols;
    mat->is_square = (rows == cols) ? 1 : 0;
    mat->data = (double*) calloc(rows*cols,sizeof(double));
    return mat;
}

void free_matrix(Matrix* mat){
    if (mat){
        free(mat->data);
        free(mat);
    }
}

Matrix* create_random_matrix(unsigned int num_rows, unsigned int num_cols, double min, double max){
    Matrix* r = create_matrix(num_rows,num_cols);
    unsigned int i, j;
    for(i = 0; i < num_rows; ++i){
        for(j = 0; j < num_cols; ++j){
            r->data[j + num_cols*i] = rand_interval(min, max);
        }
    }
    return r;
}

Matrix* create_identity_matrix(unsigned int size){
    Matrix* r = create_matrix(size,size);
    for(size_t i = 0; i < size; ++i){
        r->data[i + size*i] = 1.0;
    }
    return r;
}

Matrix* create_matrix_fromfilef(FILE* f){
    unsigned int num_rows = 0, num_cols = 0;
    fscanf(f, "%d", &num_rows);
    fscanf(f, "%d", &num_cols);
    Matrix* r = create_matrix(num_rows, num_cols);
    size_t i, j;
    for(i = 0; i < num_rows; i++) {
        for(j = 0; j < num_cols; j++) {
        fscanf(f, "%lf\t", &r->data[i*num_cols+j]);
        }
    }
    return r;
}

Matrix* copy_matrix(const Matrix* a){
    u32 num_rows = a->rows;
    u32 num_cols = a->cols;
    Matrix* r = create_matrix(num_rows,num_cols);
    memcpy(r->data,a->data,num_rows*num_cols*sizeof(r->data[0]));
    return r;
}

bool check_equality(Matrix* m1, Matrix* m2, double tolerance){
    if (m1->rows != m2->rows) return false;
    if (m1->cols != m2->cols) return false;
    size_t rows = m1->rows;
    size_t cols = m1->cols;
    size_t i, j;
    // double tolerance = 1e-9;
    for(i = 0; i < rows;++i){
        for(j = 0; j < cols;++j){
            if (fabs(m1->data[i*cols+j] - m2->data[i*cols+j]) > tolerance){
                return false;
            };
        }
    }
    return true;
}

void print_matrix(const Matrix* mat){
    int rows = mat->rows;
    int cols = mat->cols;
    for (int i=0; i < rows; i++){
        for (int j=0; j < cols; j++){
            printf("%0.4lf ", mat->data[i*cols + j]);
        }
        printf("\n");
    }
}

double get_value(const Matrix* mat, const unsigned int row, const unsigned int col){
    return mat->data[row*mat->cols + col];
}

void set_value(Matrix* mat, const unsigned int row, const unsigned int col, double value){
    mat->data[row*mat->cols + col] = value;
}

Matrix* add_matrices(const Matrix* a, const Matrix* b){
    if (a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "Error: Matrices dimensions do not match for addition.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    int rows = a->rows;
    int cols = a->cols;
    for (int i=0;i < rows*cols; i++){
        result->data[i] = a->data[i] + b->data[i];
    }
    return result;
}

Matrix* subtract_matrices(const Matrix* a, const Matrix* b){
    if (a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "Error: Matrices dimensions do not match for addition.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    int rows = a->rows;
    int cols = a->cols;
    for (int i=0;i < rows*cols; i++){
        result->data[i] = a->data[i] - b->data[i];
    }
    return result;
}

Matrix* multiply_matrices(const Matrix* a, const Matrix* b){
    return multiply_matrices1(a,b);
}

Matrix* multiply_matrices0(const Matrix* a, const Matrix* b){
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Matrices dimensions do not match for multiplication.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    int a_rows = a->rows;
    int a_cols = a->cols; // equal to b->rows
    int b_cols = b->cols;
    for (int i = 0; i < a_rows; i++) {
        for (int j = 0; j < b_cols; j++) {
            result->data[i*b_cols + j] = 0;
            for (int k = 0; k < a_cols; k++) {
                result->data[i*b_cols + j] += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
            }
        }
    }
    return result;
}

Matrix* multiply_matrices1(const Matrix* a, const Matrix* b){
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Matrices dimensions do not match for multiplication.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    int a_rows = a->rows;
    int a_cols = a->cols; // equal to b->rows
    int b_cols = b->cols;
    double sum;
    for (int i = 0; i < a_rows; i++) {
        for (int j = 0; j < b_cols; j++) {
            sum = 0; // Hoisting
            for (int k = 0; k < a_cols; k++) {
                sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
            }
            result->data[i*b_cols+j] = sum;
            // set_value(result,i,j,sum);
        }
    }
    return result;
    // TODO : Tiling + Blocking
}

Matrix* multiply_matrices2(const Matrix* a, const Matrix* b){
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Matrices dimensions do not match for multiplication.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    int a_rows = a->rows;
    int a_cols = a->cols; // equal to b->rows
    int b_cols = b->cols;
    int i, j;
    for (i = 0; i < a_rows-1; i += 2) {
        for (j = 0; j < b_cols-1; j += 2) {
            double t00 = result->data[i*b_cols+j];
            double t01 = result->data[i*b_cols+j+1];
            double t10 = result->data[(i+1)*b_cols+j];
            double t11 = result->data[(i+1)*b_cols+j+1];
            for (int k = 0; k < a_cols; ++k) {
                t00 += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
                t01 += a->data[i*a_cols+ k] * b->data[k*b_cols + j + 1];
                t10 += a->data[(i+1)*a_cols+ k] * b->data[k*b_cols + j];
                t11 += a->data[(i+1)*a_cols+ k] * b->data[k*b_cols + j + 1];
            }
            result->data[i*b_cols+j] = t00;
            result->data[i*b_cols+j + 1] = t01;
            result->data[(i+1)*b_cols+j] = t10;
            result->data[(i+1)*b_cols+j+1] = t11;
        }
    }
    int i_clean = i;
    int j_clean = j;
    double sum;
    for (i = i_clean; i < a_rows; ++i) {
        for (j = 0; j < b_cols - 1; ++j) {
            sum = 0; // Hoisting
            for (int k = 0; k < a_cols; k++) {
                sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
            }
            result->data[i*b_cols+j] = sum;
        }
    }
    for (i = 0; i < a_rows - 1; ++i) {
        for (j = j_clean; j < b_cols; ++j) {
            sum = 0; // Hoisting
            for (int k = 0; k < a_cols; k++) {
                sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
            }
            result->data[i*b_cols+j] = sum;
        }
    }
    sum = 0;
    i = i_clean;
    j = j_clean;
    if (i < a_rows && j < b_cols){
        for (int k = 0; k < a_cols; k++) {
            sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
        }
    }
    result->data[i*b_cols+j] = sum;
    return result;
    // TODO : Blocking
}

Matrix* multiply_matrices3(const Matrix* a, const Matrix* b){
    if (a->cols != b->rows) {
        fprintf(stderr, "Error: Matrices dimensions do not match for multiplication.\n");
        return NULL;
    }

    Matrix* result = create_matrix(a->rows,b->cols);
    size_t a_rows = a->rows;
    size_t a_cols = a->cols; // equal to b->rows
    size_t b_cols = b->cols;
    size_t blocksize = 32;

    assert(a_rows % blocksize == 0);
    assert(a_cols % blocksize == 0);
    assert(b_cols % blocksize == 0);

    // size_t ii, jj;
    for (size_t ii = 0; ii < a_rows; ii += blocksize) {
        for (size_t jj = 0; jj < b_cols; jj += blocksize) {
            for (size_t kk = 0; kk < a_cols; kk += blocksize) {
                size_t stop_i  = MIN(ii + blocksize, a_rows);
                size_t stop_j  = MIN(jj + blocksize, b_cols);
                size_t stop_k  = MIN(kk + blocksize, a_cols);
                for (size_t i = ii; i < stop_i; i += 2) {
                    for (size_t j = jj; j < stop_j; j += 2) {

                        // unrolling
                        double t00 = result->data[i*b_cols+j];
                        double t01 = result->data[i*b_cols+j+1];
                        double t10 = result->data[(i+1)*b_cols+j];
                        double t11 = result->data[(i+1)*b_cols+j+1];
                        for (size_t k = kk; k < stop_k; ++k) {
                            t00 += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
                            t01 += a->data[i*a_cols+ k] * b->data[k*b_cols + j + 1];
                            t10 += a->data[(i+1)*a_cols+ k] * b->data[k*b_cols + j];
                            t11 += a->data[(i+1)*a_cols+ k] * b->data[k*b_cols + j + 1];
                        }
                        result->data[i*b_cols+j] = t00;
                        result->data[i*b_cols+j + 1] = t01;
                        result->data[(i+1)*b_cols+j] = t10;
                        result->data[(i+1)*b_cols+j+1] = t11;
                    }
                }
            }
        }
    }
    // printf("i: %zu, j: %zu\n",ii,jj);
    // int i_clean = i;
    // int j_clean = j;
    // double sum;
    // for (i = i_clean; i < a_rows; ++i) {
    //     for (j = 0; j < b_cols - 1; ++j) {
    //         sum = 0; // Hoisting
    //         for (int k = 0; k < a_cols; k++) {
    //             sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
    //         }
    //         result->data[i*b_cols+j] = sum;
    //     }
    // }
    // for (i = 0; i < a_rows - 1; ++i) {
    //     for (j = j_clean; j < b_cols; ++j) {
    //         sum = 0; // Hoisting
    //         for (int k = 0; k < a_cols; k++) {
    //             sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
    //         }
    //         result->data[i*b_cols+j] = sum;
    //     }
    // }
    // sum = 0;
    // i = i_clean;
    // j = j_clean;
    // if (i < a_rows && j < b_cols){
    //     for (int k = 0; k < a_cols; k++) {
    //         sum += a->data[i*a_cols+ k] * b->data[k*b_cols + j];
    //     }
    // }
    // result->data[i*b_cols+j] = sum;
    return result;
    // TODO : Blocking
}

double f_norm(const Matrix* mat){
    double sum = 0.0;
    size_t num_rows = mat->rows;
    size_t num_cols = mat->cols;
    for (size_t i = 0; i < num_rows; ++i) {    // compute column sum
        for (size_t j = 0; j < num_cols; ++j) {    // for each column
            sum += mat->data[i * num_cols + j] * mat->data[i*num_cols + j];
        }
    }
    return sqrt(sum);
}

// double one_norm(const Matrix& A) {
//   double sum = 0.0;
//   for (size_t j = 0; j < A.num_cols(); ++j) {    // for each column
//     double tmp = 0.0;
//     for (size_t i = 0; i < A.num_rows(); ++i) {    // compute column sum
//       tmp += std::abs(A(i, j));
//     }
//     sum = std::max(sum, tmp);
//   }
//   return sum;
// }


// double inf_norm(const Matrix& A) {
//   double sum = 0.0;
//   for (size_t i = 0; i < A.num_rows(); ++i) {    // for each row
//     double tmp = 0.0;
//     for (size_t j = 0; j < A.num_cols(); ++j) {    // compute row sum
//       tmp += std::abs(A(i, j));
//     }
//     sum = std::max(sum, tmp);
//   }
//   return sum;
// }

Matrix* get_col(Matrix* mat, u32 col){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (col >= mat_cols) {
        fprintf(stderr, "col >= mat->cols \n");
        return NULL;
    }
    Matrix* r = create_matrix(mat_rows,1);
    for(u32 j=0;j<mat_rows;++j){
        r->data[j] = mat->data[mat_cols*j + col];
    }
    return r;
}

Matrix* get_row(Matrix* mat, u32 row){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (row >= mat_rows) {
        fprintf(stderr, "row >= mat->rows \n");
        return NULL;
    }
    Matrix* r = create_matrix(1,mat_cols);
    memcpy(r->data,&mat->data[mat_cols*row],mat_cols * sizeof(r->data[0]));
    return r;
}

void set_all(Matrix* mat, double value){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    for(u32 i = 0;i<mat_rows;++i){
        for(u32 j = 0;j<mat_cols;++j){
            mat->data[mat_cols*i+j] = value;
        }
    }
}

void set_diag(Matrix* mat, double value){
    if (!mat->is_square) {
        fprintf(stderr, "matrix is not square \n");
    }
    u32 mat_rows = mat->rows;
    for(u32 i = 0;i<mat_rows;++i){
        mat->data[mat_rows*i+i] = value;
    }
}

bool multiply_row_with_scalar_r(Matrix* mat, u32 row, double num){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (row >= mat_rows){
        fprintf(stderr, "row >= mat->rows \n");
        return false;
    }
    for(u32 i = 0;i < mat_cols;++i){
        mat->data[mat_cols*row+i] *= num;
    }
    return true;
}
Matrix* multiply_row_with_scalar(const Matrix* mat, u32 row, double num){
    Matrix* r = copy_matrix(mat);
    if (!multiply_row_with_scalar_r(r,row,num)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

bool multiply_col_with_scalar_r(Matrix* mat, u32 col, double num){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (col >= mat_cols){
        fprintf(stderr, "col >= mat->cols \n");
        return false;
    }
    for(u32 i = 0;i < mat_rows;++i){
        mat->data[mat_cols*i+col] *= num;
    }
    return true;
}

Matrix* multiply_col_with_scalar(const Matrix* mat, u32 col, double num){
    Matrix* r = copy_matrix(mat);
    if (!multiply_col_with_scalar_r(r,col,num)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

bool multiply_with_scalar_r(Matrix* mat, double num){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    for(u32 i = 0;i < mat_rows;++i){
        for(u32 j = 0;j < mat_cols;++j){
        mat->data[mat_cols*i+j] *= num;
        }
    }
    return true;
}

Matrix* multiply_with_scalar(const Matrix* mat, double num){
    Matrix* r = copy_matrix(mat);
    if (!multiply_with_scalar_r(r,num)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

bool add_two_rows_r(Matrix* mat, u32 where, u32 row, double multiplier){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (where >= mat_rows || row >= mat_rows){
        fprintf(stderr, "row >= mat->rows \n");
        return false;
    }
    for(u32 i = 0;i < mat_cols;++i){
        mat->data[mat_cols*where+i] += multiplier * mat->data[mat_cols*row+i];
    }
    return true;
}

Matrix* add_two_rows(Matrix* mat, u32 where, u32 row, double multiplier){
    Matrix* r = copy_matrix(mat);
    if (!add_two_rows_r(r,where,row,multiplier)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

Matrix* remove_col(Matrix* mat, u32 col){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (col >= mat_cols){
        fprintf(stderr, "col >= mat->cols \n");
        return NULL;
    }
    if (mat_cols == 1){
        fprintf(stderr, "Only 1 column is left \n");
        return NULL;
    }
    Matrix* r = create_matrix(mat_rows,mat_cols-1);
    u32 i,j,k;
    for(i = 0;i < mat_rows;++i){
        k = 0;
        for(j = 0;j < mat_cols;++j){
            if (j == col) continue;
            // printf("j: %d, k: %d, val: %lf\n",j,k,mat->data[mat_cols*i+j]);
            r->data[(mat_cols-1)*i+(k++)] = mat->data[mat_cols*i+j];
        }
    }
    return r;
}

Matrix* remove_row(Matrix* mat, u32 row){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (row >= mat_rows){
        fprintf(stderr, "col >= mat->cols \n");
        return NULL;
    }
    if (mat_rows == 1){
        fprintf(stderr, "Only 1 row is left \n");
        return NULL;
    }
    Matrix* r = create_matrix(mat_rows-1,mat_cols);
    u32 i,j,k;
    for(i = 0, k=0;i < mat_rows;++i){
        if (i == row) continue;
        for(j = 0;j < mat_cols;++j){
            r->data[mat_cols*k+j] = mat->data[mat_cols*i+j];
        }
        k++;
    }
    return r;
}

bool swap_col_r(Matrix* mat, u32 col1, u32 col2){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (col1 >= mat_cols || col2 >= mat_cols){
        fprintf(stderr, "col >= mat->cols \n");
        return false;
    }
    double tmp;
    for(u32 i = 0;i < mat_rows;++i){
        tmp = mat->data[mat_cols*i+col1];
        mat->data[mat_cols*i+col1] = mat->data[mat_cols*i+col2];
        mat->data[mat_cols*i+col2] = tmp;
    }
    return true;
}

Matrix* swap_col(Matrix* mat, u32 col1, u32 col2){
    Matrix* r = copy_matrix(mat);
    if (!swap_col_r(r,col1,col2)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

bool swap_row_r(Matrix* mat, u32 row1, u32 row2){
    u32 mat_rows = mat->rows;
    u32 mat_cols = mat->cols;
    if (row1 >= mat_rows || row2 >= mat_rows){
        fprintf(stderr, "col >= mat->cols \n");
        return false;
    }
    double tmp;
    for(u32 i = 0;i < mat_cols;++i){
        tmp = mat->data[mat_cols*row1+i];
        mat->data[mat_cols*row1+i] = mat->data[mat_cols*row2+i];
        mat->data[mat_cols*row2+i] = tmp;
    }
    return true;
}

Matrix* swap_row(Matrix* mat, u32 row1, u32 row2){
    Matrix* r = copy_matrix(mat);
    if (!swap_row_r(r,row1,row2)){
        free_matrix(r);
        return NULL;
    }
    return r;
}

Matrix* vcat_matrices(u32 mnum, Matrix** marr){
    if (mnum == 0){
        return NULL;
    }
    if (mnum == 1){
        return copy_matrix(marr[0]);
    }
    u32 i,j,k,offset;
    u32 rows, cols;
    rows = marr[0]->rows;
    cols = marr[0]->cols;
    for(k = 1; k < mnum;++k){
        if (marr[k] == NULL){
            fprintf(stderr, "INCONSISTENT ARRAY \n");
            return NULL;
        }
        if (cols != marr[k]->cols){
            fprintf(stderr, "CANNOT CONCATENATE V \n");
            return NULL;
        }
        rows += marr[k]->rows;
    }
    Matrix* r = create_matrix(rows, cols);
    for(j = 0; j < cols; ++j) {
        offset = 0;
        k = 0;
        for(i = 0;i<rows;++i){
            if (i - offset == marr[k]->rows){
                offset += marr[k]->rows;
                k++;
            }
            r->data[cols*i+j] = marr[k]->data[cols*(i-offset)+j];
        }
    }
    return r;
}

Matrix* hcat_matrices(u32 mnum, Matrix** marr){
    if (mnum == 0){
        return NULL;
    }
    if (mnum == 1){
        return copy_matrix(marr[0]);
    }
    u32 i,j,k,offset;
    u32 rows, cols;
    rows = marr[0]->rows;
    cols = marr[0]->cols;
    for(k = 1; k < mnum;++k){
        if (marr[k] == NULL){
            fprintf(stderr, "INCONSISTENT ARRAY \n");
            return NULL;
        }
        if (rows != marr[k]->rows){
            fprintf(stderr, "CANNOT CONCATENATE V \n");
            return NULL;
        }
        cols += marr[k]->cols;
    }
    Matrix* r = create_matrix(rows, cols);
    for(i = 0;i<rows;++i){
        offset = 0;
        k = 0;
        for(j = 0; j < cols; ++j) {
            if (j - offset == marr[k]->cols){
                offset += marr[k]->cols;
                k++;
            }
            r->data[cols*i+j] = marr[k]->data[marr[k]->cols*i+j-offset];
        }
    }
    return r;
}