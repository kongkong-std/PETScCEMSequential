#include "test_def.h"

PetscErrorCode DefCommandProcess(int argc, char **argv,
                                 char **path_mat, char **path_rhs, char **path_pc)
{
    PetscFunctionBeginUser;
#if 0
    printf("argc = %d\n", argc);
    for (int index = 0; index < argc; ++index)
    {
        puts(argv[index]);
    }
#endif

    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mat_path", argv[index]))
        {
            // file path to csr matrix
            *path_mat = argv[index + 1];
            puts(*path_mat);
        }
        else if (strstr("-rhs_path", argv[index]))
        {
            // file path to right-hand side vector
            *path_rhs = argv[index + 1];
            puts(*path_rhs);
        }
        else if (strstr("-pc_path", argv[index]))
        {
            // file path to preconditioners
            *path_pc = argv[index + 1];
            puts(*path_pc);
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DefMatAssemble(const char *path_mat, Mat *petsc_mat_A)
{
    PetscFunctionBeginUser;

    // file io get coefficient matrix
    int n_mat = 0, nnz_mat = 0;
    int *row_idx = NULL, *row_ptr = NULL, *column_idx = NULL;
    double *mat_value_re = NULL, *mat_value_im = NULL;
    DefMatProcess(path_mat, &n_mat, &nnz_mat, &row_idx, &row_ptr, &column_idx, &mat_value_re, &mat_value_im);

    // petsc matrix mat
    PetscCall(MatCreate(PETSC_COMM_SELF, petsc_mat_A));
    PetscCall(MatSetSizes(*petsc_mat_A, PETSC_DECIDE, PETSC_DECIDE, n_mat, n_mat));
    PetscCall(MatSetType(*petsc_mat_A, MATAIJ));
    PetscCall(MatSetUp(*petsc_mat_A));

    // assigning values for matrix
    for (int index = 0; index < n_mat; ++index)
    {
        int index_start = row_ptr[index];
        int index_end = row_ptr[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar mat_tmp = mat_value_re[index_j] + mat_value_im[index_j] * PETSC_i;
            PetscCall(MatSetValue(*petsc_mat_A, index, column_idx[index_j], mat_tmp, INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(*petsc_mat_A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*petsc_mat_A, MAT_FINAL_ASSEMBLY));

    // free memory
    free(row_idx);
    free(row_ptr);
    free(column_idx);
    free(mat_value_re);
    free(mat_value_im);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DefVecAssemble(const char *path_rhs, Vec *vec_rhs, Vec *vec_sol)
{
    PetscFunctionBeginUser;

    // file io get rhs vector
    int n_vec = 0;
    double *vec_value_re = NULL, *vec_value_im = NULL;
    DefRHSProcess(path_rhs, &n_vec, &vec_value_re, &vec_value_im);

    // petsc vector rhs, sol
    PetscCall(VecCreate(PETSC_COMM_SELF, vec_rhs));
    PetscCall(VecSetType(*vec_rhs, VECSEQ));
    PetscCall(VecSetSizes(*vec_rhs, PETSC_DECIDE, n_vec));
    PetscCall(VecDuplicate(*vec_rhs, vec_sol));

    // assigning values for rhs vector
    /*
     * petsc_vec_rhs[i] = vec_value_re[i] + vec_value_im[i] * PETSC_i
     */
    for (int index = 0; index < n_vec; ++index)
    {
        PetscScalar vec_tmp = vec_value_re[index] + vec_value_im[index] * PETSC_i;
        PetscCall(VecSetValues(*vec_rhs, 1, &index, &vec_tmp, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(*vec_rhs));
    PetscCall(VecAssemblyEnd(*vec_rhs));

    // free memory
    free(vec_value_re);
    free(vec_value_im);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DefCallAGMG(Mat *linsys_mat, Vec *linsys_rhs, Vec *linsys_sol, int ijob)
{
    PetscFunctionBeginUser;
    // zagmg function
    /*
     * void zagmg_( int * n, double * a, int * ja, int * ia, double * f, double * x,
     *              int * ijob, int * iprint, int * nrest, int * maxit, double * tol)
     * n: size of matrix
     * a: non-zero elements in csr format
     * ja: column indices in csr format
     * ia: row pointer in csr format
     * f: right-hand side vector
     * x: numerical solution
     * ijob: agmg parameter, 0 is recommended
     * iprint: agmg parameter, 6 is recommended
     * nrest: number of restart of cgs
     * maxit: maximum iterations
     * tol: tolerance, relative
     */
    int n = 0, nnz = 0;
    double *a = NULL, *f = NULL, *x = NULL;
    int *ja = NULL, *ia = NULL;

    /*
     * iprint = 6, standard output
     * iprint = 0, not print on screen
     */
    // int ijob = 2;    // use previous setup
    int iprint = 0, nrest = 50, maxit = 5;
    double tol = 0.1;

    // get size
    PetscCall(MatGetSize(*linsys_mat, &n, NULL));

    // get ia, ja
    PetscBool done = PETSC_TRUE;
    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    PetscCall(MatGetRowIJ(*linsys_mat, 0, PETSC_FALSE, PETSC_TRUE, NULL, &csr_ia, &csr_ja, &done));

    // memory allocation and assigning values for ia, ja
    if ((ia = (int *)malloc((n + 1) * sizeof(int))) == NULL ||
        (ja = (int *)malloc(csr_ia[n] * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < n + 1; ++index)
    {
        ia[index] = csr_ia[index];
    }
    for (int index = 0; index < csr_ia[n]; ++index)
    {
        ja[index] = csr_ja[index];
    }

    // memory allocation for a, x, f
    /*
     * nnz = ia[n]
     * size a = nnz x 2, complex type value
     * size x = n x 2, complex type value
     * size f = n x 2, complex type value
     */
    nnz = ia[n];
    if ((a = (double *)malloc(nnz * 2 * sizeof(double))) == NULL ||
        (x = (double *)malloc(n * 2 * sizeof(double))) == NULL ||
        (f = (double *)malloc(n * 2 * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    // assigning value for a, x, f
    for (int index = 0; index < 2 * n; ++index)
    {
        x[index] = 0;
        f[index] = 0;
    }
    for (int index = 0; index < 2 * nnz; ++index)
    {
        a[index] = 0;
    }

    // non-zero element
    for (int index = 0; index < n; ++index)
    {
        int index_start = ia[index];
        int index_end = ia[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(*linsys_mat, index, ja[index_j], &val_tmp));
            a[2 * index_j] = PetscRealPart(val_tmp);
            a[2 * index_j + 1] = PetscImaginaryPart(val_tmp);
        }
    }

    // right-hand side vector
    if (linsys_rhs != NULL)
    {
        for (int index = 0; index < n; ++index)
        {
            PetscScalar val_tmp;
            PetscCall(VecGetValues(*linsys_rhs, 1, &index, &val_tmp));
            f[2 * index] = PetscRealPart(val_tmp);
            f[2 * index + 1] = PetscImaginaryPart(val_tmp);
        }
    }
#if 0 // check mat information
    puts("size of matrix:");
    printf("n= %d\n", n);
    puts("\ncsr format:");
    printf("here lists row pointer in csr format...\n");
    for (int index = 0; index < n + 1; ++index)
    {
        printf("%d\n", ia[index]);
    }
    printf("here lists column indices in csr format...\n");
    for (int index = 0; index < csr_ia[n]; ++index)
    {
        printf("%d\n", ja[index]);
    }
    printf("here lists non-zero element...\n");
    for (int index = 0; index < nnz; ++index)
    {
        printf("%021.14le\t%021.14le\n", a[2 * index], a[2 * index + 1]);
    }
    printf("here lists right-hand side vector...\n");
    for (int index = 0; index < n; ++index)
    {
        printf("%021.14le\t%021.14le\n", f[2 * index], f[2 * index + 1]);
    }
#endif

    // updating 1-base
    for (int index = 0; index < n + 1; ++index)
    {
        ia[index] += 1;
    }
    for (int index = 0; index < nnz; ++index)
    {
        ja[index] += 1;
    }

    // call agmg solve linear system
    zagmg_(&n, a, ja, ia, f, x, &ijob, &iprint, &nrest, &maxit, &tol);

    // assigning agmg numerical solution to linsys_sol
    /*
     * linsys_sol[i] = x[2 * i] + x[2 * i + 1] * PETSC_i, i = 0, 1, ..., n - 1
     */
    if (linsys_sol != NULL)
    {
        for (int index = 0; index < n; ++index)
        {
            PetscScalar vec_tmp = x[2 * index] + x[2 * index + 1] * PETSC_i;
            PetscCall(VecSetValues(*linsys_sol, 1, &index, &vec_tmp, INSERT_VALUES));
        }
        PetscCall(VecAssemblyBegin(*linsys_sol));
        PetscCall(VecAssemblyEnd(*linsys_sol));
    }

    // free memory
    free(ia);
    free(ja);
    free(a);
    free(x);
    free(f);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DefMatProcess(const char *path_mat, int *size, int *nnz, int **row_idx, int **row_ptr,
                             int **column_idx, double **value_re, double **value_im)
{
    PetscFunctionBeginUser;
    FILE *fp = NULL;
    if ((fp = fopen(path_mat, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file \"%s\"!\n", path_mat);
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d%d", size, nnz);
    if ((*row_idx = (int *)malloc(*size * sizeof(int))) == NULL ||
        (*row_ptr = (int *)malloc((*size + 1) * sizeof(int))) == NULL ||
        (*column_idx = (int *)malloc(*nnz * sizeof(int))) == NULL ||
        (*value_re = (double *)malloc(*nnz * sizeof(double))) == NULL ||
        (*value_im = (double *)malloc(*nnz * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < *size; ++index)
    {
        fscanf(fp, "%d", *row_idx + index);
    }
    for (int index = 0; index < *size + 1; ++index)
    {
        fscanf(fp, "%d", *row_ptr + index);
    }
    for (int index = 0; index < *nnz; ++index)
    {
        fscanf(fp, "%d", *column_idx + index);
    }
    for (int index = 0; index < *nnz; ++index)
    {
        fscanf(fp, "%lf%lf", *value_re + index, *value_im + index);
    }
    fclose(fp);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DefRHSProcess(const char *path_rhs, int *size, double **value_re, double **value_im)
{
    PetscFunctionBeginUser;

    *value_re = NULL;
    *value_im = NULL;
    FILE *fp = NULL;
    if ((fp = fopen(path_rhs, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file \"%s\"!\n", path_rhs);
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d", size);
    if ((*value_re = (double *)malloc(*size * sizeof(double))) == NULL ||
        (*value_im = (double *)malloc(*size * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < *size; ++index)
    {
        fscanf(fp, "%*d%lf%lf", *value_re + index, *value_im + index);
    }
    fclose(fp);

#if 0 // check file io
    printf("address of value_re = %p, value_im = %p\n", *value_re, *value_im);
#endif
    PetscFunctionReturn(PETSC_SUCCESS);
}
