#ifndef TEST_DEF_H_
#define TEST_DEF_H_

#include <petscksp.h>

extern PetscErrorCode DefCommandProcess(int, char **, char **, char **, char **);
extern PetscErrorCode DefMatAssemble(const char *, Mat *);
extern PetscErrorCode DefVecAssemble(const char *, Vec *, Vec *);
extern PetscErrorCode DefRHSProcess(const char *, int *, double **, double **);
extern PetscErrorCode DefMatProcess(const char *, int *, int *, int **, int **, int **, double **, double **);
extern PetscErrorCode DefCallAGMG(Mat *, Vec *, Vec *, int);
void zagmg_(int *, double *, int *, int *, double *, double *, int *, int *, int *, int *, double *);
#if 0
extern "C"
{
    void zagmg_(int *, double *, int *, int *, double *, double *, int *, int *, int *, int *, double *);
}
#endif

#endif