#ifndef TEST_PCSHELL_H_
#define TEST_PCSHELL_H_

#include <petscksp.h>

/*
 * agmg solves preconditioning system
 */
#if 0
typedef struct agmg_shell_pc
{
    Vec agmg_vec;
} AGMGShellPC;
#endif

// extern PetscErrorCode AGMGShellPCCreate(AGMGShellPC **);
extern PetscErrorCode AGMGShellPCSetup(PC);
extern PetscErrorCode AGMGShellPCApply(PC, Vec, Vec);
// extern PetscErrorCode AGMGShellPCDestroy(PC);

#endif