#include "test_pcshell.h"
#include "test_def.h"

#if 0
PetscErrorCode AGMGShellPCCreate(AGMGShellPC **shell)
{
    PetscFunctionBeginUser;

    AGMGShellPC *newctx;
    PetscCall(PetscNew(&newctx));
    newctx->agmg_vec = 0;
    *shell = newctx;

    PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

PetscErrorCode AGMGShellPCApply(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    // PetscCall(VecCopy(x, y));    // test pcshell
    Mat A;

    PetscCall(PCGetOperators(pc, NULL, &A));
    int ijob = 2;    // solve with previous setup
    PetscCall(DefCallAGMG(&A, &x, &y, ijob));

    PetscFunctionReturn(PETSC_SUCCESS);
}

#if 1
PetscErrorCode AGMGShellPCSetup(PC pc)
{
    PetscFunctionBeginUser;

    Mat A;

    PetscCall(PCGetOperators(pc, NULL, &A));
    int ijob = 1;    // setup only
    PetscCall(DefCallAGMG(&A, NULL, NULL, ijob));

    PetscFunctionReturn(PETSC_SUCCESS);
}
#endif

#if 0
PetscErrorCode AGMGShellPCDestroy(PC pc)
{
    PetscFunctionBeginUser;

    AGMGShellPC *shell;
    PetscCall(PCShellGetContext(pc, &shell));
    PetscCall(VecDestroy(&shell->agmg_vec));
    PetscCall(PetscFree(shell));

    PetscFunctionReturn(PETSC_SUCCESS);
}
#endif