/*
 * test petsc with agmg lib
 */

#include <petscksp.h>
#include "test_def.h"
#include "test_pcshell.h"

int main(int argc, char **argv)
{
    PetscFunctionBeginUser;

    char *path_mat = NULL, *path_rhs = NULL, *path_pc = NULL;
    PetscCall(DefCommandProcess(argc, argv, &path_mat, &path_rhs, &path_pc));

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    // assemble linear system
    /*
     * coefficient matrix, right-hand side vector, preconditioner
     */
    Vec petsc_vec_rhs, petsc_vec_sol; // right-hand side vector, solution vector
    PetscCall(DefVecAssemble(path_rhs, &petsc_vec_rhs, &petsc_vec_sol));

    Mat petsc_mat_A, petsc_mat_pc; // coefficient matrix
    PetscCall(DefMatAssemble(path_mat, &petsc_mat_A));
    PetscCall(DefMatAssemble(path_pc, &petsc_mat_pc));

#if 0 // test complex value in petsc
    PetscScalar a = 3.01 + 1.02 * PETSC_i;
    PetscReal a_re = PetscRealPart(a), a_im = PetscImaginaryPart(a);
    printf("a_re = %.2lf, a_im = %.2lf\n", a_re, a_im);
    printf("real part of a is %.2lf, imaginary part is %.2lf\n", PetscRealPart(a), PetscImaginaryPart(a));
#endif

#if 0 // vec viewer, check vector assembling
    PetscCall(VecView(*petsc_vec_rhs, PETSC_VIEWER_STDOUT_SELF));
    puts(">>>>>>>>before calling agmg, numerical solution...");
    PetscCall(VecView(petsc_vec_sol, PETSC_VIEWER_STDOUT_SELF));
#endif

#if 0 // mat viewer, check matrix assembling
    PetscCall(MatView(petsc_mat_A, PETSC_VIEWER_STDOUT_SELF));
    puts( "============" );
    PetscCall(MatView(petsc_mat_pc, PETSC_VIEWER_STDOUT_SELF));
#endif

    // call agmg
    /*
     * input: mat, rhs, sol
     * output: sol
     */
    // PetscCall(DefCallAGMG(&petsc_mat_A, &petsc_vec_rhs, &petsc_vec_sol));
#if 0 // check numerical solution via agmg
    puts(">>>>>>>>after calling agmg, numerical solution...");
    PetscCall(VecView(petsc_vec_sol, PETSC_VIEWER_STDOUT_SELF));
#endif

    // solving linear system with ksp
    KSP ksp;
    PC pc;

    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOperators(ksp, petsc_mat_A, petsc_mat_pc));
    PetscCall(KSPGetPC(ksp, &pc));

    // user-defined preconditioner
    PetscBool user_defined_pc = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-user_defined_pc", &user_defined_pc, NULL));
    if (user_defined_pc)
    {
        PetscCall(PCSetType(pc, PCSHELL));
        PetscCall(PCShellSetSetUp(pc, AGMGShellPCSetup));
        PetscCall(PCShellSetApply(pc, AGMGShellPCApply));
        PetscCall(PCShellSetContext(pc, NULL));
    }

    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp, petsc_vec_rhs, petsc_vec_sol));

    // free petsc memory
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&petsc_mat_A));
    PetscCall(VecDestroy(&petsc_vec_rhs));
    PetscCall(VecDestroy(&petsc_vec_sol));

    PetscCall(PetscFinalize());
    return 0;
}

/*
 * command option
 * -ksp_type fgmres -ksp_gmres_restart 50 -ksp_monitor_true_residual -ksp_max_it 2000
 * -ksp_rtol 1e-8 -pc_type none -user_defined_pc
 */
