subroutine destroyhypre
use para
use MPI

! destroy hypre handles
call HYPRE_StructSMGDestroy(solver,ierr)
call HYPRE_StructVectorDestroy(bvector,ierr)
call HYPRE_StructVectorDestroy(xvector,ierr)
call HYPRE_StructMatrixDestroy(Lmatrix,ierr)
call HYPRE_StructStencilDestroy(stencil,ierr)
call HYPRE_StructGridDestroy(grid,ierr)

end subroutine destroyhypre
