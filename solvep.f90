
SUBROUTINE solvep
  USE MPI
  USE para
  USE Euler
  USE stretchxy

  INTEGER:: iters,sendid
  INTEGER:: status(MPI_STATUS_SIZE)
  INTEGER, DIMENSION(1:2):: indices
  DOUBLE PRECISION:: resid,sum1,sum_p,aver_p,summ,summ1,diff

  DOUBLE PRECISION:: maxerr,pnorm
  INTEGER :: maxxc,maxyc

  ! Step 5: setup the solver
  CALL HYPRE_StructSMGCreate(MPI_COMM_WORLD,solver,ierr)
  CALL HYPRE_StructSMGSetMemoryUse(solver, 0, ierr);
  CALL HYPRE_StructSMGSetMaxIter(solver, maxit,ierr)
  CALL HYPRE_StructSMGSetTol(solver, tol,ierr)
  CALL HYPRE_StructSMGSetRelChange(solver, relch,ierr)
  !  call HYPRE_StructSMGSetNonZeroGuess(solver,ierr)
  CALL HYPRE_StructSMGSetNumPreRelax(solver, npre,ierr)
  CALL HYPRE_StructSMGSetNumPostRelax(solver, npost,ierr)
  CALL HYPRE_StructSMGSetPrintLevel(solver, 0,ierr)
  CALL HYPRE_StructSMGSetLogging(solver, 1,ierr)

  ! clean rhs of p to make solution exist
  CALL rhs4p_clean

  ! assign values to bvector and xvector
  DO j = jpbeg,jpend
     indices(2) = j
     DO i = ipbeg,ipend
        indices(1) = i
        CALL HYPRE_StructVectorSetValues(bvector,indices,rhsp(i,j),ierr)
        CALL HYPRE_StructVectorSetValues(xvector,indices,p(i,j),ierr)
     END DO
  END DO
  CALL HYPRE_StructVectorAssemble(bvector,ierr)
  CALL HYPRE_StructVectorAssemble(xvector,ierr)

  ! solve the system: Lmatrix*xvector=bvector
  CALL HYPRE_StructSMGSetup(solver,Lmatrix,bvector,xvector,ierr)
  CALL HYPRE_StructSMGSolve(solver,Lmatrix,bvector,xvector,ierr)

  resid=0
  iters=0
  ! return the norm of the final relative residual and the number of iterations
  CALL HYPRE_StructSMGGetFinalRelative(solver,resid,ierr)
  CALL HYPRE_StructSMGGetNumIterations(solver,iters,ierr)
  IF(myid==0) THEN
     PRINT*
     PRINT*, "Iterations = ", iters
     PRINT*, "Final Relative Residual Norm = ", resid
     PRINT*
  ENDIF

  ! arrange the solution in the 2d array
  DO j = jpbeg,jpend
     indices(2) = j
     DO i = ipbeg,ipend
        indices(1) = i
        CALL HYPRE_StructVectorGetValues(xvector,indices,p(i,j),ierr)
     END DO
  END DO

  ! exchange pressure
  ! i-direction
  CALL MPI_SENDRECV(p(ipend-nxghost+1, jpbeg-nyghost),1,pitype,right,0,&
       p(ipbeg-nxghost, jpbeg-nyghost),1,pitype,left, 0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(p(ipbeg,jpbeg-nyghost),1,pitype,left, 0,&
       p(ipend+1,  jpbeg-nyghost),1,pitype,right,0,&
       comm2d,status,ierr)
  ! j-direction
  CALL MPI_SENDRECV(p(ipbeg-nxghost, jpend-nyghost+1),1,pjtype,back, 0,&
       p(ipbeg-nxghost, jpbeg-nyghost),1,pjtype,front,0,&
       comm2d,status,ierr)
  CALL MPI_SENDRECV(p(ipbeg-nxghost, jpbeg),1,pjtype,front,0,&
       p(ipbeg-nxghost, jpend+1),1,pjtype,back, 0,&
       comm2d,status,ierr)

  ! clean up pressure
  sum1=zero
  sum_p=zero
  DO i=ipbeg,ipend
     DO j=jpbeg,jpend
        sum1=sum1+p(i,j)
     ENDDO
  ENDDO
  CALL MPI_Allreduce(sum1,sum_p,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm2d,ierr)
  aver_p=sum_p/nx/ny
  PRINT*, 'aver_p',aver_p
  p=p-aver_p

  diff = p(1,1)-SIN(xc(1))*SIN(yc(1))
  CALL MPI_Bcast(diff,1,MPI_DOUBLE_PRECISION,0,comm2d,ierr)
  p=p-diff

  maxerr = -100.0d0
  maxxc = -10000
  maxyc = -10000
  DO j=jpbeg,jpend
     DO i=ipbeg,ipend
        IF(iop(i,j).NE.0) THEN
           IF(maxerr.LE.ABS(p(i,j)-EXP(-xc(i))*EXP(-yc(j)))) THEN
              maxxc = i
              maxyc = j
           ENDIF
           maxerr = MAX(maxerr,ABS(p(i,j)-EXP(-xc(i))*EXP(-yc(j))))
        ELSE
           IF(maxerr.LE.ABS(p(i,j)-SIN(xc(i))*SIN(yc(j)))) THEN
              maxxc = i
              maxyc = j
           ENDIF
           maxerr = MAX(maxerr,ABS(p(i,j)-SIN(xc(i))*SIN(yc(j))))
        ENDIF
     ENDDO
  ENDDO
  CALL MPI_Allreduce(maxerr,pnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)

  PRINT*,'pnorm = ', pnorm
  ! print*, maxxc, maxyc,xc(maxxc),yc(maxyc),sqrt(xc(maxxc)*xc(maxxc)+yc(maxyc)*yc(maxyc)),rhsp(maxxc,maxyc)
  CALL HYPRE_StructSMGDestroy(solver,ierr)

END SUBROUTINE solvep
