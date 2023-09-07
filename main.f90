PROGRAM main

  USE MPI
  USE para
  USE Euler
  USE Lagrange
  DOUBLE PRECISION:: time0, time1
  INTEGER:: nn

  CALL MPI_INIT(ierr)
  time0 = MPI_WTIME()

  CALL readpara
  CALL startMPI
  CALL Cartesiangrid

  IF(nobj>0) THEN
     CALL locateobject
     CALL loadobject
     CALL exchangeobject
  END IF
  CALL allocateEuler

  CALL setuphypre
  CALL initial

  ! march in time
  nstart=0
  nstart=nstart+1
  nend=nstart+nstep-1
  tstart=t
  DO nn=1, nstep
     CALL cfl
     IF(myid==0) THEN
        PRINT *, 'nn = ',nn, ' t =',t,' dt =',dt
     ENDIF
     CALL rk
  ENDDO
  CALL run_output
  !call data_plot
  CALL outputEuler

  time1 = MPI_WTIME()

  IF(myid.EQ.0) PRINT*, 'Computational time=', time1-time0, 'seconds'
  CALL MPI_FINALIZE(ierr)

END PROGRAM main
