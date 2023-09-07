SUBROUTINE cfl

  USE MPI
  USE para
  USE Euler
  DOUBLE PRECISION:: umax,vmax,dtc,dtv,ug,vg

  IF(icfl==0) THEN
     dt=dt0
     t=t+dt
     RETURN
  ENDIF

  umax=0.0d0
  DO j=jpbeg-1,jpend+1
     DO i=ipbeg-1,ipend
        umax=MAX(umax,ABS(u(i,j)))
     ENDDO
  ENDDO

  vmax=0.0d0
  DO j=jpbeg-1,jpend
     DO i=ipbeg-1,ipend+1
        vmax=MAX(vmax,ABS(v(i,j)))
     ENDDO
  ENDDO

  !write(*,*),'umax',umax,'vmax',vmax

  dx1=1.0d0/dxi
  dy1=1.0d0/deta
  dx2=dxi*dxi
  dy2=deta*deta
  !call MPI_Allreduce(umax,ug,1,MPI_DOUBLE,MPI_MAX,comm2d,ierr)
  !call MPI_Allreduce(vmax,vg,1,MPI_DOUBLE,MPI_MAX,comm2d,ierr)
  CALL MPI_Reduce(umax,ug,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  CALL MPI_Reduce(vmax,vg,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm2d,ierr)
  IF(myid==0) THEN
     IF(ug==0.0.AND.vg==0.0) THEN
        dtc=100.0d0
     ELSE
        dtc=cflc/(ug*dx1+vg*dy1)
     ENDIF
     dtv=cflv*Re/(dx2+dy2)

     dt=MIN(dtc,dtv,dtcfl)
     dt1=1.0d0/dt
     dt2=dt1/dt

  ENDIF

  CALL MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,comm2d,ierr)
  t=t+dt

END SUBROUTINE cfl
