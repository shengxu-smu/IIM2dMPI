SUBROUTINE Cartesiangrid
  USE MPI
  USE para
  USE Euler
  USE stretchxy

  DOUBLE PRECISION:: ksi,ita
  DOUBLE PRECISION:: xc2,xc1,xc0,xf1,xf0,yc2,yc1,yc0,yf1,yf0
  DOUBLE PRECISION, DIMENSION(1:2):: corner0,corner1
  INTEGER, DIMENSION(1:15):: gatherbuf

  ALLOCATE(xic(ipbeg-5:ipend+4))
  ALLOCATE(xicx(ipbeg-5:ipend+4))
  ALLOCATE(xicxx(ipbeg-5:ipend+4))
  ALLOCATE(xc(ipbeg-5:ipend+4))

  ALLOCATE(xif(iubeg-5:iuend+4))
  ALLOCATE(xifx(iubeg-5:iuend+4))
  ALLOCATE(xifxx(iubeg-5:iuend+4))
  ALLOCATE(xf(iubeg-5:iuend+4))

  ALLOCATE(etac(jpbeg-5:jpend+4))
  ALLOCATE(etacy(jpbeg-5:jpend+4))
  ALLOCATE(etacyy(jpbeg-5:jpend+4))
  ALLOCATE(yc(jpbeg-5:jpend+4))

  ALLOCATE(etaf(jvbeg-5:jvend+4))
  ALLOCATE(etafy(jvbeg-5:jvend+4))
  ALLOCATE(etafyy(jvbeg-5:jvend+4))
  ALLOCATE(yf(jvbeg-5:jvend+4))

  ! x=x(xi) & xi=xi(x)
  ! xi0=x2xi(x0)
  ! xi1=x2xi(x1)

  IF(ipbcx0==2.AND.ipbcx1==2) THEN
     dxi=(xi1-xi0)/DBLE(nx)
  END IF
  IF(ipbcx0==1.AND.ipbcx1==1) THEN
     dxi=(xi1-xi0)/DBLE(nx+1)
  END IF
  IF((ipbcx0==1.AND.ipbcx1==2).OR.(ipbcx0==2.AND.ipbcx1==1).OR. &
       (ipbcx0==0.AND.ipbcx1==0)) THEN
     dxi=(xi1-xi0)/DBLE(nx)
  END IF
  halfdxi=HALF*dxi
  twodxi=TWO*dxi

  DO i=ipbeg-5,ipend+4
     IF(ipbcx0==2.AND.ipbcx1==2) THEN
        ksi=xi0+(DBLE(i)-HALF)*dxi
     END IF
     IF(ipbcx0==1.AND.ipbcx1==1) THEN
        ksi=xi0+DBLE(i)*dxi
     END IF
     xic(i)=ksi
     xc(i)=xi2x(ksi)
     xicx(i)=dxidx(xc(i))
     xicxx(i)=d2xidx2(xc(i))
     !  print*,i,ksi,xc(i),xicx(i),xicxx(i)
  ENDDO
  DO i=iubeg-5,iuend+4
     IF(ipbcx0==2.AND.ipbcx1==2) THEN
        ksi=xi0+(DBLE(i))*dxi
     END IF
     IF(ipbcx0==1.AND.ipbcx1==1) THEN
        ksi=xi0+(DBLE(i)+HALF)*dxi
     END IF
     xif(i)=ksi
     xf(i)=xi2x(ksi)
     xifx(i)=dxidx(xf(i))
     xifxx(i)=d2xidx2(xf(i))
  ENDDO

  ! y=y(eta) & eta=eta(y)
  ! eta0=y2eta(y0)
  ! eta1=y2eta(y1)
  IF(jpbcy0==2.AND.jpbcy1==2) THEN
     deta=(eta1-eta0)/DBLE(ny)
  END IF
  IF(jpbcy0==1.AND.jpbcy1==1) THEN
     deta=(eta1-eta0)/DBLE(ny+1)
  END IF
  IF((jpbcy0==1.AND.jpbcy1==2).OR.(jpbcy0==2.AND.jpbcy1==1).OR. &
       (jpbcy0==0.AND.jpbcy1==0)) THEN
     deta=(eta1-eta0)/DBLE(ny)
  END IF
  halfdeta=HALF*deta
  twodeta=TWO*deta

  DO j=jpbeg-5,jpend+4
     IF(jpbcy0==2.AND.jpbcy1==2) THEN
        ita=eta0+(DBLE(j)-HALF)*deta
     END IF
     IF(jpbcy0==1.AND.jpbcy1==1) THEN
        ita=eta0+DBLE(j)*deta
     END IF
     etac(j)=ita
     yc(j)=eta2y(ita)
     etacy(j)=detady(yc(j))
     etacyy(j)=d2etady2(yc(j))
  ENDDO
  DO j=jvbeg-5,jvend+4
     IF(jpbcy0==2.AND.jpbcy1==2) THEN
        ita=eta0+(DBLE(j))*deta
     END IF
     IF(jpbcy0==1.AND.jpbcy1==1) THEN
        ita=eta0+(DBLE(j)+HALF)*deta
     END IF
     etaf(j)=ita
     yf(j)=eta2y(ita)
     etafy(j)=detady(yf(j))
     etafyy(j)=d2etady2(yf(j))
  ENDDO

  ! ratios of spatial steps
  r2dxideta=(dxi*dxi)/(deta*deta)

  ! gather corners from all processors and then broadcast all corners to all processors
  corner0=(/xf(iubeg),yf(jvbeg)/)
  corner1=(/xf(iuend),yf(jvend)/)

  ALLOCATE(allcorner0(1:2,0:nprocs-1))
  ALLOCATE(allcorner1(1:2,0:nprocs-1))
  ALLOCATE(ProcsInfo(1:15,0:nprocs-1))
  gatherbuf=(/myid,0,0,ipbeg,ipend,jpbeg,jpend,iubeg,iuend,jubeg,juend,ivbeg,ivend,jvbeg,jvend/)
  CALL MPI_ALLGATHER(corner0,2,MPI_DOUBLE_PRECISION,allcorner0,2,MPI_DOUBLE_PRECISION,comm2d,ierr)
  CALL MPI_ALLGATHER(corner1,2,MPI_DOUBLE_PRECISION,allcorner1,2,MPI_DOUBLE_PRECISION,comm2d,ierr)
  CALL MPI_ALLGATHER(gatherbuf,15,MPI_INTEGER,ProcsInfo,15,MPI_INTEGER,comm2d,ierr)

  ! set up z vector for xATz = bTz that ATz = 0
  ALLOCATE(ppp(ipbeg:ipend))
  ALLOCATE(qqq(jpbeg:jpend))

  xc2=xi2x(xi0+3.0d0/2.0d0*dxi)
  xc1=xi2x(xi0+dxi/2.0d0)
  xc0=xi2x(xi0-dxi/2.0d0)
  xf1=xi2x(xi0+dxi)
  xf0=xi2x(xi0)

  yc2=eta2y(eta0+3.0d0/2.0d0*deta)
  yc1=eta2y(eta0+deta/2.0d0)
  yc0=eta2y(eta0-deta/2.0d0)
  yf1=eta2y(eta0+deta)
  yf0=eta2y(eta0)

  DO i=ipbeg,ipend
     IF(i==1) THEN
        ppp(i)=1.0d0
     ELSEIF(i==nx) THEN
        ppp(nx)=(xf(nx)-xf(nx-1))*(xc(nx+1)-xc(nx))*(xc2-xc0)/(xc(nx+1)-xc(nx-1))/(xf1-xf0)/(xc1-xc0)
     ELSE
        ppp(i)=(xf(i)-xf(i-1))*(xc2-xc0)/(xf1-xf0)/(xc1-xc0)
     ENDIF
  ENDDO

  DO j=jpbeg,jpend
     IF(j==1) THEN
        qqq(j)=1.0d0
     ELSEIF(j==ny) THEN
        qqq(ny)=(yf(ny)-yf(ny-1))*(yc(ny+1)-yc(ny))*(yc2-yc0)/(yc(ny+1)-yc(ny-1))/(yf1-yf0)/(yc1-yc0)
     ELSE
        qqq(j)=(yf(j)-yf(j-1))*(yc2-yc0)/(yf1-yf0)/(yc1-yc0)
     ENDIF
  ENDDO

  ! output xc, yc, xf, yf
  IF(myid==0) THEN

     OPEN(unit=49,file='OUTPUT/xc.dat',status='unknown')
     DO i=1,nx
  	IF(ipbcx0==2.AND.ipbcx1==2) THEN
           ksi=xi0+(DBLE(i)-HALF)*dxi
	END IF
	IF(ipbcx0==1.AND.ipbcx1==1) THEN
           ksi=xi0+DBLE(i)*dxi
  	END IF
        WRITE(49,200)xi2x(ksi)
     ENDDO
     CLOSE(49)

     OPEN(unit=69,file='OUTPUT/xf.dat',status='unknown')
     DO i=1,nx
  	IF(ipbcx0==2.AND.ipbcx1==2) THEN
           ksi=xi0+(DBLE(i)-1)*dxi
        END IF
        IF(ipbcx0==1.AND.ipbcx1==1) THEN
           ksi=xi0+(DBLE(i)-HALF)*dxi
        END IF
        WRITE(69,200)xi2x(ksi)
     ENDDO
     CLOSE(69)

     OPEN(unit=59,file='OUTPUT/yc.dat',status='unknown')
     DO j=1,ny
        IF(jpbcy0==2.AND.jpbcy1==2) THEN
           ita=eta0+(DBLE(j)-HALF)*deta
        END IF
        IF(jpbcy0==1.AND.jpbcy1==1) THEN
           ita=eta0+DBLE(j)*deta
        END IF
        WRITE(59,200)eta2y(ita)
     ENDDO
     CLOSE(59)

     OPEN(unit=79,file='OUTPUT/yf.dat',status='unknown')
     DO j=1,ny
        IF(jpbcy0==2.AND.jpbcy1==2) THEN
           ita=eta0+(DBLE(j)-1)*deta
        END IF
        IF(jpbcy0==1.AND.jpbcy1==1) THEN
           ita=eta0+(DBLE(j)-HALF)*deta
        ENDIF
        WRITE(79,200)eta2y(ita)
     ENDDO
     CLOSE(79)

  ENDIF
200 FORMAT(1x,10e16.6)


END SUBROUTINE Cartesiangrid
