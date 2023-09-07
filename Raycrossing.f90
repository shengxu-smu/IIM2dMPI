SUBROUTINE Raycrossing

  USE MPI
  USE para
  USE Euler
  USE Lagrange
  USE stretchxy

  INTEGER::myobj,myvertex,ixc,ixf,jyc,jyf,jint,iint,counter
  DOUBLE PRECISION, DIMENSION(1:2)::A,B,corner0,corner1
  CHARACTER(5) filesuffix
  DOUBLE PRECISION:: signx,signy,yint,xint,su,sv

  INTERFACE
     FUNCTION panelinPE(A,B,corner0,corner1)
       DOUBLE PRECISION, DIMENSION(1:2), INTENT(in):: A,B,corner0,corner1
       LOGICAL panelinPE
     END FUNCTION panelinPE
  END INTERFACE

  n_xc_int = 0
  n_yc_int = 0
  n_xf_int = 0
  n_yf_int = 0

  su = 0
  sv = 0
  iop=0
  iou=0
  iov=0

  corner0=allcorner0(:,myid)
  corner1=allcorner1(:,myid)
  ! ============== Raycrossing part =====================
  ! loop over xc
  DO i=ipbeg-1,ipend+1
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           ! panel
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(1)>B(1)) THEN
              IF((A(1)>xc(i)).AND.(B(1)<=xc(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xc(i)-A(1))/(B(1)-A(1))
                 jint=INT((y2eta(yint)-eta0)/deta*2)
                 IF(MOD(jint,2)==0) THEN
                    jyc=jint/2
                    jyf=jint/2
                 ELSE
                    jyc=(jint+1)/2
                    jyf=(jint-1)/2
                 ENDIF
                 jyc=MIN(jyc,jpend+1)
                 jyf=MIN(jyf,jvend)
                 IF(yint==yc(jyc)) jyc=jyc-MAX(0,INT(signy))
                 IF(yint==yf(jyf)) jyf=jyf-MAX(0,INT(signy))
                 iop(i, jpbeg-1:jyc)=myobj-iop(i, jpbeg-1:jyc)
                 iov(i, jvbeg:jyf)=myobj-iov(i, jvbeg:jyf)
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    n_xc_int = n_xc_int+1
                 ENDIF
              ENDIF
           ELSEIF(A(1)<B(1)) THEN
              IF((A(1)<=xc(i)).AND.(B(1)>xc(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xc(i)-A(1))/(B(1)-A(1))
                 jint=INT((y2eta(yint)-eta0)/deta*2)
                 IF(MOD(jint,2)==0) THEN
                    jyc=jint/2
                    jyf=jint/2
                 ELSE
                    jyc=(jint+1)/2
                    jyf=(jint-1)/2
                 ENDIF
                 jyc=MIN(jyc,jpend+1)
                 jyf=MIN(jyf,jvend)
                 IF(yint==yc(jyc)) jyc=jyc-MAX(0,INT(signy))
                 IF(yint==yf(jyf)) jyf=jyf-MAX(0,INT(signy))
                 iop(i, jpbeg-1:jyc)=myobj-iop(i, jpbeg-1:jyc)
                 iov(i, jvbeg:jyf)=myobj-iov(i, jvbeg:jyf)
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    n_xc_int = n_xc_int+1
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over xf
  DO i=iubeg,iuend
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(1)>B(1)) THEN
              IF((A(1)>xf(i)).AND.(B(1)<=xf(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))

                 yint=A(2)+(B(2)-A(2))*(xf(i)-A(1))/(B(1)-A(1))
                 jint=INT((y2eta(yint)-eta0)/deta*2)
                 IF(MOD(jint,2)==0) THEN
                    jyc=jint/2
                 ELSE
                    jyc=(jint+1)/2
                 ENDIF
                 jyc=MIN(jyc,juend)
                 IF(yint==yc(jyc)) jyc=jyc-MAX(0,INT(signy))
                 iou(i, jubeg:jyc)=myobj-iou(i, jubeg:jyc)
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    n_xf_int = n_xf_int+1
                 ENDIF
              ENDIF
           ELSEIF(A(1)<B(1)) THEN
              IF((A(1)<=xf(i)).AND.(B(1)>xf(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xf(i)-A(1))/(B(1)-A(1))
                 jint=INT((y2eta(yint)-eta0)/deta*2)
                 IF(MOD(jint,2)==0) THEN
                    jyc=jint/2
                 ELSE
                    jyc=(jint+1)/2
                 ENDIF
                 jyc=MIN(jyc,juend)
                 IF(yint==yc(jyc)) jyc=jyc-MAX(0,INT(signy))
                 iou(i, jubeg:jyc)=myobj-iou(i, jubeg:jyc)
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    n_xf_int = n_xf_int+1
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over yc
  DO j=jpbeg-1,jpend+1
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(2)>B(2)) THEN
              IF((A(2)>yc(j)).AND.(B(2)<=yc(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yc(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xc(ipbeg) .AND. xint<=xc(ipend)) THEN
                 !  IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                 IF((xint>corner0(1)).AND.(xint<corner1(1))) THEN
                    n_yc_int = n_yc_int+1
                 ENDIF
              ENDIF
           ELSEIF(A(2)<B(2)) THEN
              IF((A(2)<=yc(j)).AND.(B(2)>yc(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yc(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xc(ipbeg) .AND. xint<=xc(ipend)) THEN
                 !  IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                 IF((xint>corner0(1)).AND.(xint<corner1(1))) THEN
                    n_yc_int = n_yc_int+1
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over yf
  DO j=jvbeg, jvend
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(2)>B(2)) THEN
              IF((A(2)>yf(j)).AND.(B(2)<=yf(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yf(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xf(ivbeg) .AND. xint<=xf(ivend)) THEN
                 IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                    n_yf_int = n_yf_int+1
                 ENDIF
              ENDIF
           ELSEIF(A(2)<B(2)) THEN
              IF((A(2)<=yf(j)).AND.(B(2)>yf(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yf(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xf(ivbeg) .AND. xint<=xf(ivend)) THEN
                 IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                    n_yf_int = n_yf_int+1
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! ============= update intersection pts information ==================
  ! xc
  ALLOCATE(xc_int_info(1:14,1:n_xc_int))
  ALLOCATE(xc_int_vertex(1:5,1:n_xc_int))
  ! yc
  ALLOCATE(yc_int_info(1:14,1:n_yc_int))
  ALLOCATE(yc_int_vertex(1:5,1:n_yc_int))
  ! xf
  ALLOCATE(xf_int_info(1:10,1:n_xf_int))
  ALLOCATE(xf_int_vertex(1:5,1:n_xf_int))
  ! yf
  ALLOCATE(yf_int_info(1:10,1:n_yf_int))
  ALLOCATE(yf_int_vertex(1:5,1:n_yf_int))
  ! loop over xc
  counter = 0
  DO i=ipbeg-1,ipend+1
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           ! panel
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(1)>B(1)) THEN
              IF((A(1)>xc(i)).AND.(B(1)<=xc(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xc(i)-A(1))/(B(1)-A(1))
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    jint=INT((y2eta(yint)-eta0)/deta*2)
                    IF(MOD(jint,2)==0) THEN
                       jyc=jint/2
                       jyf=jint/2
                    ELSE
                       jyc=(jint+1)/2
                       jyf=(jint-1)/2
                    ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(xc(i)-A(1))/(B(1)-A(1))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(xc(i)-A(1))/(B(1)-A(1))
                    counter = counter + 1
                    xc_int_vertex(1,counter) = myvertex
                    xc_int_vertex(2,counter) = i
                    xc_int_vertex(3,counter) = jyc
                    xc_int_vertex(4,counter) = jyf
                    ! print*,counter,i,jyc,jyf
                    xc_int_vertex(5,counter) = myobj

                    ! print*, myvertex, i, jyc, jyf

                    xc_int_info(1, counter) = xc(i)
                    xc_int_info(2, counter) = yint
                    xc_int_info(3, counter) = signx
                    xc_int_info(4, counter) = signy
                    xc_int_info(5, counter) = su
                    xc_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ELSEIF(A(1)<B(1)) THEN
              IF((A(1)<=xc(i)).AND.(B(1)>xc(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xc(i)-A(1))/(B(1)-A(1))
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    jint=INT((y2eta(yint)-eta0)/deta*2)
                    IF(MOD(jint,2)==0) THEN
                       jyc=jint/2
                       jyf=jint/2
                    ELSE
                       jyc=(jint+1)/2
                       jyf=(jint-1)/2
                    ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(xc(i)-A(1))/(B(1)-A(1))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(xc(i)-A(1))/(B(1)-A(1))
                    counter = counter + 1
                    xc_int_vertex(1,counter) = myvertex
                    xc_int_vertex(2,counter) = i
                    xc_int_vertex(3,counter) = jyc
                    xc_int_vertex(4,counter) = jyf
                    ! print*,counter,i,jyc,jyf
                    xc_int_vertex(5,counter) = myobj

                    ! print*, myvertex, i, jyc, jyf

                    xc_int_info(1, counter) = xc(i)
                    xc_int_info(2, counter) = yint
                    xc_int_info(3, counter) = signx
                    xc_int_info(4, counter) = signy
                    xc_int_info(5, counter) = su
                    xc_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over xf
  counter = 0
  DO i=iubeg,iuend
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(1)>B(1)) THEN
              IF((A(1)>xf(i)).AND.(B(1)<=xf(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xf(i)-A(1))/(B(1)-A(1))
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    jint=INT((y2eta(yint)-eta0)/deta*2)
                    IF(MOD(jint,2)==0) THEN
                       jyc=jint/2
                       jyf=jint/2
                    ELSE
                       jyc=(jint+1)/2
                       jyf=(jint-1)/2
                    ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(xf(i)-A(1))/(B(1)-A(1))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(xf(i)-A(1))/(B(1)-A(1))
                    counter = counter + 1
                    xf_int_vertex(1, counter) = myvertex
                    xf_int_vertex(2, counter) = i
                    xf_int_vertex(3, counter) = jyc
                    xf_int_vertex(4, counter) = jyf
                    xf_int_vertex(5, counter) = myobj


                    xf_int_info(1, counter) = xf(i)
                    xf_int_info(2, counter) = yint
                    xf_int_info(3, counter) = signx
                    xf_int_info(4, counter) = signy
                    xf_int_info(5, counter) = su
                    xf_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ELSEIF(A(1)<B(1)) THEN
              IF((A(1)<=xf(i)).AND.(B(1)>xf(i))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 yint=A(2)+(B(2)-A(2))*(xf(i)-A(1))/(B(1)-A(1))
                 IF((yint>corner0(2)).AND.(yint<corner1(2))) THEN
                    jint=INT((y2eta(yint)-eta0)/deta*2)
                    IF(MOD(jint,2)==0) THEN
                       jyc=jint/2
                       jyf=jint/2
                    ELSE
                       jyc=(jint+1)/2
                       jyf=(jint-1)/2
                    ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(xf(i)-A(1))/(B(1)-A(1))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(xf(i)-A(1))/(B(1)-A(1))
                    counter = counter + 1
                    xf_int_vertex(1, counter) = myvertex
                    xf_int_vertex(2, counter) = i
                    xf_int_vertex(3, counter) = jyc
                    xf_int_vertex(4, counter) = jyf
                    xf_int_vertex(5, counter) = myobj

                    xf_int_info(1, counter) = xf(i)
                    xf_int_info(2, counter) = yint
                    xf_int_info(3, counter) = signx
                    xf_int_info(4, counter) = signy
                    xf_int_info(5, counter) = su
                    xf_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over yc
  counter = 0
  DO j=jpbeg-1,jpend+1
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(2)>B(2)) THEN
              IF((A(2)>yc(j)).AND.(B(2)<=yc(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yc(j)-A(2))/(B(2)-A(2))
                 IF((xint>corner0(1)).AND.(xint<corner1(1))) THEN
                    iint = INT((x2xi(xint)-xi0)/dxi*2)
                    IF(MOD(iint,2)==0) THEN
                       ixc=iint/2
                       ixf=iint/2
                    ELSE
                       ixc=(iint+1)/2
                       ixf=(iint-1)/2
                    ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(yc(j)-A(2))/(B(2)-A(2))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(yc(j)-A(2))/(B(2)-A(2))
                    counter = counter + 1
                    yc_int_vertex(1, counter) = myvertex
                    yc_int_vertex(2, counter) = j
                    yc_int_vertex(3, counter) = ixc
                    yc_int_vertex(4, counter) = ixf
                    yc_int_vertex(5, counter) = myobj

                    yc_int_info(1, counter) = yc(j)
                    yc_int_info(2, counter) = xint
                    yc_int_info(3, counter) = signx
                    yc_int_info(4, counter) = signy
                    yc_int_info(5, counter) = su
                    yc_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ELSEIF(A(2)<B(2)) THEN
              IF((A(2)<=yc(j)).AND.(B(2)>yc(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yc(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xc(ipbeg) .AND. xint<=xc(ipend)) THEN
                 !  IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                 IF((xint>corner0(1)).AND.(xint<corner1(1))) THEN
                   iint = INT((x2xi(xint)-xi0)/dxi*2)
                   IF(MOD(iint,2)==0) THEN
                      ixc=iint/2
                      ixf=iint/2
                   ELSE
                      ixc=(iint+1)/2
                      ixf=(iint-1)/2
                   ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(yc(j)-A(2))/(B(2)-A(2))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(yc(j)-A(2))/(B(2)-A(2))
                    counter = counter + 1
                    yc_int_vertex(1, counter) = myvertex
                    yc_int_vertex(2, counter) = j
                    yc_int_vertex(3, counter) = ixc
                    yc_int_vertex(4, counter) = ixf
                    yc_int_vertex(5, counter) = myobj

                    yc_int_info(1, counter) = yc(j)
                    yc_int_info(2, counter) = xint
                    yc_int_info(3, counter) = signx
                    yc_int_info(4, counter) = signy
                    yc_int_info(5, counter) = su
                    yc_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! loop over yf
  counter = 0
  DO j=jvbeg, jvend
     DO myobj=1,nobj4proc
        DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
           A(1)=vertex(1,myvertex)
           A(2)=vertex(2,myvertex)
           B(1)=vertex(1,myvertex+1)
           B(2)=vertex(2,myvertex+1)
           IF(A(2)>B(2)) THEN
              IF((A(2)>yf(j)).AND.(B(2)<=yf(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yf(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xf(ivbeg) .AND. xint<=xf(ivend)) THEN
                 IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                   iint = INT((x2xi(xint)-xi0)/dxi*2)
                   IF(MOD(iint,2)==0) THEN
                      ixc=iint/2
                      ixf=iint/2
                   ELSE
                      ixc=(iint+1)/2
                      ixf=(iint-1)/2
                   ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(yf(j)-A(2))/(B(2)-A(2))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(yf(j)-A(2))/(B(2)-A(2))
                    counter = counter + 1
                    yf_int_vertex(1, counter) = myvertex
                    yf_int_vertex(2, counter) = j
                    yf_int_vertex(3, counter) = ixc
                    yf_int_vertex(4, counter) = ixf
                    yf_int_vertex(5, counter) = myobj

                    yf_int_info(1, counter) = yf(j)
                    yf_int_info(2, counter) = xint
                    yf_int_info(3, counter) = signx
                    yf_int_info(4, counter) = signy
                    yf_int_info(5, counter) = su
                    yf_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ELSEIF(A(2)<B(2)) THEN
              IF((A(2)<=yf(j)).AND.(B(2)>yf(j))) THEN
                 signx=SIGN(1.0d0,B(2)-A(2))
                 signy=SIGN(1.0d0,A(1)-B(1))
                 xint=A(1)+(B(1)-A(1))*(yf(j)-A(2))/(B(2)-A(2))
                 !  IF(xint>=xf(ivbeg) .AND. xint<=xf(ivend)) THEN
                 IF(xint>=xc(ipbeg-1) .AND. xint<=xc(ipend+1)) THEN
                   iint = INT((x2xi(xint)-xi0)/dxi*2)
                   IF(MOD(iint,2)==0) THEN
                      ixc=iint/2
                      ixf=iint/2
                   ELSE
                      ixc=(iint+1)/2
                      ixf=(iint-1)/2
                   ENDIF
                    !  su=us(myvertex)+(us(myvertex+1)-us(myvertex))*(yf(j)-A(2))/(B(2)-A(2))
                    !  sv=vs(myvertex)+(vs(myvertex+1)-vs(myvertex))*(yf(j)-A(2))/(B(2)-A(2))
                    counter = counter + 1
                    yf_int_vertex(1, counter) = myvertex
                    yf_int_vertex(2, counter) = j
                    yf_int_vertex(3, counter) = ixc
                    yf_int_vertex(4, counter) = ixf
                    yf_int_vertex(5, counter) = myobj

                    yf_int_info(1, counter) = yf(j)
                    yf_int_info(2, counter) = xint
                    yf_int_info(3, counter) = signx
                    yf_int_info(4, counter) = signy
                    yf_int_info(5, counter) = su
                    yf_int_info(6, counter) = sv
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  ! if(myid==0) then
  ! DO myobj=1,nobj
  !    object=object_list(myobj)
  !   !  write(filesuffix,'(i0)') object
  !   !  open(100,file='INPUT/OBJECT/vobj'//filesuffix)
  !
  !    open(unit=69,file='OUTPUT/object1.dat',status='unknown')
  !    do myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
  !      write(69,300)(vertex(i,myvertex),i=1,2)
  !    enddo
  !    close(69)
  ! ENDDO
  ! endif
  ! 300 FORMAT(1x,30e16.6)

END SUBROUTINE Raycrossing
