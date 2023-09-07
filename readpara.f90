SUBROUTINE readpara

  USE para
  USE Lagrange
  ! read all parameters
  NAMELIST /domainpara/xi0,xi1,eta0,eta1,nx,ny,px,py,nstep
  NAMELIST /hyprepara/tol,maxit,relch,rlxtype,npre,npost
  NAMELIST /flowpara/Re,icfl,cflv,cflc,dtcfl,dt0,tout, &
       iread,iwrite,itout,iplot,ianimation
  NAMELIST /BCs/ipbcx0,ipbcx1,jpbcy0,jpbcy1, &
       iubcx0,iubcx1,jubcy0,jubcy1, &
       ivbcx0,ivbcx1,jvbcy0,jvbcy1
  NAMELIST /objectinfo/nobj

  OPEN(100, file='INPUT/domainpara.txt')
  READ(100, domainpara)
  CLOSE(100)

  OPEN(200, file='INPUT/hyprepara.txt')
  READ(200, hyprepara)
  CLOSE(200)

  OPEN(300, file='INPUT/flowpara.txt')
  READ(300, flowpara)
  CLOSE(300)

  OPEN(400, file='INPUT/BCs.txt')
  READ(400, BCs)
  CLOSE(400)

  OPEN(500,file='INPUT/OBJECT/objlist')
  READ(500,objectinfo)
  CLOSE(500)

  ! Some setups
  ! =================================================================== !

  ! 1/Re
  Re1=1.0d0/Re

  ! number of layers of ghost points
  nxghost=4
  nyghost=4

  ! average number of grid points separating two consecutive intersections along a grid line
  !????????????????????
  nseparate=8

  IF(nobj>0) isingular = .TRUE.

  ! ================================================================== !

END SUBROUTINE readpara
