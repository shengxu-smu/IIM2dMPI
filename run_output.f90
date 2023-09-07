subroutine run_output

use mpi
use para
use Euler
use Lagrange

double precision,dimension(ipbeg:ipend,jpbeg:jpend)::rhspsubarray,psubarray,usubarray,vsubarray,dsubarray,dnsubarray
integer:: pfiletype,ufiletype,vfiletype,dfiletype,pfh,ufh,vfh,dfh,info,np,nu,nv,nd
integer:: rhspfiletype,rhspfh
integer,dimension(1:2)::psizes,psubsizes,pstarts
integer,dimension(1:2)::rhspsizes,rhspsubsizes,rhspstarts
integer,dimension(1:2)::usizes,usubsizes,ustarts
integer,dimension(1:2)::vsizes,vsubsizes,vstarts
integer,dimension(1:2)::dsizes,dsubsizes,dstarts
integer(kind=MPI_OFFSET_KIND),parameter::zero_off=0

call ubc
call vbc
call interpolateuv(0)

np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
psizes=(/nx,ny/)
psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
pstarts=(/ipbeg-1,jpbeg-1/)
psubarray=p(ipbeg:ipend,jpbeg:jpend)

nd=(ipend-ipbeg+1)*(jpend-jpbeg+1)
dsizes=(/nx,ny/)
dsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
dstarts=(/ipbeg-1,jpbeg-1/)
dsubarray=d(ipbeg:ipend,jpbeg:jpend)
dnsubarray=dn(ipbeg:ipend,jpbeg:jpend)


np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
rhspsizes=(/nx,ny/)
rhspsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
rhspstarts=(/ipbeg-1,jpbeg-1/)
rhspsubarray=rhsp(ipbeg:ipend,jpbeg:jpend)

nu=(ipend-ipbeg+1)*(jpend-jpbeg+1)
usizes=(/nx,ny/)
usubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
ustarts=(/ipbeg-1,jpbeg-1/)
usubarray=ucc(ipbeg:ipend,jpbeg:jpend)

nv=(ipend-ipbeg+1)*(jpend-jpbeg+1)
vsizes=(/nx,ny/)
vsubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
vstarts=(/ipbeg-1,jpbeg-1/)
vsubarray=vcc(ipbeg:ipend,jpbeg:jpend)


! output pressure
call MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              pfiletype,ierr)
call MPI_Type_commit(pfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/p.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
call MPI_FILE_SET_VIEW(pfh,zero_off,MPI_DOUBLE_PRECISION,pfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(pfh,psubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(pfh,ierr)

! output rhs4p
call MPI_TYPE_CREATE_SUBARRAY(2, rhspsizes, rhspsubsizes, rhspstarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              rhspfiletype,ierr)
call MPI_Type_commit(rhspfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/rhsp.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,rhspfh,ierr)
call MPI_FILE_SET_VIEW(rhspfh,zero_off,MPI_DOUBLE_PRECISION,rhspfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(rhspfh,rhspsubarray,np,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(rhspfh,ierr)

! output divergence
call MPI_TYPE_CREATE_SUBARRAY(2, dsizes, dsubsizes, dstarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              dfiletype,ierr)
call MPI_Type_commit(dfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/d.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,dfh,ierr)
call MPI_FILE_SET_VIEW(dfh,zero_off,MPI_DOUBLE_PRECISION,dfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(dfh,dsubarray,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(dfh,ierr)

! dn
call MPI_TYPE_CREATE_SUBARRAY(2, dsizes, dsubsizes, dstarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              dfiletype,ierr)
call MPI_Type_commit(dfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/dn.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,dfh,ierr)
call MPI_FILE_SET_VIEW(dfh,zero_off,MPI_DOUBLE_PRECISION,dfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(dfh,dnsubarray,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(dfh,ierr)

! output ucc
call MPI_TYPE_CREATE_SUBARRAY(2, usizes, usubsizes, ustarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              ufiletype,ierr)
call MPI_Type_commit(ufiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/uc.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,ufh,ierr)
call MPI_FILE_SET_VIEW(ufh,zero_off,MPI_DOUBLE_PRECISION,ufiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(ufh,usubarray,nu,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(ufh,ierr)


! output vcc
call MPI_TYPE_CREATE_SUBARRAY(2, vsizes, vsubsizes, vstarts,            &
                              MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,&
                              vfiletype,ierr)
call MPI_Type_commit(vfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/vc.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,vfh,ierr)
call MPI_FILE_SET_VIEW(vfh,zero_off,MPI_DOUBLE_PRECISION,vfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(vfh,vsubarray,nv,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(vfh,ierr)


np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
psizes=(/nx,ny/)
psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
pstarts=(/ipbeg-1,jpbeg-1/)

! iop
call MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
                              MPI_ORDER_FORTRAN,MPI_INTEGER2,&
                              pfiletype,ierr)
call MPI_Type_commit(pfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/iop.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
call MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(pfh,iop(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(pfh,ierr)

! iou
np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
psizes=(/nx,ny/)
psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
pstarts=(/ipbeg-1,jpbeg-1/)
call MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
                              MPI_ORDER_FORTRAN,MPI_INTEGER2,           &
                              pfiletype,ierr)
call MPI_Type_commit(pfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/iou.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
call MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(pfh,iou(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(pfh,ierr)

! iov
np=(ipend-ipbeg+1)*(jpend-jpbeg+1)
psizes=(/nx,ny/)
psubsizes=(/ipend-ipbeg+1,jpend-jpbeg+1/)
pstarts=(/ipbeg-1,jpbeg-1/)
call MPI_TYPE_CREATE_SUBARRAY(2, psizes, psubsizes, pstarts,            &
                              MPI_ORDER_FORTRAN,MPI_INTEGER2,           &
                              pfiletype,ierr)
call MPI_Type_commit(pfiletype,ierr)
call MPI_Info_create(info,ierr)
call MPI_Info_set(info,'collective_buffering','true',ierr)
call MPI_File_open(comm2d,'OUTPUT/iov.dat',MPI_MODE_WRONLY + MPI_MODE_CREATE,info,pfh,ierr)
call MPI_FILE_SET_VIEW(pfh,zero_off,MPI_INTEGER2,pfiletype,'native',MPI_INFO_NULL,ierr)
call MPI_FILE_WRITE_ALL(pfh,iov(ipbeg:ipend,jpbeg:jpend),np,MPI_INTEGER2,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(pfh,ierr)

end subroutine run_output
