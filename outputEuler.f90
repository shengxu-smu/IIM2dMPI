subroutine outputEuler
use para
use Euler
character(5) filesuffix

!write(unit=filesuffix,fmt='(i0)') myid ! (i0): using minimum field width appropriate
!open(unit=12,file='OUTPUT/newresult.'//filesuffix,status='unknown')
!write(12,*)'myid = ', myid, 'coords = ',coords
! do j=juend+nyghost,jubeg-nyghost,-1
!    write(12,200) (u(i,j),i=iubeg-nxghost,iuend+nxghost)
!  enddo

call HYPRE_StructVectorDestroy(xvector,ierr)
call HYPRE_StructVectorDestroy(bvector,ierr)
call HYPRE_StructMatrixDestroy(Lmatrix,ierr)
call HYPRE_StructGridDestroy(grid,ierr)
!call HYPRE_StructSMGDestroy(solver,ierr)

200 format(1x, 64f6.2)

deallocate(u)
deallocate(ucc)
deallocate(uce)
deallocate(uee)
deallocate(un)
deallocate(rhsu1)
deallocate(rhsu2)
deallocate(rhsu3)
deallocate(rhsu4)
deallocate(v)
deallocate(vcc)
deallocate(vee)
deallocate(vec)
deallocate(vn)
deallocate(rhsv1)
deallocate(rhsv2)
deallocate(rhsv3)
deallocate(rhsv4)
deallocate(p)
deallocate(rhsp)
deallocate(d)
deallocate(dn)
deallocate(iou)
deallocate(iov)
deallocate(iop)

deallocate(xic)
deallocate(xicx)
deallocate(xicxx)
deallocate(xc)
deallocate(xif)
deallocate(xifx)
deallocate(xifxx)
deallocate(xf)
deallocate(etac)
deallocate(etacy)
deallocate(etacyy)
deallocate(yc)
deallocate(etaf)
deallocate(etafy)
deallocate(etafyy)
deallocate(yf)
deallocate(allcorner0)
deallocate(allcorner1)
deallocate(ProcsInfo)
deallocate(ppp)
deallocate(qqq)

end subroutine outputEuler
