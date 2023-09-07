subroutine moveobject
use para
use Lagrange
integer:: myobj,n0,n1

do myobj=1,nobj4proc
  n0=nvertex4proc(myobj-1)+1
  n1=nvertex4proc(myobj)
  vertex(1,n0:n1)=vertex(1,n0:n1)-1.0d0
  vertex(2,n0:n1)=vertex(2,n0:n1)-1.0d0
!  vertex(3,n0:n1)=vertex(3,n0:n1)-1.0d0
end do

end subroutine moveobject
