function panelinPE(A,B,corner0,corner1)
implicit none
double precision, dimension(1:2), intent(in):: A,B,corner0,corner1
logical panelinPE
integer:: i
integer, dimension(1:2):: m,n,l
double precision, dimension(1:2):: A2,B2,C2,P2
double precision:: coord,y,x

panelinPE=.false.

! test if vertex A,B in box [corner0,corner1]
if((corner0(1)<=A(1).and.corner1(1)>=A(1).and.&
    corner0(2)<=A(2).and.corner1(2)>=A(2)).or.&
   (corner0(1)<=B(1).and.corner1(1)>=B(1).and.&
    corner0(2)<=B(2).and.corner1(2)>=B(2))) then
  panelinPE=.true.
  return
endif

!test if vertex point on edge of box [corner0,corner1] 
if(A(1).eq.B(1)) then
  if((A(1).eq.corner0(1)).or.(A(1).eq.corner1(1))) then
        panelinPE=.true.
        return  
  endif
endif
if(A(2).eq.B(2)) then
  if((A(2).eq.corner0(2)).or.(A(2).eq.corner1(2))) then
        panelinPE=.true.
        return
  endif
endif

!test if panel cross processor but vertex not on it
if((A(1).ne.B(1)).and.(A(2).ne.B(2))) then 
  y=(B(2)-A(2))/(B(1)-A(1))*(corner0(1)-A(1))+A(2)
  if((y.lt.corner1(2)).and.(y.gt.corner0(2))) then
    panelinPE=.true.
    return
  endif
  y=(B(2)-A(2))/(B(1)-A(1))*(corner1(1)-A(1))+A(2)
  if((y.gt.corner0(2)).and.(y.lt.corner1(2))) then
    panelinPE=.true.
    return
  endif
  x=(corner0(2)-A(2))*(B(1)-A(1))/(B(2)-A(2))+A(1)
  if((x.gt.corner0(1)).and.(x.lt.corner1(1))) then
    panelinPE=.true.
    return
  endif
  x=(corner1(2)-A(2))*(B(1)-A(1))/(B(2)-A(2))+A(1)
  if((x.gt.corner0(1)).and.(x.lt.corner1(1))) then
    panelinPE=.true.
    return
  endif
endif

end function panelinPE
