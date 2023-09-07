subroutine indexrange(n,np,id,ibeg,iend)
integer, intent(in):: n,np,id
integer, intent(out):: ibeg,iend
integer:: nloc,nr

nloc = floor(dble(n/np))
nr = mod(n,np)

ibeg = id*nloc+1+min(id,nr)
if(id < nr) then
  nloc = nloc+1
endif
iend = ibeg+nloc-1
if((id.eq.np-1).or.(iend > n)) then
  iend = n
endif

end subroutine indexrange
