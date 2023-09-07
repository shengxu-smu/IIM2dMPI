subroutine interpolate(xa,ya,n,x,y,dy)


      integer:: n
      double precision:: dy,x,y
      integer,parameter:: nmax=10
      double precision, dimension(1:n):: xa,ya
      integer:: i,m,ns
      double precision:: den,dif,dift,ho,hp,w
      double precision, dimension(1:nmax):: c,d

      if(n.gt.nmax) then
        write(*,*)'  !! interpolation polynomial order too high!'
        stop
      endif

      ns=1
      dif=abs(x-xa(1))

      do i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo

      y=ya(ns)
      ns=ns-1

      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0d0) then
            print*, 'intern',xa(i),xa(i+m)
            write(*,*)'failure in interpolation!'
            stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if(2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo

end subroutine interpolate

