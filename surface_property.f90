SUBROUTINE surface_property

  USE mpi
  USE Lagrange
  DOUBLE PRECISION:: gacobi,taoxx,taoyy
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::coord

  ALLOCATE(coord(1:2,1:nvertex4proc(nobj4proc)))
  DO myobj=1,nobj4proc
     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)
        coord(1,myvertex) = vertex(1, myvertex)
        coord(2,myvertex) = vertex(2, myvertex)
     ENDDO

     DO myvertex=nvertex4proc(myobj-1)+1,nvertex4proc(myobj)-1
        taoxx = coord(1,myvertex+1)-coord(1,myvertex)
        taoyy = coord(2,myvertex+1)-coord(2,myvertex)
        gacobi=SQRT(taoxx*taoxx+taoyy*taoyy)
        taox(myvertex) = taoxx/gacobi
        taoy(myvertex) = taoyy/gacobi
        !print*, myvertex, taox(myvertex),taoy(myvertex)
     ENDDO
     taox(nvertex4proc(myobj)) = taox(nvertex4proc(myobj-1)+1)
     taoy(nvertex4proc(myobj)) = taoy(nvertex4proc(myobj-1)+1)
  ENDDO

  DEALLOCATE(coord)


END SUBROUTINE surface_property
