SUBROUTINE allocate_jcs

  USE para
  USE Lagrange
  USE Euler

  ! principle jump conditions
  ALLOCATE(ujc(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(vjc(1:6, 1:nvertex4proc(nobj4proc)))
  ALLOCATE(pjc(1:8, 1:nvertex4proc(nobj4proc)))
  ! xf
  ALLOCATE(ujcxf(1:6,1:n_xf_int))
  ALLOCATE(vjcxf(1:6,1:n_xf_int))
  ! xc
  ALLOCATE(vjcxc(1:6,1:n_xc_int))
  ALLOCATE(ujcxc(1:6,1:n_xc_int))
  ALLOCATE(pjcxc(1:8,1:n_xc_int))
  ! yf
  ALLOCATE(ujcyf(1:6,1:n_yf_int))
  ALLOCATE(vjcyf(1:6,1:n_yf_int))
  ! yc
  ALLOCATE(vjcyc(1:6,1:n_yc_int))
  ALLOCATE(ujcyc(1:6,1:n_yc_int))
  ALLOCATE(pjcyc(1:8,1:n_yc_int))

END SUBROUTINE allocate_jcs
