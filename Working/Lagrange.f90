module Lagrange
  ! local object
  integer:: nlocobj
  integer, dimension(:), allocatable:: nvertex4obj,npanel4obj
  integer, dimension(:), allocatable:: obj4proc
  integer, dimension(:,:), allocatable:: proc4obj
  integer, dimension(:,:), allocatable:: panel
  double precision, dimension(:,:), allocatable:: vertex,curvature

  ! panel-gridline intersections
  integer:: nicjc,nifjc,nicjf,nickc,nifkc,nickf,njckc,njfkc,njckf
  integer, dimension(:), allocatable:: k4icjc,k4ifjc,k4icjf,k4ickc,k4ifkc,k4ickf,k4jckc,k4jfkc,k4jckf

  ! Cartesian jump conditions
  double precision, dimension(:), allocatable:: jumpdudx,jumpdudy,jumpdudz
  double precision, dimension(:), allocatable:: jumpdvdx,jumpdvdy,jumpdvdz
  double precision, dimension(:), allocatable:: jumpdwdx,jumpdwdy,jumpdwdz
  double precision, dimension(:), allocatable:: jumpd2udx2,jumpd2udy2,jumpd2udz2,jumpd2udxy,jumpd2udxz,jumpd2udyz
  double precision, dimension(:), allocatable:: jumpd2vdx2,jumpd2vdy2,jumpd2vdz2,jumpd2vdxy,jumpd2vdxz,jumpd2vdyz
  double precision, dimension(:), allocatable:: jumpd2wdx2,jumpd2wdy2,jumpd2wdz2,jumpd2wdxy,jumpd2wdxz,jumpd2wdyz
end module Lagrange
