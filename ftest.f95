!Double loop sum
real*8 function doubleloop(inp,n)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: inp(n,n)
  real*8 :: dum,b,temp(n)
  integer :: i,j

  dum = 0
 
  !$omp parallel do private(i,j,b) reduction(+:dum)
  do i=1,n
    temp = inp(i,:)
    do j=1,n
      dum = dum + temp(j)
    end do
  end do
  !$omp end parallel do

  doubleloop = dum

end function doubleloop


!Double loop sum
real*8 function singleloop(inp,n)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: inp(n)
  real*8 :: dum
  integer :: i

  dum = 0
 
  !$omp parallel do private(i) reduction(+:dum)
  do i=1,n
    dum = dum + inp(i)
  end do
  !$omp end parallel do

  singleloop = dum

end function singleloop
