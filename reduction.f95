!Reduce vector using some floating point operations
real*8 function reduce(vec,n)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: vec(n)
  integer :: i

  reduce = 0.
  !Loop through indices, do operations and add to sum
  !$omp parallel do reduction(+:reduce)
  do i=1,n
    reduce = reduce + vec(i)*vec(i)/2. + vec(i)/10.
  end do
  !$omp end parallel do

end function reduce
