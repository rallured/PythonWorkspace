!Merit function calculation for fmin_slsqp optimizer
real*8 function merit(voltages,distortion,ifuncs,n,m)
	implicit none
	integer, intent(in):: n,m
	real*8, intent(in) :: voltages(n), distortion(m), ifuncs(m,n)
	integer :: i,j
	real*8 :: tmp, tmp2

	tmp2 = 0.

	!Outer loop for spatial dimension
	!$omp parallel do private(i,j,tmp) reduction(+:tmp2)
	do i=1,m
		tmp = 0.
		do j=1,n
			tmp = tmp + voltages(j)*ifuncs(i,j)
		end do
		tmp2 = tmp2 + (tmp - distortion(i))**2
	end do
	!$omp end parallel do

	merit = tmp2 / m

end function merit
