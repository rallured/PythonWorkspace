!Factorial Function
real*8 function factorial(n)
    implicit none
    integer, intent(in) :: n
    real*8 :: Ans
    integer :: i

    Ans = 1
    do i=1,n
        Ans = Ans * dble(i)
    end do

    factorial = Ans
end function factorial

!Function to compute radial Zernike polynomial at radius rho of order n,m
real*8 function radialpoly(rho,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho
    integer, intent(in) :: n,m
    real*8 :: output
    real*8 :: const, factorial
    integer :: j

    !n,m are assumed to be valid (m>=n, n-abs(m) is even)
    !Following Python module convention, negative m indicates sinusoidal Zernike
    output = 0.
    if (rho<=1) then
        do j = 0,(n-abs(m))/2
            const = (-1)**j*factorial(n-j)/factorial(j)/factorial((n+m)/2-j)/factorial((n-m)/2-j)
            output = output + const*rho**(n-2*j)
        end do
    end if

    radialpoly = output

end function radialpoly

!This function makes use of the q recursive method to compute a Zernike
!Radial polynomial
recursive function fastradpoly(rho,n,m) result(output)
  !Declarations
  implicit none
  real*8, intent(in) :: rho,n,m
  !integer, intent(in) :: n,m
  real*8 :: output,r2,r4,h1,h2,h3,tempm

  !Handle rho=0
  if (rho==0) then
    if (m==0) then
      output = (-1.)**(n/2)
      return
    else
      output = 0.
      return
    end if
  end if

  !Three possibilities, m=n, m=n-2, m<n-2
  if (m==n) then
    output = rho**n
    return
  else if (m==n-2) then
    output = n*rho**n - (n-1)*rho**(n-2)
    return
  else
    tempm = n-4
    r4 = rho**n
    r2 = n*rho**n - (n-1)*rho**(n-2)
    do while(tempm >= m)
      h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
      h2 = h3*(n+tempm+4)*(n-tempm-2)/4./(tempm+3) + (tempm+2)
      h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.
      output = h1*r4 + (h2 + h3/rho**2)*r2
      if (tempm == m) then
        return
      end if
      r4 = r2
      r2 = output
      tempm = tempm - 2
    end do
  end if

end function fastradpoly

!This function makes use of the q recursive method to compute a Zernike
!Radial polynomial DERIVATIVE
recursive function fastradder(rho,n,m) result(output)
  !Declarations
  implicit none
  real*8, intent(in) :: rho
  real*8, intent(in) :: n,m
  real*8 :: output,r2,r4,h1,h2,h3,fastradpoly,tempm

  !Handle rho=0
  if (rho==0) then
    if (m==1) then
      output = (-1.)**((n-1)/2)*(n+1)/2
      return
    else
      output = 0.
      return
    end if
  end if  

  !Three possibilities, m=n, m=n-2, m<n-2
  if (m==n) then
    output = n*rho**(n-1)
    return
  else if (m==n-2) then
    output = n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
    return
  else
    tempm = n-4
    r4 = n*rho**(n-1)
    r2 = n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
    do while(tempm >= m)
      h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
      h2 = h3*(n+tempm+4)*(n-tempm-2)/4./(tempm+3) + (tempm+2)
      h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.
      !Might need to literally copy and paste fastradpoly here to save on function call...
      output = h1*r4 + (h2 + h3/rho**2)*r2 - 2*h3*fastradpoly(rho,n,tempm+2)/rho**3
      if (tempm == m) then
        return
      end if
      r4 = r2
      r2 = output
      tempm = tempm - 2
    end do
  end if

end function fastradder

!This function will return a vector of Zernike polynomial evaluations for
!a point (rho,theta) and a number N of standard Zernike polynomials
!Start at n=0, and loop up in radial order
!For each radial order, start at Znn, and recursively calculate down
!Adding azimuthal dependence is easy, just note what overall index you're at
!and add sinmt if odd and cosmt if even, where Z1 = 1, so Z2=2rcost
!Still need special handling of rho=0
subroutine zernset(rho,theta,rorder,aorder,znum,polyout,derrho,dertheta)
  !Declarations
  implicit none
  integer, intent(in) :: znum
  integer, intent(in) :: rorder(znum),aorder(znum)
  real*8, intent(in) :: rho,theta
  real*8, intent(out) :: polyout(znum),derrho(znum),dertheta(znum)
  real*8, allocatable :: rnm(:,:),rprime(:,:)
  integer :: i,j,radnum,tznum
  real*8 :: n,m,mm,h1,h2,h3,norm
  real*8 :: fastradpoly,fastradder,zernthetader,zernrhoder,zernike

  !Allocate temp output array for full sets of radial polynomials
  !Will slice off unused terms at the end
  tznum = 1
  radnum = 1
  do while(tznum < znum)
    tznum = tznum + (radnum + 1)
    radnum = radnum + 1
  end do
  !print *, tznum, " ", radnum
  !Allocate radnum by radnum arrays for radial polynomials
  !This is much larger than necessary, but for practical purposes this should be ok
  allocate(rnm(radnum,radnum))
  allocate(rprime(radnum,radnum))

  !Compute radial polynomials and derivatives
  do i=1,radnum
    n = dble(i)-1
    do j=int(n)+1,1,-2
      m = dble(j)-1
      !Handle case rho=0
      if (rho==0) then
        if (m==1) then
          rprime(i,j) = (-1.)**((n-1)/2)*(n+1)/2
          rnm(i,j) = 0.
        else if (m==0) then
          rprime(i,j) = 0.
          rnm(i,j) = (-1.)**(n/2)
        else
          rprime(i,j) = 0.
          rnm(i,j) = 0.
        end if
      !Handle general case
      else
        if (n==m) then
          rnm(i,j) = rho**n
          rprime(i,j) = n*rho**(n-1)
        else if (m==n-2) then
          rnm(i,j) = n*rnm(i,i) - (n-1)*rnm(i-2,i-2)
          rprime(i,j) = n*rprime(i,i) - (n-1)*rprime(i-2,i-2)
        else
          h3 = -4*(m+2)*(m+1)/(n+m+2)/(n-m)
          h2 = h3*(n+m+4)*(n-m-2)/4./(m+3) + (m+2)
          h1 = .5*(m+4)*(m+3) - (m+4)*h2 + h3*(n+m+6)*(n-m-4)/8.
          rnm(i,j) = h1*rnm(i,j+4) + (h2+h3/rho**2)*rnm(i,j+2)
          rprime(i,j) = h1*rprime(i,j+4) + (h2+h3/rho**2)*rprime(i,j+2) - 2*h3/rho**3*rnm(i,j+2)
        end if
      end if
      !print *, i-1, " ", j-1
      !print *, rnm(i,j), " ", fastradpoly(rho,dble(i-1),dble(j-1))
      !print *, rprime(i,j), " ", fastradder(rho,dble(i-1),dble(j-1))
    end do
  end do

  !Radial terms are computed. Now construct Zernike, Zernderrho, and Zerndertheta
  !Use standard order up to znum
  do i=1,znum
    !Rorder and Aorder are passed to this function
    !Simply construct if statement on m, and compute all three!
    n = rorder(i)
    mm = aorder(i)
    m = abs(mm)
    norm = sqrt(2*(n+1))
    if (mm<0) then
        polyout(i) = norm * rnm(int(n)+1,int(m)+1) * sin(m*theta)
        derrho(i) = norm * rprime(int(n)+1,int(m)+1) * sin(m*theta)
        dertheta(i) = norm * rnm(int(n)+1,int(m)+1) * cos(m*theta) * m
    else if (mm>0) then
        polyout(i) = norm * rnm(int(n)+1,int(m)+1) * cos(m*theta)
        derrho(i) = norm * rprime(int(n)+1,int(m)+1) * cos(m*theta)
        dertheta(i) = -norm * rnm(int(n)+1,int(m)+1) * sin(m*theta) * m
    else
        polyout(i) = norm*sqrt(0.5)*rnm(int(n)+1,int(m)+1)
        derrho(i) = norm*sqrt(0.5)*rprime(int(n)+1,int(m)+1)
        dertheta(i) = 0.
    end if
    !print *, dertheta(i), " ",zernthetader(rho,theta,int(n),int(mm))
  end do

end subroutine zernset

!Function to compute full Zernike polynomial at radius rho and angle theta of order n,m
!Following Python Zernike mod conventions
real*8 function zernike(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: fastradpoly,rnm,znm,norm

    !Compute radial polynomial
    rnm = fastradpoly(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * sin(m*theta)
    else if (m>0) then
        znm = norm * rnm * cos(m*theta)
    else
        znm = norm*sqrt(0.5)*rnm
    end if

    !Return Zernike polynomial
    zernike = znm

end function zernike

!Function to compute radial polynomial derivative
real*8 function radialder(rho,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho
    integer, intent(in) :: n,m
    real*8 :: output
    real*8 :: const, factorial
    integer :: j
    !n,m are assumed to be valid (m>=n, n-abs(m) is even)
    !Following Python module convention, negative m indicates sinusoidal Zernike
    output = 0.
    if (rho<=1) then
        do j = 0,(n-abs(m))/2
            const = (-1)**j*factorial(n-j)/factorial(j)/factorial((n+m)/2-j)/factorial((n-m)/2-j)
            output = output + const*(n-2*j)*rho**(n-2*j-1)
        end do
    end if

    radialder = output

end function radialder

!Function to compute partial Zernike derivative wrt rho
real*8 function zernrhoder(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: fastradder, rnm, norm, znm

    !Compute radial polynomial
    rnm = fastradder(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * sin(abs(m)*theta)
    else if (m>0) then
        znm = norm * rnm * cos(m*theta)
    else
        znm = norm*sqrt(0.5)*rnm
    end if

    !Return Zernike polynomial
    zernrhoder = znm

end function zernrhoder

!Function to compute partial Zernike derivative wrt theta
real*8 function zernthetader(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: rnm, znm, fastradpoly, norm

    !Compute radial polynomial
    rnm = fastradpoly(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * abs(m) * cos(abs(m)*theta)
    else if (m>0) then
        znm = -norm * rnm * m * sin(m*theta)
    else
        znm = 0.
    end if

    !Return Zernike polynomial
    zernthetader = znm

end function zernthetader

!Function to trace set of rays to Zernike surface
!Inputs are position and cosines of rays
!Array of Zernike coefficients
!Output is position of rays at surface
!and new cosines after reflection
subroutine tracezern(x,y,z,l,m,n,ux,uy,uz,num,coeff,rorder,aorder,arrsize,rad)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize, num
  real*8, intent(in) :: coeff(arrsize),rad
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer, intent(in) :: rorder(arrsize),aorder(arrsize)
  integer :: i,c
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,rho,theta,delta,t
  real*8 :: zern(arrsize),rhoder(arrsize),thetader(arrsize)
  real*8 :: zernike,zernrhoder,zernthetader
  real*8 :: dum

  !Loop through each individual ray, trace to surface, and reflect
  !Establish convergence criteria (surface function should be < 1.e-12)
  !$omp parallel do private(t,delta,rho,theta,F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,i,c,zern,rhoder,thetader)
  do i=1,num
    t = 0.
    delta = 100.
    !count = 0
    !print *, x(i), y(i), z(i)
    !print *, l(i), m(i), n(i)
    !read *, dum
    do while (abs(delta) > 1.e-10)
      !count = count + 1
      !Convert to cylindrical coordinates
      rho = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      !Compute surface function and derivatives
      F = z(i)
      Frho = 0.
      Ftheta = 0.
      call zernset(rho/rad,theta,rorder,aorder,arrsize,zern,rhoder,thetader)
      !print *, zern
      !read *, dum
      do c=1,arrsize
        !F = F - coeff(c)*zernike(rho/rad,theta,rorder(c),aorder(c))
        !Frho = Frho - coeff(c)*zernrhoder(rho/rad,theta,rorder(c),aorder(c))/rad
        !Ftheta = Ftheta - coeff(c)*zernthetader(rho/rad,theta,rorder(c),aorder(c))
        F = F - coeff(c)*zern(c)
        Frho = Frho - coeff(c)*rhoder(c)/rad
        Ftheta = Ftheta - coeff(c)*thetader(c)
        !print *, rho/rad, " ", theta, " ", rorder(c), " ", aorder(c)
        !print *, zern(c), rhoder(c)
        !print *, coeff(c)
        !print *, F
        !read *, dum
      end do
      !print *, rho/rad
      !print *, theta
      !print *, zernrhoder(rho/rad,theta,rorder(c),aorder(c))
      !print *, zernthetader(rho/rad,theta,rorder(c),aorder(c))
      !Convert back to cartesian coordinates
      !print *, Frho, Ftheta
      !read *, dum
      Frhox = (x(i)/rho) * Frho
      Frhoy = (y(i)/rho) * Frho
      Fthetax = (-y(i)/rho) * Ftheta/rho
      Fthetay = (x(i)/rho) * Ftheta/rho
      !if (rho==0) then
      !  Frhox = Frho
      !  Frhoy = 0.
      !  Fthetax = 0.
      !  Fthetay = 0.
      !end if
      Fx = Frhox + Fthetax
      Fy = Frhoy + Fthetay
      Fz = 1.
      !print *, Fx, Fy, Fz
      !print *, l(i), m(i), n(i)
      !read *, dum
      !Compute delta and new ray position
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      !print *, 'Fp:' , Fp
      !read *, dum
      delta = -F/Fp
      !print *, F, Fp, delta
      !read *, dum
      x(i) = x(i) + l(i)*delta
      y(i) = y(i) + m(i)*delta
      z(i) = z(i) + n(i)*delta
      !Keep track of total length ray is traced
      t = t + delta
    end do
    !print *, "Number of iterations: ", count
    !We have converged, do any vignetting and compute surface normal
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
  end do
  !$omp end parallel do

end subroutine tracezern

!Trace a Zernike standard phase surface as in ZEMAX
subroutine zernphase(opd,x,y,z,l,m,n,ux,uy,uz,num,coeff,rorder,aorder,arrsize,rad,wave)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize, num
  real*8, intent(in) :: coeff(arrsize),rad,wave
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num),opd(num)
  integer, intent(in) :: rorder(arrsize),aorder(arrsize)
  real*8 :: zern(arrsize),rhoder(arrsize),thetader(arrsize)
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy
  real*8 :: pi,theta,rho,dum
  integer :: i,c

  pi = acos(-1.)

  do i=1,num
    !Compute the set of Zernike polynomials
    rho = sqrt(x(i)**2+y(i)**2)
    theta = atan2(y(i),x(i))
    call zernset(rho/rad,theta,rorder,aorder,arrsize,zern,rhoder,thetader)
    !Sum derivatives
    F = 0.
    Frho = 0.
    Ftheta = 0.
    do c=1,arrsize
      F = F + coeff(c)*zern(c)
      Frho = Frho + coeff(c)*rhoder(c)/rad
      Ftheta = Ftheta + coeff(c)*thetader(c)
    end do
    !Convert to Cartesian coordinates
    Frhox = (x(i)/rho) * Frho
    Frhoy = (y(i)/rho) * Frho
    Fthetax = (-y(i)/rho) * Ftheta/rho
    Fthetay = (x(i)/rho) * Ftheta/rho
    Fx = Frhox + Fthetax
    Fy = Frhoy + Fthetay
    !print *, rho,theta,Frho, Ftheta,Fx,Fy,wave
    !read *, dum
    !Add appropriate perturbations to ray wavevector
    l(i) = l(i) + Fx*wave!/2/pi
    m(i) = m(i) + Fy*wave!/2/pi
    n(i) = sign(sqrt(1.-l(i)**2-m(i)**2),n(i))
    opd(i) = opd(i) + F*wave
  end do

end subroutine zernphase

!Function to trace set of rays to Zernike surface
!Inputs are position and cosines of rays
!Array of Zernike coefficients
!Output is position of rays at surface
!and new cosines after reflection
subroutine tracezernrot(x,y,z,l,m,n,ux,uy,uz,num,coeff1,rorder1,aorder1,arrsize1,coeff2,rorder2,aorder2,arrsize2,rad,rot)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize1, arrsize2, num
  real*8, intent(in) :: coeff1(arrsize1),coeff2(arrsize2),rad,rot
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer, intent(in) :: rorder1(arrsize1),aorder1(arrsize1)
  integer, intent(in) :: rorder2(arrsize2),aorder2(arrsize2)
  integer :: i,c
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,rho,theta,delta,t
  real*8 :: zern1(arrsize1),rhoder1(arrsize1),thetader1(arrsize1)
  real*8 :: zern2(arrsize2),rhoder2(arrsize2),thetader2(arrsize2)
  real*8 :: zernike,zernrhoder,zernthetader
  real*8 :: dum

  !Loop through each individual ray, trace to surface, and reflect
  !Establish convergence criteria (surface function should be < 1.e-12)
  !$omp parallel do private(t,delta,rho,theta,F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,i,c,zern1,zern2)
  do i=1,num
    t = 0.
    delta = 100.
    !count = 0
    !print *, x(i), y(i), z(i)
    !print *, l(i), m(i), n(i)
    !read *, dum
    do while (abs(delta) > 1.e-10)
      !count = count + 1
      !Convert to cylindrical coordinates
      rho = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      !Compute surface function and derivatives
      F = z(i)
      Frho = 0.
      Ftheta = 0.
      call zernset(rho/rad,theta,rorder1,aorder1,arrsize1,zern1,rhoder1,thetader1)
      !print *, zern1(1)
      !read *, dum
      do c=1,arrsize1
        !F = F - coeff(c)*zernike(rho/rad,theta,rorder(c),aorder(c))
        !Frho = Frho - coeff(c)*zernrhoder(rho/rad,theta,rorder(c),aorder(c))/rad
        !Ftheta = Ftheta - coeff(c)*zernthetader(rho/rad,theta,rorder(c),aorder(c))
        F = F - coeff1(c)*zern1(c)
        Frho = Frho - coeff1(c)*rhoder1(c)/rad
        Ftheta = Ftheta - coeff1(c)*thetader1(c)
        !print *, rho/rad, " ", theta, " ", rorder(c), " ", aorder(c)
        !print *, zern(c), rhoder(c)
        !print *, coeff(c)
        !print *, F
        !read *, dum
      end do
      call zernset(rho/rad,theta+rot,rorder2,aorder2,arrsize2,zern2,rhoder2,thetader2)
      !print *, zern2(1)
      do c=1,arrsize2
        F = F - coeff2(c)*zern2(c)
        Frho = Frho - coeff2(c)*rhoder2(c)/rad
        Ftheta = Ftheta - coeff2(c)*thetader2(c)
      end do
      !print *, rho/rad
      !print *, theta
      !print *, zernrhoder(rho/rad,theta,rorder(c),aorder(c))
      !print *, zernthetader(rho/rad,theta,rorder(c),aorder(c))
      !Convert back to cartesian coordinates
      !print *, Frho, Ftheta
      !read *, dum
      Frhox = (x(i)/rho) * Frho
      Frhoy = (y(i)/rho) * Frho
      Fthetax = (-y(i)/rho) * Ftheta/rho
      Fthetay = (x(i)/rho) * Ftheta/rho
      !if (rho==0) then
      !  Frhox = Frho
      !  Frhoy = 0.
      !  Fthetax = 0.
      !  Fthetay = 0.
      !end if
      Fx = Frhox + Fthetax
      Fy = Frhoy + Fthetay
      Fz = 1.
      !print *, Fx, Fy, Fz
      !print *, l(i), m(i), n(i)
      !read *, dum
      !Compute delta and new ray position
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      !print *, 'Fp:' , Fp
      !read *, dum
      delta = -F/Fp
      !print *, F, Fp, delta
      !read *, dum
      x(i) = x(i) + l(i)*delta
      y(i) = y(i) + m(i)*delta
      z(i) = z(i) + n(i)*delta
      !Keep track of total length ray is traced
      t = t + delta
    end do
    !print *, "Number of iterations: ", count
    !We have converged, do any vignetting and compute surface normal
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
  end do
  !$omp end parallel do

end subroutine tracezernrot

!This function computes a Legendre polynomial of order n
real*8 function legendre(x,n)
  !Declarations
  implicit none
  real*8, intent(in) :: x
  integer, intent(in) :: n
  integer :: i
  real*8 :: factorial,x2

  if (abs(x) > 1.) then
    x2 = x/abs(x)
  else
    x2 = x
  end if

  legendre = 0.
  if (n==0) then
    legendre = 1.
  else
    do i=0,floor(real(n)/2)
      legendre = legendre + (-1)**(i)*factorial(2*n-2*i)/factorial(i)/factorial(n-i)/factorial(n-2*i)/2**n*x2**(n-2*i)
    end do
  end if

end function legendre

!This function computes a Legendre polynomial first derivative of order n
real*8 function legendrep(x,n)
  !Declarations
  implicit none
  real*8, intent(in) :: x
  integer, intent(in) :: n
  integer :: i
  real*8 :: factorial

  legendrep = 0.
  if (n==0) then
    legendrep = 0.
  else if (n==1) then
    legendrep = 1.
  else if (x==0. .and. mod(n,2)==0) then
    legendrep = 0.
  else
    do i=0,floor(real(n)/2)
      legendrep = legendrep + (-1)**(i)*factorial(2*n-2*i)/factorial(i)/factorial(n-i)/factorial(n-2*i)/2**n*(n-2*i)*x**(n-2*i-1)
    end do
  end if  

  if (abs(x) > 1.) then
    legendrep = 0.
  end if

end function legendrep

!This function rotates a vector in a right handed fashion
!axis is 1,2,3 for x,y,z axis rotation
subroutine rotatevector(x,y,z,theta,axis)
  !Declarations
  real*8, intent(inout) :: x,y,z,theta
  integer, intent(in) :: axis
  real*8 :: output(3)
  !real*8, dimension(3) :: rotatevector

  if (axis==1) then
    output(1) = x
    output(2) = cos(theta)*y-sin(theta)*z
    output(3) = sin(theta)*y+cos(theta)*z
  else if (axis==2) then
    output(1) = cos(theta)*x+sin(theta)*z
    output(2) = y
    output(3) = -sin(theta)*x+cos(theta)*z
  else
    output(1) = cos(theta)*x-sin(theta)*y
    output(2) = sin(theta)*x+cos(theta)*y
    output(3) = z
  end if

  x = output(1)
  y = output(2)
  z = output(3)

end subroutine rotatevector

!This function rotates a vector in a right handed fashion
!axis is given by input ux,uy,uz
subroutine rotateaxis(x,y,z,theta,ux,uy,uz)
  !Declarations
  real*8, intent(inout) :: x,y,z
  real*8, intent(inout) :: theta,ux,uy,uz
  real*8 :: output(3),s,c,mag
  !real*8, dimension(3) :: rotatevector
  
  !Ensure axis is normalized
  mag = sqrt(ux**2 + uy**2 + uz**2)
  ux = ux/mag
  uy = uy/mag
  uz = uz/mag

  c = cos(theta)
  s = sin(theta)
  output(1) = (c+ux**2*(1-c))*x + (ux*uy*(1-c)-uz*s)*y + (ux*uz*(1-c)+uy*s)*z
  output(2) = (uy*ux*(1-c)+uz*s)*x + (c+uy**2*(1-c))*y + (uy*uz*(1-c)-ux*s)*z
  output(3) = (uz*ux*(1-c)-uy*s)*x + (uz*uy*(1-c)+ux*s)*y + (c+uz**2*(1-c))*z
  
  x = output(1)
  y = output(2)
  z = output(3)

end subroutine rotateaxis

!Function to reflect about local surface normal
!Assume ray cosines point into surface, surface normal points out of surface
!If i = -incident, then reflected ray is 2(i.n)n-i
subroutine reflect(l,m,n,ux,uy,uz,num)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8 :: dot
  integer :: i

  !Loop through rays and reflect about normal
  !$omp parallel do private(i,dot)
  do i=1,num
    !Compute dot of incident with normal
    dot = ux(i)*l(i) + uy(i)*m(i) + uz(i)*n(i)
    !Compute reflected direction
    l(i) = l(i) - 2*dot*ux(i)
    m(i) = m(i) - 2*dot*uy(i)
    n(i) = n(i) - 2*dot*uz(i)
  end do
  !$omp end parallel do

end subroutine reflect

!Refract from one index to another
subroutine refract(l,m,n,ux,uy,uz,num,n1,n2)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: n1,n2
  real*8, intent(inout) :: l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: dot,t1,t2,alpha,cx,cy,cz
  
  !Loop through rays
  !$omp parallel do private(dot,t1,t2,cx,cy,cz,alpha)
  do i=1,num
    !Ensure normal vector is pointing into second index (dot product should be positive)
    dot = l(i)*ux(i) + m(i)*uy(i) + n(i)*uz(i)
    if (dot < 0) then
      ux(i) = -ux(i)
      uy(i) = -uy(i)
      uz(i) = -uz(i)
      dot = -dot
    end if
    !Compute Snell's law
    t1 = acos(dot)
    t2 = asin((n1/n2)*sin(t1))
    !Compute cross product
    cx = uy(i)*n(i)-m(i)*uz(i)
    cy = l(i)*uz(i)-ux(i)*n(i)
    cz = ux(i)*m(i)-l(i)*uy(i)
    !Rotate about cross vector an angle t2-t1
    call rotateaxis(l(i),m(i),n(i),t2-t1,cx,cy,cz)
    !Normalize
    alpha = sqrt(l(i)**2 + m(i)**2 + n(i)**2)
    l(i) = l(i)/alpha
    m(i) = m(i)/alpha
    n(i) = n(i)/alpha
  end do
  !$omp end parallel do

end subroutine refract

!Coordinate system transform, translations are done first, then rotations in XYZ order
!Rotations act to rotate a surface via the right hand rule
subroutine transform(x,y,z,l,m,n,ux,uy,uz,num,tx,ty,tz,rx,ry,rz)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: tx,ty,tz,rx,ry,rz
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i

  !Loop through rays
  !$omp parallel do
  do i=1,num
    !Perform translation
    x(i) = x(i) + tx
    y(i) = y(i) + ty
    z(i) = z(i) + tz
    !Perform x rotation
    call rotatevector(x(i),y(i),z(i),rx,1)
    call rotatevector(l(i),m(i),n(i),rx,1)
    call rotatevector(ux(i),uy(i),uz(i),rx,1)
    !Perform y rotation
    call rotatevector(x(i),y(i),z(i),ry,2)
    call rotatevector(l(i),m(i),n(i),ry,2)
    call rotatevector(ux(i),uy(i),uz(i),ry,2)
    !Perform z rotation
    call rotatevector(x(i),y(i),z(i),rz,3)
    call rotatevector(l(i),m(i),n(i),rz,3)
    call rotatevector(ux(i),uy(i),uz(i),rz,3)
  end do
  !$omp end parallel do

end subroutine transform

!Coordinate system transform, rotations in ZYX order, then translations
!Rotations act to rotate a surface via the right hand rule
!This is meant to invert the transform routine
subroutine itransform(x,y,z,l,m,n,ux,uy,uz,num,tx,ty,tz,rx,ry,rz)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: tx,ty,tz,rx,ry,rz
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i

  !Loop through rays
  !$omp parallel do
  do i=1,num
    !Perform z rotation
    call rotatevector(x(i),y(i),z(i),-rz,3)
    call rotatevector(l(i),m(i),n(i),-rz,3)
    call rotatevector(ux(i),uy(i),uz(i),-rz,3)
    !Perform y rotation
    call rotatevector(x(i),y(i),z(i),-ry,2)
    call rotatevector(l(i),m(i),n(i),-ry,2)
    call rotatevector(ux(i),uy(i),uz(i),-ry,2)
    !Perform x rotation
    call rotatevector(x(i),y(i),z(i),-rx,1)
    call rotatevector(l(i),m(i),n(i),-rx,1)
    call rotatevector(ux(i),uy(i),uz(i),-rx,1)
    !Perform translation
    x(i) = x(i) - tx
    y(i) = y(i) - ty
    z(i) = z(i) - tz
  end do
  !$omp end parallel do

end subroutine itransform

!Trace rays to local XY plane
subroutine flat(x,y,z,l,m,n,ux,uy,uz,num)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: l(num),m(num),n(num)
  real*8, intent(inout) :: x(num),y(num),z(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: dummy

  !Loop through rays
  !$omp parallel do private(delta)
  do i=1,num
    delta = -z(i)/n(i)
    z(i) = 0.
    x(i) = x(i) + delta*l(i)
    y(i) = y(i) + delta*m(i)
    ux(i) = 0.
    uy(i) = 0.
    uz(i) = 1.
    !print *, x(i),y(i),z(i)
    !print *, l(i),m(i),n(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dummy
  end do
  !$omp end parallel do

end subroutine flat

!Trace rays to local XY plane
subroutine flatOPD(x,y,z,l,m,n,ux,uy,uz,opd,num,nr)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: l(num),m(num),n(num),nr
  real*8, intent(inout) :: x(num),y(num),z(num),ux(num),uy(num),uz(num),opd(num)
  integer :: i

  !Loop through rays
  !$omp parallel do private(delta)
  do i=1,num
    delta = -z(i)/n(i)
    z(i) = 0.
    x(i) = x(i) + delta*l(i)
    y(i) = y(i) + delta*m(i)
    ux(i) = 0.
    uy(i) = 0.
    uz(i) = 1.
    opd(i) = opd(i) + delta*nr
  end do
  !$omp end parallel do

end subroutine flatOPD

!Trace rays to spherical surface, center assumed to be at origin
!Intersection taken to be the closest point to ray
subroutine tracesphere(x,y,z,l,m,n,ux,uy,uz,num,rad)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: rad
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: mago, dotol, determinant, d1, d2

  !Loop through rays
  !$omp parallel do private(dotol,mago,determinant,d1,d2)
  do i=1,num
    !Compute dot product
    dotol= l(i)*x(i) + m(i)*y(i) + n(i)*z(i)
    mago = x(i)**2. + y(i)**2. + z(i)**2.
    !Compute distance to move rays
    determinant = dotol**2 - mago + rad**2
    !If ray does not intersect, set position and cosine vector to NaN
    if (determinant < 0) then
      x(i) = 0.
      y(i) = 0.
      z(i) = 0.
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      d1 = -dotol + sqrt(determinant)
      d2 = -dotol - sqrt(determinant)
      if (abs(d2) < abs(d1)) then
        d1 = d2
      end if
      x(i) = x(i) + d1*l(i)
      y(i) = y(i) + d1*m(i)
      z(i) = z(i) + d1*n(i)
    end if
    !Compute surface normal, just normalized position vector
    mago = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
    ux(i) = x(i)/mago
    uy(i) = y(i)/mago
    uz(i) = z(i)/mago
  end do
  !$omp end parallel do

end subroutine tracesphere

!Trace rays to spherical surface, center assumed to be at origin
!Intersection taken to be the closest point to ray
subroutine tracesphereOPD(opd,x,y,z,l,m,n,ux,uy,uz,num,rad,nr)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: rad,nr
  real*8 , intent(inout) :: opd(num),x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: mago, dotol, determinant, d1, d2

  !Loop through rays
  !$omp parallel do private(dotol,mago,determinant,d1,d2)
  do i=1,num
    !Compute dot product
    dotol= l(i)*x(i) + m(i)*y(i) + n(i)*z(i)
    mago = x(i)**2. + y(i)**2. + z(i)**2.
    !Compute distance to move rays
    determinant = dotol**2 - mago + rad**2
    !If ray does not intersect, set position and cosine vector to NaN
    if (determinant < 0) then
      x(i) = 0.
      y(i) = 0.
      z(i) = 0.
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      d1 = -dotol + sqrt(determinant)
      d2 = -dotol - sqrt(determinant)
      if (abs(d2) < abs(d1)) then
        d1 = d2
      end if
      x(i) = x(i) + d1*l(i)
      y(i) = y(i) + d1*m(i)
      z(i) = z(i) + d1*n(i)
      opd(i) = opd(i) + d1*nr
    end if
    !Compute surface normal, just normalized position vector
    mago = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
    ux(i) = x(i)/mago
    uy(i) = y(i)/mago
    uz(i) = z(i)/mago
  end do
  !$omp end parallel do

end subroutine tracesphereOPD


!Traces onto a cylinder
!Center is assumed to be at origin, y axis is cylindrical axis
subroutine tracecyl(x,y,z,l,m,n,ux,uy,uz,num,rad)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: rad
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: a,b,c,mag,d1,d2,det,dum

  !$omp parallel do private(a,b,c,det,mag,d1,d2)
  !Compute a,b,c terms in quadratic solution for distance to move rays
  do i=1,num
    a = l(i)**2 + n(i)**2
    b = 2*(x(i)*l(i)+z(i)*n(i))
    c = x(i)**2 + z(i)**2 - rad**2
    !Compute determinant, if < 0, set ray to 0's
    det = b**2 - 4*a*c
    if (det < 0) then
      x(i) = 0.
      y(i) = 0.
      z(i) = 0.
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      !Find smallest distance to cylinder
      d1 = (-b + sqrt(det))/2/a
      d2 = (-b - sqrt(det))/2/a
      if (abs(d2) < abs(d1)) then
        d1 = d2
      end if
      !Move ray
      x(i) = x(i) + l(i)*d1
      y(i) = y(i) + m(i)*d1
      z(i) = z(i) + n(i)*d1
    end if
    !Compute surface normal
    mag = sqrt(x(i)**2+z(i)**2)
    ux(i) = x(i)/mag
    uz(i) = z(i)/mag
    uy(i) = 0.
  end do
  !$omp end parallel do

end subroutine tracecyl

!Traces onto a cylinder
!Center is assumed to be at origin, y axis is cylindrical axis
subroutine tracecylOPD(opd,x,y,z,l,m,n,ux,uy,uz,num,rad,nr)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: rad,nr
  real*8 , intent(inout) :: opd(num),x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: a,b,c,mag,d1,d2,det,dum

  !$omp parallel do private(a,b,c,det,mag,d1,d2)
  !Compute a,b,c terms in quadratic solution for distance to move rays
  do i=1,num
    a = l(i)**2 + n(i)**2
    b = 2*(x(i)*l(i)+z(i)*n(i))
    c = x(i)**2 + z(i)**2 - rad**2
    !Compute determinant, if < 0, set ray to 0's
    det = b**2 - 4*a*c
    if (det < 0) then
      x(i) = 0.
      y(i) = 0.
      z(i) = 0.
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      !Find smallest distance to cylinder
      d1 = (-b + sqrt(det))/2/a
      d2 = (-b - sqrt(det))/2/a
      if (abs(d2) < abs(d1)) then
        d1 = d2
      end if
      !Move ray
      x(i) = x(i) + l(i)*d1
      y(i) = y(i) + m(i)*d1
      z(i) = z(i) + n(i)*d1
      opd(i) = opd(i) + d1*nr
    end if
    !Compute surface normal
    mag = sqrt(x(i)**2+z(i)**2)
    ux(i) = x(i)/mag
    uz(i) = z(i)/mag
    uy(i) = 0.
  end do
  !$omp end parallel do

end subroutine tracecylOPD

!This routine traces rays to a cylindrical conic surface
!z axis is cylindrical axis, rays have been traced to the xy
!plane, and sag is in the y direction
subroutine cylconic(x,y,z,l,m,n,ux,uy,uz,num,rad,k)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: rad,k
  real*8 :: low,high,dL,dH
  real*8 :: F,Fx,Fy,Fp,delt,dum
  integer :: i

  !Z derivative of surface function is 0 due to cylindrical symmetry
  !This is essentially a 2D problem

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fp,low,high,dL,dH)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      low = 1 + sqrt(1 - (1+k)*rad**2*x(i)**2)
      high = rad*x(i)**2
      dL = -(1+k)*rad**2*x(i) / sqrt(1 - (1+k)*rad**2*x(i)**2)
      dH = 2*rad*x(i)
      F = y(i) - high/low
      Fx = (high*dL - low*dH) / low**2
      Fy = 1.
      Fp = Fx*l(i) + Fy*m(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F, Fp
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = 0.
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine cylconic

!This function traces to a Wolter I primary mirror
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine wolterprimary(x,y,z,l,m,n,ux,uy,uz,num,r0,z0)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  Fz = 2*p
  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fp)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-8)
      F = 2*p*z(i) + p**2 + 4*e**2*p*d/(e**2-1) - x(i)**2 - y(i)**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wolterprimary

!This function traces to a Wolter I secondary mirror
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine woltersecondary(x,y,z,l,m,n,ux,uy,uz,num,r0,z0)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-8)
      F = e**2*(d+z(i))**2 - z(i)**2 - x(i)**2 - y(i)**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fz = 2*e**2*(d+z(i)) - 2*z(i)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F, Fx, Fy, Fz
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersecondary

!This function traces to a Wolter I primary mirror with sinusoidal perturbation
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine woltersine(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,amp,freq)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,amp,freq
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      rad = sqrt(x(i)**2+y(i)**2) + amp*sin(2*acos(-1.)*freq*z(i))
      F = 2*p*z(i) + p**2 + 4*e**2*p*d/(e**2-1) - rad**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fz = 2.*p - 2*rad*amp*2*acos(-1.)*freq*cos(2*acos(-1.)*freq*z(i))
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersine

!Construct a paraboloid as in Wolter but with Legendre-Legendre
!deformations. Define Legendre and Legendre derivative functions.
!Pass in coeff, axial order, and azimuthal order as in Zemax implementation
subroutine wolterprimLL(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum)
  !Declarations
  implicit none
  integer, intent(in) :: num,cnum
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,zmax,zmin,dphi,coeff(cnum)
  integer, intent(in) :: axial(cnum),az(cnum)
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad,add,addx,addy,addz,G
  integer :: i,a
  real*8 :: pi,legendre,legendrep,zarg,targ,ang

  !Compute Van Speybroeck parameters
  pi = acos(-1.)
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad,ang,zarg,targ,add,addx,addy,addz,a,G)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      ang = atan2(y(i),x(i))
      zarg = (z(i)-((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
      targ = 2*ang/dphi
      
      !Compute Legendre additive terms
      add = 0.
      addx = 0.
      addy = 0.
      addz = 0.
      do a=1,cnum
        add = add + coeff(a)*legendre(zarg,axial(a))*legendre(targ,az(a))
        addx = addx - coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(y(i)/(y(i)**2+x(i)**2))
        addy = addy + coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(x(i)/(y(i)**2+x(i)**2))
        addz = addz + coeff(a)*legendrep(zarg,axial(a))*legendre(targ,az(a))*2/(zmax-zmin)
      end do
      G = sqrt(x(i)**2+y(i)**2) + add
      F = -(G**2 - p**2 - 2*p*z(i) - 4*e**2*p*d/(e**2-1))
      Fx = -2*G*(x(i)/sqrt(x(i)**2+y(i)**2)+addx)
      Fy = -2*G*(y(i)/sqrt(x(i)**2+y(i)**2)+addy)
      Fz = 2*p - 2*G*addz
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Set rays outside mirror definition to NaN
    !if (abs(zarg)>1 .or. abs(targ) > 1) then
    !  x(i) = 0.
    !  y(i) = 0.
    !  z(i) = 0.
    !end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wolterprimLL

!Construct a hyperboloid as in Wolter but with Legendre-Legendre
!deformations. Define Legendre and Legendre derivative functions.
!Pass in coeff, axial order, and azimuthal order as in Zemax implementation
subroutine woltersecLL(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum)
  !Declarations
  implicit none
  integer, intent(in) :: num,cnum
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,zmax,zmin,dphi,coeff(cnum)
  integer, intent(in) :: axial(cnum),az(cnum)
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad,add,addx,addy,addz,G
  integer :: i,a
  real*8 :: pi,legendre,legendrep,zarg,targ,ang

  !Compute Van Speybroeck parameters
  pi = acos(-1.)
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad,ang,zarg,targ,add,addx,addy,addz,a,G)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-7)
      ang = atan2(y(i),x(i))
      zarg = (z(i)-((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
      targ = 2*ang/dphi

      !print *, x(i),y(i),z(i)
      !read *, dum
      
      !Compute Legendre additive terms
      add = 0.
      addx = 0.
      addy = 0.
      addz = 0.
      do a=1,cnum
        add = add + coeff(a)*legendre(zarg,axial(a))*legendre(targ,az(a))
        addx = addx - coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(y(i)/(y(i)**2+x(i)**2))
        addy = addy + coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(x(i)/(y(i)**2+x(i)**2))
        addz = addz + coeff(a)*legendrep(zarg,axial(a))*legendre(targ,az(a))*2/(zmax-zmin)
      end do

      G = sqrt(x(i)**2+y(i)**2) + add
      F = -(G**2 - e**2*(d+z(i))**2 + z(i)**2)
      Fx = -2*G*(x(i)/sqrt(x(i)**2+y(i)**2)+addx)
      Fy = -2*G*(y(i)/sqrt(x(i)**2+y(i)**2)+addy)
      Fz = 2*e**2*(d+z(i)) - 2*z(i) - 2*G*addz
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, G, F, Fx, Fy, Fz
      !print *, l(i), m(i), n(i), delt
      !read *, dum
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Set rays outside mirror definition to NaN
    if (abs(zarg)>1 .or. abs(targ) > 1) then
      x(i) = 0.
      y(i) = 0.
      z(i) = 0.
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersecLL

!This function traces to a Wolter-Schwarzschild primary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wsprimary(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi
  real*8 :: k,kterm,dbdx,dbdy,beta,betas,ff,g,r
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,flag,r)
  do i=1,num
    delt = 100.
    c = 0
    do while(abs(delt)>1.e-8)
      beta = asin(sqrt(x(i)**2 + y(i)**2)/ff)
      flag = 0
      if (beta<=betas) then
        beta = betas
        flag = 1
        kterm = 0.
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
      end if
      F = -z(i) - ff*sin(betas/2)**2 + &
          ff**2*sin(beta)**2/(4*ff*sin(betas/2)**2) + &
          g*cos(beta/2)**4*(kterm)**(1-k)
      Fb = ff**2*sin(beta)*cos(beta)/(2*ff*sin(betas/2)**2) - &
           2*g*cos(beta/2)**3*sin(beta/2)*(kterm)**(1-k) + &
           g*(1-k)*cos(beta/2)*sin(beta/2)*(kterm)**(-k)*(1/k)
      Fz = -1.
      if (flag==1) then
        r = sqrt(x(i)**2 + y(i)**2)
        Fb = ff**2*sin(betas)*cos(betas)/(2*ff*sin(betas/2)**2) + &
              g*(1-k)*cos(betas/2)*sin(betas/2)*(1/k)
        F = F + (r - ff*sin(betas))*z(i)/(r**2+z(i)**2)*Fb
        Fz = Fz + (r-ff*sin(betas))*(r**2-z(i)**2)/(r**2+z(i)**2)**2*Fb
        !print *, Fb, F, Fz
        !read *, dum
      end if
      dbdx = x(i)/sqrt(1-(x(i)**2+y(i)**2)/ff**2)/ff/sqrt(x(i)**2+y(i)**2)
      dbdy = y(i)/sqrt(1-(x(i)**2+y(i)**2)/ff**2)/ff/sqrt(x(i)**2+y(i)**2)
      Fx = Fb * dbdx
      Fy = Fb * dbdy
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, x(i),y(i),z(i)
      !print *, F, Fx, Fy, Fz
      !print *, kterm, Fb, flag,k,tan(beta/2)**2
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 10) then
        delt = 0.
        l(i) = 0.
        m(i) = 0.
        n(i) = 0.
      end if
      c = c + 1
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = -Fx/Fp
    uy(i) = -Fy/Fp
    uz(i) = -Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wsprimary

!This function traces to a Wolter-Schwarzschild secondary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wssecondary(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi
  real*8 :: k,kterm,dbdx,dbdy,dbdz,dadb,beta,betas,ff,g,d,a
  real*8 :: gam,dbdzs,dadbs
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,dbdz,flag,c,a,dadbs,dbdzs,gam,dadb)
  do i=1,num
    delt = 100.
    c = 0
    do while(abs(delt)>1.e-8)
      beta = atan2(sqrt(x(i)**2 + y(i)**2),z(i))
      flag = 0
      if (beta<=betas) then
        beta = betas
        kterm = 0
        a = 1/ff
        flag = 1
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
        a = (1-cos(beta))/(1-cos(betas))/ff + &
          (1+cos(beta))/(2*g)*(kterm)**(1+k)
      end if
      F = -z(i) + cos(beta)/a
      !Add corrective term to F if beta was < betas
      if (flag==1) then
        Fb = 0.
        dadbs = sin(betas)/ff/(1-cos(betas)) + &
                (k+1)*(cos(betas)+1)*tan(betas/2)/cos(betas/2)**2/2/g/k
        dbdzs = -sin(betas)**2/sqrt(x(i)**2+y(i)**2)
        gam = (-ff*sin(betas)-ff**2*cos(betas)*dadbs)*dbdzs
        F = F + gam*(z(i)-sqrt(x(i)**2+y(i)**2)/tan(betas))
        Fx = -2./tan(betas)*x(i)/sqrt(x(i)**2+y(i)**2)
        Fy = -2./tan(betas)*y(i)/sqrt(x(i)**2+y(i)**2)
        Fz = gam - 1.
        !print *, x(i), y(i), z(i)
        !print *, F, Fx, Fy, Fz
        !print *, delt
        !read *, dum
      !Otherwise, business as usual
      else
        dadb = sin(beta)/ff/(1-cos(betas)) - &
               sin(beta)/(2*g)*(kterm)**(1+k) + &
               (k+1)*(cos(beta)+1)*tan(beta/2)*kterm**k/2/g/k/(cos(beta/2)**2)
        Fb = -sin(beta)/a - cos(beta)/a**2*dadb
        dbdx = x(i)*z(i)/(x(i)**2+y(i)**2+z(i)**2)/sqrt(x(i)**2+y(i)**2)
        dbdy = y(i)*z(i)/(x(i)**2+y(i)**2+z(i)**2)/sqrt(x(i)**2+y(i)**2)
        dbdz = -sqrt(x(i)**2+y(i)**2)/(x(i)**2+y(i)**2+z(i)**2)
        Fx = Fb * dbdx
        Fy = Fb * dbdy
        Fz = -1. + Fb*dbdz
      end if
      !We have derivatives, now compute the iteration
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, x(i), y(i), z(i)
      !print *, F, Fx, Fy, Fz
      !print *, delt
      !read *, dum
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 10) then
        delt = 0.
        l(i) = 0.
        m(i) = 0.
        n(i) = 0.
      end if
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
      c = c+1
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !print *, F, Fx, Fy, Fz
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wssecondary

!Intersection with SPO Cone
subroutine spoCone(x,y,z,l,m,n,ux,uy,uz,num,R0,tg)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: R0,tg
!  real*8 :: k,kterm,dbdx,dbdy,dbdz,dadb,beta,betas,ff,g,d,a
!  real*8 :: gam,dbdzs,dadbs
!  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb
  integer :: i
  real*8 :: A,B,C,sl,det,t1,t2

  !Loop through rays and trace to mirror
  sl = tan(tg)
  !$omp parallel do private(i,A,B,C,det,t1,t2)
  do i=1,num
    !Solve quadratic equation for ray advancement distance
    A = n(i)**2*sl**2 - m(i)**2 - l(i)**2
    B = 2*n(i)*sl*R0 + 2*z(i)*sl**2*n(i) - 2*x(i)*l(i) - 2*y(i)*m(i)
    C = R0**2 + 2*sl*R0*z(i) + z(i)**2*sl**2 - x(i)**2 - y(i)**2
    det = B**2 - 4*A*C
    if (det .ge. 0) then
      t1 = (-B + sqrt(det))/(2*A)
      t2 = (-B - sqrt(det))/(2*A)
      if (abs(t2) < abs(t1)) then
        t1 = t2
      end if
      !Set new ray position
      x(i) = x(i) + t1*l(i)
      y(i) = y(i) + t1*m(i)
      z(i) = z(i) + t1*n(i)
      !Set up surface normal
      ux(i) = -x(i)/sqrt(x(i)**2+y(i)**2)*cos(tg)
      uy(i) = -y(i)/sqrt(x(i)**2+y(i)**2)*cos(tg)
      uz(i) = sin(tg)!*abs(z(i))/z(i)
    else
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !print *, F, Fx, Fy, Fz
    !read *, dum
  end do
  !$omp end parallel do

end subroutine spoCone


!Radially grooved grating diffraction
!Assumes grating in x y plane, with grooves converging at 
!hubdist in positive y direction
subroutine radgrat(x,y,l,m,n,num,hubdist,dpermm,order,wave)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num)
  real*8, intent(in) :: hubdist,dpermm,wave,order
  integer :: i
  real*8 :: d, yaw, pi, dum, det

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Compute local d spacing in nm
    d = dpermm * sqrt((hubdist-y(i))**2 + x(i)**2)
    !Compute local yaw
    yaw = pi/2 + atan(x(i)/(hubdist-y(i)))
    !print *, x(i),y(i),d,yaw
    !print *, l(i),m(i),n(i)

    !Evanescence?
    det = l(i)**2+m(i)**2
    if (det<1) then
      !Compute new direction cosines
      l(i) = l(i) + sin(yaw)*order*wave/d
      m(i) = m(i) - cos(yaw)*order*wave/d
      n(i) = sqrt(1. - l(i)**2 - m(i)**2)
    end if
  end do

end subroutine radgrat

!Radially grooved grating diffraction with wavelength vector
!Assumes grating in x y plane, with grooves converging at 
!hubdist in positive y direction
subroutine radgratW(x,y,l,m,n,wave,num,hubdist,dpermm,order)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num),wave(num)
  real*8, intent(in) :: hubdist,dpermm,order
  integer :: i
  real*8 :: d, yaw, pi, dum, det

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Compute local d spacing in nm
    d = dpermm * sqrt((hubdist-y(i))**2 + x(i)**2)
    !Compute local yaw
    yaw = pi/2 + atan(x(i)/(hubdist-y(i)))
    !print *, x(i),y(i),d,yaw
    !print *, l(i),m(i),n(i)

    !Evanescence?
    det = l(i)**2+m(i)**2
    if (det<1) then
      !Compute new direction cosines
      l(i) = l(i) + sin(yaw)*order*wave(num)/d
      m(i) = m(i) - cos(yaw)*order*wave(num)/d
      n(i) = sqrt(1. - l(i)**2 - m(i)**2)
    end if
  end do

end subroutine radgratW

!Linear grating with groove period d
!Wavelength wave
!Groove direction assumed in y direction
subroutine grat(x,y,l,m,n,num,d,order,wave)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num)
  real*8, intent(in) :: d,wave,order
  integer :: i
  real*8 :: pi, dum

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Compute new direction cosines
    l(i) = l(i) - order*wave/d
    n(i) = sqrt(1 - l(i)**2 - m(i)**2)
    !Evanescence?
    if ((l(i)**2+m(i)**2)>1) then
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    end if
  end do

end subroutine grat

!Function to trace to a general conic
!Vertex assumed at origin, opening up in the +z direction
!Radius of curvature and conic constant are required parameters
!Traces to solution closest to vertex
subroutine conic(x,y,z,l,m,n,ux,uy,uz,num,R,K)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: R,K
  real*8 :: s1,s2,s,z1,z2,denom,b,c,disc
  integer :: i

  !Loop through rays and trace to the conic
  do i=1,num
    !Compute amount to move ray s
    s = 0.
    if (K .eq. -1 .and. abs(n(i))==1.) then
      s = (x(i)**2 + y(i)**2 - 2*R*z(i)) / (2*R*n(i))
    else
      denom = l(i)**2 + m(i)**2 + (K+1)*n(i)**2
      b = x(i)*l(i) + y(i)*m(i) + ((K+1)*z(i) - R)*n(i)
      b = b/denom
      c = x(i)**2 + y(i)**2 + (K+1)*z(i)**2 - 2*R*z(i)
      c = c/denom
      disc = b**2 - c
      if (disc .ge. 0.) then
        s1 = -b + sqrt(disc)
        s2 = -b - sqrt(disc)
        !Choose smallest positive resulting Z
        z1 = z(i) + s1*n(i)
        z2 = z(i) + s2*n(i)
        if (abs(z1) .le. abs(z2)) then
          s = s1
        else
          s = s2
        end if
      end if
    end if
    !Advance ray
    if (s==0.) then
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      x(i) = x(i) + l(i)*s
      y(i) = y(i) + m(i)*s
      z(i) = z(i) + n(i)*s
      !Compute normal derivative
      denom = sqrt(R**2 - K*(x(i)**2+y(i)**2))
      ux(i) = -x(i) / denom
      uy(i) = -y(i) / denom
      uz(i) = -R/abs(R) * sqrt(R**2 - (K+1)*(x(i)**2+y(i)**2))
      uz(i) = -uz(i) / denom
    end if
  end do

end subroutine conic

!Function to trace to a general conic
!Vertex assumed at origin, opening up in the +z direction
!Radius of curvature and conic constant are required parameters
!Traces to solution closest to vertex
subroutine conicopd(opd,x,y,z,l,m,n,ux,uy,uz,num,R,K,nr)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num),opd(num)
  real*8, intent(in) :: R,K,nr
  real*8 :: s1,s2,s,z1,z2,denom,b,c,disc
  integer :: i

  !Loop through rays and trace to the conic
  do i=1,num
    !Compute amount to move ray s
    s = 0.
    if (K .eq. -1 .and. abs(n(i))==1.) then
      s = (x(i)**2 + y(i)**2 - 2*R*z(i)) / (2*R*n(i))
    else
      denom = l(i)**2 + m(i)**2 + (K+1)*n(i)**2
      b = x(i)*l(i) + y(i)*m(i) + ((K+1)*z(i) - R)*n(i)
      b = b/denom
      c = x(i)**2 + y(i)**2 + (K+1)*z(i)**2 - 2*R*z(i)
      c = c/denom
      disc = b**2 - c
      if (disc .ge. 0.) then
        s1 = -b + sqrt(disc)
        s2 = -b - sqrt(disc)
        !Choose smallest positive resulting Z
        z1 = z(i) + s1*n(i)
        z2 = z(i) + s2*n(i)
        if (abs(z1) .le. abs(z2)) then
          s = s1
        else
          s = s2
        end if
      end if
    end if
    !Advance ray
    if (s==0.) then
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    else
      x(i) = x(i) + l(i)*s
      y(i) = y(i) + m(i)*s
      z(i) = z(i) + n(i)*s
      opd(i) = opd(i) + s*nr
      !Compute normal derivative
      denom = sqrt(R**2 - K*(x(i)**2+y(i)**2))
      ux(i) = -x(i) / denom
      uy(i) = -y(i) / denom
      uz(i) = -R/abs(R) * sqrt(R**2 - (K+1)*(x(i)**2+y(i)**2))
      uz(i) = -uz(i) / denom
    end if
  end do

end subroutine conicopd
