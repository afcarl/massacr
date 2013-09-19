module initialize

use globals

save

real(8) :: x(xn), y(yn), t(tn), ytemp(xn,yn)
real(8) :: permeability(xn,yn), permx(xn,yn), permy(xn,yn), permLong((xn-2)*(yn-2))
real(8) :: rhoLong((xn-2)*(yn-2))
real(8) :: permxLong((xn-2)*(yn-2)), permyLong((xn-2)*(yn-2))
real(8) :: rho0(xn,yn)
real(8) :: bcx0(xn,2), bcy0(2,yn), bcxPsi(xn,2), bcyPsi(2,yn), ic0(xn,yn)
real(8) :: kMat(xn,yn), porosity(xn,yn), permeable(xn)

contains

! ----------------------------------------------------------------------------------%%
!
! SUBROUTINE TO INITIALIZE, JUST CALL IT
!
! ----------------------------------------------------------------------------------%%
  
subroutine init ()
use globals
integer :: m,n
real(8) :: waveIntercept




! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
  dx = ( x_max - x_min ) / real ( xn - 1, kind = 8 ) 
  x = linspace ( xn, x_min, x_max )
  dy = ( y_max - y_min ) / real ( yn - 1, kind = 8 ) 
  y = linspace ( yn, y_min, y_max )
  dt = ( t_max - t_min ) / real ( tn - 1, kind = 8 ) 
  t = linspace ( tn, t_min, t_max)
  
  
! BOUNDARY CONDITIONS
  ic0(:,:) = 0.00 ! IC
  do i=1,xn
  	!ic0(i,:) = linspace(yn, 340.0D+00, 200.0D+00)
  end do
  bcx0(:,1) = 273.0 ! bottom
  bcx0(:,2) = 273.0 ! top
  !BENCHMARK
!  bcx0(:,1) = 0.0 ! bottom
  bcx0(:,2) = 0.0 ! top
  do i =1,xn
  	!bcx0(i,1) = 295.0 + real(i,kind=4)/20.0 !sqrt(real(i,kind=4))
  	!bcx0(i,1) = 100.0 + 10.0*cos(2.0*3.14*x(i)/3000.0)
  	!bcx0(i,2) = 200.0 + 10.0*cos(2.0*3.14*x(i)/1300.0)
  end do
  
  bcy0(1,:) = 273.0 ! left
  bcy0(2,:) = 273.0 ! right
  !BENCHMARK
  bcy0(1,:) = 0.5 ! left
  bcy0(2,:) = -0.5 ! right
  
!  bcx0(1:yn,1) = linspace(yn, bcy0(1,1), bcy0(1,2)) ! left
!  bcx0(1:yn,2) = linspace(yn, bcy0(1,1), bcy0(1,2)) ! right
  !bcy0(1,1:yn) = 275.0 ! left
  !bcy0(2,1:yn) = 275.0 ! right
  
  
  bcyPsi(1,1:yn) = 0.0 ! left
  bcyPsi(2,1:yn) = 0.0 ! right
  bcxPsi(1:xn,1) = 0.0 ! bottom
  bcxPsi(1:xn,2) = 0.0 ! top
  do i=1,xn
  !bcxPsi(i,2) = +(dy/2)*xn*2.0*.078/(3.14e7)-(dy)*i*2.0*.078/(3.14e7)
 end do


! SET UP MATERIAL
	 permeability(:,:) = 1e-16
     kMat = 2.6/(1000.0*4186.0)
     porosity = 0.01
	do m = 1,xn
	waveIntercept = (y_min/25.0) * cos(5.0*3.14*x(m)/x_max)
	waveIntercept = 0.0
	do n = 1,yn
	
	    if (y(n) .gt. -900 + waveIntercept) then
		 !permeability(m,n) = 5e-17
		 kMat(m,n) = 2.6/(1000.0*4186.0)
		end if
		if (y(n) .gt. -850 + waveIntercept) then
		 !permeability(m,n) = 5e-16
		 kMat(m,n) = 2.6/(1000.0*4186.0)
		end if
		if (y(n) .gt. -800 + waveIntercept) then
		 !permeability(m,n) = 5e-15
		 kMat(m,n) = 2.6/(1000.0*4186.0)
		end if
	
	
		if (y(n) .gt. -750 + waveIntercept) then
		 permeability(m,n) = 1e-13
		 kMat(m,n) = 2.6/(1000.0*4186.0)
		end if
		
		
		if (y(n) .gt. -250 + waveIntercept) then
		 permeability(m,n) = 1e-15
		 kMat(m,n) = 1.6/(1000.0*4186.0)
		end if
		
		!!! FISHER!!!!
		
		if (y(n) .gt. -950 + waveIntercept) then
		 !permeability(m,n) = 1e-17
		!kMat(m,n) = 1.85
		porosity(m,n) = 0.05
		end if
		
		if (y(n) .gt. -750 + waveIntercept) then
		 !permeability(m,n) = 1e-17
		! kMat(m,n) = 1.75
		 porosity(m,n) = 0.08
		end if
		
		if (y(n) .gt. -450 + waveIntercept) then
		 !permeability(m,n) = 5e-15
		! kMat(m,n) = 1.6
		 porosity(m,n) = 0.1
		end if
		
		if (y(n) .gt. -350 + waveIntercept) then
		! 13 -> 15
		 !permeability(m,n) = 1e-13
		! kMat(m,n) = 1.5
		 porosity(m,n) = 0.125
		end if
		
		if (y(n) .gt. -250 + waveIntercept) then
		 !permeability(m,n) = 5.94e-17
		! kMat(m,n) = 1.11
		 porosity(m,n) = 0.69
		end if
		
		if (y(n) .gt. -150 + waveIntercept) then
		 !permeability(m,n) = 8.88e-16
		! kMat(m,n) = 0.89
		 porosity(m,n) = 0.8
		end if
		
		if (y(n) .gt. -50 + waveIntercept) then
		 !permeability(m,n) = 4.88e-15
		! kMat(m,n) = 0.86
		 porosity(m,n) = 0.85
		end if


	end do
	end do
	
!permeability = 2e-17

!waveIntercept = .02 * cos(2.0*3.14*x(m)/(5.0*x_max))
waveIntercept = 0.0
write(*,*) waveIntercept
!waveIntercept = 0.0
do i=1,yn
	!permeability(:,i) = 1e-18
	!permeability(i,:) = moving_average(permeability(i,:),yn, 37)
	
	
!	if (y(i) .ge. y_min) then
!	permeability(m,i) = (0.5+0.5*tanh((y(i)+((700.0/1300.0)+waveIntercept))/.06))*5e-14 &
!	&+ (1.0 - (0.5+0.5*tanh((y(i)+((700.0/1300.0)+waveIntercept))/.06)))*1e-18
!	end if
!	if (y(i) .gt. -.3) then
!	permeability(m,i) = (0.5+0.5*tanh((y(i)+((200.0/1300.0)+waveIntercept))/.06))*5e-18 &
!	&+ (1.0 - (0.5+0.5*tanh((y(i)+((200.0/1300.0)+waveIntercept))/.06)))*5e-14
!	end if
	
	if (y(i) .ge. y_min) then
	permeability(:,i) = (0.5+0.5*tanh((y(i)+((800.0)))/70.0))*1e-13 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+((800.0)))/70.0)))*1e-21
	end if
	if (y(i) .gt. -500.0) then
	permeability(:,i) = (0.5+0.5*tanh((y(i)+((200.0)))/70.0))*1e-21 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+((200.0)))/70.0)))*1e-13
	end if
	
	!permeability = 1e-18
	
end do

 permeability =  permeability
 kMat = 2.6/(1000.0*4186.0)
ki=2.6/(1000.0*4186.0)


return

end subroutine init












end module initialize