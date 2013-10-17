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


! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
  dx = ( x_max - x_min ) / real ( xn - 1, kind = 8 ) 
  x = linspace ( xn, x_min, x_max )
  dy = ( y_max - y_min ) / real ( yn - 1, kind = 8 ) 
  y = linspace ( yn, y_min, y_max )
  dt = ( t_max - t_min ) / real ( tn - 1, kind = 8 ) 
  t = linspace ( tn, t_min, t_max)
  
  
! BOUNDARY CONDITIONS
  ic0(:,:) = 300.0 ! IC

  bcx0(:,1) = 273.0 ! bottom
  bcx0(:,2) = 273.0 ! top

  bcy0(1,:) = 273.0 ! left
  bcy0(2,:) = 273.0 ! right


  bcyPsi(1,1:yn) = 0.0 ! left
  bcyPsi(2,1:yn) = 0.0 ! right
  bcxPsi(1:xn,1) = 0.0 ! bottom
  bcxPsi(1:xn,2) = 0.0 ! top

! PERMEABILITY
do i=1,yn

	if (y(i) .ge. y_min) then
	permeability(:,i) = (0.5+0.5*tanh((y(i)+((800.0)))/20.0))*1e-13 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+((800.0)))/20.0)))*1e-21
	end if
	if (y(i) .gt. -500.0) then
	permeability(:,i) = (0.5+0.5*tanh((y(i)+((200.0)))/20.0))*1e-21 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+((200.0)))/20.0)))*1e-13
	end if
	
end do

! HEAT TRANSFER PARAMETERS
kMat = 2.6/(1000.0*4186.0)
ki=2.6/(1000.0*4186.0)


return

end subroutine init



end module initialize