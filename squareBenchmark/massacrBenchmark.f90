
! ----------------------------------------------------------------------------------%%
!
! MAIN MASSACR METHOD
! 
! SUMMARY: Solve governing equations (conservation of momentum and thermal energy) 
!          in 2D and writes data to text (.txt) and netCDF (.nc) files. 
! PARAMETERS: dx,dy : x,y grid spacing
!             dt : timestep
!             t : computation times
!             x,y : gridpoints
!             bcx0, bcy0 : boundary conditions
!             ic0 : initial conditions
!             u,v : vector field components
!             
! RETURNS:
!
! ----------------------------------------------------------------------------------%%

PROGRAM main
!use netcdf

use globals
use initialize

implicit none

interface



  function h_next (h, psi, rho_in, flux)
  use globals
  use initialize
! integers
  integer :: i, j, n, ii, c1, c2, c3, c4, m=5
! inputs 
  real(8) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2)
! velocity stuff
  real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
  real(8) ::  velocities0(xn,2*yn)
! matrix stuff
  real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
!  real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
!  real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
  real(8) :: kMatLong((xn-2)*(yn-2))
  end function h_next
  
  function psi_next (h, rhs0, psi,top_in,rho_in)
  use globals
  use initialize
  ! integers
  integer :: i, j, n, m=5
  ! inputs
  real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong((xn-2)*(yn-2))
  real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn)
  ! matrix stuff
  real(8) :: uVec((xn-2)*(yn-2)), psiLong((xn)*(yn)), psi_nextRow((xn-2)*(yn-2))
  real(8) :: psi_next(xn,yn), dsbit, fbit, conv
  real(8) :: way1((xn-2)*(yn-2)), way2((xn-2)*(yn-2))
  real(8) :: yep, top_in(xn,1)
  ! tridiag
  real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  end function psi_next
  
  function rho_next (h_in)
  use globals
  use initialize
  integer :: i,j
  real(8) :: h_in(xn,yn), rho_next(xn,yn)
  end function rho_next

  function velocities(psi)
  use globals
  use initialize
  implicit none
  real(8) :: velocities(xn,2*yn), psi(xn,yn)
  real(8) :: u0(xn,yn), v0(xn,yn)
  end function velocities


function partial(array,rows,cols,d1,d2,dim)
use globals
use initialize
implicit none
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)
end function partial

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix

end function write_matrix

function write_vec ( n, vector, filename )
use globals
  implicit none
  integer :: n, j, output_status, unit0
  character ( len = * ) filename 
  real(4)  :: vector(n), write_vec

end function write_vec

function v_next(dPdy_in, rho_in)
use globals
use initialize
  implicit none
  
  integer :: i,j
  real(8) :: dPdy_in(xn,yn), rho_in(xn,yn), v_next(xn,yn)
end function v_next

   
end interface



  real(8) :: cfl, test(10,11)
  real(8) :: h(xn,yn), psi(xn,yn) ! xn ROWS DEEP & yn COLUMNS WIDE 
!  real(8) :: hmat(xn,(yn*tn)), psimat(xn,(yn*tn)), umat(xn,(yn*tn)), vmat(xn,(yn*tn))
  real(8) :: hmat(xn,(yn)), psimat(xn,(yn)), umat(xn,(yn)), vmat(xn,(yn))
  real(8) :: rhs0(xn,yn), velocities0(xn,2*yn)
  real(8) :: u(xn,yn), uLong(xn*yn), v(xn,yn), vLong(xn*yn)
  real(8), allocatable :: linsp(:)
  
  real(8) :: pressure0(xn,yn), pressure(xn,yn)
  real(8) :: rho(xn,yn), flux(xn,2)



  
  integer :: unit
  real(8) :: yep
  
  integer :: xInt, yInt, tInt, hInt, uInt, vInt
  integer :: ncid
  integer :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
  integer :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid


integer :: i, j, ii, m, n

real(8) :: nusseltLocalv(xn,1), nuBar



! INITIALIZE
call init()


  
! PUT BOUNDARY CONDITIONS IN
psi(1,1:yn) = bcyPsi(1,1:yn)
psi(xn,1:yn) = bcyPsi(2,1:yn)
psi(1:xn,1) = bcxPsi(1:xn,1)
psi(1:xn,yn) = bcxPsi(1:xn,2)
psi=0.0


h(1:xn,1:yn) = ic0
h(1,1:yn) = bcy0(1,1:yn)
h(xn,1:yn) = bcy0(2,1:yn)
h(1:xn,1) = bcx0(1:xn,1)
h(1:xn,yn) = bcx0(1:xn,2)

! PUT INITIAL VALUES INTO MATRIX FOR FILE
hmat(1:xn,1:yn) = h
psimat(1:xn,1:yn) = psi
umat(1:xn,1:yn) = u
vmat(1:xn,1:yn) = v



!CONDUCTION ONLY SIMULATION
! RIGHT NOW IT DOES NOTHING AND CONDUCTION FUNCTION IS OUT OF DATE
do j = 2, 2
write(*,*) "conduction"
write(*,*) j

!h = conduction(h)
psi = 0.0
u = 0.0
v = 0.0

! ADD EACH TIMESTEP TO MATRICES
!hmat(1:xn,1+yn*(j-1):1+yn*(j)) = h
!psimat(1:xn,1+yn*(j-1):1+yn*(j)) = psi
!umat(1:xn,1+yn*(j-1):1+yn*(j)) = u
!vmat(1:xn,1+yn*(j-1):1+yn*(j)) = v
hmat(1:xn,1:yn) = h
psimat(1:xn,1:yn) = psi
umat(1:xn,1:yn) = u
vmat(1:xn,1:yn) = v

end do






do j = 2, tn

write(*,*) j

! HEAT FLUX BOUNDARY CONDITIONS

  ! bottom
  do i = 1,xn
  !flux(i,1) = -h(i,3)/3.0 +4.0*h(i,2)/3.0 +(.27+.01*cos(2.0*x(i)*3.14/x_max))*2.0*dx/(lambda*3.0)
  flux(i,1) = h(i,2) +(.27)*dy/(lambda)
  end do
  
  !write(*,*) maxval(flux(:,1))
  ! top
  do i = 1,xn
  !flux(i,2) = -h(i,xn-2)/3.0 +4.0*h(i,xn-1)/3.0 -.27*2.0*dy/(lambda*3.0)
  !flux(i,2) = h(i,xn-2) -(.27)*dy/(lambda)
  !flux(i,2) = (4.0/3.0)*h(i,xn-1) - (1.0/3.0)*h(i,xn-2)
  end do
  
! ATTEMPT AT MIXED BOUNDARY CONDITION

  !do i=1,xn
  !	if (v(i,yn) .gt. 0) then
  !		h(i,yn) = (4.0/3.0)*h(i,yn-1) - (1.0/3.0)*h(i,yn-2)
  !	end if
  !	if (v(i,yn) .le. 0) then
  !		h(i,yn) = 275.0
  !	end if
  !end do

  ! PUT IN BOUNDARY CONDITIONS BETWEEN STEPS
  
  rho = rho_next(h)
  h = h_next(h, psi,rho, flux)
  
!!!!!!!!!!!! THIS !!!!!!!!!!!
!  h(1,:) = (4.0/3.0)*h(2,:) - (1.0/3.0)*h(3,:) ! left
!  h(xn,:) = (4.0/3.0)*h(xn-1,:) - (1.0/3.0)*h(xn-2,:) ! right
!  h(:,1) = flux(:,1)
!  h(:,yn) = flux(:,2)
  

  
! BENCHMARK FIXED/ADIABATIC BOUNDARY CONDITIONS
  

 flux(:,1) = (4.0/3.0)*h(:,2) - (1.0/3.0)*h(:,3) ! bottom
 flux(:,2) = (4.0/3.0)*h(:,xn-1) - (1.0/3.0)*h(:,xn-2) ! top
  h(:,1) = flux(:,1)
  h(:,yn) = flux(:,2)
  h(1,2:xn-1) = 0.5
 h(xn,2:xn-1) = -0.5
  

  

! BENCHMARK 
rhs0 = - 1000.0*partial(h,xn,yn,dx,dy,1) 
! OLD METHOD
!rhs0 = - (permeability/(viscosity))*rho_fluid*g*alpha*partial(h,xn,yn,dx,dy,1)
!!!!!!rhs0 = - (1.0/(viscosity))*g*rho_fluid*alpha*partial(h,xn,yn,dx,dy,1)
! LOTS OF COEFFICIENTS FOR SOME REASON
!rhs0 = - cp*g*alpha*200.0*1300.0*rho_fluid/(viscosity*lambda)*partial(h,xn,yn,dx,dy,1)

psi = psi_next(h, rhs0, psi, permeable, rho)

! PUT IN BOUNDARY CONDITIONS BETWEEN STEPS
psi(1,1:yn) = bcyPsi(1,1:yn) ! left
!psi(1,:) = ((4.0/3.0)*psi(2,:) - (1.0/3.0)*psi(3,:)) ! left
psi(xn,1:yn) = bcyPsi(2,1:yn) ! right
!psi(xn,:)  = ((4.0/3.0)*psi(xn-1,:) - (1.0/3.0)*psi(xn-2,:)) ! right
psi(1:xn,1) = bcxPsi(1:xn,1) ! bottom
!psi(:,1) = ((4.0/3.0)*psi(:,2) - (1.0/3.0)*psi(:,3)) !bottom
psi(1:xn,yn) = bcxPsi(1:xn,2) ! top
!psi(:,yn) = ((4.0/3.0)*psi(:,yn-1) - (1.0/3.0)*psi(:,yn-2)) 
permeable = psi(:,yn)


! BENCHMARK
!psi(:,yn) = 0.0
!psi(:,1) = 0.0
!psi(1,:) = 0.0
!psi(xn,:) = 0.0


! VELOCITIES
velocities0 = velocities(psi)
u = velocities0(1:xn,1:yn)/rho
v = velocities0(1:xn,yn+1:2*yn)/rho
   
! ADD EACH TIMESTEP TO MATRICES
!	 hmat(1:xn,1+yn*(j-1):1+yn*(j)) = h
!	 psimat(1:xn,1+yn*(j-1):1+yn*(j)) = psi
!    umat(1:xn,1+yn*(j-1):1+yn*(j)) = u
!    vmat(1:xn,1+yn*(j-1):1+yn*(j)) = v
hmat(1:xn,1:yn) = h
psimat(1:xn,1:yn) = psi
umat(1:xn,1:yn) = u
vmat(1:xn,1:yn) = v

! NUSSELT NUMBER STUFF

nuBar = 0.0
do i=1,yn
	nusseltLocalv = -partial(h(:,i),xn,1,dx,dy,1)
	nuBar = nuBar + nusseltLocalv(1,1)*dy
end do
write(*,*) "nu Bar"
write(*,*) nuBar

end do




! WRITE EVERYTHING TO TEXT FILES
yep = write_vec ( xn, real(x,kind=4), 'x.txt' )
yep = write_vec ( yn, real(y,kind=4), 'y.txt' )
yep = write_vec ( tn, real(t, kind=4), 't.txt' )
!yep = write_matrix ( xn, yn*tn, real(hmat, kind = 4), 'hT.txt' )
!yep = write_matrix ( xn, yn*tn, real(psimat,kind=4), 'psiMat.txt' )
!yep = write_matrix ( xn, yn*tn, real(umat,kind=4), 'uMat.txt' )
!yep = write_matrix ( xn, yn*tn, real(vmat,kind=4), 'vMat.txt' )
yep = write_matrix ( xn, yn, real(hmat, kind = 4), 'h.txt' )
yep = write_matrix ( xn, yn, real(psimat,kind=4), 'psiMat.txt' )
yep = write_matrix ( xn, yn, real(umat,kind=4), 'uMat.txt' )
yep = write_matrix ( xn, yn, real(vmat,kind=4), 'vMat.txt' )
yep = write_matrix ( xn, yn, real(rho,kind=4), 'rho.txt' )
yep = write_matrix ( xn, yn,real(permeability,kind=4), 'permeability.txt' )

! WRITE THINGS TO NETCDF FILES
!  call check( nf90_create('thermalNRG.nc', NF90_CLOBBER, ncid) )
!  call check( nf90_def_dim(ncid, "x", xn, x_dimid) )
!  call check( nf90_def_dim(ncid, "y", yn, y_dimid) )

!  call check(nf90_def_var(ncid, "h", NF90_FLOAT, (/ x_dimid, y_dimid /), h_varid) )
!  call check(nf90_def_var(ncid, "u", NF90_FLOAT, (/ x_dimid, y_dimid /), u_varid) )
!  call check(nf90_def_var(ncid, "v", NF90_FLOAT, (/ x_dimid, y_dimid /), v_varid) )
!  call check(nf90_enddef(ncid))

!call check( nf90_put_var(ncid, h_varid, h) )
!call check( nf90_put_var(ncid, u_varid, u) )
!call check( nf90_put_var(ncid, v_varid, v) )

!  
!  CALL check(nf90_put_att(ncid, x_varid, "units", "meter")) 
!  CALL check(nf90_put_att(ncid, y_varid, "units", "meter")) 
!  CALL check(nf90_put_att(ncid, h_varid, "units", "K"))
!  CALL check(nf90_put_att(ncid, u_varid, "units", "m/s"))
!  CALL check(nf90_put_att(ncid, v_varid, "units", "m/s"))
!  
!  CALL check(nf90_put_att(ncid, x_varid, "axis", "X")) 
!  CALL check(nf90_put_att(ncid, y_varid, "axis", "Y")) 


  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
!  call check( nf90_close(ncid) )

write(*,*) " "
write(*,*) "ALL DONE!"


write(*,*) mod(0,10)
END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: Computers the thermal profile of the next timestep
!
! INPUTS: xn,yn : number of x,y gridpoints
!         x,y : values of x,y gridpoints
!         bcx0,bcy0 : boundary conditions
!         h : thermal profile of the previous timestep
!         rhs0 : source/sink term
!         u,v : x,y velocity field
!         dx,dy : x,y grid spacing
!         k : thermal diffusivity
!
! FUNCTIONS USED: Gsselm : gaussian elimination
!                 tridiag : solves tridiagonal system
!
! RETURNS: h_next(xn,yn)
!
! ----------------------------------------------------------------------------------%%

function h_next (h, psi, rho_in, flux)
use globals
use initialize
implicit none

interface
  

  function velocities(psi)
  use globals
  implicit none
  integer :: i,ii,j,jj
  real(8) :: uu(xn+1,yn+1), vv(xn+1,yn+1), velocities(xn,2*yn), psi(xn,yn)
  real(8) :: u0(xn,yn), v0(xn,yn), dLeft, dRight, dTop, dBottom
  end function
  

  
end interface

  ! integers
  integer :: i, j, n, ii, m=3
  ! inputs 
  real(8) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2)
  ! velocity stuff
  real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
  real(8) ::  velocities0(xn,2*yn)
  ! matrix stuff
  real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
!  real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
!  real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
  real(8) :: kMatLong((xn-2)*(yn-2))
  real(8) :: mn(xn,yn)
  
  
mn = h


velocities0 = velocities(psi)
u = velocities0(1:xn,1:yn)
v = velocities0(1:xn,yn+1:2*yn)
uLong = reshape(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
vLong = reshape(transpose(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! PRINT VELOCITIES AT EACH TIMESTEP TO GIVE SCOPE
write(*,*) " "
write(*,*) "max u"
write(*,*) maxval(abs(u))
write(*,*) "max v"
write(*,*) maxval(abs(v))

!write(*,*) (2.0*dt)/(dy*dy)
write(*,*) " "
write(*,*) "velocity check"
write(*,"(F10.5)") (dt*maxval(abs(u)))/(dx)
write(*,"(F10.5)") (dt*maxval(abs(v)))/(dy)
write(*,*) "conduction check"
!write(*,"(F10.5)") (2.0*dt*lambda)/(4186.0*dy*dy)
write(*,*) (2.0*dt)/(dy*dy)
write(*,*) " "



h0 = h

  qx = dt/(dx)
  qy = dt/(dy)
  sx = (2.0*dt*lambda)/(4186.0*dx*dx)
  sy = (2.0*dt*lambda)/(4186.0*dy*dy)
  ! BENCHMARK
  sx = (2.0*dt)/(dx*dx)
  sy = (2.0*dt)/(dy*dy)

! ACCOUNT FOR BOUNDARY CONDITIONS IN THE MATRIX
 h(2,2:xn-1) = h(2,2:xn-1) + h0(1,2:xn-1)*sx/2.0  ! left
 h(yn-1,2:xn-1) = h(yn-1,2:xn-1) + h0(xn,2:xn-1)*sx/2.0  ! right
 !h(2:xn-1,2) = h(2:xn-1,2) + h0(2:xn-1,1)*sy/2.0
 !h(2:xn-1,xn-1) = h(2:xn-1,xn-1) + h0(2:xn-1,xn)*sy/2.0
 

uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))


!  aa = 0.0

  
!  do i = 1,(xn-2)*(yn-2)

  ! 4th order central difference (last edge)
!  aa(i,i) = 1.0+sx
!    if (i .gt. 1) then
!  	aa(i,i-1) = -sx/2.0 - uLong(i)*qx/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	aa(i,i+1) = -sx/2.0 + uLong(i)*qx/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	aa(i,i+2) = 0.0
!  	end if
!  	if (i .gt. 2) then
!  	aa(i,i-2) = 0.0
!  	end if
  	
  	
  	! edges

  ! fwd diff at first edge no matter what
!  if (any(mod((/i-1, i-2, i-3/),xn-2) .eq. 0.0)) then
!  aa(i,i) = 1.0+ sx/2.0 - 3.0*uLong(i)*qx/2.0
!   if (i .gt. 1) then
!  	aa(i,i-1) =  0.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	aa(i,i+1) = -2.0*sx/2.0 + 4.0*uLong(i)*qx/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	aa(i,i+2) = sx/2.0 - 1.0*uLong(i)*qx/2.0
!  	end if
!  	if (i .gt. 2) then
!  	aa(i,i-2) = 0.0
!  	end if
!  end if
  
  ! bckwrd diff at last edge no matter what
!  if (any(mod((/i, i+1, i+2/),xn-2) .eq. 0.0)) then
!  aa(i,i) = 1.0 + sx/2.0 + 3.0*uLong(i)*qx/2.0
!   if (i .gt. 1) then
!  	aa(i,i-1) = -2.0*sx/2.0 - 4.0*uLong(i)*qx/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	aa(i,i+1) = 0.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	aa(i,i+2) = 0.0
!  	end if
!  	if (i .gt. 2) then
!  	aa(i,i-2) = sx/2.0 + 1.0*uLong(i)*qx/2.0
!  	end if
!  end if

!  end do
  
  
!  do i = 1,((xn-2)-1)
!    ii = i*(xn-2)
!    aa(ii+1,ii) = 0.0
!  	aa(ii,ii+1) = 0.0
!  	aa(ii+2,ii) = 0.0
!  	aa(ii,ii+2) = 0.0
!  	aa(ii+1,ii-1) = 0.0
!  	aa(ii-1,ii+1) = 0.0
!  end do



  
  ! MAKE THE BAND
  aBand = 0.0
  do i = 1,(xn-2)*(yn-2)
  aBand(i,2) = 1.0+sx
  if (i-1 .gt. 0) then
  	aBand(i,1) = -sx/2.0 + uLong(i)*qx/2.0
  end if
  if (i+1 .le. (xn-2)*(yn-2)) then
  	aBand(i,3) = -sx/2.0 - uLong(i)*qx/2.0
  end if
  
  end do
  
  do i=1,((xn-2)-1)
  ii = i*(xn-2)
  aBand(ii,3) = 0.0
  aBand(ii+1,1) = 0.0
  end do
  
  ! SOLVING EQUATION
!  a(1:(xn-2)*(yn-2),1:(xn-2)*(yn-2)) = aa
!  a(:,(xn-2)*(yn-2)+1) = uVec
  
  !!!!!!!!!!!! THIS !!!!!!!!!!!
  h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))

!  aBand = band(aBand,m,(xn-2)*(yn-2))
!  h_nextRow = solve(aBand,uVec,m,(xn-2)*(yn-2))


  h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))
!  h(:,1) =  bcx0(:,1) ! bottom
!  h(:,yn) = bcx0(:,2) ! top

!h(:,1) = (4.0/3.0)*h(:,2) - (1.0/3.0)*h(:,3) ! bottom
!h(:,yn) = (4.0/3.0)*h(:,xn-1) - (1.0/3.0)*h(:,xn-2) ! top


! ACCOUNT FOR BOUNDARY CONDITIONS IN THE MATRIX
! h(2,2:yn-1) = h(2,2:yn-1) + h0(1,2:yn-1)*sx/2.0
! h(yn-1,2:yn-1) = h(yn-1,2:yn-1) + h0(yn,2:yn-1)*sx/2.0 
 h(:,2) = h(:,2) + h0(:,1)*sy/2.0 !- h(2:xn-1,1)*qx*u(2:xn-1,1) ! bottom
 h(:,xn-1) = h(:,xn-1) + h0(:,xn)*sy/2.0 !+ h(2:xn-1,xn)*qx*u(2:xn-1,xn) ! top



  h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

 
!  bb = 0.0
! do i = 1,(xn-2)*(yn-2)
  
    
! CENTRAL DIFFERENCE (first edge)
!-------------------------------------------------!
!  bb(i,i) = 1.0+sy
!   if (i .gt. 1) then
!  	bb(i,i-1) = -sy/2.0 - vLong(i)*qy/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	bb(i,i+1) = -sy/2.0 + vLong(i)*qy/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	bb(i,i+2) = 0.0
!  	end if
!  	if (i .gt. 2) then
!  	bb(i,i-2) = 0.0
!  	end if
!-------------------------------------------------!

  	! edges

  ! fwd diff at first edge no matter what
!  if (any(mod((/i-1, i-2, i-3/),xn-2) .eq. 0.0)) then
!  bb(i,i) = 1.0+ sy/2.0 - 3.0*vLong(i)*qy/2.0
!   if (i .gt. 1) then
!  	bb(i,i-1) =  0.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	bb(i,i+1) = -2.0*sy/2.0 + 4.0*vLong(i)*qy/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	bb(i,i+2) = sy/2.0 - 1.0*vLong(i)*qy/2.0
!  	end if
!  	if (i .gt. 2) then
!  	bb(i,i-2) = 0.0
!  	end if
!  end if
  
  ! bckwrd diff at last edge no matter what
!  if (any(mod((/i, i+1, i+2/),xn-2) .eq. 0.0)) then
!  bb(i,i) = 1.0 + sy/2.0 + 3.0*vLong(i)*qy/2.0
!   if (i .gt. 1) then
!  	bb(i,i-1) = -2.0*sy/2.0 - 4.0*vLong(i)*qy/2.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)) then
!  	bb(i,i+1) = 0.0
!  	end if
!  	if (i .lt. (xn-2)*(yn-2)-1) then
!  	bb(i,i+2) = 0.0
!  	end if
!  	if (i .gt. 2) then
!  	bb(i,i-2) = sy/2.0 + 1.0*vLong(i)*qy/2.0
!  	end if
!  end if

!end do
  
!  do i = 1,((xn-2)-1)
!    ii = i*(xn-2)
!    bb(ii+1,ii) = 0.0
!  	bb(ii,ii+1) = 0.0
!  	bb(ii+2,ii) = 0.0
!  	bb(ii,ii+2) = 0.0
!  	bb(ii+1,ii-1) = 0.0
!  	bb(ii-1,ii+1) = 0.0
!  end do
  
 


    ! MAKE THE BAND
  bBand = 0.0
  do i = 1,(xn-2)*(yn-2)
  bBand(i,2) = 1.0+sy
  if (i-1 .gt. 0) then
  	bBand(i,1) = -sy/2.0 + vLong(i)*qy/2.0
  end if
  if (i+1 .le. (xn-2)*(yn-2)) then
  	bBand(i,3) = -sy/2.0 - vLong(i)*qy/2.0
  end if
  
  end do
  
  do i=1,((xn-2)-1)
  ii = i*(xn-2)
  bBand(ii,3) = 0.0
  bBand(ii+1,1) = 0.0
  end do

  ! SOLVING EQUATION
!  b(1:(xn-2)*(yn-2),1:(xn-2)*(yn-2)) = bb
!  b(:,(xn-2)*(yn-2)+1) = h_nextRow
  !h_nextRow = tridiag(b,(xn-2)*(yn-2))
  
  h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))

!  bBand = band(bBand,m,(xn-2)*(yn-2))
!  h_nextRow = solve(bBand,h_nextRow,m,(xn-2)*(yn-2))


  h_next(2:xn-1,2:yn-1) = transpose(reshape(h_nextRow, (/xn-2, yn-2/)))


write(*,*) "deltaT"
write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-h_next(2:xn-1,2:yn-1))/h_next(2:xn-1,2:yn-1)))

  return
end function h_next








! ----------------------------------------------------------------------------------%%
!
! PSI_NEXT
!
! SUMMARY: Computers the velocity field of the next timestep
!
! INPUTS: xn,yn : number of x,y gridpoints
!         x,y : values of x,y gridpoints
!         bcx0,bcy0 : boundary conditions
!         h : thermal profile of the previous timestep
!         rhs0 : source/sink term
!         u,v : x,y velocity field
!         dx,dy : x,y grid spacing
!         k : thermal diffusivity
!
! FUNCTIONS USED: Gsselm : gaussian elimination
!                 tridiag : solves tridiagonal system
!
! RETURNS: psi_next(xn,yn)
!
! ----------------------------------------------------------------------------------%%

function psi_next (h, rhs0, psi, top_in, rho_in)
use globals
use initialize
  implicit none

interface
  
function partial(array,rows,cols,d1,d2,dim)
use globals
use initialize
implicit none
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)
end function partial

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix

end function write_matrix
   
end interface

  ! integers
  integer :: i, j, ii, n, m
  ! inputs
  real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong((xn-2)*(yn-2))
  real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn)
  ! matrix stuff
  real(8) :: uVec((xn-2)*(yn-2)), psiLong((xn)*(yn)), psi_nextRow((xn-2)*(yn-2))
  real(8) :: psi_next(xn,yn), dsbit, fbit, conv
  real(8) :: way1((xn-2)*(yn-2)), way2((xn-2)*(yn-2))
  real(8) :: yep, top_in(xn,1)
  real(8) :: mn(xn,yn)
  !multigrid
  real(8) :: restr(3,3), interp(3,3), rh(xn,yn), eh(xn,yn)
  real(8) :: r2h((xn-1)/2,(xn-1)/2), e2h((xn-1)/2,(xn-1)/2), ee2h((xn-1)/2,(xn-1)/2)
  real(8) :: perm2h((xn-1)/2,(xn-1)/2), rho2h((xn-1)/2,(xn-1)/2)
  real(8) :: r4h(((xn-1)/2-1)/2,((xn-1)/2-1)/2), e4h(((xn-1)/2-1)/2,((xn-1)/2-1)/2)
 ! back to band
  real(8) :: aBand0((xn-2)*(yn-2),2*(xn-2) + 1)


mn = psi

!permLong = reshape(permeability,(/(xn-2)*(yn-2)/))

permx = partial((1/(permeability*rho_in)),xn,yn,dx,dy,1)
permy = partial((1/(permeability*rho_in)),xn,yn,dx,dy,2)

rhoLong = reshape(rho_in,(/(xn-2)*(yn-2)/))
permLong = reshape(permeability,(/(xn-2)*(yn-2)/))
permxLong = reshape(permx,(/(xn-2)*(yn-2)/))
permyLong = reshape(permy,(/(xn-2)*(yn-2)/))

 uVec = 0.0
 !uVec = reshape(rhs0,(/ (xn-2)*(yn-2) /))

 rhs1 = rhs0

 rhs1(2,:) = rhs1(2,:)
 rhs1(yn-1,:) = rhs1(xn-1,:)
 rhs1(:,2) = rhs1(:,2)
 rhs1(:,xn-1) = rhs1(:,xn-1)
 rhs1(2:xn-1,xn-1) = rhs1(2:xn-1,xn-1) !+&
 !& permeable(2:xn-1)/((2.0/(dx*dx*(permeability(2:xn-1,xn-1)*rho_in(2:xn-1,xn-1))))+&
 !&(2.0/(dy*dy*(permeability(2:xn-1,xn-1)*rho_in(2:xn-1,xn-1))))) 

  !aa0 = 0.0
  do i = 1,(xn-2)*(yn-2)
  	!aa0(i,i) = (4.0)/(dx*dx)
  	if (i .gt. 1) then
  	!aa0(i,i-1) = (-1.0)/(dx*dx) !- permLong(i)*permxLong(i)*dt0/(2.0*dx)
  	end if
  	if (i .lt. (xn-2)*(yn-2)) then
  	!aa0(i,i+1) = (-1.0)/(dx*dx) !+ permLong(i)*permxLong(i)*dt0/(2.0*dx)
  	end if
  	if (i-(xn-2) .gt. 0) then
  	!aa0(i,i-(xn-2)) = (-1.0*dt0)/(dx*dx) + permLong(i)*permyLong(i)*dt0/(2.0*dx)
  	!aa0(i-(xn-2),i) = (-1.0*dt0)/(dx*dx) + permLong(i)*permyLong(i)*dt0/(2.0*dx)
  	end if
  	if (i+(xn-2) .le. (xn-2)*(yn-2)) then
  	!aa0(i,i+(xn-2)) = (-1.0*dt0)/(dx*dx) - (permLong(i)*permyLong(i)*dt0)/(2.0*dx)
  	!aa0(i+(xn-2),i) = (-1.0*dt0)/(dx*dx) - (permLong(i)*permyLong(i)*dt0)/(2.0*dx)
  	end if
  end do
  
  do i = 1,((xn-2)-1)
    j = i*(xn-2)
    !aa0(j+1,j) = 0.0
  	!aa0(j,j+1) = 0.0
  end do

! BANDED WAY
  uVec = reshape(rhs1(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
!  aBand = band(aBand,m,(xn-2)*(yn-2))
!  psi_nextRow = solve(aBand,uVec,m,(xn-2)*(yn-2))
!  uVec = reshape(transpose(psi_next(2:xn-1,2:yn-1)),(/(xn-2)*(yn-2)/))
!  psi_nextRow = -solve(aBand,uVec,m,(xn-2)*(yn-2)) 
!  psi_next(2:xn-1,2:yn-1) = transpose(reshape(psi_nextRow, (/xn-2, yn-2/)))

psi_next = 0.0

!uVec = reshape(rhs1(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
!a0(1:(xn-2)*(yn-2),1:(xn-2)*(yn-2)) = aa0
!a0(:,(xn-2)*(yn-2)+1) = -uVec

!----------------------------------------------------------!
! MAKE THE BAND
 ! put this in for consistency permeLong(i)*rhoLong(i)*
  aBand0 = 0.0
  m = 2*(xn-2) + 1
  do i = 1,(xn-2)*(yn-2)
  aBand0(i,(m+1)/2) = (4.0)/(dx*dx)
  if (i .gt. 1) then
  	aBand0(i,((m+1)/2)-1) = (-1.0)/(dx*dx) !(-1.0)/(dx*dx)
  end if
  if (i .lt. (xn-2)*(yn-2)) then
  	aBand0(i,((m+1)/2)+1) = (-1.0)/(dx*dx) !(-1.0)/(dx*dx)
  end if
    ! extra columns
  if (i .le. (xn-2)*(yn-2)-(xn-2)) then
  	aBand0(i,m) = (-1.0)/(dx*dx) !+ permyLong(i)/(2.0*dx)
  end if
  if (i .ge. (xn-2)) then
  	aBand0(i,1) = (-1.0)/(dx*dx) !- permyLong(i)/(2.0*dx)
  end if
  end do
  
  do i=1,((xn-2)-1)
  ii = i*(xn-2)
  aBand0(ii,((m+1)/2)+1) = 0.0
  aBand0(ii+1,((m+1)/2)-1) = 0.0
  end do
  
! THIS IS FOR SOLVING IF THE MATRIX WORKED
  aBand0 = band(aBand0,m,(xn-2)*(yn-2))
  psi_nextRow = solve(aBand0,uVec,m,(xn-2)*(yn-2))
  psi_next(2:xn-1,2:yn-1) = reshape(psi_nextRow, (/xn-2, yn-2/))

  
!----------------------------------------------------------!
!multigridz
!----------------------------------------------------------!

restr(1,:) = (/.0625/2 +(.0625/2), .125, .0625/2 +(.0625/2)/)
restr(2,:) = (/.125, .250, .125/)
restr(3,:) = (/.0625/2 +(.0625/2), .125, .0625/2 +(.0625/2)/)

interp = 4.0*restr

!write(*,"(3F6.4)") restr
!write(*,"(3F6.4)") interp
  
  
!----------------------------------------------------------!

! ENTIRELY JACOBIFIED !

!write(*,"(F30.5)") permy*dx

!BENCHMARK
rho_in = 1.0
permeability = 1.0

do n=1,1000

do i=2,xn-1
do j=2,yn-1

	!psi_next(i,j) = (1.0/4.0)*(psi_next(i+1,j)+psi_next(i-1,j)&
	!&+psi_next(i,j+1)+psi_next(i,j-1)+dx*dx*rhs1(i,j))
	
	
	!psi_next(i,j) = ((1.0/(co*dx*dx))*psi_next(i+1,j)+(1.0/(co*dx*dx))*psi_next(i-1,j)&
!	&+permeability(i,j)*permy(i,j)*psi_next(i,j+1)*dx/(dy*dy) - permeability(i,j)*permy(i,j)*psi_next(i,j-1)*dx/(dy*dy)&
	!&+(1.0/(co*dy*dy))*psi_next(i,j+1)+(1.0/(co*dy*dy))*psi_next(i,j-1)+rhs1(i,j)/co)

! this is good.
!psi_next(i,j) = ((dx*dx*permeability(i,j)*rho_in(i,j))/4.0)&
!&*(psi_next(i+1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i-1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+rhs1(i,j))

!!!!!!!!! THIS BIT !!!!!!!!!!!!!
!psi_next(i,j) = (1/((4.0/(dx*dx*(permeability(i,j)*rho_in(i,j))))))&
!&*(psi_next(i+1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i-1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dy)&
!!&-(permx(i,j)*psi_next(i+1,j) + permx(i,j)*psi_next(i-1,j))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+rhs1(i,j))


end do
end do

end do

write(*,*) "deltaPSI"
write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))


!----------------------------------------------------------!
!multigridz
!----------------------------------------------------------!

!rh = rhs1 - psi_next
!
!r2h = 0.0
!do i=1,(xn-1)/2
!do j=1,(xn-1)/2
!	r2h(i,j) = sum(restr*rh(2*i:2*i+2,2*j:2*j+2))
!end do
!end do
!
!
!
!e2h(:,:) = 0.0
!
!do n=1,3
!do i=2,(xn-1)/2
!do j=2,(xn-1)/2
!e2h(i,j) = 1/((2.0/(4.0*dx*dx*(permeability(i/2,j/2)*rho_in(i/2,j/2))))+(2.0/(4.0*dy*dy*(permeability(i/2,j/2)*rho_in(i/2,j/2)))))&
!&*(e2h(i+1,j)/(permeability(i/2,j/2)*rho_in(i/2,j/2)*4.0*dx*dx)&
!&+e2h(i-1,j)/(permeability(i/2,j/2)*rho_in(i/2,j/2)*4.0*dx*dx)&
!&+(permy(i/2,j/2)*e2h(i,j+1) - permy(i/2,j/2)*e2h(i,j-1))/(2.0*dy)&
!&-(permx(i/2,j/2)*e2h(i+1,j) + permx(i/2,j/2)*e2h(i-1,j))/(2.0*dx)&
!&+e2h(i,j+1)/(permeability(i/2,j/2)*rho_in(i/2,j/2)*4.0*dy*dy)&
!&+e2h(i,j-1)/(permeability(i/2,j/2)*rho_in(i/2,j/2)*4.0*dy*dy)&
!&+r2h(i,j))
!end do
!end do
!end do
!
!
!
!do i=1,(xn-1)/2
!do j=1,(yn-1)/2
!eh(2*i,2*j) = e2h(i,j)
!eh(2*i,2*j+1) = (dy/dx)*0.5*(e2h(i,j) + e2h(i,j+1))
!
!eh(2*i+1,2*j) = (dx/dy)*0.5*(e2h(i,j) + e2h(i+1,j))
!eh(2*i+1,2*j+1) = (dx/dy)*0.25*(e2h(i,j) + e2h(i+1,j)) + .25*(dy/dx)*(e2h(i,j+1) + e2h(i+1,j+1))
!end do
!end do
!
!
!psi_next = psi_next + eh
!
!
!do n=1,100
!do i=2,xn-1
!do j=2,yn-1
!psi_next(i,j) = 1/((2.0/(dx*dx*(permeability(i,j)*rho_in(i,j))))+(2.0/(dy*dy*(permeability(i,j)*rho_in(i,j)))))&
!&*(psi_next(i+1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i-1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dy)&
!&-(permx(i,j)*psi_next(i+1,j) + permx(i,j)*psi_next(i-1,j))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+rhs1(i,j))
!end do
!end do
!if (n .eq. 100) then
!write(*,"(F30.26)") (maxval(abs(psi_next)) &
!&- maxval(abs((mn))))/(maxval(abs((mn))))
!end if
!mn = psi_next(2:xn-1,2:yn-1)
!end do



!psi_next(1,:) = ((4.0/3.0)*psi_next(2,:) - (1.0/3.0)*psi_next(3,:))
!psi_next(xn,:) = ((4.0/3.0)*psi_next(xn-1,:) - (1.0/3.0)*psi_next(xn-2,:))

return

end function psi_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY: Computes the updated fluid densities based on temperature profile
!
! INPUTS: h : temperature matrix, which is how density is calculated
!
! RETURNS: rho_next(xn,yn)
!
! ----------------------------------------------------------------------------------%%

function rho_next(h_in)
use globals
use initialize
  implicit none
  
  integer :: i,j
  real(8) :: h_in(xn,yn), rho_next(xn,yn)
  
  do i=1,xn
  do j = 1,yn
  rho_next(i,j) = rho_fluid*(1.0 - alpha*((h_in(i,j)-273.0)))
  end do 
  end do
  
  return

end function rho_next


! ----------------------------------------------------------------------------------%%
!
! V_NEXT
!
! SUMMARY: for psi boundary condish
!
! INPUTS: 
!
! RETURNS: 
!
! ----------------------------------------------------------------------------------%%

function v_next(dPdy_in, rho_in)
use globals
use initialize
  implicit none
  
  integer :: i,j
  real(8) :: dPdy_in(xn,yn), rho_in(xn,yn), v_next(xn,yn)

  v_next = -(permeability/viscosity)*(dPdy_in-rho_in*g)

  
  return

end function v_next





! ----------------------------------------------------------------------------------%%
!
! VELOCITIES
!
! Updating this function to use velocity formulation instead of the streamfunction
! formulation.
!
! ----------------------------------------------------------------------------------%%


function velocities(psi)
use globals
use initialize
implicit none

interface

function partial(array,rows,cols,d1,d2,dim)
use globals
use initialize
implicit none
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)
end function partial

end interface

real(8) :: velocities(xn,2*yn), psi(xn,yn)
real(8) :: u0(xn,yn), v0(xn,yn)


u0 = partial(psi,xn,yn,dx,dy,2)
v0 = -partial(psi,xn,yn,dx,dy,1)

velocities(1:xn,1:yn) = u0
velocities(1:xn,yn+1:2*yn) = v0

return
end function velocities



! ----------------------------------------------------------------------------------%%
!
! CONTINUITY
!
! Updating this function to use velocity formulation instead of the streamfunction
! formulation.
!
! ----------------------------------------------------------------------------------%%

function continuity(u_in,v_in)
use globals
use initialize
implicit none

real(8) :: u_in(xn,yn), v_in(xn,yn)
real(8) :: u_out(xn,yn), v_out(xn,yn)
real(8) :: continuity(xn,2*yn)

continuity(1:xn,1:yn) = u_out
continuity(1:xn,yn+1:2*yn) = v_out


return
end function continuity


! ----------------------------------------------------------------------------------%%
!
! PARTIAL
!
! second-order accurate partial derivative of given array with respect to
! given dimension
!
! ----------------------------------------------------------------------------------%%

function partial(array,rows,cols,d1,d2,dim)
use globals
use initialize
implicit none
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)


if (dim .eq. 1) then
ii = 1
jj = 0
d = d1
partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
partial(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)

!partial(2,:) = ( -3.0*array(2,:) + 4.0*array(3,:) -array(4,:)) / (2.0*d)
!partial(rows-1,:) = ( 3.0*array(rows-1,:) - 4.0*array(rows-2,:) + array(rows-3,:) ) / (2.0*d)


end if

if (dim .eq. 2) then
ii = 0
jj = 1
d = d2
partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)

!partial(:,2) = ( -3.0*array(:,2) + 4.0*array(:,3) -array(:,4)) / (2.0*d)
!partial(:,cols-1) = ( 3.0*array(:,cols-1) - 4.0*array(:,cols-2) + array(:,cols-3) ) / (2.0*d)


end if


    do i = 2-jj,rows-1+jj
    do j = 2-ii,cols-1+ii
    partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
   
    
    ! positive derivative gets upwind forward
    !if (abs((array(i+ii,j+jj)-array(i,j))/d) .gt. abs((array(i,j)-array(i-ii,j-jj))/d)) then
    !	partial(i,j) = ( -3.0*array(i,j) + 4.0*array(i+ii,j+jj) -array(i+2*ii,j+2*jj)) / (2.0*d)
	!end if
	
	! negative derivative gets upwind backward
    !if (abs((array(i+ii,j+jj)-array(i,j))/d) .lt. abs((array(i,j)-array(i-ii,j-jj))/d)) then
    !	partial(i,j) = ( 3.0*array(i,j) - 4.0*array(i-ii,j-jj) +array(i-2*ii,j-2*jj)) / (2.0*d)
	!end if
	
	end do
    end do

return
end function partial






! ----------------------------------------------------------------------------------%%
!
! WRITE_VEC
!
! SUMMARY: Write linspace-style vector to file
!
! INPUTS: n : dimension
!         vector : vector with data
!         filename : file name
!
! RETURNS: write_vec
!
! ----------------------------------------------------------------------------------%%

function write_vec ( n, vector, filename )
use globals
  implicit none
  integer :: n, j, output_status, unit0
  character ( len = * ) filename 
  real(4)  :: vector(n), write_vec



  unit0 = get_unit ()
  open ( unit = unit0, file = filename, status = 'replace', iostat = output_status )
  if ( output_status /= 0 ) then
    write ( *, '(a,i8)' ) 'COULD NOT OPEN OUTPUT FILE "' // &
      trim ( filename ) // '" USING UNIT ', unit0
    unit0 = -1
    stop
  end if
  

  if ( 0 < n ) then
    do j = 1, n
      write ( unit0, '(2x,g24.16)' ) vector(j)
    end do

  end if


  close ( unit = unit0 )
  write_vec = 1.0
  return
end function write_vec




! ----------------------------------------------------------------------------------%%
!
! WRITE_MATRIX
!
! SUMMARY: Write 2d array to file
!
! INPUTS: m,n : 2d dimensinons
!         table : 2d array with data
!         filename : file name
!
! RETURNS: write_matrix
!
! ----------------------------------------------------------------------------------%%

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix

  unit0 = get_unit ()
  open ( unit = unit0, file = filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) 'Could not open the output file "' // &
      trim ( filename ) // '" on unit ', unit0
    unit0 = -1
    stop
  end if


	write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'

    do j = 1, n
      write ( unit0, string ) table(1:m,j)
    end do


  close ( unit = unit0 )
  write_matrix = 2.0
  return
end function write_matrix

