
! ----------------------------------------------------------------------------------%%
!
! MAIN MASSACR METHOD
! 
! SUMMARY: Solve governing equations (conservation of momentum and thermal energy) 
!          in 2D and writes output to text (.txt) and netCDF (.nc) files. 
! 
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
	real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
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

end interface



	real(8) :: cfl, test(10,11)
	real(8) :: h(xn,yn), psi(xn,yn) ! xn ROWS DEEP & yn COLUMNS WIDE 
	! real(8) :: hmat(xn,(yn*tn)), psimat(xn,(yn*tn)), umat(xn,(yn*tn)), vmat(xn,(yn*tn))
	real(8) :: umat(xn,(yn)), vmat(xn,(yn))
	! hmat(xn,(yn)), psimat(xn,(yn)),
	real(8) :: hmat(xn,(yn*tn/10)), psimat(xn,(yn*tn/10))
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
	real(8) :: hc=20.0



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



do j = 2, tn

	write(*,*) j

	! HEAT FLUX BOUNDARY CONDITIONS

		! bottom
		do i = 1,xn
		flux(i,1) = h(i,2) +((.27)+0.1*i)*dy/(lambda)
		!flux(i,1) = 400.0
		end do

		! top
		do i = 1,xn
		flux(i,2) = h(i,xn-1) -(.27)*dy/(lambda)
		!flux(i,2) = ( 273.0 + h(i,xn-1) ) / 2
		end do
  
	! SOLVE THERMAL NRG EQUATION
	rho = rho_next(h)
	h = h_next(h, psi,rho, flux)

	!!!!!!!!!!!! THIS !!!!!!!!!!!
	h(1,:) = (4.0/3.0)*h(2,:) - (1.0/3.0)*h(3,:) ! left
	h(xn,:) = (4.0/3.0)*h(xn-1,:) - (1.0/3.0)*h(xn-2,:) ! right
	h(:,1) = flux(:,1)
	h(:,yn) = flux(:,2)

	! ! SOLVE STREAMFUNCTION-VORTICITY EQUATION
	! rhs0 = - (1.0/(viscosity))*g*rho_fluid*alpha*partial(h,xn,yn,dx,dy,1)
	! psi = psi_next(h, rhs0, psi, permeable, rho)

	! PUT IN BOUNDARY CONDITIONS BETWEEN STEPS
	psi(1,1:yn) = bcyPsi(1,1:yn) ! left
	psi(xn,1:yn) = bcyPsi(2,1:yn) ! right
	psi(1:xn,1) = bcxPsi(1:xn,1) ! bottom
	psi(1:xn,yn) = bcxPsi(1:xn,2) ! top
	psi(:,yn) = ((4.0/3.0)*psi(:,yn-1) - (1.0/3.0)*psi(:,yn-2)) 
	permeable = psi(:,yn)

	! VELOCITIES
	velocities0 = velocities(psi)
	u = velocities0(1:xn,1:yn)/rho
	v = velocities0(1:xn,yn+1:2*yn)/rho

	! ADD EACH TIMESTEP TO MATRICES
	if (mod(j,10) .eq. 0) then
		 hmat(1:xn,1+yn*(j/10-1):1+yn*(j/10)) = h
		 psimat(1:xn,1+yn*(j/10-1):1+yn*(j/10)) = psi
	end if
	! umat(1:xn,1+yn*(j-1):1+yn*(j)) = u
	! vmat(1:xn,1+yn*(j-1):1+yn*(j)) = v
	! hmat(1:xn,1:yn) = h
	! psimat(1:xn,1:yn) = psi
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
yep = write_vec ( xn, real(x,kind=4), 'x1.txt' )
yep = write_vec ( yn, real(y,kind=4), 'y1.txt' )
yep = write_vec ( tn, real(t, kind=4), 't1.txt' )
yep = write_matrix ( xn, yn*tn/10, real(hmat, kind = 4), 'hMat.txt' )
yep = write_matrix ( xn, yn*tn/10, real(psimat,kind=4), 'psiMat.txt' )
!yep = write_matrix ( xn, yn*tn, real(umat,kind=4), 'uMat.txt' )
!yep = write_matrix ( xn, yn*tn, real(vmat,kind=4), 'vMat.txt' )
!yep = write_matrix ( xn, yn, real(hmat, kind = 4), 'h3.txt' )
!yep = write_matrix ( xn, yn, real(psimat,kind=4), 'psiMat3.txt' )
yep = write_matrix ( xn, yn, real(umat,kind=4), 'uMat1.txt' )
yep = write_matrix ( xn, yn, real(vmat,kind=4), 'vMat1.txt' )
yep = write_matrix ( xn, yn, real(rho,kind=4), 'rho1.txt' )
yep = write_matrix ( xn, yn,real(permeability,kind=4), 'permeability1.txt' )

! WRITE THINGS TO NETCDF FILES
! call check( nf90_create('thermalNRG.nc', NF90_CLOBBER, ncid) )
! call check( nf90_def_dim(ncid, "x", xn, x_dimid) )
! call check( nf90_def_dim(ncid, "y", yn, y_dimid) )

! call check(nf90_def_var(ncid, "h", NF90_FLOAT, (/ x_dimid, y_dimid /), h_varid) )
! call check(nf90_def_var(ncid, "u", NF90_FLOAT, (/ x_dimid, y_dimid /), u_varid) )
! call check(nf90_def_var(ncid, "v", NF90_FLOAT, (/ x_dimid, y_dimid /), v_varid) )
! call check(nf90_enddef(ncid))

! call check( nf90_put_var(ncid, h_varid, h) )
! call check( nf90_put_var(ncid, u_varid, u) )
! call check( nf90_put_var(ncid, v_varid, v) )

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
! INPUTS: h(xn,yn) : temperature profile of previous timestep
!         psi(xn,yn) : 2D streamfunction array
!         rho_in(xn,yn) : 2D density array
!         flux(xn,2) : top and bottom heat flux boundary conditions
!
! RETURNS: h_next(xn,yn) : temperature profile of current timestep
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
  ! real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
  ! real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
  real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
  real(8) :: kMatLong((xn-2)*(yn-2))
  real(8) :: mn(xn,yn)
  
mn = h

velocities0 = velocities(psi)
u = velocities0(1:xn,1:yn)
v = velocities0(1:xn,yn+1:2*yn)
uLong = reshape(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
vLong = reshape(transpose(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! PRINT CONVERGENCE CONDITIONS AT EACH TIMESTEP
write(*,*) " "
write(*,*) "max u"
write(*,*) maxval(abs(u))
write(*,*) "max v"
write(*,*) maxval(abs(v))

write(*,*) " "
write(*,*) "velocity check"
write(*,"(F10.5)") (dt*maxval(abs(u)))/(dx)
write(*,"(F10.5)") (dt*maxval(abs(v)))/(dy)
write(*,*) "conduction check"
write(*,"(F10.5)") (2.0*dt*lambda)/(4186.0*dy*dy)
write(*,*) " "

h0 = h

qx = dt/(dx)
qy = dt/(dy)
sx = (2.0*dt*lambda)/(4186.0*dx*dx)
sy = (2.0*dt*lambda)/(4186.0*dy*dy)

! VERTICAL BOUNDARY CONDITIONS
 h(2,:) = h(2,:) + h0(1,:)*sx/2.0  ! left
 h(yn-1,:) = h(yn-1,:) + h0(xn,:)*sx/2.0  ! right
 
 ! STRETCH OUT LAST TIMESTEP
uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))


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

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
	aBand(i,2) = 1.0 + sy + uLong(i)*qy
	if (i .gt. 1) then
	aBand(i,1) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	aBand(i,3) = -sy/2.0 - uLong(i)*qy
	end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
	aBand(i,2) = 1.0 + sy - uLong(i)*qy
	if (i .gt. 1) then
	aBand(i,3) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	aBand(i,1) = -sy/2.0 + uLong(i)*qy
	end if
	end if
  
end do

! PUT IN EMPTY SPACE FOR BOUNDARY CONDITIONS
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do
  
  
!!!!!!!!!!!! THIS !!!!!!!!!!!
h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))


! HORIZONTAL BOUNDARY CONDITIONS
h(:,2) = h(:,2) + flux(:,1)*sy/2.0 ! bottom
h(:,xn-1) = h(:,xn-1) + flux(:,2)*sy/2.0 ! top


h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))


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

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
	bBand(i,2) = 1.0 + sy + vLong(i)*qy
	if (i .gt. 1) then
	bBand(i,1) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	bBand(i,3) = -sy/2.0 - vLong(i)*qy
	end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
	bBand(i,2) = 1.0 + sy - vLong(i)*qy
	if (i .gt. 1) then
	bBand(i,3) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	bBand(i,1) = -sy/2.0 + vLong(i)*qy
	end if
	end if
	
 end do
  
! PUT IN EMPTY SPACE FOR BOUNDARY CONDITIONS
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	bBand(ii,3) = 0.0
	bBand(ii+1,1) = 0.0
end do

h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))
h_next(2:xn-1,2:yn-1) = transpose(reshape(h_nextRow, (/xn-2, yn-2/)))

! PRINT DELTA T
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
! INPUTS: h(xn,yn) : temperature profile
!         rhs0(xn,yn) : right hand side of streamfunction-vorticity equation
!         psi(xn,yn) : 2D streamfunction array of previous timestep
!         top_in(xn,1) : permeable upper boundary
!         rho_in(xn,yn) : 2D density array
!
! RETURNS: psi_next(xn,yn): 2D streamfunction array for current timestep
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
 ! back to band
  real(8) :: aBand0((xn-2)*(yn-2),2*(xn-2) + 1)


mn = psi

permx = partial((1/(permeability*rho_in)),xn,yn,dx,dy,1)
permy = partial((1/(permeability*rho_in)),xn,yn,dx,dy,2)

rhoLong = reshape(rho_in(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permLong = reshape(permeability(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permxLong = reshape(permx(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permyLong = reshape(permy(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))


rhs1 = rhs0

rhs1(2,:) = rhs1(2,:) 
rhs1(yn-1,:) = rhs1(xn-1,:) 
rhs1(:,2) = rhs1(:,2) 
rhs1(:,xn-1) = rhs1(:,xn-1)
rhs1(2:xn-1,xn-1) = rhs1(2:xn-1,xn-1) +&
& permeable(2:xn-1)/(4.0/(dx*dx*(permeability(2:xn-1,xn-1)*rho_in(2:xn-1,xn-1))))

uVec = reshape(rhs1(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))

psi_next = 0.0


! MAKE THE BAND
aBand0 = 0.0
m = 2*(xn-2) + 1
do i = 1,(xn-2)*(yn-2)
	  
	aBand0(i,(m+1)/2) = (4.0)/(permLong(i)*rhoLong(i)*dx*dx)
	
	if (i .gt. 1) then
	aBand0(i,((m+1)/2)-1) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx)
	end if
	
	if (i .lt. (xn-2)*(yn-2)) then
	aBand0(i,((m+1)/2)+1) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx)
	end if
	
	! extra columns
	if (i .le. (xn-2)*(yn-2)-(xn-2)) then
	aBand0(i,m) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx) - (permyLong(i))/(2.0*dx)
	end if
	
	if (i .ge. (xn-2)) then
	aBand0(i,1) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx) + (permyLong(i))/(2.0*dx)
	end if
  
end do
  
! PUT IN EMPTY SPACE FOR BOUNDARY CONDITIONS
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	aBand0(ii,((m+1)/2)+1) = 0.0
	aBand0(ii+1,((m+1)/2)-1) = 0.0
end do
  
! THIS IS FOR SOLVING
aBand0 = band(aBand0,m,(xn-2)*(yn-2))
psi_nextRow = solve(aBand0,uVec,m,(xn-2)*(yn-2))
psi_next(2:xn-1,2:yn-1) = reshape(psi_nextRow, (/xn-2, yn-2/))


  

! ENTIRELY JACOBIFIED !

!do n=1,800
!do i=2,xn-1
!do j=2,yn-1

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
!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dy)&
!&-(permx(i,j)*psi_next(i+1,j) + permx(i,j)*psi_next(i-1,j))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+rhs1(i,j))

!end do
!end do
!end do

write(*,*) "deltaPSI"
write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))

return

end function psi_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
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
! VELOCITIES
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
end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
	partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
end if


do i = 2-jj,rows-1+jj
	do j = 2-ii,cols-1+ii
		partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
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

