
! ----------------------------------------------------------------------------------%%
!
! MAIN MASSACR METHOD
! 
! SUMMARY: main method runs fluid dynamic simulation coupled to geochemical
!          simulation, writes errthing to file
! 
! gfortran -c -O3 -I/usr/local/include -L/usr/local/lib -lnetcdff -liphreeqc 
! globals.f90 initialize.f90 alteration.f90 massacr.f90
!
! gfortran -I/usr/local/include -L/usr/local/lib -lnetcdff -liphreeqc 
! globals.o initialize.o alteration.o massacr.o
!
! I HAVE A MAKEFILE NOW
!
!
! MOST BASIC PARALLELIZING
! make -f theMakeFile
! mpirun -np 1 ./massacr
!
! ----------------------------------------------------------------------------------%%

PROGRAM main
!use netcdf

use globals
use initialize
use alteration

implicit none

include 'mpif.h'



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
	real(8) :: psi_next(xn,yn)
	real(8) :: top_in(xn,1)
	! tridiag
	real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
	real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
	end function psi_next
	
	function alt_next (temp, timestep, primaryList, secondaryList, soluteList)
	use globals
	use initialize
	use alteration
	! declare yo shit
	real(8) :: temp, timestep
	real(8) :: alt_next(1,58)
	real(8) :: alter0(1,58)
	real(8) :: primaryList(5), secondaryList(16), soluteList(11)
	end function alt_next

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

! DECLARE EVERYTHING
real(8) :: h(xn,yn), psi(xn,yn) ! xn ROWS DEEP & yn COLUMNS WIDE 
! real(8) :: hmat(xn,(yn*tn)), psimat(xn,(yn*tn)) 
! real(8) :: umat(xn,(yn*tn)), vmat(xn,(yn*tn))
! real(8) :: umat(xn,(yn)), vmat(xn,(yn))
! hmat(xn,(yn)), psimat(xn,(yn)),
real(8) :: hmat(xn,(yn*tn/mstep)), psimat(xn,(yn*tn/mstep))
real(8) :: umat(xn,(yn*tn/mstep)), vmat(xn,(yn*tn/mstep))
real(8) :: rhs0(xn,yn), velocities0(xn,2*yn)
real(8) :: u(xn,yn), uLong(xn*yn), v(xn,yn), vLong(xn*yn)
real(8) :: rho(xn,yn), flux(xn,2)

integer :: unit
real(8) :: yep

integer :: xInt, yInt, tInt, hInt, uInt, vInt
integer :: ncid
integer :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
integer :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid
integer :: i, j, ii, m, n, jj

real(8) :: nusseltLocalv(xn,1), nuBar
real(8) :: alt0(1,altnum)

! ALTERATION ARRAYS
real(8) :: primary(xn/cell,yn/cell,5), secondary(xn/cell,yn/cell,16), solute(xn/cell,yn/cell,11)
real(8) :: primaryMat(xn/cell,yn*tn/(cell*mstep),5), secondaryMat(xn/cell,yn*tn/(cell*mstep),16)
real(8) :: soluteMat(xn/cell,yn*tn/(cell*mstep),11)

! DECLARE MPI STUFF
!integer :: ierr
integer, parameter :: max_rows = 10000000
integer, parameter :: send_data_tag = 2001, return_data_tag = 2002
integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
integer :: num_procs, an_id, num_rows_to_receive
integer :: avg_rows_per_process, num_rows, num_rows_to_send
integer :: end_row, sender, start_row, num_rows_received
real(8) :: vector(max_rows), vector2(max_rows), partial_sum, sum
real(8) :: local_mean, global_mean
real(8) :: hLocal((xn/cell)*(yn/cell)), dt_local


! MPI STRETCHED OUT ARRAYS
real(8) :: hLong((xn/cell)*(yn/cell))
real(8) :: priLong((xn/cell)*(yn/cell),5), secLong((xn/cell)*(yn/cell),16), solLong((xn/cell)*(yn/cell),11)
real(8) :: priLocal((xn/cell)*(yn/cell),5), secLocal((xn/cell)*(yn/cell),16), solLocal((xn/cell)*(yn/cell),11)
real(8) :: priLongBit((xn/cell)*(yn/cell)), priLocalBit((xn/cell)*(yn/cell))
real(8) :: secLongBit((xn/cell)*(yn/cell)), secLocalBit((xn/cell)*(yn/cell))
real(8) :: solLongBit((xn/cell)*(yn/cell)), solLocalBit((xn/cell)*(yn/cell))


! INITIALIZE CHEMISTRY
primary(:,:,1) = 12.96 ! feldspar
primary(:,:,2) = 6.96 ! augite
primary(:,:,3) = 1.26 ! pigeonite
primary(:,:,4) = .4 ! magnetite
primary(:,:,5) = 96.77 ! basaltic glass

secondary(:,:,:) = 0.0

solute(:,:,1) = 7.5 ! ph
solute(:,:,2) = 6.0e-4 ! Ca
solute(:,:,3) = 2.0e-5 ! Mg
solute(:,:,4) = 1.0e-3 ! Na
solute(:,:,5) = 1.0e-4 ! K
solute(:,:,6) = 1.2e-6 ! Fe
solute(:,:,7) = 1.0e-4 ! S(6)
solute(:,:,8) = 2.0e-4 ! Si
solute(:,:,9) = 3.0e-4 ! Cl
solute(:,:,10) = 1.0e-6 ! Al
solute(:,:,11) = 2.0e-3 ! Alk






!-----------------------MESSAGE PASSING-----------------------!
! process #0 is the root process
root_process = 0

! initialize a process
call MPI_INIT ( ierr )

! find out the process ID and how many processes were started so far
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

write(*,*) "my_id:", my_id
write(*,*) " "

! what to do if the process is the root process
if (my_id .eq. root_process) then
!-----------------------MESSAGE PASSING-----------------------!


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





! THIS IS THE MAIN LOOP
do j = 2, tn
write(*,*) j

	! HEAT FLUX BOUNDARY CONDITIONS
	! bottom
	do i = 1,xn
	flux(i,1) = h(i,2) +((.27))*dy/(2.6)
	end do
	! top
	do i = 1,xn
	flux(i,2) = 273.0 !+ 0.005*x(i)
	end do
  
	! SOLVE THERMAL NRG EQUATION
	rho = rho_next(h)
	h = h_next(h, psi,rho, flux)
  
	! PUT HEAT FLUX BOUNDARY CONDITIONS IN ARRAY
	h(1,:) = 273.0 !(4.0/3.0)*h(2,:) - (1.0/3.0)*h(3,:) ! left
	h(xn,:) = (4.0/3.0)*h(xn-1,:) - (1.0/3.0)*h(xn-2,:) ! right
	h(:,1) = flux(:,1)
	h(:,yn) = flux(:,2)
  

	! SOLVE STREAMFUNCTION-VORTICITY EQUATION
	rhs0 = (1.0/(viscosity))*g*rho_fluid*alpha*partial(h,xn,yn,dx,dy,1)
	psi = psi_next(h, rhs0, psi, permeable, rho)

	! PUT IN STREAMFUNCTION BOUNDARY CONDITIONS BETWEEN STEPS
	psi(1,1:yn) = bcyPsi(1,1:yn) ! left
	psi(xn,1:yn) = bcyPsi(2,1:yn) ! right
	psi(1:xn,1) = bcxPsi(1:xn,1) ! bottom
	psi(1:xn,yn) = bcxPsi(1:xn,2) ! top
	psi(:,yn) = ((4.0/3.0)*psi(:,yn-1) - (1.0/3.0)*psi(:,yn-2))/1.0
	permeable = psi(:,yn)

	! GET VELOCITIES
	velocities0 = velocities(psi)
	u = velocities0(1:xn,1:yn)/rho
	v = velocities0(1:xn,yn+1:2*yn)/rho
	
! THINGS DONE ONLY EVERY mTH TIMESTEP GO HERE
if (mod(j,mstep) .eq. 0) then
	
	
	! stretch everything out
	hLong = reshape(h(1:xn-1:cell,1:yn-1:cell), (/(xn/cell)*(yn/cell)/))
	priLong = reshape(primary, (/(xn/cell)*(yn/cell), 5/))
	secLong = reshape(secondary, (/(xn/cell)*(yn/cell), 16/))
	solLong = reshape(solute, (/(xn/cell)*(yn/cell), 11/))
	
	!-----------------------MESSAGE PASSING-----------------------!
	
	! DISTRIBUTE TO SLAVE PROCESSORS
	do an_id = 1, num_procs - 1
		
		! put number of rows in vector here for hLong
		num_rows = (xn/cell)*(yn/cell)
		avg_rows_per_process = num_rows / (num_procs-1)
        start_row = ( (an_id-1) * avg_rows_per_process) + 1
        end_row = start_row + avg_rows_per_process - 1
        if (an_id .eq. (num_procs - 1)) end_row = num_rows
        num_rows_to_send = (end_row - start_row + 1)
		
		! send size of h chunk
        call MPI_SEND( num_rows_to_send, 1, MPI_INTEGER, &
		an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
		! send timestep size
        call MPI_SEND( dt, 1, MPI_DOUBLE_PRECISION, &
		an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
		! send h chunk to the an_id
        call MPI_SEND( hLong(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
		an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
		! send primary chunk to the an_id
		do ii = 1,5
			priLongBit = priLong(:,ii)
        	call MPI_SEND( priLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send secondary chunk to the an_id
		do ii = 1,16
			secLongBit = secLong(:,ii)
        	call MPI_SEND( secLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send solute chunk to the an_id
		do ii = 1,11
			solLongBit = solLong(:,ii)
        	call MPI_SEND( solLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		end do
		write(*,*) "DONE SENDING TO PROCESSOR", an_id
		
     end do

	write(*,*) "BEFORE"
	write(*,*) priLong(:,5)
	
	! RECEIVE EVERYTHING FROM SLAVE PROCESSORS HERE
	do an_id = 1, num_procs - 1
		
		! get the size of each chunk again
		num_rows = (xn/cell)*(yn/cell)
		avg_rows_per_process = num_rows / (num_procs-1)
        start_row = ( (an_id-1) * avg_rows_per_process) + 1
        end_row = start_row + avg_rows_per_process - 1
        if (an_id .eq. (num_procs - 1)) end_row = num_rows
        num_rows_to_send = (end_row - start_row + 1)
		
		! primary chunk
		do ii = 1,5
			! receive it
			call MPI_RECV( priLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			! fill it
			priLong(start_row:end_row,ii) = priLocal(1:num_rows_to_send,ii)
		end do
		
		! secondary chunk
		do ii = 1,16
			! receive it
			call MPI_RECV( secLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			! fill it
			secLong(start_row:end_row,ii) = secLocal(1:num_rows_to_send,ii)
		end do
		
		! solute chunk
		do ii = 1,11
			! receive it
			call MPI_RECV( solLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			! fill it
			solLong(start_row:end_row,ii) = solLocal(1:num_rows_to_send,ii)
		end do

		write(*,*) "DONE RECEIVING FROM PROCESSOR", an_id
	end do
	
	!-----------------------MESSAGE PASSING-----------------------!
	



	! put stretched vectors back into 2d arrays
	primary = reshape(priLong,(/(xn/cell), (yn/cell), 5/))
	secondary = reshape(secLong,(/(xn/cell), (yn/cell), 16/))
	solute = reshape(solLong,(/(xn/cell), (yn/cell), 11/))
	
	! LOOK AT PRIMARY SEE IF IT WORKS
	write(*,*) "AFTER"
	write(*,*) priLong(:,5)

	! ADD EACH TIMESTEP TO MATRICES
	 hmat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = h
	 psimat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = psi
	 umat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = u
	 vmat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = v
	 primaryMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = primary
	 secondaryMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = secondary
	 soluteMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = solute
	 
end if ! END mTH TIMESTEP LOOP

! umat(1:xn,1+yn*(j-1):1+yn*(j)) = u
! vmat(1:xn,1+yn*(j-1):1+yn*(j)) = v
! hmat(1:xn,1:yn) = h
! psimat(1:xn,1:yn) = psi
! umat(1:xn,1:yn) = u
! vmat(1:xn,1:yn) = v


end do ! END ALL TIMESTEP LOOP





! WRITE EVERYTHING TO FILE
yep = write_vec ( xn, real(x,kind=4), 'x.txt' )
yep = write_vec ( yn, real(y,kind=4), 'y.txt' )
yep = write_vec ( tn, real(t, kind=4), 't.txt' )
yep = write_matrix ( xn, yn*tn/mstep, real(hmat, kind = 4), 'hMat.txt' )
yep = write_matrix ( xn, yn*tn/mstep, real(psimat,kind=4), 'psiMat.txt' )
yep = write_matrix ( xn, yn*tn/mstep, real(umat, kind = 4), 'uMat.txt' )
yep = write_matrix ( xn, yn*tn/mstep, real(vmat,kind=4), 'vMat.txt' )

yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,1),kind=4), 'feldsparMat.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,5),kind=4), 'glassMat.txt' )


! WRITE TO FILE FOR LAST TIMESTEP CASE
!yep = write_matrix ( xn, yn*tn, real(umat,kind=4), 'uMat.txt' )
!yep = write_matrix ( xn, yn*tn, real(vmat,kind=4), 'vMat.txt' )
!yep = write_matrix ( xn, yn, real(hmat, kind = 4), 'h3.txt' )
!yep = write_matrix ( xn, yn, real(psimat,kind=4), 'psiMat3.txt' )
!yep = write_matrix ( xn, yn, real(umat,kind=4), 'uMat1.txt' )
!yep = write_matrix ( xn, yn, real(vmat,kind=4), 'vMat1.txt' )

yep = write_matrix ( xn, yn, real(rho,kind=4), 'rho.txt' )
yep = write_matrix ( xn, yn,real(permeability,kind=4), 'permeability.txt' )

! WRITE THINGS TO NETCDF FILES
! call check( nf90_create('thermalNRG.nc', NF90_CLOBBER, ncid) )
! call check( nf90_def_dim(ncid, "x", xn, x_dimid) )
! call check( nf90_def_dim(ncid, "y", yn, y_dimid) )

! call check(nf90_def_var(ncid, "h", NF90_FLOAT, (/ x_dimid, y_dimid /), h_varid) )
! call check(nf90_def_var(ncid, "u", NF90_FLOAT, (/ x_dimid, y_dimid /), u_varid) )
! call check(nf90_def_var(ncid, "v", NF90_FLOAT, (/ x_dimid, y_dimid /), v_varid) )
! call check(nf90_enddef(ncid))

!call check( nf90_put_var(ncid, h_varid, h) )
!call check( nf90_put_var(ncid, u_varid, u) )
!call check( nf90_put_var(ncid, v_varid, v) )

!  
! CALL check(nf90_put_att(ncid, x_varid, "units", "meter")) 
! CALL check(nf90_put_att(ncid, y_varid, "units", "meter")) 
! CALL check(nf90_put_att(ncid, h_varid, "units", "K"))
! CALL check(nf90_put_att(ncid, u_varid, "units", "m/s"))
! CALL check(nf90_put_att(ncid, v_varid, "units", "m/s"))
!  
! CALL check(nf90_put_att(ncid, x_varid, "axis", "X")) 
! CALL check(nf90_put_att(ncid, y_varid, "axis", "Y")) 


! Close the file. This frees up any internal netCDF resources
! associated with the file, and flushes any buffers.
!  call check( nf90_close(ncid) )

write(*,*) " "
write(*,*) "ALL DONE!"






!-----------------------MESSAGE PASSING-----------------------!
else
	do jj = 1, tn/mstep
		! here is a slave process, each process must receive a chunk of the h array and 
		! take the local mean, print it, send it back.
		
		! receive size of chunk
		call MPI_RECV ( num_rows_to_receive, 1 , MPI_INTEGER, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive timestep
		call MPI_RECV ( dt_local, 1 , MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive h chunk, save in local hLocal
		call MPI_RECV ( hLocal, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		num_rows_received = num_rows_to_receive

		! receive primary chunk, save in local priLocal
		do ii = 1,5
			call MPI_RECV ( priLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			priLocal(:,ii) = priLocalBit
		end do
	
		! receive secondary chunk, save in local priLocal
		do ii = 1,16
			call MPI_RECV ( secLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			secLocal(:,ii) = secLocalBit
		end do
	
		! receive solute chunk, save in local priLocal
		do ii = 1,11
			call MPI_RECV ( solLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			solLocal(:,ii) = solLocalBit
		end do
	
	!-----------------------MESSAGE PASSING-----------------------!

		!  DO ALL THE SLAVEWORK HERE
	
		! slave processor goes through phreeqc loop here
		do m=1,num_rows_to_receive
			!alt0 = alt_next(hLocal(m),dt_local*mstep,priLocal(m,:),secLocal(m,:),solLocal(m,:))

			!PARSING
			solLocal(m,:) = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
			alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12) /)

			secLocal(m,:) = (/ alt0(1,13), alt0(1,15), alt0(1,17), alt0(1,19), alt0(1,21), &
			alt0(1,23), alt0(1,25), alt0(1,27), alt0(1,29), alt0(1,31), alt0(1,33), alt0(1,35), &
			alt0(1,37), alt0(1,39), alt0(1,41), alt0(1,43)/)

			priLocal(m,:) = (/ alt0(1,45), alt0(1,47), alt0(1,49), alt0(1,51), alt0(1,53)/)
		end do
	
	
	!-----------------------MESSAGE PASSING-----------------------!
	
		! send primary chunk back to root process
		do ii = 1,5
			call MPI_SEND( priLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send secondary chunk back to root process
		do ii = 1,16
			call MPI_SEND( secLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send solute chunk back to root process
		do ii = 1,11
			call MPI_SEND( solLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		write(*,*) "SLAVE PROCESS IS DONE WITH WORK"
	
	end do
	
end if ! END MESSAGE PASSING LOOP
call MPI_FINALIZE ( ierr )
!-----------------------MESSAGE PASSING-----------------------!





END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: computes the 2D temperature profile for the current timestep
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
real(8) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
  
  
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
write(*,"(F10.5)") (2.0*dt*maxval(lambdaMat))/(4186.0*dy*dy)
write(*,*) " "



h0 = h

qx = dt/(dx)
qy = dt/(dy)
sx = (2.0*dt*lambda)/(4186.0*dx*dx)
sy = (2.0*dt*lambda)/(4186.0*dy*dy)
sxMat = (2.0*dt*lambdaMat)/(4186.0*dx*dx)
syMat = (2.0*dt*lambdaMat)/(4186.0*dy*dy)

! VERTICAL BOUNDARY CONDITIONS
h(2,:) = h(2,:) + h0(1,:)*sxMat(1,:)/2.0  ! left
h(xn-1,:) = h(xn-1,:) + h0(xn,:)*sxMat(xn,:)/2.0  ! right
 
uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
sxLong = reshape(sxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
syLong = reshape(syMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

! MAKE THE BAND
aBand = 0.0
do i = 1,(xn-2)*(yn-2)
	  
	aBand(i,2) = 1.0+sxLong(i)
	if (i-1 .gt. 0) then
	aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qx/2.0
	end if
	if (i+1 .le. (xn-2)*(yn-2)) then
	aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qx/2.0
	end if

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
	aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qx
	if (i .gt. 1) then
	aBand(i,1) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qx
	end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
	aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qx
	if (i .gt. 1) then
	aBand(i,3) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	aBand(i,1) = -sxLong(i)/2.0 + uLong(i)*qx
	end if
	end if
  
end do
  
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do
  
  
!!!!!!!!!!!! THIS !!!!!!!!!!!
h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))
sxMat(2:xn-1,2:yn-1) = reshape(sxLong, (/xn-2, yn-2/))
syMat(2:xn-1,2:yn-1) = reshape(syLong, (/xn-2, yn-2/))
sxLong = reshape(transpose(sxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
syLong = reshape(transpose(syMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! HORIZONTAL BOUNDARY CONDITIONS
h(:,2) = h(:,2) + flux(:,1)*syMat(:,1)/2.0 ! bottom
h(:,xn-1) = h(:,xn-1) + flux(:,2)*syMat(:,yn)/2.0 ! top

h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))


! MAKE THE BAND
bBand = 0.0
do i = 1,(xn-2)*(yn-2)
	bBand(i,2) = 1.0+syLong(i)
	if (i-1 .gt. 0) then
	bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qy/2.0
	end if
	if (i+1 .le. (xn-2)*(yn-2)) then
	bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qy/2.0
	end if

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
	bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qy
	if (i .gt. 1) then
	bBand(i,1) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qy
	end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
	bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qy
	if (i .gt. 1) then
	bBand(i,3) =  0.0
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qy
	end if
	end if
end do
  
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
! SUMMARY: computes the 2D streamfunction array of the current timestep
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
real(8) :: psi_next(xn,yn)
real(8) :: top_in(xn,1)
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
rhs1(xn-1,:) = rhs1(xn-1,:) 
rhs1(:,2) = rhs1(:,2) 
rhs1(:,yn-1) = rhs1(:,yn-1)
rhs1(:,yn-1) = rhs1(:,yn-1) +&
& top_in(:,1)/(4.0/(dy*dy*(permeability(:,yn)*rho_in(:,yn)))) !+ top_in(:,1)*(permx(:,yn))/(2.0*dx)

uVec = reshape(rhs1(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))

psi_next = 0.0


! MAKE THE BAND
aBand0 = 0.0
m = 2*(xn-2) + 1
do i = 1,(xn-2)*(yn-2)
	! DIAGONAL
	aBand0(i,(m+1)/2) = (2.0)/(permLong(i)*rhoLong(i)*dx*dx) + (2.0)/(permLong(i)*rhoLong(i)*dy*dy)
	! OFF-DIAGONALS
	if (i .gt. 1) then
	aBand0(i,((m+1)/2)-1) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx) + (permxLong(i))/(2.0*dx)
	end if
	if (i .lt. (xn-2)*(yn-2)) then
	aBand0(i,((m+1)/2)+1) = (-1.0)/(permLong(i)*rhoLong(i)*dx*dx) - (permxLong(i))/(2.0*dx)
	end if
	! MORE OFF-DIAGONALS
	if (i .le. (xn-2)*(yn-2)-(xn-2)) then
	aBand0(i,m) = (-1.0)/(permLong(i)*rhoLong(i)*dy*dy) - (permyLong(i))/(2.0*dy)
	end if
	if (i .ge. (xn-2)) then
	aBand0(i,1) = (-1.0)/(permLong(i)*rhoLong(i)*dy*dy) + (permyLong(i))/(2.0*dy)
	end if
end do
  
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
! ALT_NEXT
!
! SUMMARY: solves for equilibrium at a single grid cell using PHREEQC
!
! INPUTS: temp : temperature of grid cell
!         timestep : time elapsed
!         primaryList(5) : amounts of primary minerals
!         secondaryList(16) : amounts of secondary minerals
!         soluteList(11) : concentrations of solutes
!
! RETURNS: alt_next(1,altnum): returns everything from PHREEQC in a big pile
!          and it gets parsed in the main method's geochem loop
!
! ----------------------------------------------------------------------------------%%

function alt_next (temp, timestep, primaryList, secondaryList, soluteList)
use globals
use initialize
use alteration
implicit none

interface
end interface

! DECLARE YO SHIT
real(8) :: temp, timestep
real(8) :: alt_next(1,58)
real(8) :: alter0(1,58)
real(8) :: primaryList(5), secondaryList(16), soluteList(11)

! GRAB EVERYTHING FROM THE ALTERATION MODULE
alter0 = alter(temp-272.9, timestep, primaryList, secondaryList, soluteList)

! RENAME IT FOR A REASON THAT I FORGET
alt_next = alter0

end function alt_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
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
! SUMMARY : computes the darcy velocity (specific discharge) from the streamfunction
!           using finite difference partial derivatives
!
! INPUTS : psi(xn,yn) : 2D streamfunction array of current timestep
!
! RETURNS : velocities(xn,2*yn) : both u and v velocities in one matrix
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
! SUMMARY : versatile function for solving for second-order accurate partial 
!           derivatives of 1D or 2D arrays with respect to specified dimension
!           
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
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

