! littleChunk.f90 
! gfortran -O3 -g -I/usr/local/include -c alteration.f90

! gfortran -O3 -g -I/usr/local/include -c littleChunk.f90
! gfortran -O3 -g -I/usr/local/include -o littleChunk littleChunk.o alteration.o -L/usr/local/lib -liphreeqc
! ./littleChunk

! revising this script to use it in preliminary experiments

PROGRAM main
	use alteration
	
	implicit none
	
	interface
		! writes 2D array to file
		function write_matrix ( m, n, table, filename )
			implicit none
			integer :: m, n, j, output_status, unit0
			character ( len = * ) filename
			character ( len = 30 ) string
			real(4)  :: table(m,n) , write_matrix
		end function write_matrix

		! writes 1D array to file
		function write_vec ( n, vector, filename )
			implicit none
			integer :: n, j, output_status, unit0
			character ( len = * ) filename 
			real(4)  :: vector(n), write_vec
		end function write_vec
	end interface
	
	
	
	! inputs
	real(8) :: timestep, primary(5), secondary(28), solute(15), medium(7), temp
	real(8) :: solute0(15), water0
	!character(len=50) :: infile
	
	! other stuff
	integer :: i, j, steps
	real(8) ::  alt0(1,85) 
	real(8) :: out(100,85)
	real(8) :: yep, mix1= 0.9, mix2=0.1


	! command line arguments
	character(len=100) :: intemp
	character(len=100) :: infile
	integer :: in
	
	in = iargc()
	call getarg(1,intemp)
	call getarg(2,infile)
	
	read(intemp,*)temp
	print *,temp
	
	
	
	! initial conditions
	!infile = "prelim.txt"
	steps = 100
	timestep = 3.14e8
	!temp = 2.0

	
	primary(1) = 12.96 ! feldspar
	primary(2) = 6.96 ! augite
	primary(3) = 1.26 ! pigeonite
	primary(4) = .4 ! magnetite
	primary(5) = 96.77 ! basaltic glass
	
	secondary(:) = 0.0

! 	! columbia river
! 	solute(1) = 7.8 ! ph
! 	solute(2) = 8.451 ! pe
! 	solute(3) = 2.3e-3 ! Alk 1.6e-3
! 	solute(4) = 2.200e-3 !1.2e-2 ! H2CO3
! 	solute(5) = 6.0e-3 ! Ca
! 	solute(6) = 2.0e-5 ! Mg
! 	solute(7) = 1.0e-3 ! Na
! 	solute(8) = 1.0e-4 ! K
! 	solute(9) = 1.2e-6 ! Fe
! 	solute(10) = 1.0e-4 ! 1.0e-4 ! S(6)
! 	solute(11) = 2.0e-4 ! Si
! 	solute(12) = 3.0e-4 ! Cl
! 	solute(13) = 1.0e-6 ! Al
! 	solute(14) = 2.200e-3 ! HCO3-
! 	solute(15) = 0.0 ! CO3-2
	
	! today
	solute(1) = 7.8 ! ph
	solute(2) = 8.4 ! pe
	solute(3) = .0023 ! Alk 1.6e-3
	solute(4) = .0023 !1.2e-2 ! H2CO3
	solute(5) = .0105 ! Ca
	solute(6) = .0533 ! Mg
	solute(7) = .468 ! Na
	solute(8) = .00997 ! K
	solute(9) = 0.0 !1.2e-6 ! Fe
	solute(10) = .0281 ! 1.0e-4 ! S(6)
	solute(11) = 0.0 !2.0e-4 ! Si
	solute(12) = .546 ! Cl
	solute(13) = 0.0 !1.0e-6 ! Al
	solute(14) = .00234 ! HCO3-
	solute(15) = 0.0 ! CO3-2

! ! timestep grab
! solute(:) = 0.0
! solute(1) = 7.882959
! solute(2) = 14.66159
! solute(3) = 5.6893e-004
! solute(4) = 5.8361e-004
! solute(5) = 4.4860e-003
! solute(6) = 2.1836e-003
! solute(7) = 1.8791e-003
! solute(8) = 2.1209e-004
! solute(9) = 3.0177e-013
! solute(10) = 1.600e-003
! solute(11) = 6.061e-005
! solute(12) = 3.0036e-004
! solute(13) = 6.8355e-009
! solute(14) = 5.2320e-004
! solute(15) = 0.0 ! CO3-2

! ! hydrothermal solute concentrations [mol/kgw]
! solute(1) = 7.8 ! ph
! solute(2) = 8.451 ! pe
! solute(3) = 2.3e-3 ! Alk 1.6e-3
! solute(4) = 2.200e-3 !1.2e-2 ! H2CO3
! solute(5) = .0103 ! Ca
! solute(6) = .0528 ! Mg
! solute(7) = .469 ! Na
! solute(8) = .0102 ! K
! solute(9) = 0.0 ! Fe
! solute(10) = 0.0 !1e-6 ! 1.0e-4 ! S(6)
! solute(11) = 0.0 ! Si
! solute(12) = .546 ! Cl
! solute(13) = 0.0 ! Al
! solute(14) = 2.200e-3 ! HCO3-
! solute(15) = 0.0 ! CO3-2

solute0 = solute

medium(:) = 0.0
medium(3) = .385 ! water_volume
	
	write(*,*) "doing something..."
	
	do i=1,steps
		write(*,*) i
		
		water0 = medium(3)
		
		write(*,*) "altering"
		
		write(*,*) medium(3)
		
		alt0 = alter(temp,timestep,primary,secondary,solute,medium)
		!write(*,*) alt0
		!PARSING
		
		write(*,*) "altered"
		
		
		
		solute = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
		alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
		alt0(1,13), alt0(1,14), alt0(1,15), 0.0D+00/)
		
		

		secondary = (/ alt0(1,16), alt0(1,18), alt0(1,20), &
		alt0(1,22), alt0(1,24), alt0(1,26), alt0(1,28), alt0(1,30), alt0(1,32), alt0(1,34), &
		alt0(1,36), alt0(1,38), alt0(1,40), alt0(1,42), alt0(1,44), alt0(1,46), alt0(1,48), &
		alt0(1,50), alt0(1,52), alt0(1,54), alt0(1,56), alt0(1,58), alt0(1,60), alt0(1,62), &
		alt0(1,64), alt0(1,66), alt0(1,68), alt0(1,70)/)
		
	
		primary = (/ alt0(1,72), alt0(1,74), alt0(1,76), alt0(1,78), alt0(1,80)/)
		
		medium(3) = alt0(1,84)
		
		!solute = solute * water0 / medium(3)
		
		!solute(1) = solute(1)*mix1 + solute0(1)*mix2
		!solute(1) = solute0(1)
		!solute(2) = solute(2)*mix1 + solute0(2)*mix2
		solute(1) = -log10(mix1*10.0**(0.0-solute(1)) + mix2*10.0**(0.0-solute0(1)))
		!solute(2) = -log10(mix1*10.0**(-solute(2)) mix2*10.0**(-solute0(2)))
		solute(3) = solute(3)*mix1 + solute0(3)*mix2
		solute(4) = solute(4)*mix1 + solute0(4)*mix2
		solute(5) = solute(5)*mix1 + solute0(5)*mix2
		solute(6) = solute(6)*mix1 + solute0(6)*mix2
		solute(7) = solute(7)*mix1 + solute0(7)*mix2
		solute(8) = solute(8)*mix1 + solute0(8)*mix2
		solute(9) = solute(9)*mix1 + solute0(9)*mix2
		solute(10) = solute(10)*mix1 + solute0(10)*mix2
		solute(11) = solute(11)*mix1 + solute0(11)*mix2
		solute(12) = solute(12)*mix1 + solute0(12)*mix2
		solute(13) = solute(13)*mix1 + solute0(13)*mix2
		solute(14) = solute(14)*mix1 + solute0(14)*mix2
		solute(15) = solute(15)*mix1 + solute0(15)*mix2
		
		
		
		write(*,*) medium(3)
		
		alt0(1,2:16) = solute
		
		out(i,:) = alt0(1,:)
		
		!secondary = 0.0

	end do
	
	yep = write_matrix ( steps, 85, real(out,kind=4), infile )


	! ALL DONE!
	write(*,*) "this is phreeqing me out"
	
END PROGRAM main

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

  implicit none
  
  interface
	  function get_unit ( )
	    integer :: i, ios, get_unit
	    logical lopen
	  end function get_unit
  end interface
  
  
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

  implicit none
  
  interface
	  function get_unit ( )
	    integer :: i, ios, get_unit
	    logical lopen
	  end function get_unit
  end interface
  
  
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

! ----------------------------------------------------------------------------------%%
!
! GIVE ME A NUMBER I AM NOT USING SO I CAN USE IT TO MAKE MYSELF SOME FILES
!
! ----------------------------------------------------------------------------------%%

function get_unit ( )

  implicit none
  integer :: i, ios, get_unit
  logical lopen
  get_unit = 0
  do i = 1, 99
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
      inquire ( unit = i, opened = lopen, iostat = ios )
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          get_unit = i
          return
        end if
      end if
    end if
  end do
  return
end function get_unit