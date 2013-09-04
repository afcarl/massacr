
! ----------------------------------------------------------------------------------%%
!
! PHREEQOUT  - SINGLE BOX GEOCHEMICAL MODEL
! 
! SUMMARY: 
! COMPILE : gfortran -c -I/usr/local/include -L/usr/local/lib -liphreeqc phreeqout.f90
!			gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc phreeqout.o
!			./a.out
! 
!
! ----------------------------------------------------------------------------------%%

program phreeqout
    INCLUDE "IPhreeqc.f90.inc"
  INTEGER(KIND=4) :: id
  INTEGER(KIND=4) :: i
  CHARACTER(LEN=400) :: line
  character(len=800) :: inputz
  real(8), allocatable :: outmat(:,:)

  ! PHREEQC INPUT SCRIPT 
  inputz = "SOLUTION 1 Pure water" //NEW_LINE('')// &
  &"    pH 7.0" //NEW_LINE('')// &
  &"    temp 25.0" //NEW_LINE('')// &
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    Gypsum 0.0 1.0" //NEW_LINE('')// &
  &"    Anhydrite 0.0 1.0" //NEW_LINE('')// &
  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
  &"    25.0 75.0 in 51 steps" //NEW_LINE('')// &
!  &"SAVE solution 1" //NEW_LINE('')// &
!  &"SAVE equilibrium_phases 1" //NEW_LINE('')// &
  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -si anhydrite gypsum" //NEW_LINE('')// &
  &"END"
  write(*,*) inputz

  id = CreateIPhreeqc()
  IF (id.LT.0) THEN
     STOP
  END IF
  
  IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  IF (LoadDatabase(id, "phreeqc.dat").NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  IF (RunString(id, inputz).NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
  
  ! OPEN FILE
  OPEN(UNIT=12, FILE="outmat.txt", ACTION="write", STATUS="replace") 
  
  ! WRITE AWAY
  allocate(outmat(GetSelectedOutputStringLineCount(id),2))
  DO i=1,GetSelectedOutputStringLineCount(id)
     call GetSelectedOutputStringLine(id, i, line)
     ! HEADER BITS YOU MAY WANT
     if (i .eq. 1) then
     	write(12,*) trim(line)
     end if
     ! MEAT
     if (i .gt. 1) then
     	read(line,*) outmat(i,:)
     	write(12,"(2F10.5)") outmat(i,1), outmat(i,2)
     end if
  END DO
  
  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
     STOP
  END IF

  
  write(*,*) "phreequm dress"
  
  
    
end program phreeqout


