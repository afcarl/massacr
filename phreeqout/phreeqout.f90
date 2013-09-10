
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
  character(len=1200) :: inputz
  character(len=1200) :: inputz0
  real(8), allocatable :: outmat(:,:)

!!! INPUT INFO
! EQUILIBRIUM_PHASES 1
! phase targetSatIndex amount
! ---default takes a single # as SI and 10.0 mols
! ---SI = 0 means in equilibrium with fluid

!!! END INFO BLOCK

  ! PHREEQC INPUT SCRIPT 
  inputz = "SOLUTION 1 Pure water" //NEW_LINE('')// &
  &"    pH 7.0" //NEW_LINE('')// &
  &"    temp 25.0" //NEW_LINE('')// &
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    Gypsum 0.0 1.0" //NEW_LINE('')// &
  &"    dolomite 0.0 1.0" //NEW_LINE('')// &
  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
  &"    0.0 200.0 in 51 steps" //NEW_LINE('')// &
  &"SAVE solution 1" //NEW_LINE('')// &
  &"SAVE equilibrium_phases 1" //NEW_LINE('')// &
  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -equilibrium_phases Gypsum dolomite" //NEW_LINE('')// &
  &"    -temp" //NEW_LINE('')// &

  &"END"
  
  inputz0 = "SOLUTION 1 # 1 mmol/l NaCl" //NEW_LINE('')// &
  &"    units   ppm" //NEW_LINE('')// &
  &"    pH 8.22" //NEW_LINE('')// &
  &"    pe 8.451" //NEW_LINE('')// &
  &"    density 1.023" //NEW_LINE('')// &
  &"    temp 25.0" //NEW_LINE('')// &
  &"    redox O(0)/O(-2)" //NEW_LINE('')// &
!  &"    Ca 412.3" //NEW_LINE('')// &
!  &"    Mg 1291.8" //NEW_LINE('')// &
  &"    Ca 500.0" //NEW_LINE('')// &
  &"    Mg 500.0" //NEW_LINE('')// &
  &"    Na 10768.0" //NEW_LINE('')// &
  &"    K 399.1" //NEW_LINE('')// &
  &"    Fe .002" //NEW_LINE('')// &
  &"    Mn .0002 pe" //NEW_LINE('')// &
  &"    Si 4.28" //NEW_LINE('')// &
  &"    Cl 19353.0" //NEW_LINE('')// &
  &"    Alkalinity 141.682 as HCO3" //NEW_LINE('')// &
  &"    S(6) 2712.0" //NEW_LINE('')// &
  &"    N(5) 0.29 gfw 62.0" //NEW_LINE('')// &
  &"    N(-3) 0.03 as NH4" //NEW_LINE('')// &
  &"    O(0) 1.0 O2(g) -0.7" //NEW_LINE('')// &
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    quartz 0.0 1.0" //NEW_LINE('')// &
!  &"    diopside 0.0 1.0" //NEW_LINE('')// &
  &"    Wollastonite 0.0 1.0" //NEW_LINE('')// &
  &"    P-Wollstanite 0.0 1.0" //NEW_LINE('')// &
  &"    Ca-Olivine 0.0 1.0" //NEW_LINE('')// &
  &"    Magnetite 0.0 1.0" //NEW_LINE('')// &
  &"    Anorthite 0.0 1.0" //NEW_LINE('')// &
!  &"    Albite(low) 0.0 0.5" //NEW_LINE('')// &
!  &"    k-feldspar 0.0 1.0" //NEW_LINE('')// &
!  &"    dolomite 0.0 0.0" //NEW_LINE('')// &
!  &"    calcite 0.0 0.0" //NEW_LINE('')// &
!  &"    kaolinite 0.0 0.0" //NEW_LINE('')// &
  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
  &"    0.0 750.0 in 51 steps" //NEW_LINE('')// &
  &"SAVE solution 1" //NEW_LINE('')// &
  &"SAVE equilibrium_phases 1" //NEW_LINE('')// &
  &"DUMP" //NEW_LINE('')// &
  &"    -solution 1" //NEW_LINE('')// &
  &"    -equilibrium_phases" //NEW_LINE('')// &
  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -molalities Ca+2 Mg+2" //NEW_LINE('')// &
  &"    -temp" //NEW_LINE('')// &
  &"END"
  
  write(*,*) inputz0

  ! INITIALIZE STUFF
  id = CreateIPhreeqc()
  IF (id.LT.0) THEN
     STOP
  END IF
  
  IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
     CALL OutputErrorString(id)
     STOP
  END IF

  IF (LoadDatabase(id, "minteq.dat").NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  ! RUN INPUT
  IF (RunString(id, inputz0).NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  ! PRINT DUMP/OUTPUT
  DO i=1,GetOutputStringLineCount(id)
     call GetOutputStringLine(id, i, line)
     write(*,*) trim(line)
  END DO
  
  WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
  
  ! OPEN FILE
  OPEN(UNIT=12, FILE="outmat.txt", ACTION="write", STATUS="replace") 
  
  ! WRITE AWAY
  allocate(outmat(GetSelectedOutputStringLineCount(id),3))
  DO i=1,GetSelectedOutputStringLineCount(id)
     call GetSelectedOutputStringLine(id, i, line)
     ! HEADER BITS YOU MAY WANT
     if (i .eq. 1) then
     	!write(12,*) trim(line)
     	write(*,*) trim(line)
     end if
     ! MEAT
     if (i .gt. 2) then
     	read(line,*) outmat(i,:)
     	write(12,"(3F12.5)") outmat(i,:)
     	write(*,*) trim(line)
     end if
  END DO
  
  ! DESTROY INSTANCE
  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
     STOP
  END IF

  ! ALL DONE!
  write(*,*) "phreequm dress"
  
  
end program phreeqout


