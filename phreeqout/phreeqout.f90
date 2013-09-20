
! ----------------------------------------------------------------------------------%%
!
! PHREEQOUT  - SINGLE BOX GEOCHEMICAL MODEL
! 
! SUMMARY: 
! COMPILE : gfortran -c -I/usr/local/include -L/usr/local/lib -liphreeqc rates.f90 phreeqout.f90
!			gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc rates.o phreeqout.o
!			./a.out
! 
!
! ----------------------------------------------------------------------------------%%

program phreeqout
use rates
    INCLUDE "IPhreeqc.f90.inc"
  INTEGER(KIND=4) :: id
  INTEGER(KIND=4) :: i, k
  CHARACTER(LEN=400) :: line
  character(len=3200) :: inputz0
  real(8), allocatable :: outmat(:,:)
  real(8) :: actH, actOH
  real(8) :: plag, aug, pig, glass, mag
  real(8) :: d_plag, d_aug, d_pig, d_glass, d_mag
  real(8) :: si_plag, si_aug, si_pig, si_glass, si_mag
  real(8) :: siderite
  real(8) :: temp
  character(len=20) :: s_siderite
  character(len=20) :: s_plag, s_aug, s_pig, s_glass, s_mag
  character(len=20) :: s_temp, s_time

!---------------------------------------!
!     INITIAL AMOUNTS OF EVERYTHING     !
!---------------------------------------!

plag = 1.0
aug = 1.0
pig = 0.5
glass = 2.0
mag = 0.5
siderite = 0.0

temp = 20.0


do k=1,10

  write(s_siderite,'(F10.5)') siderite
  write(s_plag,'(F10.5)') plag
  write(s_aug,'(F10.5)') aug
  write(s_pig,'(F10.5)') pig
  write(s_glass,'(F10.5)') glass
  write(s_mag,'(F10.5)') mag
  
  write(s_temp,'(F10.5)') temp
  write(s_time,'(F10.5)') time
  
!---------------------------------------!
!        INPUT STRING GOES HERE         !
!---------------------------------------!

  inputz0 = "SOLUTION 1 # 1 mmol/l NaCl" //NEW_LINE('')// &
  &"    units   mol/kgw" //NEW_LINE('')// &
  &"    pH 7.5" //NEW_LINE('')// &
  &"    pe 8.451" //NEW_LINE('')// &
  &"    density 1.023" //NEW_LINE('')// &
  &"    temp 80.0" //NEW_LINE('')// &
! &"    redox O(0)/O(-2)" //NEW_LINE('')// &
  &"    Ca 6.0e-4" //NEW_LINE('')// &
  &"    Mg 2.0e-3" //NEW_LINE('')// &
  &"    Na 1.0e-3" //NEW_LINE('')// &
  &"    K 1.0e-4" //NEW_LINE('')// &
  &"    Fe 1.2e-6" //NEW_LINE('')// &
  &"    S(6) 1.0e-4" //NEW_LINE('')// &
  &"    Si 2.0e-4" //NEW_LINE('')// &
  &"    Cl 3.0e-4" //NEW_LINE('')// &
  &"    Al 1.0e-6" //NEW_LINE('')// &
  &"    Alkalinity 2.0e-3 as HCO3" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  &"Use solution 1" //NEW_LINE('')// &
  
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    Siderite 0.0" // trim(s_siderite) // NEW_LINE('')// &
  &"    Plagioclase 0.0" // trim(s_plag) //NEW_LINE('')// &
  &"    Augite 0.0" // trim(s_aug) //NEW_LINE('')// &
  &"    Pigeonite 0.0" // trim(s_pig) //NEW_LINE('')// &
  &"    BasaltGlass 0.0" // trim(s_glass) //NEW_LINE('')// &
  &"    Magnetite 0.0" // trim(s_mag) //NEW_LINE('')// &
  
  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
  & trim(s_temp) //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -e Plagioclase augite pigeonite basaltglass magnetite" //NEW_LINE('')// &
  &"    -si Plagioclase augite pigeonite basaltglass magnetite" //NEW_LINE('')// &
  &"    -temp" //NEW_LINE('')// &
  &"END"
  
!---------------------------------------!
!         IPHREECQ INTERFACING          !
!---------------------------------------!

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

  IF (LoadDatabase(id, "llnl.dat").NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF

  IF (RunString(id, inputz0).NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF

  DO i=1,GetOutputStringLineCount(id)
     call GetOutputStringLine(id, i, line)
  END DO
  
  DO i=1,GetSelectedOutputStringLineCount(id)
     call GetSelectedOutputStringLine(id, i, line)
     ! HEADER
     if (i .eq. 1) then
     	write(*,*) trim(line)
     end if
     ! MEAT
     if (i .gt. 1) then
     	read(line,*) temp, plag, d_plag, aug, d_aug, pig, d_pig, glass, d_glass, &
     	& mag, d_mag, si_plag, si_aug, si_pig, si_glass, si_mag
     	write(*,'(11F16.5)') temp, plag, d_plag, aug, d_aug, pig, d_pig, glass, d_glass, &
     	& mag, d_mag, si_plag, si_aug, si_pig, si_glass, si_mag
     end if
  END DO
  
  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
     STOP
  END IF
  
  
  ! ALL DONE!
  
end do


  write(*,*) "phreequm dress"
  
  
end program phreeqout


