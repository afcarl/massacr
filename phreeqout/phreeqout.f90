
! ----------------------------------------------------------------------------------%%
!
! PHREEQOUT  - SINGLE BOX GEOCHEMICAL MODEL
! 
! SUMMARY: 
! COMPILE : gfortran -c -I/usr/local/include -L/usr/local/lib -liphreeqc phreeqoutTest.f90
!			gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc phreeqoutTest.o
!			./a.out
! 
!
! ----------------------------------------------------------------------------------%%

program phreeqout
    INCLUDE "IPhreeqc.f90.inc"
  INTEGER(KIND=4) :: id
  INTEGER(KIND=4) :: i
  CHARACTER(LEN=700) :: line
  character(len=2200) :: inputz
  character(len=3200) :: inputz0, inputSS
  real(8), allocatable :: outmat(:,:)

  !!!!!!!!!!!!!!!!!!!!!!
  !  OLD INPUT SCRIPT  !
  !!!!!!!!!!!!!!!!!!!!!!
  
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  SOLID SOLUTIONS ONLY   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  inputSS = "SOLUTION 1 # 1 mmol/l NaCl" //NEW_LINE('')// &
  &"    units   mmol/kgw" //NEW_LINE('')// &
  &"    pH 9.5" //NEW_LINE('')// &
  &"    pe 8.451" //NEW_LINE('')// &
  &"    density 1.023" //NEW_LINE('')// &
  &"    temp 60.0" //NEW_LINE('')// &
  &"    Ca 6.0e-1" //NEW_LINE('')// &
  &"    Mg 2.0e-0" //NEW_LINE('')// &
  &"    Na 1.0e-0" //NEW_LINE('')// &
  &"    K 1.0e-1" //NEW_LINE('')// &
  &"    Fe 1.2e-3" //NEW_LINE('')// &
  &"    S(6) 1.0e-1" //NEW_LINE('')// &
  &"    Si 2.0e-1" //NEW_LINE('')// &
  &"    Cl 3.0e-1" //NEW_LINE('')// &
  &"    Al 1.0e-3" //NEW_LINE('')// &
  &"    F 1.0e-3" //NEW_LINE('')// &
!  &"    Alkalinity 2.0e-3 as HCO3-" //NEW_LINE('')// &
!  &"    units   kg/kgw" //NEW_LINE('')// &
!  &"    CO2 600" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  &"CALCULATE_VALUES" //NEW_LINE('')// &
  
!  &"R(Plag)" //NEW_LINE('')// &
!  &"-start" //NEW_LINE('')// &
!  &"10 something = log10(((ACT('Ca+2')^.5)*(ACT('H2O')^3.0)*(ACT('SiO2')^2.5)*" // &
!  & "(ACT('Al+3')^1.5)*(ACT('Na+')^.5))/(ACT('H+')^6.0))" //NEW_LINE('')// &
!  &"100 SAVE something" //NEW_LINE('')// &
!  &"-end" //NEW_LINE('')// &
!  
!  &"R(Aug)" //NEW_LINE('')// &
!  &"-start" //NEW_LINE('')// &
!  &"10 something = log10(((ACT('Ca+2')^.7)*(ACT('H2O')^2.0)*(ACT('SiO2')^2.0)*" //&
!  &"(ACT('Fe+2')^.6)*(ACT('Mg+2')^.7))/(ACT('H+')^4.0))" //NEW_LINE('')// &
!  &"100 SAVE something" //NEW_LINE('')// &
!  &"-end" //NEW_LINE('')// &
!  
!  &"R(Pig)" //NEW_LINE('')// &
!  &"-start" //NEW_LINE('')// &
!  &"10 something = log10(((ACT('Ca+2')^1.14)*(ACT('H2O')^2.0)*(ACT('SiO2')^2.0)*" //&
!  &"(ACT('Fe+2')^.64)*(ACT('Mg+2')^.22))/(ACT('H+')^4.0))" //NEW_LINE('')// &
!  &"100 SAVE something" //NEW_LINE('')// &
!  &"-end" //NEW_LINE('')// &
  
  &"R(Glass)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 something = log10(((ACT('Ca+2')^.015)*(ACT('H2O')^.35)*(ACT('SiO2')^0.5)*" //&
  &"(ACT('Fe+2')^.095)*(ACT('Mg+2')^.065)*(ACT('Na+')^.025)*(ACT('K+')^.001)*" // &
  &"(ACT('Al+3')^.105))/(ACT('H+')^7.0))" //NEW_LINE('')// &
  &"100 SAVE something" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"END" //NEW_LINE('')// &
  
  &"SOLID_SOLUTIONS 1" //NEW_LINE('')// &
!  &"Plagioclase" //NEW_LINE('')// &
!  &"    -comp albite 0.5" //NEW_LINE('')// &
!  &"    -comp anorthite 0.5" //NEW_LINE('')// &
!  &"Augite" //NEW_LINE('')// &
!  &"    -comp enstatite 0.3" //NEW_LINE('')// &
!  &"    -comp ferrosilite 0.3" //NEW_LINE('')// &
!  &"    -comp wollastonite 0.35" //NEW_LINE('')// &
  &"Pigeonite" //NEW_LINE('')// &
  &"    -comp enstatite 0.57" //NEW_LINE('')// &
  &"    -comp ferrosilite 0.32" //NEW_LINE('')// &
  &"    -comp wollastonite 0.11" //NEW_LINE('')// &
!  &"BasaltGlass" //NEW_LINE('')// &
!  &"    -comp Ca 0.015" //NEW_LINE('')// &
!  &"    -comp Fe 0.095" //NEW_LINE('')// &
!  &"    -comp Mg 0.065" //NEW_LINE('')// &
!  &"    -comp Na 0.025" //NEW_LINE('')// &
!  &"    -comp K 0.01" //NEW_LINE('')// &
!  &"    -comp Al 0.105" //NEW_LINE('')// &
!!  &"    -comp S(6) 0.003" //NEW_LINE('')// &
!  &"    -comp SiO2(am) 0.5" //NEW_LINE('')// &
!!  &"    -comp SiO2(am) 0.5" //NEW_LINE('')// &
!!  &"    -comp O(s) .35" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  
  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
  &"    1.0 100.0 in 51 steps" //NEW_LINE('')// &
  
  &"Use solution 1" //NEW_LINE('')// &
  &"Use solid_solutions 1" //NEW_LINE('')// &
  
  
  
  
  
  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"-reset false" //NEW_LINE('')// &
  &"-ca R(Pig)" //NEW_LINE('')// &
  &"-temp" //NEW_LINE('')// &
  &"END"
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  BASALT DISSOLUTION EXPERIMENTS   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
  &"    pH 7.5" //NEW_LINE('')// &
  &"    units   mol/kgw" //NEW_LINE('')// &
!  &"    pressure 100" //NEW_LINE('')// &
!  &"    density 1.023" //NEW_LINE('')// &
  &"    temp 40.0" //NEW_LINE('')// &
  &"    Ca 6.0e-4" //NEW_LINE('')// &
  &"    Mg 2.0e-5" //NEW_LINE('')// &
  &"    Na 1.0e-3" //NEW_LINE('')// &
  &"    K 1.0e-4" //NEW_LINE('')// &
  &"    Fe 1.2e-6" //NEW_LINE('')// &
  &"    S(6) 1.0e-4 as SO4" //NEW_LINE('')// &
  &"    Si 2.0e-4" //NEW_LINE('')// &
  &"    Cl 3.0e-4" //NEW_LINE('')// &
  &"    Al 1.0e-6" //NEW_LINE('')// &
!  &"    pe 10.68" //NEW_LINE('')// &
  &"    Alkalinity 2.0e-3 as HCO3" //NEW_LINE('')// &
!  &"    pressure 1" //NEW_LINE('')// &
!  &"    CO2(g) 2 # 100 atm" //NEW_LINE('')// &
  &"    -water		.01	# kg" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    CO2(g) 2 1000000" //NEW_LINE('')// &

  &"    Siderite 0.0 0.0" //NEW_LINE('')// &
  &"    Kaolinite 0.0 0.0" //NEW_LINE('')// &
  &"    Goethite 0.0 0.0" //NEW_LINE('')// &
!  &"    Dolomite 0.0 0.0" //NEW_LINE('')// &
  &"    Celadonite 0.0 0.0" //NEW_LINE('')// &
  &"    SiO2(am) 0.0 0.0" //NEW_LINE('')// &
  &"    Albite 0.0 0.0" //NEW_LINE('')// &
  &"    Calcite 0.0 0.0" //NEW_LINE('')// &
  &"    Hematite 0.0 0.0" //NEW_LINE('')// &
  &"    Smectite-high-Fe-Mg 0.0 0.0" //NEW_LINE('')// &
  &"    Saponite-Mg 0.0 0.0" //NEW_LINE('')// &
  &"    Stilbite 0.0 0.0" //NEW_LINE('')// &
!  &"    Dawsonite 0.0 0.0" //NEW_LINE('')// &
!  &"    Magnesite 0.0 0.0" //NEW_LINE('')// &
  &"    Clinoptilolite-Ca 0.0 0.0" //NEW_LINE('')// &
  &"    Pyrite 0.0 0.0" //NEW_LINE('')// &
!  &"    Quartz 0.0 0.0" //NEW_LINE('')// &
  
!  &"    BasaltGlass 0.0 0.0 precipitate only" //NEW_LINE('')// &
!  &"    Plagioclase 0.0 0.0 precipitate only" //NEW_LINE('')// &
!  &"    Augite 0.0 0.0 precipitate only" //NEW_LINE('')// &
!  &"    Pigeonite 0.0 0.0 precipitate only" //NEW_LINE('')// &
!  &"    Anorthite 0.0 0.0 precipitate_only" //NEW_LINE('')// &

  
  &"KINETICS" //NEW_LINE('')// &
  &"Plagioclase" //NEW_LINE('')// &
  &"-m0 12.96" //NEW_LINE('')// &
  &"Augite" //NEW_LINE('')// &
  &"-m0 6.96" //NEW_LINE('')// &
  &"Pigeonite" //NEW_LINE('')// &
  &"-m0 1.26" //NEW_LINE('')// &
  &"Magnetite" //NEW_LINE('')// &
  &"-m0 .4" //NEW_LINE('')// &
  &"BGlass" //NEW_LINE('')// &
  &"-f Ca+2 0.015 Fe+2 0.095 Mg+2 0.065 " //&
  & "Na+ 0.025 K+ 0.01 Al+3 0.105 Si 0.5 S 0.003 H2O 1.35" //NEW_LINE('')// &
  &"-m0 96.77" //NEW_LINE('')// &

  &"    -step 1.57e20 in 100" //NEW_LINE('')// &
  &"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &
  
    &"Use solution 1" //NEW_LINE('')// &
  !&"Use equilibrium_phases 1" //NEW_LINE('')// &
  
  !&"SAVE solution 1" //NEW_LINE('')// &
  
  &"RATES" //NEW_LINE('')// &
  
  &"BGlass" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0=M*46.5*(1.52*10^-5)*0.1*(1e-10)*exp(-25.5/(.008314*TK))" // &
  &"*(((ACT('H+')^3)/(ACT('Al+3')))^.33)" //NEW_LINE('')// &
!  &"    10 rate=M*125*(1.52e-2)*0.1*(1e-10)*exp(-25.5/(.008314*TK))" // &
!  &"*(((ACT('H+')^3)/(ACT('Al+3')))^.33)" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Plagioclase" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate = M*270*(1.52*10^-5)*0.1*(((1.58e-9)"//&
  &"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
  &"(4.78e-15)*exp(-59/(.008314*TK))*(ACT('OH-'))^-0.57))"//NEW_LINE('')//&
  &"    20 save rate * time"//NEW_LINE('')//&
  &"-end" //NEW_LINE('')// &
  
  &"Augite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = M*230*(1.52e-5)*0.1*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //NEW_LINE('')// & 
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Pigeonite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = M*236*(1.52e-5)*0.1*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Magnetite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = M*231*(1.52e-5)*0.1*(((2.57e-9)*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279) +" // &
  &"(1.66e-11)*exp(-18.6/(.008314*TK))))" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &


   
  &"DUMP" //NEW_LINE('')// &
  &"    -solution 1" //NEW_LINE('')// &
  &"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -k plagioclase augite pigeonite magnetite bglass" //NEW_LINE('')// &
  &"    -ph" //NEW_LINE('')// &
!  &"    -molalities HCO3-" //NEW_LINE('')// &
  &"    -p stilbite sio2(am) kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
  &"    -p pyrite hematite goethite dolomite Smectite-high-Fe-Mg Dawsonite" //NEW_LINE('')// &
  &"    -p magnesite Clinoptilolite-Ca siderite" //NEW_LINE('')// &
  
!  &"    -activities H+ Al+3 " //NEW_LINE('')// &
!  &"    -equilibrium_phases Siderite Kaolinite Dolomite SiO2(am) Goethite" //NEW_LINE('')// &
  &"    -time" //NEW_LINE('')// &
  &"END"
  


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

!  IF (LoadDatabase(id, "minteq.dat").NE.0) THEN
!     CALL OutputErrorString(id)
!     STOP
!  END IF
  
  IF (LoadDatabase(id, "llnl.dat").NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF


  
  ! RUN INPUT
  IF (RunString(id, inputz0).NE.0) THEN
     CALL OutputErrorString(id)
     STOP
  END IF
  
  write(*,*) "it works!"
  
  ! PRINT DUMP/OUTPUT
  DO i=1,GetOutputStringLineCount(id)
     call GetOutputStringLine(id, i, line)
     write(*,*) trim(line)
  END DO
  
  WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
  
  WRITE(*,*) GetSelectedOutputStringLineCount(id)
  
  
  ! OPEN FILE
  OPEN(UNIT=12, FILE="testMat.txt", ACTION="write", STATUS="replace") 
  
  ! WRITE AWAY
  allocate(outmat(GetSelectedOutputStringLineCount(id)+1,44))
  DO i=1,GetSelectedOutputStringLineCount(id)
     call GetSelectedOutputStringLine(id, i, line)
     ! HEADER BITS YOU MAY WANT
     if (i .eq. 1) then
     	!write(12,*) trim(line)
     	write(*,*) trim(line)
     end if
     ! MEAT
     if (i .gt. 1) then
     	read(line,*) outmat(i,:)
     	!write(12,"(4F12.5)") outmat(i,:)
     	write(12,*) outmat(i,:)
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


