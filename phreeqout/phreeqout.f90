
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
  CHARACTER(LEN=900) :: line
  character(len=10200) :: inputz
  character(len=10200) :: inputz0, inputSS
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
  &"    pH 9.5" //NEW_LINE('')// &
  &"    units   mol/kgw" //NEW_LINE('')// &
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
  &"    Alkalinity 2.0e-3 as HCO3" //NEW_LINE('')// &
  &"    -water		.34	# kg" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    CO2(g) 2 1000000" //NEW_LINE('')// &

  &"    Siderite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Kaolinite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Goethite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Dolomite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Celadonite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    SiO2(am) 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Albite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Calcite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Hematite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Smectite-high-Fe-Mg 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Saponite-Mg 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Stilbite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Dawsonite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Magnesite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Clinoptilolite-Ca 0.0 0.0 precipitate_only" //NEW_LINE('')// &
  &"    Pyrite 0.0 0.0 precipitate_only" //NEW_LINE('')// &
!  &"    Quartz 0.0 0.0 precipitate_only" //NEW_LINE('')// &

  &"CALCULATE_VALUES" //NEW_LINE('')// &
  
  &"R(sum)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 sum = EQUI('Stilbite')*832.2/2.15 + EQUI('SiO2(am)')*60.0/2.62" //&
  &"+ EQUI('Kaolinite')*258.2/2.6 + EQUI('Albite')*262.3/2.62" // &
  &"+ EQUI('Saponite-Mg')*385.537/2.4 + EQUI('Celadonite')*396.8/3.0" // &
  &"+ EQUI('Clinoptilolite-Ca')*1344.49/2.62 + EQUI('Pyrite')*120.0/4.84" // &
  &"+ EQUI('Hematite')*103.8/5.3 + EQUI('Goethite')*88.8/3.8" // &
  &"+ EQUI('Dolomite')*184.3/2.84 + EQUI('Smectite-high-Fe-Mg')*425.7/2.7" // &
  &"+ EQUI('Dawsonite')*144.0/2.42 + EQUI('Magnesite')*84.3/3.0" // &
  &"+ EQUI('Siderite')*115.8/3.96 + EQUI('Calcite')*100.0/2.71" // &
  &"+ KIN('Plagioclase')*270.0/2.68 + KIN('Augite')*230.0/3.4" // &
  &"+ KIN('Pigeonite')*239.6/3.38 + KIN('Magnetite')*231.0/5.15" // &
  &"+ KIN('BGlass')*46.5/2.92" // &
  &"" //NEW_LINE('')// &
  &"100 SAVE sum" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"R(phi)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 phi = 1.0-(CALC_VALUE('R(sum)')*.001/(CALC_VALUE('R(sum)')*.001+SOLN_VOL))" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE phi" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"R(water_volume)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 water_volume = SOLN_VOL" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE water_volume" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  
  &"R(rho_s)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 rho_s = 2.93e6" //NEW_LINE('')// &
!  &"10 rho_s = EQUI('Stilbite')*2.15 + EQUI('SiO2(am)')*2.62" //&
!  &"+ EQUI('Kaolinite')*2.6 + EQUI('Albite')*2.62" // &
!  &"+ EQUI('Saponite-Mg')*2.4 + EQUI('Celadonite')*3.0" // &
!  &"+ EQUI('Clinoptilolite-Ca')*2.62 + EQUI('Pyrite')*4.84" // &
!  &"+ EQUI('Hematite')*5.3 + EQUI('Goethite')*3.8" // &
!  &"+ EQUI('Dolomite')*2.84 + EQUI('Smectite-high-Fe-Mg')*2.7" // &
!  &"+ EQUI('Dawsonite')*2.42 + EQUI('Magnesite')*3.0" // &
!  &"+ EQUI('Siderite')*3.96 + EQUI('Calcite')*2.71" // &
!  &"+ KIN('Plagioclase')*2.68 + KIN('Augite')*3.4" // &
!  &"+ KIN('Pigeonite')*3.38 + KIN('Magnetite')*5.15" // &
!  &"+ KIN('BGlass')*2.92" //NEW_LINE('')// &
!  &"20 rho_s = rho_s/ (EQUI('Stilbite') + EQUI('SiO2(am)')" //&
!  &"+ EQUI('Kaolinite') + EQUI('Albite')" // &
!  &"+ EQUI('Saponite-Mg') + EQUI('Celadonite')" // &
!  &"+ EQUI('Clinoptilolite-Ca') + EQUI('Pyrite')" // &
!  &"+ EQUI('Hematite') + EQUI('Goethite')" // &
!  &"+ EQUI('Dolomite') + EQUI('Smectite-high-Fe-Mg')" // &
!  &"+ EQUI('Dawsonite') + EQUI('Magnesite')" // &
!  &"+ EQUI('Siderite') + EQUI('Calcite')" // &
!  &"+ KIN('Plagioclase') + KIN('Augite')" // &
!  &"+ KIN('Pigeonite') + KIN('Magnetite')" // &
!  &"+ KIN('BGlass'))" //NEW_LINE('')// &
  &"100 SAVE rho_s" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  
  &"R(s_sp)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 s_sp = CALC_VALUE('R(phi)')/(1.0-CALC_VALUE('R(phi)')) * 400.0 / CALC_VALUE('R(rho_s)')" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE s_sp" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &

  
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
  & "Na+ 0.025 K+ 0.01 Al+3 0.105 SiO2 0.5 SO4-- 0.003 O .23 H+.692" //NEW_LINE('')// &
  &"-m0 96.77" //NEW_LINE('')// &

  &"    -step 1.57e20 in 100" //NEW_LINE('')// &

  &"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &
  
    &"Use solution 1" //NEW_LINE('')// &
    
    
  &"RATES" //NEW_LINE('')// &
  
  &"BGlass" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0=M*46.5*(CALC_VALUE('R(s_sp)'))*0.1*(1e-10)*exp(-25.5/(.008314*TK))" // &
  &"*(((ACT('H+')^3)/(ACT('Al+3')))^.33)" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Plagioclase" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate = (1-SR('Plagioclase'))*M*270.0*(1.52*10^-5)*0.1*(((1.58e-9)"//&
  &"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
  &"(4.78e-15)*exp(-59.0/(.008314*TK))*(ACT('H+'))^-0.57))"//NEW_LINE('')//&
  &"    20 save rate * time"//NEW_LINE('')//&
  &"-end" //NEW_LINE('')// &
  
  &"Augite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Augite'))*M*230.0*(1.52e-5)*0.1*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //NEW_LINE('')// & 
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Pigeonite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Pigeonite'))*M*236.0*(1.52e-5)*0.1*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Magnetite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Magnetite'))*M*231.0*(1.52e-5)*0.1*(((2.57e-9)" // &
  &"*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279)+(1.66e-11)*exp(-18.6/(.008314*TK))))" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &


   
  &"DUMP" //NEW_LINE('')// &
  &"    -solution 1" //NEW_LINE('')// &
  &"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -k plagioclase augite pigeonite magnetite bglass" //NEW_LINE('')// &
  &"    -ph" //NEW_LINE('')// &
  &"    -molalities HCO3-" //NEW_LINE('')// &
  &"    -p stilbite sio2(am) kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
  &"    -p pyrite hematite goethite dolomite Smectite-high-Fe-Mg Dawsonite" //NEW_LINE('')// &
  &"    -p magnesite siderite calcite" //NEW_LINE('')// &
  &"    -calculate_values R(phi) R(s_sp) R(water_volume)" //NEW_LINE('')// &
!  &"    -activities H+ Al+3 " //NEW_LINE('')// &
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
     !write(*,*) trim(line)
  END DO
  
  WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
  
  WRITE(*,*) GetSelectedOutputStringLineCount(id)
  
  
  ! OPEN FILE
  OPEN(UNIT=12, FILE="testMat.txt", ACTION="write", STATUS="replace") 
  
  ! WRITE AWAY
  allocate(outmat(GetSelectedOutputStringLineCount(id)+1,48))
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
  write(*,*) outmat(2,46)
  
end program phreeqout


