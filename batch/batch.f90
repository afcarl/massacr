
! ----------------------------------------------------------------------------------%%
!
! BATCH  - LIKE PHREEQOUT BUT LESS EXTREME
! 
! SUMMARY: 
! COMPILE : gfortran -c -I/usr/local/include -L/usr/local/lib -liphreeqc batch.f90
!			gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc batch.o
!			./a.out
! 
!
! ----------------------------------------------------------------------------------%%

program batch
    INCLUDE "IPhreeqc.f90.inc"
  INTEGER(KIND=4) :: id
  INTEGER(KIND=4) :: i
  CHARACTER(LEN=900) :: line
  character(len=10200) :: inputz
  character(len=10200) :: inputz0, inputSS
  real(8), allocatable :: outmat(:,:)
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !        BATCH INPUT SCRIPT         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

!   !water composition is grande ronde river
!   inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
!   &"    pH 7.5" //NEW_LINE('')// &
!   &"    units   mol/kgw" //NEW_LINE('')// &
!   &"    temp 20.0" //NEW_LINE('')// &
!   &"    Ca 6.0e-4" //NEW_LINE('')// &
!   &"    Mg 2.0e-5" //NEW_LINE('')// &
!   &"    Na 1.0e-3" //NEW_LINE('')// &
!   &"    K 1.0e-4" //NEW_LINE('')// &
!   &"    Fe 1.2e-6" //NEW_LINE('')// &
!   &"    S(6) 1.0e-4 as SO4" //NEW_LINE('')// &
!   &"    Si 2.0e-4" //NEW_LINE('')// &
!   &"    Cl 3.0e-4" //NEW_LINE('')// &
!   &"    Al 1.0e-6" //NEW_LINE('')// &
!   &"    Alkalinity 2.0e-3 as HCO3" //NEW_LINE('')// &
!   &"    -water		.38	# kg" //NEW_LINE('')// &
!   &"END" //NEW_LINE('')// &

 ! water composition is nordstrom et al. 1979 seawater
  inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
  &"    pH 8.22" //NEW_LINE('')// &
  &"    units   ppm" //NEW_LINE('')// &
  &"    temp 10.0" //NEW_LINE('')// &
  &"    Ca 412.3" //NEW_LINE('')// &
  &"    Mg 1291.8" //NEW_LINE('')// &
  &"    Na 10768.0" //NEW_LINE('')// &
  &"    K 399.1" //NEW_LINE('')// &
  &"    Fe .002" //NEW_LINE('')// &
  &"    Mn .0002" //NEW_LINE('')// &
  &"    S(6) 2712.0 as SO4" //NEW_LINE('')// &
  &"    Si 4.28" //NEW_LINE('')// &
  &"    Cl 19353.0" //NEW_LINE('')// &
  &"    NO3 0.29" //NEW_LINE('')// &
  &"    NH4 0.03" //NEW_LINE('')// &
  &"    Alkalinity 141.682 as HCO3" //NEW_LINE('')// &
  &"    -water		5.0	# kg" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
 
!     &"SOLUTION 2 " //NEW_LINE('')// &
!     &"    pH 8.22" //NEW_LINE('')// &
!     &"    units   ppm" //NEW_LINE('')// &
!     &"    temp 10.0" //NEW_LINE('')// &
!     &"    Ca 412.3" //NEW_LINE('')// &
!     &"    Mg 1291.8" //NEW_LINE('')// &
!     &"    Na 10768.0" //NEW_LINE('')// &
!     &"    K 399.1" //NEW_LINE('')// &
!     &"    Fe .002" //NEW_LINE('')// &
!     &"    Mn .0002" //NEW_LINE('')// &
!     &"    S(6) 2712.0 as SO4" //NEW_LINE('')// &
!     &"    Si 4.28" //NEW_LINE('')// &
!     &"    Cl 19353.0" //NEW_LINE('')// &
!     &"    NO3 0.29" //NEW_LINE('')// &
!     &"    NH4 0.03" //NEW_LINE('')// &
!     &"    Alkalinity 141.682 as HCO3" //NEW_LINE('')// &
!     &"    -water		.5	# kg" //NEW_LINE('')// &
!     &"END" //NEW_LINE('')// &

  
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
  &"    CO2(g) -3.45 24" //NEW_LINE('')// &
!  &"    O2(g) -.69897 10000" //NEW_LINE('')// &
!  &"    Siderite 0.0 0.0" //NEW_LINE('')// &
  &"    Kaolinite 0.0 0.0" //NEW_LINE('')// &
  &"    Goethite 0.0 0.0" //NEW_LINE('')// &
!  &"    Dolomite 0.0 0.0" //NEW_LINE('')// &
  &"    Celadonite 0.0 0.0" //NEW_LINE('')// &
  &"    SiO2(am) 0.0 0.0" //NEW_LINE('')// &
  &"    Albite 0.0 0.0" //NEW_LINE('')// &
  &"    Calcite 0.0 0.0" //NEW_LINE('')// &
  &"    Hematite 0.0 0.0" //NEW_LINE('')// &
!  &"    Smectite-high-Fe-Mg 0.0 0.0" //NEW_LINE('')// &
  &"    Saponite-Mg 0.0 0.0" //NEW_LINE('')// &
!  &"    Stilbite 0.0 0.0" //NEW_LINE('')// &
!  &"    Dawsonite 0.0 0.0" //NEW_LINE('')// &
!  &"    Magnesite 0.0 0.0" //NEW_LINE('')// &
  &"    Clinoptilolite-Ca 0.0 0.0" //NEW_LINE('')// &
  &"    Pyrite 0.0 0.0" //NEW_LINE('')// &
  &"    Quartz 0.0 0.0" //NEW_LINE('')// &
  &"    K-Feldspar 0.0 0.0" //NEW_LINE('')// &
  
!   &"EQUILIBRIUM_PHASES 2" //NEW_LINE('')// &
!   &"    CO2(g) -3.45 10000" //NEW_LINE('')// &
!   &"    O2(g) -.69897 10000" //NEW_LINE('')// &

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
  &"+ EQUI('Quartz')*60.0/2.62 + EQUI('K-Feldspar')*193.0/2.56" // &
  &"+ KIN('Plagioclase')*270.0/2.68 + KIN('Augite')*230.0/3.4" // &
  &"+ KIN('Pigeonite')*239.6/3.38 + KIN('Magnetite')*231.0/5.15" // &
  &"+ KIN('BGlass')*46.5/2.92" // &
  &"" //NEW_LINE('')// &
  &"100 SAVE sum" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  ! .7-.3*SIM_TIME/3.14e11
  &"R(phi)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 phi = 1.0-(CALC_VALUE('R(sum)')/(CALC_VALUE('R(sum)')+(TOT('water')*1000.0)))" //&
  !&"10 phi = 0.1" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE phi" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"R(water_volume)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 water_volume = TOT('water')" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE water_volume" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  

  &"R(rho_s)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  !&"10 rho_s = CALC_VALUE('R(sum)')" //NEW_LINE('')// &
&"10 rho_s = EQUI('Stilbite')*2.15 + EQUI('SiO2(am)')*2.62" //&
  &"+ EQUI('Kaolinite')*2.6 + EQUI('Albite')*2.62" // &
  &"+ EQUI('Saponite-Mg')*2.4 + EQUI('Celadonite')*3.0" // &
  &"+ EQUI('Clinoptilolite-Ca')*2.62 + EQUI('Pyrite')*4.84" // &
  &"+ EQUI('Hematite')*5.3 + EQUI('Goethite')*3.8" // &
  &"+ EQUI('Dolomite')*2.84 + EQUI('Smectite-high-Fe-Mg')*2.7" // &
  &"+ EQUI('Dawsonite')*2.42 + EQUI('Magnesite')*3.0" // &
  &"+ EQUI('Siderite')*3.96 + EQUI('Calcite')*2.71" // &
  &"+ EQUI('Quartz')*2.62 + EQUI('k-Feldspar')*2.56" // &
  &"+ KIN('Plagioclase')*2.68 + KIN('Augite')*3.4" // &
  &"+ KIN('Pigeonite')*3.38 + KIN('Magnetite')*5.15" // &
  &"+ KIN('BGlass')*2.92" //NEW_LINE('')// &
  &"20 rho_s = rho_s/ (EQUI('Stilbite') + EQUI('SiO2(am)')" //&
  &"+ EQUI('Kaolinite') + EQUI('Albite')" // &
  &"+ EQUI('Saponite-Mg') + EQUI('Celadonite')" // &
  &"+ EQUI('Clinoptilolite-Ca') + EQUI('Pyrite')" // &
  &"+ EQUI('Hematite') + EQUI('Goethite')" // &
  &"+ EQUI('Dolomite') + EQUI('Smectite-high-Fe-Mg')" // &
  &"+ EQUI('Dawsonite') + EQUI('Magnesite')" // &
  &"+ EQUI('Siderite') + EQUI('Calcite')" // &
  &"+ EQUI('Quartz') + EQUI('K-Feldspar')" // &
  &"+ KIN('Plagioclase') + KIN('Augite')" // &
  &"+ KIN('Pigeonite') + KIN('Magnetite')" // &
  &"+ KIN('BGlass'))" //NEW_LINE('')// &
!  &"20 rho_s = rho_s / (EQUI('Stilbite')*832.2 + EQUI('SiO2(am)')*60.0" //&
!  &"+ EQUI('Kaolinite')*258.2 + EQUI('Albite')*262.3" // &
!  &"+ EQUI('Saponite-Mg')*385.537 + EQUI('Celadonite')*396.8" // &
!  &"+ EQUI('Clinoptilolite-Ca')*1344.49 + EQUI('Pyrite')*120.0" // &
!  &"+ EQUI('Hematite')*103.8 + EQUI('Goethite')*88.8" // &
!  &"+ EQUI('Dolomite')*184.3 + EQUI('Smectite-high-Fe-Mg')*425.7" // &
!  &"+ EQUI('Dawsonite')*144.0 + EQUI('Magnesite')*84.3" // &
!  &"+ EQUI('Siderite')*115.8 + EQUI('Calcite')*100.0" // &
!  &"+ KIN('Plagioclase')*270.0 + KIN('Augite')*230.0" // &
!  &"+ KIN('Pigeonite')*239.6 + KIN('Magnetite')*231.0" // &
!  &"+ KIN('BGlass')*46.5)" //NEW_LINE('')// &
  &"30 rho_s = rho_s * 1000000.0" //NEW_LINE('')// &
  &"100 SAVE rho_s" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  
  &"R(s_sp)" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"10 s_sp = (CALC_VALUE('R(phi)')/(1.0-CALC_VALUE('R(phi)')))*400.0/CALC_VALUE('R(rho_s)')" //&
  &"" //NEW_LINE('')// &
  &"100 SAVE s_sp" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
!   &"TRANSPORT" //NEW_LINE('')// &
!   &"-cells 2" //NEW_LINE('')// &
!   &"-lengths .1" //NEW_LINE('')// &
!   &"-shifts 1" //NEW_LINE('')// &
!   &"-time_step 3.14e7" //NEW_LINE('')// &
!   &"-flow_direction diffusion_only" //NEW_LINE('')// &
!   &"-diffusion_coefficient 7.0e-2" //NEW_LINE('')// &
!   !&"-stagnant 1" //NEW_LINE('')// &
!   &"END" //NEW_LINE('')// &
!   &"" //NEW_LINE('')// &
!   &"" //NEW_LINE('')// &
  
  &"KINETICS" //NEW_LINE('')// &
  &"Plagioclase" //NEW_LINE('')// &
  &"-m0 129.6" //NEW_LINE('')// &
  &"Augite" //NEW_LINE('')// &
  &"-m0 69.6" //NEW_LINE('')// &
  &"Pigeonite" //NEW_LINE('')// &
  &"-m0 12.6" //NEW_LINE('')// &
  &"Magnetite" //NEW_LINE('')// &
  &"-m0 4.0" //NEW_LINE('')// &
  &"BGlass" //NEW_LINE('')// &
  &"-f Ca 0.015 Fe 0.095 Mg 0.065 " //&
  & "Na 0.025 K 0.01 Al 0.105 Si 0.5 S 0.003 O 1.35" //NEW_LINE('')// &
  &"-m0 967.7" //NEW_LINE('')// &
  !&"-m0 0.0" //NEW_LINE('')// &
!  &"pCO2" //NEW_LINE('')// &
!  &"-f CO2" //&
!  &"-m0 10.0" //NEW_LINE('')// &

  &"    -step 3.14e13 in 100" //NEW_LINE('')// &
  
  &"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &
  
    &"Use solution 1" //NEW_LINE('')// &
    
  &"RATES" //NEW_LINE('')// &
  
  &"BGlass" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0=M*46.5*CALC_VALUE('R(s_sp)')*1.0*(1e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
  &"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Plagioclase" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate = (1-SR('Plagioclase'))*M*270.0*CALC_VALUE('R(s_sp)')*1.0*((1.58e-9)"//&
  &"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
  &"(4.78e-15)*exp(-59.0/(.008314*TK))*(ACT('H+')^-0.57))"//NEW_LINE('')//&
  &"    20 save rate * time"//NEW_LINE('')//&
  &"-end" //NEW_LINE('')// &
  
  &"Augite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Augite'))*M*230.0*CALC_VALUE('R(s_sp)')*1.0*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //NEW_LINE('')// & 
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Pigeonite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Pigeonite'))*M*236.0*CALC_VALUE('R(s_sp)')*1.0*(((1.58e-7)" // &
  &"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &
  
  &"Magnetite" //NEW_LINE('')// &
  &"-start" //NEW_LINE('')// &
  &"    10 rate0 = (1-SR('Magnetite'))*M*231.0*CALC_VALUE('R(s_sp)')*1.0*(((2.57e-9)" // &
  &"*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279)+(1.66e-11)*exp(-18.6/(.008314*TK))))" //NEW_LINE('')// &
  &"    20 save rate0 * time" //NEW_LINE('')// &
  &"-end" //NEW_LINE('')// &

!   &"pCO2" //NEW_LINE('')// &
!   &"-start" //NEW_LINE('')// &
!   &"    10 rate = 1.0e-13" //NEW_LINE('')// &
!   &"    20 save rate * time" //NEW_LINE('')// &
!   &"-end" //NEW_LINE('')// &


   
  &"DUMP" //NEW_LINE('')// &
  &"    -solution 1" //NEW_LINE('')// &
  &"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
  &"    -k plagioclase augite pigeonite magnetite bglass" //NEW_LINE('')// &
  &"    -ph" //NEW_LINE('')// &
  &"    -alkalinity" //NEW_LINE('')// &
!  &"    -molalities HCO3-" //NEW_LINE('')// &
  &"    -p stilbite sio2(am) kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
  &"    -p pyrite hematite goethite dolomite Smectite-high-Fe-Mg Dawsonite" //NEW_LINE('')// &
  &"    -p magnesite siderite calcite quartz k-feldspar" //NEW_LINE('')// &
  &"    -calculate_values R(phi) R(s_sp) R(water_volume) R(rho_s)" //NEW_LINE('')// &
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
  allocate(outmat(GetSelectedOutputStringLineCount(id)+1,53))
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
     	!write(*,*) trim(line)
     end if
  END DO
  
  ! DESTROY INSTANCE
  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
     STOP
  END IF

  ! ALL DONE!
  write(*,*) "phreequm dress"
  !write(*,*) outmat(2,46)
  
end program batch


