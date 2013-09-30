
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
  real(8) :: actH, actOH, actAL
  real(8) :: plag, aug, pig, glass, mag
  real(8) :: d_plag, d_aug, d_pig, d_glass, d_mag
  real(8) :: si_plag, si_aug, si_pig, si_glass, si_mag
  real(8) :: m_ca, m_mg, m_fe, m_s6, m_si, m_na, m_cl, m_al, m_alk, m_k
  real(8) :: siderite
  real(8) :: temp
  character(len=20) :: s_siderite
  character(len=20) :: s_plag, s_aug, s_pig, s_glass, s_mag
  character(len=20) :: s_temp, s_time
  character(len=20) :: s_ca, s_mg, s_fe, s_s6, s_si, s_na, s_cl, s_al, s_k

!---------------------------------------!
!     INITIAL AMOUNTS OF EVERYTHING     !
!---------------------------------------!

plag = 1.0
aug = 1.0
pig = 0.5
glass = 2.0
mag = 0.5
siderite = 1.0

temp = 20.0

m_ca = 6.0e-4
m_mg = 2.0e-3
m_fe = 1.26e-6
m_s6 = 1.0e-4
m_si = 2.0e-4
m_na = 1.0e-3
m_cl = 3.0e-4
m_al = 1.0e-6
!m_alk = 2.0e-3
m_k = 1.0e-4


do k=1,10

  write(s_siderite,'(F10.5)') siderite
  write(s_plag,'(F10.5)') plag
  write(s_aug,'(F10.5)') aug
  write(s_pig,'(F10.5)') pig
  write(s_glass,'(F10.5)') glass
  write(s_mag,'(F10.5)') mag
  
  write(s_temp,'(F10.5)') temp
  write(s_time,'(F10.5)') time
  

  write(s_ca,'(F10.5)') m_ca
  write(s_mg,'(F10.5)') m_mg
  write(s_fe,'(F10.5)') m_fe
  write(s_s6,'(F10.5)')  m_s6
  write(s_si,'(F10.5)')  m_si
  write(s_na,'(F10.5)')  m_na
  write(s_cl,'(F10.5)')  m_cl
  write(s_al,'(F10.5)')  m_al
  !write(s_alk,'(F10.5)')  m_alk
  write(s_k,'(F10.5)')  m_k
  
!---------------------------------------!
!        INPUT STRING GOES HERE         !
!---------------------------------------!

  inputz0 = "SOLUTION 1 # 1 mmol/l NaCl" //NEW_LINE('')// &
  &"    units   mol/kgw" //NEW_LINE('')// &
  &"    pH 7.5" //NEW_LINE('')// &
  &"    pe 8.451" //NEW_LINE('')// &
  &"    density 1.023" //NEW_LINE('')// &
  &"    temp" // trim(s_temp) //NEW_LINE('')// &
! &"    redox O(0)/O(-2)" //NEW_LINE('')// &
  &"    Ca" // trim(s_ca) //NEW_LINE('')// &
  &"    Mg" // trim(s_mg) //NEW_LINE('')// &
  &"    Na" // trim(s_na) //NEW_LINE('')// &
  &"    K" // trim(s_k) //NEW_LINE('')// &
  &"    Fe"// trim(s_fe) //NEW_LINE('')// &
  &"    S(6)" // trim(s_s6) //NEW_LINE('')// &
  &"    Si" // trim(s_si) //NEW_LINE('')// &
  &"    Cl" // trim(s_cl) //NEW_LINE('')// &
  &"    Al" // trim(s_al) //NEW_LINE('')// &
  &"    Alkalinity 2.0e-3 as HCO3" //NEW_LINE('')// &
  &"END" //NEW_LINE('')// &
  
  &"Use solution 1" //NEW_LINE('')// &
  
  &"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
!  &"    Siderite 0.0" // trim(s_siderite) // NEW_LINE('')// &
  
  
!  &"REACTION_TEMPERATURE 1" //NEW_LINE('')// &
!  & trim(s_temp) //NEW_LINE('')// &
  
  &"Save solution 1" //NEW_LINE('')// &
  
 &"DUMP" //NEW_LINE('')// &
  &"    -solution 1" //NEW_LINE('')// &
  &"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
!  &"    -e plagioclase augite pigeonite basaltglass magnetite" //NEW_LINE('')// &
  &"    -activities H+ OH- Al+3" //NEW_LINE('')// &
  &"    -m Ca+2 Mg+2 Fe+2 SO4-2 SiO2 Na+ Cl- Al+3 K+" //NEW_LINE('')// &
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
     !write(*,*) trim(line)
  END DO
  
  DO i=1,GetSelectedOutputStringLineCount(id)
     call GetSelectedOutputStringLine(id, i, line)
     ! HEADER
     if (i .eq. 1 .and. k .eq. 1) then
     	write(*,*) trim(line)
     end if
     ! MEAT
     if (i .gt. 1) then
     	read(line,*) temp, m_ca, m_mg, m_fe, m_s6, m_si, m_na, m_cl, m_al, m_k,&
     	&actH, actOH, actAL
!     	&plag, d_plag, aug, d_aug, pig, d_pig, glass, d_glass, mag, d_mag,&
!     	&si_plag, si_aug, si_pig, si_glass, si_mag
       write(*,'(19F16.12)') temp, m_ca, m_mg, m_fe, m_s6, m_si, m_na, m_cl, m_al, m_k,&
       & actH, actOH, actAL
     end if
  END DO
  
  IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
     STOP
  END IF
  
  !write(*,*) plag, aug, pig, glass, mag

  if (si_plag .lt. 0.0) then
  	
  end if
  
  
!---------------------------------------!
!    DISSOLVE BASALT MINERALS/GLASS     !
!---------------------------------------!
! get new amount of mineral: mineral = mineral - rate*time
! get new constituents into the seawater: m_i
  
end do


  write(*,*) "phreequm dress"
  
  
end program phreeqout


