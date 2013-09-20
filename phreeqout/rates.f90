module rates
save
  
  real(8) :: R = .008314

  contains
  
! ----------------------------------------------------------------------------------%%
!
! PLAGIOCLASE
! 
! ----------------------------------------------------------------------------------%%

function r_plag(temp,timestep,actH,actOH)

implicit none
real(8) :: temp, timestep, actH, actOH
real(8) :: omega, sa, r_plag

omega = 1.0 

r_plag = sa*(1.0-omega) *( (1.58e-9)*exp(-53.5/(R*temp)) * (actH**0.541) + &
						& (3.39e-12)*exp(-57.4/(R*temp)) + &
						& (4.78e-15)*exp(-59.0/(R*temp)) * (actOH**(-0.57)))
						
  return
  
end function r_plag

! ----------------------------------------------------------------------------------%%
!
! AUGITE
! 
! ----------------------------------------------------------------------------------%%

function r_aug(temp,timestep,actH,actOH)

implicit none
real(8) :: temp, timestep, actH, actOH
real(8) :: omega, sa, r_aug

omega = 1.0 

r_aug = sa*(1.0-omega) *( (1.58e-7)*exp(-78.0/(R*temp)) * (actH**0.7) + &
						& (1.07e-12)*exp(-78.0/(R*temp)))
						
  return
  
end function r_aug

! ----------------------------------------------------------------------------------%%
!
! PIGEONITE
! 
! ----------------------------------------------------------------------------------%%

function r_pig(temp,timestep,actH,actOH)

implicit none
real(8) :: temp, timestep, actH, actOH
real(8) :: omega, sa, r_pig

omega = 1.0 

r_pig = sa*(1.0-omega) *( (1.58e-7)*exp(-78.0/(R*temp)) * (actH**0.7) + &
						& (1.07e-12)*exp(-78.0/(R*temp)))
						
  return
  
end function r_pig

! ----------------------------------------------------------------------------------%%
!
! BASALTIC GLASS
! 
! ----------------------------------------------------------------------------------%%

function r_glass(temp,timestep,actH,actAL)

implicit none
real(8) :: temp, timestep, actH, actAL
real(8) :: omega, sa, r_glass

omega = 1.0 

r_glass = sa*(1.0-omega) *( (10.0**(-5.6))*exp(-25.5/(R*temp)) * &
						& ((actH**3.0)/actAL))
						
  return
  
end function r_glass


! ----------------------------------------------------------------------------------%%
!
! MAGNETITE
! 
! ----------------------------------------------------------------------------------%%

function r_mag(temp,timestep,actH,actOH)

implicit none
real(8) :: temp, timestep, actH, actOH
real(8) :: omega, sa, r_mag

omega = 1.0 

r_mag = sa*(1.0-omega) *( (2.57e-9)*exp(-18.6/(R*temp)) * (actH**0.279) + &
						& (1.66e-11)*exp(-18.6/(R*temp)))
						
  return
  
end function r_mag



end module rates


  