module Constants  
implicit none 

   integer,  parameter :: dp = REAL128 !SELECTED_REAL_KIND(16,307)
   real(dp), parameter :: pi = 4.0e0_dp*datan(1.0e0_dp)
   real(dp), parameter :: h  = 6.62607004e-34_dp
   real(dp), parameter :: hbar  = h/(2.0e0_dp*pi)
   real(dp), parameter :: cl  = 299792458.e0_dp
   real(dp), parameter :: eps0  = 8.85418781762e-12_dp
   real(dp), parameter :: elec  = 1.60217662e-19_dp
   real(dp), parameter :: m0  = 9.10938356e-31_dp
   real(dp), parameter :: Cm_to_D  = 3.33564e-30_dp

end module Constants 
