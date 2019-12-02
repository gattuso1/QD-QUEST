module Constants_au

use omp_lib
use, intrinsic :: iso_fortran_env

implicit none 

   integer,  parameter :: dp = SELECTED_REAL_KIND(15,307)
   real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)
   real(dp), parameter :: E_au = 5.1422065211e11_dp
   real(dp), parameter :: Dip_au = 8.4783532619e-30_dp
   real(dp), parameter :: a0 = 5.291772109217e-11_dp
   real(dp), parameter :: t_au = 2.41888432650516e-17_dp
   real(dp), parameter :: Energ_au = 4.3597441775e-18_dp
   real(dp), parameter :: h  = 6.62607004e-34_dp
   real(dp), parameter :: eV_cm  = 8065.54429e0_dp
   real(dp), parameter :: hbar_au  = 1.e0_dp
   real(dp), parameter :: hbar = h/(2.0e0_dp*pi)
   real(dp), parameter :: cl  = 299792458.e0_dp
   real(dp), parameter :: eps0  = 8.85418781762e-12_dp
   real(dp), parameter :: elec  = 1.60217662e-19_dp
   real(dp), parameter :: elec_au  = 1.e0_dp
   real(dp), parameter :: m0  = 9.10938356e-31_dp
   real(dp), parameter :: Cm_to_D  = 3.33564e-30_dp
   real(dp), parameter :: D_to_au  = 2.5417462310548435e0_dp 

end module Constants_au 
