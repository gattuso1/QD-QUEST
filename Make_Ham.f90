module Make_Ham

use omp_lib
use Constants_au
use Variables_au
use Integrals
use Normal

real(dp), parameter :: a11_1d_ho = 0.0191697d0  
real(dp), parameter :: a11_2d_ho = 0.242775d0  
real(dp), parameter :: a11_3d_ho = 1.28322d0 
real(dp), parameter :: a11_1e_ho = 0.0131599d0       
real(dp), parameter :: a11_2e_ho = 0.199444d0        
real(dp), parameter :: a11_3e_ho = 1.10133d0
real(dp), parameter :: a22_1d_ho = 0.0178298d0  
real(dp), parameter :: a22_2d_ho = 0.239631d0   
real(dp), parameter :: a22_3d_ho = 1.25828d0
real(dp), parameter :: a22_1e_ho = 0.00477726d0      
real(dp), parameter :: a22_2e_ho = 0.0379395d0       
real(dp), parameter :: a22_3e_ho = 1.66235d0
real(dp), parameter :: a12_1d_ho = -0.00288509
real(dp), parameter :: a12_2d_ho = 0.0260066
real(dp), parameter :: a12_3d_ho = 0.57151
real(dp), parameter :: a12_1e_ho = -0.00651569   
real(dp), parameter :: a12_2e_ho = 0.0530662     
real(dp), parameter :: a12_3e_ho = 1.71364


contains

subroutine make_Ham_he


if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif


Ham     = 0.d0
Ham_0   = 0.d0
Ham_dir = 0.d0
Ham_ex  = 0.d0

Ham_0(1)     = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1d9)**a11_3e_ho))

Ham_0(2)     = minEe(1,n) + minEh(2,n) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1d9)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n)*1d9)**a22_3e_ho))

Ham_0(3)     = minEe(1,n+ndim) + minEh(1,n+ndim) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n+ndim)*1d9)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n+ndim)*1d9)**a11_3e_ho))

Ham_0(4)     = minEe(1,n+ndim) + minEh(2,n+ndim) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n+ndim)*1d9)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n+ndim)*1d9)**a22_3e_ho))

Ham_0(5)     = minEe(1,n+ndim) + minEh(1,n) + V0 
Ham_dir(5,5) = elec*(a55_1d_he + a55_2d_he / ((aR(n)*1d9)**a55_3d_he * (aR(n+ndim)*1d9)**a55_4d_he)) 
Ham_ex(5,5)  = elec*(a55_1e_he + a55_2e_he / ((aR(n)*1d9)**a55_3e_he * (aR(n+ndim)*1d9)**a55_4e_he))
                                  
Ham_0(6)     = minEe(1,n+ndim) + minEh(2,n) + V0 
Ham_dir(6,6) = elec*(a66_1d_he + a66_2d_he / ((aR(n)*1d9)**a66_3d_he * (aR(n+ndim)*1d9)**a66_4d_he)) 
Ham_ex(6,6)  = elec*(a66_1e_he + a66_2e_he / ((aR(n)*1d9)**a66_3e_he * (aR(n+ndim)*1d9)**a66_4e_he)) 
                                  
Ham_0(7)     = minEe(1,n) + minEh(1,n+ndim) + V0 
Ham_dir(7,7) = elec*(a55_1d_he + a55_2d_he / ((aR(n+ndim)*1d9)**a55_3d_he * (aR(n)*1d9)**a55_4d_he)) 
Ham_ex(7,7)  = elec*(a55_1e_he + a55_2e_he / ((aR(n+ndim)*1d9)**a55_3e_he * (aR(n)*1d9)**a55_4e_he))

Ham_0(8)     = minEe(1,n) + minEh(2,n+ndim) + V0 
Ham_dir(8,8) = elec*(a66_1d_he + a66_2d_he / ((aR(n+ndim)*1d9)**a66_3d_he * (aR(n)*1d9)**a66_4d_he)) 
Ham_ex(8,8)  = elec*(a66_1e_he + a66_2e_he / ((aR(n+ndim)*1d9)**a66_3e_he * (aR(n)*1d9)**a66_4e_he))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1d9)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1d9)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aR(n)*1d9)**a13_3d_he * (aR(n+ndim)*1d9)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aR(n)*1d9)**a13_3e_he * (aR(n+ndim)*1d9)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aR(n)*1d9)**a14_3d_he * (aR(n+ndim)*1d9)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aR(n)*1d9)**a14_3e_he * (aR(n+ndim)*1d9)**a14_4e_he))

Ham_dir(1,5) = elec*(a15_1d_he + a15_2d_he / ((aR(n)*1d9)**a15_3d_he * (aR(n+ndim)*1d9)**a15_4d_he))
Ham_ex(1,5)  = elec*(a15_1e_he + a15_2e_he / ((aR(n)*1d9)**a15_3e_he * (aR(n+ndim)*1d9)**a15_4e_he))

Ham_dir(1,6) = elec*(a16_1d_he + a16_2d_he / ((aR(n)*1d9)**a16_3d_he * (aR(n+ndim)*1d9)**a16_4d_he))
Ham_ex(1,6)  = elec*(a16_1e_he + a16_2e_he / ((aR(n)*1d9)**a16_3e_he * (aR(n+ndim)*1d9)**a16_4e_he))

Ham_dir(1,7) = elec*(a17_1d_he + a17_2d_he / ((aR(n)*1d9)**a17_3d_he * (aR(n+ndim)*1d9)**a17_4d_he))
Ham_ex(1,7)  = elec*(a17_1e_he + a17_2e_he / ((aR(n)*1d9)**a17_3e_he * (aR(n+ndim)*1d9)**a17_4e_he))

Ham_dir(1,8) = elec*(a18_1d_he + a18_2d_he / ((aR(n)*1d9)**a18_3d_he * (aR(n+ndim)*1d9)**a18_4d_he))
Ham_ex(1,8)  = elec*(a18_1e_he + a18_2e_he / ((aR(n)*1d9)**a18_3e_he * (aR(n+ndim)*1d9)**a18_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aR(n+ndim)*1d9)**a14_3d_he * (aR(n)*1d9)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aR(n+ndim)*1d9)**a14_3e_he * (aR(n)*1d9)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aR(n)*1d9)**a24_3d_he * (aR(n+ndim)*1d9)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aR(n)*1d9)**a24_3e_he * (aR(n+ndim)*1d9)**a24_4e_he))

Ham_dir(2,5) = elec*(a25_1d_he + a25_2d_he / ((aR(n)*1d9)**a25_3d_he * (aR(n+ndim)*1d9)**a25_4d_he))
Ham_ex(2,5)  = elec*(a25_1e_he + a25_2e_he / ((aR(n)*1d9)**a25_3e_he * (aR(n+ndim)*1d9)**a25_4e_he))

Ham_dir(2,6) = elec*(a26_1d_he + a26_2d_he / ((aR(n)*1d9)**a26_3d_he * (aR(n+ndim)*1d9)**a26_4d_he))
Ham_ex(2,6)  = elec*(a26_1e_he + a26_2e_he / ((aR(n)*1d9)**a26_3e_he * (aR(n+ndim)*1d9)**a26_4e_he))

Ham_dir(2,7) = elec*(a27_1d_he + a27_2d_he / ((aR(n)*1d9)**a27_3d_he * (aR(n+ndim)*1d9)**a27_4d_he))
Ham_ex(2,7)  = elec*(a27_1e_he + a27_2e_he / ((aR(n)*1d9)**a27_3e_he * (aR(n+ndim)*1d9)**a27_4e_he))

Ham_dir(2,8) = elec*(a28_1d_he + a28_2d_he / ((aR(n)*1d9)**a28_3d_he * (aR(n+ndim)*1d9)**a28_4d_he))
Ham_ex(2,8)  = elec*(a28_1e_he + a28_2e_he / ((aR(n)*1d9)**a28_3e_he * (aR(n+ndim)*1d9)**a28_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n+ndim)*1d9)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n+ndim)*1d9)**a12_3e_ho))

Ham_dir(3,5) = elec*(a17_1d_he + a17_2d_he / ((aR(n+ndim)*1d9)**a17_3d_he * (aR(n)*1d9)**a17_4d_he)) 
Ham_ex(3,5)  = elec*(a17_1e_he + a17_2e_he / ((aR(n+ndim)*1d9)**a17_3e_he * (aR(n)*1d9)**a17_4e_he))

Ham_dir(3,6) = elec*(a18_1d_he + a18_2d_he / ((aR(n+ndim)*1d9)**a18_3d_he * (aR(n)*1d9)**a18_4d_he)) 
Ham_ex(3,6)  = elec*(a18_1e_he + a18_2e_he / ((aR(n+ndim)*1d9)**a18_3e_he * (aR(n)*1d9)**a18_4e_he))

Ham_dir(3,7) = elec*(a15_1d_he + a15_2d_he / ((aR(n+ndim)*1d9)**a15_3d_he * (aR(n)*1d9)**a15_4d_he)) 
Ham_ex(3,7)  = elec*(a15_1e_he + a15_2e_he / ((aR(n+ndim)*1d9)**a15_3e_he * (aR(n)*1d9)**a15_4e_he))

Ham_dir(3,8) = elec*(a16_1d_he + a16_2d_he / ((aR(n+ndim)*1d9)**a16_3d_he * (aR(n)*1d9)**a16_4d_he)) 
Ham_ex(3,8)  = elec*(a16_1e_he + a16_2e_he / ((aR(n+ndim)*1d9)**a16_3e_he * (aR(n)*1d9)**a16_4e_he))

Ham_dir(4,5) = elec*(a27_1d_he + a27_2d_he / ((aR(n+ndim)*1d9)**a27_3d_he * (aR(n)*1d9)**a27_4d_he)) 
Ham_ex(4,5)  = elec*(a27_1e_he + a27_2e_he / ((aR(n+ndim)*1d9)**a27_3e_he * (aR(n)*1d9)**a27_4e_he))

Ham_dir(4,6) = elec*(a28_1d_he + a28_2d_he / ((aR(n+ndim)*1d9)**a28_3d_he * (aR(n)*1d9)**a28_4d_he)) 
Ham_ex(4,6)  = elec*(a28_1e_he + a28_2e_he / ((aR(n+ndim)*1d9)**a28_3e_he * (aR(n)*1d9)**a28_4e_he))

Ham_dir(4,7) = elec*(a25_1d_he + a25_2d_he / ((aR(n+ndim)*1d9)**a25_3d_he * (aR(n)*1d9)**a25_4d_he)) 
Ham_ex(4,7)  = elec*(a25_1e_he + a25_2e_he / ((aR(n+ndim)*1d9)**a25_3e_he * (aR(n)*1d9)**a25_4e_he))

Ham_dir(4,8) = elec*(a26_1d_he + a26_2d_he / ((aR(n+ndim)*1d9)**a26_3d_he * (aR(n)*1d9)**a26_4d_he)) 
Ham_ex(4,8)  = elec*(a26_1e_he + a26_2e_he / ((aR(n+ndim)*1d9)**a26_3e_he * (aR(n)*1d9)**a26_4e_he))

Ham_dir(5,6) = elec*(a56_1d_he + a56_2d_he / ((aR(n)*1d9)**a56_3d_he * (aR(n+ndim)*1d9)**a56_4d_he)) 
Ham_ex(5,6)  = elec*(a56_1e_he + a56_2e_he / ((aR(n)*1d9)**a56_3e_he * (aR(n+ndim)*1d9)**a56_4e_he))

Ham_dir(5,7) = elec*(a57_1d_he + a57_2d_he / ((aR(n)*1d9)**a57_3d_he * (aR(n+ndim)*1d9)**a57_4d_he))  
Ham_ex(5,7)  = elec*(a57_1e_he + a57_2e_he / ((aR(n)*1d9)**a57_3e_he * (aR(n+ndim)*1d9)**a57_4e_he))   

Ham_dir(5,8) = elec*(a58_1d_he + a58_2d_he / ((aR(n)*1d9)**a58_3d_he * (aR(n+ndim)*1d9)**a58_4d_he))  
Ham_ex(5,8)  = elec*(a58_1e_he + a58_2e_he / ((aR(n)*1d9)**a58_3e_he * (aR(n+ndim)*1d9)**a58_4e_he))   

Ham_dir(6,7) = elec*(a58_1d_he + a58_2d_he / ((aR(n+ndim)*1d9)**a58_3d_he * (aR(n)*1d9)**a58_4d_he))  
Ham_ex(6,7)  = elec*(a58_1e_he + a58_2e_he / ((aR(n+ndim)*1d9)**a58_3e_he * (aR(n)*1d9)**a58_4e_he))   

Ham_dir(6,8) = elec*(a68_1d_he + a68_2d_he / ((aR(n)*1d9)**a68_3d_he * (aR(n+ndim)*1d9)**a68_4d_he)) 
Ham_ex(6,8)  = elec*(a68_1e_he + a68_2e_he / ((aR(n)*1d9)**a68_3e_he * (aR(n+ndim)*1d9)**a68_4e_he))   

Ham_dir(7,8) = elec*(a56_1d_he + a56_2d_he / ((aR(n+ndim)*1d9)**a56_3d_he * (aR(n)*1d9)**a56_4d_he)) 
Ham_ex(7,8)  = elec*(a78_1e_he + a78_2e_he / ((aR(n+ndim)*1d9)**a78_3e_he * (aR(n)*1d9)**a78_4e_he))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.d0 * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au

if ( inbox .eq. "y" ) then

TransHam_d = 0.d0
TransHam_l = 0.d0

TransHam_l(0,1,:) = vector(TransDip_Ana_h1e(n))
TransHam_l(0,2,:) = vector(TransDip_Ana_h2e(n))
TransHam_l(0,3,:) = vector(TransDip_Ana_h1e(n+ndim))
TransHam_l(0,4,:) = vector(TransDip_Ana_h2e(n+ndim))
TransHam_l(0,5,:) = vector(TransDip_Fit_h1e_he(aR(n+ndim),aR(n)))
TransHam_l(0,6,:) = vector(TransDip_Fit_h2e_he(aR(n+ndim),aR(n)))
TransHam_l(0,7,:) = vector(TransDip_Fit_h1e_he(aR(n),aR(n+ndim)))
TransHam_l(0,8,:) = vector(TransDip_Fit_h2e_he(aR(n),aR(n+ndim)))
TransHam_l(1,2,:) = vector(TransDip_Ana_h1h2(n))
TransHam_l(1,5,:) = vector(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(1,7,:) = vector(TransDip_Fit_h1h1_he(aR(n+ndim),aR(n)))
TransHam_l(1,8,:) = vector(TransDip_Fit_h1h2_he(aR(n),aR(n+ndim)))
TransHam_l(2,6,:) = vector(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(2,7,:) = vector(TransDip_Fit_h1h2_he(aR(n+ndim),aR(n)))
TransHam_l(2,8,:) = vector(TransDip_Fit_h2h2_he(aR(n+ndim),aR(n)))
TransHam_l(3,4,:) = vector(TransDip_Ana_h1h2(n+ndim))
TransHam_l(3,5,:) = vector(TransDip_Fit_h1h1_he(aR(n+ndim),aR(n)))
TransHam_l(3,6,:) = vector(TransDip_Fit_h1h2_he(aR(n+ndim),aR(n)))
TransHam_l(3,7,:) = vector(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(4,5,:) = vector(TransDip_Fit_h1h2_he(aR(n),aR(n+ndim)))
TransHam_l(4,6,:) = vector(TransDip_Fit_h2h2_he(aR(n+ndim),aR(n)))
TransHam_l(4,8,:) = vector(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(5,6,:) = vector(TransDip_Ana_h1h2(n))
TransHam_l(7,8,:) = vector(TransDip_Ana_h1h2(n+ndim))

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

else 

TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(0,3) = TransDip_Ana_h1e(n+ndim)
TransHam(0,4) = TransDip_Ana_h2e(n+ndim)
TransHam(0,5) = TransDip_Fit_h1e_he(aR(n+ndim),aR(n))
TransHam(0,6) = TransDip_Fit_h2e_he(aR(n+ndim),aR(n))
TransHam(0,7) = TransDip_Fit_h1e_he(aR(n),aR(n+ndim))
TransHam(0,8) = TransDip_Fit_h2e_he(aR(n),aR(n+ndim))
TransHam(1,2) = TransDip_Ana_h1h2(n)
TransHam(1,5) = TransDip_Fit_ee_he(aR(n+ndim),aR(n))
TransHam(1,7) = TransDip_Fit_h1h1_he(aR(n+ndim),aR(n))
TransHam(1,8) = TransDip_Fit_h1h2_he(aR(n),aR(n+ndim))
TransHam(2,6) = TransDip_Fit_ee_he(aR(n+ndim),aR(n))
TransHam(2,7) = TransDip_Fit_h1h2_he(aR(n+ndim),aR(n))
TransHam(2,8) = TransDip_Fit_h2h2_he(aR(n+ndim),aR(n))
TransHam(3,4) = TransDip_Ana_h1h2(n+ndim)
TransHam(3,5) = TransDip_Fit_h1h1_he(aR(n+ndim),aR(n))
TransHam(3,6) = TransDip_Fit_h1h2_he(aR(n+ndim),aR(n))
TransHam(3,7) = TransDip_Fit_ee_he(aR(n+ndim),aR(n))
TransHam(4,5) = TransDip_Fit_h1h2_he(aR(n),aR(n+ndim))
TransHam(4,6) = TransDip_Fit_h2h2_he(aR(n+ndim),aR(n))
TransHam(4,8) = TransDip_Fit_ee_he(aR(n+ndim),aR(n))
TransHam(5,6) = TransDip_Ana_h1h2(n)
TransHam(7,8) = TransDip_Ana_h1h2(n+ndim)

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

endif

TransHam = TransHam/D_to_au

end subroutine make_Ham_he

subroutine make_Ham_fineSt

Ham = 0.d0

Ham(1,1)   = Eeh1(n)  - Dso1/3.d0 - 3.d0*Kas
Ham(2,2)   = Eeh1(n)  - Dso1/3.d0 - 3.d0*Kcs
Ham(3,3)   = Eeh1(n)  - Kas
Ham(4,4)   = Eeh1(n)  - 3.d0*Kas
Ham(5,5)   = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(6,6)   = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(7,7)   = Eeh1(n)  - Kcs
Ham(8,8)   = Eeh1(n)  - 3.d0*Kcs
Ham(9,9)   = Eeh1(n)  + Dso1/3.d0 - 3.d0*Kas
Ham(10,10) = Eeh1(n)  - Kbs + Dxf
Ham(11,11) = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(12,12) = Eeh1(n)  + Dso1/3.d0 - 3.d0*Kcs

Ham(3,4) = -1.d0*Dso1/3.d0
Ham(4,3) = Ham(3,4)
Ham(3,5) = Ham(3,4)
Ham(5,3) = Ham(3,4)
Ham(6,7) = Ham(3,4)
Ham(7,6) = Ham(3,4)
Ham(6,8) = Ham(3,4)
Ham(8,6) = Ham(3,4)
Ham(9,10) = Ham(3,4)
Ham(10,9) = Ham(3,4)
Ham(9,11) = Ham(3,4)
Ham(11,9) = Ham(3,4)
Ham(10,12) = Ham(3,4)
Ham(12,10) = Ham(3,4)

Ham(4,5) = Dso1/3.d0
Ham(5,4) = Ham(4,5)
Ham(7,8) = Ham(4,5)
Ham(8,7) = Ham(4,5)
Ham(11,12) = Ham(4,5)
Ham(12,11) = Ham(4,5)

Ham(13,13)   = Eeh2(n)  - Dso2/3.d0 - 3.d0*Kas
Ham(14,14)   = Eeh2(n)  - Dso2/3.d0 - 3.d0*Kcs
Ham(15,15)   = Eeh2(n)  - Kas
Ham(16,16)   = Eeh2(n)  - 3.d0*Kas
Ham(17,17)   = Eeh2(n)  - 3*Kbs + Dxf
Ham(18,18)   = Eeh2(n)  - 3*Kbs + Dxf
Ham(19,19)   = Eeh2(n)  - Kcs
Ham(20,20)   = Eeh2(n)  - 3.d0*Kcs
Ham(21,21)   = Eeh2(n)  + Dso2/3.d0 - 3.d0*Kas
Ham(22,22)   = Eeh2(n)  - Kbs + Dxf
Ham(23,23)   = Eeh2(n)  - 3.d0*Kbs + Dxf
Ham(24,24)   = Eeh2(n)  + Dso2/3.d0 - 3.d0*Kcs

Ham(15,16) = -1.d0*Dso2/3.d0
Ham(16,15) = Ham(15,16)
Ham(15,17) = Ham(15,16)
Ham(17,15) = Ham(15,16)
Ham(18,19) = Ham(15,16)
Ham(19,18) = Ham(15,16)
Ham(18,20) = Ham(15,16)
Ham(20,18) = Ham(15,16)
Ham(21,22) = Ham(15,16)
Ham(22,21) = Ham(15,16)
Ham(21,23) = Ham(15,16)
Ham(23,21) = Ham(15,16)
Ham(22,24) = Ham(15,16)
Ham(24,22) = Ham(15,16)

Ham(16,17) = Dso2/3
Ham(17,16) = Ham(16,17)
Ham(19,20) = Ham(16,17)
Ham(20,19) = Ham(16,17)
Ham(23,24) = Ham(16,17)
Ham(24,23) = Ham(16,17)

Ham = Ham/Energ_au

end subroutine make_Ham_fineSt

subroutine make_TransHam_0_fineSt

TransHam0 = 0.d0

TransHam0(0,3)  = abs(TransDip_Ana_h1e(n))
TransHam0(0,7)  = abs(TransDip_Ana_h1e(n)) 
TransHam0(0,10) = abs(TransDip_Ana_h1e(n))
TransHam0(0,15) = abs(TransDip_Ana_h2e(n))
TransHam0(0,19) = abs(TransDip_Ana_h2e(n))
TransHam0(0,22) = abs(TransDip_Ana_h2e(n))

TransHam0 = TransHam0/Dip_au

do i=0,nstates-1
TransHam0(i,0) = TransHam0(0,i)
enddo

end subroutine make_TransHam_0_fineSt

subroutine make_TransHam_ei_fineSt

TransHam = 0.d0

do i=0,nstates-1
do j=0,nstates-1
TransHam(0,i) = TransHam(0,i) +  TransHam0(0,j) * Ham_ei(j,i)
enddo
enddo

TransHam = abs(TransHam)

write(6,*) (abs(TransHam(0,i)*Dip_au/Cm_to_D), i=0,nstates-1)

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

Ham = 0.d0

do i=0,nstates-1
Ham(i,i) = lambda(i)
enddo

end subroutine make_TransHam_ei_fineSt

subroutine make_Ham_l

Ham_l = 0.e0_dp

do i=0,nstates-1
Ham_l(i,i) = lambda(i)
enddo

end subroutine make_Ham_l

subroutine make_Ham_singl

Ham     = 0.0
Ham_0   = 0.0
Ham_dir = 0.0
Ham_ex  = 0.0

Ham_0(1)   = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1d9)**a11_3e_ho))

Ham_0(2)   = minEe(1,n) + minEh(2,n) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1d9)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n)*1d9)**a22_3e_ho))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1d9)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1d9)**a12_3e_ho))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.d0 * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au

if ( inbox .eq. "y" ) then

TransHam_d = 0.d0
TransHam_l = 0.d0

if ( TDM_ee .eq. 'n') then
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo
elseif ( TDM_ee .eq. 'y') then
TransHam_l(0,1,:) = vector(TransDip_Ana_h1e(n))
TransHam_l(0,2,:) = vector(TransDip_Ana_h2e(n))
TransHam_l(1,2,:) = vector(TransDip_Ana_h1h2(n))

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo
endif

else 

if ( TDM_ee .eq. 'n') then
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo
elseif ( TDM_ee .eq. 'y') then
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(1,2) = TransDip_Ana_h1h2(n)

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo
endif

endif

TransHam = TransHam/D_to_au

end subroutine make_Ham_singl

end module Make_Ham
