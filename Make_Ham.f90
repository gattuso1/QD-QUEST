module Make_Ham

use omp_lib
use Constants_au
use Variables_au
use Integrals
use Normal

real(dp), parameter :: a11_1d_ho = 0.0191697e0_dp  
real(dp), parameter :: a11_2d_ho = 0.242775e0_dp  
real(dp), parameter :: a11_3d_ho = 1.28322e0_dp 
real(dp), parameter :: a11_1e_ho = 0.0131599e0_dp       
real(dp), parameter :: a11_2e_ho = 0.199444e0_dp        
real(dp), parameter :: a11_3e_ho = 1.10133e0_dp
real(dp), parameter :: a22_1d_ho = 0.0178298e0_dp  
real(dp), parameter :: a22_2d_ho = 0.239631e0_dp   
real(dp), parameter :: a22_3d_ho = 1.25828e0_dp
real(dp), parameter :: a22_1e_ho = 0.00477726e0_dp      
real(dp), parameter :: a22_2e_ho = 0.0379395e0_dp       
real(dp), parameter :: a22_3e_ho = 1.66235e0_dp
real(dp), parameter :: a12_1d_ho = -0.00288509e0_dp
real(dp), parameter :: a12_2d_ho = 0.0260066e0_dp
real(dp), parameter :: a12_3d_ho = 0.57151e0_dp
real(dp), parameter :: a12_1e_ho = -0.00651569e0_dp  
real(dp), parameter :: a12_2e_ho = 0.0530662e0_dp     
real(dp), parameter :: a12_3e_ho = 1.71364e0_dp


contains

subroutine make_Ham_he

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

Ham_0(1)     = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1e9_dp)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1e9_dp)**a11_3e_ho))

Ham_0(2)     = minEe(1,n) + minEh(2,n) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n)*1e9_dp)**a22_3e_ho))

Ham_0(3)     = minEe(1,n+ndim) + minEh(1,n+ndim) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n+ndim)*1e9_dp)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n+ndim)*1e9_dp)**a11_3e_ho))

Ham_0(4)     = minEe(1,n+ndim) + minEh(2,n+ndim) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n+ndim)*1e9_dp)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n+ndim)*1e9_dp)**a22_3e_ho))

Ham_0(5)     = minEe(1,n+ndim) + minEh(1,n) + V0 
Ham_dir(5,5) = elec*(a55_1d_he + a55_2d_he / ((aR(n)*1e9_dp)**a55_3d_he * (aR(n+ndim)*1d9)**a55_4d_he)) 
Ham_ex(5,5)  = elec*(a55_1e_he + a55_2e_he / ((aR(n)*1e9_dp)**a55_3e_he * (aR(n+ndim)*1d9)**a55_4e_he))
                                  
Ham_0(6)     = minEe(1,n+ndim) + minEh(2,n) + V0 
Ham_dir(6,6) = elec*(a66_1d_he + a66_2d_he / ((aR(n)*1e9_dp)**a66_3d_he * (aR(n+ndim)*1d9)**a66_4d_he)) 
Ham_ex(6,6)  = elec*(a66_1e_he + a66_2e_he / ((aR(n)*1e9_dp)**a66_3e_he * (aR(n+ndim)*1d9)**a66_4e_he)) 
                                  
Ham_0(7)     = minEe(1,n) + minEh(1,n+ndim) + V0 
Ham_dir(7,7) = elec*(a55_1d_he + a55_2d_he / ((aR(n+ndim)*1e9_dp)**a55_3d_he * (aR(n)*1d9)**a55_4d_he)) 
Ham_ex(7,7)  = elec*(a55_1e_he + a55_2e_he / ((aR(n+ndim)*1e9_dp)**a55_3e_he * (aR(n)*1d9)**a55_4e_he))

Ham_0(8)     = minEe(1,n) + minEh(2,n+ndim) + V0 
Ham_dir(8,8) = elec*(a66_1d_he + a66_2d_he / ((aR(n+ndim)*1e9_dp)**a66_3d_he * (aR(n)*1d9)**a66_4d_he)) 
Ham_ex(8,8)  = elec*(a66_1e_he + a66_2e_he / ((aR(n+ndim)*1e9_dp)**a66_3e_he * (aR(n)*1d9)**a66_4e_he))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1e9_dp)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1e9_dp)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aR(n)*1e9_dp)**a13_3d_he * (aR(n+ndim)*1d9)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aR(n)*1e9_dp)**a13_3e_he * (aR(n+ndim)*1d9)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aR(n)*1e9_dp)**a14_3d_he * (aR(n+ndim)*1d9)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aR(n)*1e9_dp)**a14_3e_he * (aR(n+ndim)*1d9)**a14_4e_he))

Ham_dir(1,5) = elec*(a15_1d_he + a15_2d_he / ((aR(n)*1e9_dp)**a15_3d_he * (aR(n+ndim)*1d9)**a15_4d_he))
Ham_ex(1,5)  = elec*(a15_1e_he + a15_2e_he / ((aR(n)*1e9_dp)**a15_3e_he * (aR(n+ndim)*1d9)**a15_4e_he))

Ham_dir(1,6) = elec*(a16_1d_he + a16_2d_he / ((aR(n)*1e9_dp)**a16_3d_he * (aR(n+ndim)*1d9)**a16_4d_he))
Ham_ex(1,6)  = elec*(a16_1e_he + a16_2e_he / ((aR(n)*1e9_dp)**a16_3e_he * (aR(n+ndim)*1d9)**a16_4e_he))

Ham_dir(1,7) = elec*(a17_1d_he + a17_2d_he / ((aR(n)*1e9_dp)**a17_3d_he * (aR(n+ndim)*1d9)**a17_4d_he))
Ham_ex(1,7)  = elec*(a17_1e_he + a17_2e_he / ((aR(n)*1e9_dp)**a17_3e_he * (aR(n+ndim)*1d9)**a17_4e_he))

Ham_dir(1,8) = elec*(a18_1d_he + a18_2d_he / ((aR(n)*1e9_dp)**a18_3d_he * (aR(n+ndim)*1d9)**a18_4d_he))
Ham_ex(1,8)  = elec*(a18_1e_he + a18_2e_he / ((aR(n)*1e9_dp)**a18_3e_he * (aR(n+ndim)*1d9)**a18_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aR(n+ndim)*1e9_dp)**a14_3d_he * (aR(n)*1d9)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aR(n+ndim)*1e9_dp)**a14_3e_he * (aR(n)*1d9)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aR(n)*1e9_dp)**a24_3d_he * (aR(n+ndim)*1d9)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aR(n)*1e9_dp)**a24_3e_he * (aR(n+ndim)*1d9)**a24_4e_he))

Ham_dir(2,5) = elec*(a25_1d_he + a25_2d_he / ((aR(n)*1e9_dp)**a25_3d_he * (aR(n+ndim)*1d9)**a25_4d_he))
Ham_ex(2,5)  = elec*(a25_1e_he + a25_2e_he / ((aR(n)*1e9_dp)**a25_3e_he * (aR(n+ndim)*1d9)**a25_4e_he))

Ham_dir(2,6) = elec*(a26_1d_he + a26_2d_he / ((aR(n)*1e9_dp)**a26_3d_he * (aR(n+ndim)*1d9)**a26_4d_he))
Ham_ex(2,6)  = elec*(a26_1e_he + a26_2e_he / ((aR(n)*1e9_dp)**a26_3e_he * (aR(n+ndim)*1d9)**a26_4e_he))

Ham_dir(2,7) = elec*(a27_1d_he + a27_2d_he / ((aR(n)*1e9_dp)**a27_3d_he * (aR(n+ndim)*1d9)**a27_4d_he))
Ham_ex(2,7)  = elec*(a27_1e_he + a27_2e_he / ((aR(n)*1e9_dp)**a27_3e_he * (aR(n+ndim)*1d9)**a27_4e_he))

Ham_dir(2,8) = elec*(a28_1d_he + a28_2d_he / ((aR(n)*1e9_dp)**a28_3d_he * (aR(n+ndim)*1d9)**a28_4d_he))
Ham_ex(2,8)  = elec*(a28_1e_he + a28_2e_he / ((aR(n)*1e9_dp)**a28_3e_he * (aR(n+ndim)*1d9)**a28_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n+ndim)*1e9_dp)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n+ndim)*1e9_dp)**a12_3e_ho))

Ham_dir(3,5) = elec*(a17_1d_he + a17_2d_he / ((aR(n+ndim)*1e9_dp)**a17_3d_he * (aR(n)*1d9)**a17_4d_he)) 
Ham_ex(3,5)  = elec*(a17_1e_he + a17_2e_he / ((aR(n+ndim)*1e9_dp)**a17_3e_he * (aR(n)*1d9)**a17_4e_he))

Ham_dir(3,6) = elec*(a18_1d_he + a18_2d_he / ((aR(n+ndim)*1e9_dp)**a18_3d_he * (aR(n)*1d9)**a18_4d_he)) 
Ham_ex(3,6)  = elec*(a18_1e_he + a18_2e_he / ((aR(n+ndim)*1e9_dp)**a18_3e_he * (aR(n)*1d9)**a18_4e_he))

Ham_dir(3,7) = elec*(a15_1d_he + a15_2d_he / ((aR(n+ndim)*1e9_dp)**a15_3d_he * (aR(n)*1d9)**a15_4d_he)) 
Ham_ex(3,7)  = elec*(a15_1e_he + a15_2e_he / ((aR(n+ndim)*1e9_dp)**a15_3e_he * (aR(n)*1d9)**a15_4e_he))

Ham_dir(3,8) = elec*(a16_1d_he + a16_2d_he / ((aR(n+ndim)*1e9_dp)**a16_3d_he * (aR(n)*1d9)**a16_4d_he)) 
Ham_ex(3,8)  = elec*(a16_1e_he + a16_2e_he / ((aR(n+ndim)*1e9_dp)**a16_3e_he * (aR(n)*1d9)**a16_4e_he))

Ham_dir(4,5) = elec*(a27_1d_he + a27_2d_he / ((aR(n+ndim)*1e9_dp)**a27_3d_he * (aR(n)*1d9)**a27_4d_he)) 
Ham_ex(4,5)  = elec*(a27_1e_he + a27_2e_he / ((aR(n+ndim)*1e9_dp)**a27_3e_he * (aR(n)*1d9)**a27_4e_he))

Ham_dir(4,6) = elec*(a28_1d_he + a28_2d_he / ((aR(n+ndim)*1e9_dp)**a28_3d_he * (aR(n)*1d9)**a28_4d_he)) 
Ham_ex(4,6)  = elec*(a28_1e_he + a28_2e_he / ((aR(n+ndim)*1e9_dp)**a28_3e_he * (aR(n)*1d9)**a28_4e_he))

Ham_dir(4,7) = elec*(a25_1d_he + a25_2d_he / ((aR(n+ndim)*1e9_dp)**a25_3d_he * (aR(n)*1d9)**a25_4d_he)) 
Ham_ex(4,7)  = elec*(a25_1e_he + a25_2e_he / ((aR(n+ndim)*1e9_dp)**a25_3e_he * (aR(n)*1d9)**a25_4e_he))

Ham_dir(4,8) = elec*(a26_1d_he + a26_2d_he / ((aR(n+ndim)*1e9_dp)**a26_3d_he * (aR(n)*1d9)**a26_4d_he)) 
Ham_ex(4,8)  = elec*(a26_1e_he + a26_2e_he / ((aR(n+ndim)*1e9_dp)**a26_3e_he * (aR(n)*1d9)**a26_4e_he))

Ham_dir(5,6) = elec*(a56_1d_he + a56_2d_he / ((aR(n)*1e9_dp)**a56_3d_he * (aR(n+ndim)*1d9)**a56_4d_he)) 
Ham_ex(5,6)  = elec*(a56_1e_he + a56_2e_he / ((aR(n)*1e9_dp)**a56_3e_he * (aR(n+ndim)*1d9)**a56_4e_he))

Ham_dir(5,7) = elec*(a57_1d_he + a57_2d_he / ((aR(n)*1e9_dp)**a57_3d_he * (aR(n+ndim)*1d9)**a57_4d_he))  
Ham_ex(5,7)  = elec*(a57_1e_he + a57_2e_he / ((aR(n)*1e9_dp)**a57_3e_he * (aR(n+ndim)*1d9)**a57_4e_he))   

Ham_dir(5,8) = elec*(a58_1d_he + a58_2d_he / ((aR(n)*1e9_dp)**a58_3d_he * (aR(n+ndim)*1d9)**a58_4d_he))  
Ham_ex(5,8)  = elec*(a58_1e_he + a58_2e_he / ((aR(n)*1e9_dp)**a58_3e_he * (aR(n+ndim)*1d9)**a58_4e_he))   

Ham_dir(6,7) = elec*(a58_1d_he + a58_2d_he / ((aR(n+ndim)*1e9_dp)**a58_3d_he * (aR(n)*1d9)**a58_4d_he))  
Ham_ex(6,7)  = elec*(a58_1e_he + a58_2e_he / ((aR(n+ndim)*1e9_dp)**a58_3e_he * (aR(n)*1d9)**a58_4e_he))   

Ham_dir(6,8) = elec*(a68_1d_he + a68_2d_he / ((aR(n)*1e9_dp)**a68_3d_he * (aR(n+ndim)*1d9)**a68_4d_he)) 
Ham_ex(6,8)  = elec*(a68_1e_he + a68_2e_he / ((aR(n)*1e9_dp)**a68_3e_he * (aR(n+ndim)*1d9)**a68_4e_he))   

Ham_dir(7,8) = elec*(a56_1d_he + a56_2d_he / ((aR(n+ndim)*1e9_dp)**a56_3d_he * (aR(n)*1d9)**a56_4d_he)) 
Ham_ex(7,8)  = elec*(a78_1e_he + a78_2e_he / ((aR(n+ndim)*1e9_dp)**a78_3e_he * (aR(n)*1d9)**a78_4e_he))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.e0_dp * Ham_dir(i,j) + Ham_ex(i,j)
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

if ( rdm_ori .eq. "y" ) then

call random_number(psirot)
call random_number(phirot)
call random_number(thetarot)

psirot   = psirot   * 2._dp*pi
phirot   = phirot   * 2._dp*pi
thetarot = thetarot * 2._dp*pi

!rotmat(1,1)= cos(psirot)   * cos(phirot)   - cos(thetarot) * sin(phirot) * sin(psirot)
!rotmat(1,2)= cos(psirot)   * sin(phirot)   + cos(thetarot) * cos(phirot) * sin(psirot)
!rotmat(1,3)= sin(psirot)   * sin(thetarot)
!rotmat(2,1)=-sin(psirot)   * cos(phirot)   - cos(thetarot) * sin(phirot) * cos(psirot)
!rotmat(2,2)=-sin(psirot)   * sin(phirot)   + cos(thetarot) * cos(phirot) * cos(psirot) 
!rotmat(2,3)= cos(psirot)   * sin(thetarot)
!rotmat(3,1)= sin(thetarot) * sin(phirot)
!rotmat(3,2)=-sin(thetarot) * cos(phirot) 
!rotmat(3,3)= cos(thetarot)

rotmat(1,1)= cos(phirot)*cos(thetarot)*cos(psirot) - sin(phirot)*sin(psirot)
rotmat(1,2)=-cos(phirot)*cos(thetarot)*sin(psirot) - sin(phirot)*cos(psirot)
rotmat(1,3)= cos(phirot)*sin(thetarot)
rotmat(2,1)= sin(phirot)*cos(thetarot)*cos(psirot) + cos(phirot)*sin(psirot)
rotmat(2,2)=-sin(phirot)*cos(thetarot)*sin(psirot) + cos(phirot)*cos(psirot)
rotmat(2,3)= sin(phirot)*sin(thetarot)
rotmat(3,1)=-sin(thetarot)*cos(psirot)
rotmat(3,2)= sin(thetarot)*sin(psirot)
rotmat(3,3)= cos(thetarot)

TransHam_l(0,1,:) = matmul(rotmat,eTDM(0,1,:))*(TransDip_Ana_h1e(n))
TransHam_l(0,2,:) = matmul(rotmat,eTDM(0,2,:))*(TransDip_Ana_h2e(n))
TransHam_l(0,3,:) = matmul(rotmat,eTDM(0,3,:))*(TransDip_Ana_h1e(n+ndim))
TransHam_l(0,4,:) = matmul(rotmat,eTDM(0,4,:))*(TransDip_Ana_h2e(n+ndim))
TransHam_l(0,5,:) = matmul(rotmat,eTDM(0,5,:))*(TransDip_Fit_h1e_he(aR(n+ndim),aR(n)))
TransHam_l(0,6,:) = matmul(rotmat,eTDM(0,6,:))*(TransDip_Fit_h2e_he(aR(n+ndim),aR(n)))
TransHam_l(0,7,:) = matmul(rotmat,eTDM(0,7,:))*(TransDip_Fit_h1e_he(aR(n),aR(n+ndim)))
TransHam_l(0,8,:) = matmul(rotmat,eTDM(0,8,:))*(TransDip_Fit_h2e_he(aR(n),aR(n+ndim)))
TransHam_l(1,2,:) = matmul(rotmat,eTDM(1,2,:))*(TransDip_Ana_h1h2(n))
TransHam_l(1,5,:) = matmul(rotmat,eTDM(1,5,:))*(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(1,7,:) = matmul(rotmat,eTDM(1,7,:))*(TransDip_Fit_h1h1_he(aR(n+ndim),aR(n)))
TransHam_l(1,8,:) = matmul(rotmat,eTDM(1,8,:))*(TransDip_Fit_h1h2_he(aR(n),aR(n+ndim)))
TransHam_l(2,6,:) = matmul(rotmat,eTDM(2,6,:))*(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(2,7,:) = matmul(rotmat,eTDM(2,7,:))*(TransDip_Fit_h1h2_he(aR(n+ndim),aR(n)))
TransHam_l(2,8,:) = matmul(rotmat,eTDM(2,8,:))*(TransDip_Fit_h2h2_he(aR(n+ndim),aR(n)))
TransHam_l(3,4,:) = matmul(rotmat,eTDM(3,4,:))*(TransDip_Ana_h1h2(n+ndim))
TransHam_l(3,5,:) = matmul(rotmat,eTDM(3,5,:))*(TransDip_Fit_h1h1_he(aR(n+ndim),aR(n)))
TransHam_l(3,6,:) = matmul(rotmat,eTDM(3,6,:))*(TransDip_Fit_h1h2_he(aR(n+ndim),aR(n)))
TransHam_l(3,7,:) = matmul(rotmat,eTDM(3,7,:))*(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(4,5,:) = matmul(rotmat,eTDM(4,5,:))*(TransDip_Fit_h1h2_he(aR(n),aR(n+ndim)))
TransHam_l(4,6,:) = matmul(rotmat,eTDM(4,6,:))*(TransDip_Fit_h2h2_he(aR(n+ndim),aR(n)))
TransHam_l(4,8,:) = matmul(rotmat,eTDM(4,8,:))*(TransDip_Fit_ee_he(aR(n+ndim),aR(n)))
TransHam_l(5,6,:) = matmul(rotmat,eTDM(5,6,:))*(TransDip_Ana_h1h2(n))
TransHam_l(7,8,:) = matmul(rotmat,eTDM(7,8,:))*(TransDip_Ana_h1h2(n+ndim))

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

TransHam_l = TransHam_l/D_to_au

elseif ( rdm_ori .eq. "n" ) then

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

TransHam = TransHam/D_to_au

endif

end subroutine make_Ham_he

subroutine make_Ham_he_FO

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

Ham_0(1)     = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1e9_dp)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1e9_dp)**a11_3e_ho))

write(6,*) Ham_0(1)/elec, Ham_dir(1,1)/elec, Ham_ex(1,1)/elec

Ham_0(2)     = minEe(1,n) + minEh(2,n) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n)*1e9_dp)**a22_3e_ho))

Ham_0(3)     = minEe(1,n+ndim) + minEh(1,n+ndim) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n+ndim)*1e9_dp)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n+ndim)*1e9_dp)**a11_3e_ho))

Ham_0(4)     = minEe(1,n+ndim) + minEh(2,n+ndim) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n+ndim)*1e9_dp)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n+ndim)*1e9_dp)**a22_3e_ho))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1e9_dp)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1e9_dp)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aR(n)*1e9_dp)**a13_3d_he * (aR(n+ndim)*1e9_dp)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aR(n)*1e9_dp)**a13_3e_he * (aR(n+ndim)*1e9_dp)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aR(n)*1e9_dp)**a14_3d_he * (aR(n+ndim)*1e9_dp)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aR(n)*1e9_dp)**a14_3e_he * (aR(n+ndim)*1e9_dp)**a14_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aR(n+ndim)*1e9_dp)**a14_3d_he * (aR(n)*1e9_dp)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aR(n+ndim)*1e9_dp)**a14_3e_he * (aR(n)*1e9_dp)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aR(n)*1e9_dp)**a24_3d_he * (aR(n+ndim)*1e9_dp)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aR(n)*1e9_dp)**a24_3e_he * (aR(n+ndim)*1e9_dp)**a24_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n+ndim)*1e9_dp)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n+ndim)*1e9_dp)**a12_3e_ho))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.e0_dp * Ham_dir(i,j) + Ham_ex(i,j)
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

!if ( inbox .eq. "y" ) then

TransHam_l(0,1,:) = vector(TransDip_Ana_h1e(n))
TransHam_l(0,2,:) = vector(TransDip_Ana_h2e(n))
TransHam_l(0,3,:) = vector(TransDip_Ana_h1e(n+ndim))
TransHam_l(0,4,:) = vector(TransDip_Ana_h2e(n+ndim))
TransHam_l(1,2,:) = vector(TransDip_Ana_h1h2(n))
TransHam_l(3,4,:) = vector(TransDip_Ana_h1h2(n+ndim))

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

!else 
!
!TransHam(0,1) = TransDip_Ana_h1e(n)
!TransHam(0,2) = TransDip_Ana_h2e(n)
!TransHam(0,3) = TransDip_Ana_h1e(n+ndim)
!TransHam(0,4) = TransDip_Ana_h2e(n+ndim)
!TransHam(1,2) = TransDip_Ana_h1h2(n)
!TransHam(3,4) = TransDip_Ana_h1h2(n+ndim)
!
!do i=0,nstates-1
!do j=i+1,nstates-1
!TransHam(j,i) = TransHam(i,j)
!enddo
!enddo
!
!endif

TransHam = TransHam/D_to_au

end subroutine make_Ham_he_FO

subroutine make_Ham_fineSt

Ham(1,1)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(Dso1/3.d0 - 3.d0*Kas)
Ham(2,2)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(Dso1/3.d0 - 3.d0*Kcs)
Ham(3,3)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(Kas)
Ham(4,4)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(3.d0*Kas)
Ham(5,5)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(3.d0*Kbs + Dxf)
Ham(6,6)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(3.d0*Kbs + Dxf)
Ham(7,7)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(Kcs)
Ham(8,8)   = minEe(1,n) + minEh(1,n)  + V0  - elec*(3.d0*Kcs)
Ham(9,9)   = minEe(1,n) + minEh(1,n)  + V0  + elec*(Dso1/3.d0 - 3.d0*Kas)
Ham(10,10) = minEe(1,n) + minEh(1,n)  + V0  - elec*(Kbs + Dxf)
Ham(11,11) = minEe(1,n) + minEh(1,n)  + V0  - elec*(3.d0*Kbs + Dxf)
Ham(12,12) = minEe(1,n) + minEh(1,n)  + V0  + elec*(Dso1/3.d0 - 3.d0*Kcs)

Ham(3,4)  = elec*(-1.d0*Dso1/3.d0)
Ham(4,3)  = Ham(3,4)
Ham(3,5)  = Ham(3,4)
Ham(5,3)  = Ham(3,4)
Ham(6,7)  = Ham(3,4)
Ham(7,6)  = Ham(3,4)
Ham(6,8)  = Ham(3,4)
Ham(8,6)  = Ham(3,4)
Ham(9,10) = Ham(3,4)
Ham(10,9) = Ham(3,4)
Ham(9,11) = Ham(3,4)
Ham(11,9) = Ham(3,4)
Ham(10,12)= Ham(3,4)
Ham(12,10)= Ham(3,4)

Ham(4,5)   = elec*(Dso1/3.d0)
Ham(5,4)   = Ham(4,5)
Ham(7,8)   = Ham(4,5)
Ham(8,7)   = Ham(4,5)
Ham(11,12) = Ham(4,5)
Ham(12,11) = Ham(4,5)

Ham(13,13)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(Dso2/3.d0 - 3.d0*Kas)
Ham(14,14)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(Dso2/3.d0 - 3.d0*Kcs)
Ham(15,15)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(Kas)
Ham(16,16)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(3.d0*Kas)
Ham(17,17)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(3*Kbs + Dxf)
Ham(18,18)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(3*Kbs + Dxf)
Ham(19,19)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(Kcs)
Ham(20,20)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(3.d0*Kcs)
Ham(21,21)   = minEe(1,n) + minEh(2,n)  + V0  + elec*(Dso2/3.d0 - 3.d0*Kas)
Ham(22,22)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(Kbs + Dxf)
Ham(23,23)   = minEe(1,n) + minEh(2,n)  + V0  - elec*(3.d0*Kbs + Dxf)
Ham(24,24)   = minEe(1,n) + minEh(2,n)  + V0  + elec*(Dso2/3.d0 - 3.d0*Kcs)

Ham(15,16) = elec*(-1.d0*Dso2/3.d0)
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

Ham(16,17) = elec*(Dso2/3)
Ham(17,16) = Ham(16,17)
Ham(19,20) = Ham(16,17)
Ham(20,19) = Ham(16,17)
Ham(23,24) = Ham(16,17)
Ham(24,23) = Ham(16,17)

Ham = Ham/Energ_au

TransHam(0,3)  = abs(TransDip_Ana_h1e(n))
TransHam(0,7)  = abs(TransDip_Ana_h1e(n)) 
TransHam(0,10) = abs(TransDip_Ana_h1e(n))
TransHam(0,15) = abs(TransDip_Ana_h2e(n))
TransHam(0,19) = abs(TransDip_Ana_h2e(n))
TransHam(0,22) = abs(TransDip_Ana_h2e(n))

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_fineSt

subroutine make_Ham_l

do i=0,nstates-1
Ham_l(i,i) = lambda(i)
enddo

end subroutine make_Ham_l

subroutine make_Ham_singl

Ham_0(1)   = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1d9)**a11_3e_ho))

!write(6,*) abs((- 1.d0*Ham_dir(1,1) + Ham_ex(1,1))/elec)
!write(6,*) Ham_0(1)/elec, Ham_dir(1,1)/elec, Ham_ex(1,1)/elec

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

if ( rdm_ori .eq. "y" ) then

TransHam_l(0,1,:) = vector(TransDip_Ana_h1e(n))
TransHam_l(0,2,:) = vector(TransDip_Ana_h2e(n))
TransHam_l(1,2,:) = vector(TransDip_Ana_h1h2(n))

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

TransHam_l = TransHam_l/D_to_au

elseif  ( rdm_ori .eq. "n" ) then

TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(1,2) = TransDip_Ana_h1h2(n)

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

!write(6,*) TransDip_Ana_h1h2(n)

TransHam = TransHam/D_to_au

endif

end subroutine make_Ham_singl

subroutine make_Ham_FS_FO

!print*, (minEe(1,n) + minEh(1,n)), V0, Eeh1(n)*elec

Eg(1) = minEe(1,n) + minEh(1,n)  + V0 - Eeh1(n)*elec
!Eg(2) = minEe(1,n) + minEh(2,n)  + V0
Eg(2) = Eg(1) + elec*(-0.523865d0 + Eg(1)/elec * 0.293665d0)
Eg(3) = minEe(1,n+ndim) + minEh(1,n+ndim)  + V0 - Eeh1(n)*elec
!Eg(4) = minEe(1,n+ndim) + minEh(2,n+ndim)  + V0
Eg(4) = Eg(3) + elec*(-0.523865d0 + Eg(3)/elec * 0.293665d0)

TDM(1) = TransDip_Ana_h1e(n)
TDM(2) = TransDip_Ana_h2e(n)
TDM(3) = TransDip_Ana_h1e(n+ndim)
TDM(4) = TransDip_Ana_h2e(n+ndim)

!print*, Eg(1), Eg(2), Eg(3), Eg(4)

j = 0
k = 1
i = 0

if ( n .le. nQDA+nQDB ) then
nbands = 2
elseif ( n .gt. nQDA+nQDB ) then
nbands = 4 
endif

do while ( k .le. nbands ) 

!i = modulo(k-1,2)+1
!nsys = n
!
!if ( k .gt. 2 ) then
!nsys = n + ndim
!endif

!print*, i,j,k,n

Ham(j+1 ,j+1 ) = Eg(k)  - elec* (Dso1/3.d0 + 3.d0*Kas)
Ham(j+2 ,j+2 ) = Eg(k)  - elec* (Dso1/3.d0 + 3.d0*Kcs)
Ham(j+3 ,j+3 ) = Eg(k)  - elec* Kas
Ham(j+4 ,j+4 ) = Eg(k)  - elec* 3.d0*Kas
Ham(j+5 ,j+5 ) = Eg(k)  - elec* (3.d0*Kbs - Dxf)
Ham(j+6 ,j+6 ) = Eg(k)  - elec* (3.d0*Kbs - Dxf)
Ham(j+7 ,j+7 ) = Eg(k)  - elec* Kcs
Ham(j+8 ,j+8 ) = Eg(k)  - elec* 3.d0*Kcs
Ham(j+9 ,j+9 ) = Eg(k)  - elec* (Dso1/(-3.d0) + 3.d0*Kas)
Ham(j+10,j+10) = Eg(k)  - elec* (Kbs - Dxf)
Ham(j+11,j+11) = Eg(k)  - elec* (3.d0*Kbs - Dxf)
Ham(j+12,j+12) = Eg(k)  - elec* (Dso1/(-3.d0) + 3.d0*Kcs)

Ham(j+3 ,j+4 ) = elec*(-1.d0*Dso1/3.d0)
Ham(j+4 ,j+3 ) = Ham(j+3,j+4)
Ham(j+3 ,j+5 ) = Ham(j+3,j+4)
Ham(j+5 ,j+3 ) = Ham(j+3,j+4)
Ham(j+6 ,j+7 ) = Ham(j+3,j+4)
Ham(j+7 ,j+6 ) = Ham(j+3,j+4)
Ham(j+6 ,j+8 ) = Ham(j+3,j+4)
Ham(j+8 ,j+6 ) = Ham(j+3,j+4)
Ham(j+9 ,j+10) = Ham(j+3,j+4)
Ham(j+10,j+9 ) = Ham(j+3,j+4)
Ham(j+9 ,j+11) = Ham(j+3,j+4)
Ham(j+11,j+9 ) = Ham(j+3,j+4)
Ham(j+10,j+12) = Ham(j+3,j+4)
Ham(j+12,j+10) = Ham(j+3,j+4)

Ham(j+4 ,j+5 ) = elec*(Dso1/3.d0)
Ham(j+5 ,j+4 ) = Ham(j+4,j+5)
Ham(j+7 ,j+8 ) = Ham(j+4,j+5)
Ham(j+8 ,j+7 ) = Ham(j+4,j+5)
Ham(j+11,j+12) = Ham(j+4,j+5)
Ham(j+12,j+11) = Ham(j+4,j+5)

TransHam(0 ,j+3 ) = abs(TDM(k))/sqrt(3.d0) 
TransHam(0 ,j+7 ) = abs(TDM(k))/sqrt(3.d0)
TransHam(0 ,j+10) = abs(TDM(k))/sqrt(3.d0)

!Interdot fine structure dipole moment non diagonal elements between 0, -1, +1, S to S, T to T

!print*, TransHam(j+0 ,j+3 )

k = k + 1
j = j + 12

enddo

k = 1
i = 0

do while ( k .le. nbands ) 

!i = modulo(k-1,2)+1
!nsys = n

!if ( k .gt. 2 ) then
!nsys = n + ndim
!endif

!print*, i,j,k,n

Ham(j+1 ,j+1 ) = 2*Eg(k)  - elec*(2.d0*Kbs + 2.d0*Kcs)
Ham(j+2 ,j+2 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2.d0*Kcs - Dxf + Dso1/3.d0 + Kpp)
Ham(j+3 ,j+3 ) = 2*Eg(k)  - elec*(2.d0*Kas + Kbs + Kcs - Dxf + Dso1/3.d0 + Kpp)
Ham(j+4 ,j+4 ) = 2*Eg(k)  - elec*(2.d0*Kas + 2.d0*Kbs)
Ham(j+5 ,j+5 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2.d0*Kcs - Dxf - Kpp)
Ham(j+6 ,j+6 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2.d0*Kcs - Dxf + Kpp)
Ham(j+7 ,j+7 ) = 2*Eg(k)  - elec*(Kas + 2.d0*Kbs + Kcs + Kpp)
Ham(j+8 ,j+8 ) = 2*Eg(k)  - elec*(Kas + 2.d0*Kbs + Kcs + Kpp)
Ham(j+9 ,j+9 ) = 2*Eg(k)  - elec*(2.d0*Kas + Kbs + Kcs - Dxf - Kpp)
Ham(j+10,j+10) = 2*Eg(k)  - elec*(2.d0*Kas + Kbs + Kcs - Dxf + Kpp)
Ham(j+11,j+11) = 2*Eg(k)  - elec*(Kas + Kbs + 2.d0*Kcs - Dxf - Dso1/3.d0 + Kpp)
Ham(j+12,j+12) = 2*Eg(k)  - elec*(2.d0*Kas + 2.d0*Kcs - 2.d0*Dxf)
Ham(j+13,j+13) = 2*Eg(k)  - elec*(Kas + 2.d0*Kbs + Kcs - Kpp)
Ham(j+14,j+14) = 2*Eg(k)  - elec*(Kas + 2.d0*Kbs + Kcs + Kpp)
Ham(j+15,j+15) = 2*Eg(k)  - elec*(2.d0*Kas + Kbs + Kcs - Dxf - Dso1/3.d0 + Kpp)

Ham(j+1 ,j+4 ) = elec*Kpp
Ham(j+4 ,j+1 ) = Ham(j+1,j+4)

Ham(j+1 ,j+12 ) = elec*Kpp
Ham(j+12 ,j+1 ) = Ham(j+1,j+12)

Ham(j+4 ,j+12 ) = elec*Kpp
Ham(j+12 ,j+4 ) = Ham(j+4,j+12)

Ham(j+1 ,j+2 ) = elec*sqrt(2.d0)*Dso1/3.d0
Ham(j+2 ,j+1 ) = Ham(j+1,j+2)
Ham(j+3 ,j+4 ) = Ham(j+1,j+2)
Ham(j+4 ,j+3 ) = Ham(j+1,j+2)
Ham(j+11,j+12) = Ham(j+1,j+2)
Ham(j+12,j+11) = Ham(j+1,j+2)
Ham(j+12,j+15) = Ham(j+1,j+2)
Ham(j+15,j+12) = Ham(j+1,j+2)

Ham(j+5 ,j+6 ) = elec*Dso1/3.d0
Ham(j+6 ,j+5 ) = Ham(j+5,j+6)
Ham(j+6 ,j+7 ) = Ham(j+5,j+6)
Ham(j+7 ,j+6 ) = Ham(j+5,j+6)
Ham(j+8 ,j+10) = Ham(j+5,j+6)
Ham(j+10,j+8 ) = Ham(j+5,j+6)
Ham(j+9 ,j+10) = Ham(j+5,j+6)
Ham(j+10,j+9 ) = Ham(j+5,j+6)
Ham(j+11,j+13) = Ham(j+5,j+6)
Ham(j+13,j+11) = Ham(j+5,j+6)
Ham(j+11,j+14) = Ham(j+5,j+6)
Ham(j+14,j+11) = Ham(j+5,j+6)
Ham(j+15,j+13) = Ham(j+5,j+6)
Ham(j+13,j+15) = Ham(j+5,j+6)
Ham(j+14,j+15) = Ham(j+5,j+6)
Ham(j+15,j+14) = Ham(j+5,j+6)

Ham(j+5 ,j+7 ) = elec*Dso1/(-3.d0)
Ham(j+7 ,j+5 ) = Ham(j+5 ,j+7 )
Ham(j+8 ,j+9 ) = Ham(j+5 ,j+7 )
Ham(j+9 ,j+8 ) = Ham(j+5 ,j+7 )

Ham(j+13,j+14) = elec*2.d0*Dso1/3.d0
Ham(j+14,j+13) = Ham(j+13,j+14)

TransHam(0 ,j+1 ) = abs(TDM(k))/sqrt(6.d0)
TransHam(0 ,j+4 ) = abs(TDM(k))/sqrt(6.d0)
TransHam(0 ,j+5 ) = abs(TDM(k))/sqrt(6.d0)
TransHam(0 ,j+9 ) = abs(TDM(k))/sqrt(6.d0)
TransHam(0 ,j+12) = abs(TDM(k))/sqrt(6.d0)
TransHam(0 ,j+13) = abs(TDM(k))/sqrt(6.d0)

k = k + 1
j = j + 15

enddo

do i=1,24
do j=i+1,24
if ( mul(i) .eq. mul(j) ) then
TransHam(i,j) = abs(TransDip_Ana_h1h2(n))/sqrt(276.)
endif
enddo
enddo

Ham = Ham/Energ_au
TransHam = TransHam/D_to_au

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

end subroutine make_Ham_FS_FO

subroutine make_Ham_FS_FO_P

Eg(1) = minEe(1,n) + minEh(1,n)  + V0 - Eeh1(n)*elec
!Eg(2) = Eg(1) + elec*(-0.523865d0 + Eg(1)/elec * 0.293665d0)
!Eg(3) = minEe(1,n+ndim) + minEh(1,n+ndim)  + V0 - Eeh1(n)*elec
!Eg(4) = Eg(3) + elec*(-0.523865d0 + Eg(3)/elec * 0.293665d0)

TDM(1) = TransDip_Ana_h1e(n)
!TDM(2) = TransDip_Ana_h2e(n)
!TDM(3) = TransDip_Ana_h1e(n+ndim)
!TDM(4) = TransDip_Ana_h2e(n+ndim)

!print*, Eg(1), Eg(2), Eg(3), Eg(4)

j = 0
k = 1
i = 0

!if ( n .le. nQDA+nQDB ) then
!nbands = 2
!elseif ( n .gt. nQDA+nQDB ) then
!nbands = 4 
!endif

!do while ( k .le. nbands ) 

!i = modulo(k-1,2)+1
!nsys = n
!
!if ( k .gt. 2 ) then
!nsys = n + ndim
!endif

!print*, i,j,k,n

Ham(j+1 ,j+1 ) = Eg(k)  - elec* (2.d0*Dsop/3.d0 + 3.d0*Kp)
Ham(j+2 ,j+2 ) = Eg(k)  - elec* (2.d0*Dsop/3.d0 + 3.d0*Kp)
Ham(j+3 ,j+3 ) = Eg(k)  - elec* Kpp
Ham(j+4 ,j+4 ) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+5 ,j+5 ) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+6 ,j+6 ) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+7 ,j+7 ) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+8 ,j+8 ) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+9 ,j+9 ) = Eg(k)  - elec* Kpp
Ham(j+10,j+10) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+11,j+11) = Eg(k)  - elec* (3.d0*Kpp - 2.d0*Dsop/3.d0)
Ham(j+12,j+12) = Eg(k)  - elec* Kpp
Ham(j+13,j+13) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+14,j+14) = Eg(k)  - elec* Kpp
Ham(j+15,j+15) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+16,j+16) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+17,j+17) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+18,j+18) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+19,j+19) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+20,j+20) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+21,j+21) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+22,j+22) = Eg(k)  - elec* Kcs
Ham(j+23,j+23) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+24,j+24) = Eg(k)  - elec* Kpp
Ham(j+25,j+25) = Eg(k)  - elec* 3.d0*Kpp
Ham(j+26,j+26) = Eg(k)  - elec* (3.d0*Kpp - 2*Dsop/3.d0)
Ham(j+27,j+27) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+28,j+28) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+29,j+29) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+30,j+30) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+31,j+31) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+32,j+32) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+33,j+33) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+34,j+34) = Eg(k)  - elec* (3.d0*Kpp - Dxf)
Ham(j+35,j+35) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)
Ham(j+36,j+36) = Eg(k)  - elec* (3.d0*Kpp - Dsop/3.d0)

Ham(j+3 ,j+4 ) = elec*(-1.d0*Dso1/3.d0)
Ham(j+4 ,j+3 ) = Ham(j+3,j+4)
Ham(j+3 ,j+5 ) = Ham(j+3,j+4)
Ham(j+5 ,j+3 ) = Ham(j+3,j+4)
Ham(j+6 ,j+7 ) = Ham(j+3,j+4)
Ham(j+7 ,j+6 ) = Ham(j+3,j+4)
Ham(j+6 ,j+8 ) = Ham(j+3,j+4)
Ham(j+8 ,j+6 ) = Ham(j+3,j+4)
Ham(j+9 ,j+10) = Ham(j+3,j+4)
Ham(j+10,j+9 ) = Ham(j+3,j+4)
Ham(j+9 ,j+11) = Ham(j+3,j+4)
Ham(j+11,j+9 ) = Ham(j+3,j+4)
Ham(j+10,j+12) = Ham(j+3,j+4)
Ham(j+12,j+10) = Ham(j+3,j+4)

Ham(j+4 ,j+5 ) = elec*(Dso1/3.d0)
Ham(j+5 ,j+4 ) = Ham(j+4,j+5)
Ham(j+7 ,j+8 ) = Ham(j+4,j+5)
Ham(j+8 ,j+7 ) = Ham(j+4,j+5)
Ham(j+11,j+12) = Ham(j+4,j+5)
Ham(j+12,j+11) = Ham(j+4,j+5)

TransHam(0 ,j+3 ) = abs(TDM(k))/sqrt(3.d0) 
TransHam(0 ,j+7 ) = abs(TDM(k))/sqrt(3.d0)
TransHam(0 ,j+10) = abs(TDM(k))/sqrt(3.d0)

!print*, TransHam(j+0 ,j+3 )

!k = k + 1
!j = j + 12

!enddo

Ham = Ham/Energ_au

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_FS_FO_P

subroutine get_phases
l1( 1 )=   1.0e0_dp ; l2( 1 )=  0.0e0_dp ; l3( 1 )=  0.0e0_dp 
l1( 2 )=  -1.0e0_dp ; l2( 2 )=  0.0e0_dp ; l3( 2 )=  0.0e0_dp 
l1( 3 )=   0.0e0_dp ; l2( 3 )=  1.0e0_dp ; l3( 3 )=  0.0e0_dp 
l1( 4 )=   0.0e0_dp ; l2( 4 )= -1.0e0_dp ; l3( 4 )=  0.0e0_dp 
l1( 5 )=   0.0e0_dp ; l2( 5 )=  0.0e0_dp ; l3( 5 )=  1.0e0_dp 
l1( 6 )=   0.0e0_dp ; l2( 6 )=  0.0e0_dp ; l3( 6 )= -1.0e0_dp 
l1( 7 )=   3.0e0_dp ; l2( 7 )=  0.0e0_dp ; l3( 7 )=  0.0e0_dp 
l1( 8 )=  -3.0e0_dp ; l2( 8 )=  0.0e0_dp ; l3( 8 )=  0.0e0_dp 
l1( 9 )=   0.0e0_dp ; l2( 9 )=  3.0e0_dp ; l3( 9 )=  0.0e0_dp 
l1(10 )=   0.0e0_dp ; l2(10 )= -3.0e0_dp ; l3(10 )=  0.0e0_dp 
l1(11 )=   0.0e0_dp ; l2(11 )=  0.0e0_dp ; l3(11 )=  3.0e0_dp 
l1(12 )=   0.0e0_dp ; l2(12 )=  0.0e0_dp ; l3(12 )= -3.0e0_dp 
l1(13 )=   2.0e0_dp ; l2(13 )=  1.0e0_dp ; l3(13 )=  0.0e0_dp 
l1(14 )=  -2.0e0_dp ; l2(14 )= -1.0e0_dp ; l3(14 )=  0.0e0_dp 
l1(15 )=   2.0e0_dp ; l2(15 )= -1.0e0_dp ; l3(15 )=  0.0e0_dp 
l1(16 )=  -2.0e0_dp ; l2(16 )=  1.0e0_dp ; l3(16 )=  0.0e0_dp 
l1(17 )=   2.0e0_dp ; l2(17 )=  0.0e0_dp ; l3(17 )=  1.0e0_dp 
l1(18 )=  -2.0e0_dp ; l2(18 )=  0.0e0_dp ; l3(18 )= -1.0e0_dp 
l1(19 )=   2.0e0_dp ; l2(19 )=  0.0e0_dp ; l3(19 )= -1.0e0_dp 
l1(20 )=  -2.0e0_dp ; l2(20 )=  0.0e0_dp ; l3(20 )=  1.0e0_dp 
l1(21 )=   1.0e0_dp ; l2(21 )=  2.0e0_dp ; l3(21 )=  0.0e0_dp 
l1(22 )=  -1.0e0_dp ; l2(22 )= -2.0e0_dp ; l3(22 )=  0.0e0_dp 
l1(23 )=   1.0e0_dp ; l2(23 )= -2.0e0_dp ; l3(23 )=  0.0e0_dp 
l1(24 )=  -1.0e0_dp ; l2(24 )=  2.0e0_dp ; l3(24 )=  0.0e0_dp 
l1(25 )=   0.0e0_dp ; l2(25 )=  2.0e0_dp ; l3(25 )=  1.0e0_dp 
l1(26 )=   0.0e0_dp ; l2(26 )= -2.0e0_dp ; l3(26 )= -1.0e0_dp 
l1(27 )=   0.0e0_dp ; l2(27 )=  2.0e0_dp ; l3(27 )= -1.0e0_dp 
l1(28 )=   0.0e0_dp ; l2(28 )= -2.0e0_dp ; l3(28 )=  1.0e0_dp 
l1(29 )=   0.0e0_dp ; l2(29 )=  1.0e0_dp ; l3(29 )=  2.0e0_dp 
l1(30 )=   0.0e0_dp ; l2(30 )= -1.0e0_dp ; l3(30 )= -2.0e0_dp 
l1(31 )=   0.0e0_dp ; l2(31 )=  1.0e0_dp ; l3(31 )= -2.0e0_dp 
l1(32 )=   0.0e0_dp ; l2(32 )= -1.0e0_dp ; l3(32 )=  2.0e0_dp 
l1(33 )=   1.0e0_dp ; l2(33 )=  0.0e0_dp ; l3(33 )=  2.0e0_dp 
l1(34 )=  -1.0e0_dp ; l2(34 )=  0.0e0_dp ; l3(34 )= -2.0e0_dp 
l1(35 )=   1.0e0_dp ; l2(35 )=  0.0e0_dp ; l3(35 )= -2.0e0_dp 
l1(36 )=  -1.0e0_dp ; l2(36 )=  0.0e0_dp ; l3(36 )=  2.0e0_dp 
l1(37 )=   1.0e0_dp ; l2(37 )=  1.0e0_dp ; l3(37 )= -1.0e0_dp 
l1(38 )=  -1.0e0_dp ; l2(38 )= -1.0e0_dp ; l3(38 )=  1.0e0_dp 
l1(39 )=   1.0e0_dp ; l2(39 )= -1.0e0_dp ; l3(39 )=  1.0e0_dp 
l1(40 )=  -1.0e0_dp ; l2(40 )=  1.0e0_dp ; l3(40 )= -1.0e0_dp 
l1(41 )=  -1.0e0_dp ; l2(41 )=  1.0e0_dp ; l3(41 )=  1.0e0_dp 
l1(42 )=   1.0e0_dp ; l2(42 )= -1.0e0_dp ; l3(42 )= -1.0e0_dp 
l1(43 )=   1.0e0_dp ; l2(43 )=  1.0e0_dp ; l3(43 )=  1.0e0_dp 
l1(44 )=  -1.0e0_dp ; l2(44 )= -1.0e0_dp ; l3(44 )= -1.0e0_dp 
l1(45 )=   5.0e0_dp ; l2(45 )=  0.0e0_dp ; l3(45 )=  0.0e0_dp 
l1(46 )=  -5.0e0_dp ; l2(46 )=  0.0e0_dp ; l3(46 )=  0.0e0_dp
l1(47 )=   0.0e0_dp ; l2(47 )=  5.0e0_dp ; l3(47 )=  0.0e0_dp
l1(48 )=   0.0e0_dp ; l2(48 )= -5.0e0_dp ; l3(48 )=  0.0e0_dp
l1(49 )=   0.0e0_dp ; l2(49 )=  0.0e0_dp ; l3(49 )=  5.0e0_dp
l1(50 )=   0.0e0_dp ; l2(50 )=  0.0e0_dp ; l3(50 )= -5.0e0_dp
l1(51 )=   3.0e0_dp ; l2(51 )=  2.0e0_dp ; l3(51 )=  0.0e0_dp
l1(52 )=  -3.0e0_dp ; l2(52 )= -2.0e0_dp ; l3(52 )=  0.0e0_dp
l1(53 )=   3.0e0_dp ; l2(53 )= -2.0e0_dp ; l3(53 )=  0.0e0_dp
l1(54 )=  -3.0e0_dp ; l2(54 )=  2.0e0_dp ; l3(54 )=  0.0e0_dp
l1(55 )=   2.0e0_dp ; l2(55 )=  3.0e0_dp ; l3(55 )=  0.0e0_dp
l1(56 )=  -2.0e0_dp ; l2(56 )= -3.0e0_dp ; l3(56 )=  0.0e0_dp
l1(57 )=   2.0e0_dp ; l2(57 )= -3.0e0_dp ; l3(57 )=  0.0e0_dp
l1(58 )=  -2.0e0_dp ; l2(58 )=  3.0e0_dp ; l3(58 )=  0.0e0_dp
l1(59 )=   0.0e0_dp ; l2(59 )=  3.0e0_dp ; l3(59 )=  2.0e0_dp
l1(60 )=   0.0e0_dp ; l2(60 )= -3.0e0_dp ; l3(60 )= -2.0e0_dp
l1(61 )=   0.0e0_dp ; l2(61 )=  3.0e0_dp ; l3(61 )= -2.0e0_dp
l1(62 )=   0.0e0_dp ; l2(62 )= -3.0e0_dp ; l3(62 )=  2.0e0_dp
l1(63 )=   0.0e0_dp ; l2(63 )=  2.0e0_dp ; l3(63 )=  3.0e0_dp
l1(64 )=   0.0e0_dp ; l2(64 )= -2.0e0_dp ; l3(64 )= -3.0e0_dp
l1(65 )=   0.0e0_dp ; l2(65 )=  2.0e0_dp ; l3(65 )= -3.0e0_dp
l1(66 )=   0.0e0_dp ; l2(66 )= -2.0e0_dp ; l3(66 )=  3.0e0_dp
l1(67 )=   3.0e0_dp ; l2(67 )=  0.0e0_dp ; l3(67 )=  2.0e0_dp
l1(68 )=  -3.0e0_dp ; l2(68 )=  0.0e0_dp ; l3(68 )= -2.0e0_dp
l1(69 )=   3.0e0_dp ; l2(69 )=  0.0e0_dp ; l3(69 )= -2.0e0_dp
l1(70 )=  -3.0e0_dp ; l2(70 )=  0.0e0_dp ; l3(70 )=  2.0e0_dp
l1(71 )=   2.0e0_dp ; l2(71 )=  0.0e0_dp ; l3(71 )=  3.0e0_dp
l1(72 )=  -2.0e0_dp ; l2(72 )=  0.0e0_dp ; l3(72 )= -3.0e0_dp
l1(73 )=   2.0e0_dp ; l2(73 )=  0.0e0_dp ; l3(73 )= -3.0e0_dp
l1(74 )=  -2.0e0_dp ; l2(74 )=  0.0e0_dp ; l3(74 )=  3.0e0_dp
l1(75 )=   3.0e0_dp ; l2(75 )=  1.0e0_dp ; l3(75 )=  1.0e0_dp
l1(76 )=  -3.0e0_dp ; l2(76 )= -1.0e0_dp ; l3(76 )= -1.0e0_dp
l1(77 )=   3.0e0_dp ; l2(77 )=  1.0e0_dp ; l3(77 )= -1.0e0_dp
l1(78 )=  -3.0e0_dp ; l2(78 )= -1.0e0_dp ; l3(78 )=  1.0e0_dp
l1(79 )=   3.0e0_dp ; l2(79 )= -1.0e0_dp ; l3(79 )=  1.0e0_dp 
l1(80 )=  -3.0e0_dp ; l2(80 )=  1.0e0_dp ; l3(80 )= -1.0e0_dp
l1(81 )=   3.0e0_dp ; l2(81 )= -1.0e0_dp ; l3(81 )= -1.0e0_dp
l1(82 )=  -3.0e0_dp ; l2(82 )=  1.0e0_dp ; l3(82 )=  1.0e0_dp
l1(83 )=   1.0e0_dp ; l2(83 )=  3.0e0_dp ; l3(83 )=  1.0e0_dp
l1(84 )=  -1.0e0_dp ; l2(84 )= -3.0e0_dp ; l3(84 )= -1.0e0_dp
l1(85 )=   1.0e0_dp ; l2(85 )=  3.0e0_dp ; l3(85 )= -1.0e0_dp
l1(86 )=  -1.0e0_dp ; l2(86 )= -3.0e0_dp ; l3(86 )=  1.0e0_dp
l1(87 )=  -1.0e0_dp ; l2(87 )=  3.0e0_dp ; l3(87 )=  1.0e0_dp
l1(88 )=   1.0e0_dp ; l2(88 )= -3.0e0_dp ; l3(88 )= -1.0e0_dp
l1(89 )=   1.0e0_dp ; l2(89 )= -3.0e0_dp ; l3(89 )=  1.0e0_dp
l1(90 )=  -1.0e0_dp ; l2(90 )=  3.0e0_dp ; l3(90 )= -1.0e0_dp
l1(91 )=   1.0e0_dp ; l2(91 )=  1.0e0_dp ; l3(91 )=  3.0e0_dp
l1(92 )=  -1.0e0_dp ; l2(92 )= -1.0e0_dp ; l3(92 )= -3.0e0_dp
l1(93 )=   1.0e0_dp ; l2(93 )= -1.0e0_dp ; l3(93 )=  3.0e0_dp
l1(94 )=  -1.0e0_dp ; l2(94 )=  1.0e0_dp ; l3(94 )= -3.0e0_dp
l1(95 )=  -1.0e0_dp ; l2(95 )=  1.0e0_dp ; l3(95 )=  3.0e0_dp
l1(96 )=   1.0e0_dp ; l2(96 )= -1.0e0_dp ; l3(96 )= -3.0e0_dp
l1(97 )=   1.0e0_dp ; l2(97 )=  1.0e0_dp ; l3(97 )= -3.0e0_dp
l1(98 )=  -1.0e0_dp ; l2(98 )= -1.0e0_dp ; l3(98 )=  3.0e0_dp
l1(99 )=   2.0e0_dp ; l2(99 )=  2.0e0_dp ; l3(99 )=  1.0e0_dp
l1(100)=  -2.0e0_dp ; l2(100)= -2.0e0_dp ; l3(100)= -1.0e0_dp
l1(101)=  -2.0e0_dp ; l2(101)=  2.0e0_dp ; l3(101)=  1.0e0_dp
l1(102)=   2.0e0_dp ; l2(102)= -2.0e0_dp ; l3(102)= -1.0e0_dp
l1(103)=   2.0e0_dp ; l2(103)= -2.0e0_dp ; l3(103)=  1.0e0_dp
l1(104)=  -2.0e0_dp ; l2(104)=  2.0e0_dp ; l3(104)= -1.0e0_dp
l1(105)=   2.0e0_dp ; l2(105)=  2.0e0_dp ; l3(105)= -1.0e0_dp
l1(106)=  -2.0e0_dp ; l2(106)= -2.0e0_dp ; l3(106)=  1.0e0_dp
l1(107)=   2.0e0_dp ; l2(107)=  1.0e0_dp ; l3(107)=  2.0e0_dp
l1(108)=  -2.0e0_dp ; l2(108)= -1.0e0_dp ; l3(108)= -2.0e0_dp
l1(109)=  -2.0e0_dp ; l2(109)=  1.0e0_dp ; l3(109)=  2.0e0_dp
l1(110)=   2.0e0_dp ; l2(110)= -1.0e0_dp ; l3(110)= -2.0e0_dp
l1(111)=   2.0e0_dp ; l2(111)=  1.0e0_dp ; l3(111)= -2.0e0_dp
l1(112)=  -2.0e0_dp ; l2(112)= -1.0e0_dp ; l3(112)=  2.0e0_dp
l1(113)=   2.0e0_dp ; l2(113)= -1.0e0_dp ; l3(113)=  2.0e0_dp
l1(114)=  -2.0e0_dp ; l2(114)=  1.0e0_dp ; l3(114)= -2.0e0_dp
l1(115)=   1.0e0_dp ; l2(115)=  2.0e0_dp ; l3(115)=  2.0e0_dp
l1(116)=  -1.0e0_dp ; l2(116)= -2.0e0_dp ; l3(116)= -2.0e0_dp
l1(117)=   1.0e0_dp ; l2(117)= -2.0e0_dp ; l3(117)=  2.0e0_dp
l1(118)=  -1.0e0_dp ; l2(118)=  2.0e0_dp ; l3(118)= -2.0e0_dp
l1(119)=   1.0e0_dp ; l2(119)=  2.0e0_dp ; l3(119)= -2.0e0_dp
l1(120)=  -1.0e0_dp ; l2(120)= -2.0e0_dp ; l3(120)=  2.0e0_dp
l1(121)=  -1.0e0_dp ; l2(121)=  2.0e0_dp ; l3(121)=  2.0e0_dp
l1(122)=   1.0e0_dp ; l2(122)= -2.0e0_dp ; l3(122)= -2.0e0_dp
l1(123)=   4.0e0_dp ; l2(123)=  1.0e0_dp ; l3(123)=  0.0e0_dp
l1(124)=  -4.0e0_dp ; l2(124)= -1.0e0_dp ; l3(124)=  0.0e0_dp 
l1(125)=   4.0e0_dp ; l2(125)= -1.0e0_dp ; l3(125)=  0.0e0_dp
l1(126)=  -4.0e0_dp ; l2(126)=  1.0e0_dp ; l3(126)=  0.0e0_dp
l1(127)=   1.0e0_dp ; l2(127)=  4.0e0_dp ; l3(127)=  0.0e0_dp
l1(128)=  -1.0e0_dp ; l2(128)= -4.0e0_dp ; l3(128)=  0.0e0_dp
l1(129)=  -1.0e0_dp ; l2(129)=  4.0e0_dp ; l3(129)=  0.0e0_dp
l1(130)=   1.0e0_dp ; l2(130)= -4.0e0_dp ; l3(130)=  0.0e0_dp
l1(131)=   1.0e0_dp ; l2(131)=  0.0e0_dp ; l3(131)=  4.0e0_dp
l1(132)=  -1.0e0_dp ; l2(132)=  0.0e0_dp ; l3(132)= -4.0e0_dp
l1(133)=  -1.0e0_dp ; l2(133)=  0.0e0_dp ; l3(133)=  4.0e0_dp
l1(134)=   1.0e0_dp ; l2(134)=  0.0e0_dp ; l3(134)= -4.0e0_dp
l1(135)=   4.0e0_dp ; l2(135)=  0.0e0_dp ; l3(135)=  1.0e0_dp
l1(136)=  -4.0e0_dp ; l2(136)=  0.0e0_dp ; l3(136)= -1.0e0_dp
l1(137)=   4.0e0_dp ; l2(137)=  0.0e0_dp ; l3(137)= -1.0e0_dp
l1(138)=  -4.0e0_dp ; l2(138)=  0.0e0_dp ; l3(138)=  1.0e0_dp
l1(139)=   0.0e0_dp ; l2(139)=  4.0e0_dp ; l3(139)=  1.0e0_dp
l1(140)=   0.0e0_dp ; l2(140)= -4.0e0_dp ; l3(140)= -1.0e0_dp
l1(141)=   0.0e0_dp ; l2(141)=  4.0e0_dp ; l3(141)= -1.0e0_dp
l1(142)=   0.0e0_dp ; l2(142)= -4.0e0_dp ; l3(142)=  1.0e0_dp
l1(143)=   0.0e0_dp ; l2(143)=  1.0e0_dp ; l3(143)=  4.0e0_dp
l1(144)=   0.0e0_dp ; l2(144)= -1.0e0_dp ; l3(144)= -4.0e0_dp
l1(145)=   0.0e0_dp ; l2(145)= -1.0e0_dp ; l3(145)=  4.0e0_dp
l1(146)=   0.0e0_dp ; l2(146)=  1.0e0_dp ; l3(146)= -4.0e0_dp
end subroutine

subroutine make_eTDM

eTDM(0,1,:) = vector(1._dp) 
eTDM(0,2,:) = vector(1._dp)
eTDM(0,3,:) = vector(1._dp)
eTDM(0,4,:) = vector(1._dp)
eTDM(0,5,:) = vector(1._dp)
eTDM(0,6,:) = vector(1._dp)
eTDM(0,7,:) = vector(1._dp)
eTDM(0,8,:) = vector(1._dp)
eTDM(1,2,:) = vector(1._dp)
eTDM(1,5,:) = vector(1._dp)
eTDM(1,7,:) = vector(1._dp)
eTDM(1,8,:) = vector(1._dp)
eTDM(2,6,:) = vector(1._dp)
eTDM(2,7,:) = vector(1._dp)
eTDM(2,8,:) = vector(1._dp)
eTDM(3,4,:) = vector(1._dp)
eTDM(3,5,:) = vector(1._dp)
eTDM(3,6,:) = vector(1._dp)
eTDM(3,7,:) = vector(1._dp)
eTDM(4,5,:) = vector(1._dp)
eTDM(4,6,:) = vector(1._dp)
eTDM(4,8,:) = vector(1._dp)
eTDM(5,6,:) = vector(1._dp)
eTDM(7,8,:) = vector(1._dp)

end subroutine

subroutine make_distMat

open(newunit=getMat,file='getDmat.sh')

write(getMat,*) "#!/bin/bash"
write(getMat,*) 
write(getMat,*) "cat > plot-distrib-tmp <<EOF"
write(getMat,*) 
write(getMat,*) "width = 0.002"
write(getMat,*) "set boxwidth width absolute"
write(getMat,*) "set style fill solid 1.0"
write(getMat,*) "bin_width = width;"
write(getMat,*) "bin_number(x) = floor(x/bin_width)"
write(getMat,*) "rounded(x) = bin_width * ( bin_number(x) )"
write(getMat,*) "set table 'test1.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$4)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test2.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$5)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test3.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$6)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test4.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$7)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test5.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$8)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test6.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$9)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) "set table 'test7.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$10)):(1) smooth frequency"
write(getMat,*) "unset table"
write(getMat,*) "set table 'test8.dat'"
write(getMat,*) "plot 'Etransitions-he_ei.dat' using (rounded(\$11)):(1) smooth frequency "
write(getMat,*) "unset table"
write(getMat,*) 
write(getMat,*) "a2=2.5"
write(getMat,*) "b2=2.5" 
write(getMat,*) "c2=2.5"
write(getMat,*) "d2=2.5"
write(getMat,*) "e2=2.5"
write(getMat,*) "f2=2.5"
write(getMat,*) "g2=2.5"
write(getMat,*) "h2=2.5"
write(getMat,*) 
write(getMat,*) "f1(x) = (a1/(a3*(2*pi)**0.5))*exp(-(x-a2)**2/(2*(a3)**2))" 
write(getMat,*) "f2(x) = (b1/(b3*(2*pi)**0.5))*exp(-(x-b2)**2/(2*(b3)**2))"
write(getMat,*) "f3(x) = (c1/(c3*(2*pi)**0.5))*exp(-(x-c2)**2/(2*(c3)**2))"
write(getMat,*) "f4(x) = (d1/(d3*(2*pi)**0.5))*exp(-(x-d2)**2/(2*(d3)**2))"
write(getMat,*) "f5(x) = (e1/(e3*(2*pi)**0.5))*exp(-(x-e2)**2/(2*(e3)**2))"
write(getMat,*) "f6(x) = (f1/(f3*(2*pi)**0.5))*exp(-(x-f2)**2/(2*(f3)**2))"
write(getMat,*) "f7(x) = (g1/(g3*(2*pi)**0.5))*exp(-(x-g2)**2/(2*(g3)**2))"
write(getMat,*) "f8(x) = (h1/(h3*(2*pi)**0.5))*exp(-(x-h2)**2/(2*(h3)**2))"
write(getMat,*) 
write(getMat,*) "fit f1(x) 'test1.dat' via a1,a2,a3"
write(getMat,*) "save fit 'fit1.dat'"
write(getMat,*) "fit f2(x) 'test2.dat' via b1,b2,b3"
write(getMat,*) "save fit 'fit2.dat'"
write(getMat,*) "fit f3(x) 'test3.dat' via c1,c2,c3"
write(getMat,*) "save fit 'fit3.dat'"
write(getMat,*) "fit f4(x) 'test4.dat' via d1,d2,d3"
write(getMat,*) "save fit 'fit4.dat'"
write(getMat,*) "fit f5(x) 'test5.dat' via e1,e2,e3"
write(getMat,*) "save fit 'fit5.dat'"
write(getMat,*) "fit f6(x) 'test6.dat' via f1,f2,f3"
write(getMat,*) "save fit 'fit6.dat'"
write(getMat,*) "fit f7(x) 'test7.dat' via g1,g2,g3"
write(getMat,*) "save fit 'fit7.dat'"
write(getMat,*) "fit f8(x) 'test8.dat' via h1,h2,h3"
write(getMat,*) "save fit 'fit8.dat'"
write(getMat,*) 
write(getMat,*) "EOF"
write(getMat,*) 
write(getMat,*) "gnuplot plot-distrib-tmp 2> tmp"
write(getMat,*) 
write(getMat,*) "rm tmp2"
write(getMat,*) 
write(getMat,*) "grep 'a3              =' fit1.dat | awk '{print $3}' >> tmp2" 
write(getMat,*) "grep 'b3              =' fit2.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'c3              =' fit3.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'd3              =' fit4.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'e3              =' fit5.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'f3              =' fit6.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'g3              =' fit7.dat | awk '{print $3}' >> tmp2"
write(getMat,*) "grep 'h3              =' fit8.dat | awk '{print $3}' >> tmp2"

allocate(sigma(8))

call system('sh getDmat.sh')

open(newunit=sigfit,file='tmp2')
open(newunit=Dmat,file='Dmat.dat')

do i=1,8
read(sigfit,*) sigma(i)
enddo

write(Dmat,'(9f16.8)') 0._dp, (sigma(j), j=1,8)
do i=1,8
write(Dmat,'(9f16.8)') sigma(i), (abs(sigma(i)-sigma(j)), j=1,8)
enddo

write(Dmat,*) 

write(Dmat,'(9f16.8)') 0._dp, (1.e15_dp*6.582119570e-16_dp/sigma(j), j=1,8)
do i=1,8
write(Dmat,'(9f16.8)') 1.e15_dp*6.582119570e-16_dp/sigma(i), (1.e15_dp*6.582119570e-16_dp/abs(sigma(i)-sigma(j)), j=1,8)
enddo

end subroutine

end module Make_Ham
