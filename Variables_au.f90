module Variables_au

use omp_lib
use Constants_au
use Normal
use Vectors

implicit none

   character*5 :: vers
   character*6 :: pgeom
   character*64 :: popc, hmti, norm, tdmM, hmt0, outputdir, norm_ei, popc_ei, Re_c_ei, Im_c_ei, integ,model, line, dummy
   character*64 :: Re_c, Im_c, syst_n, Re_c_l, Im_c_L, cov, cov2, Pop_c_L
   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, dyn, hamilt, get_ei, finest, get_sp
   character*1 :: TDM_ee, Dyn_0, Dyn_ei, inbox, Dyn_L,doFT,CEP1,CEP2,CEP3,singleFT,nofiles, singleDS, doCovar,doFT_s,doAbs
   character*1 :: rdm_ori, noMat, Dyn_avg, doDmat
   character*1, allocatable :: mul(:)
   logical :: isit
   integer :: Pulse_f,Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,TransAbs, DipSpec_conv_f
   integer :: popc_0_f,popc_ei_f,norm_0_f,norm_ei_f,Re_c_ei_f,Im_c_ei_f,Re_c_0_f,Im_c_0_f,TDip_ei_f,tmp, nbands,Liou_f, getMat
   integer :: Re_c_L_f,Im_c_L_f,H_ei_f,Etr_0_f,Etr_ei_f,Abs_imp_f, t2, DipSpec, pol, npol, P_Match_f, DipSpec_R_f,DipSpec_NR_f
   character*64 :: form_mat,form_arr,form_abs,form_pop,form_com,form_TDM,form1,form_com_L, form_DipSpec, form_pop_L
   integer :: syst, ndots, n, rmin, rmax, nsys, npulses, nstates, ntime,i,j,k,l,t,lwork, info, idlink, threads, nQDA, nQDB
   integer :: nhomoA,nhomoB,nhetero,totsys,ndim,nQD, EminID,EmaxID,Estep,nmax,io,abso, kc, kl, nstates2, DipSpec_s, pol2
   integer :: TransAbs_NR, TransAbs_R, TransAbs_s, nFT, FTpow, t_ana, f_ana, Pop_c_L_f, TransAbs_P, ierr, liworku, lworku
   integer :: scos_ssin,TransAbs_7_f,TransAbs_17_f,TransAbs_33_f,TransAbs_39_f,TransAbs_41_f,TransAbs_44_f,DipSpec_7_f,DipSpec_17_f
   integer :: DipSpec_33_f,DipSpec_39_f,DipSpec_41_f,DipSpec_44_f, sphere, DipSpec_43_f, TransAbs_43_f
   integer :: Tmat_avg_f, Tmat_avgx_f, Tmat_avgy_f, Tmat_avgz_f, Etr_avg_f, popc_avg_f, Re_c_avg_f, Im_c_avg_f, TransAbs_avg
   integer :: DipSpec_avg, P_Match_avg, TransAbs_avg_7_f, TransAbs_avg_17_f, TransAbs_avg_33_f, TransAbs_avg_39_f, TransAbs_avg_41_f
   integer :: TransAbs_avg_43_f, TransAbs_avg_44_f, DipSpec_avg_7_f, DipSpec_avg_17_f, DipSpec_avg_33_f, DipSpec_avg_39_f
   integer :: DipSpec_avg_41_f, DipSpec_avg_43_f, DipSpec_avg_44_f, Dyn_avg_flag, error, label_0, sigfit, Dmat
   integer,allocatable :: seed(:),icol(:,:),irow(:,:),iwork2(:),zero(:)
   real(dp) :: a13_1d_he,a13_2d_he,a13_3d_he,a13_4d_he,a15_1d_he,a15_2d_he,a15_3d_he,a15_4d_he,a17_1d_he,a17_2d_he,a17_3d_he,&
               a17_4d_he,a24_1d_he,a24_2d_he,a24_3d_he,a24_4d_he,a26_1d_he,a26_2d_he,a26_3d_he,a26_4d_he,a28_1d_he,a28_2d_he,&
               a28_3d_he,a28_4d_he,a35_1d_he,a35_2d_he,a35_3d_he,a35_4d_he,a37_1d_he,a37_2d_he,a37_3d_he,a37_4d_he,a46_1d_he,&
               a46_2d_he,a46_3d_he,a46_4d_he,a48_1d_he,a48_2d_he,a48_3d_he,a48_4d_he,a55_1d_he,a55_2d_he,a55_3d_he,a55_4d_he,&
               a57_1d_he,a57_2d_he,a57_3d_he,a57_4d_he,a66_1d_he,a66_2d_he,a66_3d_he,a66_4d_he,a68_1d_he,a68_2d_he,a68_3d_he,&
               a68_4d_he,a77_1d_he,a77_2d_he,a77_3d_he,a77_4d_he,a88_1d_he,a88_2d_he,a88_3d_he,a88_4d_he,a14_1d_he,a14_2d_he,&
               a14_3d_he,a14_4d_he,a16_1d_he,a16_2d_he,a16_3d_he,a16_4d_he,a18_1d_he,a18_2d_he,a18_3d_he,a18_4d_he,a23_1d_he,&
               a23_2d_he,a23_3d_he,a23_4d_he,a25_1d_he,a25_2d_he,a25_3d_he,a25_4d_he,a27_1d_he,a27_2d_he,a27_3d_he,a27_4d_he,&
               a36_1d_he,a36_2d_he,a36_3d_he,a36_4d_he,a38_1d_he,a38_2d_he,a38_3d_he,a38_4d_he,a45_1d_he,a45_2d_he,a45_3d_he,&
               a45_4d_he,a47_1d_he,a47_2d_he,a47_3d_he,a47_4d_he,a56_1d_he,a56_2d_he,a56_3d_he,a56_4d_he,a58_1d_he,a58_2d_he,&
               a58_3d_he,a58_4d_he,a67_1d_he,a67_2d_he,a67_3d_he,a67_4d_he,a78_1d_he,a78_2d_he,a78_3d_he,a78_4d_he
   real(dp) :: a13_1e_he,a13_2e_he,a13_3e_he,a13_4e_he,a15_1e_he,a15_2e_he,a15_3e_he,a15_4e_he,a17_1e_he,a17_2e_he,a17_3e_he,&
               a17_4e_he,a24_1e_he,a24_2e_he,a24_3e_he,a24_4e_he,a26_1e_he,a26_2e_he,a26_3e_he,a26_4e_he,a28_1e_he,a28_2e_he,&
               a28_3e_he,a28_4e_he,a35_1e_he,a35_2e_he,a35_3e_he,a35_4e_he,a37_1e_he,a37_2e_he,a37_3e_he,a37_4e_he,a46_1e_he,&
               a46_2e_he,a46_3e_he,a46_4e_he,a48_1e_he,a48_2e_he,a48_3e_he,a48_4e_he,a55_1e_he,a55_2e_he,a55_3e_he,a55_4e_he,&
               a57_1e_he,a57_2e_he,a57_3e_he,a57_4e_he,a66_1e_he,a66_2e_he,a66_3e_he,a66_4e_he,a68_1e_he,a68_2e_he,a68_3e_he,&
               a68_4e_he,a77_1e_he,a77_2e_he,a77_3e_he,a77_4e_he,a88_1e_he,a88_2e_he,a88_3e_he,a88_4e_he,a14_1e_he,a14_2e_he,&
               a14_3e_he,a14_4e_he,a16_1e_he,a16_2e_he,a16_3e_he,a16_4e_he,a18_1e_he,a18_2e_he,a18_3e_he,a18_4e_he,a23_1e_he,&
               a23_2e_he,a23_3e_he,a23_4e_he,a25_1e_he,a25_2e_he,a25_3e_he,a25_4e_he,a27_1e_he,a27_2e_he,a27_3e_he,a27_4e_he,&
               a36_1e_he,a36_2e_he,a36_3e_he,a36_4e_he,a38_1e_he,a38_2e_he,a38_3e_he,a38_4e_he,a45_1e_he,a45_2e_he,a45_3e_he,&
               a45_4e_he,a47_1e_he,a47_2e_he,a47_3e_he,a47_4e_he,a56_1e_he,a56_2e_he,a56_3e_he,a56_4e_he,a58_1e_he,a58_2e_he,&
               a58_3e_he,a58_4e_he,a67_1e_he,a67_2e_he,a67_3e_he,a67_4e_he,a78_1e_he,a78_2e_he,a78_3e_he,a78_4e_he
   real(dp) :: tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8, time_ana, aR_avgA, aR_avgB, start, finish
   real(dp) :: tpx1,tpx2,tpx3,tpx4,tpx5,tpx6,tpx7,tpx8
   real(dp) :: tpy1,tpy2,tpy3,tpy4,tpy5,tpy6,tpy7,tpy8
   real(dp) :: tpz1,tpz2,tpz3,tpz4,tpz5,tpz6,tpz7,tpz8
   real(dp) :: aA, aB, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, slope, V0eV, minr, maxr, rsteps, side, link
   real(dp) :: sigma_conv,Emin,Emax,x,w1,w2,w,wstep,powtemp,time_FFT, FTscale, psirot, phirot, thetarot
   real(dp) :: vertex, zbase, alphae, alphah1, alphah2, betae, betah1, betah2, rdm1, rdm2
   real(dp) :: dispQD, displink, rdmlinker, rdmQDA, rdmQDB, t01, t02, t03, timestep, totaltime, distQD, Kpp, Dsop, Kp
   real(dp) :: omega01, omega02, omega03, phase01, phase02, phase03, width01, width02, width03, Ed01, Ed02, Ed03
   real(dp) :: pulse1, pulse2, pulse3, test, time, cnorm, cnormabs, cnormconj, cnorm2, cnorm2_ei, Kas, Kbs, Kcs, Dso1, Dso2, Dxf
   real(dp),allocatable :: aR(:), aRA(:), aRB(:), epsin(:), epsR(:), V0e(:), V0h(:), linker(:), Eg(:), TDM(:), eTDM(:,:,:), sigma(:)
   real(dp),allocatable :: epsinA(:), epsinB(:), epsRA(:), epsRB(:), V0eA(:), V0eB(:), V0hA(:), V0hB(:),l1(:), l2(:), l3(:)
   real(dp),allocatable :: Cb_eh1(:), Cb_eh2(:), Norm_Ana_e(:), Norm_Ana_h1(:), Norm_Ana_h2(:), Eeh1(:), Eeh2(:)
   real(dp),allocatable :: OverlapAna_h1e(:), OverlapAna_h2e(:), Cb_Num_eh1(:), Cb_Num_eh1_eh2(:), Cb_Num_eh2(:)
   real(dp),allocatable :: minEe(:,:),minEh(:,:), TransDip_Num_h1e(:), TransDip_Num_h2e(:), work(:), lambda(:)
   real(dp),allocatable :: TransDip_Ana_h1e(:), TransDip_Ana_h2e(:), Oscillator_Ana_h1e(:), Oscillator_Ana_h2e(:), Transvec(:)
   real(dp),allocatable :: ExctCoef_h1e(:), ExctCoef_h2e(:), Ham(:,:), E0(:), c0(:), TransHam(:,:), Hamt(:,:,:), c(:,:)
   real(dp),allocatable :: TransMat_ei(:,:), TransHam0(:,:), Ham_0(:), Ham_dir(:,:), Ham_ex(:,:), Ham_ei(:,:), Ham_l(:,:), haml(:,:)
   real(dp),allocatable :: TransDip_Ana_h1h2(:), TransHam_ei(:,:), Mat(:,:), QDcoor(:,:), Dcenter(:,:), Pe1(:), Pe2(:), Pe3(:)
   real(dp),allocatable :: TransHam_d(:,:,:), TransHam_l(:,:,:), TransHam_ei_l(:,:,:), k_1(:), k_2(:), k_3(:), work1(:), work2(:)
   real(dp),allocatable :: Matx(:,:), Maty(:,:), Matz(:,:),spec(:),dipole(:,:), Ham0_avg(:,:)
   real(dp),allocatable :: pow(:),pow_gaus(:),pulses(:), pow_s(:,:), pow_gaus_s(:,:), pulses_FFT(:)
   real(dp),allocatable :: Scov(:,:), pop(:,:),merge_diag(:,:),merge_odiag(:,:), scos(:), ssin(:), lfield(:,:)
   real(dp),allocatable :: TransHam_avg(:,:),TransHam_avg_l(:,:,:),lambda_avg(:), bloblo(:), blabla(:),maxid(:), rotmat(:,:)
 complex(8) :: ct1, ct2, ct3, ct4, xt01, xt02, xt03, xhbar, im, xwidth, xomega , xEd, xh, xphase, xtime, xhbar_au
 complex(8) :: integPol, integPol_diff, integPolconv, pow_41, pow_41_conv
 complex(8),allocatable :: xHam(:,:) , xHamt(:,:,:), xTransHam(:,:), xE0(:), xHamtk2(:,:,:), xHamtk3(:,:,:), xHamtk4(:,:,:)
 complex(8),allocatable :: xc0(:), xc(:,:), xc_ei(:,:), xcnew(:,:), k1(:), k2(:), k3(:) , k4(:), xHam_ei(:,:), pow_pol_conv(:)
 complex(8),allocatable :: k1_L(:), k2_L(:), k3_L(:) , k4_L(:),  k5_L(:), k6_L(:), k7_L(:) , k8_L(:), xpow_pol(:,:),xliou(:,:,:,:)
 complex(8),allocatable :: dk1(:), dk2(:), dk3(:) , dk4(:), k5(:), k6(:), k7(:) , k8(:), pow_pol(:,:), pow_pol_gaus(:,:)
 complex(8),allocatable :: xc_ei_av(:,:), xctemp(:),xlfield(:,:),xc_L(:,:), xpow_gaus(:),xpulse(:),wft(:),wftp(:),wftf(:),xhugo(:)
 complex(8),allocatable :: wft_pol(:,:),wftf_pol(:,:), pow_pol_diff(:), wft_s(:,:), wftf_s(:,:), xpow_gaus_s(:), xpulse2(:)
 complex(8),allocatable :: k1_rho(:,:), k2_rho(:,:), k3_rho(:,:), k4_rho(:,:), k5_rho(:,:), k6_rho(:,:), k7_rho(:,:), k8_rho(:,:)
 complex(8),allocatable :: xc_rho(:,:,:), wftf_t1(:)

contains 

subroutine getVariables

NAMELIST /outputs/   inbox,rdm_ori,get_sp,get_ei,Dyn_0,Dyn_ei,Dyn_L,Dyn_avg,TDM_ee,&
                     doAbs,doFT,singleDS,doFT_s,singleFT,nofiles,noMat,doCovar,doDmat
NAMELIST /elecSt/    model,me,mh,eps,epsout,V0eV,omegaLO,slope,side
NAMELIST /fineStruc/ Kas,Kbs,Kcs,Kpp,Dso1,Dso2,Dxf
NAMELIST /pulses/    integ,npulses,t01,t02,t03,timestep,totaltime,omega01,omega02,omega03,phase01,phase02,phase03,&
                     width01,width02,width03,Ed01,Ed02,Ed03,CEP1,CEP2,CEP3,pgeom,vertex
NAMELIST /syst/      nQDA,nQDB,nhomoA,nhomoB,nhetero,dispQD,idlink,aA,aB     
NAMELIST /FT/        FTpow

open(150,file='QD_quest.def',form='formatted')
read(150,NML=outputs)
read(150,NML=elecSt)
read(150,NML=fineStruc)
read(150,NML=pulses)
read(150,NML=syst)
read(150,NML=FT)

Dyn_avg_flag = 0
im         = dcmplx(0.0e0_dp,1.0e0_dp)
me         = me*m0
mh         = mh*m0
rhoe       = 1.0e0_dp/sqrt((2.e0_dp*me*omegaLO)/hbar)
rhoh       = 1.0e0_dp/sqrt((2.e0_dp*mh*omegaLO)/hbar)
V0         = V0eV*elec
!npol       = 44
npol       = 146
!dispQD = dispQD*0.9d0

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) .or. ( Dyn_L .eq. 'y' ) ) then

timestep   =  timestep*1.e-15_dp/t_au  !timestep*1.d-15/t_au
totaltime  =  totaltime*1.e-15_dp/t_au !totaltime*1.d-15/t_au
t01        =  t01*1.e-15_dp/t_au       !t01*1.d-15/t_au
t02        =  t02*1.e-15_dp/t_au       !t02*1.d-15/t_au
t03        =  t03*1.e-15_dp/t_au       !t03*1.d-15/t_au
width01    =  width01*1.e-15_dp/t_au     !width*1.d-15/t_au
width02    =  width02*1.e-15_dp/t_au     !width*1.d-15/t_au
width03    =  width03*1.e-15_dp/t_au     !width*1.d-15/t_au
omega01    =  omega01*t_au*(2.d0*pi)      !omega*1.d15*t_au
omega02    =  omega02*t_au*(2.d0*pi)      !omega*1.d15*t_au
omega03    =  omega03*t_au*(2.d0*pi)      !omega*1.d15*t_au
Ed01       =  Ed01/E_au        !0.024 !Ed/E_au
Ed02       =  Ed02/E_au        !0.024 !Ed/E_au
Ed03       =  Ed03/E_au        !0.024 !Ed/E_au
xh         =  dcmplx(timestep,0.0e0_dp)
ntime      =  nint(totaltime/timestep)

if ( CEP1 .eq. 'r' ) then 
call random_number(phase01) 
phase01      =  phase01 * 2.e0_dp * pi
elseif ( CEP1 .eq. 'p' ) then
phase01      =  pi
elseif ( CEP1 .eq. '0' ) then
phase01      =  0.e0_dp
endif
if ( CEP2 .eq. 'r' ) then 
call random_number(phase02) 
phase02      =  phase02 * 2.e0_dp * pi
elseif ( CEP2 .eq. 'p' ) then
phase02      =  pi
elseif ( CEP2 .eq. '0' ) then
phase02      =  0.e0_dp
endif
if ( CEP3 .eq. 'r' ) then 
call random_number(phase03) 
phase03      =  phase03 * 2.e0_dp * pi
elseif ( CEP3 .eq. 'p' ) then
phase03      =  pi
elseif ( CEP3 .eq. '0' ) then
phase03      =  0.e0_dp
endif

allocate(k_1(3),k_2(3),k_3(3),Pe1(3),Pe2(3),Pe3(3),source=0.e0_dp)

if ( npulses .eq. 3) then
pulse1 = 1.e0_dp 
pulse2 = 1.e0_dp
pulse3 = 1.e0_dp
elseif ( npulses .eq. 2) then
pulse1 = 1.e0_dp
pulse2 = 1.e0_dp
pulse3 = 0.e0_dp
elseif ( npulses .eq. 1) then
pulse1 = 1.e0_dp
pulse2 = 0.e0_dp
pulse3 = 0.e0_dp
elseif ( npulses .eq. 0) then
pulse1 = 0.e0_dp
pulse2 = 0.e0_dp
pulse3 = 0.e0_dp
endif

!if ( inbox .eq. 'y' ) then

if ( pgeom .eq. 'boxcar' ) then

zbase = 1000e0_dp*sqrt(2.e0_dp)*sqrt(1-cos(pi*vertex/180.e0_dp))/sqrt(cos(pi*vertex/180.e0_dp))
k_1(1) =(-2.e0_dp * pi / (cl/(omega01/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2&
+(zbase / 2.e0_dp)**2))
k_1(2) = (2.e0_dp * pi / (cl/(omega01/t_au))) * ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_1(3) =(-2.e0_dp * pi / (cl/(omega01/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_2(1) =(-2.e0_dp * pi / (cl/(omega02/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))  
k_2(2) = (2.e0_dp * pi / (cl/(omega02/t_au))) * ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_2(3) = (2.e0_dp * pi / (cl/(omega02/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_3(1) = (2.e0_dp * pi / (cl/(omega03/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_3(2) = (2.e0_dp * pi / (cl/(omega03/t_au))) * ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))
k_3(3) = (2.e0_dp * pi / (cl/(omega03/t_au))) * ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2& 
+(zbase / 2.e0_dp)**2))

Pe1(1) =  ( -zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe1(2) =  ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe1(3) =  ( -zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe2(1) =  ( -zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))  
Pe2(2) =  ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe2(3) =  ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe3(1) =  ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe3(2) =  ( 1000.e0_dp      ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))
Pe3(3) =  ( zbase / 2.e0_dp ) / (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase / 2.e0_dp)**2))

elseif ( pgeom .eq. 'triang' ) then

zbase= 1000e0_dp*sqrt(2.e0_dp*(1.e0_dp-cos(pi*vertex/180.e0_dp))/(1.e0_dp-(1-cos(pi*vertex/180.e0_dp))/&
                (1-cos(pi*120.e0_dp/180.e0_dp)))) 
k_1(1) = (2.e0_dp*pi/(cl/(omega01/t_au)))*(zbase / 2.e0_dp ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_1(2) = (2.e0_dp*pi/(cl/(omega01/t_au)))*(1000.e0_dp      ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_1(3) = (2.e0_dp*pi/(cl/(omega01/t_au)))*(zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_2(1) =(-2.e0_dp*pi/(cl/(omega02/t_au)))*(zbase / 2.e0_dp ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_2(2) = (2.e0_dp*pi/(cl/(omega02/t_au)))*(1000.e0_dp      ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_2(3) = (2.e0_dp*pi/(cl/(omega02/t_au)))*(zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+&
(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_3(1) = ( 2.e0_dp*pi/(cl/(omega03/t_au)))*0.e0_dp 
k_3(2) = ( 2.e0_dp*pi/(cl/(omega03/t_au)))*(1000.e0_dp      ) /        (sqrt((1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
k_3(3) =(-2.e0_dp*pi/(cl/(omega03/t_au)))*(zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))

Pe1(1) =  ( zbase / 2.e0_dp ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe1(2) =  ( 1000.e0_dp      ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe1(3) =  ( zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe2(1) =  ( zbase / 2.e0_dp ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe2(2) =  ( 1000.e0_dp      ) /        (sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe2(3) =  ( zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((zbase / 2.e0_dp)**2+(1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe3(1) =  0.e0_dp 
Pe3(2) =  ( 1000.e0_dp      ) /        (sqrt((1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))
Pe3(3) =  ( zbase/(2.e0_dp*sqrt(3.e0_dp)))/(sqrt((1000.e0_dp)**2+(zbase/(2.e0_dp*sqrt(3.e0_dp)))**2))

elseif ( pgeom .eq. 'straig' ) then

k_1(1) = 0.e0_dp 
k_1(2) = 0.e0_dp
k_1(3) = 0.e0_dp
k_2(1) = 0.e0_dp
k_2(2) = 0.e0_dp
k_2(3) = 0.e0_dp
k_3(1) = 0.e0_dp
k_3(2) = 0.e0_dp
k_3(3) = 0.e0_dp

Pe1(1) = 1.e0_dp/sqrt(3.e0_dp) 
Pe1(2) = 1.e0_dp/sqrt(3.e0_dp)
Pe1(3) = 1.e0_dp/sqrt(3.e0_dp)
Pe2(1) = 1.e0_dp/sqrt(3.e0_dp)
Pe2(2) = 1.e0_dp/sqrt(3.e0_dp)
Pe2(3) = 1.e0_dp/sqrt(3.e0_dp)
Pe3(1) = 1.e0_dp/sqrt(3.e0_dp)
Pe3(2) = 1.e0_dp/sqrt(3.e0_dp)
Pe3(3) = 1.e0_dp/sqrt(3.e0_dp)

!endif

endif

endif

nQD    = nQDA+nQDB
ndim   = nhomoA+nhomoB+nhetero
nsys   = nQDA+nQDB+nhomoA+nhomoB+nhetero     
totsys = nsys+nhomoA+nhomoB+nhetero

allocate(aR(totsys),linker(nsys),epsin(totsys),epsR(totsys),V0e(totsys),V0h(totsys),source=0.d0)
allocate(QDcoor(totsys,3))
allocate(Dcenter(totsys,3))

if ( idlink .eq. 20 ) then
linker = 0.2e-9_dp
elseif ( idlink .eq. 55 ) then
linker = 0.55e-9_dp
endif
 
call random_seed(size = n)
allocate(seed(n))
call random_seed(get=seed)

if ( get_sp .eq. 'n' ) then
!   do n=1,nQDA
!   aR(n) = r8_NORMAL_AB(aA,dispQD*1e-9_dp,seed(1))
!   enddo
!   do n=nQDA+1,nQDA+nQDB
!   aR(n) = r8_NORMAL_AB(aB,dispQD*1e-9_dp,seed(2))
!   enddo
!   do n=nQDA+nQDB+1,nQDA+nQDB+nhomoA
!   aR(n) = r8_NORMAL_AB(aA,dispQD*1e-9_dp,seed(1))
!   aR(n+ndim) = r8_NORMAL_AB(aA,dispQD*1e-9_dp,seed(2))
!   enddo
!   do n=nQDA+nQDB+nhomoA+1,nQDA+nQDB+nhomoA+nhomoB
!   aR(n) = r8_NORMAL_AB(aB,dispQD*1e-9_dp,seed(1))
!   aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*1e-9_dp,seed(2))
!   enddo
!   do n=nQDA+nQDB+nhomoA+nhomoB+1,nQDA+nQDB+nhomoA+nhomoB+nhetero
!   aR(n) = r8_NORMAL_AB(aA,dispQD*1e-9_dp,seed(1))
!   aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*1e-9_dp,seed(2))
!   enddo
   do n=1,nQDA
   !aR(n) = r8_NORMAL_AB(aA,dispQD*aA/sqrt(2._dp*log(2._dp)),seed(1))
   aR(n) = r8_NORMAL_AB(aA,dispQD*aA,seed(1))
   enddo
   do n=nQDA+1,nQDA+nQDB
   !aR(n) = r8_NORMAL_AB(aB,dispQD*aB/sqrt(2._dp*log(2._dp)),seed(2))
   aR(n) = r8_NORMAL_AB(aB,dispQD*aB,seed(2))
   enddo
   do n=nQDA+nQDB+1,nQDA+nQDB+nhomoA
   !aR(n) = r8_NORMAL_AB(aA,dispQD*aA/sqrt(2._dp*log(2._dp)),seed(1))
   !aR(n+ndim) = r8_NORMAL_AB(aA,dispQD*aA/sqrt(2._dp*log(2._dp)),seed(2))
   aR(n) = r8_NORMAL_AB(aA,dispQD*aA,seed(1))
   aR(n+ndim) = r8_NORMAL_AB(aA,dispQD*aA,seed(2))
   enddo
   do n=nQDA+nQDB+nhomoA+1,nQDA+nQDB+nhomoA+nhomoB
   !aR(n) = r8_NORMAL_AB(aB,dispQD*aA/sqrt(2._dp*log(2._dp)),seed(1))
   !aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*aB/sqrt(2._dp*log(2._dp)),seed(2))
   aR(n) = r8_NORMAL_AB(aB,dispQD*aB,seed(1))
   aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*aB,seed(2))
   enddo
   do n=nQDA+nQDB+nhomoA+nhomoB+1,nQDA+nQDB+nhomoA+nhomoB+nhetero
   !aR(n) = r8_NORMAL_AB(aA,dispQD*aA/sqrt(2._dp*log(2._dp)),seed(1))
   !aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*aB/sqrt(2._dp*log(2._dp)),seed(2))
   aR(n) = r8_NORMAL_AB(aA,dispQD*aA,seed(1))
   aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*aB,seed(2))
   enddo
elseif ( get_sp .eq. 'y' ) then
   call system("mv Etransitions-he_0.dat tmp.dat ")
   open(newunit=tmp,file="tmp.dat")
   do n=1,nQDA+nQDB
   read(tmp,*) aR(n)
   aR(n) = aR(n)*1.e-9_dp
   print*, aR(n)
   enddo
   do n=nQDA+nQDB+1,nQDA+nQDB+nhomoA+nhomoB+nhetero
   read(tmp,*) aR(n), aR(n+ndim)
   aR(n) = aR(n)*1.e-9_dp
   aR(n+ndim) = aR(n+ndim)*1.e-9_dp
   print*, aR(n), aR(n+ndim)
   enddo
endif

do n=1,totsys
epsin(n) = 1.0e0_dp + (eps - 1.0e0_dp) / (1.0e0_dp + (0.75e-9_dp/(2.e0_dp*aR(n)))**1.2e0_dp)
epsR(n)  = 1.0e0_dp/((1.0e0_dp/epsin(n))-((1.0e0_dp/epsin(n))-(1.0e0_dp/(epsin(n)+3.5e0_dp)))*&
           (1.e0_dp-(exp(-(36e0_dp/35.e0_dp)*aR(n)/rhoe)+exp(-(36.e0_dp/35.e0_dp)*aR(n)/rhoh))/2.e0_dp))
V0e(n)   =-1.e0_dp*(-3.49e0_dp+2.47e0_dp*(2.e9_dp*aR(n))**(-1.32e0_dp))*elec
V0h(n)   =-1.e0_dp*(-5.23e0_dp-0.74e0_dp*(2.e9_dp*aR(n))**(-0.95e0_dp))*elec
enddo

if ( inbox .eq. 'y' ) then

open(newunit=sphere,file='sphere.xyz')
write(sphere,*) nint(totsys/2.e0_dp) 
write(sphere,*) 
do n = 1,nint(totsys/2.e0_dp)
Dcenter(n,:) = vectorin(1.e16_dp) 
Dcenter(n,:) = (Dcenter(n,:) + 1e16_dp**(1._dp/3._dp)) * 1.e-9_dp
write(sphere,*) 'H', Dcenter(n,1), Dcenter(n,2), Dcenter(n,3) 
enddo
!Dcenter = (Dcenter-(1.e16_dp)**(1._dp/3._dp) * 1.e-9_dp 

endif

end subroutine getVariables
   
end module Variables_au
