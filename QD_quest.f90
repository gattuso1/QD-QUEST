include 'specfun.f90'

program ModelQD_one

!use omp_lib
use Constants_au
use Variables_au
use Integrals
use Vectors
use Normal
use Make_Ham

implicit none

real(dp), external:: s13adf, ei, eone, nag_bessel_j0
integer :: je,jh,nsteps,r,ifail, r1, r2
integer,dimension(10) :: matrices 
integer,dimension(4)  :: matrices_avg 
real(dp) :: Ef,delta, mu, A, I1eh1, I1eh2, I2eh1, I2eh2, I3eh1, I3eh2
real(dp),allocatable :: Ae(:), Ah1(:), Ah2(:), Be(:), Bh1(:), Bh2(:)
real(dp),allocatable :: kine(:), kinh1(:), kinh2(:)
real(dp),allocatable :: koute(:), kouth1(:), kouth2(:),diffe(:), diffh(:), E(:)

call getVariables

delta=  0.00001e-18_dp
Ef=     1.28174e-18_dp
nsteps= int(Ef/delta)
ifail=  1

!CALL OMP_SET_NUM_THREADS(omp_get_max_threads())
!CALL OMP_SET_NUM_THREADS(4)

!   write(*,*) omp_get_num_procs()
!   write(*,*) omp_get_max_threads()
!   write(*,*) omp_get_num_threads()
!   write(*,*) omp_get_max_threads()
!   write(*,*) omp_get_num_threads()

allocate(Eg(1000))
allocate(TDM(1000))
allocate(E(nsteps))
allocate(diffe(nsteps))
allocate(diffh(nsteps))
allocate(minEe(nsteps,totsys))
allocate(minEh(nsteps,totsys))
allocate(Eeh1(totsys)) 
allocate(Eeh2(totsys)) 
allocate(Ae(totsys)) 
allocate(Ah1(totsys)) 
allocate(Ah2(totsys)) 
allocate(Be(totsys)) 
allocate(Bh1(totsys)) 
allocate(Bh2(totsys)) 
allocate(kine(totsys)) 
allocate(kinh1(totsys)) 
allocate(kinh2(totsys)) 
allocate(koute(totsys)) 
allocate(kouth1(totsys)) 
allocate(kouth2(totsys))
allocate(TransDip_Ana_h1e(totsys))
allocate(TransDip_Ana_h2e(totsys))
allocate(TransDip_Ana_h1h2(totsys))
allocate(pow(0:ntime+1))
allocate(pow_gaus(0:ntime+1))
allocate(pow_gaus_s(totsys,0:ntime+1))
allocate(pulses(0:ntime+1))
allocate(pulses_FFT(0:ntime+1))
allocate(pow_pol(npol,0:ntime+1))
allocate(pow_s(totsys,0:ntime+1))
allocate(pow_pol_diff(0:ntime+1))
allocate(pow_pol_conv(0:ntime+1))
allocate(pow_pol_gaus(npol,0:ntime+1))
allocate(l1(npol),source=0._dp)
allocate(l2(npol),source=0._dp)
allocate(l3(npol),source=0._dp)
allocate(scos(npol),source=0.e0_dp)
allocate(ssin(npol),source=0.e0_dp)
allocate(TransHam_avg(0:2,0:2))
allocate(TransHam_avg_l(0:2,0:2,3))
allocate(Ham0_avg(0:2,0:2))
allocate(eTDM(9,9,3),source=0._dp)
allocate(rotmat(3,3))

if ( inbox .eq. 'y' ) then
call get_phases
endif

k=1

!Computation of energies
n=0

!!$OMP PARALLEL 
!!print *, "Hello"
!!$OMP END PARALLEL

do n = 1,totsys

i=0
r=0

je=1
jh=1

diffe = 0.d0
diffh = 0.d0

do i=1,nsteps
E(i)=delta*i
diffe(i)=abs(sqrt(2.0d0*me*E(i))/hbar * aR(n) * 1.0d0/tan(sqrt(2.e0_dp*me*E(i))/hbar * aR(n)) - 1.0d0 + (me/m0) + (me*aR(n))/(hbar)&
           * sqrt((2.0d0/m0)*(V0e(n)-E(i))))
if ((diffe(0) .eq. 0.000) .and. (diffh(0) .eq. 0.000)) then 
        diffe(0)=diffe(i)
        diffh(0)=diffe(i)
endif

if (diffe(i) .le. diffe(i-1)) then
        minEe(je,n) = E(i)
elseif ( (diffe(i) .ge. diffe(i-1)) .and. (E(i-1) .eq. minEe(je,n)) ) then
        je=je+1
endif

diffh(i)=abs(sqrt(2.d0*mh*E(i))/hbar * aR(n) * 1.d0/tan(sqrt(2.e0_dp*mh*E(i))/hbar * aR(n)) - 1.d0 + (mh/m0) + (mh*aR(n))/(hbar)&
           * sqrt((2.d0/m0)*(V0h(n)-E(i))))
if (diffh(i) .le. diffh(i-1)) then
        minEh(jh,n) = E(i)
elseif ( (diffh(i) .ge. diffh(i-1)) .and. (E(i-1) .eq. minEh(jh,n)) ) then
        jh=jh+1
endif
enddo

!wave vectors in and out
kine(n)=sqrt(2.d0*me*minEe(1,n))/hbar
koute(n)=sqrt(2.d0*m0*(V0e(n)-minEe(1,n)))/hbar

kinh1(n)=sqrt(2.d0*mh*minEh(1,n))/hbar
kouth1(n)=sqrt(2.d0*m0*(V0h(n)-minEh(1,n)))/hbar

kinh2(n)=sqrt(2.d0*mh*minEh(2,n))/hbar
kouth2(n)=sqrt(2.d0*m0*(V0h(n)-minEh(2,n)))/hbar

!normalization factors
Ae(n)=1.d0/sqrt(aR(n)/2.d0-sin(2.d0*kine(n)*aR(n))/(4.d0*kine(n))+sin(kine(n)*aR(n))**2/(2.d0*koute(n)))
Be(n)=Ae(n)*sin(kine(n)*aR(n))*exp(koute(n)*aR(n)) 

Ah1(n)=1.d0/sqrt(aR(n)/2.d0-sin(2.d0*kinh1(n)*aR(n))/(4.d0*kinh1(n))+sin(kinh1(n)*aR(n))**2/(2.d0*kouth1(n)))
Bh1(n)=Ah1(n)*sin(kinh1(n)*aR(n))*exp(kouth1(n)*aR(n)) 

Ah2(n)=1.d0/sqrt(aR(n)/2.d0-sin(2.d0*kinh2(n)*aR(n))/(4.d0*kinh2(n))+sin(kinh2(n)*aR(n))**2/(2.d0*kouth2(n)))
Bh2(n)=Ah2(n)*sin(kinh2(n)*aR(n))*exp(kouth2(n)*aR(n))

alphae  = kine(n) *aR(n)
alphah1 = kinh1(n)*aR(n)
alphah2 = kinh2(n)*aR(n)

betae  = koute(n) *aR(n)
betah1 = kouth1(n)*aR(n)
betah2 = kouth2(n)*aR(n)

!Coulomb Ferreyra
I1eh1=(Ae(n)**2*Ah1(n)**2*aR(n)/4)*(1-sin(2*alphae)/(2*alphae)-s13adf(2*alphae,ifail)/(2*alphae)+&
         (s13adf(2*alphae-2*alphah1,ifail)+s13adf(2*alphae+2*alphah1,ifail))/(4*alphae)) + &
         (Ae(n)**2*Ah1(n)**2*aR(n)/4)*(1-sin(2*alphah1)/(2*alphah1)-s13adf(2*alphah1,ifail)/(2*alphah1)+&
         (s13adf(2*alphah1-2*alphae,ifail)+s13adf(2*alphah1+2*alphae,ifail))/(4*alphah1))
I2eh1=(Ae(n)**2*Bh1(n)**2*aR(n)/2)*(1-sin(2*alphae)/(2*alphae))*eone(2*betah1)+&
         (Ah1(n)**2*Be(n)**2*aR(n)/2)*(1-sin(2*alphah1)/(2*alphah1))*eone(2*betae)
I3eh1=(Be(n)**2*Bh1(n)**2*aR(n)/(2*betae)*(exp(-2*betae)*eone(2*betah1)-eone(2*betae+2*betah1)))+&
         (Be(n)**2*Bh1(n)**2*aR(n)/(2*betah1)*(exp(-2*betah1)*eone(2*betae)-eone(2*betah1+2*betae)))

I1eh2=(Ae(n)**2*Ah2(n)**2*aR(n)/4)*(1-sin(2*alphae)/(2*alphae)-s13adf(2*alphae,ifail)/(2*alphae)+&
         (s13adf(2*alphae-2*alphah2,ifail)+s13adf(2*alphae+2*alphah2,ifail))/(4*alphae)) + &
         (Ae(n)**2*Ah2(n)**2*aR(n)/4)*(1-sin(2*alphah2)/(2*alphah2)-s13adf(2*alphah2,ifail)/(2*alphah1)+&
         (s13adf(2*alphah2-2*alphae,ifail)+s13adf(2*alphah2+2*alphae,ifail))/(4*alphah2))
I2eh2=(Ae(n)**2*Bh2(n)**2*aR(n)/2)*(1-sin(2*alphae)/(2*alphae))*eone(2*betah2)+&
         (Ah2(n)**2*Be(n)**2*aR(n)/2)*(1-sin(2*alphah2)/(2*alphah2))*eone(2*betae)
I3eh2=(Be(n)**2*Bh2(n)**2*aR(n)/(2*betae)*(exp(-2*betae)*eone(2*betah2)-eone(2*betae+2*betah2)))+&
         (Be(n)**2*Bh2(n)**2*aR(n)/(2*betah2)*(exp(-2*betah2)*eone(2*betae)-eone(2*betah2+2*betae)))

Eeh1(n)=(elec/(4.d0*PI*eps0*epsin(n)))*(I1eh1+I2eh1+I3eh1)
Eeh2(n)=(elec/(4.d0*PI*eps0*epsin(n)))*(I1eh2+I2eh2+I3eh2)

TransDip_Ana_h1e(n) = abs(TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
TransDip_Ana_h1h2(n) = abs(TransDip_Ana(Ah1(n),Ah2(n),Bh1(n),Bh2(n),kinh1(n),kinh2(n),kouth1(n),kouth2(n),aR(n)))
TransDip_Ana_h2e(n) = abs(TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))

!write(6,*) aR(n), minEe(1,n)/elec, minEh(1,n)/elec, minEh(2,n)/elec

enddo

open(newunit=Pulse_f     ,file='Pulse.dat')             
open(newunit=Tmat_0_f    ,file='TransMat.dat')          
open(newunit=Tmat_ei_f   ,file='TransMat_ei.dat')       
open(newunit=H_0_f       ,file='Ham0.dat')          
open(newunit=Tmat_y_f    ,file='TransMat_ei_y.dat')     
open(newunit=Tmat_z_f    ,file='TransMat_ei_z.dat')     
open(newunit=Tmat_x_f    ,file='TransMat_ei_x.dat')     
open(newunit=H_dir_f     ,file='Ham_dir.dat')           
open(newunit=H_ex_f      ,file='Ham_ex.dat')            
open(newunit=H_JK_f      ,file='Ham_JK.dat')            
open(newunit=H_ei_f      ,file='Ham_ei.dat')            
open(newunit=Etr_0_f     ,file='Etransitions-he_0.dat') 
open(newunit=Etr_ei_f    ,file='Etransitions-he_ei.dat')
open(newunit=TDip_ei_f   ,file='TransDip_ei.dat')       
open(newunit=Abs_imp_f   ,file='Absorption-imp.dat')    
open(newunit=Liou_f      ,file='Liou.dat')
open(newunit=TransAbs    ,file='TransAbs.dat')
open(newunit=TransAbs_NR ,file='TransAbs_NR.dat')
open(newunit=TransAbs_R  ,file='TransAbs_R.dat')
open(newunit=TransAbs_P  ,file='TransAbs_P.dat')
open(newunit=DipSpec     ,file='DipSpec.dat')
open(newunit=P_Match_f   ,file='Phase_Match.dat')
open(newunit=DipSpec_NR_f,file='DipSpec-NR.dat')
open(newunit=DipSpec_R_f ,file='DipSpec-R.dat')
open(newunit=DipSpec_conv_f ,file='DipSpec-conv.dat')
open(newunit=TransAbs_7_f ,file='TransAbs_7.dat')
open(newunit=TransAbs_17_f,file='TransAbs_17.dat') 
open(newunit=TransAbs_33_f,file='TransAbs_33.dat')
open(newunit=TransAbs_39_f,file='TransAbs_39.dat')
open(newunit=TransAbs_41_f,file='TransAbs_41.dat')
open(newunit=TransAbs_43_f,file='TransAbs_43.dat')
open(newunit=TransAbs_44_f,file='TransAbs_44.dat')
open(newunit=DipSpec_7_f ,file='DipSpec_7.dat')
open(newunit=DipSpec_17_f,file='DipSpec_17.dat') 
open(newunit=DipSpec_33_f,file='DipSpec_33.dat')
open(newunit=DipSpec_39_f,file='DipSpec_39.dat')
open(newunit=DipSpec_41_f,file='DipSpec_41.dat')
open(newunit=DipSpec_43_f,file='DipSpec_43.dat')
open(newunit=DipSpec_44_f,file='DipSpec_44.dat')
open(newunit=scos_ssin   ,file='ScosSsin.dat')
open(newunit=label_0   ,file='label_0.dat')
open(60,file='CoheTEST.dat')
allocate(mul(24))
mul = (/ "T", "T", "S", "T", "T", "T", "S", "T", "T", "S", "T", "T", "T", "T", "S", "T", "T", "T", "S", "T", "T", "S", "T", "T" /)

matrices = (/ Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f /)

write(Pulse_f,'("#  time                      pulse1                    pulse2                    pulse3")')

!$OMP PARALLEL DO private(lambda,work,Ham,Mat,TransHam_ei,Ham_ei,TransHam,Ham_dir,Ham_ex,Ham_l,xc0,xc,xc_ei)

call make_eTDM

do n=1,nsys

!write(*,*) omp_get_num_threads()

if ( ( n .le. nQDA+nQDB ) .and. ( model .eq. "SB" ) ) then

nstates=3
!print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_singl
include 'Core.f90'

elseif ( ( n .le. nQDA+nQDB ) .and. ( model .eq. "FS" ) ) then

nstates=55
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_FS_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "FO" ) ) then

nstates=5

print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_he_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "FS" ) ) then

!nstates=49
nstates=109
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_FS_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "SB" ) ) then

call cpu_time(start)

nstates=9
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_he
include 'Core.f90'

call cpu_time(finish)

write(6,*) "Time", finish-start

endif

enddo !end loop number of systems

deallocate(E,diffe,diffh,minEe,minEh,Eeh1,Eeh2,Ae,Ah1,Ah2,Be,Bh1,Bh2,kine,kinh1,kinh2,koute,kouth1,kouth2,TransDip_Ana_h1e,&
TransDip_Ana_h2e,TransDip_Ana_h1h2)
close(Tmat_0_f);close(Tmat_ei_f);close(H_0_f);close(Tmat_y_f);close(Tmat_z_f);close(Pulse_f);close(Tmat_x_f);close(H_dir_f)
close(H_ex_f);close(H_JK_f);close(H_ei_f);close(Etr_0_f);close(TDip_ei_f);close(Liou_f);close(Etr_ei_f)

!!!!Performs Liouville quantum dynamics on the average eigenstates Hamiltonian 
if ( Dyn_avg .eq. "y" ) then

open(newunit=Tmat_avg_f       ,file='TransMat_avg.dat')
open(newunit=Tmat_avgx_f      ,file='TransMat_avg_x.dat')
open(newunit=Tmat_avgy_f      ,file='TransMat_avg_y.dat')
open(newunit=Tmat_avgz_f      ,file='TransMat_avg_z.dat')
open(newunit=Etr_avg_f        ,file='Etransitions-he_avg.dat')
open(newunit=popc_ei_f       ,file='Pop_c_avg.dat')
open(newunit=Re_c_ei_f         ,file='Re_c_avg.dat')
open(newunit=Im_c_ei_f         ,file='Im_c_avg.dat')
open(newunit=TransAbs_avg     ,file='TransAbs_avg.dat')
open(newunit=DipSpec_avg      ,file='DipSpec_avg.dat')
open(newunit=P_Match_avg      ,file='Phase_Match_avg.dat')
open(newunit=TransAbs_avg_7_f ,file='TransAbs_avg_7.dat')
open(newunit=TransAbs_avg_17_f,file='TransAbs_avg_17.dat')
open(newunit=TransAbs_avg_33_f,file='TransAbs_avg_33.dat')
open(newunit=TransAbs_avg_39_f,file='TransAbs_avg_39.dat')
open(newunit=TransAbs_avg_41_f,file='TransAbs_avg_41.dat')
open(newunit=TransAbs_avg_43_f,file='TransAbs_avg_43.dat')
open(newunit=TransAbs_avg_44_f,file='TransAbs_avg_44.dat')
open(newunit=DipSpec_avg_7_f  ,file='DipSpec_avg_7.dat')
open(newunit=DipSpec_avg_17_f ,file='DipSpec_avg_17.dat')
open(newunit=DipSpec_avg_33_f ,file='DipSpec_avg_33.dat')
open(newunit=DipSpec_avg_39_f ,file='DipSpec_avg_39.dat')
open(newunit=DipSpec_avg_41_f ,file='DipSpec_avg_41.dat')
open(newunit=DipSpec_avg_43_f ,file='DipSpec_avg_43.dat')
open(newunit=DipSpec_avg_44_f ,file='DipSpec_avg_44.dat')

Dyn_avg_flag = 1
nstates=3
nstates2=nstates**2

allocate(haml(0:nstates-1,0:nstates-1))
allocate(k1_L(0:nstates2-1), &
         k2_L(0:nstates2-1),k3_L(0:nstates2-1),&
         k4_L(0:nstates2-1),k5_L(0:nstates2-1),&
         k6_L(0:nstates2-1),k7_L(0:nstates2-1),&
         k8_L(0:nstates2-1))
allocate(xc_L(0:nstates2-1,0:ntime+1))
allocate(k1(0:nstates-1), &
         k2(0:nstates-1),k3(0:nstates-1),&
         k4(0:nstates-1),k5(0:nstates-1),&
         k6(0:nstates-1),k7(0:nstates-1),&
         k8(0:nstates-1))
allocate(xc_ei(0:nstates-1,0:ntime+1))
allocate(merge_diag(0:nstates2-1,0:nstates2-1),merge_odiag(0:nstates2-1,0:nstates2-1), source = 0.e0_dp )
allocate(xliou(0:nstates-1,0:nstates-1,0:nstates-1,0:nstates-1),lfield(0:nstates2-1,0:nstates2-1))
allocate(irow(0:nstates2-1,2),icol(0:nstates2-1,2),source=0)
allocate(TransHam_ei(0:nstates-1,0:nstates-1),&
         TransHam_ei_l(0:nstates-1,0:nstates-1,3),&
         Mat(0:nstates-1,0:nstates-1),&
         Matx(0:nstates-1,0:nstates-1),&
         Maty(0:nstates-1,0:nstates-1),&
         Matz(0:nstates-1,0:nstates-1),&
         Ham_l(0:nstates-1,0:nstates-1),&
         Ham_ei(0:nstates-1,0:nstates-1),source=0.e0_dp)
xliou = dcmplx(0.e0_dp,0.e0_dp)
lfield = 0.e0_dp
Ham_l = 0.e0_dp

matrices_avg = (/ Tmat_avg_f, Tmat_avgx_f, Tmat_avgy_f, Tmat_avgz_f /)

print*, "Computing dynamics on averaged Hamiltonian"
Ham0_avg = Ham0_avg/nsys
TransHam_avg_l(:,:,1) = TransHam_avg_l(:,:,1)/nsys
TransHam_avg_l(:,:,2) = TransHam_avg_l(:,:,2)/nsys
TransHam_avg_l(:,:,3) = TransHam_avg_l(:,:,3)/nsys
TransHam_avg = TransHam_avg/nsys

Ham_ei = Ham0_avg
allocate(lambda(0:nstates-1),source = 0.e0_dp)
allocate(iwork2(3+5*nstates),source=0)
allocate(work1(6*nstates),work2(1+6*nstates+2*nstates*nstates),source=0.e0_dp)
lworku=1+6*nstates+2*nstates*nstates
liworku=3+5*nstates
ierr=0
call dsyevd('v','u',nstates, Ham_ei(0:nstates-1,0:nstates-1),nstates,lambda,work2,lworku,iwork2,liworku,ierr)
deallocate(work1)
deallocate(work2)
deallocate(iwork2)

call make_Ham_l

!write(6,*) (Ham_l(j,j)*Energ_au/elec, j=0,nstates-1)
!!!Make eigenstate TDM
if ( rdm_ori .eq. "n" ) then
Mat(:,:) = matmul(TransHam_avg(:,:),Ham_ei(:,:))
TransHam_ei(:,:) = matmul(transpose(Ham_ei(:,:)),Mat(:,:))
elseif ( rdm_ori .eq. "y" ) then
Matx(:,:) = matmul(TransHam_avg_l(:,:,1),Ham_ei(:,:))
Maty(:,:) = matmul(TransHam_avg_l(:,:,2),Ham_ei(:,:))
Matz(:,:) = matmul(TransHam_avg_l(:,:,3),Ham_ei(:,:))
TransHam_ei_l(:,:,1) = matmul(transpose(Ham_ei(:,:)),Matx(:,:))
TransHam_ei_l(:,:,2) = matmul(transpose(Ham_ei(:,:)),Maty(:,:))
TransHam_ei_l(:,:,3) = matmul(transpose(Ham_ei(:,:)),Matz(:,:))
TransHam_ei = sqrt(TransHam_ei_l(:,:,1)**2 + TransHam_ei_l(:,:,2)**2 + TransHam_ei_l(:,:,1)**2)
endif

write(6,*) TransHam_ei

!write(6,*) lambda*Energ_au/elec

include 'Core_avg.f90'

if ( nofiles .eq. 'y' ) then
close(popc_avg_f   ,status="delete")
close(Re_c_avg_f   ,status="delete")
close(Im_c_avg_f   ,status="delete")
endif

endif

if ( doAbs .eq. "y" ) then
call Convolution
endif


!Write DipSpec sum of all systems + multiplied by a gaussian 
timestep  = timestep  * t_au
totaltime = totaltime * t_au

do t=0,ntime

time = t*timestep

pow_gaus(t)=exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow(t)/nsys)

if ( doFT_s .eq. "y" ) then
pow_gaus_s(:,t)=exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*pow_s(:,t)
endif

if ( inbox .eq. 'y' ) then
pow_pol_gaus(39,t)=&
   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(39,t)/nsys)
pow_pol_gaus(41,t)=&
   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(41,t)/nsys)
!pow_pol_gaus(43,t)=&
!   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(43,t)/nsys)
!pow_pol_gaus(44,t)=&
!   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(44,t)/nsys)
!if ( nofiles .eq. 'n' ) then
!write(DipSpec_7_f,*)  time, dreal(pow_pol(7,t)), dimag(pow_pol(7,t))
!write(DipSpec_17_f,*) time, dreal(pow_pol(17,t)), dimag(pow_pol(17,t))
!write(DipSpec_33_f,*) time, dreal(pow_pol(33,t)), dimag(pow_pol(33,t))
!write(DipSpec_39_f,*) time, dreal(pow_pol(39,t)), dimag(pow_pol(39,t))
!write(DipSpec_41_f,*) time, dreal(pow_pol(41,t)), dimag(pow_pol(41,t))
!write(DipSpec_39_f,*) time, dreal(pow_pol_gaus(39,t)), dimag(pow_pol_gaus(39,t))
!write(DipSpec_41_f,*) time, dreal(pow_pol_gaus(41,t)), dimag(pow_pol_gaus(41,t))
!write(DipSpec_43_f,*) time, dreal(pow_pol_gaus(43,t)), dimag(pow_pol_gaus(43,t))
!write(DipSpec_44_f,*) time, dreal(pow_pol_gaus(44,t)), dimag(pow_pol_gaus(44,t))
!endif
!write(DipSpec_conv_f,*) time, dreal(pow_pol_conv(t)), dimag(pow_pol_conv(t)) 
endif

if ( nofiles .eq. 'n' ) then
write(DipSpec,*) time, pow(t), pow_gaus(t), pulses(t)
endif

enddo 

!write(6,*) l1(41), l2(41), l3(41) 
!write(6,*) l1(41)*Pe1(:), l2(41)*Pe2(:), l3(41)*Pe3(:)

!if ( inbox .eq. 'y' ) then
!do n=1,nsys
!
!!pow_41 = pow_41 + exp(-im*dcmplx(dot_product(l1(41)*k_1(:)+l2(41)*k_2(:)+l3(41)*k_3(:),Dcenter(n,:)),0.e0_dp))
!!pow_41 = pow_41 + exp(-im*dcmplx(dot_product(l1(41)*k_1(:)+l2(41)*k_2(:)+l3(41)*k_3(:),Dcenter(n,:)),0.e0_dp))
!pow_41 = pow_41 + exp(-im*dcmplx(dot_product(0.d0*Pe1(:)+0.d0*Pe2(:)+0.d0*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!write(6,*) pow_41, exp(-im*dcmplx(dot_product(l1(41)*k_1(:)+l2(41)*k_2(:)+l3(41)*k_3(:),Dcenter(n,:)),0.e0_dp))
!
!!scos = scos + cos((1.d0/545.d-9)*(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)))
!!
!!ssin = ssin + sin((1.d0/545.d-9)*(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)))
!
!!write(6,*) n, (1.d0/545.d-9)*(((l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!!((l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)) - &
!!(1.d0/545.d-9)*(((l1(39))*Pe1(1) + (l2(pol)-l2(39))*Pe2(1) + (l3(pol)-l3(39))*Pe3(1))*Dcenter(n,1) + &
!!((l1(39))*Pe1(2) + (l2(pol)-l2(39))*Pe2(2) + (l3(pol)-l3(39))*Pe3(2))*Dcenter(n,2) + &
!!((l1(39))*Pe1(3) + (l2(pol)-l2(39))*Pe2(3) + (l3(pol)-l3(39))*Pe3(3))*Dcenter(n,3))
!
!!write(6,*) ((l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1))*Dcenter(n,1), &
!!           ((l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1))*Dcenter(n,2), &
!!           ((l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1))*Dcenter(n,3)
!!write(6,*) ((l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1)), &
!
!!write(6,*) (l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1)
!
!!n, (1.d0/545.d-9)*(((l1(27)-l1(41))*Pe1(1) + (l2(27)-l2(41))*Pe2(1) + (l3(27)-l3(41))*Pe3(1))*Dcenter(n,1)+&
!!((l1(27)-l1(41))*Pe1(2) + (l2(27)-l2(41))*Pe2(2) + (l3(27)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(27)-l1(41))*Pe1(3) + (l2(27)-l2(41))*Pe2(3) + (l3(27)-l3(41))*Pe3(3))*Dcenter(n,3)), &
!!write(6,*) n, &
!!cos((1.d0/545.d-9)*(((l1(30)-l1(41))*Pe1(1) + (l2(30)-l2(41))*Pe2(1) + (l3(30)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!!((l1(30)-l1(41))*Pe1(2) + (l2(30)-l2(41))*Pe2(2) + (l3(30)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(30)-l1(41))*Pe1(3) + (l2(30)-l2(41))*Pe2(3) + (l3(30)-l3(41))*Pe3(3))*Dcenter(n,3))), &
!!sin((1.d0/545.d-9)*(((l1(30)-l1(41))*Pe1(1) + (l2(30)-l2(41))*Pe2(1) + (l3(30)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!!((l1(30)-l1(41))*Pe1(2) + (l2(30)-l2(41))*Pe2(2) + (l3(30)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!!((l1(30)-l1(41))*Pe1(3) + (l2(30)-l2(41))*Pe2(3) + (l3(30)-l3(41))*Pe3(3))*Dcenter(n,3)))
!
!enddo
!
!!write(6,*) scos, ssin
!
!!off diagonal convergence of PM signal
!do pol=1,npol
!if  ( pol .ne. 41)  then !.and. ( pol .ne. 1) .and. ( pol .ne. 3).and. ( pol .ne. 5).and. ( pol .ne. 24).and. ( pol .ne. 27) & 
!!.and. ( pol .ne. 32).and. ( pol .ne. 36).and. ( pol .ne. 37) .and. ( pol .ne. 39) .and. ( pol .ne. 19) .and. ( pol .ne. 15) ) then
!do n=1,nsys
!!pow_41_conv = pow_41_conv + exp(im * (1.d0/545.d-9) * & 
!
!scos(pol) = scos(pol) + cos(((1.d0/545.d-9)*&
!(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)))) 
!ssin(pol) = ssin(pol) + sin(((1.d0/545.d-9)*&
!(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3))))

!write(6,*) cos((1.d0/545.d-9)*&
!((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)), &
!           cos((1.d0/545.d-9)*&
! ((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
! ((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
! ((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)))



!write(6,*) pol, n, (1.d0/545.d-9)*(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)), &
!cos((1.d0/545.d-9)*&
!(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
! ((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
! ((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3))), &
!sin((1.d0/545.d-9)*(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3)))

!write(6,*) (1.d0/545.d-9)*(((l1(pol)-l1(41))*Pe1(1) + (l2(pol)-l2(41))*Pe2(1) + (l3(pol)-l3(41))*Pe3(1))*Dcenter(n,1) + &
!((l1(pol)-l1(41))*Pe1(2) + (l2(pol)-l2(41))*Pe2(2) + (l3(pol)-l3(41))*Pe3(2))*Dcenter(n,2) + &
!((l1(pol)-l1(41))*Pe1(3) + (l2(pol)-l2(41))*Pe2(3) + (l3(pol)-l3(41))*Pe3(3))*Dcenter(n,3))

!write(6,*) n, 1.d9*Dcenter(n,:)

!enddo
!
!write(scos_ssin,*) pol, scos(pol), ssin(pol)
!
!endif
!
!enddo

!pow_41_conv = pow_41_conv / npol

!write(6,*) dreal(pow_41), dimag(pow_41), dreal(pow_41_conv), dimag(pow_41_conv)

!endif

if ( inbox .eq. 'y' ) then

do pol=1,npol
integPol = dcmplx(0.d0,0.d0)
integPolconv = dcmplx(0.d0,0.d0)
do t=0,ntime
time = t*timestep
integPol = integPol + timestep*abs(pow_pol(pol,t) + pow_pol(pol,t+1))/2.e0_dp
!integPolconv = integPolconv + (dcmplx(timestep,0.e0_dp)*(pow_pol_conv(t) + pow_pol_conv(t+1)))/2.e0_dp
enddo
write(P_Match_f,*) pol, l1(pol), l2(pol), l3(pol), integPol!, abs(integPolconv)
enddo

!integPol_diff = dcmplx(0.d0,0.d0)

!do t=0,ntime
!time = t*timestep
!integPol_diff = integPol_diff + abs(dcmplx(timestep,0.e0_dp)*(pow_pol_diff(t) + pow_pol_diff(t+1))/2.e0_dp)
!enddo
!
!write(6,*) pol, dreal(integPol_diff)

endif

FTscale = h/(elec*(2.e0_dp**FTpow)*timestep)
!gives the scale of FT vectors
t=0
do while ( t*FTscale .le. 4.d0 )
t=t+1
enddo
nFT=t
FTscale = h/(elec*(2.e0_dp**FTpow)*timestep)

if ( doFT .eq. 'y' ) then

allocate(wft(0:nFT+1))
allocate(wft_s(totsys,0:nFT+1))
allocate(wft_pol(npol,0:nFT+1))
allocate(wftp(0:nFT+1))
allocate(wftf_s(totsys,0:nFT+1))
allocate(wftf(0:nFT+1))
allocate(xpow_gaus(0:nint(2.e0_dp**FTpow)))
allocate(xpulse(0:nint(2.e0_dp**FTpow)))

do t=0,ntime
xpow_gaus(t) = dcmplx(pow_gaus(t),0.d0)
xpulse(t)  = dcmplx(pulses_FFT(t),0.d0)
enddo

do t=ntime+1,nint(2.d0**FTpow)
xpow_gaus(t) = dcmplx(0.d0,0.d0)
xpulse(t)  = dcmplx(0.d0,0.d0)
enddo

call fft(xpulse)
call fft(xpow_gaus)

t=0
do while ( t*FTscale .le. 4.d0 ) 
wftf(t)= -2.e0_dp * dimag(sqrt(dreal(xpow_gaus(t))**2+dimag(xpow_gaus(t))**2) * dconjg(xpulse(t)))
!if ( nofiles .eq. 'n' ) then
write(TransAbs,*) t*FTscale, dreal(wftf(t)), dreal(xpow_gaus(t)), dimag(xpow_gaus(t)), abs(xpulse(t))
!write(6,*) t*FTscale, dreal(xpulse(t)), dimag(xpulse(t))
!endif
t = t + 1 
enddo

if ( doFT_s .eq. "y" ) then

do n=1,nsys

allocate(xpow_gaus_s(0:nint(2.e0_dp**FTpow)))

do t=0,ntime
xpow_gaus_s(t)  = dcmplx(pow_gaus_s(n,t),0.d0)
enddo
do t=ntime+1,nint(2.d0**FTpow)
xpow_gaus_s(t)  = dcmplx(0.d0,0.d0)
enddo

call fft(xpow_gaus_s)

t=0
do while ( t*FTscale .le. 4.d0 )
wftf_s(n,t)= -2.e0_dp * dimag(sqrt(dreal(xpow_gaus_s(t))**2+dimag(xpow_gaus_s(t))**2) * dconjg(xpulse(t)))
t = t + 1
enddo

deallocate(xpow_gaus_s)

enddo

if ( singleFT .eq. 'y' ) then
do n=1,nsys
write(cov2,'(a9,i0,a4)') 'TransAbs-', n, '.dat'
open(TransAbs_s,file=cov2)
t=0
do while ( t*FTscale .le. 4.d0 )
write(TransAbs_s,*) t*FTscale,  dreal(wftf_s(n,t))
t = t + 1
enddo
close(TransAbs_s)
enddo
endif

endif

if ( inbox .eq. 'y' ) then

allocate(xpow_pol(44,0:nint(2.e0_dp**FTpow)))
allocate(wftf_pol(44,0:nFT+1))

do t=0,ntime
xpow_pol(39,t)  = pow_pol_gaus(39,t)
xpow_pol(41,t)  = pow_pol_gaus(41,t)
!xpow_pol(43,t)  = pow_pol_gaus(43,t)
!xpow_pol(44,t)  = pow_pol_gaus(44,t)
enddo
do t=ntime+1,nint(2.d0**FTpow)
xpow_pol(39,t)  = dcmplx(0.d0,0.d0)
xpow_pol(41,t)  = dcmplx(0.d0,0.d0)
!xpow_pol(43,t)  = dcmplx(0.d0,0.d0)
!xpow_pol(44,t)  = dcmplx(0.d0,0.d0)
enddo

!do pol=1,npol
call fft(xpow_pol(39,:))
call fft(xpow_pol(41,:))
!call fft(xpow_pol(43,:))
!call fft(xpow_pol(44,:))
!enddo

t=0
do while ( t*FTscale .le. 4.d0 )
!wftf_pol(:,t) = -2.d0 * dimag(sqrt(dreal(xpow_pol(:,t))**2+dimag(xpow_pol(:,t))**2) * dconjg(xpulse(t)))
!write(TransAbs_7_f,*)  t*FTscale, dreal(xpow_pol(7,t)) , dimag(xpow_pol(7,t))
!write(TransAbs_17_f,*) t*FTscale, dreal(xpow_pol(17,t)), dimag(xpow_pol(17,t))
!write(TransAbs_33_f,*) t*FTscale, dreal(xpow_pol(33,t)), dimag(xpow_pol(33,t))
write(TransAbs_39_f,*) t*FTscale, dreal(xpow_pol(39,t)), dimag(xpow_pol(39,t)), abs(xpow_pol(39,t))
write(TransAbs_41_f,*) t*FTscale, dreal(xpow_pol(41,t)), dimag(xpow_pol(41,t)), abs(xpow_pol(41,t))
!write(TransAbs_43_f,*) t*FTscale, dreal(xpow_pol(43,t)), dimag(xpow_pol(43,t))
!write(TransAbs_44_f,*) t*FTscale, dreal(xpow_pol(44,t)), dimag(xpow_pol(44,t))
t = t + 1
enddo

deallocate(xpow_pol,wftf_pol)

endif

deallocate(pow,pow_gaus,xpow_gaus,pow_gaus_s,xpulse,pulses,wft,wft_s,wftp)

endif

!
!if ( doCovar .eq. 'y' ) then
!
!allocate(Scov(0:nFT+1,0:nFT+1))
!
!open(61,file='Allcov.dat')
!
!t=0
!do while ( t*FTscale .le. 4.d0 )
!t2=0
!do while ( t2*FTscale .le. 4.d0 )
!
!do k=1,nsys
!Scov(t,t2) = Scov(t,t2) + sum(dreal(wftf_s(k,t))*dreal(wftf_s(:,t2)))
!enddo
!
!t2 = t2 + 5
!enddo
!t = t + 5
!enddo
!
!t=0
!do while ( t*FTscale .le. 4.d0 )
!t2=0
!do while ( t2*FTscale .le. 4.d0 )
!write(61,'(2f10.6,ES15.6E3)') t*h/(elec*5.24288d-12), t2*h/(elec*5.24288d-12), Scov(t,t2)
!t2 = t2 + 5
!enddo
!write(61,*)
!t = t + 5
!enddo
!
!endif

!call system('sleep 1')
!call system('find . -size 0 -delete')

!!$OMP END DO

!$OMP END PARALLEL DO

end
