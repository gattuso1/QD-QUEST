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
allocate(xpow_gaus(0:ntime+1))
allocate(pow_gaus_s(totsys,0:ntime+1))
allocate(xpow_gaus_s(totsys,0:ntime+1))
allocate(xpulse(0:ntime+1))
allocate(pulses(0:ntime+1))
allocate(wft(0:ntime+1))
allocate(wft_s(totsys,0:ntime+1))
allocate(wft_pol(npol,0:ntime+1))
allocate(wftp(0:ntime+1))
allocate(wftf_s(totsys,0:ntime+1))
allocate(wftf(0:ntime+1))
allocate(wftf_pol(npol,0:ntime+1))
allocate(pow_pol(npol,0:ntime+1))
allocate(pow_s(totsys,0:ntime+1))
allocate(pow_pol_diff(0:ntime+1))
allocate(pow_pol_gaus(npol,0:ntime+1))
allocate(l1(npol))
allocate(l2(npol))
allocate(l3(npol))

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

enddo

open(newunit=Pulse_f     ,file='Pulse.dat')             
open(newunit=Tmat_0_f    ,file='TransMat.dat')          
open(newunit=Tmat_ei_f   ,file='TransMat_ei.dat')       
open(newunit=Tmat_x_f    ,file='TransMat_ei_x.dat')     
open(newunit=Tmat_y_f    ,file='TransMat_ei_y.dat')     
open(newunit=Tmat_z_f    ,file='TransMat_ei_z.dat')     
open(newunit=H_0_f       ,file='Ham0.dat')          
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
open(newunit=DipSpec     ,file='DipSpec.dat')
open(newunit=P_Match_f   ,file='Phase_Match.dat')
open(newunit=DipSpec_NR_f,file='DipSpec-NR.dat')
open(newunit=DipSpec_R_f ,file='DipSpec-R.dat')

matrices = (/ Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f /)

write(Pulse_f,'("#  time                      pulse1                    pulse2                    pulse3")')

!$OMP PARALLEL DO private(lambda,work,Ham,Mat,TransHam_ei,Ham_ei,TransHam,Ham_dir,Ham_ex,Ham_l,xc0,xc,xc_ei)

do n=1,nsys

!write(*,*) omp_get_num_threads()

if ( ( n .le. nQDA+nQDB ) .and. ( model .eq. "SB" ) ) then

nstates=3
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
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

nstates=9
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_he
include 'Core.f90'

endif

enddo !end loop number of systems

call Convolution

!Write DipSpec sum of all systems + multiplied by a gaussian 
timestep  = timestep  * t_au
totaltime = totaltime * t_au

do t=0,ntime

time = t*timestep

pow_gaus(t)=exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow(t)/totsys)

if ( singleFT .eq. "y" ) then
do n=1,totsys
pow_gaus_s(n,t)=exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_s(n,t)/totsys)
enddo
endif

!if ( inbox .eq. 'y' ) then

!pow_pol_gaus(39,t)=&
!   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(39,t)/totsys)
!pow_pol_gaus(41,t)=&
!   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(41,t)/totsys)
!do pol=1,npol
!pow_pol_gaus(pol,t)=&
!   exp(-1.d0*((time-totaltime/2.d0)*timestep)**2.d0/(2.d0*(totaltime*timestep/15.d0)**2.d0))*(pow_pol(pol,t)/totsys)
!enddo

!write(DipSpec_NR_f,*) time, dreal(pow_pol_gaus(39,t))
!write(DipSpec_R_f,*) time, dreal(pow_pol_gaus(41,t))

!endif

write(DipSpec,*) time, pow(t), pow_gaus(t), pulses(t)

enddo 

if ( inbox .eq. 'y' ) then

do pol=1,npol
integPol = dcmplx(0.d0,0.d0)
do t=0,ntime
time = t*timestep
integPol = integPol + abs(dcmplx(timestep,0.e0_dp)*(pow_pol(pol,t) + pow_pol(pol,t+1))/2.e0_dp)
enddo

!write(P_Match_f,'(i4,2x,3f6.2,es18.7e3)') pol, l1(pol), l2(pol), l3(pol), integPol
write(P_Match_f,*) pol, l1(pol), l2(pol), l3(pol), dreal(integPol)
!if ( pol .eq. 39 ) then
!write(6,*) pol, dreal(integPol)
!endif
enddo

!integPol_diff = dcmplx(0.d0,0.d0)

!do t=0,ntime
!time = t*timestep
!integPol_diff = integPol_diff + abs(dcmplx(timestep,0.e0_dp)*(pow_pol_diff(t) + pow_pol_diff(t+1))/2.e0_dp)
!enddo
!
!write(6,*) pol, dreal(integPol_diff)

endif

if ( doFT .eq. 'y' ) then

wstep = (w2-w1)/ntime
xpow_gaus_s  = dcmplx(pow_gaus_s,0.d0)
xpow_gaus  = dcmplx(pow_gaus,0.d0)
xpulse  = dcmplx(pulses,0.d0)

w = w1 

do t=0,ntime
do t2=0,ntime
!if ( inbox .eq. 'y' ) then
!wft_pol(39,t)  = wft_pol(39,t)  + exp(-2.d0*pi*im*w*t2*timestep) * pow_pol_gaus(39,t2) 
!wft_pol(41,t)  = wft_pol(41,t)  + exp(-2.d0*pi*im*w*t2*timestep) * pow_pol_gaus(41,t2) 
!do pol=30,40
!wft_pol(pol,t)  = wft_pol(pol,t)  + exp(-2.d0*pi*im*w*t2*timestep) * pow_pol_gaus(pol,t2) 
!enddo
!endif
if ( singleFT .eq. "y" ) then
do n=1,totsys
wft_s(n,t)  = wft_s(n,t)  + exp(-2.d0*pi*im*w*t2*timestep) * xpow_gaus_s(n,t2)
enddo
endif

wft(t)  = wft(t)  + exp(-2.d0*pi*im*w*t2*timestep) * xpow_gaus(t2) 
wftp(t) = wftp(t) + exp(-1.d0*im*w*t2*timestep) * xpulse(t2) 
enddo
w = w + wstep
enddo

w = w1

do t=0,ntime

wftf(t)= -2.e0_dp * dimag(sqrt(dreal(wft(t))**2+dimag(wft(t))**2) * dconjg(wftp(t)))
write(TransAbs,*)    w*h/elec, dreal(wftf(t))

!if ( inbox .eq. 'y' ) then
!wftf_pol(39,t) = -2.d0 * dimag(sqrt(dreal(wft_pol(39,t))**2+dimag(wft_pol(39,t))**2) * dconjg(wftp(t)))
!wftf_pol(41,t) = -2.d0 * dimag(sqrt(dreal(wft_pol(41,t))**2+dimag(wft_pol(41,t))**2) * dconjg(wftp(t)))
!write(TransAbs_NR,*) w*h/elec, (abs(wft_pol(pol,t)), pol=30,40) !dreal(wftf_pol(39,t))
!write(TransAbs_R,*)  w*h/elec, dreal(wft_pol(41,t)), dimag(wft_pol(41,t)) !dreal(wftf_pol(41,t))
!endif
w = w + wstep
enddo

if ( singleFT .eq. "y" ) then 

do n=1,totsys

write(cov2,'(a9,i0,a4)') 'TransAbs-', n, '.dat'
open(TransAbs_s,file=cov2)

w = w1

do t=0,ntime
wftf_s(n,t)= -2.e0_dp * dimag(sqrt(dreal(wft_s(n,t))**2+dimag(wft_s(n,t))**2) * dconjg(wftp(t)))
write(TransAbs_s,*)    w*h/elec, dreal(wftf_s(n,t))
w = w + wstep
enddo

close(TransAbs_s)

enddo

endif

endif

call system('find . -size 0 -delete')

!!$OMP END DO

!$OMP END PARALLEL DO

end
