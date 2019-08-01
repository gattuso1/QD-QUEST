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
integer :: Pulse_f,Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f,Etr_0_f,Etr_ei_f,Abs_imp_f
integer :: popc_0_f,popc_ei_f,norm_0_f,norm_ei_f,Re_c_ei_f,Im_c_ei_f,Re_c_0_f,Im_c_0_f,TDip_ei_f
character*64 :: form_mat,form_arr,form_abs
integer :: je,jh,k,nsteps,r,ifail, r1, r2
integer,dimension(10) :: matrices
real(dp) :: Ef,delta, mu, A
real(dp),allocatable :: Ae(:), Ah1(:), Ah2(:), Be(:), Bh1(:), Bh2(:)
real(dp),allocatable :: I1eh1(:), I1eh2(:), I2eh1(:), I2eh2(:), I3eh1(:), I3eh2(:), kine(:), kinh1(:), kinh2(:)
real(dp),allocatable :: koute(:), kouth1(:), kouth2(:),diffe(:), diffh(:), E(:)

call getVariables

delta=  0.00001d-18
Ef=     1.28174d-18
nsteps= int(Ef/delta)
ifail=  1

!CALL OMP_SET_NUM_THREADS(omp_get_max_threads())
!CALL OMP_SET_NUM_THREADS(4)

!   write(*,*) omp_get_num_procs()
!   write(*,*) omp_get_max_threads()
!   write(*,*) omp_get_num_threads()
!   write(*,*) omp_get_max_threads()
!   write(*,*) omp_get_num_threads()

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
diffe(i) = abs(sqrt(2.0d0*me*E(i))/hbar * aR(n) * 1.0d0/tan(sqrt(2*me*E(i))/hbar * aR(n)) - 1.0d0 + (me/m0) + (me*aR(n))/(hbar) &
           * sqrt((2.0d0/m0)*(V0e(n)-E(i))))
if ((diffe(0) .eq. 0.000) .and. (diffh(0) .eq. 0.00)) then 
        diffe(0)=diffe(i)
        diffh(0)=diffe(i)
endif

if (diffe(i) .le. diffe(i-1)) then
        minEe(je,n) = E(i)
elseif ( (diffe(i) .ge. diffe(i-1)) .and. (E(i-1) .eq. minEe(je,n)) ) then
        je=je+1
endif

diffh(i) = abs(sqrt(2.d0*mh*E(i))/hbar * aR(n) * 1.d0/tan(sqrt(2.d0*mh*E(i))/hbar * aR(n)) - 1.d0 + (mh/m0) + (mh*aR(n))/(hbar) &
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

TransDip_Ana_h1e(n) = abs(TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
TransDip_Ana_h1h2(n) = abs(TransDip_Ana(Ah1(n),Ah2(n),Bh1(n),Bh2(n),kinh1(n),kinh2(n),kouth1(n),kouth2(n),aR(n)))
TransDip_Ana_h2e(n) = abs(TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))

enddo

open(newunit=Pulse_f  ,file='Pulse.dat')             
open(newunit=Tmat_0_f ,file='TransMat.dat')          
open(newunit=Tmat_ei_f,file='TransMat_ei.dat')       
open(newunit=Tmat_x_f ,file='TransMat_ei_x.dat')     
open(newunit=Tmat_y_f ,file='TransMat_ei_y.dat')     
open(newunit=Tmat_z_f ,file='TransMat_ei_z.dat')     
open(newunit=H_0_f    ,file='Ham0.dat')          
open(newunit=H_dir_f  ,file='Ham_dir.dat')           
open(newunit=H_ex_f   ,file='Ham_ex.dat')            
open(newunit=H_JK_f   ,file='Ham_JK.dat')            
open(newunit=H_ei_f   ,file='Ham_ei.dat')            
open(newunit=Etr_0_f  ,file='Etransitions-he_0.dat') 
open(newunit=Etr_ei_f ,file='Etransitions-he_ei.dat')
open(newunit=TDip_ei_f,file='TransDip_ei.dat')       
open(newunit=Abs_imp_f,file='Absorption-imp.dat')    

matrices = (/ Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f /)

write(Pulse_f,'("#  time                      pulse1                    pulse2                    pulse3")')

!$OMP PARALLEL DO private(lambda,work,Ham,Mat,TransHam_ei,Ham_ei,TransHam,Ham_dir,Ham_ex,Ham_l,xc0,xc,xc_ei)

do n=1,nsys

!write(*,*) omp_get_num_threads()

write(6,*) "Computing system number:    ", n

!!!Opens output files
write(popc,'(a5,i5.5,a4)') 'Popc-', n, '.dat'
write(norm,'(a5,i5.5,a4)') 'Norm-', n, '.dat'
write(norm_ei,'(a8,i5.5,a4)') 'Norm_ei-', n, '.dat'
write(popc_ei,'(a8,i5.5,a4)') 'Popc_ei-', n, '.dat'
write(Re_c,'(a5,i5.5,a4)') 'Re_c-', n, '.dat'
write(Im_c,'(a5,i5.5,a4)') 'Im_c-', n, '.dat'
write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
open(newunit=popc_0_f   ,file=popc)    !; popc_0_f = 44
open(newunit=popc_ei_f  ,file=popc_ei) !; popc_ei_f = 49
open(newunit=norm_0_f   ,file=norm)    !; norm_0_f = 46
open(newunit=norm_ei_f  ,file=norm_ei) !; norm_ei_f = 50
open(newunit=Re_c_ei_f  ,file=Re_c_ei) !; Re_c_ei_f = 52
open(newunit=Im_c_ei_f  ,file=Im_c_ei) !; Im_c_ei_f = 53
open(newunit=Re_c_0_f   ,file=Re_c)    !; Re_c_0_f = 54
open(newunit=Im_c_0_f   ,file=Im_c)    !; Im_c_0_f = 55

if ( n .le. nQDA+nQDB ) then

nstates=3

allocate(TransHam(0:nstates-1,0:nstates-1),TransHam_ei_l(0:nstates-1,0:nstates-1,3),TransHam_l(0:nstates-1,0:nstates-1,3),&
TransHam_d(0:nstates-1,0:nstates-1,3),TransHam_ei(0:nstates-1,0:nstates-1),Transvec(0:nstates-1),&
TransMat_ei(0:nstates-1,0:nstates-1),Mat(0:nstates-1,0:nstates-1),Matx(0:nstates-1,0:nstates-1),&
Maty(0:nstates-1,0:nstates-1),Matz(0:nstates-1,0:nstates-1),Ham(0:nstates-1,0:nstates-1),Ham_l(0:nstates-1,0:nstates-1),&
Ham_0(0:nstates-1),Ham_dir(0:nstates-1,0:nstates-1),Ham_ex(0:nstates-1,0:nstates-1),Ham_ei(0:nstates-1,0:nstates-1),source=0.d0)
allocate(xc(0:nstates-1,0:ntime),xc_ei(0:nstates-1,0:ntime+1),c0(0:nstates-1),k1(0:nstates-1),k2(0:nstates-1),k3(0:nstates-1),&
k4(0:nstates-1),k5(0:nstates-1),k6(0:nstates-1),k7(0:nstates-1),k8(0:nstates-1))

call make_Ham_singl

include 'Core.f90'

elseif (n .gt. nQDA+nQDB) then

nstates=9

allocate(TransHam(0:nstates-1,0:nstates-1),TransHam_ei_l(0:nstates-1,0:nstates-1,3),TransHam_l(0:nstates-1,0:nstates-1,3),&
TransHam_d(0:nstates-1,0:nstates-1,3),TransHam_ei(0:nstates-1,0:nstates-1),Transvec(0:nstates-1),&
TransMat_ei(0:nstates-1,0:nstates-1),Mat(0:nstates-1,0:nstates-1),Matx(0:nstates-1,0:nstates-1),&
Maty(0:nstates-1,0:nstates-1),Matz(0:nstates-1,0:nstates-1),Ham(0:nstates-1,0:nstates-1),Ham_l(0:nstates-1,0:nstates-1),&
Ham_0(0:nstates-1),Ham_dir(0:nstates-1,0:nstates-1),Ham_ex(0:nstates-1,0:nstates-1),Ham_ei(0:nstates-1,0:nstates-1),source=0.d0)
allocate(xc(0:nstates-1,0:ntime),xc_ei(0:nstates-1,0:ntime+1),c0(0:nstates-1),k1(0:nstates-1),k2(0:nstates-1),k3(0:nstates-1),&
k4(0:nstates-1),k5(0:nstates-1),k6(0:nstates-1),k7(0:nstates-1),k8(0:nstates-1))

call make_Ham_he

include 'Core.f90'

endif

!
!!call cpu_time(finish)
!
!write(6,*) finish - start

enddo !end loop number of systems

!!$OMP END DO

!$OMP END PARALLEL DO

!call cpu_time(finish)

!write(6,*) finish - start

end
