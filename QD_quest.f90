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
integer :: popc_0_f,popc_ei_f,norm_0_f,norm_ei_f,Re_c_ei_f,Im_c_ei_f,Re_c_0_f,Im_c_0_f 
character*64 :: string1,string2
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

open(22,file='Pulse.dat')              ; Pulse_f = 22
open(32,file='TransMat.dat')           ; Tmat_0_f = 32
open(37,file='TransMat_ei.dat')        ; Tmat_ei_f = 37
open(38,file='TransMat_ei_x.dat')      ; Tmat_x_f = 38
open(39,file='TransMat_ei_y.dat')      ; Tmat_y_f = 39
open(40,file='TransMat_ei_z.dat')      ; Tmat_z_f = 40
open(33,file='Ham0.dat')               ; H_0_f = 33
open(34,file='Ham_dir.dat')            ; H_dir_f = 34
open(35,file='Ham_ex.dat')             ; H_ex_f = 35
open(36,file='Ham_JK.dat')             ; H_JK_f = 36
open(58,file='Ham_ei.dat')             ; H_ei_f = 58
open(47,file='Etransitions-he_0.dat')  ; Etr_0_f = 47
open(48,file='Etransitions-he_ei.dat') ; Etr_ei_f = 48
open(57,file='TransDip_ei.dat')        
open(60,file='Absorption-imp.dat')     ; Abs_imp_f = 60

matrices = (/ Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f /)

do i=1,size(matrices)
write(matrices(i),'("#     Number                  QDA                       QDB                    linker")')
enddo

!write(32,'("#     Number                  QDA                       QDB                    linker")')
!write(37,'("#     Number                  QDA                       QDB                    linker")')
!write(38,'("#     Number                  QDA                       QDB                    linker")')
!write(39,'("#     Number                  QDA                       QDB                    linker")')
!write(40,'("#     Number                  QDA                       QDB                    linker")')
!write(33,'("#     Number                  QDA                       QDB                    linker")')
write(22,'("#  time                      pulse1                    pulse2                    pulse3")')
!write(58,'("#     Number                  QDA                       QDB                    linker")')

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
open(44,file=popc)    ; popc_0_f = 44
open(49,file=popc_ei) ; popc_ei_f = 49
open(46,file=norm)    ; norm_0_f = 46
open(50,file=norm_ei) ; norm_ei_f = 50
open(52,file=Re_c_ei) ; Re_c_ei_f = 52
open(53,file=Im_c_ei) ; Im_c_ei_f = 53
open(54,file=Re_c)    ; Re_c_0_f = 54
open(55,file=Im_c)    ; Im_c_0_f = 55

if ( n .le. nQDA+nQDB ) then

nstates=3

!allocate(TransHam(0:nstates-1,0:nstates-1))
allocate(TransHam(0:nstates-1,0:nstates-1),TransHam_ei_l(0:nstates-1,0:nstates-1,3), source = 0.d0)
!allocate(TransHam_ei_l(0:nstates-1,0:nstates-1,3))
allocate(TransHam_l(0:nstates-1,0:nstates-1,3))
allocate(TransHam_d(0:nstates-1,0:nstates-1,3))
allocate(TransHam_ei(0:nstates-1,0:nstates-1))
allocate(Transvec(0:nstates-1))
allocate(TransMat_ei(0:nstates-1,0:nstates-1))
allocate(Mat(0:nstates-1,0:nstates-1))
allocate(Matx(0:nstates-1,0:nstates-1))
allocate(Maty(0:nstates-1,0:nstates-1))
allocate(Matz(0:nstates-1,0:nstates-1))
allocate(Ham(0:nstates-1,0:nstates-1),Ham_l(0:nstates-1,0:nstates-1), source = 0.d0)
allocate(Ham_0(0:nstates-1))
allocate(Ham_dir(0:nstates-1,0:nstates-1))
allocate(Ham_ex(0:nstates-1,0:nstates-1))
allocate(Ham_ei(0:nstates-1,0:nstates-1))
allocate(xc(0:nstates-1,0:ntime))
allocate(xc_ei(0:nstates-1,0:ntime+1))
allocate(c0(0:nstates-1))
allocate(k1(0:nstates-1))
allocate(k2(0:nstates-1))
allocate(k3(0:nstates-1))
allocate(k4(0:nstates-1))
allocate(k5(0:nstates-1))
allocate(k6(0:nstates-1))
allocate(k7(0:nstates-1))
allocate(k8(0:nstates-1))

Ham      = 0.0d0
TransHam = 0.0d0

call make_Ham_singl

write(32,*) n , aR(n) 
write(37,*) n , aR(n) 
write(38,*) n , aR(n) 
write(39,*) n , aR(n) 
write(40,*) n , aR(n) 
write(33,*) n , aR(n) 
write(34,*) n , aR(n) 
write(35,*) n , aR(n) 
write(36,*) n , aR(n) 

do i=0,nstates-1
write(33,'(3es16.6e2)') (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(34,'(3es16.6e2)') (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(35,'(3es16.6e2)') (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(36,'(3es16.6e2)') ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
write(32,'(3es16.6e2)') (TransHam(i,j), j=0,nstates-1)
enddo

write(32,*) 
write(33,*) 


write(string2,'("(",i0,"f16.10)")') nstates+1
write(6,*) string2

write(Etr_0_f,'(4f16.10)') aR(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(6,'(4g0)') aR(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_0_f,string2) aR(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)

if ( get_ei .eq. 'y' ) then
write(58,*) n , aR(n)
Ham_ei = Ham
allocate(lambda(0:nstates-1))
allocate(work(1))
call dsyev('V','U', nstates, Ham_ei(0:nstates-1,0:nstates-1), nstates, lambda, Work, -1, info)
lwork=nint(work(1))
deallocate (work)
allocate(work(0:lwork))
call dsyev('V', 'U', nstates, Ham_ei(0:nstates-1,0:nstates-1), nstates, lambda, Work, lwork, info)
deallocate (work)

!!!Make eigenstate TDM
if ( inbox .eq. "n" ) then
Mat(:,:) = matmul(TransHam(:,:),Ham_ei(:,:))
TransHam_ei(:,:) = matmul(transpose(Ham_ei(:,:)),Mat(:,:))
elseif ( inbox .eq. "y" ) then
Matx(:,:) = matmul(TransHam_l(:,:,1),Ham_ei(:,:))
Maty(:,:) = matmul(TransHam_l(:,:,2),Ham_ei(:,:))
Matz(:,:) = matmul(TransHam_l(:,:,3),Ham_ei(:,:))
TransHam_ei_l(:,:,1) = matmul(transpose(Ham_ei(:,:)),Matx(:,:))
TransHam_ei_l(:,:,2) = matmul(transpose(Ham_ei(:,:)),Maty(:,:))
TransHam_ei_l(:,:,3) = matmul(transpose(Ham_ei(:,:)),Matz(:,:))
endif

call make_Ham_l

do i=0,nstates-1
write(58,'(3es16.6e2)') (Ham_ei(i,j), j=0,nstates-1)
write(37,'(3es16.6e2)') (TransHam_ei(i,j), j=0,nstates-1)
if ( inbox .eq. "y" ) then
write(38,'(3es16.6e2)') (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(39,'(3es16.6e2)') (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(40,'(3es16.6e2)') (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif
enddo

write(58,*) 
write(48,'(4es16.6e2)') aR(n)*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)

endif

do i=1,nstates-1
write(60,*) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
enddo

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

!!!!!INITIAL POPULATIONS
c0(0) = 1.0d0
do i=1,nstates-1
c0(i) = 0.0d0
enddo

xc0 = dcmplx(c0,0.0d0)
xc = dcmplx(0.d0,0.0d0)
xc_ei = dcmplx(0.d0,0.0d0)
xc(:,0) = xc0(:)
xc_ei(:,0) = xc0(:)

call RK_0_ei

endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei)

elseif (n .gt. nQDA+nQDB) then

nstates=9

allocate(TransHam(0:nstates-1,0:nstates-1))
allocate(TransHam_ei_l(0:nstates-1,0:nstates-1,3))
allocate(TransHam_l(0:nstates-1,0:nstates-1,3))
allocate(TransHam_d(0:nstates-1,0:nstates-1,3))
allocate(TransHam_ei(0:nstates-1,0:nstates-1))
allocate(Transvec(0:nstates-1))
allocate(TransMat_ei(0:nstates-1,0:nstates-1))
allocate(Mat(0:nstates-1,0:nstates-1))
allocate(Matx(0:nstates-1,0:nstates-1))
allocate(Maty(0:nstates-1,0:nstates-1))
allocate(Matz(0:nstates-1,0:nstates-1))
allocate(Ham(0:nstates-1,0:nstates-1))
allocate(Ham_l(0:nstates-1,0:nstates-1))
allocate(Ham_0(0:nstates-1))
allocate(Ham_dir(0:nstates-1,0:nstates-1))
allocate(Ham_ex(0:nstates-1,0:nstates-1))
allocate(Ham_ei(0:nstates-1,0:nstates-1))
allocate(xc(0:nstates-1,0:ntime))
allocate(xc_ei(0:nstates-1,0:ntime+1))
allocate(c0(0:nstates-1))
allocate(k1(0:nstates-1))
allocate(k2(0:nstates-1))
allocate(k3(0:nstates-1))
allocate(k4(0:nstates-1))
allocate(k5(0:nstates-1))
allocate(k6(0:nstates-1))
allocate(k7(0:nstates-1))
allocate(k8(0:nstates-1))

Ham      = 0.0d0
TransHam = 0.0d0

call make_Ham_he

write(32,*) n , aR(n), aR(n+ndim), linker(n)
write(37,*) n , aR(n), aR(n+ndim), linker(n)
write(38,*) n , aR(n), aR(n+ndim), linker(n)
write(39,*) n , aR(n), aR(n+ndim), linker(n)
write(40,*) n , aR(n), aR(n+ndim), linker(n)
write(33,*) n , aR(n), aR(n+ndim), linker(n)
write(34,*) n , aR(n), aR(n+ndim), linker(n)
write(35,*) n , aR(n), aR(n+ndim), linker(n)
write(36,*) n , aR(n), aR(n+ndim), linker(n)

do i=0,nstates-1
write(33,'(9es14.6e2)') (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(34,'(9es14.6e2)') (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(35,'(9es14.6e2)') (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(36,'(9es14.6e2)') ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
!do j=0,nstates-1
!write(6,*) i , j , (-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec
!enddo
write(32,'(9es14.6e2)') (TransHam(i,j), j=0,nstates-1)
enddo
write(32,*)
write(33,*)

write(47,'(11f14.10)') aR(n)*1.d9, linker(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)

if ( get_ei .eq. 'y' ) then
write(58,*) n , aR(n), aR(n+ndim), linker(n)
Ham_ei = Ham
allocate(lambda(0:nstates-1))
allocate(work(1))
call dsyev('V','U', nstates, Ham_ei(0:nstates-1,0:nstates-1), nstates, lambda, Work, -1, info)
lwork=nint(work(1))
deallocate (work)
allocate(work(0:lwork))
call dsyev('V', 'U', nstates, Ham_ei(0:nstates-1,0:nstates-1), nstates, lambda, Work, lwork, info)
deallocate (work)

!!!Make eigenstate TDM
if ( inbox .eq. "n" ) then
Mat(:,:) = matmul(TransHam(:,:),Ham_ei(:,:))
TransHam_ei(:,:) = matmul(transpose(Ham_ei(:,:)),Mat(:,:))
elseif ( inbox .eq. "y" ) then
Matx(:,:) = matmul(TransHam_l(:,:,1),Ham_ei(:,:))
Maty(:,:) = matmul(TransHam_l(:,:,2),Ham_ei(:,:))
Matz(:,:) = matmul(TransHam_l(:,:,3),Ham_ei(:,:))
TransHam_ei_l(:,:,1) = matmul(transpose(Ham_ei(:,:)),Matx(:,:))
TransHam_ei_l(:,:,2) = matmul(transpose(Ham_ei(:,:)),Maty(:,:))
TransHam_ei_l(:,:,3) = matmul(transpose(Ham_ei(:,:)),Matz(:,:))
endif

call make_Ham_l

do i=0,nstates-1
write(58,'(10f12.6)') (Ham_ei(i,j), j=0,nstates-1)
write(37,'(9es14.6e2)') (TransHam_ei(i,j), j=0,nstates-1)

if ( inbox .eq. "y" ) then
write(38,'(9es14.6e2)') (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(39,'(9es14.6e2)') (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(40,'(9es14.6e2)') (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif

enddo

write(58,*) 

write(48,'(11f16.10)') aR(n)*1.d9, aR(n+ndim)*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)

endif

do i=1,nstates-1
write(60,*) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
enddo

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

!!!!!INITIAL POPULATIONS
c0(0) = 1.0d0
do i=1,nstates-1
c0(i) = 0.0d0
enddo

xc0 = dcmplx(c0,0.0d0)
xc(:,:) = dcmplx(0.d0,0.0d0)
xc_ei(:,:) = dcmplx(0.d0,0.0d0)
xc(:,0) = xc0(:)
xc_ei(:,0) = xc0(:)

call RK_0_ei

endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei)

endif

!call OK
!deallocate(xc,k1,k2,k3,k4,k5,k6,k7,k8)
!call OK
!deallocate(c0)
!call OK
!deallocate(xc_ei)
!call OK

!!close(44)
!!close(45)
!!close(46)
!!close(48)
!!close(49)
!
!!call cpu_time(finish)
!
!write(6,*) finish - start

enddo !end loop number of systems

!!$OMP END DO

!$OMP END PARALLEL DO

!call cpu_time(finish)

!write(6,*) finish - start

!if ( vers .eq. 'singl') then
!call makeOutputSingle
!elseif ( vers .eq. 'dimer') then
!call makeOutputDimer
!elseif ( vers .eq. 'range') then
!call makeOutputRange
!elseif ( vers .eq. 'randm') then
!call makeOutputRandm
!endif

!write(outputdir,'(a7,a5,a1,i2,a1,i2)') "Output-", vers, "-", nint(aA*1d9*10), "-" , nint(aB*1d9*10) 
!
!call system("mkdir  " // outputdir)
!call system("mv *dat " // outputdir)
!call system("mv *txt " // outputdir)

!call system("mkdir output-`date +%x | sed 's/\//-/g'`-`date +%r | sed 's/:/-/g'`")
!call system("mv *dat `ls -lrth | tail -n 1 | awk '{print $9}'`")

!deallocate(E,diffe,diffh,minEe,minEh)

end
