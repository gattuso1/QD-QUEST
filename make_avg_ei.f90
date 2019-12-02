program covar

implicit none
character*64 :: cov, powspec_r, powspec_i
integer :: i,j,n,k,l,nsys,ncep,cep,nstates
real*8 :: dummy
real*8,allocatable :: S(:,:),omega(:),Scov(:,:),Si(:,:,:)
real*8,allocatable :: Im_c(:,:,:), Re_c(:,:,:)
complex*8 :: im

im = dcmplx(0.d0,1.d0)
n=10000
nsys = 4000
nstates= 9

allocate(Si(n,nsys,nstates),S(n,nstates),omega(n),Scov(n,n),Re_c(n,nsys,nstates), Im_c(n,nsys,nstates))
open(15,file='Cohe_ei_avg.dat')

do k=1,nsys

!write(powspec_l,'(a8,i5.5,a4)') 'Popc_ei-', k, '.dat'
write(powspec_r,'(a8,i5.5,a4)') 'Re_c_ei-', k, '.dat'
write(powspec_i,'(a8,i5.5,a4)') 'Im_c_ei-', k, '.dat'

open(10,file=powspec_r)
open(11,file=powspec_i)

do i=1,n
read(10,*) omega(i), (Re_c(i,k,j), j=1,nstates)
read(11,*) omega(i), (Im_c(i,k,j), j=1,nstates)
enddo

close(10)
close(11)

enddo

do i=1,n
do j=1,nstates
S(i,j) =  S(i,j) + sum(real((Re_c(i,:,j) - im*Im_c(i,:,j))*(Re_c(i,:,1) + im*Im_c(i,:,1))))
enddo
enddo

do i=1,n
write(15,*) omega(i), ( S(i,j)/nsys, j=1,nstates)
enddo




end 
