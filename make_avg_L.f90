program covar

implicit none
character*64 :: cov, powspec_l
integer :: i,j,n,k,l,nsys,ncep,cep,nstates
real*4 :: dummy
real*4,allocatable :: S(:,:),omega(:),Scov(:,:),Si(:,:,:)

n=1500
nsys = 100
nstates= 9

allocate(Si(n,nsys,nstates),S(n,nstates),omega(n),Scov(n,n))
open(15,file='Re_c_L_avg.dat')

do k=1,nsys

write(powspec_l,'(a7,i5.5,a4)') 'Re_c_L-', k, '.dat'

open(10,file=powspec_l)

do i=1,n
read(10,*) omega(i), (Si(i,k,j), j=1,nstates)
enddo

close(10)

enddo

do i=1,n
do j=1,nstates
S(i,j) = sum(Si(i,:,j))
enddo
enddo 

do i=1,n
write(15,*) omega(i), ( S(i,j)/nsys, j=1,nstates)
enddo




end 
