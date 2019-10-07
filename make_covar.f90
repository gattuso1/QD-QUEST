program covar

implicit none
character*64 :: cov, powspec_l
integer :: i,j,n,k,l,nsys,ncep,cep
real*4 :: dummy
real*4,allocatable :: S(:,:),omega(:),Scov(:,:),Si(:,:,:)

ncep=100
n=20000
nsys = 50

allocate(Si(n,nsys,ncep),S(n,nsys),omega(n),Scov(n,n))
open(15,file='Allcov.dat')

do cep=1,ncep
do k=1,nsys

write(powspec_l,'(a16,i3.3,a1,i3.3,a4)') 'powspec-dip-num-', k, '-',cep,'.dat'

open(10,file=powspec_l)

do i=1,n
read(10,*) omega(i), Si(i,k,cep)
enddo

close(10)

enddo
enddo

do k=1,nsys
do i=1,n
S(i,k) = sum(Si(i,k,:))
enddo
enddo

open(13,file='Allcov-p.dat')
open(14,file='Allcov-m.dat')

do i=1,n,10
do j=1,n,10

do k=1,nsys
do l=1,nsys
Scov(i,j) = Scov(i,j) + S(i,k)*S(j,l)
!if ( ( i .eq. j ) .and.  ( abs(S(i,k)*S(j,l)) .ge. 1.*10**7) ) then
!write(6,*) omega(i), omega(j), S(i,k)*S(j,l)
!endif
enddo
enddo

enddo
enddo

do i=1,n,10
do j=1,n,10
if ( abs(Scov(i,j)) .ge. 1.*10**15 ) then
if ( Scov(i,j) .ge. 0. ) then
write(13,*) omega(i), omega(j), Scov(i,j)
elseif ( Scov(i,j) .le. 0. ) then
write(14,*) omega(i), omega(j), Scov(i,j)
endif
endif
write(15,*) omega(i), omega(j), Scov(i,j)
enddo
write(15,*) 
enddo


end 
