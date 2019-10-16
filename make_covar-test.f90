program covar

implicit none
character*64 :: cov, powspec_l
integer :: i,j,n,k,l,nsys,ncep,cep
real*4 :: dummy
real*4,allocatable :: omega(:),Scov(:,:),Si(:)

n=5071

allocate(Si(n),omega(n),Scov(n,n))
open(15,file='Allcov-test.dat')
open(10,file='TransAbs.dat')

do i=1,n
read(10,*) omega(i), Si(i)
enddo

do k=1,n,5
do l=1,n,5
Scov(k,l) = Si(k)*Si(l)
enddo
enddo

do i=1,n,5
do j=1,n,5
write(15,*) omega(i), omega(j), Scov(i,j)
enddo
write(15,*) 
enddo


end 
