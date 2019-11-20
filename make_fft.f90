program covar

use Constants_au

implicit none
character*64 :: map
integer :: ntime,k,l,nsys,ncep,cep,nstates,FTpow,t,tp,t1
real*8 :: dummy, FTscale, timestep, time, totaltime, timestep2, FTscale2, ptime, t0
real*8,allocatable :: S(:,:),omega(:),Scov(:,:),Si(:,:,:)
real*8,allocatable :: TA(:), TA_gaus(:)
complex*16 :: im
complex*16,allocatable :: TA_n(:)

im = dcmplx(0.d0,1.d0)
ntime=100000
nsys = 510
nstates= 3
FTpow = 24
timestep = 1.d-17
totaltime = 1000.d-15
ptime = 200.d-15
t0 = (totaltime-ptime)/2.d0
tp = nint((totaltime-ptime)/timestep) 
ntime = nint(totaltime/timestep)

!write(6,*) (ntime-tp), ntime, (totaltime-ptime)/2.d0

allocate(TA(ntime),TA_gaus(ntime),omega(ntime))
allocate(TA_n(nint(2.d0**FTpow)))

open(11,file='test_fft.dat')

open(10,file='DipSpec.dat')

FTscale = h/(elec*(2.d0**FTpow)*timestep)

TA = 0.d0

do t=1,(ntime-tp)
read(10,*)
enddo

t1 = 1
do t=(ntime-tp+1),ntime
read(10,*) omega(t), TA(t1)
t1 = t1+1
enddo

do t=1,ntime
write(6,*) t, TA(t)
enddo

do t=1,ntime
time = t*timestep
TA_gaus(t)=exp(-1.d0*((time-(totaltime-ptime)/2.d0))**2.d0/(2.d0*(tp*timestep/10.d0)**2.d0))*TA(t)
enddo

!do t=1,ntime
!time = t*timestep
!write(6,*) time, TA_gaus(t) 
!enddo

TA_n = dcmplx(0.d0,0.d0)

do t=1,ntime
TA_n(t)  = dcmplx(TA_gaus(t),0.d0)
enddo
do t=ntime+1,nint(2.d0**FTpow)
TA_n(t)  = dcmplx(0.d0,0.d0)
enddo

!do t=1,ntime
!write(6,*) t, dreal(TA_n(t))
!enddo

write(6,*) 

call fft(TA_n)

t=0
do while ( t*FTscale .le. 4.d0 )
write(11,*)  t*FTscale, dreal(TA_n(t)), dimag(TA_n(t))
t = t + 10
enddo

contains

recursive subroutine fft(x)
  complex(kind=dp), dimension(:), intent(inout)  :: x
  complex(kind=dp)                               :: t
  integer                                        :: N
  integer                                        :: i
  complex(kind=dp), dimension(:), allocatable    :: even, odd

  N=size(x)

  if(N .le. 1) return

  allocate(odd((N+1)/2))
  allocate(even(N/2))

  ! divide
  odd =x(1:N:2)
  even=x(2:N:2)

  ! conquer
  call fft(odd)
  call fft(even)

  ! combine
  do i=1,N/2
     t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
     x(i)     = odd(i) + t
     x(i+N/2) = odd(i) - t
  end do

  deallocate(odd)
  deallocate(even)

end subroutine fft

end
