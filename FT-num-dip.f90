program FT

implicit none

   integer,  parameter :: dp = SELECTED_REAL_KIND(15,307)
   integer, allocatable:: zero(:)
   real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
   character*128 :: line, Re_c_ei, Im_c_ei
   integer :: t,t2, tmax,wmin,wmax,error,nn,i,j,t0ei,nsamples,nstates,n, en_id, en, k
   real(dp) :: normt, normw, time0ei, timestep, w1,w2, wstep, w, elec, h, pow_ave, dummy, en_i, en_f, en_step, tottime
   real(dp),allocatable :: time(:), cohe(:), pow(:,:), pow_gaus(:), E(:,:), raA(:), raB(:), Ediff(:,:,:), c_ei(:,:)
   real(dp),allocatable :: maxid(:), TransDip(:,:,:), TransDip_ei(:,:,:),pow_tot(:), pow_tot2(:), pulse(:)
   real(dp),allocatable :: Re_c(:,:), Im_c(:,:), sum_pow(:,:), Ham_ei(:,:,:), pow_gaus2(:), p1(:), p2(:), p3(:)
   complex(8),allocatable :: xtime(:), xcohe(:),  wft(:), xpow_gaus(:),xpow_tot(:), xpulse(:), wftp(:), wftf(:), xpow(:,:)
   complex(8) :: im, xwstep, xpow_ave

im = dcmplx(0.0d0,1.0d0)

nstates=9
!tottime = 10000.d-15
timestep = 0.10d-15
time0ei = 400.d-15
!t0ei = nint(time0ei/(timestep))
t0ei = 4001
tmax = 104001
w1 = 0.01d12
w2 = 1000.d12
h  = 6.62607004d-34
elec  = 1.60217662d-19

nn=1000

!open(10,file='Etransitions-he_ei.dat')
open(13,file='powspec-dip-num.dat')
open(14,file='pow_gaus.dat')
open(20,file='pow_gaus2.dat')
open(15,file='Ham_ei-all.dat')
open(16,file='TransMat-all.dat')
open(17,file='Pulse.dat')
open(18,file='dipspec.dat')
!open(17,file='TransDipFT2.dat')

!do
!   read(10,'(A128)',iostat=error) line
!   if (error .ne. 0) exit
!   nn = nn + 1
!enddo

allocate(zero(0:nstates-1))
allocate(maxid(0:nstates-1))
allocate(raA(nn))
allocate(raB(nn))
allocate(E(nn,0:nstates-1))
allocate(c_ei(nn,0:nstates-1))
allocate(TransDip(nn,0:nstates-1,0:nstates-1))
allocate(TransDip_ei(nn,0:nstates-1,0:nstates-1))
allocate(Re_c(nn,0:nstates-1))
allocate(Im_c(nn,0:nstates-1))
allocate(Ediff(nn,0:nstates-1,0:nstates-1))
allocate(time(0:tmax+1))
allocate(pow(nn,0:tmax+1))
allocate(xpow(nn,0:tmax+1))
allocate(pow_tot(t0ei:tmax+1))
allocate(pow_tot2(t0ei:tmax+1))
allocate(pow_gaus(0:tmax+1-t0ei))
allocate(pow_gaus2(0:tmax+1-t0ei))
allocate(xpow_gaus(0:tmax+1-t0ei))
allocate(Ham_ei(nn,0:nstates-1,0:nstates-1))
allocate(wft(0:tmax+1-t0ei))
allocate(p1(0:tmax+1))
allocate(p2(0:tmax+1))
allocate(p3(0:tmax+1))
allocate(wftp(0:tmax+1))
allocate(wftf(0:tmax+1))

do n=1,nn

zero = 0
maxid = 0.d0

read(15,'(A128)')
do i=0,nstates-1
read(15,'(A128)') line
read(line,*) ((Ham_ei(n,i,j)),j=0,nstates-1)
!write(6,'(8f8.4)') (Ham_ei(n,i,j),j=1,nstates-1)
enddo
read(15,'(A128)')

!Ham_ei = abs(Ham_ei)
!
!do i=0,nstates-1
!maxid(i) = maxval(Ham_ei(n,:,i))
!enddo
!
!do i = 0,nstates-1
!do j = 0,nstates-1
!if (Ham_ei(n,i,j) .eq. maxid(j)) then
!        zero(j) = i
!        exit
!endif
!enddo
!enddo

read(16,'(A128)')
do i=0,nstates-1
read(16,'(A128)') line
read(line,*) ((TransDip(n,i,j)),j=0,nstates-1)
enddo
read(16,'(A128)')

TransDip_ei(n,:,:) = matmul(Ham_ei(n,:,:),TransDip(n,:,:))

write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
open(11,file=Re_c_ei)
open(12,file=Im_c_ei)

do i=1,t0ei-1
read(11,'(A128)') line
read(12,'(A128)') line
enddo

do i=t0ei,tmax
read(11,*)  dummy, (Re_c(n,j), j=0,nstates-1)
read(12,*)  dummy, (Im_c(n,j), j=0,nstates-1)

do j=0,nstates-1
do k=0,nstates-1
!pow(n,i) = pow(n,i) + 2*TransDip_ei(n,j,k) * (Re_c(n,j)*Re_c(n,k) + Im_c(n,j)*Im_c(n,k))
!pow(n,i) = pow(n,i) + 2*TransDip_ei(n,j,k) * (Re_c(n,j)*Re_c(n,k) + Im_c(n,j)*Im_c(n,k))
pow(n,i) = pow(n,i) + 2*TransDip_ei(n,j,k) * dreal((Re_c(n,j) - im*Im_c(n,j))*(Re_c(n,k) + im*Im_c(n,k)))
                         !                      im*(Re_c(n,j)*Im_c(n,i)-Re_c(n,i)*Im_c(n,j)))
enddo
enddo

enddo

enddo

do i=t0ei,tmax
j = i - t0ei
pow_tot(j) = sum(pow(:,i))
enddo

!!!NUMERICAL FT

read(17,*)

do i=0,tmax+1-t0ei
read(17,*) time(i), p1(i), p2(i), p3(i)
enddo

do i=0,tmax-t0ei
   pow_gaus(i) =  exp(-1.d0*((i-(tmax-t0ei)/2.d0)*timestep)**2.d0&      
                  /(2.d0*((tmax-t0ei)*timestep/20.d0)**2.d0)) * (pow_tot(i)/nn)
write(14,*) i*timestep, pow_gaus(i), pow_tot(i)
enddo

do i=0,(tmax-t0ei)/2
j = i + (tmax-t0ei)/2
pow_gaus2(i) = pow_gaus(j)
write(20,*) i*timestep, pow_gaus2(i)
enddo

do i=(tmax-t0ei)/2,tmax-t0ei
j = i - (tmax-t0ei)/2
pow_gaus2(i) = pow_gaus(j)
write(20,*) i*timestep, pow_gaus2(i)
enddo

pulse = p1 + p2 + p3

wstep = (w2-w1)/(tmax-t0ei)
xpow_gaus  = dcmplx(pow_gaus,0.d0)
xpulse  = dcmplx(pulse,0.d0)
!xpow_gaus  = dcmplx(pow_gaus,0.d0)
w = w1 

write(6,*) w1, w2, wstep

!do t=0,tmax-t0ei
!do t2=0,tmax-t0ei
!wft(t) = wft(t) + exp(-2.d0*pi*im*w*t2*timestep) * xpow_gaus(t2) 
!wftp(t) = wftp(t) + exp(-1.d0*im*w*t2*timestep) * xpulse(t2) 
!enddo
!w = w + wstep
!enddo

w = w1

do t=0,tmax-t0ei
wftf(t) = -2.d0 * aimag(sqrt(dreal(wft(t))**2+aimag(wft(t))**2) * dconjg(wftp(t)))
!wftf(t) = -2.d0 * aimag(dreal(wft(t))**2 * dconjg(wftp(t)))
!write(13,'(7ES14.4E3)') w*h/elec, dreal(wft(t)), aimag(wft(t)),&
!                        sqrt(dreal(wft(t))**2 + aimag(wft(t))**2), &
!                        sqrt(dreal(wftp(t))**2 + aimag(wftp(t))**2), & 
!                        wftf(t) 
w = w + wstep
enddo

end
