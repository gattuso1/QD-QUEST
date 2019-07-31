module Variables_au

use omp_lib
use Constants_au
use Normal

implicit none

   character*5 :: vers
   character*6 :: pgeom
   character*64 :: popc, hmti, norm, tdmM, hmt0, outputdir, norm_ei, popc_ei, Re_c_ei, Im_c_ei, integ,model, line, dummy
   character*64 :: Re_c, Im_c, syst_n
   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, dyn, hamilt, get_ei, finest, get_sp
   character*1 :: TDM_ee, Dyn_0, Dyn_ei, inbox
   integer :: syst, ndots, n, rmin, rmax, nsys, npulses, nstates, ntime, i, j, t, lwork, info, idlink, threads, nQDA, nQDB
   integer :: nhomoA,nhomoB,nhetero,totsys,ndim,nQD
   integer,allocatable :: seed(:)
   real(dp) :: aA, aB, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, slope, V0eV, minr, maxr, rsteps, side, link
   real(dp) :: vertex, zbase, start, finish
   real(dp) :: dispQD, displink, rdmlinker, rdmQDA, rdmQDB, t01, t02, t03, timestep, totaltime, distQD
   real(dp) :: omega01, omega02, omega03, phase01, phase02, phase03, width01, width02, width03, Ed01, Ed02, Ed03
   real(dp) :: pulse1, pulse2, pulse3, test, time, cnorm, cnormabs, cnormconj, cnorm2, cnorm2_ei, Kas, Kbs, Kcs, Dso1, Dso2, Dxf
   real(dp),allocatable :: aR(:), aRA(:), aRB(:), epsin(:), epsR(:), V0e(:), V0h(:), linker(:)
   real(dp),allocatable :: epsinA(:), epsinB(:), epsRA(:), epsRB(:), V0eA(:), V0eB(:), V0hA(:), V0hB(:)
   real(dp),allocatable :: Eeh1(:), Eeh2(:), Cb_eh1(:), Cb_eh2(:), Norm_Ana_e(:), Norm_Ana_h1(:), Norm_Ana_h2(:)
   real(dp),allocatable :: OverlapAna_h1e(:), OverlapAna_h2e(:), Cb_Num_eh1(:), Cb_Num_eh1_eh2(:), Cb_Num_eh2(:)
   real(dp),allocatable :: minEe(:,:),minEh(:,:), TransDip_Num_h1e(:), TransDip_Num_h2e(:), work(:), lambda(:)
   real(dp),allocatable :: TransDip_Ana_h1e(:), TransDip_Ana_h2e(:), Oscillator_Ana_h1e(:), Oscillator_Ana_h2e(:), Transvec(:)
   real(dp),allocatable :: ExctCoef_h1e(:), ExctCoef_h2e(:), Ham(:,:), E0(:), c0(:), TransHam(:,:), Hamt(:,:,:), c(:,:)
   real(dp),allocatable :: TransMat_ei(:,:), TransHam0(:,:), Ham_0(:), Ham_dir(:,:), Ham_ex(:,:), Ham_ei(:,:), Ham_l(:,:)
   real(dp),allocatable :: TransDip_Ana_h1h2(:), TransHam_ei(:,:), Mat(:,:), QDcoor(:,:), Dcenter(:,:), Pe1(:), Pe2(:), Pe3(:)
   real(dp),allocatable :: TransHam_d(:,:,:), TransHam_l(:,:,:), TransHam_ei_l(:,:,:), k_1(:), k_2(:), k_3(:)
   real(dp),allocatable :: Matx(:,:), Maty(:,:), Matz(:,:)
   complex(kind=8) :: ct1, ct2, ct3, ct4, xt01, xt02, xt03, xhbar, im, xwidth, xomega , xEd, xh, xphase, xtime, xhbar_au
   complex(kind=8),allocatable :: xHam(:,:) , xHamt(:,:,:), xTransHam(:,:), xE0(:), xHamtk2(:,:,:), xHamtk3(:,:,:), xHamtk4(:,:,:)
   complex(kind=8),allocatable :: xc0(:), xc(:,:), xc_ei(:,:), xcnew(:,:), k1(:), k2(:), k3(:) , k4(:), xHam_ei(:,:)
   complex(kind=8),allocatable :: dk1(:), dk2(:), dk3(:) , dk4(:), k5(:), k6(:), k7(:) , k8(:) 
   complex(kind=8),allocatable :: xc_ei_av(:,:), xctemp(:)

contains 

subroutine getVariables

NAMELIST /version/    syst 
NAMELIST /outputs/    inbox,get_sp,get_ei,Dyn_0,Dyn_ei,hamilt,fineSt,TDM_ee
NAMELIST /elecSt/     model,me,mh,eps,epsout,V0eV,omegaLO,slope,side
NAMELIST /fineStruc/  Kas,Kbs,Kcs,Dso1,Dso2,Dxf
NAMELIST /pulses/     integ,npulses,t01,t02,t03,timestep,totaltime,omega01,omega02,omega03,phase01,phase02,phase03,&
                      width01,width02,width03,Ed01,Ed02,Ed03,pgeom,vertex
NAMELIST /syst/       nQDA,nQDB,nhomoA,nhomoB,nhetero,dispQD,idlink,aA,aB     


open(150,file='QD_quest.def',form='formatted')
read(150,NML=version)
read(150,NML=outputs)
read(150,NML=elecSt)
read(150,NML=pulses)
read(150,NML=syst)

im         = dcmplx(0.0d0,1.0d0)
me         = me*m0
mh         = mh*m0
rhoe       = 1.0d0/sqrt((2.d0*me*omegaLO)/hbar)
rhoh       = 1.0d0/sqrt((2.d0*mh*omegaLO)/hbar)
V0         = V0eV*elec

!if ( model .eq. "SB") then
!nstates = 9
!elseif ( model .eq. "FO") then
!nstates = 4
!elseif ( model .eq. "FS") then
!nstates = 25
!endif

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

timestep   =  timestep*1.d-15/t_au  !timestep*1.d-15/t_au
totaltime  =  totaltime*1.d-15/t_au !totaltime*1.d-15/t_au
t01        =  t01*1.d-15/t_au       !t01*1.d-15/t_au
t02        =  t02*1.d-15/t_au       !t02*1.d-15/t_au
t03        =  t03*1.d-15/t_au       !t03*1.d-15/t_au
width01      =  width01*1.d-15/t_au     !width*1.d-15/t_au
width02      =  width02*1.d-15/t_au     !width*1.d-15/t_au
width03      =  width03*1.d-15/t_au     !width*1.d-15/t_au
omega01      =  omega01*t_au      !omega*1.d15*t_au
omega02      =  omega02*t_au      !omega*1.d15*t_au
omega03      =  omega03*t_au      !omega*1.d15*t_au
Ed01         =  Ed01/E_au        !0.024 !Ed/E_au
Ed02         =  Ed02/E_au        !0.024 !Ed/E_au
Ed03         =  Ed03/E_au        !0.024 !Ed/E_au
xh         =  dcmplx(timestep,0.0d0)
ntime      =  nint(totaltime/timestep)
phase01      =  pi
phase02      =  pi
phase03      =  pi

allocate(k_1(3))
allocate(k_2(3))
allocate(k_3(3))
allocate(Pe1(3))
allocate(Pe2(3))
allocate(Pe3(3))

if ( npulses .eq. 3) then
pulse1 = 1.d0 
pulse2 = 1.d0
pulse3 = 1.d0
elseif ( npulses .eq. 2) then
pulse1 = 1.d0
pulse2 = 1.d0
pulse3 = 0.d0
elseif ( npulses .eq. 1) then
pulse1 = 1.d0
pulse2 = 0.d0
pulse3 = 0.d0
elseif ( npulses .eq. 0) then
pulse1 = 0.d0
pulse2 = 0.d0
pulse3 = 0.d0
endif

if ( inbox .eq. 'y' ) then

if ( pgeom .eq. 'boxcar' ) then

zbase = 1000*sqrt(2.d0)*sqrt(1-cos(pi*vertex/180.d0))/sqrt(cos(pi*vertex/180.d0))
k_1(1) =(-2.d0 * pi / (cl/(omega01/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_1(2) = (2.d0 * pi / (cl/(omega01/t_au))) * ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_1(3) =(-2.d0 * pi / (cl/(omega01/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_2(1) =(-2.d0 * pi / (cl/(omega02/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))  
k_2(2) = (2.d0 * pi / (cl/(omega02/t_au))) * ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_2(3) = (2.d0 * pi / (cl/(omega02/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_3(1) = (2.d0 * pi / (cl/(omega03/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_3(2) = (2.d0 * pi / (cl/(omega03/t_au))) * ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
k_3(3) = (2.d0 * pi / (cl/(omega03/t_au))) * ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe1(1) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe1(2) =  ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe1(3) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe2(1) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))  
Pe2(2) =  ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe2(3) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe3(1) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe3(2) =  ( 1000.d0      ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))
Pe3(3) =  ( zbase / 2.d0 ) / (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase / 2.d0)**2))

elseif ( pgeom .eq. 'triang' ) then

zbase = 1000*sqrt(2.d0*(1.d0-cos(pi*vertex/180.d0))/(1.d0-(1-cos(pi*vertex/180.d0))/(1-cos(pi*120.d0/180.d0)))) 
k_1(1) = (2.d0*pi/(cl/(omega01/t_au)))*(zbase / 2.d0 ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_1(2) = (2.d0*pi/(cl/(omega01/t_au)))*(1000.d0      ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_1(3) = (2.d0*pi/(cl/(omega01/t_au)))*(zbase/(2.d0*sqrt(3.d0)))/(sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_2(1) =(-2.d0*pi/(cl/(omega02/t_au)))*(zbase / 2.d0 ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_2(2) = (2.d0*pi/(cl/(omega02/t_au)))*(1000.d0      ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_2(3) = (2.d0*pi/(cl/(omega02/t_au)))*(zbase/(2.d0*sqrt(3.d0)))/(sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_3(1) = (2.d0*pi/(cl/(omega03/t_au)))*0.d0 
k_3(2) = (2.d0*pi/(cl/(omega03/t_au)))*(1000.d0      ) /        (sqrt((1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
k_3(3) =(-2.d0*pi/(cl/(omega03/t_au)))*(zbase/(2.d0*sqrt(3.d0)))/(sqrt((1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe1(1) =  ( zbase / 2.d0 ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe1(2) =  ( 1000.d0      ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe1(3) =  ( zbase/(2.d0*sqrt(3.d0)))/(sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe2(1) =  ( zbase / 2.d0 ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe2(2) =  ( 1000.d0      ) /        (sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe2(3) =  ( zbase/(2.d0*sqrt(3.d0)))/(sqrt((zbase / 2.d0)**2+(1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe3(1) =  0.d0 
Pe3(2) =  ( 1000.d0      ) /        (sqrt((1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))
Pe3(3) =  ( zbase/(2.d0*sqrt(3.d0)))/(sqrt((1000.d0)**2+(zbase/(2.d0*sqrt(3.d0)))**2))

endif

endif

endif


!write(6,*) "You are requesting me to tackle a random set of heterodimer QD"
!
nQD    = nQDA+nQDB
ndim   = nhomoA+nhomoB+nhetero
nsys   = nQDA+nQDB+nhomoA+nhomoB+nhetero     
totsys = nsys+nhomoA+nhomoB+nhetero

allocate(aR(totsys))
allocate(linker(nsys))
allocate(epsin(totsys))
allocate(epsR(totsys))
allocate(V0e(totsys))
allocate(V0h(totsys))
!allocate(QDcoor(2*nsys,3))
!allocate(Dcenter(2*nsys,3))

if ( idlink .eq. 20 ) then
linker = 0.2d-9
elseif ( idlink .eq. 55 ) then
linker = 0.55d-9
endif
!
!
call random_seed(size = n)
allocate(seed(n))
call random_seed(get=seed)
!
!if ( get_sp .eq. 'y' ) then
!call system("mv Etransitions-he_0.dat tmp.dat ")
!open(11,file="tmp.dat")
!endif
!
aR = 0.d0

do n=1,nQDA
aR(n) = r8_NORMAL_AB(aA,dispQD*1d-9,seed(1))
enddo
do n=nQDA+1,nQDA+nQDB
aR(n) = r8_NORMAL_AB(aB,dispQD*1d-9,seed(2))
enddo
do n=nQDA+nQDB+1,nQDA+nQDB+nhomoA
aR(n) = r8_NORMAL_AB(aA,dispQD*1d-9,seed(1))
aR(n+ndim) = r8_NORMAL_AB(aA,dispQD*1d-9,seed(2))
enddo
do n=nQDA+nQDB+nhomoA+1,nQDA+nQDB+nhomoA+nhomoB
aR(n) = r8_NORMAL_AB(aB,dispQD*1d-9,seed(1))
aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*1d-9,seed(2))
enddo
do n=nQDA+nQDB+nhomoA+nhomoB+1,nQDA+nQDB+nhomoA+nhomoB+nhetero
aR(n) = r8_NORMAL_AB(aA,dispQD*1d-9,seed(1))
aR(n+ndim) = r8_NORMAL_AB(aB,dispQD*1d-9,seed(2))
enddo

!elseif ( get_sp .eq. 'y' ) then
!read(11,*) aR(n), aR(n+nsys)
!write(6,*) aR(n), aR(n+nsys)
!aR(n) = aR(n)*1.e-9_dp
!aR(n+nsys) = aR(n+nsys)*1.e-9_dp
!endif

do n=1,totsys
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n)= 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36.d0/35)*aR(n)/rhoe)+exp(-(36.d0/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo


if ( inbox .eq. 'y' ) then



!open(56,file='box-dimers.xyz',form='formatted',action='read')
!read(56,*) 
!read(56,*) 
!do n = 1,nsys
!read(56,*) dummy, QDcoor(n,1), QDcoor(n,2), QDcoor(n,3)
!read(56,*) dummy, QDcoor(n+ndim,1), QDcoor(n+ndim,2), QDcoor(n+ndim,3)
!Dcenter(n,1) = (QDcoor(n,1) + QDcoor(n+ndim,1))/2.
!Dcenter(n,2) = (QDcoor(n,2) + QDcoor(n+ndim,2))/2.
!Dcenter(n,3) = (QDcoor(n,3) + QDcoor(n+ndim,3))/2.
!enddo
!QDcoor(:,:) = QDcoor(:,:) * 1.e-10_dp
!Dcenter(:,:) = Dcenter(:,:) * 1.e-10_dp

endif

end subroutine getVariables
   
end module Variables_au
