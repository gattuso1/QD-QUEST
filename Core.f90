write(form_mat,'("(",i0,"ES16.5E2)")') nstates
write(form_TDM,'("(",i0,"f16.8)")') nstates
if (n .le. nQDA+nQDB) then
write(form_arr,'("(f16.8,16x,",i0,"f14.8)")') nstates+1
elseif (n .gt. nQDA+nQDB) then
write(form_arr,'("(",i0,"f16.12)")') nstates+3
endif
write(form_abs,'("(2f16.8,2x,i0)")') 
write(form_pop,'("(",i0,"ES17.8E3,ES25.16E3)")') nstates+1
write(form_com,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+1
write(form_com_L,'("(ES12.5E3,",i0,"ES18.8E3,ES28.15)")') nstates2+1
write(form_pop_L,'("(ES12.5E3,",i0,"ES18.8E3,f28.15)")') nstates
write(form_DipSpec,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+2

!!!!!Creates sum of TDM matrices and lambda vectors for avg dynamics
if ( Dyn_avg .eq. "y" ) then
Ham0_avg     = Ham0_avg     + Ham
TransHam_avg = TransHam_avg + TransHam
endif

!write(6,*) (TransHam(i,0)*Energ_au/elec, i=0,nstates-1)

if ( get_ei .eq. 'y' ) then
Ham_ei = Ham
allocate(lambda(0:nstates-1),source = 0.e0_dp)
allocate(iwork2(3+5*nstates),source=0)
allocate(work1(6*nstates),work2(1+6*nstates+2*nstates*nstates),source=0.e0_dp)
lworku=1+6*nstates+2*nstates*nstates
liworku=3+5*nstates
ierr=0
call dsyevd('v','u',nstates, Ham_ei(0:nstates-1,0:nstates-1),nstates,lambda,work2,lworku,iwork2,liworku,ierr)
deallocate(work1)
deallocate(work2)
deallocate(iwork2)

call make_Ham_l

!!!Make eigenstate TDM
if ( rdm_ori .eq. "n" ) then
Mat(:,:) = matmul(TransHam(:,:),Ham_ei(:,:))
TransHam_ei(:,:) = matmul(transpose(Ham_ei(:,:)),Mat(:,:))
elseif ( rdm_ori .eq. "y" ) then
Matx(:,:) = matmul(TransHam_l(:,:,1),Ham_ei(:,:))
Maty(:,:) = matmul(TransHam_l(:,:,2),Ham_ei(:,:))
Matz(:,:) = matmul(TransHam_l(:,:,3),Ham_ei(:,:))
TransHam_ei_l(:,:,1) = matmul(transpose(Ham_ei(:,:)),Matx(:,:))
TransHam_ei_l(:,:,2) = matmul(transpose(Ham_ei(:,:)),Maty(:,:))
TransHam_ei_l(:,:,3) = matmul(transpose(Ham_ei(:,:)),Matz(:,:))
TransHam_ei = sqrt(TransHam_ei_l(:,:,1)**2 + TransHam_ei_l(:,:,2)**2 + TransHam_ei_l(:,:,1)**2)
endif

if ( noMat .eq. "n" ) then
do i=1,size(matrices)
if (n .le. nQDA+nQDB) then
write(matrices(i),'(i0,f16.8)') n , aR(n)*1.d9
elseif (n .gt. nQDA+nQDB) then
write(matrices(i),'(i0,2f16.8)') n , aR(n)*1.d9, aR(n+ndim)*1.d9
endif
enddo
endif

do i=0,nstates-1
if ( noMat .eq. "n" ) then
write(H_0_f    ,form_mat) (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_dir_f  ,form_mat) (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_ex_f   ,form_mat) (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_JK_f   ,form_mat) ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
write(Tmat_0_f ,form_TDM) (TransHam(i,j), j=0,nstates-1)
write(H_ei_f   ,form_mat) (Ham_ei(i,j), j=0,nstates-1)
write(Tmat_ei_f,form_TDM) (TransHam_ei(i,j), j=0,nstates-1)
!if ( inbox .eq. "y" ) then
write(Tmat_x_f,form_TDM) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_y_f,form_TDM) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_z_f,form_TDM) (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif
!write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, &
!                         sqrt((TransHam_ei_l(0,i,1))**2+(TransHam_ei_l(0,i,2))**2+(TransHam_ei_l(0,i,3))**2) ,i
!elseif ( inbox .eq. "n" ) then
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
!endif
enddo

!!!!!labelling of eigenstates
Ham_ei = abs(Ham_ei)
do i=0,nstates-1
maxid(i) = maxval(Ham_ei(:,i))
enddo
do i = 0,nstates-1
do j = 0,nstates-1
if (Ham_ei(i,j) .eq. maxid(j)) then
        zero(j) = i
        exit
endif
enddo
enddo
write(label_0,*) n, (zero(j), j=0,nstates-1)

!!!!!!!!!!!!

if ( noMat .eq. "n" ) then
do i=1,size(matrices)
write(matrices(i),*)
enddo
endif

if (n .le. nQDA+nQDB) then
write(Etr_0_f,form_arr) aR(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_ei_f,form_arr) aR(n)*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)
elseif (n .gt. nQDA+nQDB) then
write(Etr_0_f,form_arr) aR(n)*1.d9, aR(n+ndim)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_ei_f,form_arr) aR(n)*1.d9, aR(n+ndim)*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)
endif

endif

if ( Dyn_L .eq. 'y' ) then

write(form1,'("(i4,i4,i4,i4,1x,100(e24.16,1x))")')
!write(form2,'("(100(f7.3,1x))")')

do i=0,nstates2-1
do j=0,nstates2-1
merge_diag(i,j)  = real(merge(1,0,i.eq.j),kind=dp)
merge_odiag(i,j) = real(merge(0,1,i.eq.j),kind=dp)
enddo
enddo

do i=0,nstates-1
do j=i+1,nstates-1
haml(i,j)=TransHam_ei(i,j)
haml(j,i)=haml(i,j)
enddo
enddo

do i=0,nstates-1
haml(i,i)=lambda(i)
!write(6,*) i, (haml(i,j),j=0,nstates-1)
enddo

!!!!!Liouvillian
do k=0,nstates-1
do l=0,nstates-1
do i=0,nstates-1
!!!!!Commutator [H,Eij]
xliou(k,l,i,l) = xliou(k,l,i,l) + dcmplx(haml(i,k),0._dp)
!print*, k,l,i,l,Ham_l(i,k),xliou(k,l,i,l)
enddo
do j=0,nstates-1
xliou(k,l,k,j) = xliou(k,l,k,j) - dcmplx(haml(l,j),0._dp)
enddo
enddo
enddo

!!!!!Print xLiouvillian
do i=0,nstates-1
do j=0,nstates-1
do k=0,nstates-1
do l=0,nstates-1
write(Liou_f,form1) i,j,k,l,xliou(i,j,k,l)
!print*, i,j,k,l,xliou(i,j,k,l)
enddo
enddo
enddo
enddo

!do i=0,nstates-1
!do j=0,nstates-1
!write(9,form2) ((xliou(i,j,k,l),l=0,nstates-1),k=0,nstates-1)
!write(*,form2) ((xliou(i,j,k,l),l=0,nstates-1),k=0,nstates-1)
!enddo
!enddo

kl=0;kc=0

!!!!Renumber xLiou
kl=-1
do i=0,nstates-1
do j=0,nstates-1
kl=kl+1
irow(kl,1)=i
irow(kl,2)=j
kc=-1
do k=0,nstates-1
do l=0,nstates-1
kc=kc+1
icol(kc,1)=k
icol(kc,2)=l
lfield(kl,kc)=xliou(i,j,k,l)
enddo
enddo
enddo
enddo

!lfield(0,6) = -1.e0_dp*lfield(0,6)
!lfield(6,0) = lfield(0,6)

!do l=0,nstates2-1
!write(6,'(i2, 9f18.12)') i, (lfield(l,j), j=0,nstates2-1)
!enddo

!xlfield = dcmplx(lfield,0.e0_dp)

endif

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) .or. ( Dyn_L .eq. 'y' ) ) then

!!!Opens output files
if ( nofiles .eq. 'n' ) then
if ( ( Dyn_0 .eq. 'y' ) ) then
write(popc,'(a5,i5.5,a4)') 'Popc-', n, '.dat'
write(Re_c,'(a5,i5.5,a4)') 'Re_c-', n, '.dat'
write(Im_c,'(a5,i5.5,a4)') 'Im_c-', n, '.dat'
open(newunit=popc_0_f   ,file=popc)    !; popc_0_f = 44
open(newunit=Re_c_0_f   ,file=Re_c)    !; Re_c_0_f = 54
open(newunit=Im_c_0_f   ,file=Im_c)    !; Im_c_0_f = 55
endif
if ( ( Dyn_ei .eq. 'y' ) ) then
write(popc_ei,'(a8,i5.5,a4)') 'Popc_ei-', n, '.dat'
write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
open(newunit=popc_ei_f  ,file=popc_ei) !; popc_ei_f = 49
open(newunit=Re_c_ei_f  ,file=Re_c_ei) !; Re_c_ei_f = 52
open(newunit=Im_c_ei_f  ,file=Im_c_ei) !; Im_c_ei_f = 53
endif
if ( ( Dyn_L .eq. 'y' ) ) then
write(Re_c_L,'(a7,i5.5,a4)') 'Re_c_L-', n, '.dat'
write(Im_c_L,'(a7,i5.5,a4)') 'Im_c_L-', n, '.dat'
write(Pop_c_L,'(a8,i5.5,a4)') 'Pop_c_L-', n, '.dat'
open(newunit=Re_c_L_f   ,file=Re_c_L)    !; Im_c_0_f = 55
open(newunit=Im_c_L_f   ,file=Im_c_L)    !; Im_c_0_f = 55
open(newunit=Pop_c_L_f  ,file=Pop_c_L)    !; Im_c_0_f = 55
endif
endif

!!!!!INITIAL POPULATIONS
c0    = 0.e0_dp
c0(0) = 1.e0_dp
xc0 = dcmplx(c0,0.0e0_dp)
xc = dcmplx(0.e0_dp,0.0e0_dp)
xc_ei = dcmplx(0.e0_dp,0.0e0_dp)
xc_L = dcmplx(0.e0_dp,0.0e0_dp)
xc(:,0) = xc0(:)
xc_ei(:,0) = xc0(:)
xc_L(0,0) = dcmplx(1.e0_dp,0.0e0_dp)

!do i=0,nstates2-1
!write(6,*) i, xc_L(i,0), xc_L(i,1)
!enddo
!xc_rho(0,0,0) = 1.e0_dp

call RK_0_ei

if ( nofiles .eq. 'y' ) then
if ( ( Dyn_0 .eq. 'y' ) ) then
close(popc_0_f   ,status="delete")
close(Re_c_0_f   ,status="delete")
close(Im_c_0_f   ,status="delete")
endif
if ( ( Dyn_ei .eq. 'y' ) ) then
close(popc_ei_f  ,status="delete")
close(Re_c_ei_f  ,status="delete")
close(Im_c_ei_f  ,status="delete")
endif
if ( ( Dyn_L .eq. 'y' ) ) then
close(Re_c_L_f   ,status="delete")
close(Im_c_L_f   ,status="delete")
endif
elseif ( nofiles .eq. 'n' ) then
if ( ( Dyn_0 .eq. 'y' ) ) then
close(popc_0_f )
close(Re_c_0_f )
close(Im_c_0_f )
endif
if ( ( Dyn_ei .eq. 'y' ) ) then
close(Re_c_ei_f)
close(Im_c_ei_f)
close(popc_ei_f)
endif
if ( ( Dyn_L .eq. 'y' ) ) then
close(Re_c_L_f )
close(Im_c_L_f )
endif
endif

endif

if (singleDS .eq. 'y' ) then
write(cov,'(a8,i0,a4)') 'DipSpec-', n, '.dat'
open(DipSpec_s,file=cov) 
endif

do t=0,ntime

time = t*timestep

powtemp = 0._dp
!pow_pol_conv = dcmplx(0.d0,0.d0)

do j=0,nstates-1
pow_s(n,t) = pow_s(n,t) + 2._dp * sum(TransHam_ei(j,:) * dreal(dconjg(xc_ei(j,t))*xc_ei(:,t)))
enddo

pow(t) = pow(t) + pow_s(n,t)

if ( inbox .eq. 'y' ) then
!do pol=1,npol
!pow_pol(pol,t) = pow_pol(pol,t) + dcmplx(pow_s(n,t),0.d0)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(pol)*Pe1(:)+l2(pol)*Pe2(:)+l3(pol)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!enddo

!Pow_pol(7,t) = pow_pol(7,t) + dcmplx(pow_s(n,t),0._dp)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(7)*Pe1(:)+l2(7)*Pe2(:)+l3(7)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!Pow_pol(8,t) = pow_pol(8,t) + dcmplx(pow_s(n,t),0._dp)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(8)*Pe1(:)+l2(8)*Pe2(:)+l3(8)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!Pow_pol(33,t) = pow_pol(33,t) + dcmplx(pow_s(n,t),0._dp)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(33)*Pe1(:)+l2(33)*Pe2(:)+l3(33)*Pe3(:),Dcenter(n,:)),0.e0_dp))

pow_pol(39,t) = pow_pol(39,t) + dcmplx(pow_s(n,t),0._dp)*&
   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(39)*Pe1(:)+l2(39)*Pe2(:)+l3(39)*Pe3(:),Dcenter(n,:)),0.e0_dp))

pow_pol(41,t) = pow_pol(41,t) + dcmplx(pow_s(n,t),0._dp)*&
   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(41)*Pe1(:)+l2(41)*Pe2(:)+l3(41)*Pe3(:),Dcenter(n,:)),0.e0_dp))

pow_pol(43,t) = pow_pol(43,t) + dcmplx(pow_s(n,t),0._dp)*&
   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))

pow_pol(44,t) = pow_pol(44,t) + dcmplx(pow_s(n,t),0._dp)*&
   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))

!off diagonal convergence of PM signal
!do pol=1,npol
!if ( pol .ne. 41) then
!!pow_pol_conv(t) = pow_pol_conv(t) + dcmplx(pow_s(n,t),0.d0)*&
!pow_pol_conv(t) = pow_pol_conv(t) + & !exp(-im*dcmplx(dot_product(l1(pol)*k_1(:)+l2(pol)*k_2(:)+l3(pol)*k_3(:),Dcenter(n,:)),0.e0_dp))*&
!                  exp(-im*dcmplx(dot_product((l1(pol)*Pe1(:)+l2(pol)*Pe2(:)+l3(pol)*Pe3(:))-&
!                                             (l1(41) *Pe1(:)+l2(41) *Pe2(:)+l3(41) *Pe3(:)),Dcenter(n,:)),0.e0_dp))
!endif
!enddo

!write(6,*) pow_pol(41,t), pow_pol_conv(t)

endif

if (singleDS .eq. 'y' ) then
write(DipSpec_s,form_DipSpec) time*t_au, pow_s(n,t) 
endif

enddo

!write(6,'(11a8)') 'Re', 'Im', 'Dx' , 'Dy', 'Dz', 'Pex', 'Pey', 'Pez', 'l1', 'l2', 'l3'
!write(6,'(11f8.4,2x,i5)') &
!dreal(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!dimag(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!Dcenter(n,:)*10**5, Pe1(:), l1(43), l2(43), l3(43), n
!write(6,'(11f8.4,2x,i5)') &
!dreal(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!dimag(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!Dcenter(n,:)*10**5, Pe1(:), l1(44), l2(44), l3(44), n

if (singleDS .eq. 'y' ) then
close(DipSpec_s)
endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei,haml)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei,xc_L,xc0,pop)
!deallocate(k1_L,k2_L,k3_L,k4_L,k5_L,k6_L,k7_L,k8_L,xliou,lfield)
deallocate(merge_diag,merge_odiag,icol,irow,maxid,zero)
!deallocate(xc_rho,k1_rho,k2_rho,k3_rho,k4_rho,k5_rho,k6_rho,k7_rho,k8_rho)
