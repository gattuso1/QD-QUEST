write(form_mat,'("(",i0,"ES16.5E2)")') nstates
write(form_TDM,'("(",i0,"f16.8)")') nstates
if (n .le. nQDA+nQDB) then
write(form_arr,'("(f16.8,16x,",i0,"f16.8)")') nstates+1
elseif (n .gt. nQDA+nQDB) then
write(form_arr,'("(",i0,"f16.8)")') nstates+3
endif
write(form_abs,'("(2f16.8,2x,i0)")') 
write(form_pop,'("(",i0,"ES17.8E3,ES25.16E3)")') nstates+1
write(form_com,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+1
write(form_com_L,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates2+2
write(form_DipSpec,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+1

if ( get_ei .eq. 'y' ) then
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

do i=1,size(matrices)
if (n .le. nQDA+nQDB) then
write(matrices(i),'(i0,f16.8)') n , aR(n)*1.d9
elseif (n .gt. nQDA+nQDB) then
write(matrices(i),'(i0,2f16.8)') n , aR(n)*1.d9, aR(n+ndim)*1.d9
endif
enddo

do i=0,nstates-1
write(H_0_f    ,form_mat) (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_dir_f  ,form_mat) (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_ex_f   ,form_mat) (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_JK_f   ,form_mat) ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
write(Tmat_0_f ,form_TDM) (TransHam(i,j), j=0,nstates-1)
write(H_ei_f   ,form_mat) (Ham_ei(i,j), j=0,nstates-1)
write(Tmat_ei_f,form_TDM) (TransHam_ei(i,j), j=0,nstates-1)
if ( inbox .eq. "y" ) then
write(Tmat_x_f,form_TDM) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_y_f,form_TDM) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_z_f,form_TDM) (TransHam_ei_l(i,j,3), j=0,nstates-1)
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, &
                         sqrt((TransHam_ei_l(0,i,1))**2+(TransHam_ei_l(0,i,2))**2+(TransHam_ei_l(0,i,3))**2) ,i
elseif ( inbox .eq. "n" ) then
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
endif
enddo

!call Convolution

do i=1,size(matrices)
write(matrices(i),*)
enddo

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
merge_diag(i,j)  = merge(1,0,i.eq.j)
merge_odiag(i,j) = merge(0,1,i.eq.j)
enddo
enddo

do i=0,nstates-1
do j=i+1,nstates-1
haml(i,j)=TransHam_ei(i,j)
haml(j,i)=TransHam_ei(j,i)
enddo
enddo

do i=0,nstates-1
haml(i,i)=lambda(i) 
enddo

!!!!!Liouvillian
do k=0,nstates-1
do l=0,nstates-1
do i=0,nstates-1
!!!!!Commutator [H,Eij]
xliou(k,l,i,l) = xliou(k,l,i,l) + haml(i,k)
!print*, k,l,i,l,Ham_l(i,k),xliou(k,l,i,l)
enddo
do j=0,nstates-1
xliou(k,l,k,j) = xliou(k,l,k,j) - haml(l,j)
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

xlfield(:,:) = dcmplx(lfield,0.e0_dp)

endif

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

!!!Opens output files
write(popc,'(a5,i5.5,a4)') 'Popc-', n, '.dat'
write(popc_ei,'(a8,i5.5,a4)') 'Popc_ei-', n, '.dat'
write(Re_c,'(a5,i5.5,a4)') 'Re_c-', n, '.dat'
write(Im_c,'(a5,i5.5,a4)') 'Im_c-', n, '.dat'
write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
write(Re_c_L,'(a7,i5.5,a4)') 'Re_c_L-', n, '.dat'
write(Im_c_L,'(a7,i5.5,a4)') 'Im_c_L-', n, '.dat'
open(newunit=popc_0_f   ,file=popc)    !; popc_0_f = 44
open(newunit=popc_ei_f  ,file=popc_ei) !; popc_ei_f = 49
open(newunit=Re_c_ei_f  ,file=Re_c_ei) !; Re_c_ei_f = 52
open(newunit=Im_c_ei_f  ,file=Im_c_ei) !; Im_c_ei_f = 53
open(newunit=Re_c_0_f   ,file=Re_c)    !; Re_c_0_f = 54
open(newunit=Im_c_0_f   ,file=Im_c)    !; Im_c_0_f = 55
open(newunit=Re_c_L_f   ,file=Re_c_L)    !; Im_c_0_f = 55
open(newunit=Im_c_L_f   ,file=Im_c_L)    !; Im_c_0_f = 55

!!!!!INITIAL POPULATIONS
c0    = 0.e0_dp
c0(0) = 1.e0_dp
xc0 = dcmplx(c0,0.0e0_dp)
xc = dcmplx(0.e0_dp,0.0e0_dp)
xc_ei = dcmplx(0.e0_dp,0.0e0_dp)
xc(:,0) = xc0(:)
xc_ei(:,0) = xc0(:)
xc_L(:,0) = xc0(:)

call RK_0_ei

if ( nofiles .eq. 'y' ) then
close(popc_0_f   ,status="delete")
close(popc_ei_f  ,status="delete")
close(norm_0_f   ,status="delete")
close(norm_ei_f  ,status="delete")
close(Re_c_ei_f  ,status="delete")
close(Im_c_ei_f  ,status="delete")
close(Re_c_0_f   ,status="delete")
close(Im_c_0_f   ,status="delete")
close(Re_c_L_f   ,status="delete")
close(Im_c_L_f   ,status="delete")
endif

endif

if ( doFT .eq. 'y' ) then

if (singleFT .eq. 'y' ) then
write(cov,'(a8,i0,a4)') 'DipSpec-', n, '.dat'
open(DipSpec_s,file=cov) 
endif

do t=0,ntime

time = t*timestep

do j=0,nstates-1
powtemp = 2._dp * sum(TransHam_ei(j,:) * dreal(dconjg(xc_ei(j,t))*xc_ei(:,t)))
enddo

pow(t) = pow(t) + powtemp

if (singleFT .eq. 'y' ) then
write(DipSpec_s,form_DipSpec) time, powtemp 
endif

enddo

if (singleFT .eq. 'y' ) then
close(DipSpec_s)
endif

endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei,haml)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei,xc_L,xc0)
deallocate(k1_L,k2_L,k3_L,k4_L,k5_L,k6_L,k7_L,k8_L)
deallocate(merge_diag,merge_odiag,icol,irow,xliou,lfield,xlfield)
