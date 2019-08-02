write(form_mat,'("(",i0,"ES16.5E2)")') nstates
write(form_TDM,'("(",i0,"f16.8)")') nstates
if (n .le. nQDA+nQDB) then
write(form_arr,'("(f16.8,16x,",i0,"f16.8)")') nstates+1
elseif (n .gt. nQDA+nQDB) then
write(form_arr,'("(",i0,"f16.8)")') nstates+3
endif
write(form_abs,'("(2f16.8,2x,i0)")') 
write(form_pop,'("(",i0,"ES16.8E3,ES20.12E2)")') nstates+1
write(form_com,'("(ES12.5E3,",i0,"ES16.8E3)")') nstates+1

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
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
if ( inbox .eq. "y" ) then
write(Tmat_x_f,form_TDM) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_y_f,form_TDM) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_z_f,form_TDM) (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif
enddo

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

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

!!!Opens output files
write(popc,'(a5,i5.5,a4)') 'Popc-', n, '.dat'
write(norm,'(a5,i5.5,a4)') 'Norm-', n, '.dat'
write(norm_ei,'(a8,i5.5,a4)') 'Norm_ei-', n, '.dat'
write(popc_ei,'(a8,i5.5,a4)') 'Popc_ei-', n, '.dat'
write(Re_c,'(a5,i5.5,a4)') 'Re_c-', n, '.dat'
write(Im_c,'(a5,i5.5,a4)') 'Im_c-', n, '.dat'
write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
open(newunit=popc_0_f   ,file=popc)    !; popc_0_f = 44
open(newunit=popc_ei_f  ,file=popc_ei) !; popc_ei_f = 49
open(newunit=norm_0_f   ,file=norm)    !; norm_0_f = 46
open(newunit=norm_ei_f  ,file=norm_ei) !; norm_ei_f = 50
open(newunit=Re_c_ei_f  ,file=Re_c_ei) !; Re_c_ei_f = 52
open(newunit=Im_c_ei_f  ,file=Im_c_ei) !; Im_c_ei_f = 53
open(newunit=Re_c_0_f   ,file=Re_c)    !; Re_c_0_f = 54
open(newunit=Im_c_0_f   ,file=Im_c)    !; Im_c_0_f = 55

!!!!!INITIAL POPULATIONS
c0    = 0.d0
c0(0) = 1.d0
xc0 = dcmplx(c0,0.0d0)
xc = dcmplx(0.d0,0.0d0)
xc_ei = dcmplx(0.d0,0.0d0)
xc(:,0) = xc0(:)
xc_ei(:,0) = xc0(:)

call RK_0_ei

endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei)
