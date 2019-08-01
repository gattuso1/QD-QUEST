write(form_mat,'("(",i0,"f16.8)")') nstates
write(form_arr,'("(",i0,"f16.8)")') nstates+1
write(form_abs,'("(2f16.8,2x,i0)")') 

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
write(matrices(i),*) n , aR(n) 
enddo

do i=0,nstates-1
write(H_0_f    ,form_mat) (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_dir_f  ,form_mat) (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_ex_f   ,form_mat) (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_JK_f   ,form_mat) ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
write(Tmat_0_f ,form_mat) (TransHam(i,j), j=0,nstates-1)
write(H_ei_f   ,form_mat) (Ham_ei(i,j), j=0,nstates-1)
write(Tmat_ei_f,form_mat) (TransHam_ei(i,j), j=0,nstates-1)
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
if ( inbox .eq. "y" ) then
write(Tmat_x_f,form_mat) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_y_f,form_mat) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_z_f,form_mat) (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif
enddo

write(Tmat_0_f,*) 
write(H_0_f,*) 
write(H_ei_f,*) 
write(Etr_0_f,form_arr) aR(n)*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_ei_f,form_arr) aR(n)*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)

endif

if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) ) then

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
