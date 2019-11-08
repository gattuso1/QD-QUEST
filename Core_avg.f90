write(form_TDM,'("(",i0,"f16.8)")') nstates
write(form_arr,'("(",i0,"f16.8)")') nstates+3
write(form_com_L,'("(ES12.5E3,",i0,"ES18.8E3,ES28.15)")') nstates2+1
write(form_pop_L,'("(ES12.5E3,",i0,"ES18.8E3,f28.15)")') nstates
write(form_DipSpec,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+2
write(form1,'("(i4,i4,i4,i4,1x,100(e24.16,1x))")')

if ( noMat .eq. "n" ) then
do i=1,size(matrices_avg)
write(matrices_avg(i),*)  aR_avgA*1.d9, aR_avgB*1.d9
enddo

do i=0,nstates-1
write( Tmat_avg_f,form_TDM) (TransHam_ei(i,j)    , j=0,nstates-1)
write(Tmat_avgx_f,form_TDM) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_avgy_f,form_TDM) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_avgz_f,form_TDM) (TransHam_ei_l(i,j,3), j=0,nstates-1)
enddo

do i=1,size(matrices_avg)
write(matrices_avg(i),*)
enddo
endif

write(Etr_avg_f,form_arr) aR_avgA*1.d9, aR_avgB*1.d9, (Ham_ei(i,i)*Energ_au/elec, i=0,nstates-1)

!do i=0,nstates2-1
!do j=0,nstates2-1
!merge_diag(i,j)  = real(merge(1,0,i.eq.j),kind=dp)
!merge_odiag(i,j) = real(merge(0,1,i.eq.j),kind=dp)
!enddo
!enddo
!
!do i=0,nstates-1
!do j=i+1,nstates-1
!haml(i,j)=TransHam_avg(i,j)
!haml(j,i)=haml(i,j)
!enddo
!enddo
!
!do i=0,nstates-1
!haml(i,i)=lambda_avg(i)
!enddo
!
!!!!!!Liouvillian
!do k=0,nstates-1
!do l=0,nstates-1
!do i=0,nstates-1
!!!!!!Commutator [H,Eij]
!xliou(k,l,i,l) = xliou(k,l,i,l) + dcmplx(haml(i,k),0._dp)
!enddo
!do j=0,nstates-1
!xliou(k,l,k,j) = xliou(k,l,k,j) - dcmplx(haml(l,j),0._dp)
!enddo
!enddo
!enddo
!
!!!!!!!Print xLiouvillian
!!do i=0,nstates-1
!!do j=0,nstates-1
!!do k=0,nstates-1
!!do l=0,nstates-1
!!write(Liou_f,form1) i,j,k,l,xliou(i,j,k,l)
!!!print*, i,j,k,l,xliou(i,j,k,l)
!!enddo
!!enddo
!!enddo
!!enddo
!
!kl=0;kc=0
!
!!!!!Renumber xLiou
!kl=-1
!do i=0,nstates-1
!do j=0,nstates-1
!kl=kl+1
!irow(kl,1)=i
!irow(kl,2)=j
!kc=-1
!do k=0,nstates-1
!do l=0,nstates-1
!kc=kc+1
!icol(kc,1)=k
!icol(kc,2)=l
!lfield(kl,kc)=xliou(i,j,k,l)
!enddo
!enddo
!enddo
!enddo

!!!!!INITIAL POPULATIONS
!xc_L = dcmplx(0.e0_dp,0.0e0_dp)
!xc_L(0,0) = dcmplx(1.e0_dp,0.0e0_dp)
xc_ei = dcmplx(0.e0_dp,0.0e0_dp)
xc_ei(0,0) = dcmplx(1.e0_dp,0.0e0_dp)

call RK_0_ei

!if (singleDS .eq. 'y' ) then
!write(cov,'(a8,i0,a4)') 'DipSpec-', n, '.dat'
!open(DipSpec_s,file=cov) 
!endif
!
!do t=0,ntime
!
!time = t*timestep
!
!powtemp = 0._dp
!!pow_pol_conv = dcmplx(0.d0,0.d0)
!
!do j=0,nstates-1
!pow_s(n,t) = pow_s(n,t) + 2._dp * sum(TransHam_ei(j,:) * dreal(dconjg(xc_ei(j,t))*xc_ei(:,t)))
!enddo
!
!pow(t) = pow(t) + pow_s(n,t)
!
!if ( inbox .eq. 'y' ) then
!!do pol=1,npol
!!pow_pol(pol,t) = pow_pol(pol,t) + dcmplx(pow_s(n,t),0.d0)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(pol)*Pe1(:)+l2(pol)*Pe2(:)+l3(pol)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!enddo
!
!!pow_pol(7,t) = pow_pol(7,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(7)*Pe1(:)+l2(7)*Pe2(:)+l3(7)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(8,t) = pow_pol(8,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(8)*Pe1(:)+l2(8)*Pe2(:)+l3(8)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(33,t) = pow_pol(33,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(33)*Pe1(:)+l2(33)*Pe2(:)+l3(33)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(39,t) = pow_pol(39,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(39)*Pe1(:)+l2(39)*Pe2(:)+l3(39)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(41,t) = pow_pol(41,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(41)*Pe1(:)+l2(41)*Pe2(:)+l3(41)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(43,t) = pow_pol(43,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!!pow_pol(44,t) = pow_pol(44,t) + dcmplx(pow_s(n,t),0._dp)*&
!!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!
!!off diagonal convergence of PM signal
!!do pol=1,npol
!!if ( pol .ne. 41) then
!!!pow_pol_conv(t) = pow_pol_conv(t) + dcmplx(pow_s(n,t),0.d0)*&
!!pow_pol_conv(t) = pow_pol_conv(t) + & !exp(-im*dcmplx(dot_product(l1(pol)*k_1(:)+l2(pol)*k_2(:)+l3(pol)*k_3(:),Dcenter(n,:)),0.e0_dp))*&
!!                  exp(-im*dcmplx(dot_product((l1(pol)*Pe1(:)+l2(pol)*Pe2(:)+l3(pol)*Pe3(:))-&
!!                                             (l1(41) *Pe1(:)+l2(41) *Pe2(:)+l3(41) *Pe3(:)),Dcenter(n,:)),0.e0_dp))
!!endif
!!enddo
!
!!write(6,*) pow_pol(41,t), pow_pol_conv(t)
!
!endif
!
!if (singleDS .eq. 'y' ) then
!write(DipSpec_s,form_DipSpec) time*t_au, pow_s(n,t) 
!endif
!
!enddo
!
!!write(6,'(11a8)') 'Re', 'Im', 'Dx' , 'Dy', 'Dz', 'Pex', 'Pey', 'Pez', 'l1', 'l2', 'l3'
!!write(6,'(11f8.4,2x,i5)') &
!!dreal(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!!dimag(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!!Dcenter(n,:)*10**5, Pe1(:), l1(43), l2(43), l3(43), n
!!write(6,'(11f8.4,2x,i5)') &
!!dreal(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!!dimag(exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0.e0_dp))),&
!!Dcenter(n,:)*10**5, Pe1(:), l1(44), l2(44), l3(44), n
!
!if (singleDS .eq. 'y' ) then
!close(DipSpec_s)
!endif






deallocate(haml,xc_L,k1_L,k2_L,k3_L,k4_L,k5_L,k6_L,k7_L,k8_L,merge_diag,merge_odiag,xliou,lfield,irow,icol)

