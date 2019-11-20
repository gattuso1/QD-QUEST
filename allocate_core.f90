nstates2=nstates**2

allocate(TransHam(0:nstates-1,0:nstates-1),&
         TransHam_ei_l(0:nstates-1,0:nstates-1,3),&
         TransHam_l(0:nstates-1,0:nstates-1,3),&
         TransHam_d(0:nstates-1,0:nstates-1,3),&
         TransHam_ei(0:nstates-1,0:nstates-1),&
         Transvec(0:nstates-1),&
         TransMat_ei(0:nstates-1,0:nstates-1),&
         Mat(0:nstates-1,0:nstates-1),&
         Matx(0:nstates-1,0:nstates-1),&
         Maty(0:nstates-1,0:nstates-1),&
         Matz(0:nstates-1,0:nstates-1),&
         Ham(0:nstates-1,0:nstates-1),&
         Ham_l(0:nstates-1,0:nstates-1),&
         haml(0:nstates-1,0:nstates-1),&
         pop(0:nstates-1,0:ntime+1),&
         Ham_0(0:nstates-1),&
         Ham_dir(0:nstates-1,0:nstates-1),&
         Ham_ex(0:nstates-1,0:nstates-1),&
         Ham_ei(0:nstates-1,0:nstates-1),&
         maxid(0:nstates-1),source=0.e0_dp)

allocate(xc(0:nstates-1,0:ntime+1),&
         xc0(0:nstates-1),&
         xc_ei(0:nstates-1,0:ntime+1),&
         xc_L(0:nstates2-1,0:ntime+1),&
         c0(0:nstates-1),k1(0:nstates-1),&
         k2(0:nstates-1),k3(0:nstates-1),&
         k4(0:nstates-1),k5(0:nstates-1),&
         k6(0:nstates-1),k7(0:nstates-1),&
         k8(0:nstates-1),&
         k1_L(0:nstates2-1),&
         k2_L(0:nstates2-1),k3_L(0:nstates2-1),&
         k4_L(0:nstates2-1),k5_L(0:nstates2-1),&
         k6_L(0:nstates2-1),k7_L(0:nstates2-1),&
         k8_L(0:nstates2-1))

allocate(xc_rho(0:nstates-1,0:nstates-1,0:ntime+1),&
         k1_rho(0:nstates-1,0:nstates-1),&
         k2_rho(0:nstates-1,0:nstates-1),k3_rho(0:nstates-1,0:nstates-1),&
         k4_rho(0:nstates-1,0:nstates-1),k5_rho(0:nstates-1,0:nstates-1),&
         k6_rho(0:nstates-1,0:nstates-1),k7_rho(0:nstates-1,0:nstates-1),&
         k8_rho(0:nstates-1,0:nstates-1))

allocate(merge_diag(0:nstates2-1,0:nstates2-1),&
         merge_odiag(0:nstates2-1,0:nstates2-1), source = 0.e0_dp )

allocate(xliou(0:nstates-1,0:nstates-1,0:nstates-1,0:nstates-1),lfield(0:nstates2-1,0:nstates2-1))

xliou = dcmplx(0.e0_dp,0.e0_dp)
lfield = dcmplx(0.e0_dp,0.e0_dp)


allocate(irow(0:nstates2-1,2),icol(0:nstates2-1,2),zero(0:nstates-1),source=0)

