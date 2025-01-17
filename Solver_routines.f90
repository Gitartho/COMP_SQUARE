!********************* SOLVER ROUTINES *****************************************************
      SUBROUTINE UNSTEADY(stepl)
      use declare_variables
      implicit none	  
	  
      real Ucont,Vcont,Wcont
      real rhl, ul, vl, wl, pl, Tl, El, mul, vol
      real ixl, jxl, kxl, iyl, jyl, kyl, izl, jzl, kzl	
      real uil, ujl, ukl, vil, vjl, vkl, wil, wjl, wkl, Til, Tjl, Tkl
      real bx, by, bz, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, T_x, T_y, T_z
      real Txx, Txy, Txz, Tyy, Tyz, Tzz
      real div2b3, facprM, Sterm
      integer stepl, var	 	  
	  
	  do nbl = 1,nblocks
	  do k = 1,NK(nbl)
	  do j = 1,NJ(nbl)
	  do i = 1,NI(nbl) 
	  
	  ixl = ix(i,j,k,nbl)
	  jxl = jx(i,j,k,nbl)
	  kxl = kx(i,j,k,nbl)
	  iyl = iy(i,j,k,nbl)
	  jyl = jy(i,j,k,nbl)
	  kyl = ky(i,j,k,nbl)
	  izl = iz(i,j,k,nbl)
	  jzl = jz(i,j,k,nbl)
	  kzl = kz(i,j,k,nbl)
	  
	  rhl = Qp(i,j,k,nbl,1)
	  ul =  Qp(i,j,k,nbl,2)
	  vl =  Qp(i,j,k,nbl,3)
	  wl =  Qp(i,j,k,nbl,4)
      pl =  Qp(i,j,k,nbl,5)
	  Tl =  Qp(i,j,k,nbl,6)
	  
	  Ucont = ul*ixl + vl*iyl + wl*izl
	  Vcont = ul*jxl + vl*jyl + wl*jzl
	  Wcont = ul*kxl + vl*kyl + wl*kzl
	  
	  El = pl/(rhl*(gamma-1.d0))+0.5d0*(ul**2+vl**2+wl**2)
	  vol = 1.d0/Jac(i,j,k,nbl)	  

	  !*********** Inviscid flux estimation ***********************
	  Fflux(i,j,k,nbl,1) = - rhl*Ucont*vol
	  Fflux(i,j,k,nbl,2) = -(rhl*ul*Ucont+ixl*pl)*vol
	  Fflux(i,j,k,nbl,3) = -(rhl*vl*Ucont+iyl*pl)*vol
	  Fflux(i,j,k,nbl,4) = -(rhl*wl*Ucont+izl*pl)*vol
	  Fflux(i,j,k,nbl,5) = -(rhl*El*Ucont+Ucont*pl)*vol
	  
	  Gflux(i,j,k,nbl,1) = - rhl*Vcont*vol
	  Gflux(i,j,k,nbl,2) = -(rhl*ul*Vcont+jxl*pl)*vol
	  Gflux(i,j,k,nbl,3) = -(rhl*vl*Vcont+jyl*pl)*vol
	  Gflux(i,j,k,nbl,4) = -(rhl*wl*Vcont+jzl*pl)*vol
	  Gflux(i,j,k,nbl,5) = -(rhl*El*Vcont+Vcont*pl)*vol
	  
	  Hflux(i,j,k,nbl,1) = - rhl*Wcont*vol
	  Hflux(i,j,k,nbl,2) = -(rhl*ul*Wcont+kxl*pl)*vol
	  Hflux(i,j,k,nbl,3) = -(rhl*vl*Wcont+kyl*pl)*vol
	  Hflux(i,j,k,nbl,4) = -(rhl*wl*Wcont+kzl*pl)*vol
	  Hflux(i,j,k,nbl,5) = -(rhl*El*Wcont+Wcont*pl)*vol
	
	Enddo
	Enddo
	Enddo
	Enddo
	
  !*********** Viscous flux estimation ************************
	if(euler.eq.0) then
	
	if(dscheme.eq.1.or.dscheme.eq.2) then
	call DISCRETIZATION_I_EXP(Qp,Qpi,nprims)
	call DISCRETIZATION_J_EXP(Qp,Qpj,nprims)
		if(grid2d.ne.1) call DISCRETIZATION_K_EXP(Qp,Qpk,nprims)
		if(grid2d.eq.1) Qpk = 0.d0
	elseif(dscheme.eq.3.or.dscheme.eq.4) then
	call DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
	call DISCRETIZATION_J_COMP(Qp,Qpj,nprims)
		if(grid2d.ne.1) call DISCRETIZATION_K_COMP(Qp,Qpk,nprims)
		if(grid2d.eq.1) Qpk = 0.d0
	endif

	
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
	ixl = ix(i,j,k,nbl)
	jxl = jx(i,j,k,nbl)
	kxl = kx(i,j,k,nbl)
	iyl = iy(i,j,k,nbl)
	jyl = jy(i,j,k,nbl)
	kyl = ky(i,j,k,nbl)
	izl = iz(i,j,k,nbl)
	jzl = jz(i,j,k,nbl)
	kzl = kz(i,j,k,nbl)
	
	vol =1.d0/Jac(i,j,k,nbl)
	
	rhl = Qp(i,j,k,nbl,1)
	ul  = Qp(i,j,k,nbl,2)
	vl  = Qp(i,j,k,nbl,3)
	wl  = Qp(i,j,k,nbl,4)
	pl  = Qp(i,j,k,nbl,5)
	Tl  = Qp(i,j,k,nbl,6)
	
	uil  = Qpi(i,j,k,nbl,2)
	vil  = Qpi(i,j,k,nbl,3)
	wil  = Qpi(i,j,k,nbl,4)
	Til  = Qpi(i,j,k,nbl,6)
	
	ujl  = Qpj(i,j,k,nbl,2)
	vjl  = Qpj(i,j,k,nbl,3)
	wjl  = Qpj(i,j,k,nbl,4)
	Tjl  = Qpj(i,j,k,nbl,6)
	
	ukl  = Qpk(i,j,k,nbl,2)
	vkl  = Qpk(i,j,k,nbl,3)
	wkl  = Qpk(i,j,k,nbl,4)
	Tkl  = Qpk(i,j,k,nbl,6)
	
	u_x = uil*ixl + ujl*jxl +ukl*kxl
	u_y = uil*iyl + ujl*jyl +ukl*kyl
	u_z = uil*izl + ujl*jzl +ukl*kzl
	
	v_x = vil*ixl + vjl*jxl +vkl*kxl
	v_y = vil*iyl + vjl*jyl +vkl*kyl
	v_z = vil*izl + vjl*jzl +vkl*kzl
	
	w_x = wil*ixl + wjl*jxl +wkl*kxl
	w_y = wil*iyl + wjl*jyl +wkl*kyl
	w_z = wil*izl + wjl*jzl +wkl*kzl
	
	T_x = Til*ixl + Tjl*jxl +Tkl*kxl
	T_y = Til*iyl + Tjl*jyl +Tkl*kyl
	T_z = Til*izl + Tjl*jzl +Tkl*kzl
	
	mul = Tl**1.5*(1.d0 + 110.4d0/T_ref)/(Tl + 110.4d0/T_ref)
	div2b3 = 2.d0/3.d0*(u_x + v_y + w_z)
	facprM = mul/((gamma-1)*prandtl*Mach**2)
	
	Txx = mul*(2.d0*u_x - div2b3)    										!Reynolds stress tensor
	Tyy = mul*(2.d0*v_y - div2b3)
	Tzz = mul*(2.d0*w_z - div2b3)
	Txy = mul*(u_y + v_x)
	Txz = mul*(u_z + w_x)
	Tyz = mul*(v_z + w_y)
	
	bx = ul*Txx + vl*Txy + wl*Txz + facprM*T_x								!Heat Flux
	by = ul*Txy + vl*Tyy + wl*Tyz + facprM*T_y
	bz = ul*Txz + vl*Tyz + wl*Tzz + facprM*T_z
	
!*********** Viscous flux estimation ***********************
	!Fflux(i,j,k,nbl,1) = Fflux(i,j,k,nbl,1)
	Fflux(i,j,k,nbl,2) = Fflux(i,j,k,nbl,2) + vol/Re*(ixl*Txx + iyl*Txy + izl*Txz)
	Fflux(i,j,k,nbl,3) = Fflux(i,j,k,nbl,3) + vol/Re*(ixl*Txy + iyl*Tyy + izl*Tyz)
	Fflux(i,j,k,nbl,4) = Fflux(i,j,k,nbl,4) + vol/Re*(ixl*Txz + iyl*Tyz + izl*Tzz)
	Fflux(i,j,k,nbl,5) = Fflux(i,j,k,nbl,5) + vol/Re*(ixl*bx + iyl*by + izl*bz)
	
	!Gflux(i,j,k,nbl,1) = Gflux(i,j,k,nbl,1)
	Gflux(i,j,k,nbl,2) = Gflux(i,j,k,nbl,2) + vol/Re*(jxl*Txx + jyl*Txy + jzl*Txz)
	Gflux(i,j,k,nbl,3) = Gflux(i,j,k,nbl,3) + vol/Re*(jxl*Txy + jyl*Tyy + jzl*Tyz)
	Gflux(i,j,k,nbl,4) = Gflux(i,j,k,nbl,4) + vol/Re*(jxl*Txz + jyl*Tyz + jzl*Tzz)
	Gflux(i,j,k,nbl,5) = Gflux(i,j,k,nbl,5) + vol/Re*(jxl*bx +  jyl*by  + jzl*bz)

	!Hflux(i,j,k,nbl,1) = Hflux(i,j,k,nbl,1)
	Hflux(i,j,k,nbl,2) = Hflux(i,j,k,nbl,2) + vol/Re*(kxl*Txx + kyl*Txy + kzl*Txz)
	Hflux(i,j,k,nbl,3) = Hflux(i,j,k,nbl,3) + vol/Re*(kxl*Txy + kyl*Tyy + kzl*Tyz)
	Hflux(i,j,k,nbl,4) = Hflux(i,j,k,nbl,4) + vol/Re*(kxl*Txz + kyl*Tyz + kzl*Tzz)
	Hflux(i,j,k,nbl,5) = Hflux(i,j,k,nbl,5) + vol/Re*(kxl*bx  + kyl*by +  kzl*bz)
		

 !********* Post processing steps (if any) ************************************

	!if(stepl.eq.1)then
	if(stepl.eq.1.and.i.lt.NI(nbl).and.j.lt.NJ(nbl).and.k.lt.NK(nbl))then
	tke = tke + 0.5d0*(ul**2 + vl**2 + wl**2)*rhl
	enstpt = enstpt + (1.d0/2.d0)*((w_y - v_z)**2 + (u_z - w_x)**2 + (v_x - u_y)**2)*rhl
	endif
	
	
	Enddo
	Enddo
	Enddo
	Enddo
	
	endif
 
  
  
  !*************** ESTIMATING FLUX DERIVATIVES ***********************************		  	  
	if(dscheme.eq.1.or.dscheme.eq.2) then
	call DISCRETIZATION_I_EXP(Fflux,fluxD,nconserv)
	net_flux = fluxD															! Residual
	call DISCRETIZATION_J_EXP(Gflux,fluxD,nconserv)
	net_flux = net_flux + fluxD	 
		if(grid2d.ne.1) call DISCRETIZATION_K_EXP(Hflux,fluxD,nconserv)
		if(grid2d.eq.1) fluxD=0.d0
	net_flux = net_flux + fluxD
	
	elseif(dscheme.eq.3.or.dscheme.eq.4) then
	call DISCRETIZATION_I_COMP(Fflux,fluxD,nconserv)
	net_flux = fluxD
	call DISCRETIZATION_J_COMP(Gflux,fluxD,nconserv)
	net_flux = net_flux + fluxD	  
		if(grid2d.ne.1) call DISCRETIZATION_K_COMP(Hflux,fluxD,nconserv)
		if(grid2d.eq.1) fluxD=0.d0
	net_flux = net_flux + fluxD
	
	end if
	
	do var = 1,nconserv
	do nbl = 1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
	!*************** ESTIMATING PART OF NEW CONSERVATIVE VARIABLES ***********************************	
	vol = 1.d0/Jac(i,j,k,nbl)
	Qcnew(i,j,k,nbl,var) = Qcnew(i,j,k,nbl,var) + fac_RK(stepl)*net_flux(i,j,k,nbl,var)*time_step/vol
	!*************** ESTIMATING PART OF SUB-RK_STAGE CONSERVATIVE VARIABLE to estimate new net flux ***********************************
	if(stepl.lt.rk_steps) then
	Qc(i,j,k,nbl,var) = Qcini(i,j,k,nbl,var) + fac_qini(stepl+1)*net_flux(i,j,k,nbl,var)*time_step/vol
	else
	!*************** UPDATING CONSERVATIVE VARIABLES ***********************************
	Qc(i,j,k,nbl,var) = Qcnew(i,j,k,nbl,var)
	!*************** ESTIMATING MAXIMUM RESIDUAL ***********************************
	res(var) = max(res(var),abs(Qcnew(i,j,k,nbl,var)-Qcini(i,j,k,nbl,var)))   !max(loc(res(var),Qcnew(i,j,k,nbl,var)-Qcini(i,j,k,nbl,var)))-- to get location of max Value
	
	endif
	
	
	  
	
	Enddo
	Enddo
	Enddo
	Enddo
	Enddo
	END
	  
	SUBROUTINE SET_PRIMITIVES()
	use declare_variables	 
	implicit none	 
	
	real rhl, ul, vl, wl, pl, Tl, El	
		
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
		rhl = Qc(i,j,k,nbl,1)
		ul  = Qc(i,j,k,nbl,2)/rhl
		vl  = Qc(i,j,k,nbl,3)/rhl
		wl  = Qc(i,j,k,nbl,4)/rhl
		pl  = (Qc(i,j,k,nbl,5)/rhl - 0.5d0*(ul**2 + vl**2 + wl**2)) *(rhl*(gamma-1.d0))
		Tl  = gamma*mach**2*pl/rhl
		
		Qp(i,j,k,nbl,1) = rhl
		Qp(i,j,k,nbl,2) = ul
		Qp(i,j,k,nbl,3) = vl
		Qp(i,j,k,nbl,4) = wl
		Qp(i,j,k,nbl,5) = pl
		Qp(i,j,k,nbl,6) = Tl
	enddo
	enddo
	enddo
	enddo	
	END	 
	  
!******************************************************************************************************************	  
!*********************************** PERIODIC DISCRETIZATION ROUTINES FOR FLOW VARIABLES *******************************
!******************************************************************************************************************	  

	SUBROUTINE DISCRETIZATION_I_EXP(PHI,PHID,nvars)     !Discretization in i to get derivative using Explicit scheme
	use declare_variables	 
	implicit none	  
	
	integer var,ip1,ip2,im1,im2,nvars
	real bb4, ab2, LLx	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do var=1,nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 3,NI(nbl)-2
		
		ip2 = i+2
		im2 = i-2
		ip1 = i+1
		im1 = i-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
	
	enddo
	
		i=1
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i=2
		ip2 = 4
		im2 = NI(nbl)-1
		ip1 = 3
		im1 = NI(nbl)
		PHID(i,j,k,nbl,var) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i= NI(nbl)-1
		ip2 = 2
		im2 = NI(nbl)-3
		ip1 = 1
		im1 = NI(nbl)-2
		PHID(i,j,k,nbl,var) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i=NI(nbl)
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
	enddo
	enddo
	enddo
	enddo
	END	

	SUBROUTINE DISCRETIZATION_J_EXP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	 
	
	integer var,jp1,jp2,jm1,jm2,nvars,coeff,jl
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do var=1,nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do i = 1,NI(nbl)
	do j = 3,NJ(nbl)-2
		
		jp2 = j+2
		jm2 = j-2
		jp1 = j+1
		jm1 = j-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
	
	enddo
		
		j=1
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j=2
		jp2 = 4
		jm2 = NJ(nbl)-1
		jp1 = 3
		jm1 = NJ(nbl)
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j= NJ(nbl)-1
		jp2 = 2
		jm2 = NJ(nbl)-3
		jp1 = 1
		jm1 = NJ(nbl)-2
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j=NJ(nbl)
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
	enddo
	enddo
	enddo
	enddo
	END	
	
	SUBROUTINE DISCRETIZATION_K_EXP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,kp2,km1,km2,nvars
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	
	
	 bb4= bdisc/4.d0
	 ab2= adisc/2.d0
	 
	do var=1,nvars
	do nbl =1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	do k = 3,NK(nbl)-2
	 	 
		kp2 = k+2
		km2 = k-2
		kp1 = k+1
		km1 = k-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
	enddo
	 	 
		k=1
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k=2
		kp2 = 4
		km2 = NK(nbl)-1
		kp1 = 3
		km1 = NK(nbl)
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k   = NK(nbl)-1
		kp2 = 2
		km2 = NK(nbl)-3
		kp1 = 1
		km1 = NK(nbl)-2
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k   = NK(nbl)
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		PHID(i,j,k,nbl,var) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
	enddo
	enddo
	enddo
	enddo
	END	
	
	SUBROUTINE DISCRETIZATION_K2D_EXP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,km1,nvars
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	
	END		  
	
!****************************** COMPACT SCHEMES *******************************************************************

	SUBROUTINE DISCRETIZATION_I_COMP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	  
	
	integer var,ip1,ip2,im1,im2,nvars,nm1
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NImax) :: RHS
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do var=1,nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 3,NI(nbl)-2
		
		ip2 = i+2
		im2 = i-2
		ip1 = i+1
		im1 = i-1
		RHS(i) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
	
	enddo
	
		i=1
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		RHS(i) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i=2
		ip2 = 4
		im2 = NI(nbl)-1
		ip1 = 3
		im1 = NI(nbl)
		RHS(i) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i= NI(nbl)-1
		ip2 = 2
		im2 = NI(nbl)-3
		ip1 = 1
		im1 = NI(nbl)-2
		RHS(i) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		i=NI(nbl)
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		RHS(i) = bb4*(PHI(ip2,j,k,nbl,var)-PHI(im2,j,k,nbl,var)) + ab2*(PHI(ip1,j,k,nbl,var)-PHI(im1,j,k,nbl,var))
		
		!Do i = 1,NI(nbl)-1
		!print*,RHS(1:nm1)
		!enddo
		
		nm1 = NI(nbl)-1
		call TDMAP(1,NI(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(1:nm1,j,k,nbl,var) = RHS(1:nm1)
		PHID(NI(nbl),j,k,nbl,var) = RHS(1)
		!print*, PHI(i,j,k,nbl,var),APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1)
		!Do i = 1,NI(nbl)-1
		!print*,RHS(1:nm1)
		!enddo
		!pause
		
	enddo
	enddo
	enddo
	enddo
	END	
	
	SUBROUTINE DISCRETIZATION_J_COMP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	 

	
	integer var,jp1,jp2,jm1,jm2,nvars,nm1
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NJmax) :: RHS	  
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do var=1,nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do i = 1,NI(nbl)
	do j = 3,NJ(nbl)-2
		
		jp2 = j+2
		jm2 = j-2
		jp1 = j+1
		jm1 = j-1
		RHS(j) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
	
	enddo
		
		j=1
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		RHS(j)= bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j=2
		jp2 = 4
		jm2 = NJ(nbl)-1
		jp1 = 3
		jm1 = NJ(nbl)
		RHS(j)= bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j= NJ(nbl)-1
		jp2 = 2
		jm2 = NJ(nbl)-3
		jp1 = 1
		jm1 = NJ(nbl)-2
		RHS(j) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		j=NJ(nbl)
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		RHS(j) = bb4*(PHI(i,jp2,k,nbl,var)-PHI(i,jm2,k,nbl,var)) + ab2*(PHI(i,jp1,k,nbl,var)-PHI(i,jm1,k,nbl,var))
		
		nm1 = NJ(nbl)-1
		call TDMAP(1,NJ(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(i,1:nm1,k,nbl,var) = RHS(1:nm1)
		PHID(i,NJ(nbl),k,nbl,var) = RHS(1)
		
	enddo
	enddo
	enddo
	enddo
	END	
	
	SUBROUTINE DISCRETIZATION_K_COMP(PHI,PHID,nvars)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,kp2,km1,km2,nvars,nm1
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI,PHID
	real,dimension(NKmax) :: RHS	  
	
	 bb4= bdisc/4.d0
	 ab2= adisc/2.d0
	 
	do var=1,nvars
	do nbl =1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	do k = 3,NK(nbl)-2
	 	 
		kp2 = k+2
		km2 = k-2
		kp1 = k+1
		km1 = k-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
	enddo
	 	 
		k=1
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k=2
		kp2 = 4
		km2 = NK(nbl)-1
		kp1 = 3
		km1 = NK(nbl)
		RHS(k) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k   = NK(nbl)-1
		kp2 = 2
		km2 = NK(nbl)-3
		kp1 = 1
		km1 = NK(nbl)-2
		RHS(k) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		k   = NK(nbl)
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl,var)-PHI(i,j,km2,nbl,var)) + ab2*(PHI(i,j,kp1,nbl,var)-PHI(i,j,km1,nbl,var))
		
		nm1 = NK(nbl)-1
		call TDMAP(1,NK(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(i,j,1:nm1,nbl,var) = RHS(1:nm1)
		PHID(i,j,NK(nbl),nbl,var) = RHS(1)
		
	enddo
	enddo
	enddo
	enddo	
	END		  
	  
!******************************************************************************************************************	  
!***********************************PERIODIC DISCRETIZATION ROUTINES FOR GRID COORDINATES *************************
!******************************************************************************************************************

	SUBROUTINE DISCRETIZATION_I_EXP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	  
	
	integer var,ip1,ip2,im1,im2
	real bb4, ab2, LLx  
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	  
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 3,NI(nbl)-2
		
		ip2 = i+2
		im2 = i-2
		ip1 = i+1
		im1 = i-1
		PHID(i,j,k,nbl) = bb4*(PHI(ip2,j,k,nbl)-PHI(im2,j,k,nbl)) + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl))
	
	enddo
	
		LLx = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)
		
		i=1
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i=2
		ip2 = 4
		im2 = NI(nbl)-1
		ip1 = 3
		im1 = NI(nbl)
		PHID(i,j,k,nbl)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i= NI(nbl)-1
		ip2 = 2
		im2 = NI(nbl)-3
		ip1 = 1
		im1 = NI(nbl)-2
		PHID(i,j,k,nbl)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i=NI(nbl)
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		!PHID(1,j,k,nbl) 		= bb4*(PHI(3,j,k,nbl)-(PHI(NI(nbl)-2,j,k,nbl)- LLx))  + ab2*(PHI(2,j,k,nbl)-PHI(NI(nbl)-1,j,k,nbl)+LLx)
		!PHID(2,j,k,nbl) 		= bb4*(PHI(4,j,k,nbl)-PHI(NI(nbl)-1,j,k,nbl) + LLx)   + ab2*(PHI(3,j,k,nbl)-PHI(1,j,k,nbl))
		!PHID(NI(nbl)-1,j,k,nbl) = bb4*(PHI(2,j,k,nbl)-PHI(NI(nbl)-3,j,k,nbl) + LLx)   + ab2*(PHI(1,j,k,nbl)-PHI(NI(nbl)-2,j,k,nbl))
		!PHID(NI(nbl),j,k,nbl) 	= bb4*(PHI(3,j,k,nbl)-PHI(NI(nbl)-2,j,k,nbl) + LLx)   + ab2*(PHI(2,j,k,nbl)-PHI(NI(nbl)-1,j,k,nbl)+LLx)
	enddo
	enddo
	enddo	
	END	

	SUBROUTINE DISCRETIZATION_J_EXP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	 
	
	integer var,jp1,jp2,jm1,jm2
	real bb4, ab2, LLy 
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl = 1,nblocks
	do k   = 1,NK(nbl)
	do i   = 1,NI(nbl)
	do j   = 3,NJ(nbl)-2
		
		jp2 = j+2
		jm2 = j-2
		jp1 = j+1
		jm1 = j-1
		PHID(i,j,k,nbl) = bb4*(PHI(i,jp2,k,nbl)-PHI(i,jm2,k,nbl)) + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl))
		
	enddo
		
		LLy = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)
		
		j   = 1
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		j   = 2
		jp2 = 4
		jm2 = NJ(nbl)-1
		jp1 = 3
		jm1 = NJ(nbl)
		PHID(i,j,k,nbl)= bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		j   = NJ(nbl)-1
		jp2 = 2
		jm2 = NJ(nbl)-3
		jp1 = 1
		jm1 = NJ(nbl)-2
		PHID(i,j,k,nbl)= bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		j   = NJ(nbl)
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
	enddo
	enddo
	enddo
	END	

	SUBROUTINE DISCRETIZATION_K_EXP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,kp2,km1,km2
	real bb4, ab2, LLz
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID


	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl =1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	do k = 3,NK(nbl)-2
		
		kp2 = k+2
		km2 = k-2
		kp1 = k+1
		km1 = k-1
		PHID(i,j,k,nbl) = bb4*(PHI(i,j,kp2,nbl)-PHI(i,j,km2,nbl)) + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl))
	
	enddo
	
		LLz = PHI(i,J,NK(nbl),nbl) - PHI(i,j,1,nbl)
		
		k=1
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		k=2
		kp2 = 4
		km2 = NK(nbl)-1
		kp1 = 3
		km1 = NK(nbl)
		PHID(i,j,k,nbl)= bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		k   = NK(nbl)-1
		kp2 = 2
		km2 = NK(nbl)-3
		kp1 = 1
		km1 = NK(nbl)-2
		PHID(i,j,k,nbl)= bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		k   = NK(nbl)
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		PHID(i,j,k,nbl)= bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
	enddo
	enddo
	enddo	  
	END	

	SUBROUTINE DISCRETIZATION_K2D_EXP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,km1
	real bb4, ab2	  
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	
	
	END	

!****************************** COMPACT SCHEMES FOR GRID *******************************************************************

	SUBROUTINE DISCRETIZATION_I_COMP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	  
	
	integer var,ip1,ip2,im1,im2,nm1
	real bb4, ab2, LLx	  
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	real,dimension(NImax) :: RHS	  

	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
	do i = 3,NI(nbl)-2
		
		ip2 = i+2
		im2 = i-2
		ip1 = i+1
		im1 = i-1
		RHS(i) = bb4*(PHI(ip2,j,k,nbl)-PHI(im2,j,k,nbl)) + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl))
	
	enddo
		LLx = PHI(NI(nbl),j,k,nbl) - PHI(1,j,k,nbl)

		i=1
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		RHS(i)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i=2
		ip2 = 4
		im2 = NI(nbl)-1
		ip1 = 3
		im1 = NI(nbl)
		RHS(i)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i= NI(nbl)-1
		ip2 = 2
		im2 = NI(nbl)-3
		ip1 = 1
		im1 = NI(nbl)-2
		RHS(i)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		i=NI(nbl)
		ip2 = 3
		im2 = NI(nbl)-2
		ip1 = 2
		im1 = NI(nbl)-1
		RHS(i)= bb4*(PHI(ip2,j,k,nbl)-(PHI(im2,j,k,nbl)- LLx))  + ab2*(PHI(ip1,j,k,nbl)-PHI(im1,j,k,nbl)+LLx)
		
		nm1 = NI(nbl)-1
		call TDMAP(1,NI(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(1:nm1,j,k,nbl) = RHS(1:nm1)
		PHID(NI(nbl),j,k,nbl) = RHS(1)
	enddo
	enddo
	enddo
	END	

	SUBROUTINE DISCRETIZATION_J_COMP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	 
	
	integer var,jp1,jp2,jm1,jm2,nm1
	real bb4, ab2, LLy	  
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	real,dimension(NJmax) :: RHS	  
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do i = 1,NI(nbl)
	do j = 3,NJ(nbl)-2
		
		jp2 = j+2
		jm2 = j-2
		jp1 = j+1
		jm1 = j-1
		RHS(j) = bb4*(PHI(i,jp2,k,nbl)-PHI(i,jm2,k,nbl)) + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl))
	
	enddo

		LLy = PHI(i,NJ(nbl),k,nbl) - PHI(i,1,k,nbl)

		j=1
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		RHS(j) = bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		
		j=2
		jp2 = 4
		jm2 = NJ(nbl)-1
		jp1 = 3
		jm1 = NJ(nbl)
		RHS(j) = bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		
		j= NJ(nbl)-1
		jp2 = 2
		jm2 = NJ(nbl)-3
		jp1 = 1
		jm1 = NJ(nbl)-2
		RHS(j) = bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)
		
		j=NJ(nbl)
		jp2 = 3
		jm2 = NJ(nbl)-2
		jp1 = 2
		jm1 = NJ(nbl)-1
		RHS(j) = bb4*(PHI(i,jp2,k,nbl)-(PHI(i,jm2,k,nbl)- LLy))  + ab2*(PHI(i,jp1,k,nbl)-PHI(i,jm1,k,nbl)+LLy)	
		
		nm1 = NJ(nbl)-1
		call TDMAP(1,NJ(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(i,1:nm1,k,nbl) = RHS(1:nm1)
		PHID(i,NJ(nbl),k,nbl) = RHS(1)
	enddo
	enddo
	enddo	
	END	
	
	SUBROUTINE DISCRETIZATION_K_COMP_GRID(PHI,PHID)
	use declare_variables	 
	implicit none	  
	
	integer var,kp1,kp2,km1,km2,nvars,nm1
	real bb4, ab2, LLz
	real,dimension(NImax,NJmax,NKmax,nblocks) :: PHI,PHID
	real,dimension(NKmax) :: RHS	  
	
	bb4= bdisc/4.d0
	ab2= adisc/2.d0
	
	do nbl =1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	do k = 3,NK(nbl)-2
	
		kp2 = k+2
		km2 = k-2
		kp1 = k+1
		km1 = k-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl)-PHI(i,j,km2,nbl)) + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl))
	
	enddo
		LLz = PHI(i,j,NK(nbl),nbl) - PHI(i,j,1,nbl)

		k=1
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		
		k=2
		kp2 = 4
		km2 = NK(nbl)-1
		kp1 = 3
		km1 = NK(nbl)
		RHS(k) = bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		
		k= NK(nbl)-1
		kp2 = 2
		km2 = NK(nbl)-3
		kp1 = 1
		km1 = NK(nbl)-2
		RHS(k) = bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)
		
		k=NK(nbl)
		kp2 = 3
		km2 = NK(nbl)-2
		kp1 = 2
		km1 = NK(nbl)-1
		RHS(k) = bb4*(PHI(i,j,kp2,nbl)-(PHI(i,j,km2,nbl)- LLz))  + ab2*(PHI(i,j,kp1,nbl)-PHI(i,j,km1,nbl)+LLz)	

		nm1 = NK(nbl)-1
		call TDMAP(1,NK(nbl)-1,APD(1:nm1),ACD(1:nm1),AMD(1:nm1),RHS(1:nm1),nm1)
		PHID(i,j,1:nm1,nbl) = RHS(1:nm1)
		PHID(i,j,NK(nbl),nbl) = RHS(1)
	enddo
	enddo
	enddo
	END		  
	  	  
!*************************************************************************************
!**************************  FILTERING SUBROUTINES ***********************************
!*************************************************************************************

    SUBROUTINE FILTERING_I(PHI,nvars)
    use declare_variables	 
    implicit none	  
	
    integer var,ipc,imc,coeff,nvars,npts,nm1  
    real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
    real,dimension(NImax)	:: RHS
	real pert
	
	npts = fscheme/2
	
	do var = 1, nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do j = 1,NJ(nbl)
		
		!!Testing for filtering
		!do i = 1,NI(nbl)
		!	pert = 0.d0
		!	if (i.ge.2.and.i.le.(NI(nbl)-1)) then 
		!	if (mod(i,2).eq.0) then 
		!	pert = -1.d0
		!	else
		!	pert = 1.d0
		!	endif
		!	endif
		!	PHI(i,j,k,nbl,var) = sin(xgrid(i,j,k,nbl)) + 0.1d0*pert
		!	if(i.eq.NI(nbl)) PHI(i,j,k,nbl,var)=0.d0
		!	print*, PHI(i,j,k,nbl,var)
		!enddo
		!!!!
		
	do i = 6,NI(nbl)-5						! For all interior points
		RHS(i) = 0.d0
		do coeff = 0,npts
		ipc = i+coeff
		imc = i-coeff	
		RHS(i) = RHS(i) + 0.5d0*fcoeff(coeff+1)*(PHI(ipc,j,k,nbl,var)+PHI(imc,j,k,nbl,var))
		enddo	
	enddo

	do i = 1,5						! For aboundary points from 1 to 5
		RHS(i) = 0.d0
		do coeff = 0,npts
		ipc = i+coeff
		imc = i-coeff	
		if (imc.le.0) imc = imc + NI(nbl) - 1
		RHS(i) = RHS(i) + 0.5d0*fcoeff(coeff+1)*(PHI(ipc,j,k,nbl,var)+PHI(imc,j,k,nbl,var))
		enddo	
	enddo
	
	
	do i = NI(nbl)-4,NI(nbl)						! For aboundary points from NI(nbl)-4 to NI(nbl)
		RHS(i) = 0.d0
		do coeff = 0,npts
		ipc = i+coeff
		imc = i-coeff	
		if (ipc.gt.NI(nbl)) ipc = ipc - NI(nbl) + 1
		RHS(i) = RHS(i) + 0.5d0*fcoeff(coeff+1)*(PHI(ipc,j,k,nbl,var)+PHI(imc,j,k,nbl,var))
		enddo	
	enddo
		
		nm1 = NI(nbl)-1
		call TDMAP(1,NI(nbl)-1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
		
			!Testing for filtering
			!do i = 1,NI(nbl)-1
			!print*, i, PHI(i,j,k,nbl,var), RHS(i)
			!enddo
			!pause
			!!!
			
		PHI(1:nm1,j,k,nbl,var) = RHS(1:nm1)
		PHI(NI(nbl),j,k,nbl,var) = RHS(1)
		
	enddo
	enddo
	enddo
	enddo
    END		  


    SUBROUTINE FILTERING_J(PHI,nvars)
    use declare_variables	 
    implicit none	  
	
    integer var,jpc,jmc,coeff,nvars,npts,nm1   
    real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
    real,dimension(NJmax)	:: RHS
	
	npts = fscheme/2
	
	do var=1,nvars
	do nbl =1,nblocks
	do k = 1,NK(nbl)
	do i = 1,NI(nbl)

	
	do j = 6,NI(nbl)-5						! For all interior points
		RHS(j) = 0.d0
		do coeff = 0,npts
		jpc = j+coeff
		jmc = j-coeff	
		RHS(j) = RHS(j) + 0.5d0*fcoeff(coeff+1)*(PHI(i,jpc,k,nbl,var)+PHI(i,jmc,k,nbl,var))
		enddo	
	enddo
	
	do j = 1,5								! For a boundary points from 1 to 5
		RHS(j) = 0.d0
		do coeff = 0,npts
		jpc = j+coeff
		jmc = j-coeff	
		if (jmc.le.0) jmc = jmc + NJ(nbl) - 1
		RHS(j) = RHS(j) + 0.5d0*fcoeff(coeff+1)*(PHI(i,jpc,k,nbl,var)+PHI(i,jmc,k,nbl,var))
		enddo	
	enddo
		
	do j = NJ(nbl)-4,NJ(nbl)				! For a boundary points from NI(nbl)-4 to NI(nbl)
		RHS(j) = 0.d0
		do coeff = 0,npts
		jpc = j+coeff
		jmc = j-coeff	
		if (jpc.gt.NJ(nbl)) jpc = jpc - NJ(nbl) + 1
		RHS(j) = RHS(j) + 0.5d0*fcoeff(coeff+1)*(PHI(i,jpc,k,nbl,var)+PHI(i,jmc,k,nbl,var))
		enddo	
	enddo
	
		nm1 = NJ(nbl)-1
		call TDMAP(1,NJ(nbl)-1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
		PHI(i,1:nm1,k,nbl,var) = RHS(1:nm1)
		PHI(i,NJ(nbl),k,nbl,var) = RHS(1)
		
	enddo
	enddo
	enddo
	enddo
	END	


	SUBROUTINE FILTERING_K(PHI,nvars)
	use declare_variables	 
	implicit none	  
	
	integer var,kpc,kmc,coeff,nvars,npts,nm1   
	real,dimension(NImax,NJmax,NKmax,nblocks,nvars) :: PHI
	real,dimension(NKmax)	:: RHS

	npts = fscheme/2
	
	do var=1,nvars
	do nbl =1,nblocks
	do j = 1,NJ(nbl)
	do i = 1,NI(nbl)
	
	do k = 6,NK(nbl)-5						! For all interior points
		RHS(k) = 0.d0
		do coeff = 0,npts
		kpc = k+coeff
		kmc = k-coeff	
		RHS(k) = RHS(k) + 0.5d0*fcoeff(coeff+1)*(PHI(i,j,kpc,nbl,var)+PHI(i,j,kmc,nbl,var))
		enddo	
	enddo
	 	 
	do k = 1,5								! For a boundary points from 1 to 5
		RHS(k) = 0.d0
		do coeff = 0,npts
		kpc = k+coeff
		kmc = k-coeff	
		if (kmc.le.0) kmc = kmc + NK(nbl) - 1
		RHS(k) = RHS(k) + 0.5d0*fcoeff(coeff+1)*(PHI(i,j,kpc,nbl,var)+PHI(i,j,kmc,nbl,var))
		enddo	
	enddo
		
	do k = NK(nbl)-4,Nk(nbl)				! For a boundary points from NI(nbl)-4 to NI(nbl)
		RHS(k) = 0.d0
		do coeff = 0,npts
		kpc = k+coeff
		kmc = k-coeff	
		if (kpc.gt.NK(nbl)) kpc = kpc - NK(nbl) + 1
		RHS(k) = RHS(k) + 0.5d0*fcoeff(coeff+1)*(PHI(i,j,kpc,nbl,var)+PHI(i,j,kmc,nbl,var))
		enddo	
	enddo


		nm1 = NK(nbl)-1
		call TDMAP(1,NK(nbl)-1,APF(1:nm1),ACF(1:nm1),AMF(1:nm1),RHS(1:nm1),nm1)
		PHI(i,j,1:nm1,nbl,var) = RHS(1:nm1)
		PHI(i,j,NK(nbl),nbl,var) = RHS(1)
		
	enddo
	enddo
	enddo
	enddo	
	END		

!************************** CYCLIC TDMA SOLVER *************************************	  

subroutine TDMAP(ji,jf,ap,ac,am,fi,NMAXL)

!ap - super diagonal
!ac - diagonal
!am - sub diagonal

implicit none

! -----------------------
!  Input/Output variables
! -----------------------
   integer:: ji, jf, NMAXL
   real:: ap(NMAXL), ac(NMAXL), am(NMAXL), fi(NMAXL)

! -------------------
!  Internal variables
! -------------------
   integer:: i, j, ja, jj
   real:: fnn, pp
   real:: qq(NMAXL), ss(NMAXL), fei(NMAXL)

   ja=ji+1
   jj=ji+jf

   qq(ji)=-ap(ji)/ac(ji)
   ss(ji)=-am(ji)/ac(ji)
   fnn=fi(jf)
   fi(ji)=fi(ji)/ac(ji)

!  forward elimination sweep
!----------------------------
   do j=ja,jf
      pp=1.0d0/(ac(j)+am(j)*qq(j-1))
          qq(j)=-ap(j)*pp
          ss(j)=-am(j)*ss(j-1)*pp
          fi(j)=(fi(j)-am(j)*fi(j-1))*pp
   enddo
   
!  backward pass
!----------------

   ss(jf)=1.0d0
   fei(jf)=0.0d0
 
   do i=ja,jf
      j=jj-i
          ss(j)=ss(j)+qq(j)*ss(j+1)
          fei(j)=fi(j)+qq(j)*fei(j+1)
   enddo
   
   fi(jf)=(fnn-ap(jf)*fei(ji)-am(jf)*fei(jf-1))/    &
   &      (ap(jf)*ss(ji)+am(jf)*ss(jf-1)+ac(jf)) 
   
!  backward substitution
!------------------------
   do i=ja,jf
      j=jj-i
          fi(j)=fi(jf)*ss(j)+fei(j)
   enddo

end subroutine TDMAP