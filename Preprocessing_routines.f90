!********************* READ INPUT ***********************************************************
      SUBROUTINE READ_INPUT()
      use declare_variables
      implicit none	  
      
	  open(unit=finput, file='input.dat', form='formatted')
	  
	  read(finput, *) Lx, Ly, Lz
	  read(finput, *) nblocks
	  
	  allocate(NI(nblocks))
	  allocate(NJ(nblocks))
	  allocate(NK(nblocks))
	  
	  read(finput,*) (NI(nbl), NJ(nbl), NK(nbl), nbl=1,nblocks) 
	  
	  read(finput,*) nconserv, nprims
	  
	  read(finput,*) Re, Mach, gamma, prandtl, T_ref
	  
	  read(finput,*) grid2D, dscheme
	  
	  read(finput,*) fscheme, alpha_f
	  
	  read(finput,*) euler
	  
	  read(finput,*) nsteps, time_step, animfreq
	  
	  !print*, 'print the following:', Lx, Ly, Lz, nblocks, (NI(nbl), NJ(nbl), NK(nbl), nbl=1,nblocks) 
	  
	  close(finput)
	  
      END
!********************************************************************************************

!********************* ALLOCATE_ROUTINE *****************************************************
      SUBROUTINE ALLOCATE_ROUTINE()
      use declare_variables
      implicit none	  
	  
	  NImax = maxval(NI)
	  NJmax = maxval(NJ)
	  NKmax = maxval(NK)
	  
	  print*, NImax, NJmax, NKmax
	  allocate(xgrid(NImax, NJmax, NKmax, nblocks))
	  allocate(ygrid(NImax, NJmax, NKmax, nblocks))
	  allocate(zgrid(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(xi(NImax, NJmax, NKmax, nblocks))
	  allocate(yi(NImax, NJmax, NKmax, nblocks))
	  allocate(zi(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(xj(NImax, NJmax, NKmax, nblocks))
	  allocate(yj(NImax, NJmax, NKmax, nblocks))
	  allocate(zj(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(xk(NImax, NJmax, NKmax, nblocks))
	  allocate(yk(NImax, NJmax, NKmax, nblocks))
	  allocate(zk(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(kx(NImax, NJmax, NKmax, nblocks))
	  allocate(ky(NImax, NJmax, NKmax, nblocks))
	  allocate(kz(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(ix(NImax, NJmax, NKmax, nblocks))
	  allocate(iy(NImax, NJmax, NKmax, nblocks))
	  allocate(iz(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(jx(NImax, NJmax, NKmax, nblocks))
	  allocate(jy(NImax, NJmax, NKmax, nblocks))
	  allocate(jz(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(Jac(NImax, NJmax, NKmax, nblocks))
	  
	  allocate(Qc(NImax, NJmax, NKmax, nblocks, nconserv))
	  allocate(Qcini(NImax, NJmax, NKmax, nblocks, nconserv))
	  allocate(Qp(NImax, NJmax, NKmax, nblocks, nprims))
	  
	  allocate(Qpi(NImax, NJmax, NKmax, nblocks, nprims))
	  allocate(Qpj(NImax, NJmax, NKmax, nblocks, nprims))
	  allocate(Qpk(NImax, NJmax, NKmax, nblocks, nprims))
      
	  ptsmax = max(NImax, NJmax, NKmax)
	  allocate(AMD(ptsmax))
	  allocate(ACD(ptsmax))
	  allocate(APD(ptsmax)) !TDMA diagonals
	  allocate(AMF(ptsmax))
	  allocate(ACF(ptsmax))
	  allocate(APF(ptsmax))
	  
	  allocate(FFlux(NImax, NJmax, NKmax, nblocks, nconserv))
	  allocate(GFlux(NImax, NJmax, NKmax, nblocks, nconserv))
	  allocate(HFlux(NImax, NJmax, NKmax, nblocks, nconserv))
	  
	  allocate(Net_flux(NImax, NJmax, NKmax, nblocks, nconserv))
	  allocate(fluxD(NImax, NJmax, NKmax, nblocks, nconserv))
	  
	  allocate(res(nconserv))
	  allocate(fcoeff(6))
	  
      END 
!********************************************************************************************

!********************* GENERATE_GRID *****************************************************
      SUBROUTINE GENERATE_GRID()
      use declare_variables 
      implicit none	  
	  
	  do nbl = 1,nblocks
	  do k = 1,NK(nbl)
	  do j = 1,NJ(nbl)
	  do i = 1,NI(nbl)
	  
	  xgrid(i,j,k,nbl) = (i-1.d0)/(NI(nbl)-1.d0)*Lx
	  ygrid(i,j,k,nbl) = (j-1.d0)/(NJ(nbl)-1.d0)*Ly
	  zgrid(i,j,k,nbl) = (k-1.d0)/(NK(nbl)-1.d0)*Lz
	  !print*, i,j,k, xgrid(i,j,k,nbl), ygrid(i,j,k,nbl), zgrid(i,j,k,nbl)
	  
	  enddo
	  enddo
	  enddo
	  enddo
	  
      END 
!********************************************************************************************

!********************* INITIALIZE_NON_DIMENSIONALIZE ****************************************
	SUBROUTINE INITIALIZE_NON_DIMENSIONALIZE()
	use declare_variables 
	implicit none	  
	
    real xl,yl,zl, Etotal
	real rhl,ul,vl,wl,pl
    integer nvars, n
	  
	DO nbl =1,nblocks
	DO k = 1,NK(nbl)
		DO j = 1,NJ(nbl)
			DO i = 1,NI(nbl)
			
			xl = xgrid(i,j,k,nbl)
			yl = ygrid(i,j,k,nbl)
			zl = zgrid(i,j,k,nbl)
			
			Qp(i,j,k,nbl,2) = sin(xl)*cos(yl)*cos(zl)    !!non dimensional u
			Qp(i,j,k,nbl,3) = -cos(xl)*sin(yl)*cos(zl)   !!non dimensional v
			Qp(i,j,k,nbl,4) = 0.d0                       !!non dimensional w
			Qp(i,j,k,nbl,5) = 1.d0/(mach**2*gamma) + (1.d0/16.d0)*((cos(2.d0*xl) + cos(2.d0*yl))*(cos(2.d0*zl) + 2.d0 )) !! non dim pressure
			Qp(i,j,k,nbl,6) = 1.d0                       !!non dimensional T
			Qp(i,j,k,nbl,1) = gamma*mach**2*Qp(i,j,k,nbl,5)/Qp(i,j,k,nbl,6) 				 ! check this a non dimensional variable will come here.
			
			
			rhl = Qp(i,j,k,nbl,1)
			ul  = Qp(i,j,k,nbl,2)
			vl  = Qp(i,j,k,nbl,3)
			wl  = Qp(i,j,k,nbl,4)
			pl  = Qp(i,j,k,nbl,5)
			
			Qc(i,j,k,nbl,1) = rhl
			Qc(i,j,k,nbl,2) = rhl*ul
			Qc(i,j,k,nbl,3) = rhl*vl
			Qc(i,j,k,nbl,4) = rhl*wl
			Etotal = pl/(rhl*(gamma-1.d0)) + 0.5d0*(ul**2 + vl**2 + wl**2)
			Qc(i,j,k,nbl,5) = rhl*Etotal 													! this expression will also change here.

			ENDDO
		ENDDO	
	ENDDO
	ENDDO
	END 
!********************************************************************************************

!********************* DISCRETIZATION_VALS **************************************************
      SUBROUTINE DISCRETIZATION_FILTER_RK_VALS()
      use declare_variables 
      implicit none
	  
	  !Assign Discretization Coefficients -- Interior Points
	  if (dscheme.eq.1) then
	  alpha = 0.d0
	  adisc = 1.d0
	  bdisc = 0.d0
	  elseif (dscheme.eq.2) then
	  alpha = 0.d0
	  adisc = 4.d0/3.d0
	  bdisc = -1.d0/3.d0
	  elseif (dscheme.eq.3) then
	  alpha = 1.d0/4.d0
	  adisc = 3.d0/2.d0
	  bdisc = 0.d0
	  
	  AMD(1:ptsmax)= alpha
	  ACD(1:ptsmax)= 1.d0
	  APD(1:ptsmax)= alpha

	  elseif (dscheme.eq.4) then
	  alpha = 1.d0/3.d0
	  adisc = 14.d0/9.d0
	  bdisc = 1.d0/9.d0
	  
	  AMD(1:ptsmax)= alpha
	  ACD(1:ptsmax)= 1.d0
	  APD(1:ptsmax)= alpha
	  endif
	  
	  !**************** Assign Filter coefficient- For Interior Points ***************
	  
	  if(fscheme.eq.2) then
	  
      fcoeff(1) = 0.5d0 + alpha_f
      fcoeff(2) = 0.5d0 + alpha_f
      fcoeff(3) = 0.d0
      fcoeff(4) = 0.d0
      fcoeff(5) = 0.d0
      fcoeff(6) = 0.d0
	  
      elseif(fscheme.eq.4) then
      fcoeff(1) = 5.d0/8.d0 + 3.d0*alpha_f/4.d0
      fcoeff(2) = 0.5d0 + alpha_f
      fcoeff(3) = -1.d0/8.d0 + 1.d0*alpha_f/4.d0
      fcoeff(4) = 0.d0
      fcoeff(5) = 0.d0
      fcoeff(6) = 0.d0
	  
      elseif(fscheme.eq.6) then
      fcoeff(1) = 11.d0/16.d0 + 5.d0*alpha_f/8.d0
      fcoeff(2) = 15.d0/32.d0 + 17.d0*alpha_f/16.d0
      fcoeff(3) = -3.d0/16.d0 + 3.d0*alpha_f/8.d0
      fcoeff(4) = 1.d0/32.d0 - 1.d0*alpha_f/16.d0
      fcoeff(5) = 0.d0
      fcoeff(6) = 0.d0
	  
      elseif(fscheme.eq.8) then
      fcoeff(1) = 93.d0/128.d0 + 70.d0*alpha_f/128.d0
      fcoeff(2) = 7.d0/16.d0 + 18.d0*alpha_f/16.d0
      fcoeff(3) = -7.d0/32.d0 + 14.d0*alpha_f/32.d0
      fcoeff(4) = 1.d0/16.d0 - 1.d0*alpha_f/8.d0
      fcoeff(5) = -1.d0/128.d0 + 1.d0*alpha_f/64.d0
      fcoeff(6) = 0.d0
	  
      elseif(fscheme.eq.10) then
      fcoeff(1) = 193.d0/256.d0 + 126.d0*alpha_f/256.d0
      fcoeff(2) = 105.d0/256.d0 + 302.d0*alpha_f/256.d0
      fcoeff(3) = -15.d0/64.d0 + 30.d0*alpha_f/64.d0
      fcoeff(4) = 45.d0/512.d0 - 90.d0*alpha_f/512.d0
      fcoeff(5) = -5.d0/256.d0 + 10.d0*alpha_f/256.d0
      fcoeff(6) = 1.d0/512.d0 - 2.d0*alpha_f/512.d0
      endif
	  
	  AMF(1:ptsmax)= alpha_f
	  ACF(1:ptsmax)= 1.d0
	  APF(1:ptsmax)= alpha_f
	  
	  !**************** Assign RK coefficient ***************
	  rk_steps = 4   					! later have to update to 4.
	  allocate(fac_qini(rk_steps))
	  allocate(fac_RK(rk_steps))
	  
	  fac_qini(1) = 0.d0
	  fac_qini(2) = 0.5d0
	  fac_qini(3) = 0.5d0
	  fac_qini(4) = 1.d0
	  
	  fac_RK(1)   = 1.d0/6.d0
	  fac_RK(2)   = 2.d0/6.d0
	  fac_RK(3)   = 2.d0/6.d0
	  fac_RK(4)   = 1.d0/6.d0
	  
      END 
!********************************************************************************************

!********************* METRICS **************************************************************
    SUBROUTINE METRICS()
    use declare_variables 
    implicit none	  
	
    real xil,yil,zil,xjl,yjl,zjl,xkl,ykl,zkl
    real vol
	
	if(dscheme.eq.1.or.dscheme.eq.2) then
	call DISCRETIZATION_I_EXP_GRID(xgrid,xi)
	call DISCRETIZATION_I_EXP_GRID(ygrid,yi)
	call DISCRETIZATION_I_EXP_GRID(zgrid,zi)
	
	call DISCRETIZATION_J_EXP_GRID(xgrid,xj)
	call DISCRETIZATION_J_EXP_GRID(ygrid,yj)
	call DISCRETIZATION_J_EXP_GRID(zgrid,zj)
	
	if(grid2d.ne.1) then
		call DISCRETIZATION_K_EXP_GRID(xgrid,xk)
		call DISCRETIZATION_K_EXP_GRID(ygrid,yk)
		call DISCRETIZATION_K_EXP_GRID(zgrid,zk)
	else
		xk = 0.d0
		yk = 0.d0
		zk = 1.d0
	endif
	
	
	elseif(dscheme.eq.3.or.dscheme.eq.4) then
	
	call DISCRETIZATION_I_COMP_GRID(xgrid,xi)
	call DISCRETIZATION_I_COMP_GRID(ygrid,yi)
	call DISCRETIZATION_I_COMP_GRID(zgrid,zi)
						  
	call DISCRETIZATION_J_COMP_GRID(xgrid,xj)
	call DISCRETIZATION_J_COMP_GRID(ygrid,yj)
	call DISCRETIZATION_J_COMP_GRID(zgrid,zj)
		
	if(grid2d.ne.1) then
		call DISCRETIZATION_K_COMP_GRID(xgrid,xk)
		call DISCRETIZATION_K_COMP_GRID(ygrid,yk)
		call DISCRETIZATION_K_COMP_GRID(zgrid,zk)
	else
		xk = 0.d0
		yk = 0.d0
		zk = 1.d0
	endif
	
	endif
	  !print*, 'abcd3'
	  
	  print*, 'discretization of grid done'
	  
	  DO nbl = 1,nblocks
		DO k = 1,NK(nbl)
			DO j = 1,NJ(nbl)
				DO i = 1,NI(nbl)
					xil = xi(i, j, k, nbl)
					yil = yi(i, j, k, nbl)
					zil = zi(i, j, k, nbl)
					
					xjl = xj(i, j, k, nbl)
					yjl = yj(i, j, k, nbl)
					zjl = zj(i, j, k, nbl)
					
					xkl = xk(i, j, k, nbl)
					ykl = yk(i, j, k, nbl)
					zkl = zk(i, j, k, nbl)
					
					vol = xil * (yjl * zkl - ykl * zjl) - xjl * (yil * zkl - ykl * zil) + xkl * (yil * zjl - yjl * zil)
					Jac(i, j, k, nbl) = 1.d0 / vol
					
					ix(i, j, k, nbl) = (yjl * zkl - zjl * ykl) / vol
					iy(i, j, k, nbl) = (zjl * xkl - xjl * zkl) / vol
					iz(i, j, k, nbl) = (xjl * ykl - yjl * xkl) / vol
					
					jx(i, j, k, nbl) = (ykl * zil - zkl * yil) / vol 
					jy(i, j, k, nbl) = (zkl * xil - xkl * zil) / vol
					jz(i, j, k, nbl) = (xkl * yil - ykl * xil) / vol
					
					kx(i, j, k, nbl) = (yil * zjl - zil * yjl) / vol 
					ky(i, j, k, nbl) = (zil * xjl - xil * zjl) / vol
					kz(i, j, k, nbl) = (xil * yjl - yil * xjl) / vol
					
					!print*, '(i,j,k,nbl), iz, jz, kz',i,j,k,nbl,iz(i, j, k, nbl),jz(i, j, k, nbl),kz(i, j, k, nbl)
					
				ENDDO
			ENDDO
		ENDDO
	  ENDDO
    END 
!********************************************************************************************