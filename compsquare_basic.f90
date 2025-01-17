!***************************************************************************************
!****************** WRITTEN BY DR. NAGABHUSHANA RAO VADLAMANI **************************
!***************BASED ON THE HIGH ORDER COMPSQUARE SOLVER DEVELOPED BY DR. NRV**********
!***************DISTRIBUTED AS A PART OF AS6041 COURSE ON ADVANCED CFD *****************
!***************************************************************************************

      PROGRAM MISSION_AS6041
	  	  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
!@@@@@@@@@@@@@ PRE PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	

!****** Declare Arrays, Integers, Real numbers *********************************	
      use declare_variables
      implicit none
	  
      integer step,var
	  real tke1, enstpt1 
	                    
      print*, 'Declared Variables'

!******** READ INPUT ***********************************************************
      print*, 'Reading input file....'
      call READ_INPUT()
      print*, 'Reading input file Done, Discretization scheme:', dscheme		  
!****** INITIALIZE ARRAYS ******************************************************
      call ALLOCATE_ROUTINE()
      print*, 'Allocated Integers, Real numbers and Arrays Done'	
!****** GENERATE GRID ******************************************************
      call GENERATE_GRID()
      print*, 'Grid Generation Done'	  
!******** Initialize and Non-dimensionalize arrays *****************************
      print*, 'Initializing and non-dimensionalizing arrays....'
      call INITIALIZE_NON_DIMENSIONALIZE()
      print*, 'Initializing and non-dimensionalizing arrays done'
	  call output(1)
	  !pause
!******** Read Discretization coeffs for flow !*********************************
      print*, 'Obtaining Discretization Filter RK coefficients for flow....'
      call DISCRETIZATION_FILTER_RK_VALS()
      print*, 'Obtaining Discretization Filter RK coefficients for flow Done'

	  !call DISCRETIZATION_I_EXP(Qp,Qpi,nprims)
	  !call DISCRETIZATION_I_EXP_GRID(xgrid,xi)
	  !call DISCRETIZATION_I_EXP_GRID(ygrid,yi)
	  !call DISCRETIZATION_I_EXP_GRID(zgrid,zi)
!******** Compute Metric terms *************************************************
      print*, 'Computing Metrics....'
      call METRICS()
      print*, 'Computing Metrics Done'
	  !CALL DISCRETIZATION_I_COMP(Qp,Qpi,nprims)
	  
!****************** Monitor residuals ******************************************
      if(restart.eq.1) OPEN(fresidual,file='Monitor.out',form='formatted',access='append')	  	  
      if(restart.eq.0) OPEN(fresidual,file='Monitor.out',form='formatted')	  	  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	  time = 0.d0           !Initialize Time
	  
      DO iter = 1,nsteps   !EXPLICIT TIME STEPPING
      
	    Qcini = Qc         !Conservative variables before entering RK time integration              -- initial
		Qcnew = Qc         !Conservative variable for the new time step, initialized to Old value   -- new
		
		res = 0.d0 		   ! Initializeresidual at every time step
		tke = 0.d0
		enstpt = 0.d0
		
       Do step = 1,rk_steps
       
        !Estimate Inviscid and Viscous fluxes + Time integration
        call UNSTEADY(step)
        
        !Filter the solution in all three directions for stability        
        if(step.eq.4) then
         call FILTERING_I(Qc,nconserv)		
         call FILTERING_J(Qc,nconserv)		
         if(grid2d.ne.1) call FILTERING_K(Qc,nconserv)				
        endif
        
        !Estimate Primitive variables from Conservative variables
        call SET_PRIMITIVES()
        
       Enddo
		
        !Any additional post processing specific to test case
        !if(taylor.eq.1) call volume_integral()
		
		tke = tke/real((NI(1)-1.d0)*(NJ(1)-1.d0)*(NK(1)-1.d0))
		enstpt = enstpt/real((NI(1)-1.d0)*(NJ(1)-1.d0)*(NK(1)-1.d0)) 
		
        print*, time, tke, enstpt  !divide tke,enstpt by (NI*NJ*NK) 
        write(fresidual,*) time, tke, enstpt
		!res = 0.d0 !Reinitializing

        !Animation		
        !if(mod(iter,animfreq).eq.0) call OUTPUT(1)		
        
        !Time increment
		time = time + time_step		
		
      ENDDO	

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ END OF MAIN TIME LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      CLOSE(fresidual)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ POST PROCESSING STEPS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
	!****** Write the output ***************
101	  print*, 'Writing Output....'
      call SET_PRIMITIVES()
      !call OUTPUT(0)
      print*, 'Writing output done'
 
	!****** Deallocate Arrays **************
      call DEALLOCATE_ROUTINE()
      print*, 'Deallocated Arrays'	  

      END PROGRAM

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!@@@@@@@@@@@@@@@@@@@@@@@ SUBROUTINES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!********* Subroutines in Pre-processing, Solver, Postprocessing.f90 files ****************
