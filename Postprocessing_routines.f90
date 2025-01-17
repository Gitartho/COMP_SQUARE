!********************* OUTPUT GRID & FLOW *****************************************************
      SUBROUTINE OUTPUT(flag)	  
      use declare_variables
      implicit none	  
	  
      integer n,nvars,flag
	  character(30)::filename
	  
	  print*, 'Writing Output....'
	  if(flag.eq.0) then
	  filename= 'grid.xyz'
		
	  !writing in plot3d file format to visualize in paraview, follow the link
	  open(unit=fgrid, form='unformatted', file=filename)
	  write(fgrid) nblocks
	  write(fgrid) (NI(nbl), NJ(nbl), NK(nbl), nbl=1,nblocks)
	  
	  do nbl = 1,nblocks
	  write(fgrid) (((xgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)) &
		,		   (((ygrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)) &
		,		   (((zgrid(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)) 
		
		!! all of these is being written into a single line, in formatted file time there might be newlines just for representation
		!! by formatting we make is readable to us, but also increase the file size and comutation time
		!! Writing order in the plot3d file (blk1) x coord -> y coord -> z coord (newline) (next block)
      enddo
	  close(fgrid)
		
	  elseif(flag.eq.1) then 
	  write(filename,'(a,i5.5,a)') 'flow',iter,'.xyz'
	  
	  open(unit=fflow, form='unformatted', file=filename)
	  
	  nvars = nprims + nconserv + 3*3 + 1
	  write(fflow) nblocks
	  write(fflow) (NI(nbl), NJ(nbl), NK(nbl), nvars, nbl=1,nblocks)
	  
	  do nbl = 1,nblocks
	  write(fflow)     (((( Qc(i,j,k,nbl,n), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), n = 1, nconserv)&
		,			   (((( Qp(i,j,k,nbl,n), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl)), n = 1, nprims)  &
		,			   ((( ix(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( iy(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))   				   &
		,			   ((( iz(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( jx(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( jy(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))   				   &
		,			   ((( jz(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( kx(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( ky(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))   				   &
		,			   ((( kz(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))					   &
		,			   ((( Jac(i,j,k,nbl), i=1,NI(nbl)), j=1,NJ(nbl)), k=1,NK(nbl))
	  enddo
	  close(fflow)
	  ENDIF
	  
	  END
	  
      SUBROUTINE DEALLOCATE_ROUTINE()
      use declare_variables
      implicit none	  
	  
	  deallocate(NI,NJ,NK,xgrid,ygrid,zgrid)
	  
      END	

      SUBROUTINE volume_integral()
      use declare_variables
      implicit none	  
      real,dimension(NJmax,NKmax,nblocks,2) :: plane
      real,dimension(NKmax,nblocks,2) :: line	  
      real muavg, dcell

	  
      END	  