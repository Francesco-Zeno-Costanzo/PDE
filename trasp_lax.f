      program trasplax
C=========================================================================
C Code to solve the transport equation: du/dt + v*du/dx = 0. Where v is
C the speed of the propagation and is constant. We use the lax method:
C u^{n+1}_j = (u^n_{j+1} + u^n_{j-1})/2 - vdt/dx(u^n_{j+1} - u^n_{j-1}) 
C The parameters are written to a file so that they are readable by both
C the simulation code and the plot code
c=========================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: sol_v, sol_n
	
	open(1, file='input_trasp.txt', status='old')  !parameter file
	open(2, file='tra_lax.dat', status='unknown')  !results file
	
	read(1,*) N	      !nuber of ponits on spatial  axis
	read(1,*) i_T	!nuber of ponits on temporal axis
      read(1,*) v	      !speed of the propagation
      read(1,*) dt	!temporal step
      read(1,*) dx	!spatial  step 
      read(1,*) i_P	!how often to save the solution
        
	allocate(sol_v(0:N), sol_n(0:N))
     	
      alpha = v*dt/dx
      print*, alpha
      	
      !Initial condition
	q = 2.d0 * 4.d0 * datan(1.d0) !2 pi 
      do i = 0, N
          x = i*dx
          sol_v(i) = 15.d0*dsin(q*0.5d0*x)
      enddo
      	
      !Save the initial condition
      do i = 0, N
          write(2,*) sol_v(i)
	enddo
	
      !Time evolution
	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        do i = 1, N - 1
		      sol_n(i) = 0.5*(sol_v(i+1)*(1 - alpha))+                         
     &		           0.5*(sol_v(i-1)*(1 + alpha))
		  enddo
		    	
		  !poriodic boundary conditions
	        sol_n(0) = sol_n(N - 1)
		  sol_n(N) = sol_n(1)
		  !update solutions
	        sol_v = sol_n
		
	    enddo
	    !save the solution on file
          do i = 0, N
        	  write(2,*) sol_v(i)
	    enddo
	    
	enddo
	
	deallocate(sol_v, sol_n)
	stop
	
	end program trasplax
