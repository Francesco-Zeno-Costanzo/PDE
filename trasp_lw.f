	program laxwendroff
C=============================================================================
C Code to solve the transport equation: du/dt + v*du/dx = 0. Where v is
C the speed of the propagation and is constant. Using the lax wendroff method:
C a first half step with lax and then leap frog between lax and previous solution
C
C u^{n+1/2}_{j+1/2} = (u^n_{j+1} + u^n_{j-1})/2 - (alpha/2)(u^n_{j+1} - u^n_j)
C u^{n+1}_j = u^n_j - alpha(u^{n+1/2}_{j+1/2} - u^{n+1/2}_{j+1/2})
C
C Where alpha = v*dt/dx
C The parameters are written to a file so that they are readable by both
C the simulation code and the plot code
c=============================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: sol_v, sol_n, sol_m
	
	open(1, file='input_trasp.txt', status='old')  !parameter file
	open(2, file='tra_lxw.dat', status='unknown')  !results file
	
	read(1,*) N	      !nuber of ponits on spatial  axis
	read(1,*) i_T	!nuber of ponits on temporal axis
      read(1,*) v	      !speed of the propagation
      read(1,*) dt	!temporal step
      read(1,*) dx	!spatial  step 
      read(1,*) i_P	!how often to save the solution
        
	allocate(sol_v(0:N), sol_n(0:N), sol_m(0:N))
     	
      alpha = v*dt/dx
      print*, alpha
        
	!Initial condition
	q = 2.d0 * 4.d0 * datan(1.d0) !2 pi 
      do i = 0, N
          x = i*dx
          sol_v(i) = 15.d0*dsin(q*2.d0*x)*dexp(-(x-5)**2)
      enddo
      
      !Save the initial condition
      do i = 0, N
          write(2,*) sol_v(i)
	enddo
      	
      do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        !solution with lax
	        do i = 0, N - 1
		      sol_m(i) = 0.5*(sol_v(i+1) + sol_v(i)) +
     &                 0.5*alpha*(sol_v(i+1) - sol_v(i))
		  enddo
		
		  !solution with leap frog
		  do i = 1, N - 1
                  sol_n(i) = sol_v(i)+alpha*(sol_m(i) - sol_m(i-1))
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
	
	deallocate(sol_v, sol_n, sol_m)
	stop
	
	end program
