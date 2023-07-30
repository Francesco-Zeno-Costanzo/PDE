	program traspnonlinearlw
C=============================================================================
C Code to solve the transport equation: du/dt + df(u)/dx = 0. Where f can be
C an arbritary function of u, so also non linear case
C The code use the lax wendroff method.
c=============================================================================
	implicit real*8 (a-h,p-z)
	
	f(u) = v*u**2 ! non linear equations, wavefront break
	
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
     	
      c = dt/dx
        
	!Initial condition
	q = 2.d0 * 4.d0 * datan(1.d0) !2 pi 
      do i = 0, N
          x = i*dx
          sol_v(i) = dsin(q*0.2d0*x)
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
     &                     0.5*c*(f(sol_v(i+1))-f(sol_v(i)))
		  enddo
		
		  !solution with leap frog
		  do i = 1, N - 1
                  sol_n(i) = sol_v(i)+c*(f(sol_m(i))-f(sol_m(i-1)))
		  enddo
		
		  !poriodic boundary conditions
	        sol_n(0) = sol_n(N - 1)
		  sol_n(N) = sol_n(1)
		  !update solutions
	        sol_v = sol_n
		
	    enddo
	    
	    !scrivo sul file
          do i = 0, N
        	  write(2,*) sol_v(i)
	    enddo
	    
      enddo  
	
	deallocate(sol_v, sol_n, sol_m)
	stop
	
	end program 
