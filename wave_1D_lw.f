	program wavelaxwendroff
C=============================================================================
C Code to solve the wave equation: d^2u/dt^2 + vd^2u/dx^2 = 0
C The code uses lax wendroff's method indeed, if we call r = vdu/dx, s = du/dt
C we can rewrite the equation as a transport equation:
C dU/dt + dF(U)/dx = 0 where: U = (r, s) and F(U) = ((0, -v),(-v, 0))U
c=============================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: r_v, r_n, s_v, s_n
	real*8, dimension(:), allocatable :: sol_v, sol_n
	
	open(1, file='input_wave.txt', status='old')  !parameter file
	open(2, file='wave_lw.dat', status='unknown') !results file
	
	read(1,*) N	      !nuber of ponits on spatial  axis
	read(1,*) i_T	!nuber of ponits on temporal axis
      read(1,*) v	      !speed of the propagation
      read(1,*) dt	!temporal step
      read(1,*) dx	!spatial  step 
      read(1,*) i_P	!how often to save the solution
        
      allocate(r_v(0:N), r_n(0:N), s_v(0:N), s_n(0:N))
	allocate(sol_v(0:N), sol_n(0:N))
     	
      alpha = v*dt/dx
      print*, alpha
        
	!Initial condition
	q = 2.d0 * 4.d0 * datan(1.d0) !2 pi 
      do i = 0, N
          x = i*dx
          !u(x, t=0)
          sol_v(i) = dsin(q*0.3d0*x) !exp(-((x-0.5)/0.1)**2)
          !v*du/dx
          r_v(i) = v*q*0.3d0*dcos(q*0.3d0*x) !-v*exp(-((x-0.5)/0.1)**2)*200*(x-0.5)
          !du/dt
          s_v(i) = 0.d0
      enddo
      	
      !Save the initial condition
      do i = 0, N
          write(2,*) sol_v(i)
	enddo
      	
      do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        !solutions with lax
		  do i = 1, N - 1
                  r_n(i) = r_v(i) + 0.5*alpha*(s_v(i+1) - s_v(i-1))+
     &		    0.5*alpha**2*(r_v(i+1) - 2.0*r_v(i) + r_v(i-1))
     
     		      s_n(i) = s_v(i) + 0.5*alpha*(r_v(i+1) - r_v(i-1))+
     &		    0.5*alpha**2*(s_v(i+1) - 2.0*s_v(i) + s_v(i-1))
		  enddo
		
		  !poriodic boundary conditions
		  r_n(0) = r_n(N - 1)
		  r_n(N) = r_n(1)
		
		  s_n(0) = s_n(N - 1)
		  s_n(N) = s_n(1)
		
		  !solution with leap frog
		  do i = 0, N !- 1
                  sol_n(i) = sol_v(i) + 0.5*dt*(s_n(i) + s_v(i))
		  enddo
		
		  !update solutions 
		  r_v = r_n
	 	  s_v = s_n
		  sol_v = sol_n
		
	    enddo
	    
	    !scrivo sul file
          do i = 0, N
        	  write(2,*) sol_v(i)
	    enddo
	    
      enddo  
	
	deallocate(r_v, r_n, s_v, s_n, sol_v, sol_n)
	stop
	
	end program
