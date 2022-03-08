	program ondeconlaxwendroff
C=============================================================================
C Codice per risolvere l'equazione delle onde: d^2u/dt^2 + vd^2u/dx^2 = 0
C Il codice utilizza il metodo di lax wendroff, infatti eseguendo se chiamiamo
C r = vdu/dx e s = du/dt possiamo riscrivere l'equazione come un' equazione di
C trasporto: dU/dt + dF(U)/dx = 0 dove: U = (r, s) e F(U) = ((0, -v),(-v, 0))U
c=============================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: r_v, r_n, s_v, s_n
	real*8, dimension(:), allocatable :: sol_v, sol_n
	
	open(1, file='input_onde.txt', status='old')	!file coi parametri
	open(2, file='onde_lw.dat', status='unknown')	!file coi risultati
	
	read(1,*) N	!numeri punti sulle x
	read(1,*) i_T	!numeri di punti nel tempo
        read(1,*) v	!velocità di propagazione
        read(1,*) dt	!passo temporale
        read(1,*) dx	!passo spaziale 	 
        read(1,*) i_P	!ogni quanto salvare i dati su file
        
        allocate(r_v(0:N), r_n(0:N), s_v(0:N), s_n(0:N))
	allocate(sol_v(0:N), sol_n(0:N))
     	
      	alpha = v*dt/dx
      	print*, alpha
        
	q = 2*4*ATAN(1.) !2pigreco
        do i = 0, N
            x = i*dx
            !u(x, t=0)
            sol_v(i) = sin(q*x) !exp(-((x-0.5)/0.1)**2)
            !v*du/dx
            r_v(i) = v*q*cos(q*x) !-v*exp(-((x-0.5)/0.1)**2)*200*(x-0.5)
            !du/dt
            s_v(i) = 0.0
      	enddo
      	
      	!scrivo sul file la condizioone iniziale
        do i = 0, N
            write(2,*) sol_v(i)
	enddo
      	
      	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        !soluzioni con lax wendroff
		do i = 1, N - 1
                    r_n(i) = r_v(i) + 0.5*alpha*(s_v(i+1) - s_v(i-1))+
     &		    0.5*alpha**2*(r_v(i+1) - 2.0*r_v(i) + r_v(i-1))
     
     		    s_n(i) = s_v(i) + 0.5*alpha*(r_v(i+1) - r_v(i-1))+
     &		    0.5*alpha**2*(s_v(i+1) - 2.0*s_v(i) + s_v(i-1))
		enddo
		
		!condizioni al bordo periodico
		r_n(0) = r_n(N - 1)
		r_n(N) = r_n(1)
		
		s_n(0) = s_n(N - 1)
		s_n(N) = s_n(1)
		    	
		do i = 0, N !- 1
                    sol_n(i) = sol_v(i) + 0.5*dt*(s_n(i) + s_v(i))
		enddo
		
		!aggiorno le soluzioni
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
	
	end program ondeconlaxwendroff
