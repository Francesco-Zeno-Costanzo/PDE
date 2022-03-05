	program traspnonllw
C=============================================================================
C Codice per risolvere l'equazione del trasporto: du/dt + df(u)/dx = 0
C Il codice utilizza il metodo di lax wendroff:
c=============================================================================
	implicit real*8 (a-h,p-z)
	f(u) = v*u**2
	
	real*8, dimension(:), allocatable :: sol_v, sol_n, sol_m
	
	open(1, file='input_lw.txt', status='old')	!file coi parametri
	open(2, file='tra_lw1.dat', status='unknown')	!file coi risultati
	
	read(1,*) N	!numeri punti sulle x
	read(1,*) i_T	!numeri di punti nel tempo
        read(1,*) v	!velocità di propagazione
        read(1,*) dt	!passo temporale
        read(1,*) dx	!passo spaziale 	 
        read(1,*) i_P	!ogni quanto salvare i dati su file
        
	allocate(sol_v(0:N), sol_n(0:N), sol_m(0:N))
     	
      	c = dt/dx
        
	q = 2*4*ATAN(1.) !2pigreco
        do i = 0, N
            x = i*dx
            sol_v(i) = sin(q*2*x)
      	enddo
      	
      	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        !soluzione intermedia con metodo di lax
	        do i = 0, N - 1
		    sol_m(i) = 0.5*(sol_v(i+1) + sol_v(i)) +
     &                         0.5*c*(f(sol_v(i+1))-f(sol_v(i)))
		enddo
		
		!soluzione finale con leap frog
		do i = 1, N - 1
                    sol_n(i) = sol_v(i)+c*(f(sol_m(i))-f(sol_m(i-1)))
		enddo
		
		!condizioni periodiche al bordo
		sol_n(0) = sol_n(N - 1)
		sol_n(N) = sol_n(1)
		    	
		!aggiorno la soluzione
		sol_v = sol_n
		
	    enddo
	    
	    !scrivo sul file
            do i = 0, N
        	write(2,*) sol_v(i)
	    enddo
	    
        enddo  
	deallocate(sol_v, sol_n, sol_m)
	stop
	
	end program traspnonllw
