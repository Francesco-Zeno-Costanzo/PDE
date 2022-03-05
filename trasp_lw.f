	program trasportoconlaxwendroff
C=============================================================================
C Codice per risolvere l'equazione del trasporto: du/dt + vdu/dx = 0
C Il codice utilizza il metodo di lax wendroff: (alpha = vdt/dx)
C un primo mezzo passo con lax e poi leap frog tra lax e sol precedente 
C u^{n+1/2}_{j+1/2} = (u^n_{j+1} + u^n_{j-1})/2 - (alpha/2)(u^n_{j+1} - u^n_j)
C u^{n+1}_j = u^n_j - alpha(u^{n+1/2}_{j+1/2} - u^{n+1/2}_{j+1/2})
c=============================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: sol_v, sol_n, sol_m
	
	open(1, file='input_lw.txt', status='old')	!file coi parametri
	open(2, file='tra_lw.dat', status='unknown')	!file coi risultati
	
	read(1,*) N	!numeri punti sulle x
	read(1,*) i_T	!numeri di punti nel tempo
        read(1,*) v	!velocità di propagazione
        read(1,*) dt	!passo temporale
        read(1,*) dx	!passo spaziale 	 
        read(1,*) i_P	!ogni quanto salvare i dati su file
        
	allocate(sol_v(0:N), sol_n(0:N), sol_m(0:N))
     	
      	alpha = v*dt/dx
      	print*, alpha
        
	q = 2*4*ATAN(1.) !2pigreco
        do i = 0, N
            x = i*dx
            sol_v(i) = sin(q*2*x) + sin(q*10*x)
      	enddo
      	
      	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        !soluzione intermedia con metodo di lax
	        do i = 0, N - 1
		    sol_m(i) = 0.5*(sol_v(i+1) + sol_v(i)) +
     &                         0.5*alpha*(sol_v(i+1) - sol_v(i))
		enddo
		
		!soluzione finale con leap frog
		do i = 1, N - 1
                    sol_n(i) = sol_v(i)+alpha*(sol_m(i) - sol_m(i-1))
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
	
	end program trasportoconlaxwendroff
