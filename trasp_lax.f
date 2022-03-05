	program trasportoconlax
C=========================================================================
C Codice per risolvere l'equazione del trasporto: du/dt + v*du/dx = 0
C dove v è la velocità della propagazione ed è costante con metodo di lax:
C u^{n+1}_j = (u^n_{j+1} + u^n_{j-1})/2 - vdt/dx(u^n_{j+1} - u^n_{j-1}) 
c=========================================================================
	implicit real*8 (a-h,p-z)
	
	real*8, dimension(:), allocatable :: sol_v, sol_n
	
	open(1, file='input_l.txt', status='old')	!file coi parametri
	open(2, file='tra_lax.dat', status='unknown')	!file coi risultati
	
	read(1,*) N	!numeri punti sulle x
	read(1,*) i_T	!numeri di punti nel tempo
        read(1,*) v	!velocità di propagazione
        read(1,*) dt	!passo temporale
        read(1,*) dx	!passo spaziale 
        read(1,*) i_P	!ogni quanto salvare i dati su file
        
	allocate(sol_v(0:N), sol_n(0:N))
     	
      	alpha = v*dt/dx
      	print*, alpha
      	
      	!condizione iniziale
	q = 2*4*ATAN(1.) !2pigreco
        do i = 0, N
            x = i*dx
            sol_v(i) = 10*sin(q*1*x) + 15*sin(q*5*x)
      	enddo
      	
      	!evoluzione temporale con lax
	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	        do i = 1, N - 1
		    sol_n(i) = 0.5*(sol_v(i+1)*(1 - alpha))+                         
     &		               0.5*(sol_v(i-1)*(1 + alpha))
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
	
	deallocate(sol_v, sol_n)
	stop
	
	end program trasportoconlax
