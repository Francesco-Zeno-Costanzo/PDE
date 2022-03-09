	program calore
	implicit real*8(a-h,o-z)
	
	real*8, dimension(:), allocatable :: sol_v, sol_n, a, b, c
	
	open(1, file='input_c.txt', status='old')	!file coi parametri
	open(2, file='calorekn.dat', status='unknown')	!file coi risultati
	
	read(1,*) N	!numeri punti sulle x
	read(1,*) i_T	!numeri di punti nel tempo
        read(1,*) D	!parammetro diffusione
        read(1,*) dt	!passo temporale
        read(1,*) dx	!passo spaziale 	 
        read(1,*) i_P	!ogni quanto salvare i dati su file
        
	allocate(sol_v(N), sol_n(N), b(N))
	allocate(a(N), c(N)) !solo due array perche due diagonali sono uguali
	
	r = D*dt/(2*dx*dx)
	print*,r
	
	!condizione iniziale
	l = dx*N
	do i = 1, N
            x = i*dx
            sol_v(i) = 50*exp(-((x - l/2)/0.2)**2)
      	enddo
      	
      	!scrivo sul file la condizione iniziale
         do i = 1, N
       	     write(2,*) sol_v(i)
	 enddo
	    
	do i_time = 1, i_T
	    do i_tpass = 1, i_P
	    	!creo la matrice che mi da l'evoluzione del sistema
	    	!devo fare questa operazione ogni volta perchè la
	    	!subroutine modifica gli elementi di matrice 
            	do i = 1, N
                    a(i) = 1.0 + 2.0*r
               	    c(i) = -r
            	enddo
            	!parte della soluzione esplicita
            	do i = 2, N-1
                    b(i) = r*sol_v(i+1)+(1-2*r)*sol_v(i)+r*sol_v(i-1)    
                enddo
                b(1) = r*sol_v(1+1)+(1-2*r)*sol_v(1)
                b(N) = (1-2*r)*sol_v(N)+r*sol_v(N-1)
                
            	!chiamo la soubroutine che risolve il sistema
            	!per la parte implicita
            	call solve(N, a, c, c, b, sol_n)
            	!aggiorno la soluzione
            	sol_v = sol_n
	    enddo
	    
            !scrivo sul file
            do i = 1, N
        	write(2,*) sol_v(i)
	    enddo
         
      	enddo
      	
      	deallocate(sol_v, sol_n, a, c)
      	end program calore
	
C=============================================================================
C soubroutine per la risoluzione di sistemi con matrici tridiagonali
C=============================================================================
	
	subroutine solve(N, d_ii, d_sup, d_inf, b, x)
	implicit real*8 (a-h,o-z)
	real*8, dimension(N) :: d_ii, d_sup, d_inf, b, x
	
	!d_ii è la diagonale della matrice
	!d_sup è la diagonale superiore
	!d_inf è la diagonale inferiore
	
	!rendo la matrice una triangolare superiore
	!d_inf viene quindi annulato
	do i = 2, N
	    if(d_ii(i-1) == 0.0) then
                print *,"nemmeno dio può dividere per zero"
                return
            endif
            a = d_inf(i-1)/d_ii(i-1)
            d_ii(i) = d_ii(i) - a*d_sup(i-1)
            b(i) = b(i) - a*b(i-1)
	enddo

	!a questo punto è facile trovare x_n dato che
	!nell'ultima riga della matrice è rimasto solo
	!l'elemento sulla diagonale ad essere non nullo
	if(d_ii(N) == 0.0) then
       	    print *,"nemmeno dio può dividere per zero"
       	    return
       	endif
      	x(N) = b(N)/d_ii(N)
      	
      	!e adesso andando a ritroso posso trovare tutti gli altri x_i
      	do i = N-1, 1, -1
            b(i) = b(i) - d_sup(i)*x(i+1)
            if(d_ii(i) == 0.0) then
                print *,"nemmeno dio può dividere per zero"
                return
            endif
            x(i) = b(i)/d_ii(i)
      	enddo

	return
	end
