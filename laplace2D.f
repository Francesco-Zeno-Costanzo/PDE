      program laplace
      implicit real*8(a-h,o-z)
      real*8, dimension(:,:), allocatable :: a, rho
      real*8, dimension(:), allocatable :: x, b
      
      open(unit=0, file='laplace.dat', status='unknown')
      
      N = 50
      M = N**2 
      write(0, *) N
      write(0, *) M
      allocate(a(M,M), x(M), b(M), rho(N,N))
      
      omega = 1.5	!fattore di rilassamento
      tol = 1E-6	!tolleranza
      
      !matrice creata secondo l'ordinamento lessicografico
      do i = 1, M
          a(i, i) = -4
      enddo
      
      do i = 1, M-1
          a(i, i+1) = 1
          a(i+1, i) = 1
      enddo
      
      do i = 1, N-1
          a(i*N+1, i*N) = 0
          a(i*N, i*N+1) = 0
      enddo
      
      do i = 1, M-N
          a(i, i+N) = 1
          a(i+N, i) = 1
      enddo
      
      !condizioni al contorno
      !e funzione  sorgente

      rho(N/4, N/2) = -1
      rho(3*N/4, N/2) = -1

      rho(:,1) = 0
      rho(:,N) = 0
      rho(1,:) = 1
      rho(N,:) = -1
      
      !trasformo matrice in array
      do i = 0, N-1
          b(i*N + 1: (i+1)*N) = rho(i+1,:) 
      enddo
      !scrivo la distribuzione di carica sul file
      do i = 1, M
          write(0,*) b(i)
      enddo
      
      call SOR(M, a, b, x, tol, omega)
      !scrivo la soluzione sul file     
      do i = 1, M
          write(0,*) x(i)
      enddo

      end program laplace

C========================================================================
C soubroutine algoritmo Test successive over-relaxation
C========================================================================
   
      subroutine SOR(N, A, b, x, tol, omega)
      real*8, dimension(N, N) :: A
      real*8, dimension(N) :: b, x, prod
      real*8 :: tol, omega, sigma, res
      
      !calcolo residuo iniziale
      res = 0    
      do i = 1, N
          prod(i) = sum(A(i,:)*x)
      enddo
      res = sqrt(sum((prod-b)**2))
      
      !iter = 0
      do while (res>tol)
          do i = 1, N
              sigma = 0
              do j = 1, N
                  if (j/=i) then
                      sigma = sigma + A(i, j)*x(j)
                  endif 
              enddo
              x(i) = (1-omega)*x(i) + (omega/A(i, i))*(b(i)-sigma)
          enddo
          
          do i = 1, N
              prod(i) = sum(A(i,:)*x)
          enddo
          res = sqrt(sum((prod-b)**2))
          !iter = iter + 1
          !print*, iter, res
      enddo
      
      
      return
      end
