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
      
      omega = 1.5	!over relaxation factor
      tol = 1E-6	!tollerance
      
      !matrix created according to the lexicographic ordering
      do i = 1, M
          a(i, i) = -4.d0
      enddo
      
      do i = 1, M-1
          a(i, i+1) = 1.d0
          a(i+1, i) = 1.d0
      enddo
      
      do i = 1, N-1
          a(i*N+1, i*N) = 0.d0
          a(i*N, i*N+1) = 0.d0
      enddo
      
      do i = 1, M-N
          a(i, i+N) = 1.d0
          a(i+N, i) = 1.d0
      enddo
      
      !boundary conditions
      !and source function

      rho(N/4, N/2) = -1.d0
      rho(3*N/4, N/2) = -1.d0

      rho(:,1) = 0.d0
      rho(:,N) = 0.d0
      rho(1,:) = 1.d0
      rho(N,:) = -1.d0
      
      !transform matrix into array
      do i = 0, N-1
          b(i*N + 1: (i+1)*N) = rho(i+1,:) 
      enddo
      !I write the charge distribution to the file
      do i = 1, M
          write(0,*) b(i)
      enddo
      
      call SOR(M, a, b, x, tol, omega)
      !I write the solution on the file
      do i = 1, M
          write(0,*) x(i)
      enddo

      end program laplace

C========================================================================
C subroutine successive over-relaxation
C========================================================================
   
      subroutine SOR(N, A, b, x, tol, omega)
C==================================================================
C     Subroutine for the solutions of linear system with SOR method
C     
C     Parameters
C     ----------
C     N : int
C         dimension of matrix
C     A : NxN matrix
C         matrix of the system
C     b : 1d array
C         known term
C     x : 1d array
C         initial guess, will be the solution ot the end
C     tol : flot
C         tollerance
C     omega : float
C         parameter of algorithm 
C==================================================================
      real*8, dimension(N, N) :: A
      real*8, dimension(N) :: b, x, prod
      real*8 :: tol, omega, sigma, res
      
      !initial residual
      res = 0    
      do i = 1, N
          prod(i) = sum(A(i,:)*x)
      enddo
      res = sqrt(sum((prod-b)**2.d0))
      
      !iter = 0
      do while (res>tol)
          do i = 1, N
              sigma = 0
              do j = 1, N
                  if (j/=i) then
                      sigma = sigma + A(i, j)*x(j)
                  endif 
              enddo
              x(i) = (1.d0-omega)*x(i) + (omega/A(i, i))*(b(i)-sigma)
          enddo
          
          do i = 1, N
              prod(i) = sum(A(i,:)*x)
          enddo
          res = sqrt(sum((prod-b)**2.d0))
          !iter = iter + 1
          !print*, iter, res
      enddo
      
      
      return
      end
