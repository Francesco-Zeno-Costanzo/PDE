      program calore
C==================================================================================================
C     Code to solve the heat equation: du/dt - D*d^2u/dx^2 = 0. Where D is
C     the diffusion parameter and is constant. We use crank nicolson method:
C
C     (u^{n+1}_j - u^n_j)/dt = D/(2dx^2)(
C                              (u^{n+1}_{j+1} - 2u^{n+1}_j + u^{n+1}_{j-1}) +
C                            + (u^n_{j+1}     - 2u^n_j     + u^n_{j-1}) )
C
C     if r = Ddt/(2dx^2) we have:
C
C     -ru^{n+1}_{j-1} + (1+2r)u^{n+1}_j - ru^{n+1}_{j+1} = ru^n_{j+1} + (1-2r)u^n_j + ru^n_{j-1}
C
C     That is in the form of A*u = b; with A that is a tridiagonal matrix
C
C     The parameters are written to a file so that they are readable by both
C     the simulation code and the plot code
C=================================================================================================
      implicit real*8(a-h,o-z)

      real*8, dimension(:), allocatable :: sol_v, sol_n, a, b, c

      open(1, file='input_ht.txt', status='old')  !parameter file
      open(2, file='heat.dat', status='unknown')  !results file

      read(1,*) N       ! Number of ponits on spatial  axis
      read(1,*) i_T     ! Number of ponits on temporal axis
      read(1,*) D       ! Diffusion parameter
      read(1,*) dt      ! Temporal step
      read(1,*) dx      ! Spatial  step
      read(1,*) i_P     ! How often to save the solution
        
      allocate(sol_v(N), sol_n(N), b(N))
      allocate(a(N), c(N)) !only two array beacuse two diagonals are equals

      r = D*dt/(2.d0*dx**2.d0)
      print*,r

      !initial condition
      l = dx*N
      do i = 1, N
          x = i*dx
          sol_v(i) = 500.d0*dexp(-((x - l/2.d0)/0.1d0)**2.d0)
      enddo

      !Save the initial condition
      do i = 1, N
          write(2,*) sol_v(i)
      enddo

      ! Initalize matrix A
      do i = 1, N
          a(i) = 1.d0 + 2.d0*r
          c(i) = -r
      enddo

      do i_time = 1, i_T
          do i_tpass = 1, i_P
              !Explicit solutions
              do i = 2, N-1
                  b(i) = r*sol_v(i+1)+(1-2*r)*sol_v(i)+r*sol_v(i-1)    
              enddo
              b(1) = r*sol_v(1+1)+(1-2*r)*sol_v(1)
              b(N) = (1-2*r)*sol_v(N)+r*sol_v(N-1)
                
              !Solution of implicit part
              call solve(N, a, c, c, b, sol_n)
              !update solution
              sol_v = sol_n
          enddo

          !Save on file
          do i = 1, N
              write(2,*) sol_v(i)
          enddo
         
      enddo

      deallocate(sol_v, sol_n, a, b, c)
      close(1)
      close(2)
      end program calore

C=============================================================================
C Subroutine for the solution of system with tridiagonal matrix
C=============================================================================

      subroutine solve(N, diag, sup_diag, inf_diag, b, x)
C=============================================================================
C     Subroutine for the solution of system with tridiagonal matrix
C     via gauss elimination. Computational order O(N)
C
C     Parameters
C     ----------
C     N : init
C         size of matrix A (N x N)
C     diag : array
C         A(i, i) = diag
C     sup_diag : array
C         A(i, i + 1) = sup_diag
C     inf_diag : array
C         A(i - 1, i) = inf_diag
C     b : array
C         RHS of the equation
C     x : array
C         solution of the system: x = A^{-1} b
C=============================================================================
      implicit real*8 (a-h,o-z)
      real*8, dimension(N) :: diag, sup_diag, inf_diag
      real*8, dimension(N) :: d_ii, d_up, d_lo, b, x

      ! We write the information on these arrays
      ! so that the subroutine doesn't modify the old ones
      d_ii = diag       !d_ii is the diagonal of the matrix
      d_up = sup_diag   !d_up is the upper diagonal
      d_lo = inf_diag   !d_lo is the lower diagonal

      ! make the matrix an upper triangular by canceling d_lo
      do i = 2, N
         if(d_ii(i-1) == 0.0) then
              print *,"Division by zero, non invertible matrix"
              return
          endif
          a = d_lo(i-1)/d_ii(i-1)
          d_ii(i) = d_ii(i) - a*d_up(i-1)
          b(i) = b(i) - a*b(i-1)
      enddo

      ! at this point it is easy to find x_n given that
      ! in the last row of the matrix, there is only
      ! the element left on the diagonal which is non-zero
      if(d_ii(N) == 0.0) then
          print *,"Division by zero, non invertible matrix"
          return
      endif

      x(N) = b(N)/d_ii(N)

      !and now going backwards I can find all the other x_i
      do i = N-1, 1, -1
          b(i) = b(i) - d_up(i)*x(i+1)
          if(d_ii(i) == 0.0) then
              print *,"Division by zero, non invertible matrix"
              return
          endif
          x(i) = b(i)/d_ii(i)
      enddo

      return
      end
