      program wavelaxwendroff
C==================================================================================
C     Code to solve the wave equation: d^2u/dt^2 + vd^2u/dx^2 = 0
C     The code uses lax wendroff's method indeed, if we call r = vdu/dx, s = du/dt
C     we can rewrite the equation as a transport equation:
C
C     dU    dF(U)
C     --  + -----  = 0     where :
C     dt      dx
C
C         | r |               | 0    -v |
C     U = |   |    and F(U) = |         | U
C         | s |               | -v    0 |
C
C     To plot the result you can use wave_1D_plot.py
c==================================================================================
      implicit real*8 (a-h,p-z)

      real*8, dimension(:), allocatable :: r_v, r_n, s_v, s_n
      real*8, dimension(:), allocatable :: sol_v, sol_n

      open(1, file='input_wave.txt', status='old')  !parameter file
      open(2, file='wave_lw.dat', status='unknown') !results file

      read(1,*) N       !Number of ponits on spatial  axis
      read(1,*) i_T     !Number of ponits on temporal axis
      read(1,*) v       !Speed of the propagation
      read(1,*) dt      !Temporal step
      read(1,*) dx      !Spatial  step
      read(1,*) i_P     !How often to save the solution

      allocate(r_v(0:N), r_n(0:N), s_v(0:N), s_n(0:N))
      allocate(sol_v(0:N), sol_n(0:N))

      alpha = v*dt/dx
      print*, alpha
        
      !Initial condition
      q = 2.d0 * 4.d0 * datan(1.d0) !2 pi
      do i = 0, N
          x = i*dx
          !u(x, t=0)
          sol_v(i) = dexp(-((x-0.2d0)/0.05d0)**2.d0)!dsin(q*0.3d0*x)
          !v*du/dx
          r_v(i) =-v*dexp(-((x-0.2d0)/0.05d0)**2.d0)*800.d0*(x-0.2d0)!v*q*0.3d0*dcos(q*0.3d0*x)
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
     &            0.5*alpha**2*(r_v(i+1) - 2.0*r_v(i) + r_v(i-1))

                  s_n(i) = s_v(i) + 0.5*alpha*(r_v(i+1) - r_v(i-1))+
     &            0.5*alpha**2*(s_v(i+1) - 2.0*s_v(i) + s_v(i-1))
              enddo

              ! Periodic boundary conditions
              r_n(0) = r_n(N - 1)
              r_n(N) = r_n(1)
              s_n(0) = s_n(N - 1)
              s_n(N) = s_n(1)

              ! For absorbing boundary condition (zero derivative):
              !r_n(0) = r_n(1)
              !r_n(N) = r_n(N-1)
              !s_n(0) = s_n(1)
              !s_n(N) = s_n(N-1)

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
      close(1)
      close(2)
      stop

      end program
