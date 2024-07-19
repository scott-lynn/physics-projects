!Duffing Oscillator v5
module subprograms
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)
    real(kind=dp), parameter :: pi = 4.D0*DATAN(1.D0) !Standard Fortran implementation of pi
    !Module contains all numerical methods and differential equations required.


    contains

        function dydt(tn,yn,zn,params)
        !Evaluates the dy/dt part of the Duffing Oscillator equation for the given input values

            real(kind=dp) :: dydt
            real(kind=dp), intent(in) ::  tn, yn, zn, params

            !Calculate derivative dy/dt and return value
            dydt = zn

        end function

        function dzdt(tn,yn,zn,params)
        !Evaluates the dz/dt part of the Duffing Oscillator equation for the given input values

            real(kind=dp) :: dzdt
            real(kind=dp), intent(in) :: tn, yn, zn
            real(kind=dp), dimension(5) :: params
            real(kind=dp) :: A, B, C, D, omega

            A = params(1)
            B = params(2)
            C = params(3)
            D = params(4)
            omega = params(5)

            !Calculate derivative dz/dt and return value
            dzdt = -A*yn**3.0_dp + B*yn - C*zn + D*sin(omega*tn)

        end function

        function potential_energy(params,yn)
        !Calculates the external double well potential at current position yn

            real(kind=dp) :: potential_energy
            real(kind=dp), dimension(5) :: params
            real(kind=dp), intent(in) :: yn
            real(kind=dp) :: A, B, V

            A = params(1)
            B = params(2)

            !Calculate external potential V
            V = (A/4.0_dp)*(yn**4) - (B/2.0_dp)*(yn**2.0_dp)

            !Return value of V
            potential_energy = V

        end function

        function kinetic_energy(m,v)
        !Calculates the kinetic energy of the particle for its current value of velocity

            real(kind=dp) :: kinetic_energy
            real(kind=dp), intent(in) :: m, v
            real(kind=dp) :: K

            !Calculate kinetic energy of particle
            K = 0.5_dp*m*(v**2)

            kinetic_energy = K

        end function

        function total_energy(K,V)
        !Calculates the total energy given kinetic energy K and potential energy V

            real(kind=dp) :: total_energy
            real(kind=dp), intent(in) :: K,V
            real(kind=dp) :: E

            !Calculate total energy
            E = K + V

            total_energy = E

        end function

        function shm_y(tn,yn,zn)
        !Evaluates the dydt part of the shm equation

            real(kind=dp) :: shm_y
            real(kind=dp), intent(in) :: tn, yn, zn
            real(kind=dp) :: dydt

            !Calculate derivative dy/dt
            dydt = zn

            !Return value
            shm_y = dydt

        end function

        function shm_z(tn,yn,zn)
        !Evaluates the dz/dt part of the shm equation

            real(kind=dp) :: shm_z
            real(kind=dp), intent(in) :: tn, yn, zn
            real(kind=dp) :: omega
            real(kind=dp) :: dzdt

            omega = 2*pi

            !Calculate derivative dz/dt
            dzdt = -omega**2.0_dp*yn

            !Return value
            shm_z = dzdt

        end function

        subroutine rk2(dt,tn,yn,zn,fy,fz,y_next,z_next,params)
        !Calculates yn+1 and zn+1 using the RK2 method given a second order ODE split into two functions fy and fz. Requires inputs of time tn and time_step dt, alongside yn and zn.
        !Takes optional argument params, which is a 1D array of length 5, containing A,B,C,D,omega values for the Duffing Oscillator equation
            real(kind=dp), intent(in) :: tn, yn, zn, dt
            real(kind=dp), dimension(5), intent(in), optional :: params
            real(kind=dp) :: fy, fz
            real(kind=dp) :: y_next, z_next
            real(kind=dp) :: k1y, k2y
            real(kind=dp) :: k1z, k2z

            !Determine if constant is required
            if(present(params)) then

                !Calculate k1 values for y and z
                k1y = fy(tn,yn,zn,params) * dt
                k1z = fz(tn,yn,zn,params) * dt

                !Calculate k2 values for y and z
                k2y = fy(tn + dt/2, yn + k1y/2, zn + k1z/2, params) * dt
                k2z = fz(tn + dt/2, yn + k1y/2, zn + k1z/2, params) * dt

                !Calculate yn+1 and zn+1 - return these in subroutine output
                y_next = yn + k2y
                z_next = zn + k2z
            else
                !Calculate k1 values for y and z
                k1y = fy(tn,yn,zn) * dt
                k1z = fz(tn,yn,zn) * dt

                !Calculate k2 values for y and z
                k2y = fy(tn + dt/2, yn + k1y/2, zn + k1z/2) * dt
                k2z = fz(tn + dt/2, yn + k1y/2, zn + k1z/2) * dt

                !Calculate yn+1 and zn+1 - return these in subroutine output
                y_next = yn + k2y
                z_next = zn + k2z
            end if




        end subroutine

        subroutine rk4(dt,tn,yn,zn,fy,fz,y_next,z_next, params)
        !Calculates yn+1 and zn+1 using the RK2 method given a second order ODE split into two functions fy and fz. Requires inputs of time tn and time_step dt, alongside yn and zn.
        !Takes optional argument params, which is a 1D array of length 5, containing A,B,C,D,omega values for the Duffing Oscillator equation

            real(kind=dp), intent(in) :: tn, yn, zn, dt
            real(kind=dp), dimension(5), intent(in), optional :: params
            real(kind=dp) :: fy, fz
            real(kind=dp) :: y_next, z_next
            real(kind=dp) :: k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z

            !Determine if constant is required
            if(present(params)) then

                !Calculate k1 values for y and z
                k1y = fy(tn,yn,zn,params) * dt
                k1z = fz(tn,yn,zn,params) * dt

                !Calculate k2 values for y and z
                k2y = fy(tn + dt/2, yn + k1y/2, zn + k1z/2, params) * dt
                k2z = fz(tn + dt/2, yn + k1y/2, zn + k1z/2, params) * dt

                !Calculate k3 values for y and z
                k3y = fy(tn + dt/2, yn + k2y/2, zn + k2z/2, params) * dt
                k3z = fz(tn + dt/2, yn + k2y/2, zn + k2z/2, params) * dt

                !Calulate k4 values for y and z
                k4y = fy(tn + dt, yn + k3y, zn + k3z, params) * dt
                k4z = fz(tn + dt, yn + k3y, zn + k3z, params) * dt

                !Calculate yn+1 and zn+1
                y_next = yn + k1y/6 + k2y/3 + k3y/3 + k4y/6
                z_next = zn + k1z/6 + k2z/3 + k3z/3 + k4z/6

            else

                !Calculate k1 values for y and z
                k1y = fy(tn,yn,zn) * dt
                k1z = fz(tn,yn,zn) * dt

                !Calculate k2 values for y and z
                k2y = fy(tn + dt/2, yn + k1y/2, zn + k1z/2) * dt
                k2z = fz(tn + dt/2, yn + k1y/2, zn + k1z/2) * dt

                !Calculate k3 values for y and z
                k3y = fy(tn + dt/2, yn + k2y/2, zn + k2z/2) * dt
                k3z = fz(tn + dt/2, yn + k2y/2, zn + k2z/2) * dt

                !Calulate k4 values for y and z
                k4y = fy(tn + dt, yn + k3y, zn + k3z) * dt
                k4z = fz(tn + dt, yn + k3y, zn + k3z) * dt

                !Calculate yn+1 and zn+1
                y_next = yn + k1y/6 + k2y/3 + k3y/3 + k4y/6
                z_next = zn + k1z/6 + k2z/3 + k3z/3 + k4z/6

            end if

        end subroutine

        function shm_analytical(t,y0,v0)
        !This function calculates the analytical solution to the SHM equation at time t
        !Takes input parameters, time tn, initial position y0, and initial velocity v0

            real(kind=dp) :: shm_analytical
            real(kind=dp), intent(in) :: y0, t, v0
            real(kind=dp) :: y, omega

            !Set value of omega
            omega = 2*pi

            !Calculate y(t)
            y = y0 * cos(omega*t)+(v0/omega)*sin(omega*t)

            !Return
            shm_analytical = y

        end function

end module subprograms
!-------------------------------------------------------------------------------------------
program duffosc
    use subprograms
    implicit none
    !This program calculates the solution to the Duffing Oscillator using the RK4 method, by default it only outputs the poincare section for the specified parameters.
    !Other outputs are commented out but can be re-enabled through the use of an editor.

    !I/O and system
    integer :: ios, out_a, out_b, out_c, out_d, out_e, out_f, out_g, out_h, out_i, out_j, out_k, out_l
    integer :: in_a
    integer :: switch
    real(kind=dp) :: start_rk2, finish_rk2, comptime_rk2
    real(kind=dp) :: start_rk4, finish_rk4, comptime_rk4

    !Declare Variables for rk2 and rk4
    real(kind=dp) :: t0, tn, dt, t_max, n, period
    real(kind=dp) :: y0, z0, v0
    real(kind=dp) :: yn, zn
    real(kind=dp) :: y_true
    real(kind=dp) :: rk2_diff, rk4_diff

    !Duffing Oscillator Variables
    real(kind=dp) :: K, V, E, m
    real(kind=dp) :: A, B, C, D
    real(kind=dp) :: omega
    real(kind=dp), dimension(5) :: params

    !Define I/O paths
    out_a = 11
    out_b = 12
    out_c = 13
    out_d = 14
    out_e = 15
    out_f = 16
    out_g = 17
    out_h = 18
    out_i = 19
    out_j = 20
    out_k = 21
    out_l = 22
    in_a = 23

    !Time variables
    dt = 1.0e-4_dp !Time-step value
    t0 = 0.0_dp
    tn = t0
    t_max = 100.0_dp !Max time
    
    !SHM Variables
    y0 = 100.0_dp
    z0 = 0.0_dp
    v0 = 0.0_dp
    yn = y0
    zn = z0

    !Read if user wants to run Runge-Kutta testing
    open(unit=in_a, file='switch.dat', status='old', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file switch.dat'
    read(unit=in_a, fmt=*, iostat=ios) switch
    if(ios .ne. 0) stop 'Error in reading file switch.dat'
    close(unit=out_a, iostat=ios)
    if(ios .ne. 0) stop 'Error closing file switch.dat'

    !Toggle RK2 and RK4 test (speeds up program for duffosc)
    if (switch .ne. 0) then 
    !Analytical Solution to SHM ODE ------------------------------------------------------
    open(unit=out_a, file='shm_test_true.dat', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file shm_true.dat'

    do while (tn <= t_max)

        !Calculate analytical solution
        y_true = shm_analytical(tn, y0, v0)

        !Write out data
        write(unit=out_a, fmt=*, iostat=ios) tn, y_true
        if(ios .ne. 0) stop 'Error writing to file shm_true.dat'

        tn = tn + dt
    end do

    close(unit=out_a, iostat=ios)
    if(ios .ne. 0) stop 'Error closing file shm_true.dat'

    !RK2 test ----------------------------------------------------------------------------
    open(unit=out_a, file='shm_rk2.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_rk2.dat'

    !Initialise parameters for RK2
    tn=t0
    yn=y0
    zn=z0

    !Start recording computation time
    call cpu_time(start_rk2)

    !RK2 loop
    do while (tn <= t_max)

        !Calculate yn+1 and zn+1 using RK2
        call rk2(dt, tn, yn, zn, shm_y, shm_z, yn, zn)

        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) tn, yn
        if(ios.ne.0) stop 'Error writing to file shm_rk2.dat'

        !Go to next time step tn+1
        tn = tn + dt
    end do

    !Finish recording computation time
    call cpu_time(finish_rk2)
    comptime_rk2 = finish_rk2 - start_rk2

    !Calculate difference in solutions
    rk2_diff = abs(y_true - yn)

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_rk2.dat'

    print *, 'RK2 SHM ----------------------------'
    print *, 'y(t) =', y_true
    print *, 'y_rk2 =', yn
    print *, 'RK2 diff =', rk2_diff
    print *, 'RK2 time:', comptime_rk2

    !Write out dt against computational time
    open(unit=out_a, file='rk2_shm_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rk2_shm_comptime_v_dt.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, comptime_rk2
    if(ios.ne.0) stop 'Error writing to file rk2_shm_comptime_v_dt.dat'
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rk2_shm_comptime_v_dt.dat'

    !Write out dt against accuracy
    open(unit=out_b, file='rk2_shm_acc_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rk2_shm_acc_v_dt.dat'
    write(unit=out_b, fmt=*, iostat=ios) dt, rk2_diff
    if(ios.ne.0) stop 'Error writing to file rk2_shm_acc_v_dt.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rk2_shm_acc_v_dt.dat'

    !RK4 test -----------------------------------------------------------------------
    open(unit=out_a, file='shm_rk4.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_rk4.dat'

    !Initialise parameters for RK4
    tn=t0
    yn=y0
    zn=z0

    !Start recording computation time for RK4
    call cpu_time(start_rk4)

    !RK4 loop
    do while (tn <= t_max)

        !Calculate yn+1 and zn+1 using RK4 method
        call rk4(dt, tn, yn, zn, shm_y, shm_z, yn, zn)

        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) tn, yn
        if(ios.ne.0) stop 'Error writing to file shm_rk4.dat'

        !Go to next time step tn+1
        tn = tn + dt

    end do

    !Finish recording computation time for RK4
    call cpu_time(finish_rk4)

    comptime_rk4 = finish_rk4 - start_rk4
    rk4_diff = abs(y_true - yn)

    print *, 'RK4 SHM----------------------------'
    print *, 'y(t) =', y_true
    print *, 'y_rk4 =', yn
    print *, 'RK4 diff =', rk4_diff
    print *, 'RK4 time:', comptime_rk4
    print *, ''

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_rk4.dat'

    !Write out dt against computational time
    open(unit=out_a, file='rk4_shm_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rk4_shm_comptime_v_dt.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, comptime_rk4
    if(ios.ne.0) stop 'Error writing to file rk4_shm_comptime_v_dt.dat'
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rk4_shm_comptime_v_dt.dat'

    !Write out dt against accuracy
    open(unit=out_b, file='rk4_shm_acc_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rk4_shm_acc_v_dt.dat'
    write(unit=out_b, fmt=*, iostat=ios) dt, rk4_diff
    if(ios.ne.0) stop 'Error writing to file rk4_shm_acc_v_dt.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rk4_shm_acc_v_dt.dat'

end if 

!Duffing Oscillator
!---------------------------------------------------------------------------------
    print *, 'Running Duffing Oscillator Simulation'
    print *, ''

    !Initialise Variables
    n = 0.0_dp
    tn = 0.0_dp
    t_max = 10000.0_dp !Max time

    open(unit=in_a, file='params.dat', status='old', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file params.dat'

    !Read in parameter values from file
    read(unit=in_a, fmt=*, iostat=ios) A       
    if(ios .ne. 0) stop 'Error in reading file params.dat'
    read(unit=in_a, fmt=*, iostat=ios) B        
    if(ios .ne. 0) stop 'Error in reading file params.dat'
    read(unit=in_a, fmt=*, iostat=ios) C       
    if(ios .ne. 0) stop 'Error in reading file params.dat'
    read(unit=in_a, fmt=*, iostat=ios) D       
    if(ios .ne. 0) stop 'Error in reading file params.dat'
    read(unit=in_a, fmt=*, iostat=ios) omega    
    if(ios .ne. 0) stop 'Error in reading file params.dat'
    
    close(unit=out_a, iostat=ios)
    if(ios .ne. 0) stop 'Error closing file params.dat'

    !Set up parameter array
    params(1) = A       !Restoring force coefficient
    params(2) = B       !Pushes particle from origin
    params(3) = C       !Damping coefficient
    params(4) = D       !Magnitude of driving force
    params(5) = omega   !Angular frequency of driving force

    yn = 0.0_dp     !Initial position
    zn = 0.0_dp     !Initial velocity
    m = 1.0_dp      !Mass of particle =1, dimensionless

    period = ((2.0_dp*pi)*n)/omega !Period, T

    !Print to console
    print *, 'Input parameters:  ', 'A=', A, 'B=', B, 'C=', C, 'D=', D, 'omega=', omega
    print *, 'Initial values:    ', 'y0 =', yn, 'v0 =', zn
    print *, ''

    !Open files (Comment to enable)
    ! open(unit=out_a, file='trajectory.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file trajectory.dat'
    ! open(unit=out_b, file='phase_space.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file phase_space.dat'
    ! open(unit=out_e, file='potential_v_time.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file potential_v_time.dat'
    ! open(unit=out_f, file='kinetic_energy.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file kinetic_energy.dat'
    ! open(unit=out_g, file='total_energy.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file total_energy.dat'
    ! open(unit=out_h, file='potential_v_displacement.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file potential_v_displacement.dat'
    ! open(unit=out_i, file='energy_v_displacement.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file energy_v_displacement.dat'
    ! open(unit=out_j, file='kinetic_v_displacement.dat', iostat=ios)
    ! if(ios.ne.0) stop 'Error in opening file kinetic_v_displacement.dat'
    open(unit=out_k, file='poincare.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file poincare.dat'
    !open(unit=out_l, file='poincare_position_time.dat', iostat=ios)
    !if(ios.ne.0) stop 'Error in opening file poincare_position_time.dat'

    !Start recording computation time
    print *, ''
    call cpu_time(start_rk4)

        !Duffing Oscillator RK4 loop
        do while (tn <= t_max)

            !Calculate yn+1 and zn+1 using RK4
            call rk4(dt, tn, yn, zn, dydt, dzdt, yn, zn, params)

            !Calculate potential energy
            V = potential_energy(params,yn)

            !Calculate kinetic energy
            K = kinetic_energy(m,zn)

            !Calculate total energy
            E = total_energy(K,V)

            !poincare section
            if (period - tn <= 1.0e-10_dp) then
                
                !Check if 10% of total time has passed and write poincare data
                if (tn >= 0.10*t_max) then
                    write(unit=out_k, fmt=*, iostat=ios) yn, zn, tn
                    if(ios.ne.0) stop 'Error writing to file poincare.dat'
                    !write(unit=out_l, fmt=*, iostat=ios) tn, yn
                    !if(ios.ne.0) stop 'Error writing to file poincare_position_time.dat'
                end if
                !Calculate time of next period
                n = n + 1.0_dp
                period = ((2.0_dp*pi)*n)/omega
            end if

            !Write output to file (Comment to enable)
            ! write(unit=out_a, fmt=*, iostat=ios) tn, yn
            ! if(ios.ne.0) stop 'Error writing to file trajectory.dat'
            ! write(unit=out_b, fmt=*, iostat=ios) yn, zn, tn
            ! if(ios.ne.0) stop 'Error writing to file phase_space.dat'
            ! write(unit=out_c, fmt=*, iostat=ios) tn, zn
            ! if(ios.ne.0) stop 'Error writing to file potential_v_time.dat'
            ! write(unit=out_f, fmt=*, iostat=ios) tn, K
            ! if(ios.ne.0) stop 'Error writing to file kinetic_energy.dat'
            ! write(unit=out_g, fmt=*, iostat=ios) tn, E
            ! if(ios.ne.0) stop 'Error writing to file total_energy.dat'
            ! write(unit=out_h, fmt=*, iostat=ios) yn, V
            ! if(ios.ne.0) stop 'Error writing to file potential_v_displacement.dat'
            ! write(unit=out_i, fmt=*, iostat=ios) yn, E
            ! if(ios.ne.0) stop 'Error writing to file energy_v_displacement.dat'
            ! write(unit=out_j, fmt=*, iostat=ios) yn, K
            ! if(ios.ne.0) stop 'Error writing to file kinetic_v_displacement.dat'

            !Go to next time step tn+1
            tn = tn + dt
        end do

    !Finish recording computation time
    call cpu_time(finish_rk4)
    comptime_rk4 = finish_rk4 - start_rk4
    print *, 'Comp time =', comptime_rk4

    !Close files (Comment to enable)
    ! close(unit=out_a, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file trajectory.dat'
    ! close(unit=out_b, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file phase_space.dat'
    ! close(unit=out_c, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file potential_v_time.dat'
    ! close(unit=out_f, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file kinetic_energy.dat'
    ! close(unit=out_g, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file total_energy.dat'
    ! close(unit=out_h, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file potential_v_displacement.dat'
    ! close(unit=out_i, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file energy_v_displacement.dat'
    ! close(unit=out_j, iostat=ios)
    ! if(ios.ne.0) stop 'Error closing file kinetic_v_displacement.dat'
    close(unit=out_k, iostat=ios)
    if(ios.ne.0) stop 'Error closing file poincare.dat'
    !close(unit=out_l, iostat=ios)
    !if(ios.ne.0) stop 'Error closing file poincare_position_time.dat'
end program duffosc
