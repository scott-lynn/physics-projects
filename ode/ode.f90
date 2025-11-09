! ODE v1.8
! This program solves both the radioactive decay ODE and the SHM ODE
! Contains Euler and Leapfrog methods
module methods
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)
    contains

        function rad_decay_ODE(N,t,lambda)
        !This function calculates the derivative dN/dt for the radioactive decay ODE.
        !It takes two inputs, the value of N and the time t, which is a dummy argument so that the function is of the form that can be used with the numerical methods later within the program
            real(kind=dp) :: rad_decay_ODE
            real(kind=dp), intent(in) :: N, t, lambda
            real(kind=dp) :: dNdt

            !Calculate the differential
            dNdt = (-lambda*N)

            !Return value of dNdt
            rad_decay_ODE = dNdt

        end function rad_decay_ODE

        function rad_decay(N0,t,lambda)
        !This function gives the analytical solution for the radioactive decay equation given input values of N0, t and lambda

            real(kind=dp) :: rad_decay
            real(kind=dp), intent(in) :: N0, t, lambda
            real(kind=dp) :: N

            !Calculate N(t)
            N = N0*exp(-lambda*t)
            rad_decay = N

        end function rad_decay

        function SHM_ODE(y,t,omega)
        ! Contains part of the coupled ODE for the SHM SODE of dz/dt
            real(kind=dp) :: SHM_ODE
            real(kind=dp), intent(in) :: y, t, omega
            real(kind=dp) :: dzdt

            !Calculate the differential
            dzdt = -omega**2.0_dp*y

            !Return value of dzdt
            SHM_ODE = dzdt

        end function SHM_ODE

        function SHM(y0,t,v0,omega)
        !This function calculates the analytical solution to the SHM equation at time t

        real(kind=dp) :: SHM
        real(kind=dp), intent(in) :: y0, t, v0, omega
        real(kind=dp) :: y

        !Calculate y(t)
        y = y0 * cos(omega*t)+(v0/omega)*sin(omega*t)
        SHM = y
        end function SHM

        function euler_method(dt,f,y0,t,c_in)
        !This function uses the Euler numerical scheme to calculate the value of an input function f evaluated at input values y0, t and increment dt.
        !Can pass an optional argument c_in, which represents a constant in the equation, if this is not needed it defaults to zero

            real(kind=dp) :: euler_method

            real(kind=dp), intent(in) :: y0, dt
            real(kind=dp), intent(in), optional :: c_in
            real(kind=dp) ::  t, y, y_next, c
            real(kind=dp) :: f !input function f(y,t,c)

            !Determine if constant is required
            if(present(c_in)) then
                c = c_in
            else
                c = 0.0_dp
            end if

            !Set value of y to input value
            y = y0

            !Calculate value of function at the next step
            y_next = y + f(y,t,c)*dt

            !Return this
            euler_method = y_next

        end function euler_method

        function leapfrog(dt,f,y,y_prev,t,c_in)
        !This function uses the Leapfrog numerical scheme to calculate the value of an input function f evalued at input values y_prev, y, t, and increment dt.
        !Can pass an optional argument c_in which represents a constant in the equation, if this is not needed it defaults to zero.
        !Two starting values are required, y and y_prev, it is recommended that the function euler_method is used to generate these values.
            real(kind=dp) :: leapfrog

            real(kind=dp), intent(in) :: y, y_prev, dt
            real(kind=dp), intent(in), optional :: c_in
            real(kind=dp) :: t, y_next, c
            real(kind=dp) :: f !input function f(y,t,c)

            !Determine if constant is required
            if(present(c_in)) then
                c = c_in
            else
                c = 0.0_dp
            end if

            !Calculate value of function at next step
            y_next = y_prev +2*f(y,t,c)*dt

            !Return value
            leapfrog = y_next

        end function leapfrog
end module methods
!----------------------------------------------------------------------
program ode
    use methods
    implicit none

    !I/O
    integer :: ios, out_a, out_b
    real(kind=dp) :: start, finish

    !Declare Variables
    real(kind=dp) :: temp
    real(kind=dp) :: dt, t, t_max, t0, diff

    !Variables for Radioactive Decay equation
    real(kind=dp) :: N, N_true, N0, N_prev, lambda

    !Variabes for SHM equation
    real(kind=dp), parameter :: pi = 4.D0*DATAN(1.D0)
    real(kind=dp) :: omega !Angular frequency
    real(kind=dp) :: v0 !Initial velocity: dy/dt at t=0
    real(kind=dp) :: z, z_prev, z_temp, z0
    real(kind=dp) :: y, y_prev, y0, y_true, y_temp

    !-----------------------------------------------------------------------------
    !I/O Settings
    out_a = 11 !Channel 1
    out_b = 12 !Channel 2

    !------------------------------------------------------------------------------
    !Time values for all ODE solvers

    !Time values
    dt = 0.001_dp  !Time step
    t0 = 0.0_dp !Start time
    t_max = 10.0_dp !Stop time

    !--------------------------------------------------------------------------------
    !Radioactive Decay Settings

    !Initial values
    N0 = 100.0_dp   !Initial value of N
    N = N0
    lambda = 1.0_dp !Radioactive decay constant

    !--------------------------------------------------------------------------------
    !SHM Settings

    !Initial values
    y0 = 1.0_dp
    z0 = 0.0_dp
    omega = 2*pi
    v0 = 0.0_dp

    !--------------------------------------------------------------------------------
    !Radioactive Decay ODE with Euler Method

    !Open file to write euler results
    open(unit=out_a, file='rd_euler_true.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file rd_euler_true.dat'

    open(unit=out_b, file='rd_euler_approx.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file rd_euler_approx.dat'

    !Start recording computational time
    call cpu_time(start)

    !Calculate approximation for the radioactive decay equation evaluated at time point t_max
    do while (t <= t_max)

        !Calculate analytical solution
        N_true = rad_decay(N0,t,lambda)

        !Approximate solution with Euler method
        N = euler_method(dt,rad_decay_ODE,N,t,lambda)

        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) t, N_true
        if(ios.ne.0) stop 'Error writing to file rd_euler_true.dat'

        write(unit=out_b, fmt=*, iostat=ios) t, N
        if(ios.ne.0) stop 'Error writing to file rd_euler_approx.dat'

        !Update time
        t = t + dt

    end do

    !Finish recording computation time
    call cpu_time(finish)

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_euler_true.dat'

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_euler_approx.dat'

    diff = abs(N_true - N)

    !Print output
    print *, '--------------------------------------------------'
    print *, 'Radioactive Decay ODE with Euler Method:'
    print '("Computation time:",f6.3," seconds.")',finish-start
    print *, 'Evaluated at time t=', t_max
    print *, 'dt:', dt
    print *, 'N_approx(t):', N
    print *, 'N_true(t):', N_true
    print *, 'Difference in solutions:', diff
    print *, '--------------------------------------------------'

    !Write out dt against computational time
    open(unit=out_a, file='rd_euler_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rd_euler_dt_v_comptime.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, finish-start
    if(ios.ne.0) stop 'Error writing to file rd_euler_dt_v_comptime.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_euler_dt_v_comptime.dat.dat'

    !Write out dt against difference in solutions
    open(unit=out_a, file='rd_euler_dt_v_diff.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rd_euler_dt_v_diff.dat.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, diff
    if(ios.ne.0) stop 'Error writing to file rd_euler_dt_v_diff.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_euler_dt_v_diff.dat.dat'

    !------------------------------------------------------------------------------------------------------------
    !Radioactive Decay ODE with Leapfrog Method

    !Open file to write results to
    open(unit=out_a, file='rd_leapfrog_true.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file rd_leapfrog_true.dat'

    open(unit=out_b, file='rd_leapfrog_approx.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file rd_leapfrog_approx.dat'

    !Reset values
    N = N0
    t = 0.0_dp

    !Generate y and y_prev using one step of Euler Method
    N = euler_method(dt,rad_decay_ODE,N,t,lambda)

    N_prev = N0

    !Start recording computational time
    call cpu_time(start)

    !Calculate approximation for the radioactive decay equation at time t_max
    do while (t <= t_max)

        !Calculate analytical solution
        N_true = rad_decay(N0,t,lambda)

        temp = N
        !Calculate next step using Leapfrog Method
        N = leapfrog(dt,rad_decay_ODE,N,N_prev,t,lambda)

        !Pass held value to update N_prev
        N_prev = temp

        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) t, N_true
        if(ios.ne.0) stop 'Error writing to file rd_leapfrog_true.dat'

        write(unit=out_b, fmt=*, iostat=ios) t, N
        if(ios.ne.0) stop 'Error writing to file rd_leapfrog_approx.dat'

        !Update time
        t = t + dt

    end do

    !Finish recording computation time
    call cpu_time(finish)

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_leapfrog_true.dat'

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_leapfrog_approx.dat'

    diff = abs(N_true - N)

    !Print output
    print *, '--------------------------------------------------'
    print *, 'Radioactive Decay ODE with Leapfrog Method:'
    print '("Computation time:",f6.3," seconds.")',finish-start
    print *, 'Evaluated at time t=', t_max
    print *, 'dt:', dt
    print *, 'N_approx(t):', N
    print *, 'N_true(t):', N_true
    print *, 'Difference in solutions:', diff
    print *, '--------------------------------------------------'

    !Write out dt against computational time
    open(unit=out_a, file='rd_leapfrog_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rd_leapfrog_dt_v_comptime.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, finish-start
    if(ios.ne.0) stop 'Error writing to file rd_leapfrog_dt_v_comptime.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_leapfrog_dt_v_comptime.dat.dat'

    !Write out dt against difference in solutions
    open(unit=out_a, file='rd_leapfrog_dt_v_diff.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file rd_leapfrog_dt_v_diff.dat.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, diff
    if(ios.ne.0) stop 'Error writing to file rd_leapfrog_dt_v_diff.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_leapfrog_dt_v_diff.dat.dat'

    !--------------------------------------------------------------------------
    !SHM with Euler Method

    !Open file to write results to
    open(unit=out_a, file='shm_euler_true.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_euler_true.dat'

    open(unit=out_b, file='shm_euler_approx.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_euler_approx.dat'

    !Reset Values
    z=z0
    y=y0
    t = t0

    !Start recording computational time
    call cpu_time(start)

    !Calculate Euler approximation for SHM equation
    do while (t < t_max)

        !Calculate analytical solution
        y_true = SHM(y0,t,v0,omega)

        !Hold this value for use in calculating yn+1
        z_temp = z

        !Calculate next steps using the Euler method
        z = z + SHM_ODE(y,t,omega)*dt    !zn+1
        y = y + z_temp*dt               !yn+1


        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) t, y_true
        if(ios.ne.0) stop 'Error writing to file shm_euler_true.dat'

        write(unit=out_b, fmt=*, iostat=ios) t, y
        if(ios.ne.0) stop 'Error writing to file shm_euler_approx.dat'

        t = t + dt !Go to next time step

        end do

    !Finish recording computation time
    call cpu_time(finish)

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_euler_true.dat'

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file rd_euler_approx.dat'

    diff = abs(y_true - y)

    !Print output
    print *, '--------------------------------------------------'
    print *, 'SHM SODE with Euler Method:'
    print '("Computation time:",f6.3," seconds.")',finish-start
    print *, 'Evaluated at time t=', t_max
    print *, 'dt:', dt
    print *, 'y_approx(t):', y
    print *, 'y_true(t):', y_true
    print *, 'Difference in solutions:', diff
    print *, '--------------------------------------------------'

    !Write out dt against computational time
    open(unit=out_a, file='shm_euler_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file shm_euler_dt_v_comptime.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, finish-start
    if(ios.ne.0) stop 'Error writing to file shm_euler_dt_v_comptime.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_euler_dt_v_comptime.dat.dat'

    !Write out dt against difference in solutions
    open(unit=out_a, file='shm_euler_dt_v_diff.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file shm_euler_dt_v_diff.dat.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, diff
    if(ios.ne.0) stop 'Error opening file shm_euler_dt_v_diff.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_euler_dt_v_diff.dat.dat'

    !------------------------------------------------------------------------------
    !SHM with Leapfrog Method

    !Open file to write results to
    open(unit=out_a, file='shm_leapfrog_true.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_leapfrog_true.dat'

    open(unit=out_b, file='shm_leapfrog_approx.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file shm_leapfrog_approx.dat'

    !Reset Values
    z=z0
    y=y0
    t = t0

    !Start recording computational time
    call cpu_time(start)

    !Generate z, z_prev and y, y_prev using one step of euler method

    z_prev = z !zn-1
    y_prev = y !yn-1

    z = z - SHM_ODE(y,t,omega)*(dt/2)   !zn
    y = y + z_prev*(dt/2)               !yn

    !Calculate Leapfrog approximation for SHM equation
    do while (t < t_max)

        !Calculate analytical solution
        y_true = SHM(y0,t,v0,omega)

        !Hold this value
        z_temp = z !zn
        y_temp = y !yn

        !Calculate next steps using the Leapfrog method
        z = z_prev + 2*SHM_ODE(y,t,omega)*dt    !zn+1
        y = y_prev + 2*z_temp*dt               !yn+1

        !n becomes n-1
        z_prev = z_temp
        y_prev = y_temp

        !Write output to file
        write(unit=out_a, fmt=*, iostat=ios) t, y_true
        if(ios.ne.0) stop 'Error writing to file shm_leapfrog_true.dat'

        write(unit=out_b, fmt=*, iostat=ios) t, y
        if(ios.ne.0) stop 'Error writing to file shm_leapfrog_approx.dat'

        !Go to next time step
        t = t + dt

        end do

    !Finish recording computation time
    call cpu_time(finish)

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_leapfrog_true.dat'

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_leapfrog_approx.dat'

    diff = abs(y_true - y)

    !Print output
    print *, '--------------------------------------------------'
    print *, 'SHM SODE with Leapfrog Method:'
    print '("Computation time:",f6.3," seconds.")',finish-start
    print *, 'Evaluated at time t=', t_max
    print *, 'dt:', dt
    print *, 'y_approx(t):', y
    print *, 'y_true(t):', y_true
    print *, 'Difference in solutions:', diff
    print *, '--------------------------------------------------'

    !Write out dt against computational time
    open(unit=out_a, file='shm_leapfrog_comptime_v_dt.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file shm_leapfrog_dt_v_comptime.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, finish-start
    if(ios.ne.0) stop 'Error writing to file shm_leapfrog_dt_v_comptime.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_leapfrog_dt_v_comptime.dat.dat'

    !Write out dt against difference in solutions
    open(unit=out_a, file='shm_leapfrog_dt_v_diff.dat', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error opening file shm_leapfrog_dt_v_diff.dat.dat'
    write(unit=out_a, fmt=*, iostat=ios) dt, diff
    if(ios.ne.0) stop 'Error opening file shm_leapfrog_dt_v_diff.dat.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file shm_leapfrog_dt_v_diff.dat.dat'

end program
