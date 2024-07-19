module double_precision
    implicit none 
    integer, parameter :: dp = selected_real_kind(15,300)
end module double_precision

module constants 
    !Contains constants used within Physics.
    use double_precision
    implicit none 
        
    !Constants
    real(kind=dp), parameter :: pi = 4.D0*DATAN(1.D0)       !Pi (Standard Fortran Implementation)
    real(kind=dp), parameter :: eV = 1.602176634e-19        !Electron Volt [J]
    real(kind=dp), parameter :: h = 6.62607015e-4_dp        !Planck's constant [J/Hz ]
    real(kind=dp), parameter :: h_bar = h/(2.0_dp*pi)       !Reduced Planck's constant [J]
    real(kind=dp), parameter :: e_mass = 9.1093837015e-31   !Electron mass [kg]
    real(kind=dp), parameter :: c = 299792458               !Speed of Light [m/s]
    
end module constants

module IO 
    !Module which contains I/O paths and the iostat error checking variable
    implicit none 
    integer :: ios  !Error checking

    !IO channels    
    integer, parameter :: out_a = 11                
    integer, parameter :: out_b = 12    
    integer, parameter :: out_c = 13
    integer, parameter :: out_d = 14 
    integer, parameter :: out_e = 15
end module IO

module global
    !Defines global parameters so that they may be passed to functions/subroutines without requiring an argument
    use double_precision 
    implicit none 

    !Potential Well values
    real(kind=dp) :: L, L_min, L_max 
    real(kind=dp) :: V_inf, V_none, V0 
    integer :: well_type     

    !Energy value
    real(kind=dp) :: E 
    
end module global

module equations 
    !Contains equations required for the experiment
    use double_precision 
    use global
    implicit none 

    contains 

        function potential(x)
            !Calculates the potential of the well as a function of x
            !Makes use of the global variable, well_type which defines if the program should model the potential for the infinite or finite well
            !well_type = 1 -> infinite square well
            !well_type = 2 -> finite square well
            real(kind=dp) :: potential
            real(kind=dp), intent(in) :: x 

            select case(well_type)
            case(1)
                if (abs(x) < (L/2.0_dp)) then
                    potential = 0.0_dp
                else 
                    potential = 1000.0_dp
                end if

            case(2)
                if (abs(x) < (L/2.0_dp)) then 
                    potential = V_none 
                else 
                    potential = V0
                end if 
            end select 
        end function potential

        function dydx(x, y, z) 
            !Returns psi when integrated, first part of the coupled ODEs for the RK4 method
            real(kind=dp) :: dydx
            real(kind=dp), intent(in) :: x, y, z

            !Gives psi when integrated
            dydx = z

        end function dydx

        function dzdx(x, y, z)
            !Returns dpsi/dx when integrated, second part of the coupled ODEs for the RK4 method
            implicit none 
            real(kind=dp) :: dzdx
            real(kind=dp), intent(in) :: x, y, z

            !Gives dpsi when integerated
            dzdx = (potential(x) - E) * y

        end function dzdx

        function k_squared(x) 
            !Used for the Numerov method so that the TISE SODE can be integrated in the form of two functions f(x)y(x), this is the f(x) part
            real(kind=dp) :: k_squared 
            real(kind=dp), intent(in) :: x 

            !returns the kinetic energy value of the TISE
            k_squared = potential(x) - E
        end function k_squared

    end module equations 

module methods 
    !Contains numerical methods
    use double_precision 
    use equations 
    implicit none 

    contains 

        subroutine rk4_coupled(dx, x, y, z, fy, fz, y_next, z_next)
            !Contains the RK4 method for a SODE split into two coupled first order ODEs.

            !dx is the step size input, this can be inputted as positive to step forwards or negative to step backwards
            !x is the current value of the independent variable, required as an input
            !y is the solution to the SODE at y(x)
            !z is the first order derivative integrated from the SODE at z(x)
            !fy and fz are the coupled input ODEs to be integrated 
            !y_next and z_next are the output variables at y(x+dx) and z(x+dx)

            real(kind=dp), intent(in) :: dx                     !Step size
            real(kind=dp), intent(in) :: x                      !Initial value of independent variable
            real(kind=dp), intent(in) :: y, z                   !Input values of y and z
            real(kind=dp) :: fy, fz                             !Input functions, fy = dy/dx, fz = dz/dx
            real(kind=dp) :: y_next, z_next                     !Output variables: y_next = yn+1, z_next=zn+1
            real(kind=dp) :: k1y, k2y, k3y, k4y 
            real(kind=dp) :: k1z, k2z, k3z, k4z
            
            !Calculate values at the beginning of the slope
            k1y = fy(x, y, z) * dx
            k1z = fz(x, y, z) * dx 

            !Calculate values at the midpoint of the slope
            k2y = fy(x + dx/2.0_dp, y + k1y/2.0_dp, z + k1z/2.0_dp) * dx 
            k2z = fz(x + dx/2.0_dp, y + k1y/2.0_dp, z + k1z/2.0_dp) * dx 

            !Calculate values at the midpoint of the slope using the previous step
            k3y = fy(x + dx/2.0_dp, y + k2y/2.0_dp, z + k2z/2.0_dp) * dx 
            k3z = fz(x + dx/2.0_dp, y + k2y/2.0_dp, z + k2z/2.0_dp) * dx 

            !Calculate values at the end of the slope
            k4y = fy(x + dx, y + k3y, z + k3z) * dx 
            k4z = fz(x + dx, y + k3y, z + k3z) * dx 

            !Calculate yn+1 and zn+1 values at point x+dx 
            y_next = y + k1y/6.0_dp + k2y/3.0_dp + k3y/3.0_dp + k4y/6.0_dp 
            z_next = z + k1z/6.0_dp + k2z/3.0_dp + k3z/3.0_dp + k4z/6.0_dp
             
        end subroutine rk4_coupled

        function G(h, x, y, f)
            !Calculates G values for the Numerov method

            !h is the value of the step size
            !x is the current value of the indepdendent variable 
            !y is the current value of the solution, the y(x) part in the numerov method
            !f is an input function of the form f(x) required for the numerov method

            real(kind=dp) :: G
            real(kind=dp), intent(in) :: x, y
            real(kind=dp) :: f 
            real(kind=dp) :: h

            !Calculate G at point x
            G = (1.0_dp - (h**2/12.0_dp) * f(x)) * y
        end function G

        subroutine numerov(h, x, y, G_n, G_prev, f, y_next, y_new_prev)
            !A recursive function which calculates the yn+1 value using the Numerov Method 
            
            !h is the value of the step size
            !x is the indepdendent variable 
            !y is the solution to the SODE at y(x)
            !G_n is the current value of G(x)
            !G_prev is the value of G(x-h)
            !f is the function f(x) required for the Numerov method 
            !y_next is the value of the solution to the SODE at yn+1 to be returned
            !y_new_prev is the new value of yn to be returned

            real(kind=dp), intent(inout) :: G_n, G_prev 
            real(kind=dp), intent(in) :: y 
            real(kind=dp) :: G_next, G_new_prev
            real(kind=dp) :: f 
            real(kind=dp) :: h, x
            real(kind=dp) :: y_next, y_new_prev
            
            !Hold values for output as n becomes n-1
            y_new_prev = y
            G_new_prev = G_n

            !Calculate G for the next step n+1
            G_next = 2*G_n - G_prev + h**2 * f(x) * y

            !Gn -> Gn+1
            G_n = G_next 

            !Gn-1 -> Gn
            G_prev = G_new_prev

            !Calculate and return value of yn+1
            y_next = G_next/(1.0_dp - ((h**2)/12.0_dp)*f(x))
        end subroutine numerov 
    end module methods

    module testing 
        !Contains routines for testing areas within the program
        use double_precision
        use IO 
        use global
        use methods 
        use equations 

        contains 
        subroutine well_check(dx, L) 
            !Tests that the infinite/finite square well is modelled correctly, plot the .dat files using a graphing program to check visually
            real(kind=dp), intent(in) :: dx, L
            real(kind=dp) :: V, x

            !Check potential well is modelled correctly for infinite square well
            x = -L
            well_type = 1
            open(unit=out_a, file='potential_test_infinite.dat', iostat=ios)
            if(ios.ne.0) stop 'Error in opening file potential_test_infinite.dat'

            do while (x <= L)
                V = potential(x)
                write(unit=out_a, fmt=*, iostat=ios) x, V 
                if(ios.ne.0) stop 'Error writing to file potential_test_infinite.dat'
                x = x + dx 
            end do 

            close(unit=out_a, iostat=ios)
            if(ios.ne.0) stop 'Error closing file potential_test_infinite.dat'

            !Check potential well is modelled correctly for finite square well 
            x = -L
            well_type = 2
            open(unit=out_b, file='potential_test_finite.dat', iostat=ios)
            if(ios.ne.0) stop 'Error in opening file potential_test_finite.dat'

            do while (x <= L)
                V = potential(x)
                write(unit=out_b, fmt=*, iostat=ios) x, V 
                if(ios.ne.0) stop 'Error writing to file potential_test_finite.dat'
                x = x + dx 
            end do 

            close(unit=out_b, iostat=ios)
            if(ios.ne.0) stop 'Error closing file potential_test_finite.dat'

        end subroutine well_check 
    end module testing 

program schrodinger
!Complab Experiment 2.5. integrates the TISE to determine the energy eigenvalues using the shooting method and DLS.
    use double_precision            !Allows use of real_kind_dp(15,300)
    use constants, only : pi        !Used for calculation of analytical value
    use IO                          !Allows use of I/O paths
    use global                      !Global variables
    use equations                   !Allows use of coupled ODEs
    use methods                     !Allows use of numerical methods
    use testing                     !Testing harness 
    implicit none 
     
    real(kind=dp) :: x, dx, x_prev                              !Position values
    real(kind=dp) :: E_min, E_max, dE, true_E, E_numerical      !Energy Values
    real(kind=dp) :: psi_L, dpsi_L                              !Values for LHS of shooting method
    real(kind=dp) :: psi_R, dpsi_R                              !Values for RHS of shooting method
    real(kind=dp) :: psi_initial, dpsi_initial                  !Initial values for psi and dpsi
    real(kind=dp) :: DLS_L, DLS_R                               !DLS variables for LHS and RHS
    real(kind=dp) :: psi_prev_L, psi_prev_R                     !Holding variables for psi n-1
    real(kind=dp) :: dpsi_prev_L, dpsi_prev_R                   !Holding variables for dpsi n-1
    real(kind=dp) :: G_n, G_prev                                !G values to pass into Numerov Method subroutine

    !Settings for the potential well 
    L = 4.0_dp              !Width of the well
    L_min = -L/2.0_dp       !LHS boundary position
    L_max = L/2.0_dp        !RHS boundary position
    V_none = 0.0_dp         !Potential V=0
    V_inf = 1000.0_dp       !Potential outside of infinite square well
    V0 = 50.0_dp            !Potential outside of finite square well
    
    dx = 1.0e-4_dp          !Step size for position
    dE = 1.0e-4_dp          !Energy step size
    psi_initial = 0.0_dp    !Initial values for psi
    dpsi_initial = 1.0_dp   !Initial values for dpsi

    E_min = 0.5977_dp       !Minimum Energy boundary to scan from
    E_max = 0.5978_dp       !Maximum Energy boundary to stop at
    well_type = 1           !Sets potential well to the infinite square well

    !Test that the potential wells are functioning correctly
    !call well_check(dx, L)

    !-----------------------------------------------------------------------------------------------------------------
    !Infinite potential well - Shooting Method 
    well_type = 1    !Set case to infinite potential well
    E = E_min        !Initial value of trial energy

    !Open files
    open(unit=out_c, file='DLS_L.dat', status='replace', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file DLS_L.dat'
    open(unit=out_d, file='DLS_R.dat', status='replace', position='append', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file DLS_R.dat'

    !Scan from E_min to E_max to determine energy value through a DLS plot
    do while (E <= E_max)

        !Integrate from the left side -----------------------------------
        psi_L = psi_initial
        dpsi_L = dpsi_initial
        x = L_min - 0.2_dp  

        !RK4 loop for LHS
        do while (x <= 0.0_dp)
            call rk4_coupled(dx, x, psi_L, dpsi_L, dydx, dzdx, psi_L, dpsi_L)
            x = x + dx 
        end do

        !integrate from the right side ------------------------------------
        psi_R = psi_initial
        dpsi_R = -dpsi_initial      !Reverse partity to dpsi_L
        x = L_max + 0.2_dp

        !RK4 loop for RHS
        do while (x >= 0.0_dp)
            call rk4_coupled((-dx), x, psi_R, dpsi_R, dydx, dzdx, psi_R, dpsi_R)
            x = x - dx 
        end do 

        !Calculate the value of DLS for both the LHS and RHS at the fixed point x=0
        DLS_L = dpsi_L/psi_L 
        DLS_R = dpsi_R/psi_R

        !Write the DLS values out to file against the Energy value
        write(unit=out_c, fmt=*, iostat=ios) DLS_L, E
        if(ios.ne.0) stop 'Error writing to file DLS_L.dat'
        write(unit=out_d, fmt=*, iostat=ios) DLS_R, E
        if(ios.ne.0) stop 'Error writing to file DLS_R.dat'

        !Iterate to the next value of trial energy
        E = E + dE
    end do 

    close(unit=out_c, iostat=ios)
    if(ios.ne.0) stop 'Error closing file DLS_L.dat'
    close(unit=out_d, iostat=ios)
    if(ios.ne.0) stop 'Error closing file DLS_R.dat'
!-------------------------------------------------------------------------------------------
    !Once a 'ideal' value for E is discovered from the DLS plot, it can be entered here to plot psi and dpsi against x
    true_E = ((1.0_dp**2.0_dp) * (pi ** 2.0_dp))/L**2.0_dp  !Expected analytical solution
    E_numerical = 0.597785                                  !Value of E1 from DLS plot
    E = E_numerical                                         !Set value for plotting psi and dpsi

    !Integrate from the left side -------
    psi_L = psi_initial
    dpsi_L = dpsi_initial
    x = L_min - 0.2_dp

    !Open files
    open(unit=out_a, file='psi_L.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file psi_L.dat'
    open(unit=out_b, file='dpsi_L.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file dpsi_L.dat'

    !RK4 Loop
    do while (x <= 0.0_dp)
        call rk4_coupled(dx, x, psi_L, dpsi_L, dydx, dzdx, psi_L, dpsi_L)

        !Write out values of psi and dpsi against x
        write(unit=out_a, fmt=*, iostat=ios) x, psi_L
        if(ios.ne.0) stop 'Error writing to file psi_L.dat'
        write(unit=out_b, fmt=*, iostat=ios) x, dpsi_L
        if(ios.ne.0) stop 'Error writing to file dpsi_L.dat'

        !Iterate x by dx
        x = x + dx 
    end do

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file psi_L.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file dpsi_L.dat'

    !integrate from the right side ------------------
    psi_R = psi_initial
    dpsi_R = -dpsi_initial      !Opposite parity to dpsi_L
    x = L_max + 0.2_dp

    open(unit=out_a, file='psi_R.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file psi_R.dat'
    open(unit=out_b, file='dpsi_R.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file dpsi_R.dat'

    !RK4 loop
    do while (x >= 0.0_dp)
        call rk4_coupled((-dx), x, psi_R, dpsi_R, dydx, dzdx, psi_R, dpsi_R)

        !Write out psi and dpsi against x
        write(unit=out_a, fmt=*, iostat=ios) x, psi_R
        if(ios.ne.0) stop 'Error writing to file psi_R.dat'
        write(unit=out_b, fmt=*, iostat=ios) x, dpsi_R
        if(ios.ne.0) stop 'Error writing to file dpsi_R.dat'

        !Iterate x by -dx
        x = x - dx 
    end do 
    
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file psi_R.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file dpsi_R.dat'

    !--------------------------------------------------------------------------------------------------
    !NUMEROV METHOD 
    E = E_numerical                     !Set energy value

    !Do left hand side ----------------------------------
    x_prev = L_min - 2.0_dp             !xn-1  
    x = x_prev + dx                     !xn
    psi_prev_L = psi_initial 
    dpsi_prev_L = dpsi_initial  

    open(unit=out_a, file='psi_L_numerov.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file psi_L_numerov.dat'

    !Call RK4 to generate psi n and psi n-1
    call rk4_coupled(dx, x_prev, psi_prev_L, dpsi_prev_L, dydx, dzdx, psi_L, dpsi_L)
    
    write(unit=out_a, fmt=*, iostat=ios) x_prev, psi_prev_L
    if(ios.ne.0) stop 'Error writing to file psi_L_numerov.dat'
    write(unit=out_a, fmt=*, iostat=ios) x, psi_L
    if(ios.ne.0) stop 'Error writing to file psi_L_numerov.dat'
    
    !calculate Gn-1
    G_prev = G(dx, x_prev, psi_prev_L, k_squared)

    !Calculate Gn
    G_n = G(dx, x, psi_L, k_squared)

    !Numerov loop
    do while (x <= 0.0_dp)

        !Call Numerov to calculate G values and psi n+1
        call numerov(dx, x, psi_L, G_n, G_prev, k_squared, psi_L, psi_prev_L)
        
        !Iterate x by dx
        x = x + dx 

        !Write out psi against x
        write(unit=out_a, fmt=*, iostat=ios) x, psi_L
        if(ios.ne.0) stop 'Error writing to file psi_L_numerov.dat'
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file psi_L_numerov.dat'

    !Do right hand side ----------------------------------------------------
    x_prev = L_max + 2.0_dp
    x = x_prev - dx                 !Step to the left towards x=0
    psi_prev_R = psi_initial 
    dpsi_prev_R = -dpsi_initial     !Opposite parity to LHS

    open(unit=out_a, file='psi_R_numerov.dat', status='replace', iostat=ios)
    if(ios .ne. 0) stop 'Error in opening file psi_R_numerov.dat'

    !Call RK4 to generate psi n and psi n-1
    call rk4_coupled((-dx), x_prev, psi_prev_R, dpsi_prev_R, dydx, dzdx, psi_R, dpsi_R)
    
    !Write initial values to file
    write(unit=out_a, fmt=*, iostat=ios) x_prev, psi_prev_R
    if(ios.ne.0) stop 'Error writing to file psi_L_numerov.dat'
    write(unit=out_a, fmt=*, iostat=ios) x, psi_R
    if(ios.ne.0) stop 'Error writing to file psi_L_numerov.dat'
    
    !calculate Gn-1
    G_prev = G((-dx), x_prev, psi_prev_R, k_squared)

    !Calculate Gn
    G_n = G((-dx), x, psi_R, k_squared)

    !Numerov loop
    do while (x >= 0.0_dp)

        !Call Numerov to calculate G values and psi n+1
        call numerov(-dx, x, psi_R, G_n, G_prev, k_squared, psi_R, psi_prev_R)
        
        !Iterate x by -dx
        x = x - dx 
        write(unit=out_a, fmt=*, iostat=ios) x, psi_R
        if(ios.ne.0) stop 'Error writing to file psi_R_numerov.dat'
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file psi_R_numerov.dat'
end program 

