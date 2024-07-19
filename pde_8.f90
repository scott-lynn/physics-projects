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

module methods 
    use double_precision
    use IO, only: ios, out_a, out_b
    implicit none 

    contains 

        function uniform_rng()
            !This function generates a uniform random number in the range of 0 to 1 with a type of real
                implicit none
                !Function type
                real(kind=dp) :: uniform_rng
            
                !Local variables
                integer :: A, B, M
                integer :: int_x = 12345 !Initial Seed
            
                !Define magic numbers
                A = 137
                B = 150887
                M = 714025
            
                !Generate uniform random integer
                int_x = mod(A*int_x+B,M)
                !print *, 'random num:', int_x
            
                !Normalise random number for output
                uniform_rng = real(int_x,kind=dp)/real(M,kind=dp)
            
                !Catch for wrap around
                if (uniform_rng<0.0_dp) uniform_rng=uniform_rng+1.0_dp
            
                return
            
        end function uniform_rng
                
        subroutine initialise(integ_select, P, lambda, coeff_matrix)
            !Initialises for either FTCS (case 1) or BTCS (case 2)
            implicit none 
            real(kind=dp), intent(in) :: lambda     !Value of lambda, affects stability
            integer, intent(in) :: P               !Size of P x P square matrix
            integer, intent(in) :: integ_select    !Case select
            real(kind=dp), dimension(:,:), allocatable :: coeff_matrix

            if (allocated(coeff_matrix)) then 
                deallocate(coeff_matrix, stat=ios)
                if(ios.ne.0) stop 'Error deallocating coefficient matrix'
            end if 

            allocate(coeff_matrix(P,P), stat=ios)
            if(ios.ne.0) stop 'Error allocating coefficient matrix'

            !Are we creating matrix A or B?
            select case(integ_select)
                case(1) !FTCS
                    !Initialise coefficient matrix
                    call FTCS_coeff(coeff_matrix, lambda)
                case(2) !BTCS
                    !Initialise coefficient matrix
                    call BTCS_coeff(coeff_matrix, lambda)
                case(3) !FTCS Dirichlet Boundary Conditions
                    ! Initialise coefficient matrix
                    call FTCS_coeff_fixed(coeff_matrix,lambda)
                case(4) !BTCS Dirichlet Boundary Conditions 
                    !Initialise coefficient matrix
                    call BTCS_coeff_fixed(coeff_matrix,lambda)
            end select 

        end subroutine initialise 

        subroutine FTCS_coeff(coeff_matrix, lambda)
            !Initialises the FTCS coefficient matrix for a given value of lambda
            implicit none 
            real(kind=dp), dimension(:,:) :: coeff_matrix 
            real(kind=dp) :: lambda 
            integer :: i, P 

            !Check coeff_matrix is square
            if (size(coeff_matrix,1) .ne. size(coeff_matrix,2)) then
                stop 'coeff_matrix is not square'
            end if 
            
            !Get size
            P = size(coeff_matrix,1)

            coeff_matrix = 0.0_dp !Initialise matrix values to zero 
                        
            occupy: do i=1,(P-1),1
                coeff_matrix(i,i) = 1.0_dp - 2.0_dp * lambda    !Diagonals
                coeff_matrix(i,i+1) = lambda                    !Slot to the right
                coeff_matrix(i+1,i) = lambda                    !Slot below
            end do occupy

            !Define final slot (P,P)
            coeff_matrix(P,P) = 1.0_dp - 2.0_dp * lambda 

            !Define periodic boundary conditions 
            coeff_matrix(1,P) = lambda 
            coeff_matrix(P,1) = lambda 
        end subroutine FTCS_coeff 

        subroutine BTCS_coeff(coeff_matrix, lambda)
            implicit none
            real(kind=dp), dimension(:,:) :: coeff_matrix 
            real(kind=dp) :: lambda 
            integer :: i, P 

            ! Check coeff is square
            if (size(coeff_matrix,1) .ne. size(coeff_matrix,2)) then
                stop 'coeff_matrix is not square'
            end if 
            
            P = size(coeff_matrix,1)

            coeff_matrix = 0.0_dp 

            !Assign diagonals 
            occupy: do i=1,(P-1),1
                coeff_matrix(i,i) = 1.0_dp + (2.0_dp * lambda)  !Diagonals
                coeff_matrix(i,i+1) = -lambda                   !Slot to the right
                coeff_matrix(i+1,i) = -lambda                   !Slot below
            end do occupy

            !Define final slot (P,P)
            coeff_matrix(P,P) = 1.0_dp + 2.0_dp * lambda 

            !Define periodic boundary conditions 
            coeff_matrix(1,P) = -lambda 
            coeff_matrix(P,1) = -lambda 

            !Invert the matrix using LAPACK
            call invert_matrix(coeff_matrix)
        end subroutine BTCS_coeff 

        subroutine FTCS_coeff_fixed(coeff_matrix, lambda)
            !Initialises the FTCS coefficient matrix for a given value of lambda with Dirichlet (fixed) boundary conditions
            implicit none 
            real(kind=dp), dimension(:,:) :: coeff_matrix 
            real(kind=dp) :: lambda 
            integer :: i, P 

            !Check coeff_matrix is square
            if (size(coeff_matrix,1) .ne. size(coeff_matrix,2)) then
                stop 'coeff_matrix is not square'
            end if 
            
            !Get size
            P = size(coeff_matrix,1)

            coeff_matrix = 0.0_dp !Initialise matrix values to zero 
                        
            occupy: do i=1,(P-1),1
                coeff_matrix(i,i) = 1.0_dp - 2.0_dp * lambda    !Diagonals
                coeff_matrix(i,i+1) = lambda                    !Slot to the right
                coeff_matrix(i+1,i) = lambda                    !Slot below
            end do occupy

            !Define final slot (P,P)
            coeff_matrix(P,P) = 1.0_dp - 2.0_dp * lambda 
        end subroutine FTCS_coeff_fixed

        subroutine BTCS_coeff_fixed(coeff_matrix, lambda)
            implicit none
            real(kind=dp), dimension(:,:) :: coeff_matrix 
            real(kind=dp) :: lambda 
            integer :: i, P 

            ! Check coeff is square
            if (size(coeff_matrix,1) .ne. size(coeff_matrix,2)) then
                stop 'coeff_matrix is not square'
            end if 
            
            P = size(coeff_matrix,1)

            coeff_matrix = 0.0_dp 

            !Assign diagonals 
            occupy: do i=1,(P-1),1
                coeff_matrix(i,i) = 1.0_dp + (2.0_dp * lambda)  !Diagonals
                coeff_matrix(i,i+1) = -lambda                   !Slot to the right
                coeff_matrix(i+1,i) = -lambda                   !Slot below
            end do occupy

            !Define final slot (P,P)
            coeff_matrix(P,P) = 1.0_dp + 2.0_dp * lambda 

            !Invert the matrix using LAPACK
            call invert_matrix(coeff_matrix)
        end subroutine BTCS_coeff_fixed

        subroutine initialise_heat(P, u_t)
            !Initialises the heat array values to a random series of numbers
            implicit none 
            integer, intent(in) :: P 
            integer :: i 
            real(kind=dp), dimension(:,:), allocatable :: u_t 

            !Check for allocation of u_t
            if (allocated(u_t)) then 
                deallocate(u_t, stat=ios)
                if(ios.ne.0) stop 'Error deallocating u_t matrix'
            end if 

            !Allocate u_t
            allocate(u_t(P,1), stat=ios)
            if(ios.ne.0) stop 'Error allocating u_t matrix'

            !Initialise the heat values 
            randomize: do i=1,P,1 
                 u_t(i,1) =  uniform_rng() * 10.0_dp
            end do randomize

        end subroutine initialise_heat 

        subroutine integrator(integ_select, dx, lambda, coeff_matrix, u_t, total_heat, bound_cond)
            implicit none 
            integer, intent(in) :: integ_select 
            real(kind=dp), intent(in) :: dx
            real(kind=dp), dimension(:,:), allocatable :: coeff_matrix, u_t
            real(kind=dp), dimension(:), allocatable :: total_heat 
            real(kind=dp), dimension(2), intent(in), optional :: bound_cond
            real(kind=dp), dimension(:,:), allocatable :: c_matrix
            real(kind=dp) :: VL, VR, lambda 
            integer :: i 
            integer :: P
            integer :: N
            real(kind=dp) :: t 

            !Check coeff_matrix is square
            if (size(coeff_matrix,1) .ne. size(coeff_matrix,2)) then
                stop 'coeff_matrix is not square'
            end if 

            !Get size
            P = size(coeff_matrix,1)

            !Get no of iterations 
            N = size(u_t,2)

            !Check if fixed boundary conditions (Dirichlet) are present
            if (present(bound_cond)) then 
                !Assign left (VL) and right (VR) boundary conditions
                VL = bound_cond(1)
                VR = bound_cond(2)

                !Check for allocation
                if (allocated(c_matrix)) then 
                    deallocate(c_matrix, stat=ios)
                    if(ios.ne.0) stop 'Error deallocating c_matrix'
                end if 
                
                !Allocate size of c_matrix
                allocate(c_matrix(P,1), stat=ios)
                if(ios.ne.0) stop 'Error allocating c_matrix'

                !Initialise to zero
                c_matrix = 0.0_dp 

                !Apply boundary conditions
                c_matrix(1,1) = VL 
                c_matrix(P,1) = VR 

            end if

            !Check allocation of heat array
            if (allocated(total_heat)) then 
                deallocate(total_heat, stat=ios)
                if(ios.ne.0) stop 'Error deallocating total_heat matrix'   
            end if 

            !Allocate size of total heat array 
            allocate(total_heat(N), stat=ios)
            if(ios.ne.0) stop 'Error allocating total_heat matrix'
            
            !Initialise total_heat
            total_heat = 0.0_dp

            !Initialise time 
            t = 0.0_dp

            !Determine case 
            select case(integ_select)
            case(1) !FTCS PBC

                !Integrate over length of rod 
                do i=1,(N-1),1
                    !Calculate sum of heat for time t 
                    total_heat(i) = sum(u_t(:,i)) * dx

                    !Integrate using FTCS 
                    u_t(:,i+1) = matmul(coeff_matrix, u_t(:,i))

                    ! !Next timestep 
                    ! t = t + dt 
                end do
                
                !Calculate total heat at P
                total_heat(N) = sum(u_t(:,N)) * dx
            
            case(2) !BTCS PBC

                !Integrate over length of rod 
                do i=1,(N-1),1
                    !Calculate sum of heat for time t 
                    total_heat(i) = sum(u_t(:,i)) * dx

                    !Integrate using BTCS
                    u_t(:,i+1) = matmul(coeff_matrix, u_t(:,i))

                    ! !Next timestep 
                    ! t = t + dt 
                end do 
                
                !Calculate total heat at P
                total_heat(N) = sum(u_t(:,N)) * dx 

            case(3) !FTCS Dirichlet fixed boundary conditions

                !Integrate over length of rod 
                do i=1,(N-1),1
                    !Calculate sum of heat for time t 
                    total_heat(i) = sum(u_t(:,i)) * dx 

                    !Integrate using FTCS 
                    u_t(:,i+1) = matmul(coeff_matrix, (u_t(:,i))) + (lambda * c_matrix(:,1))

                    ! !Next timestep 
                    ! t = t + dt 
                end do
                
                !Calculate total heat at N
                total_heat(N) = sum(u_t(:,N)) * dx 
            
            case(4) !BTCS Dirichlet fixed boundary conditions 

                !Integrate over length of rod 
                do i=1,(N-1),1
                    !Calculate sum of heat for time t 
                    total_heat(i) = sum(u_t(:,i)) * dx 
                    !Integrate using BTCS 
                    u_t(:,i+1) = matmul(coeff_matrix, (u_t(:,i) + (lambda * c_matrix(:,1))))

                    ! !Next timestep 
                    ! t = t + dt 
                end do 
                
                !Calculate total heat at N
                total_heat(N) = sum(u_t(:,N)) * dx 

            end select 
        end subroutine integrator

        subroutine invert_matrix(matrix)
            ! Invert the supplied matrix using LAPACK
            implicit none
            real(kind=dp), dimension(:,:), intent(inout) :: matrix
            integer :: N, LWORK, IERR
            integer, dimension(:), allocatable :: IPIV
            real(kind=dp), dimension(:), allocatable :: WORK
            if (size(matrix,1) /= size(matrix,2)) STOP "Matrix is not square"
            N = size(matrix,1)
            allocate(IPIV(N),stat=IERR)
            if (IERR/=0) STOP "Failed to allocate IPIV"
            LWORK = N**2
            allocate(WORK(LWORK),stat=IERR)
            if (IERR/=0) STOP "Failed to allocate WORK"
            call dgetrf(N,N,matrix,N,IPIV,IERR)
            if (IERR/=0) STOP "Error in dgetrf: Matrix is singular"
            call dgetri(N,matrix,N,IPIV,WORK,LWORK,IERR)
            if (IERR/=0) STOP "Error in dgetri: Matrix is singular"
        end subroutine
                
end module methods 

!*******************************************************************************
!Main program
program pde 
!Declare modules
use double_precision
use IO, only: out_a, out_b, ios 
use methods, only: uniform_rng, initialise, integrator, initialise_heat 

implicit none

    integer :: P
    integer :: N
    integer :: i
    integer :: integ_select 
    real(kind=dp) :: L
    real(kind=dp) :: alpha 
    real(kind=dp) :: lambda
    real(kind=dp) :: t  
    real(kind=dp) :: x 
    real(kind=dp) :: dt, dx 
    real(kind=dp) :: ddt 
    real(kind=dp) :: VL, VR 
    real(kind=dp), dimension(:,:), allocatable :: matrix_A, matrix_B 
    real(kind=dp), dimension(:,:), allocatable :: u_t1, u_t2, u_t_initial
    real(kind=dp), dimension(:), allocatable :: total_heat1, total_heat2
    real(kind=dp), dimension(:), allocatable :: FTCS_diff, BTCS_diff 
    real(kind=dp), dimension(2) :: bound_cond 
    !*****************************************************************************
    !Part 1: Solve the heat equation for a 1D ring with periodic boundary conditions

    P = 100                             !No of points along ring of length l
    N = 100                             !No of iterations
    L = 1.0_dp                          !Length of ring in m
    alpha = 1.0e-4_dp                   !Thermal conductivity
    t = 0.0_dp                          !Initial time 
    dt = 0.1_dp                         !Time step
    dx = L/100.0_dp                     !Change in position
    lambda = (alpha * dt)/((dx)**2)     !Calculate lambda for use in coefficient matrix

    !*****************************************************************************
    !Part 1b: FTCS 
    integ_select = 1 !FTCS selected

    !Initialise the coeff matrix for FTCS
    call initialise(integ_select, P, lambda, matrix_A)

    !Randomise the initial heat values at time t0
    call initialise_heat(P, u_t_initial)

    !Assign initial heat values to matrix
    allocate(u_t1(P,N))
    u_t1(:,1) = u_t_initial(:,1)

    !Integrate using matrix multiplication FTCS method
    call integrator(integ_select, dx, lambda, matrix_A, u_t1, total_heat1)

    open(unit=out_a, file='FTCS1.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file FTCS1.dat'

    !Write out total heat against time 
    t = 0.0_dp 
    do i=1,N,1
        write(unit=out_a, fmt=*, iostat=ios) t, total_heat1(i)
        if(ios.ne.0) stop 'Error writing to file FTCS1.dat'
        t = t + dt
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file FTCS1.dat'

    !***************************************************************************
    !Part 1c: BTCS
    integ_select = 2 !BTCS
    !Initialise coefficient matrix for BTCS
    call initialise(integ_select, P, lambda, matrix_B)

    !We want the initial values to be the same as the FTCS test from part 1b
    allocate(u_t2(P,N))
    u_t2(:,1) = u_t_initial(:,1)

    !Integrate using matrix multiplication BTCS method
    call integrator(integ_select, dx, lambda, matrix_B, u_t2, total_heat2)

    open(unit=out_a, file='BTCS1.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file BTCS1.dat'

    !Write out total heat against time
    t = 0.0_dp 
    do i=1,N,1
        write(unit=out_a, fmt=*, iostat=ios) t, total_heat2(i)
        if(ios.ne.0) stop 'Error writing to file BTCS1.dat'
        t = t + dt 
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file BTCS1.dat'

    !***********************************************************************
    !Part 1d: Conservation of heat with varying time step dt

    dt = 0.1_dp                         !Initial time step value
    ddt = 0.01_dp                       !Amount to increment dx by
    lambda = (alpha * dt)/((dx)**2)     !Initial value of lambda
    N = 100

    !Allocate arrays for difference values
    allocate(FTCS_diff(N))
    allocate(BTCS_diff(N))

    !Open files
    open(unit=out_a, file='FTCS_diff.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening FTCS_diff.dat'
    open(unit=out_b, file='BTCS_diff.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening BTCS_diff.dat'

    !Calculate difference in heat conservation for FTCS and BTCS methods for a range of dt values
    do i=1,N - nint(dt/ddt) + 1,1

        !Reset heat arrays 
        deallocate(u_t1, u_t2, total_heat1, total_heat2)

        !FTCS ****************************************
        integ_select = 1 !FTCS selected

        !Initialise
        call initialise(integ_select, P, lambda, matrix_A)
        allocate(u_t1(P,N))
        u_t1(:,1) = u_t_initial(:,1)

        !Integrate using FTCS
        call integrator(integ_select, dx, lambda, matrix_A, u_t1, total_heat1)

        !How well is heat conserved throughout the ring? Quantify heat conservation
        FTCS_diff = total_heat1(N) - total_heat1(1)

        !Write out FTCS against dt
        write(unit=out_a, fmt=*, iostat=ios) dt, FTCS_diff(i)
        if(ios.ne.0) stop 'Error writing to file FTCS_diff.dat'

        !BTCS ******************************************
        integ_select = 2 !BTCS selected

        !Initialise
        call initialise(integ_select, P, lambda, matrix_B)

        !We want the initial values to be the same as the FTCS test from part 1b
        allocate(u_t2(P,N))
        u_t2(:,1) = u_t_initial(:,1)

        !Integrate using BTCS
        call integrator(integ_select, dx, lambda, matrix_B, u_t2, total_heat2)

        !How well is heat conserved throught the ring? Quantify heat conservation
        BTCS_diff(i) = total_heat2(N) - total_heat2(1)

        !Write out BTCS against dt
        write(unit=out_b, fmt=*, iostat=ios) dt, BTCS_diff(i)
        if(ios.ne.0) stop 'Error writing to file BTCS_diff.dat'

        dt = dt + ddt                       !Increment dt
        lambda = (alpha * dt)/((dx)**2)     !Recalculate lambda
    end do 

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error closing file FTCS_diff.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error closing file BTCS_diff.dat'

    

!**********************************************************************************
!Part 2 Heat equation equilibrium and Dirichlet boundary conditions
!Part 2b Dirichlet Boundary Conditions
    
    L = 1.0_dp                          !Length of bar in m
    P = 100                             !No of points along bar of length l
    N = 100000                          !No of iterations
    alpha = 1.0e-3_dp                   !Thermal conductivity
    t = 0.0_dp                          !Initial time 
    x = 0.0_dp                          !Initial position
    dt = 0.49_dp                         !Time step
    dx = L/100.0_dp                     !Change in position
    lambda = (alpha * dt)/((dx)**2)     !Calculate lambda for use in coefficient matrix
    VL = 0.0_dp                         !Left boundary condition
    VR = 100.0_dp                       !Right boundary condition

    !Set up boundary condition array for passing to subroutines
    bound_cond(1) = VL 
    bound_cond(2) = VR 

    !Reset heat arrays 
    deallocate(u_t1, u_t2, total_heat1, total_heat2)
    print *, 'lambda:', lambda 
    print *, 'dt:', dt 

    !********************************************************
    !FTCS
    integ_select = 3 !FTCS with Dirichlet conditions selected

    !Initialise the coeff matrix for FTCS
    call initialise(integ_select, P, lambda, matrix_A)

    !Assign initial heat values to matrix
    allocate(u_t1(P,N))
    !Multiply to get initial values between 0 and 100
    u_t_initial = u_t_initial * 10.0_dp 
    u_t1(:,1) = u_t_initial(:,1)

    !Integrate using matrix multiplication FTCS method
    call integrator(integ_select, dx, lambda, matrix_A, u_t1, total_heat1, bound_cond)

    open(unit=out_a, file='FTCS2.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file FTCS2.dat'

    !Write out total heat against time 
    t = 0.0_dp 
    do i=1,N,1
        write(unit=out_a, fmt=*, iostat=ios) t, total_heat1(i)
        if(ios.ne.0) stop 'Error writing to file FTCS2.dat'
        t = t + dt
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file FTCS2.dat'

    open(unit=out_b, file='dirichlet_FTCS.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening file dirichlet_FTCS.dat'

    !Write out equilibrium solution for each point on the bar  
    do i=1,P,1
        write(unit=out_b, fmt=*, iostat=ios) x, u_t1(i,N)
        if(ios.ne.0) stop 'Error writing to file dirichlet_FTCS.dat'
        x = x + dx 
    end do 

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file dirichlet_FTCS.dat'

    !********************************************************
    !BTCS
    integ_select = 4 !BTCS with Dirichlet conditions selected

    !Initialise coefficient matrix for BTCS
    call initialise(integ_select, P, lambda, matrix_B)

    !We want the initial values to be the same as the FTCS
    allocate(u_t2(P,N))
    u_t2(:,1) = u_t_initial(:,1)

    !Integrate using matrix multiplication BTCS method
    call integrator(integ_select, dx, lambda, matrix_B, u_t2, total_heat2, bound_cond)

    open(unit=out_a, file='BTCS2.dat', iostat=ios)
    if(ios.ne.0) stop 'Error in opening file BTCS2.dat'

    !Write out total heat against time
    t = 0.0_dp 
    do i=1,N-1,1
        ! write(unit=out_a, fmt=*, iostat=ios) t, total_heat2(i)
        write(unit=out_a, fmt=*, iostat=ios) t, total_heat2(i)
        if(ios.ne.0) stop 'Error writing to file BTCS2.dat'
        t = t + dt 
    end do 

    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file BTCS2.dat'

    open(unit=out_b, file='dirichlet_BTCS.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening file dirichlet_BTCS.dat'

    
    !Write out equilibrium solution for each point on the bar  
    x = 0.0_dp 
    do i=1,P,1
        write(unit=out_b, fmt=*, iostat=ios) x, u_t2(i,N)
        if(ios.ne.0) stop 'Error writing to file dirichlet_BTCS.dat'
        x = x + dx 
    end do 

    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file dirichlet_BTCS.dat'
    

end program pde