!DLA Model Version FINAL

program dla
    implicit none
    integer, parameter :: dp = selected_real_kind(15,300)
    real(kind=dp), parameter :: pi = 3.14159265359_dp

    !lattice variables
    integer :: L, S, i
    integer, allocatable, dimension(:,:) :: lattice
    integer, dimension(2) :: origin, location

    ! cartesian variables
    integer :: x, y

    !polar variables
    real(kind=dp) :: r, theta
    real(kind=dp) :: r_start, r_outer, r_cluster, r_max

    !Logical operators
    logical :: found_cluster = .false. !Is particle adjacent to cluster?
    logical :: exit_circle = .false. !Is particle at max radius?
    logical :: max_bounds = .false. !Is r_outer at max bounds?

    !Mean variables
    real(kind=dp) :: sum_r2, mean_R, log_S, log_R

    !Open files for writing values
    open(file='sv.dat',unit=11)
    open(file='log_sv.dat',unit=12)

    print *, 'Simulation initiated..'

    !Set lattice size and initialise the grid
    L = 800 !Lattice size
    allocate(lattice(-L/2 +1:L/2, -L/2 +1:L/2))
    lattice = 0
    !Set maximum radius of outer circle
    r_max = L/2

    !Set nucleus at origin of lattice
    origin = 0
    lattice(origin(1),origin(2)) = 1
    !Initialise the size of the cluster
    S = 1
    !Initialise radius of the cluster
    r_cluster = 1.0_dp
    sum_r2 = r_cluster**2
    !Set inner circle and outer circle sizes
    r_start = r_cluster + 5.0_dp !Radius to generate particle at
    r_outer = r_start*1.5_dp !Maximum radius particle can drift

    !Start the simulation
    generate_loop : do while (max_bounds .eqv. .false.)
        !Check if the loop condition is valid
        call check_bounds()

        !Generate a particle at a point on the circle
        call new_particle(r_start)

        i = 1
        walk_loop : do
            !print *, 'step:', i
            !Walk the particle
            call random_walk()
            call check_particle()
            !If particle is outside r_outer generate a new one
            if (exit_circle .eqv. .true.) then
            exit_circle = .false.
                cycle generate_loop
            !If particle is touching cluster then grow the cluster and generate a new particle
            else if (found_cluster .eqv. .true.) then
                !Occupy location in lattice
                lattice(location(1),location(2)) = 1
                !Update the cluster
                call cluster_growth()
                !Reset logical operator
                found_cluster = .false.
                cycle generate_loop
            end if
            i = i + 1
            !if (i >= 1000) cycle generate_loop
        end do walk_loop
    end do generate_loop

print *, 'Max Bounds:', max_bounds
print *, 'Number of Particles in cluster:', S
call cluster_plot()

contains
function uniform_rng()
!This function generates a uniform random number in the range of 0 to 1 with a type of real

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
subroutine new_particle(r_start)
!This subroutine generates a new particle at a random point on a circle
!The circle has the radius r_start and a random angle theta
!Theta is generated by calling a uniform random number generator which generates a real number in the range of 0,1, this is then scaled by a factor of 2pi
!The position of the particle on the circle is converted into a location on the underlying lattice

    !Radius of circle to generate particle on
    real(kind=dp), intent(in) :: r_start

    !Update current radius
    r = r_start

    !Generate a random value of theta (in radians)
    theta = uniform_rng()*(2.0_dp*pi)

    !Convert the current polar coordinates into cartesian coordinates
    x = nint(r*cos(theta))
    y = nint(r*sin(theta))

    !Convert the current cartesian coordinates into a location on the lattice map
    location(1) = x
    location(2) = y

    !Print position of new particle to console
    !print *, 'New particle generated'
    !print *, 'New particle location: ', 'x=', location(1), 'y=', location(2)
    !print *, 'New particle Polar position: ', 'radius=', r, 'theta=', theta

end subroutine new_particle
subroutine random_walk()
!This subroutine simulates a random walk for a particle.
!The particle randomly moves in one of the four cardinal directions
!After taking a step, the location of the particle in the lattice is updated and the coordinates subroutine is called to calculate the new values of r and theta

    !Determines direction to take step
    integer :: direction
    !Direction tags
    integer, parameter :: up = 0, down = 1, left = 2, right = 3

    !Generate a random number between 0 and 3 (int function rounds down)
    direction = int(uniform_rng()*4)

    !Simulate the random walk
    if (direction == up) then
        !Move up by one step
        location(2) = location(2) + 1
        !print *, 'Particle has moved up one step'
        !Update coordinates
        call coordinates(location)

    else if (direction == down) then
        !Move down by one step
        location(2) = location(2) - 1
        !print *, 'Particle has moved down one step'
        !Update coordinates
        call coordinates(location)

    else if (direction == left) then
        !Move left by one step
        location(1) = location(1) - 1
        !print *, 'Particle has moved left one step'
        !update coordinates
        call coordinates(location)

    else if (direction == right) then
        !Move right by one step
        location(1) = location(1) + 1
        !print *, 'Particle has moved right one step'
        !Update coordinates
        call coordinates(location)
    end if


end subroutine random_walk
subroutine coordinates(location)
!This subroutine converts the current location of the particle on the lattice into cartesian and polar coordinates

    integer, dimension(2), intent(in) :: location
    !Convert to cartesian coordinates
    x = location(1)
    y = location(2)
    !Convert to polar coordinates
    r = sqrt(real(x,kind=dp)**2 + real(y,kind=dp)**2) !Current radius
    theta = atan2(real(y,kind=dp),real(x,kind=dp)) !Current value of theta

    !print *, 'Current Location: ', 'x=', location(1), 'y=', location(2)
    !print *, 'Current Polar coordinates: ', 'r=',r, 'theta=',theta
end subroutine coordinates
subroutine check_particle()
!This subroutine checks the current location of the particle
!If the distance of the particle from the nucleus is equal to or greater than the maximum radius then the particle is lost
!If the particle touches the cluster (is adjacent to it on the grid), it sticks to its current position and the grid location becomes occupied. The size of the cluster is incremented, and the radius of the cluster is recalculated.


    !Check if particle is outside r_outer

    !If particle is outside the outer circle it is lost
    if (nint(r) >= nint(r_outer)) then
        exit_circle = .true.
        !print *, 'exit_circle:', exit_circle
    end if

    !Check if particle is adjacent to the cluster
    !Above
    if (lattice(location(1),location(2) + 1) == 1) found_cluster = .true.
    !Below
    if (lattice(location(1),location(2) - 1) == 1) found_cluster = .true.
    !Left
    if (lattice(location(1) - 1,location(2)) == 1) found_cluster = .true.
    !Right
    if (lattice(location(1) + 1,location(2)) == 1) found_cluster = .true.


    !print *, 'found_cluster:', found_cluster


end subroutine check_particle
subroutine cluster_growth()
!This subroutine deals with the growth of the cluster and its conditions
!If a particle is added to the cluster then the cluster size must increase and the max distance r_cluster  of any particle in the cluster to the centre of the lattice must be recalculated

    !Check if particle has been added to the cluster
    if (found_cluster .eqv. .true.) then
        S = S + 1
        !print *, 'cluster_size:', S
        if (r > r_cluster) r_cluster = r
            !print *, 'r_cluster:', r_cluster
        sum_r2 = sum_r2 + abs(r)**2
        mean_R = sqrt(1/real(S,kind=dp) * sum_r2)
        log_R = log(mean_R)
        log_S = log(real(S,kind=dp))
        write (unit=11, fmt=*) mean_R, S
        write (unit=12, fmt=*) log_R, log_S
        call check_circles()
    end if

end subroutine cluster_growth
subroutine check_bounds()
!This subroutine checks if the outer radius is on the edge of the grid

    !Check if outer_radius is on edge of grid
    if (nint(r_outer) >= int(r_max-1)) max_bounds = .true.
    !print *, 'max_bounds:', max_bounds

end subroutine check_bounds
subroutine check_circles()
!This subroutine handles the scaling of the particle generation circle and the outer circle
!If the cluster grows to some factor of the particle generation circle, then the size of r_start must increase
!If r_start increases, so must the size of the outer circle

    !Increase size of r_start in relation to cluster
    if (r_cluster >= r_start - 5.0_dp) then
        r_start = r_start + 5.0_dp
        !increase r_outer in relation to r_start to preserve ratio
        r_outer = r_start*1.5_dp
        call check_bounds()
    end if

    !print *, 'r_start:', r_start
    !print *, 'r_outer:', r_outer

end subroutine check_circles
subroutine cluster_plot()
!This subroutine outputs the cluster as a PGM file 'cluster.pgm'

    integer :: i, j, Nx, Ny, max_greys, out_unit
    Nx = L
    Ny = L
    max_greys = 1

    out_unit=10
    open(file='cluster.pgm',unit=out_unit,status='unknown')
    write(out_unit,11) 'P2' !pgm magic number
    write(out_unit,12) Nx,Ny !width, height
    write(out_unit,13) max_greys !max gray value

    !Test by printing one pixel per line
    !Cycle through array
    do j=-L/2 + 1,L/2
        do i=-L/2 + 1,L/2
            write(out_unit,*) lattice(i,j)
        end do
    end do

    close(unit=out_unit)

    11 format(a2) !write out text in a 2 character wide field
    12 format(i3,1x,i3)
    13 format (i5)

end subroutine cluster_plot
end program dla