module double_precision
    implicit none 
    integer, parameter :: dp = selected_real_kind(15,300)
end module double_precision

module constants 
    !Contains constants used within Physics.
    use double_precision
    implicit none 
        
    !Constants
    real(kind=dp), parameter :: pi = 4.D0*DATAN(1.D0)           !Pi (Standard Fortran Implementation)
    real(kind=dp), parameter :: two_pi = 2.0_dp * pi            !2pi
    real(kind=dp), parameter :: eV = 1.602176634e-19_dp         !Electron Volt [J]
    real(kind=dp), parameter :: h = 6.62607015e-4_dp            !Planck's constant [J/Hz ]
    real(kind=dp), parameter :: h_bar = h/(two_pi)              !Reduced Planck's constant [J]
    real(kind=dp), parameter :: e_mass = 9.1093837015e-31_dp    !Electron mass [kg]
    real(kind=dp), parameter :: c = 299792458_dp                !Speed of Light [m/s]
    complex(kind=dp), parameter :: i = (0.0_dp,1.0_dp)          !Imaginary number i = sqrt(-1)
    
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

module DFT 
    !Module containing functions/subroutines required to perform a discrete fourier transform (DFT)
    use double_precision
    use IO, only : ios 
    implicit none 

    contains 
        subroutine DFT_coeff(N, W_matrix)
            !Subroutine which takes matrix W_matrix and initialises it as the N,N square coefficient matrix required for the DFT 
            use constants, only : two_pi, i 
            implicit none 

            !Inputs
            integer, intent(in) :: N 
            complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: W_matrix
            
            !Local variables 
            complex(kind=dp) :: W
            integer :: j, k

            !Check if input matrix is already allocated
            if (allocated(W_matrix)) then 
                deallocate(W_matrix, stat=ios)
                if(ios.ne.0) stop 'Error deallocating W'
            end if 

            !Allocate W_matrix to size N,N (square matrix)
            allocate(W_matrix(0:N-1,0:N-1), stat=ios)
            if(ios.ne.0) stop 'Error allocating W'

            !Determine constant 
            W = exp((two_pi * i)/ real(N, kind=dp))

            !Occupy the matrix
            do j = 0, N-1, 1          !Rows
                do k = 0, N-1, 1      !Columns 
                    W_matrix(j,k) = W**(real(j,kind=dp)*real(k,kind=dp))
                end do 
            end do 
        end subroutine DFT_coeff 

    subroutine do_DFT (N, t_data, W_matrix, f_data)
        !Subroutine which performs the Discrete Fourier Transform (DFT) given an input vector t_data of length (N) and the coefficient matrix: W_matrix of size (N,N)
        !Returns a vector f_data of length (N), containing the values of the DFT result
        !The t_data vector represents an input signal in the time domain, the f_data vector represents the corresponding fourier transform of the input signal in the frequency domain
        !The f_data vector is of type complex, and contains a real and imaginary part 

        !Input 
        integer, intent(in) :: N 
        complex(kind=dp), dimension(:), allocatable :: t_data !Time domain array
        complex(kind=dp), dimension(:,:), allocatable :: W_matrix !Coeff matrix

        !Inout 
        complex(kind=dp), dimension(:), allocatable, intent(inout) :: f_data !Frequency domain array


        !Check W_matrix is preallocated 
        if (allocated(W_matrix).neqv..true.) then 
            stop 'W_matrix is not preallocated, allocate W_matrix to size (N,N) before input'
        end if 

        !Check W_matrix is square
        if (size(W_matrix,1) .ne. size(W_matrix,2)) then
            stop 'W_matrix is not square, ensure that W_matrix is of size (N,N)'
        end if 

        !Check that t_data is preallocated 
        if (allocated(t_data).neqv..true.) then 
            stop 't_data is not preallocated, allocate t_data to size (N) and initialise'
        end if 

        !Check if f_data is allocated
        if (allocated(f_data)) then 
            deallocate(f_data, stat = ios)
            if(ios.ne.0) stop 'Error deallocating f_data'
        end if 

        !Allocate f_data 
        allocate(f_data(0:N-1), stat=ios)
        if(ios.ne.0) stop 'Error allocating f_data'

        !Calculate DFT
        f_data = matmul(W_matrix, t_data)
        !Normalise 
        f_data = f_data(:)/real(N,kind=dp)

    end subroutine do_DFT 

    subroutine triangle_window(N, signal)
        !Applies a triangle windowing function to a 1D input signal of size N 

        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal 
        integer, intent(in) :: N 
        integer :: iter 
        real(kind=dp) :: weight

        !Check that signal is allocated
        if (allocated(signal).neqv..true.) then 
            stop 'signal is not preallocated, allocate signal to size (N) and initialise'
        end if 

        !Apply triangle window 
        do iter = 0, N-1, 1
            weight = 1.0_dp - abs(2.0 * (real(iter,kind=dp) - 0.5_dp * (real(N,kind=dp) + 1.0_dp))) / (real(N,kind=dp) + 1.0_dp)
            signal(iter) = signal(iter) * weight 
        end do 

    end subroutine triangle_window 

    subroutine cosine_bell_window(N, signal)
        !Applies a cosine bell windowing function to a 1D input signal of size N 

        use constants, only: two_pi 
        implicit none 

        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal 
        integer, intent(in) :: N 
        integer :: iter 
        real(kind=dp) :: weight, val 

        !Apply cosine bell window
         do iter = 0, N-1, 1
            val =  two_pi * real(iter, kind=dp) / (real(N,kind=dp) -1.0_dp )
            weight = 0.5_dp - 0.5_dp * cos(val)
            signal(iter) = signal(iter) * weight 
        end do 
        

    end subroutine cosine_bell_window

    subroutine gaussian_bell_window(N, signal)
        !Applies a gaussian bell windowing function to a 1D input signal of size N 

        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal 
        integer, intent(in) :: N 
        integer :: iter 
        real(kind=dp) :: weight, sigma 

        !Adjust sigma as required 
        sigma = real(N,kind=dp) / 6.0_dp 

        !Apply gaussian bell window 
        do iter = 0, N-1, 1
            weight = exp(-0.5_dp * (((real(iter,kind=dp) - 0.5_dp * real(N,kind=dp)+1.0_dp))/sigma)**2.0_dp )
            signal(iter) = signal(iter) * weight
        end do

    end subroutine gaussian_bell_window

    subroutine rectangular_pulse(N, signal, amp, L, sample_rate, start_pos)
        
    !---------------------------------------------------------------------------------------------------!
    !                                       RECTANGULAR PULSE                                           !
    !                                                                                                   !
    !       Generates a rectangular pulse starting at t=0, ending at t=L                                !
    !                                                                                                   !
    !---------------------------------------------------------------------------------------------------!

        use IO, only: out_a, ios 
        implicit none 

        !Inputs
        integer, intent(in) :: N                                                    !Number of samples
        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal        !Signal vector
        real(kind=dp), intent(in) :: amp                                            !Amplitude of pulse
        real(kind=dp), intent(in) :: L                                              !Length of pulse
        real(kind=dp), intent(in) :: sample_rate                                    !Sample rate of signal
        real(kind=dp), intent(in) :: start_pos                                      !Start position                             
        
        !Local variables
        integer :: iter                                                             !Iteration loop counter
        real(kind=dp) :: t                                                          !Current time
        real(kind=dp) :: dt                                                         !Timestep size

        !Calculate timestep
        dt = 1.0_dp / sample_rate

        open(unit=out_a, file='rectangle_pulse.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file rectangle_pulse.dat'

        do iter = 0, N-1, 1
            !Calculate time
            t = real(iter,kind=dp) * dt 

            if (t <= L .and. t >= start_pos) then 
                !Generate rectangle pulse
                signal(iter) = cmplx(amp, 0.0_dp, kind=dp)
            else   
                !Zero signal
                signal(iter) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            end if 

            !Write data out to file
            write(unit=out_a, fmt=*, iostat=ios) t, signal(iter)%re 
            if(ios.ne.0) stop 'Error writing to rectangle_pulse.dat'
        end do 

        close(unit=out_a, iostat=ios)
        if(ios.ne.0) stop 'Error writing to file rectangle_pulse.dat'

    end subroutine rectangular_pulse

    subroutine triangle_pulse(N, signal, amp, L, sample_rate, start_pos)
        !Generates a symmetric triangular pulse signal

    !---------------------------------------------------------------------------------------------------!
    !                                       TRIANGLE PULSE                                              !
    !                                                                                                   !
    !       Generates a symmetric triangle pulse starting at t=0, ending at t=L, with its peak at L/2   !                                 !
    !                                                                                                   !
    !---------------------------------------------------------------------------------------------------!
        use IO, only: out_a, ios 
        implicit none 

        !Inputs
        integer, intent(in) :: N                                                    !Number of samples
        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal        !Signal vector
        real(kind=dp), intent(in) :: amp                                            !Amplitude of pulse
        real(kind=dp), intent(in) :: L                                              !Length of pulse
        real(kind=dp), intent(in) :: sample_rate                                    !Sample rate of signal
        real(kind=dp), intent(in) :: start_pos 

        !Local variables
        integer :: iter                                                             !Iteration loop counter
        real(kind=dp) :: t                                                          !Current time
        real(kind=dp) :: dt                                                         !Timestep size
        real(kind=dp) :: peak                                                       !Peak position

        !Calculate peak position
        peak = (start_pos + L )/2.0_dp 

        !Calculate timestep
        dt = 1.0_dp / sample_rate

        open(unit=out_a, file='triangle_pulse.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file triangle_pulse.dat'

        signal = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do iter = 0, N-1, 1
            !Calculate time
            t = real(iter,kind=dp) * dt 
        
            if (t < start_pos) then 
                !Zero signal before start_pos
                signal(iter) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else if (t <= peak) then
                !Linear increase towards peak
                signal(iter) = cmplx(amp * (t - start_pos) / (peak - start_pos), 0.0_dp, kind=dp)
            else if (t <= L) then 
                !Linear decrease from peak to zero
                signal(iter) = cmplx(amp * (L -t) / (L-peak), 0.0_dp, kind=dp)
            else
                !Zero signal after L
                signal(iter) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            end if

            !Write data out to file
            write(unit=out_a, fmt=*, iostat=ios) t, signal(iter)%re 
            if(ios.ne.0) stop 'Error writing to triangle_pulse.dat'
        end do 

        close(unit=out_a, iostat=ios)
        if(ios.ne.0) stop 'Error writing to file triangle_pulse.dat'

    end subroutine triangle_pulse

    subroutine gaussian_pulse(N, signal, amp, L, sigma, sample_rate, start_pos)

    !---------------------------------------------------------------------------------------------------!
    !                                       GAUSSIAN PULSE                                              !        
    !                                                                                                   !
    !       Generates a gaussian pulse with its peak at L/2                                             !
    !       Larger values of sigma may cause the pulse width to become larger than the sample window,   !
    !       it is recommended to keep sigma in the range of 0 <= sigma <= 1.0                           !
    !                                                                                                   !
    !---------------------------------------------------------------------------------------------------!
        
        use IO, only: out_a, ios
        use constants, only: two_pi
        implicit none
        
        !Inputs
        integer, intent(in) :: N                                                !Number of samples
        complex(kind=dp), dimension(:), allocatable, intent(inout) :: signal    !Signal vector
        real(kind=dp), intent(in) :: amp                                        !Amplitude of pulse
        real(kind=dp), intent(in) :: L                                          !Length of pulse
        real(kind=dp), intent(in) :: sigma                                      !Width of pulse / standard deviation of gaussian curve
        real(kind=dp), intent(in) :: sample_rate                                !Sample rate of signal
        real(kind=dp), intent(in) :: start_pos                                  !Start position of signal
        
        !Local variables
        integer :: iter                                                         !Iteration loop counter
        real(kind=dp) :: t                                                      !Current time
        real(kind=dp) :: dt                                                     !Timestep size
        real(kind=dp) :: peak                                                   !Peak of pulse
        
        !Calculate time step
        dt = 1.0_dp / sample_rate

        !Calculate peak of pulse
        peak = (start_pos + L)/ 2.0_dp 
    
        open(unit=out_a, file='gaussian_pulse.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file gaussian_pulse.dat'
        

        do iter = 0, N-1, 1
            !Calculate time
            t = real(iter, kind=dp) * dt
                
            if (t < start_pos) then
                !Before the start position, the signal is zero
                signal(iter) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else if (t > L) then
                !After the end position, the signal is also zero
                signal(iter) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            else
                !Gaussian pulse centered at peak with standard deviation sigma
                signal(iter) = cmplx(amp * exp(-0.5_dp * ((t - peak) / sigma)**2.0_dp), 0.0_dp, kind=dp)
            end if
            
            ! Write data out to file
            write(unit=out_a, fmt=*, iostat=ios) t, signal(iter)%re
            if(ios.ne.0) stop 'Error writing to gaussian_pulse.dat'
        end do
    
        close(unit=out_a, iostat=ios)
        if(ios.ne.0) stop 'Error writing to file gaussian_pulse.dat'
    
    end subroutine gaussian_pulse
    
end module DFT 

program signal_processing 
    use double_precision 
    use constants, only: two_pi, i 
    use IO, only : ios, out_a, out_b, out_c, out_d 
    use DFT 
    implicit none 

    !Define variables
    integer :: N                        !Number of samples
    integer :: iter                     !Iteration variable
    real(kind=dp) :: period             !Period of input signal
    real(kind=dp) :: t                  !Time
    real(kind=dp) :: dt                 !Timestep size
    real(kind=dp) :: amp                !Amplitude of input signal 
    real(kind=dp) :: freq               !Frequency of input signal
    real(kind=dp) :: nyquist_freq       !Nyquist critical frequency
    real(kind=dp) :: sample_rate        !Sampling rate (No of samples recorded per second)
    real(kind=dp) :: freq_resolution    !Frequency resolution: Width of each frequency bin in DFT
    real(kind=dp) :: duration           !Total time of sample window
    real(kind=dp) :: real_freq          !Frequency shift for x-axis of DFT plot
    real(kind=dp) :: L                  !Length of finite pulse
    real(kind=dp) :: sigma              !Width input for gaussian pulse subroutine
    real(kind=dp) :: start_pos          !Start position of rectangle pulse

    !Define matrices
    complex(kind=dp), dimension(:,:), allocatable :: W_matrix   !The DFT matrix
    complex(kind=dp), dimension(:), allocatable :: t_data       !Time data matrix
    complex(kind=dp), dimension(:), allocatable :: f_data       !Freq data matrix

    !Case selectors
    character(len=3) :: input       !Input signal selector
    character(len=3) :: window      !Window function selector
    
    !---------------------------------------------------------------------------------------------------!
    !                                       PROGRAM CONTROLS                                            !
    !---------------------------------------------------------------------------------------------------!
    !                                                                                                   !
    !       INPUT SIGNAL SELECTION                                                                      !
    !                                                                                                   !
    !       input = 'sin' -> sine wave                                                                  !
    !       input = 'rec' -> finite rectangle pulse                                                     !
    !       input = 'tri' -> finite symmetric triangle pulse                                            !
    !       input = 'gau' -> finite gaussian pulse                                                      !
    !                                                                                                   !
    !       WINDOW FUNCTION SELECTION                                                                   !
    !                                                                                                   !
    !       window = 'non' -> no window function                                                        !
    !       window = 'tri' -> triangle window function                                                  !
    !       window = 'cos' -> cosine bell window function                                               !
    !       window = 'gau' -> gaussian bell window function                                             !
    !                                                                                                   !
    !---------------------------------------------------------------------------------------------------!

    !Selectors
    input = 'sin'
    window = 'non'
    sigma = 1.0_dp !For gaussian pulse input (max = 1.0_dp)

    !Input signal settings
    L = 10.0_dp      
    amp = 1.0_dp 
    !Only apply to sine wave
    freq = 1.0_dp
    period = 1.0_dp / freq

    !Sampling settings
    duration = 10.0_dp 
    sample_rate = 400.0_dp 
    !Scale N with sample rate and duration, max function ensures that N >= 100 
    N = max(100, nint(sample_rate * duration))
    freq_resolution = sample_rate / real(N, kind=dp)
    dt = 1.0_dp/sample_rate 
    nyquist_freq = 1.0_dp / (2.0_dp * dt)

    !---------------------------------------------------------------------------------------------------

    !Print values to console

    print *, 'INPUT SIGNAL:'
    print *, 'input=', input 
    print *, 'window=', window 
    print *, 'amp=', amp
    !Print freq and period if sine wave
    if (input == 'sin') then 
        print *, 'freq=', freq 
        print *, 'period=', period
    end if 
    print *, 'L=', L 
    !Print sigma if gaussian pulse is selected
    if (input == 'gau') then 
        print *, 'sigma=', sigma 
    end if 

    print *, '' !space

    print *, 'SAMPLE WINDOW:'
    print *, 'N=', N
    print *, 'sample_rate=', sample_rate
    print *, 'nyquist freq=',  nyquist_freq
    print *, 'freq_resolution=', freq_resolution
    print *, 'duration=', duration 
    
    !---------------------------------------------------------------------------------------------------!
    !                                       INPUT SIGNAL GENERATION                                     !
    !---------------------------------------------------------------------------------------------------!
    
    !Define t_data from 0 to N-1 for a size N array
    allocate(t_data(0:N-1), stat=ios)
    if(ios.ne.0) stop ' Error allocating t_data'

    select case (input)
        case('sin')

            open(unit=out_a, file='sin.dat', iostat=ios)
            if(ios.ne.0) stop 'Error opening file sin.dat'

            !Occupy sine wave data
            t = 0.0_dp
            do iter = 0, N-1, 1
                t = iter * dt 
                t_data(iter) = cmplx(amp*sin(two_pi*freq*t), 0.0_dp,kind=dp)
                write(unit=out_a, fmt=*, iostat=ios) t, t_data(iter)%re 
                if(ios.ne.0) stop 'Error writing to sin.dat'
            end do 

            close(unit=out_a, iostat=ios)
            if(ios.ne.0) stop 'Error writing to file sin.dat'
        
        case('rec')
            start_pos = duration - L 
            call rectangular_pulse(N, t_data, amp, L, sample_rate, start_pos)
        case('tri')
            start_pos = duration - L 
            call triangle_pulse(N, t_data, amp, L, sample_rate, start_pos)
        case('gau')
            start_pos = duration - L 
            call gaussian_pulse(N, t_data, amp, L, sigma, sample_rate, start_pos)
        end select 

    !---------------------------------------------------------------------------------------------------!
    !                                       DFT GENERATION                                              !
    !---------------------------------------------------------------------------------------------------!
    
    !Calculate DFT
    call DFT_coeff(N, W_matrix)                 !Generate DFT coefficient matrix (W_matrix)
    call do_DFT(N, t_data, W_matrix, f_data)    !Calculate DFT        
    
    !Open files
    open(unit=out_a, file='dft_re.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening file dft_re.dat'
    open(unit=out_b, file='dft_im.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening file dft_im.dat'
    open(unit=out_c, file='dft_mag.dat', iostat=ios)
    if(ios.ne.0) stop 'Error opening file dft_mag.dat'
    
    !Write to file 
    do iter = 0, N-1, 1 
        !Calculate frequencies respective to positive and negative freq spectrum
        if (iter <= N/2) then
            real_freq = real(iter, kind=dp) * freq_resolution
        else
            real_freq = real(iter - N, kind=dp) * freq_resolution
        end if

        !Real part
        write(unit=out_a, fmt=*, iostat=ios) real_freq, f_data(iter)%re
        if(ios.ne.0) stop 'Error writing to dft_re.dat'
        !Imaginary part
        write(unit=out_b, fmt=*, iostat=ios) real_freq, f_data(iter)%im
        if(ios.ne.0) stop 'Error writing to dft_im.dat'
        !Magnitude Normalised from -sample_rate/2 to +sample_rate/2
        write(unit=out_c, fmt=*, iostat=ios) real_freq , abs(f_data(iter))
        if(ios.ne.0) stop 'Error writing to dft_mag.dat'

    end do 

    !Close files
    close(unit=out_a, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file dft_re.dat'
    close(unit=out_b, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file dft_im.dat'
    close(unit=out_c, iostat=ios)
    if(ios.ne.0) stop 'Error writing to file dft_mag.dat'

    !---------------------------------------------------------------------------------------------------!
    !                                       WINDOWING FUNCTION                                          !
    !---------------------------------------------------------------------------------------------------!
    
    ! window = non -> No windowing applied 
    ! window = tri -> triangle window function applied 
    ! window = cos -> cosine bell window function applied
    ! window = gau -> gaussian bell window function applied 

    select case (window)
    case('tri')
        call triangle_window(N, t_data)
    case('cos')
        call cosine_bell_window(N, t_data)
    case ('gau')
        call gaussian_bell_window(N, t_data)
    end select 

    !If window function is applied then output data as a dat file
    if (window == 'tri' .or. window == 'cos' .or. window == 'gau') then 
    
        t = 0.0_dp 
        open(unit=out_a, file='windowed_signal.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file windowed_signal.dat'

        do iter = 0, N-1, 1 
            t = iter * dt 
            !Write out input signal with window function applied
            write(unit=out_a, fmt=*, iostat=ios) t, t_data(iter)%re 
            if(ios.ne.0) stop 'Error writing to file windowed_signal.dat'
        end do 
        
        close(unit=out_a, iostat=ios)
        if(ios.ne.0) stop 'Error writing to file windowed_signal.dat'

        !Calculate DFT of windowed signal
        call DFT_coeff(N, W_matrix)
        call do_DFT(N, t_data, W_matrix, f_data)
        
        !Open files
        open(unit=out_a, file='dft_re_windowed.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file dft_re_windowed.dat'
        open(unit=out_b, file='dft_im_windowed.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file dft_im_windowed.dat'
        open(unit=out_c, file='dft_mag_windowed.dat', iostat=ios)
        if(ios.ne.0) stop 'Error opening file dft_mag_windowed.dat'
        
        !Write to file 
        do iter = 0, N-1, 1 
            !Calculate frequencies respective to positive and negative freq spectrum
            if (iter <= N/2) then
                real_freq = real(iter, kind=dp) * freq_resolution
            else
                real_freq = real(iter - N, kind=dp) * freq_resolution
            end if

            !Real part
            write(unit=out_a, fmt=*, iostat=ios) real_freq, f_data(iter)%re
            if(ios.ne.0) stop 'Error writing to dft_re_windowed.dat'
            !Imaginary part
            write(unit=out_b, fmt=*, iostat=ios) real_freq, f_data(iter)%im
            if(ios.ne.0) stop 'Error writing to dft_im_windowed.dat'
            !Magnitude Normalised from -sample_rate/2 to +sample_rate/2
            write(unit=out_c, fmt=*, iostat=ios) real_freq , abs(f_data(iter))
            if(ios.ne.0) stop 'Error writing to dft_mag_windowed.dat'

        end do 
    end if 

end program signal_processing