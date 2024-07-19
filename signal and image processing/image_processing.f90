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

module pgm_manager 
    implicit none 

    contains   
    subroutine read_PGM(filename, grid)
        !reads in a pgm file and converts the image into an array of size Nx, Ny
        integer :: max_greys, i, j
        integer, dimension(:,:), allocatable :: grid 
        integer :: Nx, Ny
        character(len = *), intent(in) :: filename
        character(len=2) :: magic
        integer :: in_unit
        integer :: ios 

        in_unit = 10 
    
        11 format(a2)
        12 format(i3,1x,i3)
        13 format (i5)
    
        if (allocated(grid)) then 
            deallocate(grid, stat=ios)
            if(ios.ne.0) stop 'Error deallocating input grid'
        end if 
    
        !Open existing file of 'filename'
        open(file=filename, unit=in_unit, status='old', iostat=ios)
        if(ios.ne.0) stop 'Error opening .pgm file' 
    
        read(in_unit, 11, iostat=ios) magic
        if (ios.ne.0) stop 'Error reading magic number'
        if (magic /= 'P2') stop 'Error - file should be P2'
    
        read(in_unit, 12, iostat=ios) Nx, Ny 
        if(ios.ne.0) stop 'Error reading image dimensions'
        if (Nx<1 .or. Ny<1) stop 'Error - invalid grid size in .pgm file'

        allocate(grid(1:Nx,1:Ny), stat=ios)
        if(ios.ne.0) stop 'Error allocating grid'
    
        read(in_unit, 13, iostat=ios) max_greys
        if (max_greys<=0 .or. max_greys>256) stop 'Error - invalid max greys'
        if (ios /= 0) stop 'Error reading maximum grayscale value'
        
        do j=1,Ny 
            do i=1,Nx-17,17
                read(in_unit,*, iostat=ios) grid(i:i+16,j)
                if(ios.ne.0) stop 'Error reading i index'
            end do 
            read(in_unit,*, iostat=ios) grid(i:Nx,j)
            if(ios.ne.0) stop 'Error reading j index'
        end do 

        close(unit=in_unit, iostat=ios)
        if(ios.ne.0) stop 'Error closing .pgm file'

        ! print *, 'read_PGM - .pgm read successful'
        
    end subroutine read_PGM
    
    subroutine write_PGM(filename, Nx, Ny, grid)
        !Takes input of an integer array of size Nx,Ny containing the grey value at each pixel
        character(len = *), intent(in) :: filename 
        integer, dimension(:,:), allocatable:: grid
        integer, intent(in) :: Nx, Ny
        integer :: i, j, max_greys
        integer:: out_unit
        integer :: ios 

        out_unit = 10 
        max_greys = 255

        !Open file to write out to
        open(file=filename, unit=out_unit, status='replace', action='write', iostat=ios)
        if(ios.ne.0) stop 'Error writing to .pgm file'

        write(out_unit, 11) 'P2'        !PGM Magic Number
        write(out_unit, 12) Nx,Ny       !Width, Height
        write(out_unit,13) max_greys    !Max grey value

        do j=1,Ny
            do i=1,Nx-17,17
                write(out_unit,14) grid(i:i+16,j)
            end do 
            write(out_unit,14) grid(i:Nx,j)
        end do 
        close(unit=out_unit)

        11 format(a2)           !write out text in a 2 character wide field
        12 format(i3,1x,i3)
        13 format (i5)
        14 format (17(i3.3,1x))

        ! print *, 'write_PGM - .pgm write successful'

    end subroutine write_PGM
end module pgm_manager 

module sobel_filter 
    use double_precision, only: dp 
    use IO, only: ios 
    implicit none 

    contains 
        subroutine sobel_x(image_in, image_out)
            !Applies the sobel filter in the x-direction
            integer, dimension(:,:), allocatable, intent(in) :: image_in           !Image input
            integer, dimension(:,:), allocatable, intent(out) :: image_out         !Image output
            integer, dimension(3,3) :: GX                                          !GX operator
            integer :: sx     
            integer :: Nx, Ny                                                      !Image input size x,y                                      
            integer :: i, j                                                        !Iteration variables for image
            integer :: kx, ky                                                      !Iteration variables for sobel operator

            !Check allocation
            if (.not.allocated(image_in)) stop 'image_in is not preallocated, ensure image is mapped to matrix'

            !Get image size
            Nx = size(image_in,1)
            Ny = size(image_in,2)

            !Allocate image_out array to size of image input
            allocate(image_out(Nx-2, Ny-2), stat=ios)
            if(ios.ne.0) stop 'Error allocating image_out'

            !Define GX matrix
            GX(1,1) = -1
            GX(1,2) = 0
            GX(1,3) = 1
            GX(2,1) = -2
            GX(2,2) = 0
            GX(2,3) = 2 
            GX(3,1) = -1
            GX(3,2) = 0
            GX(3,3) = 1

            do j = 2, Ny-1, 1 
                do i = 2, Nx-1, 1 
            
                    sx = 0 

                    do ky = 1, 3, 1
                        do kx = 1, 3, 1

                            !Calculate sx
                            sx = sx + image_in(i+kx-2, j+ky-2)*GX(kx,ky)

                        end do 
                    end do 
                    !Assign value to image_out
                    image_out(i-1,j-1) = sx 
                end do 
            end do 

        end subroutine sobel_x 

        subroutine sobel_y(image_in, image_out)
            !Applies the sobel filter in the y-direction

            integer, dimension(:,:), allocatable, intent(in) :: image_in           !Image input
            integer, dimension(:,:), allocatable, intent(out) :: image_out         !Image output
            integer, dimension(3,3) :: GY                                          !GY operator
            integer :: sy      
            integer :: Nx, Ny                                                      !Input image size x,y                                     
            integer :: i, j                                                        !Iteration variables for image
            integer :: kx, ky                                                      !Iteration variables for sobel operator
            
            !Get image size
            Nx = size(image_in,1)
            Ny = size(image_in,2)

            !Check allocation
            if (.not.allocated(image_in)) stop 'image_in is not preallocated, ensure image is mapped to matrix'

            !Allocate image_out array to size of image input
            allocate(image_out(Nx-2, Ny-2), stat=ios)
            if(ios.ne.0) stop 'Error allocating image_out'

            !Define G matrix
            GY(1,1) = 1
            GY(1,2) = 2
            GY(1,3) = 1
            GY(2,1) = 0
            GY(2,2) = 0
            GY(2,3) = 0 
            GY(3,1) = -1
            GY(3,2) = -2
            GY(3,3) = -1

            do j = 2, Ny-1, 1 
                do i = 2, Nx-1, 1 
                    
                    sy = 0 

                    do ky = 1, 3, 1
                        do kx = 1, 3, 1

                            !Calculate sy
                            sy = sy + image_in(i+kx-2, j+ky-2)*GY(kx,ky)

                        end do 
                    end do 
                    !Assign value to image_out
                    image_out(i-1,j-1) = sy
                end do 
            end do 

        end subroutine sobel_y 

        subroutine sobel_magnitude(image_x, image_y, image_out)
            !Calculates the magnitude of the gradient at each pixel

            integer, dimension(:,:), allocatable :: image_x, image_y
            integer, dimension(:,:), allocatable :: image_out
            integer :: Nx, Ny
            real(kind=dp) :: G 
            integer :: i, j

            !Check allocation
            if (.not.allocated(image_x)) stop 'image_x not preallocated, run through sobel_x subroutine first'
            if (.not.allocated(image_x)) stop 'image_y not preallocated, run through sobel_y subroutine first'

            !Get image_out size
            Nx = size(image_x,1)
            Ny = size(image_x,2)

            !Allocate image_out size
            allocate(image_out(Nx,Ny), stat=ios)
            if(ios.ne.0) stop 'Error allocating image_out'

            do i = 1, Ny, 1
                do j = 1, Ny, 1

                    !Calculate magnitude of gradient at each pixel
                    G = sqrt(real(image_x(i,j)**2, kind=dp) + real(image_y(i,j)**2, kind=dp))

                    !Assign to image out
                    image_out(i,j) = nint(G) 
                end do 
            end do 

        end subroutine sobel_magnitude

end module sobel_filter 

module fft 
    use double_precision, only: dp 
    use IO, only: ios 
    implicit none 

    contains 
        subroutine create_mask(image, radius)
            !Subroutine which creates and applies mask of a given radius to the image

            !Inputs
            complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: image 
            integer, intent(in) :: radius 

            !Local variables 
            complex(kind=dp), dimension(:,:), allocatable :: mask 
            integer :: Nx, Ny 
            integer :: i, j

            !Check allocation
            if (.not.allocated(image)) stop 'create_mask - image input not preallocated'

            !Get size 
            Nx = size(image,1)
            Ny = size(image,2)


            allocate(mask(Nx,Ny), stat=ios)
            if(ios.ne.0) stop 'Error allocating mask'

            !Create circular mask of size radius
            do j = 1, Ny, 1
                do i = 1, Nx, 1
                    if (((i**2 + j**2) < radius**2) .or. ((i**2 + (ny - j)**2) < radius**2)) then
                        mask(i,j) = 0.0_dp 
                    else 
                        mask(i,j) = 1.0_dp 
                    end if 
                end do 
            end do 

            !Apply the circular mask to the image 
            do i = 1, Nx 
                do j = 1, Ny 
                    image(i,j) = image(i,j) * mask(i,j)
                end do 
            end do 

            deallocate(mask, stat=ios)
            if(ios.ne.0) stop 'Error deallocating mask'
            
        end subroutine create_mask 

        subroutine blur_func(blur, blur_width, start_x, start_y)
            !Subroutine which creates a blur function in real space

            real(kind=dp), dimension(:,:), allocatable, intent(inout) :: blur 
            integer, intent(in) :: blur_width 
            integer, intent(in) :: start_x, start_y
            integer :: i
            integer :: Nx, Ny 

            if(.not.allocated(blur)) stop 'blur_func - blur matrix not preallocated'
            
            !Get size
            Nx = size(blur,1)
            Ny = size(blur,2)

            !Initialise to zero
            blur = 0.0_dp 

            !Ensure start is within bounds
            if (start_x < 1 .or. start_x > Nx) stop 'Invalid start_x position for blur line'
            if (start_y < 1 .or. start_y > Ny) stop 'Invalid start_y position for blur line'

            !Occupy matrix with blur rectangle
            do i = 1, blur_width
                blur(start_x + i - 1, start_y) = 1.0/real(blur_width, kind=dp) 
            end do 

        end subroutine blur_func
end module fft 

program image_processing 
    use double_precision
    use constants, only: pi 
    use IO, only: ios 
    use pgm_manager, only : read_PGM, write_PGM
    use sobel_filter, only: sobel_x, sobel_y, sobel_magnitude 
    use fft, only: create_mask, blur_func
    implicit none 

    integer :: iter, iter2                                          !Iteration counters
    integer, parameter :: i64 = selected_int_kind(18)               !Decimal range 10 ^18
    integer, parameter :: fftw_estimate=64                          !Used for FFTW3
    integer(kind=i64) :: plan_forward, plan_backward                !FFTW3 plans
    real(kind=dp) :: min_val, max_val                               !Used for normalisation

    !Images
    integer :: Nx, Ny                                               !Image size x,y
    integer, dimension(:,:), allocatable :: image                   !Original image
    real(kind=dp), dimension(:,:), allocatable :: image_real        !Image in real space
    integer, dimension(:,:), allocatable :: image_sobel             !Image after sobel edge detection
    integer, dimension (:,:), allocatable :: image_x, image_y       !Sobel filter x,y components
    complex(kind=dp), dimension(:,:), allocatable :: image_fft      !FFT of the image
    complex(kind=dp), dimension(:,:), allocatable :: blur_fft       !FFT of blur function
    
    !High-pass filter
    integer :: mask_radius                                          !Radius of high-pass filter mask

    !Blur function
    real(kind=dp), dimension(:,:), allocatable :: blur              !Blur kernel
    integer :: blur_width                                           !Width of blur rectangle function
    integer :: start_x, start_y                                     !Start position of blur rectangle

    !---------------------------------------------------------------------------------------------------!
    !                                       PROGRAM CONTROLS                                            !
    !---------------------------------------------------------------------------------------------------!
    !                                                                                                   !
    !       FFT EDGE DETECTION                                                                          !
    !                                                                                                   !
    !       mask_radius - controls radius of high pass filter mask                                      !
    !                                                                                                   !
    !       DECONVOLUTION                                                                               !
    !                                                                                                   !
    !       blur_width - controls width of the blur rectangle in horizontal direction                   !
    !       start_x - controls x start position of blur rectangle                                       !
    !       start_y - controls y start position of blur rectangle                                       !
    !                                                                                                   !
    !---------------------------------------------------------------------------------------------------!
    
    !FFT EDGE DETECTION
    mask_radius = 10

    !DECONVOLUTION
    blur_width = 21
    start_x = 1
    start_y = 1

    !---------------------------------------------------------------------------------------------------!
    !                                       PART 2A - SOBEL EDGE DETECTION                              !
    !---------------------------------------------------------------------------------------------------!

    print *, 'Sobel Edge Detection - start'

    !Read PGM image file into memory 
    call read_PGM('clown.pgm', image)
    print *, 'read_PGM successful - clown.pgm'

    !Apply Sobel filter in the x-direction 
    call sobel_x(image, image_x)

    !Apply Sobel filter in the y-direction
    call sobel_y(image, image_y)

    !Combine Sobel filter x/y to make the magnitude of the gradient at each pixel 
    call sobel_magnitude(image_x, image_y, image_sobel)

    !Write out the gradient of the image to disk as a PGM file         
    call write_PGM('clown_sobel.pgm', size(image_sobel,1), size(image_sobel,2), image_sobel)
    print *, 'write_PGM succesful - clown_sobel.pgm'

    deallocate(image, stat=ios)
    if(ios.ne.0) stop 'Part 2A - Error deallocating image'
    deallocate(image_x, stat=ios) 
    if(ios.ne.0) stop 'Part 2A - Error deallocating image_x'
    deallocate(image_y, stat=ios) 
    if(ios.ne.0) stop 'Part 2A - Error deallocating image_y'
    deallocate(image_sobel, stat=ios)
    if(ios.ne.0) stop 'Part 2A - Error deallocating image_sobel'

    print *, 'Sobel Edge Detection - completed'
    print *, ''
    !-----------------------------------------------------------------------------------------------------------!
    !                                       PART 2B - FFT EDGE DETECTION                                        !
    !-----------------------------------------------------------------------------------------------------------!

    print *, 'FFT Edge Detection - start'

    !Read PGM image file into memory 
    call read_PGM('clown.pgm', image)
    print *, 'read_PGM succesful - clown.pgm'

    !Get size
    Nx = size(image,1)
    Ny = size(image,2)

    !Allocate image_real
    allocate(image_real(Nx,Ny), stat=ios)
    if(ios.ne.0) stop 'Part 2A - Error allocating image_real'

    !Allocate image_fft
    allocate(image_fft((Nx / 2 + 1), Ny), stat=ios)
    if (ios /= 0) stop "Error allocating image_fft"

    !Convert to real space
    image_real = real(image, kind=dp)

    !Deallocate image
    deallocate(image, stat=ios)
    if (ios.ne.0) stop 'Part 2B - Error deallocating image'

    !-----------------------------------------------------------------------------------------------------------
    !Use FFTW3 to do a 2D transform of the image

    !Define FFT forward plan
    call dfftw_plan_dft_r2c_2d(plan_forward, Nx, Ny, image_real, image_fft, fftw_estimate)
    if (plan_forward == 0) stop 'Part 2B - Error creating FFTW forward plan'

    !Perform FFT
    call dfftw_execute(plan_forward)

    !Destroy FFT plan
    call dfftw_destroy_plan(plan_forward)

    !Deallocate image_real
    deallocate(image_real, stat=ios)
    if (ios /= 0) stop 'Error deallocating image_real'
    
    !-----------------------------------------------------------------------------------------------------------
    !Apply the high pass filter to each element of the FFT image
    call create_mask(image_fft, mask_radius)

    !-----------------------------------------------------------------------------------------------------------
    !Back transform the image into real-space, normalise, write out to disk as PGM file 

    !Allocate image matrix for FFT output
    allocate(image_real(Nx,Ny))
    if(ios.ne.0) stop 'Part 2B - Error allocating image_real'
    allocate(image(Nx,Ny))
    if(ios.ne.0) stop 'Part 2B - Error allocating image'

    !Inverse FFT 
    call dfftw_plan_dft_c2r_2d(plan_backward, Nx, Ny, image_fft, image_real, fftw_estimate)
    if (plan_backward == 0) stop 'Part 2B - Error creating FFTW backward plan'
    call dfftw_execute(plan_backward)
    call dfftw_destroy_plan(plan_backward)

    !Deallocate image_fft
    deallocate(image_fft, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image_fft'

    image_real = abs(image_real)
    min_val = minval(image_real)
    max_val = maxval(image_real)

    do iter = 1, size(image_real,1)
        do iter2 = 1, size(image_real,2)
            !Normalise between zero and 1
            image_real(iter,iter2) = ((image_real(iter,iter2) - min_val) / (max_val - min_val))
            !Scale between 0-255
            image(iter,iter2) = nint(image_real(iter,iter2) * 255)
        end do 
    end do 


    deallocate(image_real, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image_real'

    Nx = size(image,1)
    Ny = size(image,2)

    !Write out as .pgm file 
    call write_PGM('clown_edges.pgm', Nx, Ny, image)
    print *, 'write_PGM succesful - clown_edges.pgm'

    !Deallocate
    deallocate(image, stat=ios)
    if(ios.ne.0) stop 'Part 2B - Error deallocating image'
   
    print *, 'FFT Edge Detection - completed'
    print *, ''
    !-----------------------------------------------------------------------------------------------------------!
    !                                       PART 3 - DECONVOLUTION                                              !
    !-----------------------------------------------------------------------------------------------------------!
    
    print *, 'Deconvolution - start'

    !Read PGM image file into memory 
    call read_PGM('blurred.pgm', image)
    print *, 'read_PGM succesful - blurred.pgm'

    !Get size
    Nx = size(image,1)
    Ny = size(image,2)

    !Allocate image_real
    allocate(image_real(Nx,Ny), stat=ios)
    if(ios.ne.0) stop 'Part 2A - Error allocating image_real'

    !Convert to real space 
    image_real = real(image,kind=dp)

    !Deallocate image
    deallocate(image, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image'

    allocate(image_fft((Nx / 2 + 1), Ny), stat=ios)
    if(ios.ne.0) stop 'Error allocating image_fft'

    !Use FFTW3 to do a 2D transform of the image 
    call dfftw_plan_dft_r2c_2d(plan_forward, Nx, Ny, image_real, image_fft, fftw_estimate)
    if (plan_forward == 0) stop 'Error creating FFTW forward plan'
    call dfftw_execute(plan_forward)
    call dfftw_destroy_plan(plan_forward)

    deallocate(image_real, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image_real'

    !-------------------------------------------------------
    !Define a trial blur function in real space 
    allocate(blur(Nx,Ny), stat=ios)
    if(ios.ne.0) stop 'Error allocating blur'

    call blur_func(blur,blur_width, start_x, start_y)

    !Display blur to check
    image = int(blur * 255 * blur_width)
    call write_PGM('blur.pgm', size(blur,1), size(blur,2), image)
    print *, 'write_PGM succesful - blur.pgm'

    deallocate(image,stat=ios)
    if(ios.ne.0) stop 'Error deallocating image'

    allocate(blur_fft(Nx / 2 + 1 ,Ny), stat=ios)
    if(ios.ne.0) stop 'Error allocating blur_fft'

    !Use FFTW3 to do a 2D transform of the blur function
    call dfftw_plan_dft_r2c_2d(plan_forward, Nx, Ny, blur, blur_fft, fftw_estimate)
    if (plan_forward == 0) stop 'Error creating FFTW forward plan'
    call dfftw_execute(plan_forward)
    call dfftw_destroy_plan(plan_forward)

    !Perform deconvolution of image_fft and blur_fft
    do iter = 1, size(image_fft,1)
        do iter2 = 1, size(image_fft,2)
            !Avoid dividing by small values or zero
            if (abs(blur_fft(iter,iter2)) > 1.0e-4_dp) then 
                image_fft(iter,iter2) = image_fft(iter,iter2) / blur_fft(iter,iter2)
            else 
                image_fft(iter,iter2) = 0.0_dp
            end if 
        end do 
    end do 

    !Back transform the image, normalise and write out to disk as .PGM file 

    !Allocate matrix for FFT output
    allocate(image_real(Nx,Ny), stat=ios)
    if(ios.ne.0) stop 'Error allocating image real'

    call dfftw_plan_dft_c2r_2d(plan_backward, Nx, Ny, image_fft, image_real, fftw_estimate)
    if (plan_backward == 0) stop 'Error creating FFTW backward plan'
    call dfftw_execute(plan_backward)
    call dfftw_destroy_plan(plan_backward)

    deallocate(image_fft, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image_fft'

    !Normalise the image after FFT 

    allocate(image(size(image_real,1),size(image_real,2)))

    image_real = abs(image_real)
    min_val = minval(image_real)
    max_val = maxval(image_real)

    do iter = 1, size(image_real,1)
        do iter2 = 1, size(image_real,2)
            !Normalise between zero and 1
            image_real(iter,iter2) = ((image_real(iter,iter2) - min_val) / (max_val - min_val))
            !Scale between 0-255
            image(iter,iter2) = nint(image_real(iter,iter2) * 255)
        end do 
    end do 


    deallocate(image_real, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image_real'

    !Get size
    Nx = size(image,1)
    Ny = size(image,2)

    !Write out as .pgm file 
    call write_PGM('deblur.pgm', Nx, Ny, image)
    print *, 'write_PGM succesful - deblur.pgm'

    deallocate(image, stat=ios)
    if(ios.ne.0) stop 'Error deallocating image'

    print *, 'Deconvolution - completed'

end program 
