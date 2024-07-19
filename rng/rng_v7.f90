program rng_gen
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer :: N, i, count, interval
  real(kind=dp) :: seed, A, B, M, uniform_rand_gen
  real(kind=dp) :: uni_min, uni_max, norm_max, norm_min
  real(kind=dp),allocatable,dimension(:) :: uniform_array, norm_array, correlation_x, correlation_y, hist_x, hist_y
  !test
  real(kind=dp) :: mean, std_error, var, var_error, std_dev

  !initial values
  A = 100_dp
  B = 104001_dp
  M = 714025_dp

  !amount of random numbers to generate
  N = 10000_dp

  !set interval to split range of distribution into
  interval = 100
  
  !allocate array to size N
  allocate(uniform_array(N))
  allocate(norm_array(N))

  !generate seed from clock count
  call system_clock(count)
  seed = real(count,kind=dp)
  
  !generate N uniformly distributed random numbers
  uniform_array(1) = uniform_rand_gen(seed,A,B,M)
  do i=1,N-1
     !call function to generate normalized uniform random number in range of 0-1
     uniform_array(i+1) =  uniform_rand_gen(uniform_array(i)*M,A,B,M)
  end do

  !calculate statistics for uniform distribution
  open(unit=10,file='uniform_statistics.csv')
  write(unit=10,fmt=*) 'Statistics for each set of magic numbers'
  write(unit=10,fmt=*) 'MN Set',',','Mean',',','Standard Error',',','Variance',',','Error in Variance',',','Standard Deviation'
  call statistics(uniform_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=10,fmt=*) 'Default',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !uniform distribution histogram
  allocate(hist_x(interval))
  allocate(hist_y(interval))
  uni_max = 1.0_dp
  uni_min = 0.0_dp
  call distribute(uniform_array, N, interval,uni_max,uni_min,hist_x,hist_y)
  open(unit=12,file='histogram_uniform.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',' ,hist_y(i)
  end do
  close(unit=12)

  !generate N random numbers with a gaussian probability distribution
  !reset seed
  call system_clock(count)
  seed=real(count,kind=dp)
  !generate two normal distribution random numbers
  do i=1,N,2
  call normal_rand_gen(seed,A,B,M,norm_array(i),norm_array(i+1))
  seed = uniform_rand_gen(seed,A,B,M)*M
  end do

  !calculate statistics for gaussian distribution
  open(unit=11,file='gaussian_statistics.csv')
  write(unit=11,fmt=*) 'Statistics for each set of magic numbers'
  write(unit=11,fmt=*) 'MN Set',',','Mean',',','Standard Error',',','Variance',',','Error in Variance',',','Standard Deviation'
  call statistics(norm_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=11,fmt=*) 'Default',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !gaussian distribution histogram
  norm_max = mean + 3 * std_dev
  norm_min = mean - 3 * std_dev
  call distribute(norm_array,N,interval,norm_max,norm_min,hist_x,hist_y)
  open(unit=12,file='histogram_gaussian.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in Bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',', hist_y(i)
  end do
  close(unit=12)

  !3c correlations --------------------------------------------------------------------------
  !uniform correlations
  open(unit=13,file='uniform_correlation.csv')
  allocate(correlation_x(N/2))
  allocate(correlation_y(N/2))
  do i=1,N/2
    correlation_x(i) = uniform_rand_gen(seed,A,B,M)
    seed = correlation_x(i)*M
    correlation_y(i) = uniform_rand_gen(seed,A,B,M)
    write(unit=13,fmt=*) correlation_x(i), ',' ,correlation_y(i)
  end do
  close(unit=13)

  correlation_x = 0
  correlation_y = 0

  !gaussian correlations
  open(unit=13,file='normal_correlation.csv')
  do i=1,N/2
    call normal_rand_gen(seed,A,B,M,correlation_x(i),correlation_y(i))
    write(unit=13,fmt=*) correlation_x(i), ',' ,correlation_y(i)
  end do
  close(unit=13)

  !Part 3d - periodicity
  !PERIODICITY SUBROUTINE HERE


  !PART 4 NEW MAGIC NUMBERS
  !MAGIC NUMBERS 1 ---------------------------------------------------------------------------
  A = 1001_dp
  B = 100000_dp
  M = 714025_dp

  !generate N uniformly distributed random numbers
  uniform_array(1) = uniform_array(N)*M
  do i=1,N-1
     !call function to generate normalized uniform random number in range of 0-1
     uniform_array(i+1) =  uniform_rand_gen(uniform_array(i)*M,A,B,M)
     !set seed to non-normalized random number
  end do

  !calculate statistics for uniform distribution
  open(unit=10,file='uniform_statistics.csv')
  call statistics(uniform_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=10,fmt=*) 'Magic 1',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !distribute uniform rng
  uni_max = 1.0_dp
  uni_min = 0.0_dp
  call distribute(uniform_array, N, interval,uni_max,uni_min,hist_x,hist_y)
  open(unit=12,file='histogram_uniform_magic_1.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',' ,hist_y(i)
  end do
  close(unit=12)

  !Gaussian rng
  seed = uniform_array(N)
  do i=1,N,2
  call normal_rand_gen(seed,A,B,M,norm_array(i),norm_array(i+1))
  seed = uniform_rand_gen(seed,A,B,M)*M
  end do

  !calculate statistics for gaussian rng
  open(unit=11,file='gaussian_statistics.csv')
  call statistics(norm_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=11,fmt=*) 'Magic 1',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !Distribute gaussian rng
  norm_max = mean + 3 * std_dev
  norm_min = mean - 3 * std_dev
  call distribute(norm_array,N,interval,norm_max,norm_min,hist_x,hist_y)
  open(unit=12,file='histogram_gaussian_magic_1.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in Bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',', hist_y(i)
  end do
  close(unit=12)


  !MAGIC NUMBERS 2 ---------------------------------------------------------------------------
  A = 137_dp
  B = 150887_dp
  M = 714025_dp

  !generate uniform rng
  uniform_array(1) = uniform_array(N)
  do i=1,N-1
     uniform_array(i+1) =  uniform_rand_gen(uniform_array(i)*M,A,B,M)
  end do

  !statistics uniform rng
  open(unit=10,file='uniform_statistics.csv')
  call statistics(uniform_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=10,fmt=*) 'Magic 2',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !distribute uniform rng
  uni_max = 1.0_dp
  uni_min = 0.0_dp
  call distribute(uniform_array, N, interval,uni_max,uni_min,hist_x,hist_y)
  open(unit=12,file='histogram_uniform_magic_2.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',' ,hist_y(i)
  end do
  close(unit=12)

  !Gaussian rng
  seed = uniform_array(N)*M
  do i=1,N,2
  call normal_rand_gen(seed,A,B,M,norm_array(i),norm_array(i+1))
  seed = uniform_rand_gen(seed,A,B,M)*M
  end do

  !calculate statistics for gaussian rng
  open(unit=11,file='gaussian_statistics.csv')
  call statistics(norm_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=11,fmt=*) 'Magic 2',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !Distribute gaussian rng
  norm_max = mean + 3 * std_dev
  norm_min = mean - 3 * std_dev
  call distribute(norm_array,N,interval,norm_max,norm_min,hist_x,hist_y)
  open(unit=12,file='histogram_gaussian_magic_2.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in Bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',', hist_y(i)
  end do
  close(unit=12)

  !MAGIC NUMBERS 3 ---------------------------------------------------------------------------
  A = 1103515245_dp
  B = 12345_dp
  M = 2_dp**(32_dp)-1_dp

!generate uniform rng
  uniform_array(1) = uniform_array(N)
  do i=1,N-1
     uniform_array(i+1) =  uniform_rand_gen(uniform_array(i)*M,A,B,M)
  end do

  !statistics uniform rng
  open(unit=10,file='uniform_statistics.csv')
  call statistics(uniform_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=10,fmt=*) 'Magic 3',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !distribute uniform rng
  uni_max = 1.0_dp
  uni_min = 0.0_dp
  call distribute(uniform_array, N, interval,uni_max,uni_min,hist_x,hist_y)
  open(unit=12,file='histogram_uniform_magic_3.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',' ,hist_y(i)
  end do
  close(unit=12)

  !Gaussian rng
  seed = uniform_array(N)*M
  do i=1,N,2
  call normal_rand_gen(seed,A,B,M,norm_array(i),norm_array(i+1))
  seed = uniform_rand_gen(seed,A,B,M)*M
  end do

  !calculate statistics for gaussian rng
  open(unit=11,file='gaussian_statistics.csv')
  call statistics(norm_array,N,mean,std_error,var,var_error,std_dev)
  write(unit=11,fmt=*) 'Magic 3',',',mean,',',std_error,',',var,',',var_error,',',std_dev

  !Distribute gaussian rng
  norm_max = mean + 3 * std_dev
  norm_min = mean - 3 * std_dev
  call distribute(norm_array,N,interval,norm_max,norm_min,hist_x,hist_y)
  open(unit=12,file='histogram_gaussian_magic_3.csv')
  write(unit=12,fmt=*) 'Interval', ',', 'Number in Bin'
  do i = 1,interval
  write(unit=12,fmt=*) hist_x(i), ',', hist_y(i)
  end do
  close(unit=12)

end program rng_gen

!----------------------------------------------------------
!function that outputs a single uniform random number in range of [0,1]
function uniform_rand_gen(seed,A,B,M)
  implicit none
  integer,parameter :: dp = selected_real_kind(15,300)

  real(kind=dp) :: uniform_rand_gen, rand_num, seed, A, B, M

  !calculate random number
  rand_num = mod(A*seed+B,M)
  uniform_rand_gen = rand_num/M
end function uniform_rand_gen
!--------------------------------------------------------------
!subroutine that generates a random number with a normal distribution
subroutine normal_rand_gen(seed,A,B,M,out1,out2)
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  !real(kind=dp) :: normal_rand_gen
  
  !define variables
  real(kind=dp),intent(in) :: A, B, M
  real(kind=dp),intent(inout) :: seed
  real(kind=dp) :: y1,y2, y_mean, std_dev, uniform_rand_gen
  real(kind=dp),intent(out) :: out1, out2
  real(kind=dp),dimension(2) :: x
  real(kind=dp),parameter :: pi = 3.14159265358979323846264338327950288419716939937510
  !real(kind=dp),dimension(2),intent(out) :: output
  integer :: i
  
  y_mean = 0
  std_dev = sqrt(1.0)

  !generate 2 uniformly distributed random numbers x1 and x2
  do i=1,2
     !call function to generate normalized uniform random number in range of 0-1
     x(i) =  uniform_rand_gen(seed,A,B,M)
     !set seed to non-normalized random number
     seed = x(i)*M
  end do
 
  !calculate y1 and y2
  y1 = std_dev*sqrt(-2*log(x(1)))*cos(2*pi*x(2)) + y_mean
  y2 = sqrt(-2*log(x(1)))*sin(2*pi*x(2)) + y_mean
  out1 = y1
  out2 = y2
end subroutine normal_rand_gen
!-------------------------------------------------------------
!subroutine that calculates mean of a number set
subroutine mean_calc(array,N,mean)
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, intent(in) :: N
  real(kind=dp),dimension(N),intent(in) :: array
  real(kind=dp),intent(out) :: mean
  real(kind=dp) :: sum
  integer :: i

  sum = 0_dp
  
  do i=1,N
     sum = sum + array(i)
  end do
  mean = sum/N

end subroutine mean_calc
!-----------------------------------------------------------
!subroutine that calculates mean of squares of a number set
subroutine square_mean_calc(array,N,square_mean)
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,300)
  integer, intent(in) :: N
  real(kind=dp),dimension(N),intent(in) :: array
  real(kind=dp),intent(out) ::square_mean
  real(kind=dp) :: sum
  integer :: i

  sum = 0_dp
  
  !calculate sum of squares
  do i=1,N
     sum = sum + array(i)**2
  end do
  square_mean = sum/N
  
end subroutine square_mean_calc
!----------------------------------------------------
!subroutine that calculates variance of a number set
subroutine variance_calc(mean,square_mean,variance)
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  real(kind=dp),intent(in) :: mean,square_mean
  real(kind=dp),intent(out) :: variance

  variance = 0_dp
  
  !calculate variance of random number set
  variance = square_mean - mean**2
end subroutine variance_calc
!--------------------------------------------------------
!subroutine that calculates standard error in mean
subroutine std_error_calc(N,variance,output)
  implicit none
  integer,parameter :: dp = selected_real_kind(15,300)
  integer,intent(in) :: N
  real(kind=dp),intent(in) :: variance
  real(kind=dp),intent(out) :: output
  real(kind=dp) :: N_real
  N_real = real(N,kind=dp)

  !std error
  output = sqrt(variance)/sqrt(N_real)

  !analytic std error
  !output = sqrt(1.0)/100
end subroutine std_error_calc
!------------------------------------------------------
!subroutine that calculates the error in the variance
subroutine variance_error_calc(variance,N,output)
  implicit none
  integer,parameter :: dp = selected_real_kind(15,300)
  integer,intent(in) :: N
  real(kind=dp),intent(in) :: variance
  real(kind=dp),intent(out) :: output
  real(kind=dp) :: N_real

  N_real = real(N,kind=dp)
  !calculate error in variance
  output = variance*sqrt(2/(N_real-1))
  
end subroutine variance_error_calc
!----------------------------------------------------------------
!subroutine that outputs x and y coordinates for a histogram distribution
subroutine distribute(array,N,interval,x_max,x_min,x_out,y_out)
  implicit none
  integer,parameter :: dp = selected_real_kind(15,300)
  integer, intent(in) :: N, interval
  integer :: i, slot
  real(kind=dp),dimension(N) :: array
  real(kind=dp),intent(in) :: x_max, x_min
  real(kind=dp),dimension(interval),intent(out) :: x_out, y_out
  real(kind=dp) :: dx, x

  !initialize output arrays
  x_out = 0
  Y_out = 0

  !set bin width
  dx = (x_max-x_min)/interval

  !if random number is in range then place it in a bin
  do i=1,N
    x = array(i)
    if (x >= x_min .and. x <= x_max) then
    !round down to nearest integer to get index of bin
      slot = int((x-x_min)/dx) +1
      !insert item in bin
      y_out(slot) = y_out(slot)+1
    end if
  end do

  !map x values to array
  x_out(1) = x_min
  do i=2,interval
    x_out(i) = x_out(i-1) + dx
  end do

end subroutine distribute
!---------------------------------------------------------------
!subroutine that takes care of all statistic subroutines in a single call
subroutine statistics(array,N,mean,std_error,var,var_error,std_dev)
  implicit none
  integer,parameter :: dp = selected_real_kind(15,300)
  integer,intent(in) :: N
  real(kind=dp),dimension(N),intent(in) :: array
  real(kind=dp),intent(out) :: mean,std_error,var,var_error,std_dev
  real(kind=dp) :: square_mean

  call mean_calc(array,N,mean)
  call square_mean_calc(array,N,square_mean)
  call variance_calc(mean,square_mean,var)
  std_dev = sqrt(var)
  call std_error_calc(N,var,std_error)
  call variance_error_calc(var,N,var_error)


end subroutine statistics
