!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module houses the functions that execute the volume integrals. 
!! The first integral method is Boole's quadrature. The second method
!! is Monte Carlo quadrature. The latter utilizes random sampling, 
!! whereas the former utilizes 5 points of equal spacing to approximate
!! an integral in a particular direction.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! 
!! monte_carlo_quad
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
        ! This is the interface we saw in class that allows callbacks
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function receives incominning arrays and breaks them into pieces
!! in order to be sent to the booles_rule function, when then formulates
!! the and evaluates expansion approximation given for each chunk. The 
!! function then sums over the results as well as provides a check that 
!! the incoming array is of correct size.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i

    fx_size = size(fx)

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. Implement a test 
    ! to make sure that the number of evaluated points in the fx array
    ! is the correct one

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4


    if (modulo(fx_size - 1, 4) /= 0) then
         print *, 'fx array size plus 1 has to be divisible by 4'
         stop
     endif

    ! We could implement the full integration here, however to make a cleaner,
    ! easy to read (and debug or maintain) code we will define a smaller
    ! function that returns Boole's five point rule and pass slices (1:5), (5:9),
    ! (9:13), ... of fx to such function to then add all the results. 

    s = 0._dp

     do i = 1, ((fx_size-1)/4)
         s = s + booles_rule(fx((4*i - 3):(4*i + 1)), delta_x)
     enddo
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function receives the chunks sent from booles_quadrature and
!! evaluates each location by utilizing the 5 equally spaced point formula.
!! Thus, it is evaluating f(x) for each x to then be sent back to 
!! booles_quadrature and summed over all points to give the full integral.
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Let's make an additional test to make sure that the array
    ! received has 5 and only 5 points 

     if (fx_size /= 5) then
         print*, 'Array fx must have 5 points. Currently has', fx_size
         stop
    endif
    
     fx0 = fx(1)
     fx1 = fx(2)
     fx2 = fx(3)
     fx3 = fx(4)
     fx4 = fx(5)
    
     s = (delta_x*2/45._dp)*(7*fx0 + 32*fx1 +12*fx2 + 32*fx3 + 7*fx4)

end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine generates random points and utilizes probability
!! to evaluate the integral. More specifically, random points are generated
!! within a volume and those locations are stored in x_vector. The function 
!! is evaluated at that location and stored in fx. Then, the
!! integral is computed by summing the values of the function fx over the 
!! total number of samples, and then multiplying by the volume of the region
!! while dividing by the total number of sample points.
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:)
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s
    real(dp) :: avg_square, square_avg, var_f, volume

    integer :: i, vector_size, b_size
    real(dp), allocatable :: x_vector(:), fx(:), delta_ab(:)


    vector_size = size(a)
    b_size = size(b)

    ! We're defining a Monte Carlo routine that works for an arbitrary number of 
    ! dimensions in the integral (Remember, that's the advantage of Monte Carlo integration,
    ! it's very efficient for high dimensional integrals)

    ! Since a and b give the lower and upper limits they need to have the same size.
    ! Make a check to see if they do have the same size

    if (vector_size /= b_size) then
          print *, 'a and b arrays in monte_carlo_quad have to be the same size'
         stop       
     endif

    ! Here we allocate memory for the vector containing the sample points and 
    ! for a vector that contains the evaluated function

    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples))
    allocate(delta_ab(1:vector_size))

    ! First find volume of the region

    delta_ab = b - a 
    volume = 1

    do i = 1, vector_size
        volume = volume*(delta_ab(i))
    enddo

    s = 0
    square_avg = 0
    avg_square = 0

    do i=1,n_samples

        !generates an array with random numbers in the [0,1) interval

        call random_number(x_vector) 

        !rescaling to the integration volume [a,b)

        x_vector = a + x_vector*delta_ab

        ! Sampling point was chosen
        ! Now evaluate at sampling point
        fx(i) = f(x_vector,data)

        ! Compute integral

        s = s + (volume/(float(n_samples)))*fx(i)

        ! Compute values to be used error estimation

        square_avg = square_avg + (1/float(n_samples))*fx(i)**2

        avg_square = avg_square + (1/float(n_samples))*fx(i)
    enddo

    ! Variance calculation

    var_f = square_avg - avg_square**2

    ! STD calculation

    sigma_s = volume * sqrt(var_f)/sqrt(float(n_samples))

end subroutine monte_carlo_quad

end module quadrature