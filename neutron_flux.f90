!-----------------------------------------------------------------------
!Module: neutron_flux
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module contains subroutines that calculate the flux for the 
!! box reactor both with and without a hollow spherical chamber of specif-
!! ied radius = "radius". More detail is given within the individual
!! descriptions of subroutines.
!! 
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! box_flux_monte_carlo
!! hollow_box_flux_mc
!!----------------------------------------------------------------------
!! Included functions:
!!
!! box_flux_booles
!! sphere_flux_booles
!! total_flux_booles
!! sphere_flux_kernel
!! flux_kernel
!! flux_kernel_vector
!! hollow_box_flux_kernel
!! large_x0_flux
!-----------------------------------------------------------------------
module neutron_flux
use types
use quadrature, only : booles_quadrature, monte_carlo_quad

implicit none

private
public :: box_flux_booles, large_x0_flux, box_flux_monte_carlo, hollow_box_flux_mc, total_flux_booles

contains

! Let's do the Boole's integration first

!-----------------------------------------------------------------------
!! Function: box_flux_booles
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine computes the 3D flux integral over the volume of the box.
!! Each particular do statement performs a loop which fills an array and
!! is sent to the quadrature module to perform the integration in the 
!! particular direction (x,y,z). The total flux is calculated by sending
!! the final 3D array to the booles quadrature. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_grid
    
    real(dp) :: delta_x, delta_y, delta_z
    real(dp), allocatable :: f_x(:), g_xy(:), h_xyz(:)
    integer :: n_bins, i_x, i_y, i_z
    real(dp) :: x, y, z

    ! I'll give you some help here. 

    ! First we need to determine the distance between
    ! the lattice points at which the function to integrate
    ! will be evaluated 

    ! Hopefully the little diagram bellow will help you
    ! figure it out

    ! bins:              1   2   3   4   5   6   7   8
    ! x interval:      |---|---|---|---|---|---|---|---|
    ! grid points:     1   2   3   4   5   6   7   8   9 
    
    ! interval length: |-------------------------------|
    !                  0                               depth
    ! delta x length:  |---|
    !                  0   delta_x

    n_bins = n_grid -1

    ! n_bins is the number of intervals we have.

    delta_x = depth/(float(n_bins))
    delta_y = width/(float(n_bins))
    delta_z = height/(float(n_bins))

    ! We also have step sizes for x,y,z defined above. Achieved by simply
    ! dividing the distance by the number of bins.


    ! Now we need to allocate memory for the arrays that will contain
    ! the evaluated function to integrate

    allocate(  f_x(1:n_grid))
    allocate( g_xy(1:n_grid))
    allocate(h_xyz(1:n_grid))

    ! Now we need to implement the do loop in the README file

     x = -delta_x
    ! x is defined so that integration starts at x=0

    do i_x = 1, n_grid

        ! First loop in x-direction (depth). Terminated if x reaches
        ! maximum distance = "depth" user input.

        if (x < depth) then
            x = x + delta_x
        else
            exit
        endif

        y = -delta_y
        ! y is defined so that integration starts at y=0

        do i_y = 1, n_grid

        ! Second loop in y-direction (width). Terminated if y reaches
        ! maximum distance = "width" user input.

        if (y < width) then
            y = y + delta_y
        else
            exit
        endif

        z = -delta_z
        ! z is defined so that integration starts at z=0

            do i_z = 1, n_grid

                if (z < height) then
                    z = z + delta_z
                else
                    exit
                endif
                ! Now you can fill the h_xyz array with an evaluation of the function
                ! to integrate. To make things cleaner will define a function 
                ! below that returns the value we want

                h_xyz(i_z) = flux_kernel(x, y, z, x_zero, y_zero)

                ! flux_kernel contains 3D integrand. The array h_xyz is filled and 
                ! sent to booles_quadrature to be calculated across x,y, and z.
            enddo
            
            ! We fill an array by sending the previous array (the integrand) 
            ! and using booles_quadrature to evaluate integral over z.

            g_xy(i_y) = booles_quadrature(h_xyz, delta_z)

        enddo

        ! Filling a new array by sending previous evaluated array over z
        ! and using booles_quadrature to evaluate integral over y.

        f_x(i_x) = booles_quadrature(g_xy, delta_y)

    enddo

    ! Filling the final array by sending previous array evaluated over
    ! z and y and using booles_quadrature to evaluate integral over x. 
    ! this should then provide the flux for the solid box at a specified
    ! detector position.

    flux = booles_quadrature(f_x, delta_x)

end function box_flux_booles

!-----------------------------------------------------------------------
!! Function: flux_kernel
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function takes in the inputs and defines the integrand that is
!! to be computed by the function box_flux_booles.
!!
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x            real        x coordinate of the small integration volume
!! y            real        y coordinate of the small integration volume
!! y            real        z coordinate of the small integration volume
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function flux_kernel(x, y, z, x0, y0) result(k)
    implicit none
    real(dp), intent(in) :: x, y, z, x0, y0

    ! The function to be integrated goes here.
    ! The value of pi = 3.141592... is already available to you 
    ! because it's defined in the types module (take a look!)
    
    k = (1/(4._dp*pi*((x+x0)**2 + (y-y0)**2 + (z**2))))

end function flux_kernel

!-----------------------------------------------------------------------
!! Function: large_x0_flux
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function computes the flux integrand assuming a point source detector. 
!! Thus, x_0 is assumed to be large.
!! 
!!
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! d            real        Depth of the rectangular nuclear reactor
!! w            real        Width of the rectangular nuclear reactor
!! h            real        Height of the rectangular nuclear reactor
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!!----------------------------------------------------------------------
real(dp) function large_x0_flux(d, w, h, x0, y0) result(flux)
    implicit none
    real(dp), intent(in) :: d, w, h, x0, y0

    flux = (d*w*h)/(4._dp*pi*((x0+d/2._dp)**2 + (y0-w/2._dp)**2 + (h/2._dp)**2))

end function large_x0_flux

!-----------------------------------------------------------------------
!! Subroutine: box_flux_monte_carlo
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine calls upoon Monte Carlo integration to achieve the same ends
!! as the Newton-Coates integration --- that is, calculating the flux
!! of the box reactor. More specifically, the subroutine fills the arrays
!! a(:), b(:), and data(:) and sends them to monte_carlo_quad to compute
!! the flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine box_flux_monte_carlo(depth, width, height, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: flux, sigma_f
    
    real(dp) :: a(1:3), b(1:3), data(1:2)

    ! This I'll give you for free!

    ! a is the lower integration limit in the x, y, z coordinates. 
    ! since the origin was placed at the corner of the nuclear reactor the
    ! lower limit is zero in all coordinates

    a = 0._dp

    ! b is the upper integration limit in the x, y, z coordinates.

    b(1) = depth
    b(2) = width
    b(3) = height

    ! This is the 'work array' we saw in class and contains parameters
    ! (other than the sample point) needed to evaluate the function to
    ! integrate

    data(1) = x_zero
    data(2) = y_zero

    ! We're sending a function called flux_kernel_vector to
    ! the Monte Carlo subroutine to calculate the flux. 
    ! We define that function below. 

    call monte_carlo_quad(flux_kernel_vector, a, b, data, n_samples, flux, sigma_f)
end subroutine box_flux_monte_carlo

!-----------------------------------------------------------------------
!! Function: flux_kernel_vector
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function reads the both the array containing the x, y, z 
!! coordinates of the integration volume and the array containing
!! the x,y coordinates of the detector. Essentially, this function
!! organizes the two arrays so that flux_kernel can be called to fill
!! in the relevant data points for the integrand k in order to be 
!! integrated by the MC quadrature.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
! Because of the interface defined in the quadrature module the 
! Monte Carlo subroutine expects a kernel function that receives two
! arrays, the first one contains the sampling point, the second one
! contains the parameters needed to calculate the kernel. 
real(dp) function flux_kernel_vector(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0

    x = x_vector(1)
    y = x_vector(2)
    z = x_vector(3)
    x0 = data(1)
    y0 = data(2)

    ! We're going to use the function we already defined for the 
    ! Boole's integration.

    k = flux_kernel(x, y, z, x0, y0)

end function flux_kernel_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ADVANCED PART STARTS HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! As explained in the README, for Boole's method we need to calculate
! the flux of the solid box and subtract the flux from a solid sphere.
! Let's start defining the function that calculates the flux of a 
! solid sphere

!-----------------------------------------------------------------------
!! Function: sphere_flux_booles
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function uses the booles quadrature method in order to compute
!! the flux produced by a spherical region as if it is the reactor itself.
!! This is useful as we will subtract this flux from the box flux in order
!! to describe a hollow spherical shell within the box reactor. B/c the
!! volume is spherical, notice we use spherical coordinates.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! distance     real        Distance from the center of the reactor to the detector
!! radius       real        Radius of the spherical reactor
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function sphere_flux_booles(distance, radius, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: distance, radius
    integer, intent(in) :: n_grid

    real(dp) :: delta_r, delta_theta
    real(dp), allocatable :: f_r(:), g_rtheta(:)
    integer :: n_bins, i_r, i_theta
    real(dp) :: r, theta

    ! Base the rest of the function on the one for the solid box.
    ! Here we're integrating only two variables (r and theta).
    ! r is integrated from 0 to radius while theta is integrated
    ! from 0 to pi 

    ! distance is the distance from the center of the sphere 
    ! to the position of the detector (capital R in the README).
    ! It will be given as a input to this function and passed
    ! to the sphere_flux_kernel function defined below

    ! We establish the number of bins and grid points in order to
    ! find correct dr and dtheta intervals.

     n_bins = n_grid - 1

     delta_r = radius/(float(n_bins))
     delta_theta = pi/(float(n_bins))

    ! We now allocate arrays for the do loops

     allocate(g_rtheta(1:n_grid))
     allocate(f_r(1:n_grid))

    ! Next, we perform the loops
    ! Again, need r to start at -delta_r to be considered at 0

    r = -delta_r

    do i_r = 1, n_grid

        if (r < radius) then
            r = r + delta_r
            !print*, 'r=', r
        else
            exit
        endif

        ! again, need theta to start at -delta_theta

             theta = -delta_theta

            do i_theta = 1, n_grid

            if (theta < pi) then
                theta = theta + delta_theta
                !print*, 'theta=', theta
            else
                exit
            endif

                ! Now, g_rtheta calls sphere_flux_kernel to produce integrand
                ! and fill it with the theta array

                g_rtheta(i_theta) = sphere_flux_kernel(r, theta, distance)

            enddo

            ! Since integrand was established, we now integrate
            ! over theta by calling booles_quadrature

             f_r(i_r) = booles_quadrature(g_rtheta, delta_theta)

        enddo

        ! Now compute the flux by using booles_quadrature by
        ! integrating over r

        flux = booles_quadrature(f_r, delta_r)

end function sphere_flux_booles

!-----------------------------------------------------------------------
!! Function: sphere_flux_kernel
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function creates the integrand for the spherical volume in order
!! to be evaluated by sphere_flux_booles.
!!----------------------------------------------------------------------
!! Input:
!!
!! r_prime      real        r coordinate of the small integration volume
!! theta        real        theta coordinate of the small integration volume
!! big_r        integer     distance from the center of the sphere to the detector
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function sphere_flux_kernel(r_prime, theta, big_r) result(k)
    implicit none
    real(dp), intent(in) :: r_prime, theta, big_r
    
    k = (sin(theta)*r_prime**2)/(2*(r_prime**2 + big_r**2 - 2*r_prime*big_r*cos(theta)))
end function sphere_flux_kernel

!-----------------------------------------------------------------------
!! Function: total_flux_booles
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function computes the total flux by subtracting a spherical
!! flux reactor from the box flux reactor. It evaluates the distance
!! from the center of the spherical region to the detector. It calls
!! upon the box_flux and sphere_flux functions in order to subtract
!! the two and find total flux.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_grid

    real(dp) distance, box_flux, sphere_flux

    ! Now that we have a function to calculate the flux of the solid box and
    ! another one for the solid sphere we just need to use both functions 
    ! and calculate the difference.

    ! distance is the distance between the position of the detector (x_zero, y_zero)
    ! and the center of the sphere (which is also the center of the box)

    distance = sqrt((x_zero + depth/2._dp)**2 + (y_zero-width/2._dp)**2 + (height/2._dp)**2)

    ! Calculate flux from solid box at detector position

    box_flux = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)

    ! Calculate flux from solid sphere

    sphere_flux = sphere_flux_booles(distance, radius, n_grid)

    ! Compute total flux by subtracting the two

    flux = box_flux - sphere_flux


end function total_flux_booles

! As explained in the README the Monte Carlo approach is simpler.
! We just need to define a new kernel function that is zero if the 
! sampling point is inside the sphere and the original kernel if
! the sampling point is outside of the sphere.

! Again, this new kernel will take to arrays, one with the coordinates
! of the sampling point and one with all the other needed parameters
! this time it's more than just the position of the detector.

!-----------------------------------------------------------------------
!! Function: hollow_box_flux_kernel
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This function utilizes the random sampling points, the detector position,
!! and the location of the center of the sphere. The function then
!! calculates the distance from a random point within the box and tests
!! for whether it is inside or outside the sphere. If inside the sphere,
!! the location is not accounted. Otherwise, the location is counted and
!! sent to flux_kernel to write the integrand.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the sphere's radius and x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function hollow_box_flux_kernel(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0, radius, x_start, y_start, z_start
    real(dp) :: distance_to_center

    ! Random point locations inside box

    x = x_vector(1)
    y = x_vector(2)
    z = x_vector(3)

    ! Detector location

    x0 = data(1)
    y0 = data(2)
    radius = data(3)


    ! We need to determine whether or not the sampling point is inside the 
    ! sphere. For that you can calculate the distance from the sampling point
    ! (THIS IS NOT THE POSITION OF THE DETECTOR) and the center of the sphere 
    ! and compare it with the sphere's radius
    ! Remember that the origin was located at one corner of the box.

    ! Need to establish location of center of sphere

    x_start = data(4)/2._dp
    y_start = data(5)/2._dp
    z_start = data(6)/2._dp

    ! Establish distance to center

    distance_to_center = sqrt((x - x_start)**2 + (y - y_start)**2 + (z - z_start)**2)

    ! Need to distinguish inside/outside sphere 

    if (distance_to_center <= radius) then
        ! do not count for flux
        k = 0._dp

    else
        ! do count for flux
        k = flux_kernel(x, y, z, x0, y0)

    endif


end function hollow_box_flux_kernel


!-----------------------------------------------------------------------
!! Subroutine: hollow_box_flux_mc
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine organizes the data to be sent to monte_carlo_quad
!! for MC integration.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) ::  flux, sigma_f

    real(dp) :: a(1:3), b(1:3), data(1:6), deltaAB(1:3)

    ! Base the rest of the function on the one from the basic part of the project
    ! Remember that the 'work array' data now contains more information 
    
    ! Empty 'a' array for scaling

    a = 0._dp

    ! Fill b array with dimensionality

    b(1) = depth
    b(2) = width
    b(3) = height

    ! 3D array with scaled values (b-a) in x,y,z directions

    deltaAB = b - a

    ! fill data array with detector location, radius of chamber, and 
    ! box dimensions

    data(1) = x_zero
    data(2) = y_zero
    data(3) = radius
    data(4) = deltaAB(1)
    data(5) = deltaAB(2)
    data(6) = deltaAB(3)

    ! Monte_carlo_quad uses the information from the function hollow_box_flux_kernel
    ! to ensure integration happens within solid portion of box and excludes hollow sphere
    ! 
    
    call monte_carlo_quad(hollow_box_flux_kernel, a, b, data, n_samples, flux, sigma_f)
end subroutine hollow_box_flux_mc

end module neutron_flux