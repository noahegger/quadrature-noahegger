!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This module receives user specific inputs in order to establish
!! the reactor dimensions and calculate the flux for specified
!! location, as well as specifying the size of the spherical
!! hollow region.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use neutron_flux, only : box_flux_booles, large_x0_flux, box_flux_monte_carlo, total_flux_booles, hollow_box_flux_mc

implicit none

private
public :: read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Receives user inputs to specify the program. For example,
!! depth, width, and height of the reactor. Also allows user to specify
!! the detector location as well as the min and max x value of position.
!! The user also supplies the number of grid points and sample points.
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(out) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(out) ::  n_grid, m_samples

    ! use the print statement to give a message to the user describing
    ! what the program does and what input should the user give to the
    ! program  

    print*,'This program calculates the neutron flux from a box reactor.'
    print*,'You will be prompted to specify the sizing of the reactor.'
    print*,'The flux is measured a variable distance away from the reactor.'
    print*,'You will define the y-coordinate of the detector.'
    print*,'Also, you will define the range of x distances away from reactor as well as the size increment.'
  

    depth = read_real('depth D')
    width = read_real('width W')
    height = read_real('height H')
    y_zero = read_real('y coordinate of detector')
    x_min = read_real('min x detector location')
    x_max = read_real('max x detector location')
    x_step = read_real('size increment for x coordinates')

    ! The function read_real used above returns a double precision real,
    ! however n_grin and m_samples are integers, that means that we need
    ! another function to get integers. For that we'll define read_integer
    ! as well

    n_grid = read_integer('lattice points N')
    m_samples = read_integer('Monte Carlo samples M')

end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Verifies user specified inputs are of correct type, i.e. real.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print*, 'Provide a nonzero, positive value for the '//trim(name)//':'

    ! Using the do loop to get an input from the user.



     do
         read(*,'(a)',iostat=ierror) string
         ! Input is nonempty, proceed
         if(string.ne.'') then
            read(string,*,iostat=ierror) x
            ! If input can be formulated as number, proceed
            if (ierror == 0) then
                ! if number is positive, exit
                if (x > 0) exit
                print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive number'
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif
        else
            print *, 'that was an empty input, please a provide positive number'
        endif

     enddo
end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Verifies inputs are the correct type, i.e integers.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print*, 'Provide a nonzero positive value for the '//trim(name)//':'
  
     do
        read(*,'(a)',iostat=ierror) string
        if(string.ne.'') then
            read(string,*,iostat=ierror) x
            if (ierror == 0) then
                if (x > 0) exit       
                    print *, "'"//trim(string)//"'"// 'cannot be negative or zero, please provide a positive number' 
            else
                print *, "'"//trim(string)//"'"//' is not a number, please provide a number'
            endif          
        else
            print *, 'that was an empty input, please provide a positive number'
        
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!Subroutine: write_neutron_flux
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Writes the results of computations for basic box flux into a file named 
!! 'results_basic.dat'. Also writes results for given detector location.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: x_zero, box_booles, box_mc, box_large_x0, sigma_box
    character(len=*), parameter :: file_name = 'results_basic.dat'
    integer :: unit

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'x_0', 'booles', 'large x_0', 'monte carlo', 'MC uncertainty'
    x_zero = x_min
    do 
        if(x_zero > x_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        box_large_x0 = large_x0_flux(depth, width, height, x_zero, y_zero)
        call box_flux_monte_carlo(depth, width, height, x_zero, y_zero, m_samples, box_mc, sigma_box)
        write(unit,'(5e28.16)') x_zero, box_booles, box_large_x0, box_mc, sigma_box
        x_zero = x_zero + x_step
    enddo
    close(unit)
    print*, 'The fluxes were written in the '//file_name//' file'

end subroutine write_neutron_flux

! Below are the read and write subroutines for the advanced part the project.
! Remember to make them public at the top of the module. And to also `use` them
! in the main program

!-----------------------------------------------------------------------
!! Subroutine: read_advanced_input
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! This subroutine now takes user inputs to specify the hollow spherical
!! shell region. Similarly to the box reactor subroutines earlier,
!! the inputs must be verified to be correct type by using the read_real
!! routine.
!!----------------------------------------------------------------------
!! Output:
!!
!! x_zero       real        x coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!-----------------------------------------------------------------------
 subroutine read_advanced_input(x_zero, r_min, r_max, r_step)
     implicit none
     real(dp), intent(out) :: x_zero, r_min, r_max, r_step

print*, 'This portion of the program calculates flux of a reactor with a hollow spherical chamber inside.'
print*, 'You will specify the size of the hollow region inside the reactor.'
print*, 'You will also specify the range of radial sizes, along with the increment size.'

x_zero = read_real('x-coordinate of detector position.')
r_min = read_real('minumum radius of hollow sphere.')
r_max = read_real('maximum redius of hollow sphere.')
r_step = read_real('increment size for the radius')



 end subroutine read_advanced_input

!-----------------------------------------------------------------------
!Subroutine: write_advanced_flux
!-----------------------------------------------------------------------
!! By: Noah Egger
!!
!! Similar to before, this subroutine writes the flux calculations into
!! a file named resuls_advanced.dat. However, the difference is this
!! is for fluxes with a hollow spherical region. Fluxes are calculated
!! for between variable radius min and max with user defined increment size.
!!----------------------------------------------------------------------
!! Input:
!! 
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
 subroutine write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)
     implicit none
     real(dp), intent(in) :: depth, width, height, x_zero, y_zero, r_min, r_max, r_step
     integer, intent(in) :: n_grid, m_samples

     real(dp) :: radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
     character(len=*), parameter :: file_name = 'results_advanced.dat'
     integer :: unit

     open(newunit=unit, file=file_name)
     write(unit,'(5a28)') 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'
     radius = r_min

     do
        if (radius > r_max) exit

        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        hollow_booles = total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid)
        call hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, m_samples, hollow_mc, sigma_hollow)
        write(unit,'(5e28.16)') radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
      radius = radius + r_step
    enddo
    close(unit)
    print*, 'The fluxes were written in the '//file_name//' file'

 end subroutine write_advanced_flux

!     ! The goal is to compare Boole's and Monte Carlo integration when there's a hollow 
!     ! sphere inside the reactor to the calculation of the solid reactor using Boole's method
!     ! Base the rest of the subroutine on the one from the basic part of the project.
! end subroutine write_advanced_flux

end module read_write
