! Program: nuclear_reactor
! By: Noah Egger
!-----------------------------------------------------------------------------
! This program calculates the neutron flux from a box reactor of variable 
! size via a detector placed in the z=0 plane at a distance x = x_0 and 
! y=y_0. The program uses Newton-Coates and Monte Carlo quadrature
! to achieve this calculation.
!
! The additional program (advanced project) similarly caclculates the neutron 
! flux from the box reactor except with a hollow spherical inside of radius
! specified by the user. It is achieved by subtracting the flux of the sphere
! from the flux of the box reactor. The advanced project uses both Newton-Coates
! and MC quadrature to achieve this calculation.
!
!

! You're free to give different values when running the code, the ones 
! bellow are just a suggestion that worked fine for me.

! Basic part of the project
! depth = 40.0_dp
! width = 100.0_dp
! height = 60.0_dp
! y_zero = 30._dp
! x_min = 5_dp
! x_max = 200._dp
! x_step = 5._dp
! n_grid = 25
! m_samples = 10000

! Advanced part of the project
! x_zero = 80._dp
! r_min = 1.0_dp
! r_max = 19._dp
! r_step = 1.0_dp
!-----------------------------------------------------------------------------
program nuclear_reactor
use types 

use read_write, only : read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

implicit none

real(dp) :: depth, width, height, y_zero, x_min, x_max, x_step
integer :: n_grid, m_samples
real(dp) :: x_zero, r_min, r_max, r_step


! Basic part of the project
call read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
call write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)

! Advanced part of the project
call read_advanced_input(x_zero, r_min, r_max, r_step)
call write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)

end program nuclear_reactor