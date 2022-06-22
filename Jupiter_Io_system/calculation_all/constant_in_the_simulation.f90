module constant_in_the_simulation
    use constant_parameter

    implicit none

    !------------
    ! grid number
    !------------

    integer, parameter :: real_grid_number = 233    !for real space
    integer, parameter :: adiabatic_invariant_grid_number = 150    !for adiabatic invariant


    !-------------------
    ! planet's constants
    !-------------------

    double precision, parameter :: planet_dipole_moment = 1.16d27  ![m2 A]
    double precision, parameter :: planet_mass = 1.8982d27  ![kg]
    double precision, parameter :: planet_radius = 7.1492d7     ![m]
    double precision, parameter :: planet_rotation = 2d0 * pi / 9.9258d0 / 3600d0   ![rad s-1]
    double precision, parameter :: planet_l_shell = 5.84760d0
    double precision, parameter :: planet_mlat_1 = - acos(sqrt((1d0 + 25d5 / planet_radius) / planet_l_shell))  ![rad]
    double precision, parameter :: planet_mlat_2 = - planet_mlat_1  ![rad]

    
    !--------------------------------------------------------------------------------------------
    ! satellite's constants (if not set the satellite at the equator, set these parameters to 0.)
    !--------------------------------------------------------------------------------------------

    double precision, parameter :: satellite_mass = 8.931938d22     ![kg]
    double precision, parameter :: satellite_l_shell = 5.89856d0

end module constant_in_the_simulation