module main_variables
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    !------------------
    ! calculation field
    !------------------
    
    double precision, dimension(real_grid_number) :: mlat, length2planet, length2satellite, coordinate_FA
    double precision, dimension(real_grid_number - 1) :: diff_coordinate_FA
    double precision, dimension(real_grid_number) :: magnetic_flux_density


    !------------------
    ! initial condition
    !------------------

    double precision, dimension(real_grid_number) :: initial_electrostatic_potential


    !-------------------
    ! boundary condition
    !-------------------

    double precision, dimension(boundary_series_number) :: boundary_number_density, initial_boundary_number_density
    double precision, dimension(boundary_series_number) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(boundary_series_number) :: charge_number
    double precision, dimension(boundary_series_number) :: particle_mass
    integer, dimension(boundary_series_number) :: injection_grid_number


    !-----------------------------
    ! variables (count_h: 0, +, -)
    !-----------------------------

    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number) :: adiabatic_invariant
    double precision, dimension(3, real_grid_number) :: electrostatic_potential_diff
    double precision, dimension(3, boundary_series_number) :: boundary_number_density_diff
    double precision, dimension(3, boundary_series_number, real_grid_number) :: potential_energy_diff
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, amax, alim
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: &
        & potential_plus_Bmu_diff
    double precision, dimension(3, boundary_series_number, real_grid_number) :: number_density_diff
    double precision, dimension(3, real_grid_number) :: charge_density_diff, charge_density_plus_diff, charge_density_minus_diff
    double precision, dimension(3, real_grid_number) :: charge_density_poisson_diff


    !--------------
    ! for iteration
    !--------------

    double precision, dimension(real_grid_number) :: electrostatic_potential
    double precision, dimension(3, real_grid_number) :: convergence_number_diff
    double precision :: convergence_number_sum
    double precision :: convergence_number_sum_min


    !--------
    ! counter
    !--------

    integer :: count_iteration, count_iteration_timer, count_h, count_s, count_i, count_mu


    !----------
    ! file name
    !----------

    character(len = 188) :: result_file, format_character


    !------
    ! dummy
    !------

    character(len = 128) :: dummy

end module main_variables