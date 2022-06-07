module variables
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    !------------------
    ! calculation field
    !------------------
    
    double precision, dimension(real_grid_number) :: mlat, length2planet, length2satellite, coordinate_FA
    double precision, dimension(real_grid_number) :: magnetic_flux_density


    !------------------
    ! initial condition
    !------------------

    double precision, dimension(real_grid_number) :: initial_electrostatic_potential


    !-------------------
    ! boundary condition
    !-------------------

    double precision, dimension(boundary_series_number) :: boundary_number_density
    double precision, dimension(boundary_series_number) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(boundary_series_number) :: charge_number
    double precision, dimension(boundary_series_number) :: particle_mass
    integer, dimension(boundary_series_number) :: injection_grid_number


    !-----------------------
    ! variables (3: 0, +, -)
    !-----------------------

    double precision, dimension(3, real_grid_number) :: electrostatic_potential
    double precision, dimension(3, boundary_series_number) :: boundary_number_density_diff
    double precision, dimension(3, boundary_series_number, real_grid_number) :: potential_energy
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, amax
    double precision, dimension(3, boundary_series_number, real_grid_number) :: number_density
    double precision, dimension(3, real_grid_number) :: charge_density, charge_density_plus, charge_density_minus
    double precision, dimension(3, real_grid_number) :: charge_density_poisson


    !--------------
    ! for iteration
    !--------------

    double precision, dimension(real_grid_number) :: new_electrostatic_potential
    double precision, dimension(boundary_series_number) :: new_boundary_number_density
    double precision, dimension(3, real_grid_number) :: convergence_number
    double precision :: convergence_number_sum
    double precision :: convergence_number_sum_min

end module variables