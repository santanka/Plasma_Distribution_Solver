subroutine make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none

    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(out) :: adiabatic_invariant

    integer :: count_s, count_mu
    double precision :: max_mu

    do count_s = 1, boundary_series_number
        
        do count_mu = 1, adiabatic_invariant_grid_number

            max_mu = particle_mass(count_s) * speed_of_light**2d0 / 2d0 / magnetic_flux_density(injection_grid_number(count_s))

            if ( count_mu == 1 ) then
                adiabatic_invariant(count_s, count_mu) = 0d0
            else if ( count_mu /= 1 ) then
                adiabatic_invariant(count_s, count_mu) &
                    & = 1d-30 * (max_mu / 1d-30)**(dble(count_mu - 2) / dble(adiabatic_invariant_grid_number - 2))
            end if

        end do  !count_mu

    end do  !count_s
    
end subroutine make_adiabatic_invariant
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_electrostatic_potential_diff(electrostatic_potential,  electrostatic_potential_diff)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: electrostatic_potential
    double precision, dimension(3, real_grid_number), intent(out) ::  electrostatic_potential_diff

    electrostatic_potential_diff(1, :) = electrostatic_potential
    electrostatic_potential_diff(2, :) = electrostatic_potential + 1d-7
    electrostatic_potential_diff(3, :) = electrostatic_potential - 1d-7
        
end subroutine make_electrostatic_potential_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_boundary_number_density_diff(boundary_number_density,  boundary_number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(3, boundary_series_number), intent(out) ::  boundary_number_density_diff

    boundary_number_density_diff(1, :) = boundary_number_density

    boundary_number_density_diff(2, boundary_ionosphere_1_variable_species) &
        & = boundary_number_density(boundary_ionosphere_1_variable_species) * (1d0 + 1d-7)
    boundary_number_density_diff(3, boundary_ionosphere_1_variable_species) &
        & = boundary_number_density(boundary_ionosphere_1_variable_species) * (1d0 - 1d-7)
    
    boundary_number_density_diff(2, boundary_ionosphere_2_variable_species) &
        & = boundary_number_density(boundary_ionosphere_2_variable_species) * (1d0 + 1d-7)
    boundary_number_density_diff(3, boundary_ionosphere_2_variable_species) &
        & = boundary_number_density(boundary_ionosphere_2_variable_species) * (1d0 - 1d-7)

    boundary_number_density_diff(2, boundary_magnetosphere_variable_species) &
        & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 + 1d-7)
    boundary_number_density_diff(3, boundary_magnetosphere_variable_species) &
        & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 - 1d-7)

end subroutine make_boundary_number_density_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_energy_diff(mlat, length2planet, length2satellite, charge_number, particle_mass, &
    & electrostatic_potential_diff, potential_energy_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: mlat, length2planet, length2satellite
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) ::  potential_energy_diff

    integer :: count_h, count_i

    do count_h = 1, 3
        
        do count_i = 1, real_grid_number

            !gravity of planet
            potential_energy_diff(count_h, :, count_i) = - constant_of_gravitation * planet_mass * particle_mass &
                & / length2planet(count_i)
            
            !centrifugal force of planet
            potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                & - particle_mass * (planet_rotation * length2planet(count_i) * cos(mlat(count_i)))**2d0 / 2d0
            
            !gravity of satellite
            if ( satellite_mass /= 0d0 ) then
                potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                    & - constant_of_gravitation * satellite_mass * particle_mass / length2satellite(count_i)
            end if

            !Coulomb force
            potential_energy_diff(count_h, :, count_i) = potential_energy_diff(count_h, :, count_i) &
                & + charge_number * electrostatic_potential_diff(count_h, count_i)
        
        end do  !count_i
    
    end do  !count_h
    
end subroutine make_potential_energy_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!