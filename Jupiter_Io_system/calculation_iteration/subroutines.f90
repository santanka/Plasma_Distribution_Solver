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
subroutine make_amin(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, injection_grid_number, amin)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) ::  magnetic_flux_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: total_energy

    !total energy
    do count_h = 1, 3

        do count_s = 1, boundary_series_number

            do count_i = 1, real_grid_number

                total_energy(count_h, count_s, count_i, :) = potential_energy_diff(count_h, count_s, count_i) &
                    & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

            end do  !count_i

        end do  !count_s

    end do  !count_h

    !amin
    do count_s = 1, boundary_series_number

        do count_i = 1, real_grid_number

            if ( count_i /= injection_grid_number(count_s) ) then
                do count_h = 1, 3
                    
                    do count_mu = 1, adiabatic_invariant_grid_number
                        
                        if ( count_i < injection_grid_number(count_s) ) then
                            Emax_grid = count_i + 1
                            do count4max = count_i + 1, injection_grid_number(count_s)

                                if ( total_energy(1, count_s, count4max, count_mu) &
                                    & > total_energy(1, count_s, Emax_grid, count_mu) ) then
                                    Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        
                        else if ( count_i > injection_grid_number(count_s) ) then
                            Emax_grid = injection_grid_number(count_s)
                            do count4max = injection_grid_number(count_s), count_i - 1

                                if ( total_energy(1, count_s, count4max, count_mu) &
                                    & > total_energy(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        end if

                        if ( total_energy(count_h, count_s, count_i, count_mu) > total_energy(1, count_s, Emax_grid, count_mu) )then
                            amin(count_h, count_s, count_i, count_mu) = sqrt(total_energy(count_h, count_s, count_i, count_mu) &
                                & - total_energy(1, count_s, injection_grid_number(count_s), count_mu))
                        else
                            amin(count_h, count_s, count_i, count_mu) = sqrt(total_energy(1, count_s, Emax_grid, count_mu) &
                                & - total_energy(1, count_s, injection_grid_number(count_s), count_mu))
                        end if

                    end do  !count_mu

                end do  !count_h

            else if ( count_i == injection_grid_number(count_s) ) then
                amin(:, count_s, count_i, :) = 0d0

            end if
        end do  !count_i
        
    end do  !count_s

    
end subroutine make_amin
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, injection_grid_number, particle_mass, &
    & amin, amax)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) ::  magnetic_flux_density
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: total_energy
    double precision, dimension(3) :: energy_difference

    !total energy
    do count_h = 1, 3

        do count_s = 1, boundary_series_number

            do count_i = 1, real_grid_number

                total_energy(count_h, count_s, count_i, :) = potential_energy_diff(count_h, count_s, count_i) &
                    & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

            end do  !count_i

        end do  !count_s

    end do  !count_h

    !amax
    do count_s = 1, boundary_series_number

        do count_i = 1, real_grid_number

            if ( (count_i == 1 .and. injection_grid_number(count_s) /= 1) &
                & .or. (count_i == real_grid_number .and. injection_grid_number(count_s) /= real_grid_number) ) then
                amax(:, count_s, count_i, :) = 0d0

            else if ( count_i == injection_grid_number(count_s) .and. injection_grid_number(count_s) /= 1 &
                & .and. injection_grid_number(count_s) /= real_grid_number ) then
                amax(:, count_s, count_i, :) = 1d100
            
            else
                do count_mu = 1, adiabatic_invariant_grid_number

                    if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                        Emax_grid = 1
                        do count4max = 1, count_i

                            if ( total_energy(1, count_s, count4max, count_mu) &
                                & > total_energy(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    
                    else if (count_i > injection_grid_number(count_s) .or. injection_grid_number(count_s) == 1) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, real_grid_number
                            
                            if ( total_energy(1, count_s, count4max, count_mu) &
                            & > total_energy(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    end if

                    energy_difference = total_energy(1, count_s, Emax_grid, count_mu) - total_energy(:, count_s, count_i, count_mu)
                    do count_h = 1, 3
                        
                        if ( energy_difference(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = 0d0

                        else if ( amin(count_h, count_s, count_i, count_mu)**2d0 < energy_difference(count_h) &
                            & .and. energy_difference(count_h) < 5d-1 * particle_mass(count_s) *speed_of_light**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference(count_h))

                        else if ( energy_difference(count_h) >= 5d-1 * particle_mass(count_s) *speed_of_light**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(particle_mass(count_s) / 2d0) * speed_of_light
                        end if

                    end do  !count_h
                    
                end do  !count_mu

            end if
            
        end do  !count_i

    end do  !count_s

    
end subroutine make_amax
!
!-----------------------------------------------------------------------------------------------------------------------------------
!