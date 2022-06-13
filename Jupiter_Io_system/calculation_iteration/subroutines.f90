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
    boundary_number_density_diff(2, :) = boundary_number_density
    boundary_number_density_diff(3, :) = boundary_number_density

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
                            & .and. energy_difference(count_h) < 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference(count_h))

                        else if ( energy_difference(count_h) >= 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
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
subroutine make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
    & potential_energy_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, particle_mass, amin, amax, &
    & number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, amax
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) :: number_density_diff

    integer :: count_h, count_s, count_i
    double precision :: integral, num
    double precision, dimension(adiabatic_invariant_grid_number) :: alpha_mu, amax_mu, amin_mu

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number
            
            if ( count_i == injection_grid_number(count_s) ) then
                alpha_mu = magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :) / boundary_temperature_perp(count_s)
                amax_mu = amax(1, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
                amin_mu = 0d0

                call calculation_integral_exp_erf(adiabatic_invariant(count_s, :), alpha_mu, amax_mu, amin_mu, integral)

                num = 1d0 - exp(- particle_mass(count_s) * speed_of_light**2d0 / 2d0 / boundary_temperature_perp(count_s)) &
                    & + magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) * integral
                
                number_density_diff(:, count_s, count_i) = boundary_number_density_diff(:, count_s) / 2d0 * num
                
            else if ( count_i /= injection_grid_number(count_s) ) then
                do count_h = 1, 3

                    alpha_mu = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
                        & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
                        & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)
                    amax_mu = amax(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
                    amin_mu = amin(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))

                    call calculation_integral_exp_erf(adiabatic_invariant(count_s, :), alpha_mu, amax_mu, amin_mu, integral)

                    num = 1d0 - exp(- ((boundary_temperature_para(count_s) - boundary_temperature_perp(count_s)) &
                        & * magnetic_flux_density(injection_grid_number(count_s)) + boundary_temperature_perp(count_s) &
                        & * magnetic_flux_density(count_i)) / boundary_temperature_para(count_s) / magnetic_flux_density(count_i) &
                        & * particle_mass(count_s) * speed_of_light**2d0 / 2d0 / boundary_temperature_perp(count_s))
                    
                    num = num + magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) * integral

                    num = boundary_number_density_diff(1, count_s) / 2d0 &
                        & * exp(- (potential_energy_diff(count_h, count_s, count_i) - &
                        & potential_energy_diff(1, count_s, injection_grid_number(count_s)))) * num

                end do  !count_h
                
            end if

        end do  !count_i

    end do  !count_s
    
end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_integral_exp_erf(mu, alpha_mu, amax_mu, amin_mu, integral_result)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: mu, alpha_mu, amax_mu, amin_mu
    double precision, intent(out) :: integral_result

    integer :: count_mu
    double precision :: integral_former, integral_latter

    integral_result = 0d0

    do count_mu = 1, adiabatic_invariant_grid_number - 1

        if ( amax_mu(count_mu) /= 0d0 ) then
            integral_former = exp(- alpha_mu(count_mu)) * (erf(amax_mu(count_mu)) - 2d0 * erf(amin_mu(count_mu)))

        else if ( amax_mu(count_mu) == 0d0 ) then
            integral_former = exp(- alpha_mu(count_mu)) * (- erf(amin_mu(count_mu)))

        end if

        if ( amax_mu(count_mu + 1) /= 0d0 ) then
            integral_former = exp(- alpha_mu(count_mu + 1)) * (erf(amax_mu(count_mu + 1)) - 2d0 * erf(amin_mu(count_mu + 1)))

        else if ( amax_mu(count_mu + 1) == 0d0 ) then
            integral_former = exp(- alpha_mu(count_mu + 1)) * (- erf(amin_mu(count_mu + 1)))

        end if

        integral_result = integral_result + (integral_former + integral_latter) / 2d0 * (mu(count_mu + 1) - mu(count_mu))

    end do  !count_mu
  
end subroutine calculation_integral_exp_erf
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_charge_density_from_number_density(number_density_diff, charge_number, charge_density_diff, &
    & charge_density_plus_diff, charge_density_minus_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: charge_number
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_diff, charge_density_plus_diff
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_minus_diff

    integer :: count_s

    charge_density_diff = 0d0
    charge_density_plus_diff = 0d0
    charge_density_minus_diff = 0d0

    do count_s = 1, boundary_series_number
        
        charge_density_diff = charge_density_diff + charge_number(count_s) * number_density_diff(:, count_s, :)

        if ( charge_number(count_s) > 0d0 ) then
            charge_density_plus_diff = charge_density_plus_diff + charge_number(count_s) * number_density_diff(:, count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            charge_density_minus_diff = charge_density_minus_diff + charge_number(count_s) * number_density_diff(:, count_s, :)
            
        end if

    end do  !count_s

end subroutine make_charge_density_from_number_density