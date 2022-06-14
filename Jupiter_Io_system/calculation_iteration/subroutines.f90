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
subroutine make_potential_plus_Bmu_diff(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, &
    & potential_plus_Bmu_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) ::  magnetic_flux_density
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: &
        & potential_plus_Bmu_diff

    integer :: count_h, count_s, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            do count_i = 1, real_grid_number
                
                potential_plus_Bmu_diff(count_h, count_s, count_i, :) = potential_energy_diff(count_h, count_s, count_i) &
                    & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

            end do  !count_i
            
        end do  !count_s

    end do  !count_h

end subroutine make_potential_plus_Bmu_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amin(potential_plus_Bmu_diff, injection_grid_number, amin)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max

    do count_s = 1, boundary_series_number

        do count_i = 1, real_grid_number

            if ( count_i /= injection_grid_number(count_s) ) then
                do count_h = 1, 3
                    
                    do count_mu = 1, adiabatic_invariant_grid_number
                        
                        if ( count_i < injection_grid_number(count_s) ) then
                            Emax_grid = count_i + 1
                            do count4max = count_i + 1, injection_grid_number(count_s)

                                if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                    & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                    Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        
                        else if ( count_i > injection_grid_number(count_s) ) then
                            Emax_grid = injection_grid_number(count_s)
                            do count4max = injection_grid_number(count_s), count_i - 1

                                if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                    & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                                end if
                            
                            end do  !count4max
                        end if

                        if ( potential_plus_Bmu_diff(count_h, count_s, count_i, count_mu) &
                            & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) )then
                            amin(count_h, count_s, count_i, count_mu) = 0d0

                        else
                            amin(count_h, count_s, count_i, count_mu) &
                                & = sqrt(potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                                & - potential_plus_Bmu_diff(count_h, count_s, count_i, count_mu))

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
subroutine make_alim(potential_plus_Bmu_diff, injection_grid_number, particle_mass, amin, alim)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: alim

    integer :: count_h, count_s, count_mu, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number
            
            do count_mu = 1, adiabatic_invariant_grid_number
                
                alim(count_h, count_s, :, count_mu) = sqrt(5d-1 * particle_mass(count_s) * speed_of_light**2d0 &
                    & + potential_plus_Bmu_diff(count_h, count_s, injection_grid_number(count_s), count_mu) &
                    & - potential_plus_Bmu_diff(count_h, count_s, :, count_mu))

                do count_i = 1, real_grid_number

                    if ( alim(count_h, count_s, count_i, count_mu) <= amin(count_h, count_s, count_i, count_mu) ) then
                        alim(count_h, count_s, count_i, count_mu) = 0d0
                    end if
                
                end do  !count_i

            end do  !count_mu

        end do  !count_i

    end do  !count_h
    
end subroutine make_alim
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_plus_Bmu_diff, injection_grid_number, particle_mass, amin, amax)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max
    double precision, dimension(3) :: energy_difference, energy_difference_boundary

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

                            if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                                & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    
                    else if (count_i > injection_grid_number(count_s) .or. injection_grid_number(count_s) == 1) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, real_grid_number
                            
                            if ( potential_plus_Bmu_diff(1, count_s, count4max, count_mu) &
                            & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    end if

                    energy_difference = potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu_diff(:, count_s, count_i, count_mu)
                    
                    energy_difference_boundary = potential_plus_Bmu_diff(:, count_s, injection_grid_number(count_s), count_mu) &
                    & - potential_plus_Bmu_diff(:, count_s, count_i, count_mu) + 5d-1 * particle_mass(count_s) * speed_of_light**2d0

                    do count_h = 1, 3
                        
                        if ( energy_difference(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 .or. &
                            & energy_difference_boundary(count_h) <= amin(count_h, count_s, count_i, count_mu)**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = 0d0

                        else if ( amin(count_h, count_s, count_i, count_mu)**2d0 < energy_difference(count_h) .and. &
                            & energy_difference(count_h) <= energy_difference_boundary(count_h) ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference(count_h))
                        
                        else if ( energy_difference(count_h) > energy_difference_boundary(count_h) ) then
                            amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference_boundary(count_h))

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
    & potential_plus_Bmu_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
    & number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, alim, amax
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) :: number_density_diff

    integer :: count_h, count_s, count_i
    double precision :: integral
    double precision, dimension(adiabatic_invariant_grid_number) :: coefficient4integral, xmin, xlim, xmax

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number
            
            do count_i = 1, real_grid_number
                
                coefficient4integral = (potential_plus_Bmu_diff(count_h, count_s, count_i, :) &
                    & - potential_plus_Bmu_diff(count_h, count_s, injection_grid_number(count_s), :)) &
                    & / boundary_temperature_para(count_s) + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :) &
                    & / boundary_temperature_perp(count_s)
                
                xmin = amin(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
                xlim = alim(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
                xmax = amax(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))

                call calculation_x_mu_integral(xmin, xlim, xmax, coefficient4integral, adiabatic_invariant(count_s, :), integral)

                if ( count_i == injection_grid_number(count_s) ) then
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(count_h, count_s) &
                        & * magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) / sqrt(pi) * integral
                
                else
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(1, count_s) &
                        & * magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) / sqrt(pi) * integral
                
                end if

            end do  !count_i

        end do  !count_s

    end do  !count_h
    
end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_x_mu_integral(xmin, xlim, xmax, coefficient4integral, adiabatic_invariant, integral_result)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: coefficient4integral, xmin, xlim, xmax
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, intent(out) :: integral_result

    integer :: count_mu, count_x
    double precision, dimension(adiabatic_invariant_grid_number) :: result_halfway, xmin2xlim, xmin2xmax
    double precision :: integral_former_x, integral_latter_x

    integral_result = 0d0
    result_halfway = 0d0

    do count_mu = 1, adiabatic_invariant_grid_number

        if ( xlim(count_mu) > xmin(count_mu) ) then
            if ( xmin(count_mu) == 0d0 ) then
                xmin2xlim(1) = 0d0
                if ( xmax(count_mu) <= xmin(count_mu) ) then
                    do count_x = 2, adiabatic_invariant_grid_number

                        xmin2xlim(count_x) = 1d-30 * (xlim(count_mu) / 1d-30)**(dble(count_x - 2) &
                            & / dble(adiabatic_invariant_grid_number - 2))
                    
                    end do  !count_x
                    xmin2xmax = 0d0

                else if ( xmax(count_mu) > xmin(count_mu) ) then
                    xmin2xmax(1) = 0d0
                    do count_x = 2, adiabatic_invariant_grid_number

                        xmin2xlim(count_x) = 1d-30 * (xlim(count_mu) / 1d-30)**(dble(count_x - 2) &
                            & / dble(adiabatic_invariant_grid_number - 2))

                        xmin2xmax(count_x) = 1d-30 * (xmax(count_mu) / 1d-30)**(dble(count_x - 2) &
                            & / dble(adiabatic_invariant_grid_number - 2))
                    
                    end do  !count_x
                end if

            else if ( xmin(count_mu) /= 0d0 ) then
                if ( xmax(count_mu) <= xmin(count_mu) ) then
                    do count_x = 1, adiabatic_invariant_grid_number

                        xmin2xlim(count_x) = xmin(count_mu) * (xlim(count_mu) / xmin(count_mu))**(dble(count_x - 1) &
                            & / dble(adiabatic_invariant_grid_number - 1))
                    
                    end do  !count_x
                    xmin2xmax = 0d0

                else if ( xmax(count_mu) > xmin(count_mu) ) then
                    do count_x = 1, adiabatic_invariant_grid_number

                        xmin2xlim(count_x) = xmin(count_mu) * (xlim(count_mu) / xmin(count_mu))**(dble(count_x - 1) &
                            & / dble(adiabatic_invariant_grid_number - 1))

                        xmin2xmax(count_x) = xmin(count_mu) * (xmax(count_mu) / xmin(count_mu))**(dble(count_x - 1) &
                            & / dble(adiabatic_invariant_grid_number - 1))
                    
                    end do  !count_x
                end if

            end if

            do count_x = 1, adiabatic_invariant_grid_number - 1

                integral_former_x = exp(- (xmin2xlim(count_x))**2d0 - coefficient4integral(count_mu))
                integral_latter_x = exp(- (xmin2xlim(count_x + 1))**2d0 - coefficient4integral(count_mu))

                result_halfway(count_mu) = result_halfway(count_mu) + (integral_former_x + integral_latter_x) / 2d0 &
                    & * (xmin2xlim(count_x + 1) - xmin2xlim(count_x))

                if ( xmax(count_mu) > xmin(count_mu) ) then
                    integral_former_x = exp(- (xmin2xmax(count_x))**2d0 - coefficient4integral(count_mu))
                    integral_latter_x = exp(- (xmin2xmax(count_x + 1))**2d0 - coefficient4integral(count_mu))

                    result_halfway(count_mu) = result_halfway(count_mu) + (integral_former_x + integral_latter_x) / 2d0 &
                        & * (xmin2xlim(count_x + 1) - xmin2xlim(count_x))

                end if

                !print *, xmin2xlim(count_x), xmin2xmax(count_x), integral_former_x, integral_latter_x

            end do  !count_x
        
        end if
        
    end do  !count_mu

    do count_mu = 2, adiabatic_invariant_grid_number
        
        integral_result = integral_result + (result_halfway(count_mu - 1) + result_halfway(count_mu)) / 2d0 &
                & * (adiabatic_invariant(count_mu) - adiabatic_invariant(count_mu - 1))

    end do  !count_mu

end subroutine calculation_x_mu_integral
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
!subroutine make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
!    & potential_energy_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, amax, number_density_diff)
!    use constant_parameter
!    use constant_in_the_simulation
!    use boundary_and_initial_conditions
!
!    implicit none
!    
!    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
!    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
!    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
!    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
!    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
!    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
!    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, amax
!    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) :: number_density_diff
!
!    integer :: count_h, count_s, count_i
!    double precision :: integral, beta
!    double precision, dimension(adiabatic_invariant_grid_number) :: alpha_mu, amax_mu, amin_mu
!
!    do count_s = 1, boundary_series_number
!        
!        do count_i = 1, real_grid_number
!            
!            if ( count_i == injection_grid_number(count_s) ) then
!                do count_h = 1, 3
!                    alpha_mu = magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :) / boundary_temperature_perp(count_s)
!                    beta = 0d0
!                    amax_mu = amax(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
!                    amin_mu = 0d0
!
!                    call calculation_integral_exp_erf(adiabatic_invariant(count_s, :), alpha_mu, beta, amax_mu, amin_mu, integral)
!                
!                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(count_h, count_s) &
!                    & * magnetic_flux_density(count_i) / 2d0 / boundary_temperature_perp(count_s) * integral
!                end do  !count_h
!
!            else if ( count_i /= injection_grid_number(count_s) ) then
!                do count_h = 1, 3
!
!                    alpha_mu = (magnetic_flux_density(injection_grid_number(count_s)) / boundary_temperature_perp(count_s) &
!                        & + (magnetic_flux_density(count_i) - magnetic_flux_density(injection_grid_number(count_s))) &
!                        & / boundary_temperature_para(count_s)) * adiabatic_invariant(count_s, :)
!                    
!                    beta = (potential_energy_diff(count_h, count_s, count_i) &
!                        & - potential_energy_diff(1, count_s, injection_grid_number(count_s))) / boundary_temperature_para(count_s)
!
!                    amax_mu = amax(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
!                    amin_mu = amin(count_h, count_s, count_i, :) / sqrt(boundary_temperature_para(count_s))
!
!                    call calculation_integral_exp_erf(adiabatic_invariant(count_s, :), alpha_mu, beta, amax_mu, amin_mu, integral)
!
!                    number_density_diff(count_h, count_s, count_i) = integral
!                    !number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(count_h, count_s) &
!                    !    & * magnetic_flux_density(count_i) / 2d0 / boundary_temperature_perp(count_s) * integral
!
!                end do  !count_h
!                
!            end if
!
!        end do  !count_i
!
!    end do  !count_s
!    
!end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
!subroutine calculation_integral_exp_erf(mu, alpha_mu, beta, amax_mu, amin_mu, integral_result)
!    use constant_parameter
!    use constant_in_the_simulation
!
!    implicit none
!    
!    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: mu, alpha_mu, amax_mu, amin_mu
!    double precision, intent(in) :: beta
!    double precision, intent(out) :: integral_result
!
!    integer :: count_mu
!    double precision :: integral_former, integral_latter
!
!    integral_result = 0d0
!
!    do count_mu = 1, adiabatic_invariant_grid_number - 1
!
!        if ( amax_mu(count_mu) > amin_mu(count_mu) ) then
!            integral_former = exp(- alpha_mu(count_mu) - beta) * (1d0 + erf(amax_mu(count_mu)) - 2d0 * erf(amin_mu(count_mu)))
!
!        else if ( amax_mu(count_mu) <= amin_mu(count_mu) ) then
!            integral_former = exp(- alpha_mu(count_mu) - beta) * (1d0 - erf(amin_mu(count_mu)))
!
!        end if
!
!        if ( amax_mu(count_mu + 1) > amin_mu(count_mu + 1) ) then
!            integral_latter = exp(- alpha_mu(count_mu + 1) - beta) * (1d0 + erf(amax_mu(count_mu + 1)) &
!                & - 2d0 * erf(amin_mu(count_mu + 1)))
!
!        else if ( amax_mu(count_mu + 1) <= amin_mu(count_mu + 1) ) then
!            integral_latter = exp(- alpha_mu(count_mu + 1) - beta) * (1d0 - erf(amin_mu(count_mu + 1)))
!
!        end if
!
!        integral_result = integral_result + (integral_former + integral_latter) / 2d0 * (mu(count_mu + 1) - mu(count_mu))
!
!    end do  !count_mu
!  
!end subroutine calculation_integral_exp_erf
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine cannot_reach_check(number_density_diff, injection_grid_number)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(inout) :: number_density_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number

    integer :: cannot_reach_point
    integer :: count_h, count_s, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number
            
            if ( injection_grid_number(count_s) /= real_grid_number ) then
                cannot_reach_point = 0
                do count_i = injection_grid_number(count_s), real_grid_number

                    if( number_density_diff(count_h, count_s, count_i) == 0d0 .and. cannot_reach_point == 0) then
                        cannot_reach_point = count_i
                    end if

                    if ( cannot_reach_point /= 0 .and. count_i >= cannot_reach_point ) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                    end if
                    
                end do  !count_i
            end if

            if ( injection_grid_number(count_s) /= 1 ) then
                cannot_reach_point = 0
                do count_i = injection_grid_number(count_s), 1, -1

                    if( number_density_diff(count_h, count_s, count_i) == 0d0 .and. cannot_reach_point == 0) then
                        cannot_reach_point = count_i
                    end if

                    if ( cannot_reach_point /= 0 .and. count_i <= cannot_reach_point ) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                    end if
                    
                end do  !count_i
            end if

        end do  !count_s

    end do  !count_h

end subroutine cannot_reach_check
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
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_charge_density_from_Poisson_eq(electrostatic_potential_diff, diff_coordinate_FA, charge_density_poisson_diff)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(real_grid_number - 1), intent(in) :: diff_coordinate_FA
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_poisson_diff

    integer :: count_i

    do count_i = 1, real_grid_number
        
        if ( count_i == 1 .or. count_i == real_grid_number ) then
            charge_density_poisson_diff(:, count_i) = 0d0

        else
            charge_density_poisson_diff(:, count_i) = electrostatic_potential_diff(1, count_i - 1) &
                & / diff_coordinate_FA(count_i - 1) / (diff_coordinate_FA(count_i - 1) + diff_coordinate_FA(count_i))
            
            charge_density_poisson_diff(:, count_i) = charge_density_poisson_diff(:, count_i) &
                & + electrostatic_potential_diff(1, count_i + 1) / diff_coordinate_FA(count_i) &
                & / (diff_coordinate_FA(count_i - 1) + diff_coordinate_FA(count_i))
            
            charge_density_poisson_diff(:, count_i) = charge_density_poisson_diff(:, count_i) &
                & - electrostatic_potential_diff(:, count_i) / diff_coordinate_FA(count_i - 1) / diff_coordinate_FA(count_i)

            charge_density_poisson_diff(:, count_i) = - 2d0 * electric_constant * charge_density_poisson_diff(:, count_i)

        end if

    end do  !count_i

end subroutine make_charge_density_from_Poisson_eq
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_convergence_number(charge_density_diff, charge_density_plus_diff, charge_density_minus_diff, &
    & charge_density_poisson_diff, convergence_number_diff, convergence_number_sum)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: charge_density_diff, charge_density_plus_diff
    double precision, dimension(3, real_grid_number), intent(in) :: charge_density_minus_diff, charge_density_poisson_diff
    double precision, dimension(3, real_grid_number), intent(out) :: convergence_number_diff
    double precision, intent(out) :: convergence_number_sum

    convergence_number_diff = (charge_density_diff - charge_density_poisson_diff)**2d0
    convergence_number_diff = convergence_number_diff / (- charge_density_minus_diff) / charge_density_plus_diff
    
    convergence_number_sum = sqrt(sum(convergence_number_diff(1, :)) / real_grid_number)
   
end subroutine make_convergence_number
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
