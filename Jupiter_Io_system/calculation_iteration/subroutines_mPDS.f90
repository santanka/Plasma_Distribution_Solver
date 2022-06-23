subroutine make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none

    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(out) :: adiabatic_invariant

    integer :: count_s, count_mu
    double precision :: max_mu

    do count_s = 1, boundary_series_number

        max_mu = particle_mass(count_s) * speed_of_light**2d0 / 2d0 / magnetic_flux_density(injection_grid_number(count_s))

        !$omp parallel
        !$omp do
        do count_mu = 1, adiabatic_invariant_grid_number

            if ( count_mu == 1 ) then
                adiabatic_invariant(count_s, count_mu) = 0d0
            else if ( count_mu /= 1 ) then
                adiabatic_invariant(count_s, count_mu) &
                    & = 1d-20 * (max_mu / 1d-20)**(dble(count_mu - 2) / dble(adiabatic_invariant_grid_number - 2))
            end if

        end do  !count_mu
        !$omp end do
        !$omp end parallel

    end do  !count_s
    
end subroutine make_adiabatic_invariant
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_electrostatic_potential_diff(electrostatic_potential,  electrostatic_potential_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: electrostatic_potential
    double precision, dimension(3, real_grid_number), intent(out) ::  electrostatic_potential_diff

    electrostatic_potential_diff(1, :) = electrostatic_potential
    electrostatic_potential_diff(2, :) = electrostatic_potential + 1d-6
    electrostatic_potential_diff(3, :) = electrostatic_potential - 1d-6
    electrostatic_potential_diff(2, 1) = electrostatic_potential(1)
    electrostatic_potential_diff(2, real_grid_number) = electrostatic_potential(real_grid_number)
    electrostatic_potential_diff(3, 1) = electrostatic_potential(1)
    electrostatic_potential_diff(3, real_grid_number) = electrostatic_potential(real_grid_number)
    if ( initial_fix_grid > initial_min_grid_1 .and. initial_fix_grid < initial_min_grid_2 ) then
        electrostatic_potential_diff(2, initial_fix_grid) = electrostatic_potential(initial_fix_grid)
        electrostatic_potential_diff(3, initial_fix_grid) = electrostatic_potential(initial_fix_grid)
    end if
        
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

    if ( 1 <= boundary_magnetosphere_variable_species .and. boundary_magnetosphere_variable_species <= boundary_series_number ) then
        boundary_number_density_diff(2, boundary_magnetosphere_variable_species) &
            & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 + 1d-7)
        boundary_number_density_diff(3, boundary_magnetosphere_variable_species) &
            & = boundary_number_density(boundary_magnetosphere_variable_species) * (1d0 - 1d-7)
    end if

end subroutine make_boundary_number_density_diff
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_energy_diff(mlat, length2planet, length2satellite, charge_number, particle_mass, &
    & electrostatic_potential_diff, potential_energy_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: mlat, length2planet, length2satellite
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) ::  potential_energy_diff

    integer :: count_h, count_i

    do count_h = 1, 3
        
        !$omp parallel
        !$omp do
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
        !$omp end do
        !$omp end parallel
    
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
    use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(in) :: potential_energy_diff
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) ::  magnetic_flux_density
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: &
        & potential_plus_Bmu_diff

    integer :: count_h, count_s, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            !$omp parallel
            !$omp do
            do count_i = 1, real_grid_number
                
                potential_plus_Bmu_diff(count_h, count_s, count_i, :) = potential_energy_diff(count_h, count_s, count_i) &
                    & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

            end do  !count_i
            !$omp end do
            !$omp end parallel
            
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
    use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max

    do count_s = 1, boundary_series_number

        !$omp parallel private(Emax_grid, count_h, count_mu, count4max)
        !$omp do
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
                            & > potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) ) then
                            amin(count_h, count_s, count_i, count_mu) &
                                & = sqrt(potential_plus_Bmu_diff(count_h, count_s, count_i, count_mu) &
                                & - potential_plus_Bmu_diff(1, count_s, injection_grid_number(count_s), count_mu))

                        else
                            amin(count_h, count_s, count_i, count_mu) &
                                & = sqrt(potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                                & - potential_plus_Bmu_diff(1, count_s, injection_grid_number(count_s), count_mu))

                        end if

                    end do  !count_mu

                end do  !count_h

            else if ( count_i == injection_grid_number(count_s) ) then
                amin(:, count_s, count_i, :) = 0d0

            end if
        end do  !count_i
        !$omp end do
        !$omp end parallel
        
    end do  !count_s

    
end subroutine make_amin
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_alim(particle_mass, amin, alim)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: alim

    integer :: count_h, count_s, count_mu, count_i

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            !$omp parallel private(count_i)
            !$omp do
            do count_mu = 1, adiabatic_invariant_grid_number
                
                alim(count_h, count_s, :, count_mu) = sqrt(5d-1 * particle_mass(count_s)) * speed_of_light

                do count_i = 1, real_grid_number

                    if ( alim(count_h, count_s, count_i, count_mu) <= amin(count_h, count_s, count_i, count_mu) ) then
                        alim(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)
                    end if
                
                end do  !count_i

            end do  !count_mu
            !$omp end do
            !$omp end parallel

        end do  !count_s

    end do  !count_h
    
end subroutine make_alim
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_plus_Bmu_diff, injection_grid_number, particle_mass, amin, amax)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_h, count_s, count_i, count_mu, Emax_grid, count4max, h_checker
    double precision :: energy_difference

    do count_s = 1, boundary_series_number

        !$omp parallel private(count_mu, count4max, Emax_grid, energy_difference, h_checker)
        !$omp do
        do count_i = 1, real_grid_number

            if ( (count_i == 1 .and. injection_grid_number(count_s) /= 1) &
                & .or. (count_i == real_grid_number .and. injection_grid_number(count_s) /= real_grid_number) ) then
                amax(:, count_s, count_i, :) = amin(:, count_s, count_i, :)

            else if ( count_i == injection_grid_number(count_s) .and. injection_grid_number(count_s) /= 1 &
                & .and. injection_grid_number(count_s) /= real_grid_number ) then
                amax(:, count_s, count_i, :) = sqrt(particle_mass(count_s) / 2d0) * speed_of_light
            
            else
                do count_mu = 1, adiabatic_invariant_grid_number

                    if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                        Emax_grid = 1
                        do count4max = 1, count_i - 1

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

                    do count_h = 1, 3

                        if ( count_i == injection_grid_number(count_s) .and. count_i == initial_fix_grid ) then
                            h_checker = count_h
                        else
                            h_checker = 1
                        end if

                        energy_difference = potential_plus_Bmu_diff(1, count_s, Emax_grid, count_mu) &
                            & - potential_plus_Bmu_diff(h_checker, count_s, injection_grid_number(count_s), count_mu)
                        
                        if ( energy_difference <= amin(count_h, count_s, count_i, count_mu)**2d0 ) then
                            amax(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)

                        else
                            if ( energy_difference < 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
                                amax(count_h, count_s, count_i, count_mu) = sqrt(energy_difference)

                            else
                                amax(count_h, count_s, count_i, count_mu) = sqrt(5d-1 * particle_mass(count_s)) * speed_of_light

                            end if
                            
                        end if

                        if ( amax(count_h, count_s, count_i, count_mu) < amin(count_h, count_s, count_i, count_mu) ) then
                            amax(count_h, count_s, count_i, count_mu) = amin(count_h, count_s, count_i, count_mu)
                        end if

                    end do  !count_h
                    
                end do  !count_mu

            end if
            
        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s

    
end subroutine make_amax
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
    & magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
    & number_density_diff)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(3, boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, alim, amax
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(out) :: number_density_diff

    integer :: count_h, count_s, count_i, injection_grid
    double precision :: integral, sqrt_temperature_para
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax, coefficient4integral

    do count_h = 1, 3
        
        do count_s = 1, boundary_series_number

            !$omp parallel private(integral, sqrt_temperature_para, xmin, xlim, xmax, injection_grid, coefficient4integral)
            !$omp do
            do count_i = 1, real_grid_number
                
                integral = 0d0

                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))
                
                xmin = amin(count_h, count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_h, count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_h, count_s, count_i, :) / sqrt_temperature_para

                injection_grid = injection_grid_number(count_s)
                coefficient4integral = magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, :) &
                    & / boundary_temperature_perp(count_s)

                call calculation_exp_erf(xmin, xlim, xmax, coefficient4integral, adiabatic_invariant(count_s, :), integral)

                if ( count_i == injection_grid_number(count_s) .and. (count_i == 1 .or. count_i == initial_fix_grid &
                    & .or. count_i == real_grid_number) ) then
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(count_h, count_s) &
                        & * magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) / 2d0 * integral
                
                else
                    number_density_diff(count_h, count_s, count_i) = boundary_number_density_diff(1, count_s) &
                        & * magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) / 2d0 * integral
                
                end if

            end do  !count_i
            !$omp end do
            !$omp end parallel

        end do  !count_s

    end do  !count_h
    
end subroutine make_number_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_exp_erf(xmin, xlim, xmax, coefficient4integral, adiabatic_invariant, integral_result)
    use constant_parameter
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: coefficient4integral, xmin, xlim, xmax
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, intent(out) :: integral_result

    integer :: count_mu
    double precision, dimension(adiabatic_invariant_grid_number) :: exp_erf

    integral_result = 0d0

    exp_erf = exp(- coefficient4integral) * (erf(xlim) + erf(xmax) - 2d0 * erf(xmin))

    do count_mu = 2, adiabatic_invariant_grid_number
        
        integral_result = integral_result + (exp_erf(count_mu - 1) + exp_erf(count_mu)) / 2d0 &
            & * (adiabatic_invariant(count_mu) - adiabatic_invariant(count_mu - 1))

    end do  !count_mu

end subroutine calculation_exp_erf
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine cannot_reach_check(number_density_diff, injection_grid_number)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(3, boundary_series_number, real_grid_number), intent(inout) :: number_density_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number

    integer :: cannot_reach_point
    integer :: count_h, count_s, count_i

    do count_h = 1, 3

        !$omp parallel private(count_i, cannot_reach_point)
        !$omp do
        do count_s = 1, boundary_series_number
            
            if ( injection_grid_number(count_s) /= real_grid_number ) then
                cannot_reach_point = 0
                do count_i = injection_grid_number(count_s), real_grid_number

                    if( number_density_diff(count_h, count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
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

                    if( number_density_diff(count_h, count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                        cannot_reach_point = count_i
                    end if

                    if ( cannot_reach_point /= 0 .and. count_i <= cannot_reach_point ) then
                        number_density_diff(count_h, count_s, count_i) = 0d0
                    end if
                    
                end do  !count_i
            end if

        end do  !count_s
        !$omp end do
        !$omp end parallel

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
    use omp_lib

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff
    double precision, dimension(real_grid_number - 1), intent(in) :: diff_coordinate_FA
    double precision, dimension(3, real_grid_number), intent(out) :: charge_density_poisson_diff

    integer :: count_i

    !$omp parallel private(count_i)
    !$omp do
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
    !$omp end do
    !$omp end parallel

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
    
    convergence_number_sum = sqrt(sum(convergence_number_diff(1, :)**2d0) / real_grid_number)
   
end subroutine make_convergence_number
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_name(result_file)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    character(len = 59), intent(out) ::  result_file
    character(len = 3) :: grid_front, grid_back, file_number, min_grid

    write(grid_front, "(I3.3)") initial_grid_ionophere_middle_1
    write(grid_back, "(I3.3)") initial_grid_middle_magnetosphere_1
    write(file_number, "(I3.3)") boundary_file_number
    write(min_grid, "(I3.3)") initial_min_grid_1

    result_file = result_file_front // "_" // grid_front // "_" // grid_back // "_BC_" // file_number &
        & // "_min_" // min_grid //".csv"
    
end subroutine make_result_file_name
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_format(format_character)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    character(len = 33), intent(out) ::  format_character
    character(len = 2) :: series_number

    write(series_number, "(I2.1)") boundary_series_number + 9

    format_character = "(1PE25.15E3, " // series_number // "(',', 1PE25.15E3))"

end subroutine make_result_file_format
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine Newton_method_for_electrostatic_potential(electrostatic_potential_diff, convergence_number_diff, electrostatic_potential)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions
    use omp_lib

    implicit none
    
    double precision, dimension(3, real_grid_number), intent(in) :: electrostatic_potential_diff, convergence_number_diff
    double precision, dimension(real_grid_number), intent(out) :: electrostatic_potential

    integer :: count_i
    double precision :: update
    double precision, dimension(initial_min_grid_1 - 1) :: sorting_potential_1, sorting_potential_1_reverse
    double precision, dimension(real_grid_number - initial_min_grid_2) :: sorting_potential_2
    double precision, dimension(initial_fix_grid - initial_min_grid_1 - 1) :: sorting_potential_3
    double precision, dimension(initial_min_grid_2 - initial_fix_grid - 1) :: sorting_potential_4, sorting_potential_4_reverse

    !----------
    ! iteration
    !----------

    electrostatic_potential = electrostatic_potential_diff(1, :)

    !$omp parallel private(update)
    !$omp do
    do count_i = 2, real_grid_number - 1
        
        if ( count_i /= initial_fix_grid .and. (convergence_number_diff(1, count_i) > convergence_number_diff(2, count_i) &
            & .or. convergence_number_diff(1, count_i) > convergence_number_diff(3, count_i)) ) then
            
            if ( convergence_number_diff(1, count_i) == convergence_number_diff(2, count_i) ) then
                update = (electrostatic_potential_diff(1, count_i) - electrostatic_potential_diff(3, count_i)) &
                    & / (convergence_number_diff(1, count_i) - convergence_number_diff(3, count_i)) &
                    & * convergence_number_diff(1, count_i)

            else if ( convergence_number_diff(1, count_i) == convergence_number_diff(3, count_i) ) then
                update = (electrostatic_potential_diff(2, count_i) - electrostatic_potential_diff(1, count_i)) &
                    & / (convergence_number_diff(2, count_i) - convergence_number_diff(1, count_i)) &
                    & * convergence_number_diff(1, count_i)

            else
                update = ((electrostatic_potential_diff(1, count_i) - electrostatic_potential_diff(3, count_i)) &
                    & / (convergence_number_diff(1, count_i) - convergence_number_diff(3, count_i)) &
                    & + (electrostatic_potential_diff(2, count_i) - electrostatic_potential_diff(1, count_i)) &
                    & / (convergence_number_diff(2, count_i) - convergence_number_diff(1, count_i))) / 2d0 &
                    & * convergence_number_diff(1, count_i)
                
            end if

            if ( abs(update) <= 1d1 ) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update
            
            else if ( abs(update) > 1d1 .and. sqrt(convergence_number_diff(1, count_i)) <= 2d1 ) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update / abs(update) &
                    & * sqrt(convergence_number_diff(1, count_i))

            else if ( abs(update) > 1d1 .and. sqrt(convergence_number_diff(1, count_i)) > 2d1 ) then
                electrostatic_potential(count_i) = electrostatic_potential_diff(1, count_i) - update / abs(update) * 2d1
                
            end if
            
        end if

    end do  !count_i
    !$omp end do
    !$omp end parallel


    !--------
    ! sorting
    !--------

    sorting_potential_1 = electrostatic_potential(2:initial_min_grid_1)
    do count_i = 1, initial_min_grid_1 - 1

        if ( sorting_potential_1(count_i) > electrostatic_potential(1) ) then
            sorting_potential_1(count_i) = 2d0 * electrostatic_potential(1) &
                & - sorting_potential_1(count_i) 
        end if

    end do  !count_i    
    call heapsort(initial_min_grid_1 - 1, sorting_potential_1)


    sorting_potential_2 = electrostatic_potential(initial_min_grid_2:real_grid_number-1)
    do count_i = 1, real_grid_number - initial_min_grid_2

        
        if ( sorting_potential_2(count_i) > electrostatic_potential(real_grid_number) ) then
            sorting_potential_2(count_i) = 2d0 * electrostatic_potential(real_grid_number) &
                & - sorting_potential_2(count_i) 
        end if

    end do

    call heapsort(real_grid_number - initial_min_grid_2, sorting_potential_2)

    do count_i = 1, initial_min_grid_1 - 1
        
        sorting_potential_1_reverse(count_i) = sorting_potential_1(initial_min_grid_1 - count_i)

    end do  !count_i

    if ( initial_fix_grid > initial_min_grid_1 .and. initial_fix_grid < initial_min_grid_2 ) then
        sorting_potential_3 = electrostatic_potential(initial_min_grid_1 + 1:initial_fix_grid - 1)
        sorting_potential_4 = electrostatic_potential(initial_fix_grid + 1:initial_min_grid_2 - 1)

        do count_i = 1, initial_fix_grid - initial_min_grid_1 - 1

            if ( sorting_potential_3(count_i) > electrostatic_potential(initial_fix_grid) ) then
                sorting_potential_3(count_i) = 2d0 * electrostatic_potential(initial_fix_grid) &
                    & - sorting_potential_3(count_i) 
            end if
            
        end do  !count_i

        do count_i = 1, initial_min_grid_2 - initial_fix_grid - 1

            if ( sorting_potential_4(count_i) > electrostatic_potential(initial_fix_grid) ) then
                sorting_potential_4(count_i) = 2d0 * electrostatic_potential(initial_fix_grid) &
                    & - sorting_potential_4(count_i) 
            end if
            
        end do  !count_i

        call heapsort(initial_fix_grid - initial_min_grid_1, sorting_potential_3)
        call heapsort(initial_min_grid_2 - initial_fix_grid, sorting_potential_4)
        
        do count_i = 1, initial_min_grid_2 - initial_fix_grid - 1
        
            sorting_potential_4_reverse(count_i) = sorting_potential_4(initial_min_grid_2 - initial_fix_grid - count_i)
    
        end do  !count_i

    end if

    !$omp parallel
    !$omp do
    do count_i = 1, real_grid_number
        
        if ( count_i > 1 .and. count_i <= initial_min_grid_1) then
            electrostatic_potential(count_i) = sorting_potential_1_reverse(count_i - 1)
        
        else if ( initial_min_grid_1 < count_i .and. count_i < initial_fix_grid .and. initial_fix_grid > initial_min_grid_1 &
            & .and. initial_fix_grid < initial_min_grid_2 ) then
            electrostatic_potential(count_i) = sorting_potential_3(count_i - initial_min_grid_1)
        
        else if ( initial_fix_grid < count_i .and. count_i < initial_min_grid_2 .and. initial_fix_grid > initial_min_grid_1 &
            & .and. initial_fix_grid < initial_min_grid_2 ) then
            electrostatic_potential(count_i) = sorting_potential_4_reverse(count_i - initial_fix_grid)
        
        else if ( count_i >= initial_min_grid_2 .and. count_i < real_grid_number ) then
            electrostatic_potential(count_i) = sorting_potential_2(count_i - initial_min_grid_2 + 1)

        end if

    end do  !count_i
    !$omp end do
    !$omp end parallel

end subroutine Newton_method_for_electrostatic_potential
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine heapsort(n,array)
    ! https://slpr.sakura.ne.jp/qp/sortf90/
    implicit none
    integer,intent(in) :: n
    double precision,intent(inout) :: array(1:n)
   
    integer ::i,k,j,l
    double precision :: t
   
    if(n.le.0)then
       write(6,*)"Error, at heapsort"; stop
    endif
    if(n.eq.1)return
  
    l=n/2+1
    k=n
    do while(k.ne.1)
       if(l.gt.1)then
          l=l-1
          t=array(L)
       else
          t=array(k)
          array(k)=array(1)
          k=k-1
          if(k.eq.1) then
             array(1)=t
             exit
          endif
       endif
       i=l
       j=l+l
       do while(j.le.k)
          if(j.lt.k)then
             if(array(j).lt.array(j+1))j=j+1
          endif
          if (t.lt.array(j))then
             array(i)=array(j)
             i=j
             j=j+j
          else
             j=k+1
          endif
       enddo
       array(i)=t
    enddo
  
    return
end subroutine heapsort
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine Newton_method_for_boundary_number_density(boundary_number_density_diff, convergence_number_diff, injection_grid_number, &
    & boundary_number_density)
    use constant_parameter
    use constant_in_the_simulation
    use boundary_and_initial_conditions

    implicit none
    
    double precision, dimension(3, boundary_series_number), intent(in) :: boundary_number_density_diff
    double precision, dimension(3, real_grid_number), intent(in) :: convergence_number_diff
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(out) :: boundary_number_density

    integer :: count_s, check_s, check_i
    double precision :: update

    boundary_number_density = boundary_number_density_diff(1, :)

    do count_s = 1, boundary_series_number

        if ( count_s == boundary_ionosphere_1_variable_species .or. count_s == boundary_ionosphere_2_variable_species &
            & .or. count_s == boundary_magnetosphere_variable_species ) then
            check_s = count_s
            check_i = injection_grid_number(check_s)
            
            if ( convergence_number_diff(1, check_i) > convergence_number_diff(2, check_i) &
            & .or. convergence_number_diff(1, check_i) > convergence_number_diff(3, check_i) ) then

                update = ((boundary_number_density_diff(1, check_s) - boundary_number_density_diff(3, check_s)) &
                    & / (convergence_number_diff(1, check_i) - convergence_number_diff(3, check_i)) &
                    & + (boundary_number_density_diff(2, check_s) - boundary_number_density_diff(1, check_s)) &
                    & / (convergence_number_diff(2, check_i) - convergence_number_diff(1, check_i))) / 2d0 &
                    & * convergence_number_diff(1, check_i)

                boundary_number_density(check_s) = boundary_number_density_diff(1, check_s) - update

                if ( boundary_number_density(check_s) < 0d0 ) then
                    print *, "boundary_number_density < 0d0"
                    stop
                end if
                
            end if

        end if

    end do  !count_i

    
end subroutine Newton_method_for_boundary_number_density