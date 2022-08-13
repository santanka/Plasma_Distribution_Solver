subroutine make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

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
subroutine make_potential_energy(mlat_rad, length2planet, length2satellite, charge_number, particle_mass, &
    & electrostatic_potential, potential_energy)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: mlat_rad, length2planet, length2satellite
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(real_grid_number), intent(in) :: electrostatic_potential
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) ::  potential_energy

    integer :: count_i

    !$omp parallel
    !$omp do
    do count_i = 1, real_grid_number
        
        !gravity of planet
        potential_energy(:, count_i) = - constant_of_gravitation * planet_mass * particle_mass / length2planet(count_i)

        !centrifugal force of planet
        potential_energy(:, count_i) = potential_energy(:, count_i) &
            & - particle_mass * (planet_rotation * length2planet(count_i) * cos(mlat_rad(count_i)))**2d0 / 2d0

        !gravity of satellite
        if ( satellite_mass /= 0d0 ) then
            potential_energy(:, count_i) = potential_energy(:, count_i) &
                & - constant_of_gravitation * satellite_mass * particle_mass / length2satellite(count_i)
        end if
        
        !Coulomb force
        potential_energy(:, count_i) = potential_energy(:, count_i) &
            & + charge_number * electrostatic_potential(count_i)

    end do  !count_i
    !$omp end do
    !$omp end parallel

end subroutine make_potential_energy
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_potential_plus_Bmu(potential_energy, adiabatic_invariant, magnetic_flux_density, potential_plus_Bmu)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) ::  potential_energy
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: &
        & potential_plus_Bmu

    integer :: count_s, count_i

    do count_s = 1, boundary_series_number

        !$omp parallel
        !$omp do
        do count_i = 1, real_grid_number
            
            potential_plus_Bmu(count_s, count_i, :) = potential_energy(count_s, count_i) &
                & + magnetic_flux_density(count_i) * adiabatic_invariant(count_s, :)

        end do  !count_i
        !$omp end do
        !$omp end parallel
        
    end do  !count_s

end subroutine make_potential_plus_Bmu
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amin(potential_plus_Bmu, injection_grid_number, amin)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amin

    integer :: count_s, count_i, count_mu, Emax_grid, count4max

    do count_s = 1, boundary_series_number
        
        !$omp parallel private(Emax_grid, count_mu, count4max)
        !$omp do
        do count_i = 1, real_grid_number
            
            do count_mu = 1, adiabatic_invariant_grid_number

                if ( count_i <= injection_grid_number(count_s) ) then
                    Emax_grid = count_i
                    do count4max = count_i, injection_grid_number(count_s)
                        
                        if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                            & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                            Emax_grid = count4max
                        end if
                    
                    end do  !count4max
                
                else if ( count_i > injection_grid_number(count_s) ) then
                    Emax_grid = injection_grid_number(count_s)
                    do count4max = injection_grid_number(count_s), count_i

                        if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                            & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                        Emax_grid = count4max
                        end if
                    
                    end do  !count4max

                end if

                amin(count_s, count_i, count_mu) = sqrt(potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                    & - potential_plus_Bmu(count_s, injection_grid_number(count_s), count_mu))
                
            end do  !count_mu

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
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: alim

    integer :: count_s, count_i, count_mu

    do count_s = 1, boundary_series_number

        alim(count_s, :, :) = sqrt(5d-1 * particle_mass(count_s)) * speed_of_light
        
        !$omp parallel private(count_mu)
        !$omp do
        do count_i = 1, real_grid_number
            
            do count_mu = 1, adiabatic_invariant_grid_number
                
                if ( alim(count_s, count_i, count_mu) <= amin(count_s, count_i, count_mu) ) then
                    alim(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)
                end if

            end do  !count_mu

        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s

end subroutine make_alim
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_amax(potential_plus_Bmu, injection_grid_number, particle_mass, amin, amax)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(out) :: amax

    integer :: count_s, count_i, count_mu, Emax_grid, count4max, injection_grid
    double precision :: energy_difference

    do count_s = 1, boundary_series_number
        
        !$omp parallel private(count_mu, count4max, Emax_grid, energy_difference, injection_grid)
        !$omp do
        do count_i = 1, real_grid_number
            
            if ( (count_i == 1 .and. injection_grid_number(count_s) /= 1) &
                & .or. (count_i == real_grid_number .and. injection_grid_number(count_s) /= real_grid_number) ) then
                amax(count_s, count_i, :) = amin(count_s, count_i, :)

            else if ( count_i == injection_grid_number(count_s) .and. injection_grid_number(count_s) /= 1 &
                & .and. injection_grid_number(count_s) /= real_grid_number ) then
                amax(count_s, count_i, :) = sqrt(particle_mass(count_s) / 2d0) * speed_of_light
            
            else
                do count_mu = 1, adiabatic_invariant_grid_number
                    
                    if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                        Emax_grid = 1
                        do count4max = 1, count_i - 1

                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                                & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max
                    
                    else if (count_i > injection_grid_number(count_s) .or. injection_grid_number(count_s) == 1) then
                        Emax_grid = count_i + 1
                        do count4max = count_i + 1, real_grid_number
                            
                            if ( potential_plus_Bmu(count_s, count4max, count_mu) &
                            & > potential_plus_Bmu(count_s, Emax_grid, count_mu) ) then
                                Emax_grid = count4max
                            end if

                        end do  !count4max

                    end if

                    injection_grid = injection_grid_number(count_s)
                    energy_difference = potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu(count_s, injection_grid, count_mu)

                    if ( energy_difference <= amin(count_s, count_i, count_mu)**2d0 ) then
                        amax(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)

                    else if ( energy_difference > amin(count_s, count_i, count_mu)**2d0 ) then

                        if ( energy_difference < 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
                            amax(count_s, count_i, count_mu) = sqrt(energy_difference)

                        else if ( energy_difference >= 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
                            amax(count_s, count_i, count_mu) = sqrt(5d-1 * particle_mass(count_s)) * speed_of_light

                        end if

                    end if

                    if ( amax(count_s, count_i, count_mu) < amin(count_s, count_i, count_mu) ) then
                        amax(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)
                    end if

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
subroutine make_number_density(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & magnetic_flux_density, adiabatic_invariant, injection_grid_number, amin, alim, amax, &
    & number_density)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number) :: amin, alim, amax
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: number_density

    integer :: count_s, count_i, injection_grid
    double precision :: integral, sqrt_temperature_para
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax, coefficient4integral
        
    do count_s = 1, boundary_series_number

        !$omp parallel private(integral, sqrt_temperature_para, xmin, xlim, xmax, injection_grid, coefficient4integral)
        !$omp do
        do count_i = 1, real_grid_number
            
            integral = 0d0
            sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))
            
            xmin = amin(count_s, count_i, :) / sqrt_temperature_para
            xlim = alim(count_s, count_i, :) / sqrt_temperature_para
            xmax = amax(count_s, count_i, :) / sqrt_temperature_para

            injection_grid = injection_grid_number(count_s)

            coefficient4integral = magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, :) &
                & / boundary_temperature_perp(count_s)

            call calculation_exp_erf(xmin, xlim, xmax, coefficient4integral, adiabatic_invariant(count_s, :), integral)

            number_density(count_s, count_i) = boundary_number_density(count_s) &
                & * magnetic_flux_density(count_i) / boundary_temperature_perp(count_s) / 2d0 * integral

        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s
    
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
subroutine cannot_reach_check(number_density, injection_grid_number)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(inout) :: number_density
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number

    integer :: cannot_reach_point
    integer :: count_s, count_i

    !$omp parallel private(count_i, cannot_reach_point)
    !$omp do
    do count_s = 1, boundary_series_number
            
        if ( injection_grid_number(count_s) /= real_grid_number ) then
            cannot_reach_point = 0
            do count_i = injection_grid_number(count_s), real_grid_number

                if( number_density(count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                    number_density(count_s, count_i) = 0d0
                    cannot_reach_point = count_i
                end if

                if ( cannot_reach_point /= 0 .and. count_i >= cannot_reach_point ) then
                    number_density(count_s, count_i) = 0d0
                end if
                    
            end do  !count_i
        end if

        if ( injection_grid_number(count_s) /= 1 ) then
            cannot_reach_point = 0
            do count_i = injection_grid_number(count_s), 1, -1

                if( number_density(count_s, count_i) < 1d-5 .and. cannot_reach_point == 0) then
                    number_density(count_s, count_i) = 0d0
                    cannot_reach_point = count_i
                end if

                if ( cannot_reach_point /= 0 .and. count_i <= cannot_reach_point ) then
                    number_density(count_s, count_i) = 0d0
                end if
                    
            end do  !count_i
        end if

    end do  !count_s
    !$omp end do
    !$omp end parallel

end subroutine cannot_reach_check
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_particle_flux_density(boundary_number_density, magnetic_flux_density, injection_grid_number, &
    & boundary_temperature_perp, boundary_temperature_para, particle_mass, adiabatic_invariant, potential_plus_Bmu, &
    & alim, amax, number_density, particle_flux_density, parallel_mean_velocity)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density, particle_mass
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu, alim, amax
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: particle_flux_density
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: parallel_mean_velocity

    integer :: count_s, count_i, injection_grid
    double precision :: integral, sqrt_temperature_para, flux_direction
    double precision, dimension(adiabatic_invariant_grid_number) :: xlim, xmax, coefficient4integral, y_mu

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number

            if ( number_density(count_s, count_i) /= 0d0 ) then
                injection_grid = injection_grid_number(count_s)

                if ( count_i < injection_grid_number(count_s) .or. injection_grid_number(count_s) == real_grid_number ) then
                    flux_direction = - 1d0

                else
                    flux_direction = 1d0

                end if
            
                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))

                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                y_mu = magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, :) &
                    & / boundary_temperature_perp(count_s)

                coefficient4integral = (potential_plus_Bmu(count_s, injection_grid, :) &
                    & - potential_plus_Bmu(count_s, count_i, :)) / boundary_temperature_para(count_s)

                call calculation_sqrt_exp(xlim, xmax, y_mu, coefficient4integral, integral)

                particle_flux_density(count_s, count_i) = boundary_number_density(count_s) / sqrt(pi) &
                    & * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid) &
                    & * sqrt(2d0 * boundary_temperature_para(count_s) / particle_mass(count_s)) * integral * flux_direction
                
                parallel_mean_velocity(count_s, count_i) = particle_flux_density(count_s, count_i) /number_density(count_s, count_i)
                
            else if ( number_density(count_s, count_i) /= 0d0 ) then
                particle_flux_density(count_s, count_i) = 0d0
                parallel_mean_velocity(count_s, count_i) = 0d0

            end if

        end do  !count_i

    end do  !count_s

end subroutine make_particle_flux_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_sqrt_exp(xlim, xmax, y_mu, coefficient4integral, integral_result)
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xlim, xmax, coefficient4integral, y_mu
    double precision, intent(out) :: integral_result

    integer :: count_y, count_x
    double precision, dimension(adiabatic_invariant_grid_number) :: result_halfway, xmax2xlim
    double precision :: integral_former_x, integral_latter_x

    integral_result = 0d0
    result_halfway = 0d0

    !$omp parallel private(xmax2xlim, count_x, integral_former_x, integral_latter_x)
    !$omp do
    do count_y = 1, adiabatic_invariant_grid_number
        
        if ( xlim(count_y) > xmax(count_y) ) then
            if ( xmax(count_y) == 0d0 ) then
                xmax2xlim(1) = 0d0
                do count_x = 2, adiabatic_invariant_grid_number
                    
                    xmax2xlim(count_x) = xlim(count_y)*1d-5 * (xlim(count_y) / (xlim(count_y)*1d-5)) &
                        & **(dble(count_x - 2) / dble(adiabatic_invariant_grid_number -2))

                end do  !count_x
            
            else if ( xmax(count_y) /= 0d0 ) then
                do count_x = 1, adiabatic_invariant_grid_number
                    
                    xmax2xlim(count_x) = xmax(count_y) * (xlim(count_y) / xmax(count_y)) &
                        & **(dble(count_x - 1) / dble(adiabatic_invariant_grid_number - 1))

                end do  !count_x

            end if

            do count_x = 1, adiabatic_invariant_grid_number - 1

                integral_former_x = sqrt(abs(xmax2xlim(count_x)**2d0 + coefficient4integral(count_y))) &
                    & * exp(- xmax2xlim(count_x)**2d0) * exp(- y_mu(count_y))
                    
                integral_latter_x = sqrt(abs(xmax2xlim(count_x + 1)**2d0 + coefficient4integral(count_y))) &
                    & * exp(- xmax2xlim(count_x + 1)**2d0) * exp(- y_mu(count_y))
                
                result_halfway(count_y) = result_halfway(count_y) + (integral_former_x + integral_latter_x) / 2d0 &
                    & * (xmax2xlim(count_x + 1) - xmax2xlim(count_x))
                
            end do  !count_x

        end if

    end do  !count_y
    !$omp end do
    !$omp end parallel

    do count_y = 1, adiabatic_invariant_grid_number - 1
        
        integral_result = integral_result + (result_halfway(count_y) &
            & + result_halfway(count_y + 1)) / 2d0 * (y_mu(count_y + 1) - y_mu(count_y))

    end do  !count_y

end subroutine calculation_sqrt_exp
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_perpendicular(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & injection_grid_number, magnetic_flux_density, adiabatic_invariant, number_density, amin, alim, amax, pressure_perp, &
    & temperature_perp)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number
    !$use omp_lib

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin, alim
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amax
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_perp, temperature_perp

    integer :: count_s, count_i, injection_grid
    double precision :: sqrt_temperature_para, integral
    double precision, dimension(adiabatic_invariant_grid_number) :: xmin, xlim, xmax, y_mu

    do count_s = 1, boundary_series_number
        
        !$omp parallel private(sqrt_temperature_para, xmin, xlim, xmax, y_mu, injection_grid, integral)
        !$omp do
        do count_i = 1, real_grid_number
            
            if ( number_density(count_s, count_i) /= 0d0 ) then

                integral = 0d0
            
                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))

                xmin = amin(count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                injection_grid = injection_grid_number(count_s)

                y_mu = magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, :) / boundary_temperature_perp(count_s)

                call calculation_y_exp_erf(xmin, xlim, xmax, y_mu, integral)

                pressure_perp(count_s, count_i) = boundary_number_density(count_s) * boundary_temperature_perp(count_s) / 2d0 &
                    & * (magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid))**2d0 * integral

                temperature_perp(count_s, count_i) = pressure_perp(count_s, count_i) / number_density(count_s, count_i)
                
            else
                pressure_perp(count_s, count_i) = 0d0
                temperature_perp(count_s, count_i) = 0d0

            end if

        end do  !count_i
        !$omp end do
        !$omp end parallel

    end do  !count_s

end subroutine make_pressure_perpendicular
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_y_exp_erf(xmin, xlim, xmax, y_mu, integral_result)
    use constant_in_the_simulation

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xmin, xlim, xmax, y_mu
    double precision, intent(out) :: integral_result

    double precision, dimension(adiabatic_invariant_grid_number) :: y_exp_erf
    integer :: count_y

    y_exp_erf = y_mu * exp(- y_mu) * (erf(xlim) + erf(xmax) - 2d0 * erf(xmin))

    integral_result = 0d0

    do count_y = 1, adiabatic_invariant_grid_number - 1
        
        integral_result = integral_result + (y_exp_erf(count_y) + y_exp_erf(count_y + 1)) / 2d0 &
            & * (y_mu(count_y + 1) - y_mu(count_y))
        
    end do  !count_y
    
end subroutine calculation_y_exp_erf
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_parallel(boundary_number_density, boundary_temperature_perp, boundary_temperature_para, &
    & injection_grid_number, magnetic_flux_density, adiabatic_invariant, potential_plus_Bmu, particle_mass, number_density, &
    & parallel_mean_velocity, amin, alim, amax, pressure_para, temperature_para)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: boundary_number_density, particle_mass
    double precision, dimension(boundary_series_number), intent(in) :: boundary_temperature_perp, boundary_temperature_para
    integer, dimension(boundary_series_number), intent(in) :: injection_grid_number
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number, adiabatic_invariant_grid_number), intent(in) :: adiabatic_invariant
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: &
        & potential_plus_Bmu
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, parallel_mean_velocity
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amin, alim
    double precision, dimension(boundary_series_number, real_grid_number, adiabatic_invariant_grid_number), intent(in) :: amax
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_para, temperature_para

    integer :: count_s, count_i, injection_grid
    double precision :: sqrt_temperature_para, integral
    double precision, dimension(adiabatic_invariant_grid_number) :: y_mu, coefficient_1, coefficient_2, xlim, xmax, xmin

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number
            
            if ( number_density(count_s, count_i) /= 0d0 ) then

                integral = 0d0
            
                sqrt_temperature_para = sqrt(boundary_temperature_para(count_s))

                xmin = amin(count_s, count_i, :) / sqrt_temperature_para
                xlim = alim(count_s, count_i, :) / sqrt_temperature_para
                xmax = amax(count_s, count_i, :) / sqrt_temperature_para

                injection_grid = injection_grid_number(count_s)

                y_mu = magnetic_flux_density(injection_grid) * adiabatic_invariant(count_s, :) / boundary_temperature_perp(count_s)

                coefficient_1 = (potential_plus_Bmu(count_s, injection_grid, :) - potential_plus_Bmu(count_s, count_i, :)) &
                    & / boundary_temperature_para(count_s)

                coefficient_2 = sqrt(particle_mass(count_s) / 2d0 / boundary_temperature_para(count_s)) &
                    & * abs(parallel_mean_velocity(count_s, count_i))

                call calculation_integral_for_parallel_pressure(xmin, xlim, xmax, y_mu, coefficient_1, coefficient_2, integral)

                pressure_para(count_s, count_i) = 2d0 / sqrt(pi) * boundary_number_density(count_s) &
                    & * boundary_temperature_para(count_s) * magnetic_flux_density(count_i) / magnetic_flux_density(injection_grid)&
                    & * integral

                temperature_para(count_s, count_i) = pressure_para(count_s, count_i) / number_density(count_s, count_i)
                
            else
                pressure_para(count_s, count_i) = 0d0
                temperature_para(count_s, count_i) = 0d0

            end if

        end do  !count_i
        
    end do  !count_s

end subroutine make_pressure_parallel
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine calculation_integral_for_parallel_pressure(xmin, xlim, xmax, y_mu, coefficient_1, coefficient_2, integral_result)
    use constant_in_the_simulation
    !$use omp_lib

    implicit none
    
    double precision, dimension(adiabatic_invariant_grid_number), intent(in) :: xmin, xlim, xmax, y_mu, coefficient_1, coefficient_2
    double precision, intent(out) :: integral_result

    integer :: count_y, count_x
    double precision, dimension(adiabatic_invariant_grid_number) :: result_halfway, xmin2xlim, xmin2xmax
    double precision :: integral_former_x, integral_latter_x

    integral_result = 0d0
    result_halfway = 0d0

    !$omp parallel private(xmin2xmax, count_x, integral_former_x, integral_latter_x, xmin2xlim)
    !$omp do
    do count_y = 1, adiabatic_invariant_grid_number

        if ( xmax(count_y) > xmin(count_y) ) then
            if ( xmin(count_y) == 0d0 ) then
                xmin2xmax(1) = 0d0
                do count_x = 2, adiabatic_invariant_grid_number
                    
                    xmin2xmax(count_x) = xmax(count_y)*1d-5 * (xmax(count_y) / (xmax(count_y)*1d-5)) &
                        & **(dble(count_x - 2) / dble(adiabatic_invariant_grid_number -2))

                end do  !count_x
            
            else if ( xmin(count_y) /= 0d0 ) then
                do count_x = 1, adiabatic_invariant_grid_number
                    
                    xmin2xmax(count_x) = xmin(count_y) * (xmax(count_y) / xmin(count_y)) &
                        & **(dble(count_x - 1) / dble(adiabatic_invariant_grid_number - 1))

                end do  !count_x

            end if

            do count_x = 1, adiabatic_invariant_grid_number - 1

                integral_former_x = (sqrt(abs(xmin2xmax(count_x)**2d0 + coefficient_1(count_y))) + coefficient_2(count_y))**2d0 &
                    & * exp(- xmin2xmax(count_x)**2d0 - y_mu(count_y))
                    
                integral_latter_x = (sqrt(abs(xmin2xmax(count_x + 1)**2d0 + coefficient_1(count_y))) + coefficient_2(count_y))**2d0&
                & * exp(- xmin2xmax(count_x + 1)**2d0 - y_mu(count_y))
                
                result_halfway(count_y) = result_halfway(count_y) + (integral_former_x + integral_latter_x) / 2d0 &
                    & * (xmin2xmax(count_x + 1) - xmin2xmax(count_x))
                
            end do  !count_x

        end if

        if ( xlim(count_y) > xmin(count_y) ) then
            if ( xmin(count_y) == 0d0 ) then
                xmin2xlim(1) = 0d0
                do count_x = 2, adiabatic_invariant_grid_number
                    
                    xmin2xlim(count_x) = xlim(count_y)*1d-5 * (xlim(count_y) / (xlim(count_y)*1d-5)) &
                        & **(dble(count_x - 2) / dble(adiabatic_invariant_grid_number -2))
                    if ( xmin2xlim(count_x) /= xmin2xlim(count_x) ) print *, count_y, count_x

                end do  !count_x
            
            else if ( xmin(count_y) /= 0d0 ) then
                do count_x = 1, adiabatic_invariant_grid_number
                    
                    xmin2xlim(count_x) = xmin(count_y) * (xlim(count_y) / xmin(count_y)) &
                        & **(dble(count_x - 1) / dble(adiabatic_invariant_grid_number - 1))
                    if ( xmin2xlim(count_x) /= xmin2xlim(count_x) ) print *, count_y, count_x

                end do  !count_x

            end if

            do count_x = 1, adiabatic_invariant_grid_number - 1

                integral_former_x = (sqrt(abs(xmin2xlim(count_x)**2d0 + coefficient_1(count_y))) - coefficient_2(count_y))**2d0 &
                    & * exp(- xmin2xlim(count_x)**2d0 - y_mu(count_y))
                    
                integral_latter_x = (sqrt(abs(xmin2xlim(count_x + 1)**2d0 + coefficient_1(count_y))) - coefficient_2(count_y))**2d0&
                    & * exp(- xmin2xlim(count_x + 1)**2d0 - y_mu(count_y))
                
                result_halfway(count_y) = result_halfway(count_y) + (integral_former_x + integral_latter_x) / 2d0 &
                    & * (xmin2xlim(count_x + 1) - xmin2xlim(count_x))
                
            end do  !count_x

        end if
        
    end do  !count_y
    !$omp end do
    !$omp end parallel

    do count_y = 1, adiabatic_invariant_grid_number - 1
        
        integral_result = integral_result + (result_halfway(count_y) &
            & + result_halfway(count_y + 1)) / 2d0 * (y_mu(count_y + 1) - y_mu(count_y))

    end do  !count_y

end subroutine calculation_integral_for_parallel_pressure
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_pressure_dynamic(number_density, parallel_mean_velocity, particle_mass, pressure_dynamic)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, parallel_mean_velocity
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: pressure_dynamic

    integer :: count_s

    do count_s = 1, boundary_series_number
        
        pressure_dynamic(count_s, :) = 5d-1 * number_density(count_s, :) * particle_mass(count_s) &
            & * parallel_mean_velocity(count_s, :)**2d0

    end do  !count_s

    
end subroutine make_pressure_dynamic
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_Alfven_speed(magnetic_flux_density, particle_mass, number_density, Alfven_speed, Alfven_speed_per_lightspeed)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number), intent(in) :: particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(real_grid_number), intent(out) :: Alfven_speed, Alfven_speed_per_lightspeed

    double precision, dimension(real_grid_number) :: mass_density
    integer :: count_i

    do count_i = 1, real_grid_number
        
        mass_density(count_i) = sum(particle_mass * number_density(:, count_i))

    end do  !count_i

    Alfven_speed_per_lightspeed = magnetic_flux_density / sqrt(magnetic_flux_density**2d0 + mass_density / electric_constant)
    Alfven_speed = Alfven_speed_per_lightspeed * speed_of_light
    
end subroutine make_Alfven_speed
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_inertial_length(charge_number, particle_mass, number_density, ion_inertial_length, electron_inertial_length)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none

    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density
    double precision, dimension(real_grid_number), intent(out) :: ion_inertial_length, electron_inertial_length

    integer :: count_s
    double precision, dimension(real_grid_number) :: charge_density_ion, mass_density_ion, number_density_electron

    charge_density_ion = 0d0
    mass_density_ion = 0d0
    number_density_electron = 0d0

    do count_s = 1, boundary_series_number
        
        if ( charge_number(count_s) > 0d0 ) then
            charge_density_ion = charge_density_ion + charge_number(count_s) * number_density(count_s, :)
            mass_density_ion = mass_density_ion + particle_mass(count_s) * number_density(count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            number_density_electron = number_density_electron + number_density(count_s, :)

        end if

    end do  !count_s

    ion_inertial_length = sqrt(mass_density_ion / magnetic_constant) / charge_density_ion
    electron_inertial_length = sqrt(electron_mass / magnetic_constant / number_density_electron) / elementary_charge
    
end subroutine make_inertial_length
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_Larmor_radius(magnetic_flux_density, number_density, pressure_perp, charge_number, particle_mass, &
    & ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius)
    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none

    double precision, dimension(real_grid_number), intent(in) :: magnetic_flux_density
    double precision, dimension(boundary_series_number), intent(in) :: charge_number, particle_mass
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: number_density, pressure_perp
    double precision, dimension(real_grid_number), intent(out) :: ion_Larmor_radius, ion_acoustic_radius, electron_Larmor_radius

    integer :: count_s
    double precision, dimension(real_grid_number) :: pressure_perp_ion, pressure_perp_electron, mass_density_ion
    double precision, dimension(real_grid_number) :: number_density_electron, number_density_ion, charge_density_ion

    pressure_perp_ion = 0d0
    pressure_perp_electron = 0d0
    mass_density_ion = 0d0
    number_density_ion = 0d0
    number_density_electron = 0d0
    charge_density_ion = 0d0

    do count_s = 1, boundary_series_number
        
        if ( charge_number(count_s) > 0d0 ) then
            pressure_perp_ion = pressure_perp_ion + pressure_perp(count_s, :)
            mass_density_ion = mass_density_ion + particle_mass(count_s) * number_density(count_s, :)
            number_density_ion = number_density_ion + number_density(count_s, :)
            charge_density_ion = charge_density_ion + charge_number(count_s) * number_density (count_s, :)

        else if ( charge_number(count_s) < 0d0 ) then
            pressure_perp_electron = pressure_perp_electron + pressure_perp(count_s, :)
            number_density_electron = number_density_electron + number_density(count_s, :)

        end if

    end do  !count_s

    ion_Larmor_radius = sqrt(2d0 * pressure_perp_ion * mass_density_ion) / magnetic_flux_density / charge_density_ion
    ion_acoustic_radius = sqrt(2d0 * number_density_ion / number_density_electron * mass_density_ion * pressure_perp_electron) &
        & / magnetic_flux_density / charge_density_ion
    electron_Larmor_radius = sqrt(2d0 * pressure_perp_electron * number_density_electron * electron_mass) / magnetic_flux_density &
        & / number_density_electron / elementary_charge

end subroutine make_Larmor_radius
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_current_density(charge_number, particle_flux_density, current_density)
    use constant_in_the_simulation
    use reference_results_setting, only: boundary_series_number

    implicit none
    
    double precision, dimension(boundary_series_number), intent(in) :: charge_number
    double precision, dimension(boundary_series_number, real_grid_number), intent(in) :: particle_flux_density
    double precision, dimension(real_grid_number), intent(out) :: current_density
    
    integer :: count_s

    current_density = 0d0

    do count_s = 1, boundary_series_number
        
        current_density = current_density + charge_number(count_s) * particle_flux_density(count_s, :)

    end do  !count_s

end subroutine make_current_density
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
subroutine make_result_file_format(format_character)
    use reference_results_setting

    implicit none
    
    character(len = 34), intent(out) ::  format_character
    character(len = 3) :: series_number

    write(series_number, "(I3)") 8 * boundary_series_number + 17

    format_character = "(1PE25.15E3, " // series_number // "(',', 1PE25.15E3))"

end subroutine make_result_file_format