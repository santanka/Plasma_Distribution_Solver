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

    integer :: count_s, count_i, count_mu, Emax_grid, count4max
    double precision :: energy_difference

    do count_s = 2, boundary_series_number
        
        !$omp parallel private(count_mu, count4max, Emax_grid, energy_difference)
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

                    energy_difference = potential_plus_Bmu(count_s, Emax_grid, count_mu) &
                        & - potential_plus_Bmu(count_s, injection_grid_number(count_s), count_mu)

                    if ( energy_difference <= amin(count_s, count_i, count_mu)**2d0 ) then
                        amax(count_s, count_i, count_mu) = amin(count_s, count_i, count_mu)

                    else
                        if ( energy_difference < 5d-1 * particle_mass(count_s) * speed_of_light**2d0 ) then
                            amax(count_s, count_i, count_mu) = sqrt(energy_difference)

                        else
                            amax(count_s, count_i, count_mu) = sqrt(5d-1 * particle_mass(count_s)) * speed_of_light

                        end if

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
subroutine make_particle_flux_density(boundary_number_density, magnetic_flux_density, injection_grid_number, &
    & boundary_temperature_perp, boundary_temperature_para, particle_mass, adiabatic_invariant, potential_plus_Bmu, &
    & alim, amax, particle_flux_density)
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
    double precision, dimension(boundary_series_number, real_grid_number), intent(out) :: particle_flux_density

    integer :: count_s, count_i, injection_grid
    double precision :: integral, sqrt_temperature_para, flux_direction
    double precision, dimension(adiabatic_invariant_grid_number) :: xlim, xmax, coefficient4integral, y_mu

    do count_s = 1, boundary_series_number
        
        do count_i = 1, real_grid_number

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

            particle_flux_density(count_s, count_i) = boundary_number_density(count_s) / sqrt(pi) * magnetic_flux_density(count_i) &
                & / magnetic_flux_density(injection_grid) * sqrt(2d0 * boundary_temperature_para(count_s) &
                & / particle_mass(count_s)) * integral * flux_direction

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