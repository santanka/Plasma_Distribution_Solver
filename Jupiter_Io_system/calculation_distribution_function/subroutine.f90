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
subroutine make_result_file_name(count_s, result_file_name)
    use reference_results_setting

    implicit none
    
    integer, intent(in) :: count_s
    character(len = 106), intent(out) :: result_file_name

    character(len = 2) :: count_s_character
    character(len = 3) :: plot_grid_number_character

    write(count_s_character, "(I2.2)") count_s
    write(plot_grid_number_character, "(I3.3)") plot_grid_number

    result_file_name = result_file_front // result_file_back // '_grid_' // plot_grid_number_character // '_series_' &
        & // count_s_character // '.csv'
    
end subroutine make_result_file_name
