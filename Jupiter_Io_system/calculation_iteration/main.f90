program main

    use constant_parameter
    use constant_in_the_simulation
    use main_variables

    implicit none

    !-------------------------
    ! initial setting of field
    !-------------------------

    !MLAT
    do count_i = 1, real_grid_number

        if ( count_i <= 41 ) then
            mlat(count_i) = planet_mlat_1 + (planet_mlat_2 - planet_mlat_1) * dble(count_i - 1) / dble(1600)

        else if ( count_i >= 42 .and. count_i <= 193 ) then
            mlat(count_i) = planet_mlat_1 + (planet_mlat_2 - planet_mlat_1) * dble(count_i - 41 + 4) / dble(160)

        else if ( count_i >= 194 ) then
            mlat(count_i) = planet_mlat_1 + (planet_mlat_2 - planet_mlat_1) * dble(count_i - 193 + 1560) / dble(1600)

        end if

    end do  !count_i

    !length2planet
    length2planet = planet_l_shell * planet_radius * cos(mlat)**2d0

    !length2satellite
    if ( satellite_mass /= 0d0 ) then
        length2satellite =  sqrt(length2planet**2d0 + (planet_radius * planet_l_shell)**2d0 &
            & - 2d0 * length2planet * planet_radius * planet_l_shell * cos(mlat))
    end if

    !coordinate_FA
    coordinate_FA = planet_l_shell * planet_radius * (sin(mlat) * sqrt(1d0 + 3d0 * sin(mlat)**2d0) / 2d0 &
        & + asinh(sqrt(3d0) * sin(mlat)) / 2d0 / sqrt(3d0))
    
    !diff_coordinate_FA
    do count_i = 1, real_grid_number - 1
        diff_coordinate_FA(count_i) = coordinate_FA(count_i + 1) - coordinate_FA(count_i)
    end do  !count_i

    !magnetic_flux_density
    magnetic_flux_density = magnetic_constant * planet_dipole_moment / 4d0 / pi / (planet_l_shell * planet_radius)**3d0 &
        & * sqrt(1d0 + 3d0 * sin(mlat)**2d0) / cos(mlat)**6d0

    
    !------------------
    ! initial condition
    !------------------
    
    do count_i = 1, real_grid_number

        if ( (count_i <= initial_grid_ionophere_middle_1 - 1) &
            & .or. (initial_grid_ionophere_middle_2 <= count_i .and. count_i <= real_grid_number) ) then
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_ionosphere

        else if ( (initial_grid_ionophere_middle_1 <= count_i .and. count_i <= initial_grid_middle_magnetosphere_1 - 1) &
            & .or. (initial_grid_middle_magnetosphere_2 <= count_i .and. count_i <= initial_grid_ionophere_middle_2 - 1) ) then
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_middle
        
        else if ( initial_grid_middle_magnetosphere_1 <= count_i .and. count_i <= initial_grid_middle_magnetosphere_2 - 1 ) then
            initial_electrostatic_potential(count_i) = initial_electrostatic_potential_magnetosphere
        
        end if

    end do  !count_i

    electrostatic_potential = initial_electrostatic_potential


    !--------------------
    ! boundary conditions
    !--------------------

    open(30, file = boundary_file, action = 'read', status = 'old')
    read(30, *)
    do count_s = 1, boundary_series_number

        read(30, *) dummy, boundary_number_density(count_s), boundary_temperature_perp(count_s), &
            & boundary_temperature_para(count_s), charge_number(count_s), particle_mass(count_s), injection_grid_number(count_s)
    
    end do  !count_s
    close(30)

    charge_number = charge_number * elementary_charge
    boundary_temperature_perp = boundary_temperature_perp * elementary_charge
    boundary_temperature_para = boundary_temperature_para * elementary_charge


    !------------------------
    ! 1st adiabatic invariant
    !------------------------

    call make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)


    !----------------
    ! iteration start
    !----------------

    count_iteration = 0

    do  !count_iteration
        count_iteration = count_iteration + 1


        !-----------
        ! make _diff
        !-----------

        call make_electrostatic_potential_diff(electrostatic_potential, electrostatic_potential_diff)

        call make_boundary_number_density_diff(boundary_number_density, boundary_number_density_diff)

        call make_potential_energy_diff(mlat, length2planet, length2satellite, charge_number, particle_mass, &
            & electrostatic_potential_diff, potential_energy_diff)

        
        !--------------
        ! accessibility
        !--------------
        
        call make_amin(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, injection_grid_number, amin)

        call make_amax(potential_energy_diff, adiabatic_invariant, magnetic_flux_density, injection_grid_number, particle_mass, &
        & amin, amax)

        !--------
        ! density
        !--------

        call make_number_density(boundary_number_density_diff, boundary_temperature_perp, boundary_temperature_para, &
            & potential_energy_diff, magnetic_flux_density, adiabatic_invariant, injection_grid_number, particle_mass, amin, &
            & amax, number_density_diff)

        call make_charge_density_from_number_density(number_density_diff, charge_number, charge_density_diff, &
            & charge_density_plus_diff, charge_density_minus_diff)

        call make_charge_density_from_Poisson_eq(electrostatic_potential_diff, diff_coordinate_FA, charge_density_poisson_diff)

        
        !-------
        ! finish
        !-------

        !NaN check
        do count_h = 1, 3
            do count_s = 1, boundary_series_number
                do count_i = 1, real_grid_number
                    do count_mu = 1, adiabatic_invariant_grid_number
                        if ( charge_density_poisson_diff(count_h, count_i) /= charge_density_poisson_diff(count_h, count_i) ) then
                        print *, "NaN", count_h, count_s, count_i, count_mu
                        end if
                    end do
                end do
            end do
        end do

        print *, count_iteration

        if (count_iteration == 3) then
            exit
        end if 

    end do !count_iteration





























end program main