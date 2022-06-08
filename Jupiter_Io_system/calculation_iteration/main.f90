program main

    use constant_parameter
    use constant_in_the_simulation
    use variables

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


    !NaN check
    do count_s = 1, boundary_series_number
        do count_mu = 1, adiabatic_invariant_grid_number
            if ( adiabatic_invariant(count_s, count_mu) /= adiabatic_invariant(count_s, count_mu) ) then
                print *, "NaN", count_s, count_mu
            end if
        end do
    end do





























end program main