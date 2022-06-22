program main

    use constant_parameter
    use constant_in_the_simulation
    use reference_results_setting
    use main_variables

    implicit none

    !-------------
    ! data loading
    !-------------

    open(20, file = reference_file, action = 'read', status = 'old')
    read(20, *) convergence_number_sum
    do count_i = 1, real_grid_number
        
        read(20, *) coordinate_FA(count_i), length2planet(count_i), mlat_rad(count_i), mlat_degree(count_i), &
            & magnetic_flux_density(count_i), initial_electrostatic_potential(count_i), electrostatic_potential(count_i), &
            & number_density(:, count_i), charge_density(count_i), charge_density_poisson(count_i), convergence_number(count_i)

    end do  !count_i

    do count_s = 1, boundary_series_number
        
        read(20, *) boundary_number_density(count_s), boundary_temperature_perp(count_s), boundary_temperature_para(count_s), &
            & charge_number(count_s), particle_mass(count_s), injection_grid_number(count_s)

    end do  !count_s
    close(20)

    charge_number = charge_number * elementary_charge
    boundary_temperature_perp = boundary_temperature_perp * elementary_charge
    boundary_temperature_para = boundary_temperature_para * elementary_charge


    !----------
    ! satellite
    !----------

    if ( satellite_mass /= 0d0 ) then
        length2satellite =  sqrt(length2planet**2d0 + (planet_radius * satellite_l_shell)**2d0 &
            & - 2d0 * length2planet * planet_radius * satellite_l_shell * cos(mlat_rad))
    end if


    !------------------------
    ! 1st adiabatic invariant
    !------------------------

    call make_adiabatic_invariant(particle_mass, magnetic_flux_density, injection_grid_number, adiabatic_invariant)

    !do count_i = 1, adiabatic_invariant_grid_number
!
    !    print *, exp(- magnetic_flux_density(count_i)*adiabatic_invariant(:, count_i)/boundary_temperature_perp)
    !    
    !end do



    !----------
    ! potential
    !----------

    call make_potential_energy(mlat_rad, length2planet, length2satellite, charge_number, particle_mass, electrostatic_potential, &
        & potential_energy)

    call make_potential_plus_Bmu(potential_energy, adiabatic_invariant, magnetic_flux_density, potential_plus_Bmu)


    !--------------
    ! accessibility
    !--------------

    call make_amin(potential_plus_Bmu, injection_grid_number, amin)
    call make_alim(particle_mass, amin, alim)
    call make_amax(potential_plus_Bmu, injection_grid_number, particle_mass, amin, amax)


    !------------
    ! calculation
    !------------

    call make_particle_flux_density(boundary_number_density, magnetic_flux_density, injection_grid_number, &
        & boundary_temperature_perp, boundary_temperature_para, particle_mass, adiabatic_invariant, potential_plus_Bmu, alim, &
        & amax, particle_flux_density)

    !次は流速


    do count_i = 1, real_grid_number

        print *, particle_flux_density(:, count_i), particle_flux_density(:, count_i)/number_density(:, count_i)
        !, (alim(:, count_i, 1)-amax(:, count_i, 1))/sqrt(boundary_temperature_para)
        
    end do



    
end program main