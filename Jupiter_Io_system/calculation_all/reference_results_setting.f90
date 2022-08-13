module reference_results_setting

    implicit none
    
    !-----------------------
    ! reference results file
    !-----------------------

    character(len=59), parameter :: reference_file = '../results/result_number_density_027_050_BC_001_min_088.csv'
    integer, parameter :: boundary_series_number = 10


    !------------
    ! result file
    !------------

    character(len=25), parameter :: result_file_front = '../results_all/result_all'
    character(len=27), parameter :: result_file_back = reference_file(33:59)
    character(len=52), parameter :: result_file = result_file_front // result_file_back

    ! bsn -> boundary_series_number
    ! coordinate_FA(1), length2planet(2), mlat_rad(3), mlat_degree(4), magnetic_flux_density(5), initial_electrostatic_potential(6),
    ! electrostatic_potential(7), number_density(8:bsn+7), charge_density(bsn+8), charge_density_Poisson(bsn+9),
    ! convergence_number(bsn+10), particle_flux_density(bsn+11:2*bsn+10), parallel_mean_velocity(2*bsn+11:3*bsn+10), 
    ! pressure_perp(3*bsn+11:4*bsn+10), pressure_para(4*bsn+11:5*bsn+10), pressure_dynamic(5*bsn+11:6*bsn+10), 
    ! temperature_perp(6*bsn+11:7*bsn+10), temperature_para(7*bsn+11:8*bsn+10), Alfven_speed(8*bsn+11), 
    ! Alfven_speed_per_lightspeed(8*bsn+12), ion_inertial_length(8*bsn+13), electron_inertial_length(8*bsn+14),
    ! ion_Larmor_radius(8*bsn+15), ion_acoustic_radius(8*bsn+16), electron_Larmor_radius(8*bsn+17), current_density(8*bsn+18)


end module reference_results_setting