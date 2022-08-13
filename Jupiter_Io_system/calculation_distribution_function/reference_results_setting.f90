module reference_results_setting

    implicit none
    
    !-----------------------
    ! reference results file
    !-----------------------

    character(len=59), parameter :: reference_file = '../results/result_number_density_027_027_BC_002_min_063.csv'
    integer, parameter :: boundary_series_number = 10
    integer, parameter :: plot_grid_number = 50


    !------------
    ! result file
    !------------

    character(len=60), parameter :: result_file_front = '../result_distribution_function/result_distribution_function'
    character(len=23), parameter :: result_file_back = reference_file(33:55)

    ! mlat_degree(1), v_perp_i(2), v_para_i(3), distribution_function_i(4), v_perp_b(5), v_para_b(6), distribution_function_b(7)


end module reference_results_setting