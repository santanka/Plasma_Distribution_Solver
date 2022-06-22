module reference_results_setting

    implicit none
    
    !-----------------------
    ! reference results file
    !-----------------------

    character(len=59), parameter :: reference_file = '../results/result_number_density_020_020_BC_001_min_075.csv'
    integer, parameter :: boundary_series_number = 10


    !------------
    ! result file
    !------------

    character(len=25), parameter :: result_file_front = '../results_all/result_all'
    character(len=27), parameter :: result_file_back = reference_file(33:59)
    character(len=52), parameter :: result_file = result_file_front // result_file_back


end module reference_results_setting