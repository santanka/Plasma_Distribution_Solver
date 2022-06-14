module boundary_and_initial_conditions
    use constant_in_the_simulation

    implicit none

    !------------------
    ! initial condition
    !------------------

    integer, parameter :: initial_fix_grid = 117
    integer, parameter :: initial_min_grid_1 = 74
    integer, parameter :: initial_min_grid_2 = real_grid_number + 1 - initial_min_grid_1

    integer, parameter :: initial_grid_ionophere_middle_1 = 5
    ! 1 ~ initial_grid_ionophere_middle_1 - 1   initial_electrostatic_potential_ionosphere
    
    integer, parameter :: initial_grid_middle_magnetosphere_1 = 30
    ! initial_grid_ionosphere_middle_1 ~ initial_grid_middle_magnetosphere_1 - 1    initial_electrostatic_potential_middle

    integer, parameter :: initial_grid_middle_magnetosphere_2 = real_grid_number + 2 - initial_grid_middle_magnetosphere_1
    ! initial_grid_middle_magnetosphere_1 ~ initial_grid_middle_magnetosphere_2 - 1    initial_electrostatic_potential_magnetosphere
  
    integer, parameter :: initial_grid_ionophere_middle_2 = real_grid_number + 2 - initial_grid_ionophere_middle_1
    ! initial_grid_middle_magnetosphere_2 ~ initial_grid_ionophere_middle_2 - 1     initial_electrostatic_potential_middle
    ! initial_grid_ionophere_middle_2 ~ real_grid_number    initial_electrostatic_potential_ionosphere

    double precision, parameter :: initial_electrostatic_potential_ionosphere = 3d4     ![V] = [kg m2 s-3 A-1]
    double precision, parameter :: initial_electrostatic_potential_middle = 2.6d4       ![V] = [kg m2 s-3 A-1]
    double precision, parameter :: initial_electrostatic_potential_magnetosphere = 0d0  ![V] = [kg m2 s-3 A-1]


    !-------------------
    ! boundary condition
    !-------------------

    character(len=128) :: boundary_file = '../conditions/boundary_conditions_1.csv'
    integer, parameter :: boundary_file_number = 1  !boundary_file
    integer, parameter :: boundary_series_number = 10
    integer, parameter :: boundary_ionosphere_1_variable_species = 2
    integer, parameter :: boundary_ionosphere_2_variable_species = 4
    integer, parameter :: boundary_magnetosphere_variable_species = 6


    !------------
    ! result file
    !------------

    character(len = 32) :: result_file_front = '../results/result_number_density' 

end module boundary_and_initial_conditions    