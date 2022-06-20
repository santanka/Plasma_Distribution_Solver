program main

    implicit none

    !固定値
    double precision, parameter :: max_B_ratio = 91.1d0
    double precision, parameter :: max_velocity_ratio = 500
    integer, parameter :: B_grid_number = 1000

    !speed grid
    integer, parameter :: velocity_space_grid = 1000

    !配列
    double precision, dimension(B_grid_number) :: B_ratio, mlat_degree
    double precision, dimension(B_grid_number) :: number_density_svc, number_density_pds, number_density_mpds
    double precision, dimension(B_grid_number) :: particle_flux_svc, particle_flux_pds, particle_flux_mpds
    double precision, dimension(B_grid_number) :: pressure_perp_svc, pressure_perp_pds, pressure_perp_mpds
    double precision, dimension(B_grid_number) :: pressure_para_svc, pressure_para_pds, pressure_para_mpds
    !double precision, dimension(B_grid_number) :: magnetic_pressure

    !保存ファイル
    character(len=128) :: file_result = "Jupiter_Io_calculation_result_SVC_PDS_mPDS.csv"
    integer :: file_i
    !file contents
    !B_ratio(1), mlat_degree(2), number_density_svc(3), number_density_pds(4), number_density_mpds(5), particle_flux_svc(6), 
    !particle_flux_pds(7), particle_flux_mpds(8), pressure_perp_svc(9), pressure_perp_pds(10), pressure_perp_mpds(11),
    !pressure_para_svc(12), pressure_para_pds(13), pressure_para_mpds(14)
    90 format(1PE25.15E3, 13(',', 1PE25.15E3))

    
    !calculation - field
    call make_B_ratio(max_B_ratio, B_grid_number, B_ratio)

    call make_MLAT(B_grid_number, B_ratio, mlat_degree)

    !calculation - number density
    call make_number_density_SVC(B_grid_number, B_ratio, number_density_svc)
    call make_number_density_PDS(B_grid_number, B_ratio, number_density_pds)
    call make_number_density_mPDS(B_grid_number, B_ratio, number_density_mpds)

    !calculation - particle flux
    call make_particle_flux_SVC(B_grid_number, B_ratio, particle_flux_svc)
    call make_particle_flux_PDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, particle_flux_pds)
    call make_particle_flux_mPDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, particle_flux_mpds)

    !calculation - pressure perp
    call make_pressure_perp_SVC(B_grid_number, B_ratio, pressure_perp_svc)
    call make_pressure_perp_PDS(B_grid_number, B_ratio, pressure_perp_pds)
    call make_pressure_perp_mPDS(B_grid_number, B_ratio, pressure_perp_mpds)

    !calculation - pressure para
    call make_pressure_para_SVC(B_grid_number, B_ratio, pressure_para_svc)
    call make_pressure_para_PDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, pressure_para_pds)
    call make_pressure_para_mPDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, pressure_para_mpds)

    !calculation - magnetic pressure
    !call make_magnetic_pressure(B_grid_number, B_ratio, magnetic_pressure)


    !file output
    open(50, file=file_result)
    do file_i = 1, B_grid_number

        write(50, 90) B_ratio(file_i), mlat_degree(file_i), number_density_svc(file_i), number_density_pds(file_i), &
            & number_density_mpds(file_i), particle_flux_svc(file_i), particle_flux_pds(file_i), particle_flux_mpds(file_i), &
            & pressure_perp_svc(file_i), pressure_perp_pds(file_i), pressure_perp_mpds(file_i), pressure_para_svc(file_i), &
            & pressure_para_pds(file_i), pressure_para_mpds(file_i)

    end do  !file_i


    
end program main

subroutine make_B_ratio(max_B_ratio, B_grid_number, B_ratio)

    implicit none
    
    double precision, intent(in) :: max_B_ratio
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number),intent(out) :: B_ratio

    integer :: grid_i

    do grid_i = 1, B_grid_number

        B_ratio(grid_i) = 1d0 * (max_B_ratio / 1d0)**(dble(grid_i - 1) / dble(B_grid_number - 1))

    end do  !grid_i

end subroutine make_B_ratio

subroutine Newton_method(initial, option, answer)

    implicit none

    double precision :: beta_func, beta_func_differential    !function
    double precision, intent(in) :: initial, option
    double precision, intent(out) :: answer
    
    double precision :: old_number, new_number
    integer :: iteration_number
    double precision, parameter :: pi = 4.d0 * atan(1.d0)

    iteration_number = 0
    new_number = initial

    do
        iteration_number = iteration_number + 1
        old_number = new_number
        new_number = old_number - beta_func(option, old_number)/beta_func_differential(old_number)
        new_number = mod(new_number, pi/2d0)
        
        if (abs(old_number - new_number) < 1E-7) then
            exit
        end if
        
        if (iteration_number == 1000000) then
            print *, "Error!: solution is not found. initial =", initial
            exit
        end if
    end do
    
    answer = new_number

end subroutine Newton_method

double precision function beta_func(beta, mlat)

    implicit none
    
    double precision, intent(in) :: beta, mlat

    beta_func = beta - sqrt(1d0 + 3d0*sin(mlat)**2d0) / cos(mlat)**6d0

end function beta_func

double precision function beta_func_differential(mlat)

    implicit none
    
    double precision, intent(in) :: mlat

    beta_func_differential = - 3d0 * tan(mlat) * (5d0 * sin(mlat)**2d0 + 3d0) / cos(mlat)**6d0 / sqrt(1d0 + 3d0*sin(mlat)**2d0)

end function beta_func_differential

subroutine make_MLAT(B_grid_number, B_ratio, mlat_degree)

    implicit none

    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: mlat_degree

    integer :: grid_i
    double precision :: mlat_initial, mlat_rad

    double precision, parameter :: pi = 4.d0 * atan(1.d0)

    mlat_degree(1) = 0d0

    do grid_i = 2, B_grid_number
        mlat_initial = 1d-5
        call Newton_method(mlat_initial, B_ratio(grid_i), mlat_rad)
        mlat_degree(grid_i) = mlat_rad / pi * 180d0
    end do !grid_i
    
end subroutine make_MLAT

subroutine make_number_density_SVC(B_grid_number, B_ratio, number_density_svc)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: number_density_svc

    number_density_svc = 5d-1 * (1d0 + sqrt(1d0 - B_ratio / B_ratio(B_grid_number)))
    
end subroutine make_number_density_SVC

subroutine make_number_density_PDS(B_grid_number, B_ratio, number_density_pds)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: number_density_pds

    number_density_pds = 5d-1 * (1d0 + sqrt(1d0 - 1d0 / B_ratio(B_grid_number)) - 2d0 * sqrt(1d0 - 1d0 / B_ratio))
    
end subroutine make_number_density_PDS

subroutine make_number_density_mPDS(B_grid_number, B_ratio, number_density_mpds)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: number_density_mpds

    number_density_mpds = 5d-1 * B_ratio * (1d0 + sqrt(1d0 - 1d0 / B_ratio(B_grid_number)) - 2d0 * sqrt(1d0 - 1d0 / B_ratio))
    
end subroutine make_number_density_mPDS

subroutine make_particle_flux_SVC(B_grid_number, B_ratio, particle_flux_svc)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: particle_flux_svc

    double precision, parameter :: pi = 4.d0 * atan(1.d0)

    particle_flux_svc = 1d0 / 2d0 / sqrt(pi) * B_ratio / B_ratio(B_grid_number)
    
end subroutine make_particle_flux_SVC

subroutine make_particle_flux_PDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, particle_flux_pds)

    implicit none
    
    integer, intent(in) :: B_grid_number, velocity_space_grid
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, intent(in) :: max_velocity_ratio
    double precision, dimension(B_grid_number), intent(out) :: particle_flux_pds

    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    integer :: grid_i
    double precision :: integral

    do grid_i = 1, B_grid_number
        
        call zeta_function_for_PDS(max_velocity_ratio, velocity_space_grid, B_ratio(B_grid_number)/B_ratio(grid_i), &
            & 1d0/B_ratio(grid_i), integral)
        
        particle_flux_pds(grid_i) = 1d0 / sqrt(pi) / B_ratio(grid_i) * integral

    end do  !grid_i
    
end subroutine make_particle_flux_PDS

subroutine make_particle_flux_mPDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, particle_flux_mpds)

    implicit none
    
    integer, intent(in) :: B_grid_number, velocity_space_grid
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, intent(in) :: max_velocity_ratio
    double precision, dimension(B_grid_number), intent(out) :: particle_flux_mpds

    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    integer :: grid_i
    double precision :: integral

    do grid_i = 1, B_grid_number
        
        call zeta_function_for_PDS(max_velocity_ratio, velocity_space_grid, B_ratio(B_grid_number)/B_ratio(grid_i), &
            & 1d0/B_ratio(grid_i), integral)
        
        particle_flux_mpds(grid_i) = 1d0 / sqrt(pi) * integral

    end do  !grid_i
    
end subroutine make_particle_flux_mPDS

subroutine zeta_function_for_PDS(max_velocity_ratio, velocity_space_grid, xx, yy, integral_result)
    !$use omp_lib

    implicit none
    
    double precision, intent(in) :: max_velocity_ratio, xx, yy
    integer, intent(in) :: velocity_space_grid
    double precision, intent(out) :: integral_result

    integer :: grid_u, grid_v
    double precision, dimension(velocity_space_grid) :: matrix_u, matrix_v, result_halfway
    double precision :: integral_former, integral_latter

    matrix_u(1) = 1d-30

    !$omp parallel
    !$omp do
    do grid_u = 2, velocity_space_grid

        matrix_u(grid_u) = 1d-5 * (max_velocity_ratio / 1d-5)**(dble(grid_u - 2) / dble(velocity_space_grid - 2))

    end do  !grid_u
    !$omp end do
    !$omp end parallel

    result_halfway = 0d0

    !$omp parallel private(matrix_v, grid_v, integral_former, integral_latter)
    !$omp do
    do grid_u = 1, velocity_space_grid
        
        if ( (xx - 1d0) * matrix_u(grid_u) <= 1d-25 ) then
            matrix_v(1) = 1d-30
            do grid_v = 2, velocity_space_grid
                
                matrix_v(grid_v) = 1d-5 * (max_velocity_ratio / 1d-5)**(dble(grid_v - 2) / dble(velocity_space_grid - 2))

            end do  !grid_v

        else if ( (xx - 1d0) * matrix_u(grid_u) > 1d-25 ) then
            do grid_v = 1, velocity_space_grid
                
                matrix_v(grid_v) = sqrt((xx - 1d0) * matrix_u(grid_u)) * &
                    & (max_velocity_ratio / sqrt((xx - 1d0) * matrix_u(grid_u)))**(dble(grid_v - 1) / dble(velocity_space_grid - 1))

            end do  !grid_v
        end if

        do grid_v = 1, velocity_space_grid - 1

            integral_former = matrix_v(grid_v)**2d0 / sqrt(matrix_v(grid_v)**2d0 + (1d0 - yy) * matrix_u(grid_u)) &
                & * exp(- matrix_v(grid_v)**2d0 - matrix_u(grid_u))

            integral_latter = matrix_v(grid_v + 1)**2d0 / sqrt(matrix_v(grid_v + 1)**2d0 + (1d0 - yy) * matrix_u(grid_u)) &
                & * exp(- matrix_v(grid_v + 1)**2d0 - matrix_u(grid_u))

            result_halfway(grid_u) = result_halfway(grid_u) + (integral_former + integral_latter) / 2d0 &
                & * (matrix_v(grid_v + 1) - matrix_v(grid_v))
            
        end do  !grid_v

    end do  !grid_u
    !$omp end do
    !$omp end parallel

    integral_result = 0d0

    do grid_u = 1, velocity_space_grid - 1
        
        integral_result = integral_result + (result_halfway(grid_u) + result_halfway(grid_u + 1)) / 2d0 &
            & * (matrix_u(grid_u + 1) - matrix_u(grid_u))

    end do  !grid_u
    
end subroutine zeta_function_for_PDS

subroutine make_pressure_perp_SVC(B_grid_number, B_ratio, pressure_perp_svc)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_perp_svc

    pressure_perp_svc = 2.5d-1 * (2d0 + (2d0 + B_ratio / B_ratio(B_grid_number)) * sqrt(1d0 - B_ratio / B_ratio(B_grid_number)))
    
end subroutine make_pressure_perp_SVC

subroutine make_pressure_perp_PDS(B_grid_number, B_ratio, pressure_perp_pds)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_perp_pds

    pressure_perp_pds = 2.5d-1 * B_ratio * (2d0 + (2d0 + 1 / B_ratio(B_grid_number)) * sqrt(1d0 - 1d0 / B_ratio(B_grid_number)) &
        & - 2d0 * (2d0 + 1d0 / B_ratio) * sqrt(1d0 - 1d0 / B_ratio))
    
end subroutine make_pressure_perp_PDS

subroutine make_pressure_perp_mPDS(B_grid_number, B_ratio, pressure_perp_mpds)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_perp_mpds

    pressure_perp_mpds = 2.5d-1 * B_ratio**2d0 * (2d0 + (2d0 + 1 / B_ratio(B_grid_number)) &
    & * sqrt(1d0 - 1d0 / B_ratio(B_grid_number)) - 2d0 * (2d0 + 1d0 / B_ratio) * sqrt(1d0 - 1d0 / B_ratio))
    
end subroutine make_pressure_perp_mPDS

subroutine make_pressure_para_SVC(B_grid_number, B_ratio, pressure_para_svc)

    implicit none
    
    integer, intent(in) :: B_grid_number
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_para_svc

    double precision, parameter :: pi = 4.d0 * atan(1.d0)

    pressure_para_svc = 5d-1 * (1d0 + (1d0 - B_ratio / B_ratio(B_grid_number))**1.5d0 &
        & - 2d0 / pi * (B_ratio / B_ratio(B_grid_number))**2d0 / (1d0 + sqrt(1d0 - B_ratio / B_ratio(B_grid_number))))
    
end subroutine make_pressure_para_SVC

subroutine make_pressure_para_PDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, pressure_para_pds)

    implicit none
    
    integer, intent(in) :: B_grid_number, velocity_space_grid
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, intent(in) :: max_velocity_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_para_pds

    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    integer :: grid_i
    double precision :: integral, bb_bi, be_bi, bb_be

    do grid_i = 1, B_grid_number

        bb_bi = 1d0/B_ratio(grid_i)
        be_bi = B_ratio(B_grid_number)/B_ratio(grid_i)
        bb_be = 1d0/B_ratio(B_grid_number)
        
        
        call zeta_function_for_PDS(max_velocity_ratio, velocity_space_grid, be_bi, bb_bi, integral)
        
        pressure_para_pds(grid_i) = 5d-1 * ((1d0 + sqrt(1d0 - bb_be) - 2d0 * sqrt(1d0 - bb_bi)) &
            & - 8d0 / pi * bb_bi**2d0 / (1d0 + sqrt(1d0 - bb_be) - 2d0 * sqrt(1d0 - bb_bi)) * integral**2d0 &
            & + (1d0 - 1d0 / bb_bi) * (2d0 + sqrt(1d0 - bb_be) * (2d0 + bb_be) - 2d0 * sqrt(1d0 - bb_bi) * (2d0 + bb_bi)) &
            & - bb_be * sqrt(1d0 - bb_be) + 2d0 * bb_bi * sqrt(1d0 - bb_bi))

    end do  !grid_i
    
end subroutine make_pressure_para_PDS

subroutine make_pressure_para_mPDS(B_grid_number, B_ratio, max_velocity_ratio, velocity_space_grid, pressure_para_mpds)

    implicit none
    
    integer, intent(in) :: B_grid_number, velocity_space_grid
    double precision, dimension(B_grid_number), intent(in) :: B_ratio
    double precision, intent(in) :: max_velocity_ratio
    double precision, dimension(B_grid_number), intent(out) :: pressure_para_mpds

    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    integer :: grid_i
    double precision :: integral, bb_bi, be_bi, bb_be

    do grid_i = 1, B_grid_number
        
        bb_bi = 1d0/B_ratio(grid_i)
        be_bi = B_ratio(B_grid_number)/B_ratio(grid_i)
        bb_be = 1d0/B_ratio(B_grid_number)
        
        
        call zeta_function_for_PDS(max_velocity_ratio, velocity_space_grid, be_bi, bb_bi, integral)
        
        pressure_para_mpds(grid_i) = 5d-1 / bb_bi * ((1d0 + sqrt(1d0 - bb_be) - 2d0 * sqrt(1d0 - bb_bi)) &
        & - 8d0 / pi * bb_bi**2d0 / (1d0 + sqrt(1d0 - bb_be) - 2d0 * sqrt(1d0 - bb_bi)) * integral**2d0 &
        & + (1d0 - 1d0 / bb_bi) * (2d0 + sqrt(1d0 - bb_be) * (2d0 + bb_be) - 2d0 * sqrt(1d0 - bb_bi) * (2d0 + bb_bi)) &
        & - bb_be * sqrt(1d0 - bb_be) + 2d0 * bb_bi * sqrt(1d0 - bb_bi))

    end do  !grid_i
    
end subroutine make_pressure_para_mPDS

!subroutine make_magnetic_pressure(B_grid_number, B_ratio, magnetic_pressure)
!    
!    implicit none
!    
!    integer, intent(in) :: B_grid_number
!    double precision, dimension(B_grid_number), intent(in) :: B_ratio
!    double precision, dimension(B_grid_number), intent(out) :: magnetic_pressure
!
!    double precision, parameter :: magnetic_constant = 1.25663706212d-6
!
!    magnetic_pressure = B_ratio**2d0 / 2d0 / magnetic_constant
!    
!end subroutine make_magnetic_pressure