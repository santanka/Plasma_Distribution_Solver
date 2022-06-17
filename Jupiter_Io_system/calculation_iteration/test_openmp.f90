program test_openmp
    use omp_lib

    implicit none

    integer :: i

    !$omp parallel do private(i)
    do i = 0, 15
        print *, "i and thread num", i, omp_get_thread_num()
    end do
    
end program test_openmp