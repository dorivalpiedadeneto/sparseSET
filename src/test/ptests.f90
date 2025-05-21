! This file contains general performance tests used to define which options
! are best for being used in the design of the library (even though the 
! performance depends a lot on the machine characteristics, some testing
! will hopefully help to understand some performance issues)

program ptest
    implicit none
    integer, parameter::spip = 4, spdp = 8
    integer(spip),dimension(:), allocatable::indexes

    call test_indexes_creation()

    contains

    subroutine test_indexes_creation()
        implicit none
        integer(spip):: i, n, nt
        real(spdp):: ti, tf
        !Allocating the global variable indexes
        allocate(indexes(1000))
        do i=1,1000
            indexes(i) = i
        enddo
        nt = 1000
        do n = 100,1000,100
            ! Test creating index array n times
            write(*,'(a,i4)')"For an index array of size ",n
            call cpu_time(ti)
            do i = 1, nt
                call create_index_array(n)
            enddo
            call cpu_time(tf)
            write(*,'(a,i4,a,f16.8,a)') " - Time to create index array (",&
            nt," times):", tf-ti, " seconds"
            call cpu_time(ti)
            do i = 1, nt
                call copy_index_array(n, indexes)
            enddo
            call cpu_time(tf)
            write(*,'(a,i4,a,f16.8,a)') " - Time to copy index array (",&
            nt," times):", tf-ti, " seconds"
        enddo

    end subroutine test_indexes_creation

    subroutine create_index_array(n)
        implicit none
        integer(spip), intent(in):: n
        integer(spip), dimension(:), allocatable:: ind_arr
        integer(spip):: i
        allocate(ind_arr(n))
        do i = 1, n
            ind_arr(i)  = i
        enddo
    end subroutine create_index_array

    subroutine copy_index_array(n, ind_src)
        implicit none
        integer(spip), intent(in):: n
        integer(spip), dimension(:), intent(in):: ind_src
        integer(spip), dimension(:), allocatable:: ind_arr
        allocate(ind_arr(n))
        ind_arr(1:n) = ind_src(1:n)
    end subroutine

end program ptest
