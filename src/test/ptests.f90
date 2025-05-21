! This file contains general performance tests used to define which options
! are best for being used in the design of the library

program ptest
    implicit none
    integer, parameter::spip = 4
    integer(spip),dimension(:), allocatable::indexes

    call test_indexes_creation()

    contains

    subroutine test_indexes_creation()
        implicit none
        write(*,*)"Testing time to create indexes (for sorting)"

    end subroutine test_indexes_creation



end program ptest
