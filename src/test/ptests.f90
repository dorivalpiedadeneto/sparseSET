! This file contains general performance tests used to define which options
! are best for being used in the design of the library (even though the 
! performance depends a lot on the machine characteristics, some testing
! will hopefully help to understand some performance issues)

program ptest
    implicit none
    integer, parameter::spip = 4, spdp = 8
    integer(spip),dimension(:), allocatable::indexes

    call test_indexes_creation()
    call test_quicksort()

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

    ! ChatGPT suggested algorithm (to help finding error in the adapted
    ! Numerical Recipes algorithm)
    recursive subroutine quicksort_indices_sub(values, idx, low, high)
        implicit none
        integer, intent(in) :: low, high
        integer, dimension(:), intent(in) :: values
        integer, dimension(:), intent(inout) :: idx
        integer :: pivotIndex

        if (low < high) then
            call partition_indices(values, idx, low, high, pivotIndex)
            call quicksort_indices_sub(values, idx, low, pivotIndex - 1)
            call quicksort_indices_sub(values, idx, pivotIndex + 1, high)
        end if
    end subroutine quicksort_indices_sub

    !------------------------------------------------------------
    subroutine partition_indices(values, idx, low, high, pivotIndex)
        implicit none
        integer, intent(in) :: low, high
        integer, dimension(:), intent(in) :: values
        integer, dimension(:), intent(inout) :: idx
        integer, intent(out) :: pivotIndex
        integer :: i, j, temp, pivotValue

        pivotValue = values(idx(high))
        i = low - 1

        do j = low, high - 1
            if (values(idx(j)) <= pivotValue) then
                i = i + 1
                temp = idx(i)
                idx(i) = idx(j)
                idx(j) = temp
            end if
        end do

        temp = idx(i + 1)
        idx(i + 1) = idx(high)
        idx(high) = temp

        pivotIndex = i + 1
    end subroutine partition_indices


    subroutine test_quicksort()
        implicit none
        integer, dimension(:), allocatable::in_, out_, exp_, ind_
        integer::i, n
        real(8)::ti, tf

        n = 10000

        allocate(in_(n),out_(n), exp_(n))
        do i=1,n
            in_(i) = n + 1 - i
            exp_(i) = n + 1 - i ! expected result
            out_(i) = i
        enddo

        call cpu_time(ti)
        call quicksort_indices_sub(in_, out_, 1, n)
        call cpu_time(tf)

        write(*,*)"Testing quicksort: worked? (T/F?) ->",&
        (all(out_.eq.exp_))
        write(*,*)"Time to sort ",n," terms:",(tf-ti)," (s)"


    end subroutine test_quicksort

end program ptest
