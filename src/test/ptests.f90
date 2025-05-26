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

    ! Corrected NR version of quicksort
    function sorted_indexes(indexes) result(sorted)
        implicit none
        integer(spip), dimension(:), intent(in):: indexes
        integer(spip), dimension(size(indexes)):: sorted
        ! Variables (used for performing quicksort *1)
        integer(spip), dimension(size(indexes))::inds
        integer(spip), parameter:: NN = 15, NSTACK = 50
        integer(spip):: ia, sa ! pivots
        integer(spip):: temp ! for swaping
        integer(spip):: n, k, i, j, jstack, l, r
        integer(spip), dimension(NSTACK):: istack
        ! Code
        n = size(indexes)
        inds = indexes ! Copy indexes to inds to preserve indexes order
        ! Creating array to return the position of indexes ordered
        do i = 1, n
            sorted(i) = i
        enddo
        jstack = 0
        l = 1
        r = n
        do ! outer loop
            if (r-l.lt.NN) then ! Insertion sort when subarray is small enough
                do j = l+1, r
                    ia = inds(j)
                    sa = sorted(j)
                    do i = j-1, l, -1
                        if (inds(i).le.ia) exit
                        inds(i+1) = inds(i)
                        sorted(i+1) = sorted(i)
                    enddo
                    inds(i+1) = ia
                    sorted(i+1) = sa
                enddo
                if (jstack.eq.0) return
                r = istack(jstack)      ! Pop stack and begin a new round of
                l = istack(jstack-1)    ! partitioning.
                jstack = jstack - 2
            else ! Choose median of left, center and right elements as
                 ! partition element ia. Also rearrange so that a(1) <= a(l+1)
                 ! <= a(r)
                k = (l+r)/2
                ! swap inds(k) with inds(l+1)
                temp = inds(k)
                inds(k) = inds(l+1)
                inds(l+1) = temp
                temp = sorted(k)
                sorted(k) = sorted(l+1)
                sorted(l+1) = temp
                ! swap inds(l) with inds(r) if inds(l) > inds(r)
                if (inds(l).gt.inds(r)) then
                    temp = inds(l)
                    inds(l) = inds(r)
                    inds(r) = temp
                    temp = sorted(l)
                    sorted(l) = sorted(r)
                    sorted(r) = temp
                endif
                ! swap inds(l+1) with inds(r) if inds(l+1) > inds(r)
                if (inds(l+1).gt.inds(r)) then
                    temp = inds(l+1)
                    inds(l+1) = inds(r)
                    inds(r) = temp
                    temp = sorted(l+1)
                    sorted(l+1) = sorted(r)
                    sorted(r) = temp
                endif
                ! swap inds(l) with inds(l+1) if inds(l) > inds(l+1)
                if (inds(l).gt.inds(l+1)) then
                    temp = inds(l)
                    inds(l) = inds(l+1)
                    inds(l+1) = temp
                    temp = sorted(l)
                    sorted(l) = sorted(l+1)
                    sorted(l+1) = temp
                endif
                i = l + 1 ! Initialize pointers for partitioning
                j = r
                ia = inds(l+1) ! Partitioning element
                sa = sorted(l+1)
                do ! inner loop
                    do  ! Scan up to find element >= ia
                        i = i + 1
                        if (inds(i).ge.ia) exit
                    enddo
                    do ! Scan dou to find element <= ia
                        j = j - 1
                        if (inds(j).le.ia) exit
                    enddo
                    if (j.lt.i) exit ! Pointers crossed. Exit with partition
                                     ! complete.
                    ! swap ind(i) with ind(j) (exchange elements)
                    temp = inds(i)
                    inds(i) = inds(j)
                    inds(j) = temp
                    temp = sorted(i)
                    sorted(i) = sorted(j)
                    sorted(j) = temp
                enddo ! end of inner loop
                inds(l+1) = inds(j)
                sorted(l+1) = sorted(j)
                inds(j) = ia
                sorted(j) = sa
                jstack = jstack + 2
                ! Push pointers to large subarray on stack; process smaller
                ! subarray immediately
                if (jstack.gt. NSTACK) then   !NSTACK is too small
                    sorted = -1 ! return with error
                    return
                endif
                if ((r-i+1).ge.(j-l)) then
                    istack(jstack) = r
                    istack(jstack-1) = i
                    r = j - 1
                else
                    istack(jstack) = j - 1
                    istack(jstack-1) = l
                    l = i
                endif
            endif
        enddo ! end of outer loop
    end function sorted_indexes
        ! *1: Based on the implementation presented in 
        ! Numerical recipes in Fortran 90: The art of
        ! parallel scientific computing (ISBN 0-521-57439-0)
        ! (page 1169-1170)



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

        write(*,*)"Recursive version"
        write(*,*)"Testing quicksort: worked? (T/F?) ->",&
        (all(out_.eq.exp_))
        write(*,*)"Time to sort",n," terms:",(tf-ti)," (s)"
        write(*,*)
        write(*,*)"Numerical Recipes version"
        
        out_ = 0
    
        call cpu_time(ti)
        out_ = sorted_indexes(in_)
        call cpu_time(tf)
        write(*,*)"Testing quicksort: worked? (T/F?) ->",&
        (all(out_.eq.exp_))
        write(*,*)"Time to sort",n," terms:",(tf-ti)," (s)"
 


    end subroutine test_quicksort

end program ptest
