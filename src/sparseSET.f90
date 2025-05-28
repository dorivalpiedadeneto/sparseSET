! sparseSET - A library for assembling sparse matrices in fortran
! Authors: Dorival Piedade Neto
!          Rodrigo Ribeiro Paccola
!          Rogerio Carrazedo
!          Version: 2.0 (2024)
! License: BSD 2-Clause License (FreeBSD/Simplified)
module sparseset
    implicit none
    
    ! Types of variables (hardcoded; user must change if needed)
    integer, parameter::spip = 4 ! integer precision for sparseSET library
    integer, parameter::spdp = 8 ! real precision in sparseSET library
    ! The defaults used for the library are: spip for integers, spdp for reals
    ! (the sp in the variables name is used to avoid naming conflicts - this  
    !  library will be used by others, so it is better avoiding name conflicts)
    !  ip 'stands' for integer precision, dp 'stands' for double precision of
    !  floating point numbers

    ! Remark: we are using types and 'structured' programming instead of OOP
    ! (so that it can be used in older compilers, that don't support OOP)

    ! A sparse line is used to store rows or columns of a sparse matrix
    type sparse_line
        integer(spip)::lsize = 0
        integer(spip)::lcount = 0
        integer(spip), dimension(:), allocatable::lindex
        real(spdp), dimension(:), allocatable::lvalue
        integer(spip)::rpstage = 0
        logical::assembled=.false.
        integer(spip)::length
    end type sparse_line

    ! About the sparse_line type
    ! Sparse lines are data structures to hold rows or column data while a 
    ! sparse matrix is being assembled. Terms are pushed to the end, and when
    ! all  contributions to the sparse matrix are done, them the terms are 
    ! ordered and terms of the same index are summed up.
    ! The line data (index and value) are stored in one dimensional arrays.
    ! The initial space allocated to store these values depend on the system
    ! of equation characteristics. In some cases, terms of same index are
    ! pushed several times, resulting in excessive space usage and reallocation.
    ! In these cases, sometimes a better approach is to sum all terms of same
    ! index to find more space in the line without reallocating.
    ! In short, to deal with all possible cenarios, the library defines a type
    ! resize_policy to be followed for each sparse_matrix entity.
    !
    ! Variables in the line type
    !
    ! lsize - the current size of lindex and lvalue arrays
    ! count - the number of current terms pushed in the row
    ! lindex - array of index, to hold the term index
    ! lvalue - array of values, to the the term value
    ! rpstate - an integer representing the current state
    !           of the line considering the sparse matrix
    !           resize policy (see more below)
    ! assembled - a logical representing the current state
    !             of the line; every time it is 'assembled',
    !             it is set to .true.; every time terms are
    !             pushed, it is set to .false.
    ! length - the 'real' length of the row/column

    type sparse_matrix
        character(3)::mtype
        character(5)::storage
        integer(spip)::isize
        integer(spip),dimension(:), allocatable::resize_policy
        integer(spip)::nlines
        integer(spip)::nrows
        integer(spip)::ncols
        type(sparse_line), dimension(:), allocatable::line
    end type sparse_matrix

    ! About sparse_matrix type
    ! The sparse_matrix type is used to hold data in lines while the sparse
    ! matrix is being assembled.
    !
    ! Variables in the sparse_matrix type
    ! mtype - indicates if the lines represent rows ('row') or columns
    ! ('col')
    ! storage - indicates which parts of the matrix is stored: 'full', 'upper'
    ! or 'lower'
    ! resize_policy - an array of integers indicanting the policy to resize
    ! the matrix; the integers in such a matrix have the following meaning
    ! -> positive number - the size to reallocate the line, considering its
    ! initial size. For instance, 2 indicates to resize to 2*isize (the initial
    ! line size); 3 indicates to resize to 3*isize, and so on.
    ! -> zero - instead of reallocating, try to sum up equal terms without 
    ! reallocating to find more space in the line
    ! -> negative number - the size to reallocate the line, considering its
    ! current size (ignoring the sign). For istance, -2 means: resize to
    ! 2*current_line_size; -3 means: 3 * current_line_size
    ! (It would be possible to use a more fancy data structure to represent
    ! the policy; at least for now, lets use this simple and naive approach)
    ! nlines - number of lines representing the sparse_matrix
    ! nrows - number of rows in the sparse matrix
    ! ncols - number of columns in the sparse matrix
    ! line - an array of lines

    ! A note on the resize policy:
    ! Initially resize_policy was a fixed array with 16 positions:
    !
    !   integer(spip),dimension(16)::resize_policy = (\ &
    !   2, 3, 4, 8, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\)
    !
    ! The current approach is to set is as allocatable; if the
    ! user doesn't define it when allocating the matrix (optional
    ! parameter, adopt it as:
    ! (\2, 3, 4, 8, 16, 0\)
    ! (in this case for rp >= 6, always try to assemble; if the available space
    ! is not sufficient, quit the program)

    type CSC
        integer(spip)::msize
        integer(spip), dimension(:), allocatable::nz
        integer(spip), dimension(:), allocatable::row
        integer(spdp), dimension(:), allocatable::mvalue
    end type CSC

    type CSR
        integer(spip)::msize
        integer(spip), dimension(:), allocatable::nz
        integer(spip), dimension(:), allocatable::col
        integer(spdp), dimension(:), allocatable::mvalue
    end type CSR

    type triplet
        integer(spip)::msize
        integer(spip), dimension(:), allocatable::nnz
        integer(spip), dimension(:), allocatable::row
        integer(spip), dimension(:), allocatable::col
        integer(spdp), dimension(:), allocatable::mvalue
    end type triplet

    contains

    subroutine allocate_sparse_line(sp_line, line_size, length, stat)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip), intent(in):: line_size
        integer(spip), intent(in):: length
        integer(spip), intent(inout),optional:: stat
        ! Testing if informed line_size and lenght are valid
        if (present(stat)) then
            if (.not.((line_size.gt.0).and.(length.gt.0))) then
                stat = 1
                return
            else
                stat = 0
            endif
        endif
        call deallocate_sparse_line(sp_line)
        sp_line%lsize = line_size
        sp_line%length = length
        allocate(sp_line%lindex(line_size))
        allocate(sp_line%lvalue(line_size))
    end subroutine allocate_sparse_line

    ! About the allocate_sparse_line subroutine:
    ! line_size is the size to be allocated in the sparse_line variables
    ! lindex and lvalue. Length is the actual size of the line represented
    ! by the sparse_line (row or column of a matriz, or even a vector).
    ! While pushing terms, several duplicate terms occur in these arrays.
    ! Therefore, line_size can be greater than length (to hold duplicates).
    ! Both line_size and length must be greater than zero. The optional
    ! parameter 'stat' returns 1 if one tries to use invalid values for them,
    ! and returns 0 when they are suitable and it is possible to perform the
    ! allocation.

    subroutine deallocate_sparse_line(sp_line)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        if (allocated(sp_line%lindex)) deallocate(sp_line%lindex)
        if (allocated(sp_line%lvalue)) deallocate(sp_line%lvalue)
        sp_line%lsize = 0
        sp_line%lcount = 0
        sp_line%rpstage = 0
        sp_line%assembled = .false.
        sp_line%length = 0
    end subroutine deallocate_sparse_line

    function available_space(sp_line) result(n)
        implicit none
        type(sparse_line),intent(in)::sp_line
        integer(spip)::n
        n = sp_line%lsize - sp_line%lcount
        return
    end function available_space

    subroutine push_term_to_line(sp_line, index, value, stat)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip), intent(in):: index
        real(spdp), intent(in):: value
        integer(spip), intent(out), optional::stat
        if (present(stat)) then
            if ((sp_line%lsize - sp_line%lcount).gt.0) then
                stat = 0
            else
                stat = 1
                return
            endif
        endif
        ! if not present stat, test is not performed (the program will
        ! crash if user tries to push when there is no space)
        sp_line%lcount = sp_line%lcount + 1
        sp_line%lindex(sp_line%lcount) = index
        sp_line%lvalue(sp_line%lcount) = value
        sp_line%assembled = .false.
    end subroutine push_term_to_line

    subroutine push_terms_to_line(sp_line, indexes, values, stat)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip), dimension(:), intent(in):: indexes
        real(spdp), dimension(:), intent(in):: values
        integer(spip), intent(out), optional::stat
        integer(spip)::nt,av,sb
        av = sp_line%lsize - sp_line%lcount
        ! av - available space to push new terms
        nt = size(indexes)
        ! nt - number of terms to push
        if (present(stat)) then
            if (size(indexes).ne.size(values)) then
                stat = 1
                return
            else
                if (nt.gt.av) then
                    stat = 1
                    return
                endif
            endif
            if (nt.gt.av) then
                stat = 1
                return
            else
                stat = 0
            endif
        endif
        ! if not present stat, test is not performed (the program will
        ! crash if user tries to push when there is no space)
        sb = sp_line%lcount + 1
        ! First position where terms will be pushed (slice begin)
        sp_line%lindex(sb:sb+nt) = indexes
        sp_line%lvalue(sb:sb+nt) = values
        sp_line%lcount = sp_line%lcount + nt !lcount + number of terms pushed
        sp_line%assembled = .false.
    end subroutine push_terms_to_line
    
    subroutine clear_sparse_line(sp_line)
        implicit none
        type(sparse_line), intent(inout)::sp_line
        ! Clear = reset all data in line keeping the allocated space
        sp_line%lcount = 0
        sp_line%assembled = .false.
    end subroutine clear_sparse_line

    subroutine copy_sparse_line_terms(origin, destination, stat)
        implicit none
        type(sparse_line), intent(in):: origin
        type(sparse_line), intent(inout):: destination
        integer(spip), intent(out), optional:: stat
        integer(spip)::lcount
        if (present(stat)) then
            if (destination%lsize.lt.origin%lcount) then
                stat = 1
                return
            else
                stat = 0
            endif
        endif
        ! If stat is not present and space is insufficient,
        ! the program will crash)
        lcount = origin%lcount
        destination%lcount = origin%lcount
        destination%lindex(1:lcount) = origin%lindex(1:lcount)
        destination%lvalue(1:lcount) = origin%lvalue(1:lcount)
        destination%assembled = origin%assembled
    end subroutine copy_sparse_line_terms

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

    subroutine assemble_sparse_line(spline)
        implicit none
        type(sparse_line), intent(inout):: spline
        integer(spip), dimension(spline%lcount):: indexes,  ind
        real(spdp), dimension(spline%lcount):: values
        integer(spip):: i, pos
        ind = sorted_indexes(spline%lindex(1:spline%lcount))
        ! Copying indexes in ascending order (and respective values)
        do i = 1, spline%lcount
            indexes(i) = spline%lindex(ind(i))
            values(i) = spline%lvalue(ind(i))
        enddo
        ! Summing equal terms
        pos = 1
        do i = 2, spline%lcount
            if (indexes(i).eq.indexes(pos)) then
                values(pos) = values(pos) + values(i)
            else
                pos = pos + 1
                indexes(pos) = indexes(i)
                values(pos) = values(i)
            endif
        enddo
        spline%lindex(1:pos) = indexes(1:pos)
        spline%lvalue(1:pos) = values(1:pos)
        spline%lcount = pos
        spline%assembled = .true.
    end subroutine assemble_sparse_line

end module sparseset

