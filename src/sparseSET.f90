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

    ! Default values
    integer(spip):: default_isize = 100
    character(3):: default_mtype = 'row'
    character(5):: default_storage = 'full'
    integer(spip), dimension(6):: default_resize_policy = (/2,3,4,8,16,0/)
    ! End of default values

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

    ! This version used sorted_indexes; the used approach is much slower
    ! than the new version, that uses sort_indexes_and_values (keeping
    ! the old version here, for now)
    subroutine assemble_sparse_line_old_version(sp_line)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip), dimension(sp_line%lcount):: indexes,  ind
        real(spdp), dimension(sp_line%lcount):: values
        integer(spip):: i, pos
        ind = sorted_indexes(sp_line%lindex(1:sp_line%lcount))
        ! Copying indexes in ascending order (and respective values)
        do i = 1, sp_line%lcount
            indexes(i) = sp_line%lindex(ind(i))
            values(i) = sp_line%lvalue(ind(i))
        enddo
        ! Summing equal terms
        pos = 1
        do i = 2, sp_line%lcount
            if (indexes(i).eq.indexes(pos)) then
                values(pos) = values(pos) + values(i)
            else
                pos = pos + 1
                indexes(pos) = indexes(i)
                values(pos) = values(i)
            endif
        enddo
        sp_line%lindex(1:pos) = indexes(1:pos)
        sp_line%lvalue(1:pos) = values(1:pos)
        sp_line%lcount = pos
        sp_line%assembled = .true.
    end subroutine assemble_sparse_line_old_version

    subroutine sort_indexes_and_values(indexes, values, length, stat)
        implicit none
        integer(spip), dimension(:), intent(inout):: indexes
        real(spdp), dimension(:), intent(inout):: values
        integer(spip), intent(in):: length
        integer(spip), intent(out), optional:: stat
        ! Variables (used for performing quicksort *1)
        integer(spip), parameter:: NN = 15, NSTACK = 50
        integer(spip):: ia, itemp ! pivot and swap var.
        real(spdp)::ra, rtemp ! pivot and swap var.
        integer(spip):: n, k, i, j, jstack, l, r
        integer(spip), dimension(NSTACK):: istack
        if (present(stat)) stat = 0
        ! Code
        n = size(indexes)
        jstack = 0
        l = 1
        r = n
        do ! outer loop
            if (r-l.lt.NN) then ! Insertion sort when subarray is small enough
                do j = l+1, r
                    ia = indexes(j)
                    ra = values(j)
                    do i = j-1, l, -1
                        if (indexes(i).le.ia) exit
                        indexes(i+1) = indexes(i)
                        values(i+1) = values(i)
                    enddo
                    indexes(i+1) = ia
                    values(i+1) = ra
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
                itemp = indexes(k)
                indexes(k) = indexes(l+1)
                indexes(l+1) = itemp
                rtemp = values(k)
                values(k) = values(l+1)
                values(l+1) = rtemp
                ! swap inds(l) with inds(r) if inds(l) > inds(r)
                if (indexes(l).gt.indexes(r)) then
                    itemp = indexes(l)
                    indexes(l) = indexes(r)
                    indexes(r) = itemp
                    rtemp = values(l)
                    values(l) = values(r)
                    values(r) = rtemp
                endif
                ! swap inds(l+1) with inds(r) if inds(l+1) > inds(r)
                if (indexes(l+1).gt.indexes(r)) then
                    itemp = indexes(l+1)
                    indexes(l+1) = indexes(r)
                    indexes(r) = itemp
                    rtemp = values(l+1)
                    values(l+1) = values(r)
                    values(r) = rtemp
                endif
                ! swap inds(l) with inds(l+1) if inds(l) > inds(l+1)
                if (indexes(l).gt.indexes(l+1)) then
                    itemp = indexes(l)
                    indexes(l) = indexes(l+1)
                    indexes(l+1) = itemp
                    rtemp = values(l)
                    values(l) = values(l+1)
                    values(l+1) = rtemp
                endif
                i = l + 1 ! Initialize pointers for partitioning
                j = r
                ia = indexes(l+1) ! Partitioning element
                ra = values(l+1)
                do ! inner loop
                    do  ! Scan up to find element >= ia
                        i = i + 1
                        if (indexes(i).ge.ia) exit
                    enddo
                    do ! Scan dou to find element <= ia
                        j = j - 1
                        if (indexes(j).le.ia) exit
                    enddo
                    if (j.lt.i) exit ! Pointers crossed. Exit with partition
                                     ! complete.
                    ! swap ind(i) with ind(j) (exchange elements)
                    itemp = indexes(i)
                    indexes(i) = indexes(j)
                    indexes(j) = itemp
                    rtemp = values(i)
                    values(i) = values(j)
                    values(j) = rtemp
                enddo ! end of inner loop
                indexes(l+1) = indexes(j)
                values(l+1) = values(j)
                indexes(j) = ia
                values(j) = ra
                jstack = jstack + 2
                ! Push pointers to large subarray on stack; process smaller
                ! subarray immediately
                if (jstack.gt. NSTACK) then   !NSTACK is too small
                    if (present(stat)) stat = -1 ! return with error
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
    end subroutine sort_indexes_and_values

    subroutine assemble_sparse_line(sp_line)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip)::i, pos
        call sort_indexes_and_values(sp_line%lindex(1:sp_line%lcount),&
             sp_line%lvalue(1:sp_line%lcount), sp_line%lcount)
        pos = 1
        do i = 2, sp_line%lcount
            if (sp_line%lindex(i).eq.sp_line%lindex(pos)) then
                sp_line%lvalue(pos) = sp_line%lvalue(pos) + sp_line%lvalue(i)
            else
                pos = pos + 1
                sp_line%lindex(pos) = sp_line%lindex(i)
                sp_line%lvalue(pos) = sp_line%lvalue(i)
            endif
        enddo
        sp_line%lcount = pos
        sp_line%assembled = .true.
    end subroutine assemble_sparse_line

    ! If index is found, returns its position in line;
    ! if not, returns -1.
    function search_line(sp_line, ind) result(pos)
        implicit none
        type(sparse_line), intent(inout):: sp_line
        integer(spip), intent(in):: ind
        integer(spip):: pos
        integer(spip):: ib, ie, im !b -> begin; e -> end; m -> middle

        if (sp_line%lcount.ne.0) then
            ! Sparse line must be assembled (if it is not, assemble it)
            if (.not.sp_line%assembled) call assemble_sparse_line(sp_line)
            ib = 1
            ie = sp_line%lcount
            if ((sp_line%lindex(ib).gt.ind).or.(sp_line%lindex(ie).lt.ind)) then
                pos = -1
                return
            endif
            if (sp_line%lindex(ib).eq.ind) then
                pos = ib
                return
            endif
            if (sp_line%lindex(ie).eq.ind) then
                pos = ie
                return
            endif

            do while ((ie - ib).gt.1)
                if (sp_line%lindex(ib).eq.ind) then
                    pos = ib
                    return
                endif
                if (sp_line%lindex(ie).eq.ind) then
                    pos = ie
                    return
                endif
                im = ib + (ie - ib) / 2
                if (sp_line%lindex(im).eq.ind) then
                    pos = im
                    return
                endif
                if (ind.lt.sp_line%lindex(im)) then
                    ie = im
                else
                    ib = im
                endif
            enddo
            ! If it gets here, index was not found
            pos = -1
            return
        endif
        ! Line has zero terms
        pos = -1
        return
    end function search_line

    subroutine sparse_line_to_array(sp_line, array)
        implicit none
        type(sparse_line), intent(in):: sp_line
        real(spdp), dimension(:), allocatable, intent(inout):: array
        integer(spip):: i
        if (size(array).ne.sp_line%length) then
           if (allocated(array)) deallocate(array)
           allocate(array(sp_line%length))
        endif
        array = 0.0_spdp
        do i=1, sp_line%lcount
            array(sp_line%lindex(i)) = array(sp_line%lindex(i)) + & 
            sp_line%lvalue(i)
        enddo
    end subroutine sparse_line_to_array

    subroutine array_to_sparse_line(array, sp_line, tolerance)
        implicit none
        real(spdp), dimension(:), intent(in):: array
        type(sparse_line), intent(inout):: sp_line
        real(spdp), optional, intent(in):: tolerance
        real(spdp):: tol
        integer(spip):: i, length, lsize, lsze, cnt
        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 1.0e-8_spdp
        endif
        length = size(array)
        lsize = count((dabs(array).gt.tol), dim = 1, kind = spip)
        if (sp_line%lsize.lt.lsize) then
            call deallocate_sparse_line(sp_line)
            ! An 'heuristic' for defininf the new lsize value
            lsze = (lsize / 50 + 1) * 50
            call allocate_sparse_line(sp_line, lsze, length)
        endif
        cnt = 0
        do i = 1, length
            if (dabs(array(i)).gt.tol) then
                cnt = cnt + 1
                sp_line%lindex(cnt) = i
                sp_line%lvalue(cnt) = array(i)
            endif
        enddo
        sp_line%lcount = cnt
        sp_line%length = length
        sp_line%assembled = .true.
    end subroutine array_to_sparse_line

    ! Sparse Matrix Subroutines

    subroutine allocate_sparse_matrix(sp_matrix, nrows, ncols, isize, &
                mtype, storage, resize_policy, stat)
        implicit none
        type(sparse_matrix), intent(inout)::sp_matrix
        integer(spip), intent(in)::nrows, ncols
        integer(spip), optional, intent(in)::isize
        character(3), optional, intent(in)::mtype
        character(5), optional, intent(in)::storage
        integer(spip), optional, dimension(:), intent(in):: resize_policy
        integer(spip), optional, intent(out):: stat
        ! Local variables
        integer(spip)::i, nlines, length
        ! Testing values if stat is present (else ignore and code may crash)
        if (present(stat)) then
            stat = 0
            if ((nrows.le.0).or.(ncols.le.0)) then
                stat = 1
                return
            else
                sp_matrix%nrows = nrows
                sp_matrix%ncols = ncols
            endif
            if (present(isize)) then
                if (isize.le.0) then
                    stat = 1
                    return
                else
                    sp_matrix%isize = isize
                endif
            else
                !Default value for isize
                sp_matrix%isize = default_isize
            endif
            if (present(mtype)) then
                if (.not.((mtype.eq.'col').or.(mtype.eq.'row'))) then
                    stat = 1
                    return
                else
                    sp_matrix%mtype = mtype
                endif
            else
                sp_matrix%mtype = default_mtype
            endif
            if (present(storage)) then
                if (.not.(storage.eq.'upper').or.(storage.eq.'lower')&
                    .or.(storage.eq.'full')) then
                    stat = 1
                    return
                else
                    sp_matrix%storage = default_storage
                endif
            else
                sp_matrix%storage = default_storage
            endif
        endif
        if (allocated(sp_matrix%line)) then
            do i = 1, size(sp_matrix%line)
                call deallocate_sparse_line(sp_matrix%line(i))
            enddo
            deallocate(sp_matrix%line)
        endif
        if (allocated(sp_matrix%resize_policy)) then
            deallocate(sp_matrix%resize_policy)
        endif
        ! Allocating lines and resize_policy
        if (present(resize_policy)) then
            allocate(sp_matrix%resize_policy(size(resize_policy)))
            sp_matrix%resize_policy = resize_policy
        else
            allocate(sp_matrix%resize_policy(size(default_resize_policy)))
            sp_matrix%resize_policy = default_resize_policy
        endif
        if (sp_matrix%mtype.eq.'row') then
            nlines = nrows
            length = ncols
        else
            nlines = ncols
            length = nrows
        endif
        allocate(sp_matrix%line(nlines))
        do i = 1, nlines
            call allocate_sparse_line(sp_matrix%line(i), &
                 line_size=sp_matrix%isize, length=length)
        enddo
        sp_matrix%nlines = nlines
    end subroutine allocate_sparse_matrix

    subroutine deallocate_sparse_matrix(sp_matrix)
        implicit none
        type(sparse_matrix), intent(inout)::sp_matrix
        integer(spip)::i
        sp_matrix%nrows = 0
        sp_matrix%ncols = 0
        sp_matrix%isize = 0
        sp_matrix%mtype = ""
        sp_matrix%storage = ""
        if (allocated(sp_matrix%resize_policy)) then
            deallocate(sp_matrix%resize_policy)
        endif
        if (allocated(sp_matrix%line)) then
            do i = 1, size(sp_matrix%line)
                call deallocate_sparse_line(sp_matrix%line(i))
            enddo
            deallocate(sp_matrix%line)
        endif
        sp_matrix%nlines = 0
    end subroutine deallocate_sparse_matrix

end module sparseset

