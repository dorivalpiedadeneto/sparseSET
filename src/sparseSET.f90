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
        logical::sym
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
    ! sym - indicates if the matrix os symmetric (.true.) or not (.false.)
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
            if (.not.((line_size.gt.0).and.(length.gt.0).and.&
                (line_size.le.length))) then
                stat = -1
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

end module sparseset

