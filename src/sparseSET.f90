! sparseSET - A library for assembling sparse matrices in fortran
! Authors: Dorival Piedade Neto
!          Rodrigo Ribeiro Paccola
!          Rogerio Carrazedo
!          Version: 2.0 (2024)
! License: BSD 2-Clause License (FreeBSD/Simplified)
module sparseset
    implicit none
    
    ! Types of variables (hardcoded; user must change if needed)
    integer, parameter::smip = 4 ! integer precision for sparseSET library
    integer, parameter::smrp = 8 ! real precision in sparseSET library
    ! The defaults used for the library are: spip for integers, spdp for reals

    ! A line is used to store rows or columns of a sparse matrix
    type line
        integer(smip)::lsize = 0
        integer(smip)::lcount = 0
        integer(smip), dimension(:), allocatable::lindex
        real(smrp), dimension(:), allocatable::lvalue
        integer(smip)::rpstage = 0
        logical::assembled=.false.
    end type line

    ! About the line type
    ! Lines are data structures to hold rows or column data while the sparse
    ! matrix is being assembled. Terms are pushed to the end, and when all 
    ! contributions to the sparse matrix are done, them the terms are ordered
    ! and terms of the same index are summed up.
    ! The line data (index and value) are stored in one dimensional arrays.
    ! The initial space allocated to store these values depend on the system
    ! of equation characteristics. In some cases, terms of same index are
    ! pushed several times, resulting in excessive space usage e reallocation.
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

    type sparse_matrix
        character(3)::mtype = 'row'
        logical::sym = .true.
        integer(smip)::isize
        integer(smip),dimension(16)::resize_policy = (\ &
        0, 2, 0, 3, 0, -2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0\)
        integer(smip)::nlines
        integer(smpi)::nrows
        integer(smpi)::ncols
        type(lines), dimension(:), allocatable::line
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
    ! line size); 3 indites to resize to 3*isize, and so on.
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

    type CSC
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nz
        integer(smip), dimension(:), allocatable::row
        integer(smrp), dimension(:), allocatable::mvalue
    end type CSC

    type CSR
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nz
        integer(smip), dimension(:), allocatable::col
        integer(smrp), dimension(:), allocatable::mvalue
    end type CSR

    type triplet
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nnz
        integer(smip), dimension(:), allocatable::row
        integer(smip), dimension(:), allocatable::col
        integer(smrp), dimension(:), allocatable::mvalue
    end type triplet

end module sparseset

