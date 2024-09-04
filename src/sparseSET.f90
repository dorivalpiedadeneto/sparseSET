! sparseSET - A library for assembling sparse matrices in fortran
! Authors: Dorival Piedade Neto
!          Rodrigo Ribeiro Paccola
!          Rogerio Carrazedo
!          Version: 2.0 (2024)
! License: BSD 2-Clause License (FreeBSD/Simplified)
module sparseset
    implicit none
    
    ! Types of variables (hardcoded; user must change if needed)
    integer, parameter::smip = 4 !integer precision for sparseSET library
    integer, parameter::smlp = 8 !long integer precision for sparseSET library
    integer, parameter::smsp = 4 !single precision real for sparseSET library
    integer, parameter::smdp = 8 !double precision real for sparseSET library
    ! The defaults used for the library are: spip for integers, spdp for reals

    ! A line is used to store rows or columns of a sparse matrix
    type line
        integer(smip)::lsize = 0
        integer(smip)::lcount = 0
        integer(smip), dimension(:), allocatable::lindex
        real(smdp), dimension(:), allocatable::lvalue
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
        character(3)::mtype
        logical::sym
        integer(smip)::nlines
        type(lines), dimension(:), allocatable::line
    end type sparse_matrix

    type CSC
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nz
        integer(smip), dimension(:), allocatable::row
        integer(smdp), dimension(:), allocatable::mvalue
    endtype CSC

    type CSR
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nz
        integer(smip), dimension(:), allocatable::col
        integer(smdp), dimension(:), allocatable::mvalue
    end type CSR

    type triplet
        integer(smip)::msize
        integer(smip), dimension(:), allocatable::nnz
        integer(smip), dimension(:), allocatable::row
        integer(smip), dimension(:), allocatable::col
        integer(smdp), dimension(:), allocatable::mvalue
    end type triplet

end module sparseset

