! First tests (during development)

module dev_tests
    use sparseset
    implicit none


contains

    subroutine test_allocate_line()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, right = 0

        call allocate_sparse_line(sp_line, 10, 100)
        write(*,"(a)",advance="no")"Testing line allocation:"
        ! Testing if line size was correctly set
        tested = tested + 1
        if (sp_line%lsize.eq.10) right = right + 1
        ! Testing if length was correctly set
        tested = tested + 1
        if (sp_line%length.eq.100) right = right + 1
        ! Testing if lcount was correctly set to zero
        tested = tested + 1
        if (sp_line%lcount.eq.0) right = right + 1
        ! Testing if rpstage was correctly set to zero
        tested = tested + 1
        if (sp_line%rpstage.eq.0) right = right + 1
        ! Testing if assembed was correctly set to .false.
        tested = tested + 1
        if (.not.sp_line%assembled) right = right + 1
        ! Testing if lindex was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lindex).and.&
            (size(sp_line%lindex).eq.10)) right = right + 1
        ! Testing if lvalue was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lvalue).and.&
            (size(sp_line%lvalue).eq.10)) right = right + 1


        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"
        


    end subroutine test_allocate_line

end module dev_tests

