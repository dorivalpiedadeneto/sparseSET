! First tests (during development)

module dev_tests
    use sparseset
    implicit none


contains

    subroutine test_allocate_deallocate_line()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, right = 0, err_stat

        call allocate_sparse_line(sp_line, 10, 100)
        write(*,"(a)",advance="no")"Testing line allocation (from scratch):"
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
        ! Now messing values and reallocating it
        sp_line%lcount = 8
        sp_line%rpstage = 3
        sp_line%assembled = .true.
        tested = 0; right = 0
        call allocate_sparse_line(sp_line, 30, 300, err_stat)
        write(*,"(a)",advance="no")"Testing line allocation (reallocation):"
        ! Testing if stat is returning 0
        tested = tested + 1
        if (err_stat.eq.0) right = right + 1
        ! Testing if line size was correctly set
        tested = tested + 1
        if (sp_line%lsize.eq.30) right = right + 1
        ! Testing if length was correctly set
        tested = tested + 1
        if (sp_line%length.eq.300) right = right + 1
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
            (size(sp_line%lindex).eq.30)) right = right + 1
        ! Testing if lvalue was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lvalue).and.&
            (size(sp_line%lvalue).eq.30)) right = right + 1
        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"
        ! Testing deallocation
        write(*,"(a)",advance="no")"Testing line deallocation:"
        tested = 0; right = 0
        call deallocate_sparse_line(sp_line)
        ! Testing if line size was correctly set to zero
        tested = tested + 1
        if (sp_line%lsize.eq.0) right = right + 1
        ! Testing if length was correctly set to zero
        tested = tested + 1
        if (sp_line%length.eq.0) right = right + 1
        ! Testing if lcount was correctly set to zero
        tested = tested + 1
        if (sp_line%lcount.eq.0) right = right + 1
        ! Testing if rpstage was correctly set to zero
        tested = tested + 1
        if (sp_line%rpstage.eq.0) right = right + 1
        ! Testing if assembed was correctly set to .false.
        tested = tested + 1
        if (.not.sp_line%assembled) right = right + 1
        ! Testing if lindex was correcttly deallocated
        tested = tested + 1
        if (.not.allocated(sp_line%lindex)) right = right + 1
        ! Testing if lvalue was correcttly deallocated
        tested = tested + 1
        if (.not.allocated(sp_line%lvalue)) right = right + 1
        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"
        ! Testing if the allocation subroutine can identify wrong parameters
        write(*,"(a)",advance="no")"Testing attempts to use wrong parameters:"
        tested = 0; right = 0
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,-10,100, err_stat)
        if (err_stat.eq.1) right = right + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,10,-100, err_stat)
        if (err_stat.eq.1) right = right + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,-10,-100, err_stat)
        if (err_stat.eq.1) right = right + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,0,100, err_stat)
        if (err_stat.eq.1) right = right + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"

    end subroutine test_allocate_deallocate_line

    subroutine test_available_space()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, right = 0, err_stat, i, num
        write(*,"(a)",advance="no")"Testing function to get available space in&
        a sparse line:"
        num = 10
        ! First 'test': allocation is OK!
        tested = tested + 1
        call allocate_sparse_line(sp_line,num,100,err_stat)
        if (err_stat.eq.0) right = right + 1
        ! First test: all terms are available to store values
        tested = tested + 1
        if (available_space(sp_line).eq.num) right = right + 1 
        ! Test availabe space while pushing terms one by one
        do i = 1,num
            tested = tested + 1
            sp_line%lcount = sp_line%lcount + 1
            sp_line%lindex(i) = i
            sp_line%lvalue(i) = dble(i)
            if (available_space(sp_line).eq.(num-i)) right = right + 1 
        enddo
        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"

    end subroutine test_available_space

    subroutine test_pushing_to_line()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, right = 0, err_stat, i, num
        write(*,"(a)",advance="no")"Testing push_term_to_line sub:"
        num = 10
        ! First 'test': allocation is OK!
        tested = tested + 1
        call allocate_sparse_line(sp_line,num,100,err_stat)
        if (err_stat.eq.0) right = right + 1
        ! First test: all terms are available to store values
        tested = tested + 1
        if (available_space(sp_line).eq.num) right = right + 1 
        ! Test availabe space while pushing terms one by one
        do i = 1,num
            if (mod(i,2).eq.0) then
                tested = tested + 1
                call push_term_to_line(sp_line,i,dble(i))
                if (available_space(sp_line).eq.(num-i)) right = right + 1
            else
                ! Test pushing using stat variable
                tested = tested + 1
                err_stat = 1
                call push_term_to_line(sp_line,i,dble(i),err_stat)
                if (available_space(sp_line).eq.(num-i)) right = right + 1
                tested = tested + 1
                if (err_stat.eq.0) right = right + 1
            endif
        enddo
        ! Testing if the subroutine detects error correctly
        err_stat = 0
        tested = tested + 1
        call push_term_to_line(sp_line,i,dble(i),err_stat)
        if (err_stat.eq.1) right = right + 1
        ! Clear line and test pushing terms
        sp_line%assembled = .true.
        call clear_sparse_line(sp_line)
        tested = tested + 1
        if (sp_line%lcount.eq.0) right = right + 1
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) right = right + 1
        ! Test pushing 3 terms at once in sp_line
        call push_terms_to_line(sp_line,(/1, 2, 3/),(/0.0_spdp,&
        -1.0_spdp, 2.0_spdp/))
        tested = tested + 1
        if (sp_line%lcount.eq.3) right = right + 1
        tested = tested + 1
        if (available_space(sp_line).eq.7) right = right + 1     
        ! Test pushing 2 terms at once in sp_line
        err_stat = 1
        call push_terms_to_line(sp_line,(/-1, -2/),(/0.1_spdp,&
        -1.5_spdp/),err_stat)
        tested = tested + 1
        if (sp_line%lcount.eq.5) right = right + 1
        tested = tested + 1
        if (available_space(sp_line).eq.5) right = right + 1     
        tested = tested + 1
        if (err_stat.eq.0) right = right + 1
        ! Testing if sub correctly identify invalid data size
        err_stat = 0
        call push_terms_to_line(sp_line,(/-1, -2, 0/),(/0.1_spdp,&
        -1.5_spdp/),err_stat)
        tested = tested + 1
        if (err_stat.eq.1) right = right + 1
        ! Testing pushing 1 term using push_terms sub
        err_stat = 1
        call push_terms_to_line(sp_line,(/5/),(/0.1_spdp/),&
        err_stat)
        tested = tested + 1
        if (sp_line%lcount.eq.6) right = right + 1
        tested = tested + 1
        if (available_space(sp_line).eq.4) right = right + 1     
        tested = tested + 1
        if (err_stat.eq.0) right = right + 1

        ! Testing if sub correctly identify no space for pushing data
        err_stat = 0
        call push_terms_to_line(sp_line,(/-1, -2, 0, 0, 0/),(/0.1_spdp,&
        -1.5_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp/),err_stat)
        tested = tested + 1
        if (err_stat.eq.1) right = right + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",right,"/",tested,"]"
    end subroutine test_pushing_to_line

    subroutine perform_all_dev_tests()
        implicit none
        call test_allocate_deallocate_line()
        call test_available_space()
        call test_pushing_to_line()
    end subroutine perform_all_dev_tests

end module dev_tests

