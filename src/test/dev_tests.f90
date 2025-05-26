! First tests (during development)

module dev_tests
    use sparseset
    implicit none


contains

    subroutine test_allocate_deallocate_line()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, correct = 0, err_stat

        call allocate_sparse_line(sp_line, 10, 100)
        write(*,"(a)",advance="no")"Testing line allocation (from scratch):"
        ! Testing if line size was correctly set
        tested = tested + 1
        if (sp_line%lsize.eq.10) correct = correct + 1
        ! Testing if length was correctly set
        tested = tested + 1
        if (sp_line%length.eq.100) correct = correct + 1
        ! Testing if lcount was correctly set to zero
        tested = tested + 1
        if (sp_line%lcount.eq.0) correct = correct + 1
        ! Testing if rpstage was correctly set to zero
        tested = tested + 1
        if (sp_line%rpstage.eq.0) correct = correct + 1
        ! Testing if assembed was correctly set to .false.
        tested = tested + 1
        if (.not.sp_line%assembled) correct = correct + 1
        ! Testing if lindex was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lindex).and.&
            (size(sp_line%lindex).eq.10)) correct = correct + 1
        ! Testing if lvalue was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lvalue).and.&
            (size(sp_line%lvalue).eq.10)) correct = correct + 1
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
        ! Now messing values and reallocating it
        sp_line%lcount = 8
        sp_line%rpstage = 3
        sp_line%assembled = .true.
        tested = 0; correct = 0
        call allocate_sparse_line(sp_line, 30, 300, err_stat)
        write(*,"(a)",advance="no")"Testing line allocation (reallocation):"
        ! Testing if stat is returning 0
        tested = tested + 1
        if (err_stat.eq.0) correct = correct + 1
        ! Testing if line size was correctly set
        tested = tested + 1
        if (sp_line%lsize.eq.30) correct = correct + 1
        ! Testing if length was correctly set
        tested = tested + 1
        if (sp_line%length.eq.300) correct = correct + 1
        ! Testing if lcount was correctly set to zero
        tested = tested + 1
        if (sp_line%lcount.eq.0) correct = correct + 1
        ! Testing if rpstage was correctly set to zero
        tested = tested + 1
        if (sp_line%rpstage.eq.0) correct = correct + 1
        ! Testing if assembed was correctly set to .false.
        tested = tested + 1
        if (.not.sp_line%assembled) correct = correct + 1
        ! Testing if lindex was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lindex).and.&
            (size(sp_line%lindex).eq.30)) correct = correct + 1
        ! Testing if lvalue was correcttly allocated
        tested = tested + 1
        if (allocated(sp_line%lvalue).and.&
            (size(sp_line%lvalue).eq.30)) correct = correct + 1
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
        ! Testing deallocation
        write(*,"(a)",advance="no")"Testing line deallocation:"
        tested = 0; correct = 0
        call deallocate_sparse_line(sp_line)
        ! Testing if line size was correctly set to zero
        tested = tested + 1
        if (sp_line%lsize.eq.0) correct = correct + 1
        ! Testing if length was correctly set to zero
        tested = tested + 1
        if (sp_line%length.eq.0) correct = correct + 1
        ! Testing if lcount was correctly set to zero
        tested = tested + 1
        if (sp_line%lcount.eq.0) correct = correct + 1
        ! Testing if rpstage was correctly set to zero
        tested = tested + 1
        if (sp_line%rpstage.eq.0) correct = correct + 1
        ! Testing if assembed was correctly set to .false.
        tested = tested + 1
        if (.not.sp_line%assembled) correct = correct + 1
        ! Testing if lindex was correcttly deallocated
        tested = tested + 1
        if (.not.allocated(sp_line%lindex)) correct = correct + 1
        ! Testing if lvalue was correcttly deallocated
        tested = tested + 1
        if (.not.allocated(sp_line%lvalue)) correct = correct + 1
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
        ! Testing if the allocation subroutine can identify wrong parameters
        write(*,"(a)",advance="no")"Testing attempts to use wrong parameters:"
        tested = 0; correct = 0
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,-10,100, err_stat)
        if (err_stat.eq.1) correct = correct + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,10,-100, err_stat)
        if (err_stat.eq.1) correct = correct + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,-10,-100, err_stat)
        if (err_stat.eq.1) correct = correct + 1
        err_stat = 0
        tested = tested + 1
        call allocate_sparse_line(sp_line,0,100, err_stat)
        if (err_stat.eq.1) correct = correct + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"

    end subroutine test_allocate_deallocate_line

    subroutine test_available_space()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, correct = 0, err_stat, i, num
        write(*,"(a)",advance="no")"Testing function to get available space in&
        a sparse line:"
        num = 10
        ! First 'test': allocation is OK!
        tested = tested + 1
        call allocate_sparse_line(sp_line,num,100,err_stat)
        if (err_stat.eq.0) correct = correct + 1
        ! First test: all terms are available to store values
        tested = tested + 1
        if (available_space(sp_line).eq.num) correct = correct + 1 
        ! Test availabe space while pushing terms one by one
        do i = 1,num
            tested = tested + 1
            sp_line%lcount = sp_line%lcount + 1
            sp_line%lindex(i) = i
            sp_line%lvalue(i) = dble(i)
            if (available_space(sp_line).eq.(num-i)) correct = correct + 1 
        enddo
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"

    end subroutine test_available_space

    subroutine test_pushing_to_line()
        implicit none
        type(sparse_line)::sp_line
        integer::tested = 0, correct = 0, err_stat, i, num
        write(*,"(a)",advance="no")"Testing push_term_to_line sub:"
        num = 10
        ! First 'test': allocation is OK!
        tested = tested + 1
        call allocate_sparse_line(sp_line,num,100,err_stat)
        if (err_stat.eq.0) correct = correct + 1
        ! First test: all terms are available to store values
        tested = tested + 1
        if (available_space(sp_line).eq.num) correct = correct + 1 
        ! Test availabe space while pushing terms one by one
        do i = 1,num
            if (mod(i,2).eq.0) then
                tested = tested + 1
                call push_term_to_line(sp_line,i,dble(i))
                if (available_space(sp_line).eq.(num-i)) correct = correct + 1
            else
                ! Test pushing using stat variable
                tested = tested + 1
                err_stat = 1
                call push_term_to_line(sp_line,i,dble(i),err_stat)
                if (available_space(sp_line).eq.(num-i)) correct = correct + 1
                tested = tested + 1
                if (err_stat.eq.0) correct = correct + 1
            endif
        enddo
        ! Testing if the subroutine detects error correctly
        err_stat = 0
        tested = tested + 1
        call push_term_to_line(sp_line,i,dble(i),err_stat)
        if (err_stat.eq.1) correct = correct + 1
        ! Clear line and test pushing terms
        sp_line%assembled = .true.
        call clear_sparse_line(sp_line)
        tested = tested + 1
        if (sp_line%lcount.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) correct = correct + 1
        ! Test pushing 3 terms at once in sp_line
        call push_terms_to_line(sp_line,(/1, 2, 3/),(/0.0_spdp,&
        -1.0_spdp, 2.0_spdp/))
        tested = tested + 1
        if (sp_line%lcount.eq.3) correct = correct + 1
        tested = tested + 1
        if (available_space(sp_line).eq.7) correct = correct + 1     
        ! Test pushing 2 terms at once in sp_line
        err_stat = 1
        call push_terms_to_line(sp_line,(/-1, -2/),(/0.1_spdp,&
        -1.5_spdp/),err_stat)
        tested = tested + 1
        if (sp_line%lcount.eq.5) correct = correct + 1
        tested = tested + 1
        if (available_space(sp_line).eq.5) correct = correct + 1     
        tested = tested + 1
        if (err_stat.eq.0) correct = correct + 1
        ! Testing if sub correctly identify invalid data size
        err_stat = 0
        call push_terms_to_line(sp_line,(/-1, -2, 0/),(/0.1_spdp,&
        -1.5_spdp/),err_stat)
        tested = tested + 1
        if (err_stat.eq.1) correct = correct + 1
        ! Testing pushing 1 term using push_terms sub
        err_stat = 1
        call push_terms_to_line(sp_line,(/5/),(/0.1_spdp/),&
        err_stat)
        tested = tested + 1
        if (sp_line%lcount.eq.6) correct = correct + 1
        tested = tested + 1
        if (available_space(sp_line).eq.4) correct = correct + 1     
        tested = tested + 1
        if (err_stat.eq.0) correct = correct + 1

        ! Testing if sub correctly identify no space for pushing data
        err_stat = 0
        call push_terms_to_line(sp_line,(/-1, -2, 0, 0, 0/),(/0.1_spdp,&
        -1.5_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp/),err_stat)
        tested = tested + 1
        if (err_stat.eq.1) correct = correct + 1
        tested = tested + 1
        if (available_space(sp_line).eq.4) correct = correct + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_pushing_to_line

    subroutine test_copy_line_terms()
        implicit none
        type(sparse_line):: oline, dline ! origin and destination lines
        integer:: correct, tested, err_stat, i
        write(*,"(a)",advance="no")"Testing copy_line_terms sub:"
        tested = 0; correct = 0
        call allocate_sparse_line(oline,10,20)
        call allocate_sparse_line(dline,5,10)
        call push_terms_to_line(oline,(/1, 2, 3, 4, 5, 6/),&
        (/1.0_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp/))
        ! Test if copy_line status identify error
        err_stat = 0
        call copy_sparse_line_terms(oline, dline, err_stat)
        tested = tested + 1
        if (err_stat.eq.1) correct = correct + 1
        ! Test if clear data is working
        ! Are lcount correct before clear?
        tested = tested + 1
        if (oline%lcount.eq.6) correct = correct + 1
        tested = tested + 1
        call clear_sparse_line(oline)
        if (oline%lcount.eq.0) correct = correct + 1
        ! Now the original line will have a legth so that it can
        ! be copied to dline
        call push_terms_to_line(oline,(/1, 2, 3, 4, 5/),&
        (/1.0_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp, 1.0_spdp/))
        ! Test if copy_line status error i zero now
        err_stat = 1
        call copy_sparse_line_terms(oline, dline, err_stat)
        tested = tested + 1
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (dline%lcount.eq.5) correct = correct + 1
        do i = 1, oline%lcount
            tested = tested + 1
            if (oline%lindex(i).eq.dline%lindex(i)) correct = correct + 1
            tested = tested + 1
            if (oline%lvalue(i).eq.dline%lvalue(i)) correct = correct + 1
        enddo

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_copy_line_terms

    subroutine test_quicksort()
        implicit none
        integer(spip), dimension(:), allocatable::in_arr, out_arr, exp_arr
        integer(spip)::tested, correct, i
        write(*,"(a)",advance="no")"Testing sorting function:"
        tested = 0
        correct = 0
        allocate(in_arr(5), out_arr(5), exp_arr(5))
        ! Reverse order, 5 terms
        in_arr = (/5, 4, 3, 2, 1/)
        exp_arr = (/5, 4, 3, 2, 1/)
        out_arr = sorted_indexes(in_arr)
        tested = tested + 1
        ! Already sorted, 5 terms
        if (all(exp_arr.eq.out_arr)) correct = correct + 1
        in_arr = (/1, 2, 3, 4, 5/)
        exp_arr = (/1, 2, 3, 4, 5/)
        out_arr = sorted_indexes(in_arr)
        tested = tested + 1
        if (all(exp_arr.eq.out_arr)) correct = correct + 1
        !  5 terms, no repeted value
        in_arr = (/3, 2, 5, 1, 4/)
        exp_arr = (/4, 2, 1, 5, 3/)
        out_arr = sorted_indexes(in_arr)
        tested = tested + 1
        if (all(exp_arr.eq.out_arr)) correct = correct + 1
        ! 5 terms with repetition
        in_arr = (/3, 5, 5, 1, 5/)
        exp_arr = (/4, 1, 2, 3, 5/)
        out_arr = sorted_indexes(in_arr)
        tested = tested + 1
        if (all(exp_arr.eq.out_arr)) correct = correct + 1
        ! Reverse order, 100 terms
        deallocate(in_arr, out_arr, exp_arr)
        allocate(in_arr(100), out_arr(100), exp_arr(100))
        do i=1,100
            in_arr(i) = 101 - i
            exp_arr(i) = 101 - i
        enddo
        out_arr = sorted_indexes(in_arr)
        tested = tested + 1
        if (all(exp_arr.eq.out_arr)) correct = correct + 1
        write(*,*)'in: ',in_arr
        write(*,*)'out: ',out_arr
        write(*,*)'exp: ',exp_arr
        write(*,*)(out_arr.eq.exp_arr)

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_quicksort

    subroutine perform_all_dev_tests()
        implicit none
        call test_allocate_deallocate_line()
        call test_available_space()
        call test_pushing_to_line()
        call test_copy_line_terms()
        call test_quicksort()
    end subroutine perform_all_dev_tests

end module dev_tests

