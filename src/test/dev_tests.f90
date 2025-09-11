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
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_quicksort

    subroutine test_assemble_line()
        implicit none
        type(sparse_line)::sp_line
        integer:: correct, tested, i
        real(spdp)::tol=1.0e-6
        real(spdp),dimension(10)::ex_vals
        call allocate_sparse_line(sp_line, 10, 10)
        write(*,"(a)",advance="no")"Testing assemble line sub:"
        tested = 0; correct = 0
        call push_terms_to_line(sp_line, (/3,1,2,1,2,4,5,2/), &
        (/0.5_spdp, 1.0_spdp, 0.2_spdp, 1.0_spdp, 0.2_spdp,-2.0_spdp, &
        3.0_spdp, 0.2_spdp/))
        tested = tested + 1
        if (sp_line%lcount.eq.8) correct = correct + 1
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) correct = correct + 1
        call assemble_sparse_line(sp_line)
        tested = tested + 1
        if (sp_line%lcount.eq.5) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:sp_line%lcount).eq.(/1,2,3,4,5/)))&
            correct = correct + 1
        ! For testing the floating point values, use tolerance
        ! (testing one by one, at least for now)
        ex_vals(1:5) = (/2.0_spdp, 0.6_spdp, 0.5_spdp, -2.0_spdp, 3.0_spdp/)
        do i = 1, 5
            tested = tested + 1
            if ((sp_line%lvalue(i) - ex_vals(i))**2.lt.tol) &
                correct = correct + 1
        enddo
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        ! Testing if pushing after assembling is correctly done
        call push_terms_to_line(sp_line,(/3,2,9,8/),&
        (/0.5_spdp, 0.4_spdp, -1.0_spdp, 2.3_spdp/))
        tested = tested + 1
        if (sp_line%lcount.eq.9) correct = correct + 1
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) correct = correct + 1
        call assemble_sparse_line(sp_line)
        tested = tested + 1
        if (sp_line%lcount.eq.7) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:sp_line%lcount).eq.(/1,2,3,4,5,8,9/)))&
            correct = correct + 1
        ! For testing the floating point values, use tolerance
        ! (testing one by one, at least for now)
        ex_vals(1:7) = (/2.0_spdp, 1.0_spdp, 1.0_spdp, -2.0_spdp,&
        3.0_spdp, 2.3_spdp, -1.0_spdp/)
        do i = 1, 7
            tested = tested + 1
            if ((sp_line%lvalue(i) - ex_vals(i))**2.lt.tol) &
                correct = correct + 1
        enddo
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        ! Now test adding one term
        call push_term_to_line(sp_line,10,-10.0_spdp)
        tested = tested + 1
        if (sp_line%lcount.eq.8) correct = correct + 1
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) correct = correct + 1
        call assemble_sparse_line(sp_line)
        tested = tested + 1
        if (sp_line%lcount.eq.8) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:sp_line%lcount).eq.(/1,2,3,4,5,8,9,10/)))&
            correct = correct + 1
        ! For testing the floating point values, use tolerance
        ! (testing one by one, at least for now)
        ex_vals(1:8) = (/2.0_spdp, 1.0_spdp, 1.0_spdp, -2.0_spdp,&
        3.0_spdp, 2.3_spdp, -1.0_spdp, -10.0_spdp/)
        do i = 1, 8
            tested = tested + 1
            if ((sp_line%lvalue(i) - ex_vals(i))**2.lt.tol) &
                correct = correct + 1
        enddo
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_assemble_line

    subroutine test_search_line()
        implicit none
        type(sparse_line)::sp_line
        integer(spip)::tested, correct, pos, ind, err_stat, i
        real(spdp)::tol=1.0e-6_spdp
        write(*,"(a)",advance="no")"Testing search_line function:"
        tested = 0
        correct = 0
        tested = tested + 1
        pos = search_line(sp_line,1)
        if (pos.eq.-1) correct = correct + 1
        call allocate_sparse_line(sp_line, 20, 100)
        tested = tested + 1
        pos = search_line(sp_line,1)
        if (pos.eq.-1) correct = correct + 1
        call push_terms_to_line(sp_line,(/17, 21, 12, 15, 12, 25, 27, 17/),&
        (/9.0_spdp, 21.0_spdp, 6.0_spdp, 15.0_spdp, 6.0_spdp, 25.0_spdp,&
        27.0_spdp, 8.0_spdp/), err_stat)
        ! Does the pushing went right?
        tested = tested + 1
        if (err_stat.eq.0) correct = correct + 1
        ! Sparse line must be not assembled yet
        tested = tested + 1
        if (sp_line%assembled.eqv..false.) correct = correct + 1
        ! Searching must assemble automatically (even if ind is not found)
        tested = tested + 1
        pos = search_line(sp_line,1)
        if (pos.eq.-1) correct = correct + 1
        ! Is it assembled now (it must be!)  
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        do i=1,11
            tested = tested + 1
            pos = search_line(sp_line,i)
            if (pos.eq.-1) correct = correct + 1
        enddo
        ! The first term must be 12
        tested = tested + 1
        pos = search_line(sp_line, 12)
        if (pos.eq.1) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-12.0_spdp).lt.tol) correct = correct + 1
        ! 13 and 14 are not in the line
        tested = tested + 1
        pos = search_line(sp_line, 13)
        if (pos.eq.-1) correct = correct + 1
        tested = tested + 1
        pos = search_line(sp_line, 14)
        if (pos.eq.-1) correct = correct + 1
        ! The second term must be 15
        tested = tested + 1
        pos = search_line(sp_line, 15)
        if (pos.eq.2) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-15.0_spdp).lt.tol) correct = correct + 1
        ! 16 is not in the line
        tested = tested + 1
        pos = search_line(sp_line, 16)
        if (pos.eq.-1) correct = correct + 1
        ! The third term must be 17
        tested = tested + 1
        pos = search_line(sp_line, 17)
        if (pos.eq.3) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-17.0_spdp).lt.tol) correct = correct + 1
        ! 18 to 20 and 22 to 24 are not in the line 
        do i=18, 20
            tested = tested + 1
            pos = search_line(sp_line,i)
            if (pos.eq.-1) correct = correct + 1
        enddo
        do i=22, 24
            tested = tested + 1
            pos = search_line(sp_line,i)
            if (pos.eq.-1) correct = correct + 1
        enddo
        ! The 4th term must be 21
        tested = tested + 1
        pos = search_line(sp_line, 21)
        if (pos.eq.4) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-21.0_spdp).lt.tol) correct = correct + 1
        ! The 5th term must be 25
        tested = tested + 1
        pos = search_line(sp_line, 25)
        if (pos.eq.5) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-25.0_spdp).lt.tol) correct = correct + 1
        ! The 6th term must be 27
        tested = tested + 1
        pos = search_line(sp_line, 27)
        if (pos.eq.6) correct = correct + 1
        tested = tested + 1
        if (dabs(sp_line%lvalue(pos)-27.0_spdp).lt.tol) correct = correct + 1
        ! 26 is not in the line
        tested = tested + 1
        pos = search_line(sp_line, 26)
        if (pos.eq.-1) correct = correct + 1
        ! 28 and more are not in the line 
        do i=28, 31
            tested = tested + 1
            pos = search_line(sp_line,i)
            if (pos.eq.-1) correct = correct + 1
        enddo
 
 
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_search_line

    subroutine test_sparse_line_to_array()
        implicit none
        type(sparse_line)::sp_line
        real(spdp), dimension(:), allocatable:: arr, expected
        integer(spip)::err_stat, tested, correct
        real(spdp)::tol=1.0e-8_spdp
        tested = 0; correct = 0
        write(*,"(a)",advance="no")"Testing sparse_line_to_array:"
        tested = tested + 1
        call allocate_sparse_line(sp_line, 20, 10, err_stat)
        if (err_stat.eq.0) correct = correct + 1
        ! Testing for  a 'unassebled' sparse_line
        call push_terms_to_line(sp_line,(/2, 5, 3/),(/2.0_spdp, -5.0_spdp, &
        3.0_spdp/))
        call sparse_line_to_array(sp_line, arr)
        allocate(expected(10)); expected = 0.0_spdp
        expected(2) = 2.0_spdp; expected(3) = 3.0_spdp; expected(5) = -5.0_spdp
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1
        ! Testing if the result is correct (even if the line is not assembled)
        call push_terms_to_line(sp_line,(/2, 5, 3/),(/2.0_spdp, -5.0_spdp, &
        3.0_spdp/))
        call push_terms_to_line(sp_line,(/2, 5, 3/),(/2.0_spdp, -5.0_spdp, &
        3.0_spdp/))
        call push_terms_to_line(sp_line,(/2, 5, 3/),(/2.0_spdp, -5.0_spdp, &
        3.0_spdp/))
        expected = 4 * expected
        call sparse_line_to_array(sp_line, arr)
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1
        ! More terms...
        call push_terms_to_line(sp_line,(/1, 5, 8/),(/-1.0_spdp, -5.0_spdp, &
        8.8_spdp/))
        expected(1) = -1.0_spdp; expected(5) = expected(5) - 5.0_spdp
        expected(8) = 8.8_spdp
        call sparse_line_to_array(sp_line, arr)
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1
        ! Assembled must return the same result!
        call assemble_sparse_line(sp_line)
        call sparse_line_to_array(sp_line, arr)
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1
        ! Now, testing if changing the array size it still work correctly
        call deallocate_sparse_line(sp_line)
        call allocate_sparse_line(sp_line, 20, 100, err_stat)
        deallocate(expected)
        allocate(expected(100))
        expected = 0.0_spdp
        call push_terms_to_line(sp_line,(/8, 12, 33, 78, 99/),(/1.0_spdp,&
        1.2_spdp, -3.3e-2_spdp, 7.8e3_spdp, 9.99_spdp/))
        expected(8) = 1.0_spdp; expected(12) = 1.2_spdp
        expected(33) = -3.3e-2_spdp; expected(78) = 7.8e3_spdp
        expected(99) = 9.99_spdp
        call sparse_line_to_array(sp_line, arr)
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1
        call assemble_sparse_line(sp_line)
        tested = tested + 1
        if (all(dabs(arr-expected).lt.tol)) correct = correct + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_sparse_line_to_array

    subroutine test_array_to_sparse_line()
        implicit none
        type(sparse_line):: sp_line
        real(spdp), dimension(:), allocatable:: arr
        real(spdp)::tol
        integer(spdp)::tested, correct
        allocate(arr(10))
        arr = 0.0_spdp
        arr(1:10:2) = 1.0_spdp
        call array_to_sparse_line(arr, sp_line)
        write(*,"(a)",advance="no")"Testing array_to_sparse_line:"
        tested = 0; correct = 0
        ! Basic tests
        tested = tested + 1
        if (sp_line%lcount.eq.5) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:5).eq.(/1,3,5,7,9/))) correct = correct + 1
        tested = tested + 1
        tol = 1.0e-8_spdp
        if (all(dabs(sp_line%lvalue(1:5)-1.0_spdp).lt.tol))&
            correct = correct+ 1
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        tested = tested + 1
        if (sp_line%length.eq.10) correct = correct + 1
        ! Testing tolerance
        arr = 1.0e-9_spdp
        arr(1:10:2) = 1.0_spdp
        call array_to_sparse_line(arr, sp_line, tol)
        tested = tested + 1
        if (sp_line%lcount.eq.5) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:5).eq.(/1,3,5,7,9/))) correct = correct + 1
        tested = tested + 1
        if (all(dabs(sp_line%lvalue(1:5)-1.0_spdp).lt.tol))&
            correct = correct+ 1
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        tested = tested + 1
        if (sp_line%length.eq.10) correct = correct + 1
        ! Still testing tolerance
        arr = 1.0e-7_spdp
        arr(1:10:2) = 1.0_spdp
        call array_to_sparse_line(arr, sp_line, tol)
        tested = tested + 1
        if (sp_line%lcount.eq.10) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:10).eq.(/1,2,3,4,5,6,7,8,9,10/)))&
            correct = correct + 1
        tested = tested + 1
        if (all(dabs(sp_line%lvalue(1:10:2)-1.0_spdp).lt.tol))&
            correct = correct+ 1
        tested = tested + 1
        if (all(dabs(sp_line%lvalue(2:10:2)-1.0e-7_spdp).lt.tol**2))&
            correct = correct+ 1
        tested = tested + 1
        if (sp_line%assembled.eqv..true.) correct = correct + 1
        tested = tested + 1
        if (sp_line%length.eq.10) correct = correct + 1
        
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_array_to_sparse_line

    subroutine test_allocate_deallocate_sparse_matrix()
        implicit none
        type(sparse_matrix):: sp_matrix
        integer::tested = 0, correct = 0, err_stat
        write(*,"(a)",advance="no")"Testing allocate_sparse_matrix:"
        ! Testing if subroutine correctly identifies bad nrows and ncols values
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = -1, ncols = 1,&
        stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 1, ncols = -1,&
        stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = -1, ncols = -1,&
        stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        isize = -1, stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        mtype='cow', stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        mtype='rol', stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        storage='uppeR', stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        storage='lowEr', stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 10,&
        storage='FulL', stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        ! Test default values
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 20,&
        stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%isize.eq.default_isize) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%mtype.eq.default_mtype) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%nlines.eq.10) correct = correct + 1
        tested = tested + 1
        if (allocated(sp_matrix%line)) correct = correct + 1
        tested = tested + 1
        if (size(sp_matrix%line).eq.10) correct = correct + 1
        ! Testing if first line has the right length
        tested = tested + 1
        if (sp_matrix%line(1)%length.eq.20) correct = correct + 1
        tested = tested + 1
        if (all(sp_matrix%resize_policy.eq.default_resize_policy)) &
            correct = correct + 1
        !Same matrix, col mtype
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = 10, ncols = 20,&
        mtype='col',stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%mtype.eq.'col') correct = correct + 1
        tested = tested + 1
        if (sp_matrix%nlines.eq.20) correct = correct + 1
        tested = tested + 1
        if (allocated(sp_matrix%line)) correct = correct + 1
        tested = tested + 1
        if (size(sp_matrix%line).eq.20) correct = correct + 1
        ! Testing if first line has the right length
        tested = tested + 1
        if (sp_matrix%line(1)%length.eq.10) correct = correct + 1
        ! Testing deallocate
        call deallocate_sparse_matrix(sp_matrix)
        tested = tested + 1
        if (sp_matrix%nrows.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%ncols.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%nlines.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%mtype.eq."") correct = correct + 1
        tested = tested + 1
        if (sp_matrix%storage.eq."") correct = correct + 1
        tested = tested + 1
        if (.not.allocated(sp_matrix%resize_policy)) correct = correct + 1
        tested = tested + 1
        if (.not.allocated(sp_matrix%line)) correct = correct + 1

        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_allocate_deallocate_sparse_matrix

    subroutine test_resize_sparse_line()
        implicit none
        type(sparse_line):: sp_line
        integer(spip):: tested = 0, correct = 0, err_stat, as

        ! Allocating and pushing terms
        write(*,"(a)",advance="no")"Testing resize sparse line:"
        ! Testing if line size was correctly set
        tested = tested + 1
        call allocate_sparse_line(sp_line, 10, 100)
        call push_terms_to_line(sp_line,(/1, 2, 3, 4, 5, 6/),&
        (/1.0_spdp, 2.0_spdp, 3.0_spdp, 4.0_spdp, 5.0_spdp, 6.0_spdp/),&
        err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        call resize_sparse_line(sp_line,(/2, -2, 0/), 10, as, err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (as.eq.14) correct = correct + 1
        tested = tested + 1
        if (sp_line%rpstage.eq.1) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:6).eq.(/1, 2, 3, 4, 5, 6/))) then
            correct = correct + 1
        endif
        tested = tested + 1
        if (all(sp_line%lvalue(1:6).eq.(/1.0_spdp, 2.0_spdp, 3.0_spdp,&
            4.0_spdp, 5.0_spdp, 6.0_spdp/))) then
            correct = correct + 1
        endif
        tested = tested + 1
        if (size(sp_line%lindex).eq.20) correct = correct + 1
        tested = tested + 1
        if (size(sp_line%lvalue).eq.20) correct = correct + 1
        ! Pushing again
        tested = tested + 1
        call push_terms_to_line(sp_line,(/1, 2, 3, 4, 5, 6/),&
        (/1.0_spdp, 2.0_spdp, 3.0_spdp, 4.0_spdp, 5.0_spdp, 6.0_spdp/),&
        err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_line%lsize.eq.20) correct = correct + 1
        tested = tested + 1
        if (sp_line%lcount.eq.12) correct = correct + 1
        tested = tested + 1
        if (available_space(sp_line).eq.8) correct = correct + 1
        ! Resizing again (using current size as reference to resize)
        tested = tested + 1
        call resize_sparse_line(sp_line,(/2, -2, 0/), 10, as, err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (as.eq.28) correct = correct + 1
        tested = tested + 1
        if (sp_line%lsize.eq.40) correct = correct + 1
        tested = tested + 1
        if (sp_line%lcount.eq.12) correct = correct + 1
        tested = tested + 1
        if (sp_line%rpstage.eq.2) correct = correct + 1
        ! Now trying to resize (resize_policy defines assembling instead)
        tested = tested + 1
        call resize_sparse_line(sp_line,(/2, -2, 0/), 10, as, err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (as.eq.34) correct = correct + 1
        tested = tested + 1
        if (sp_line%lsize.eq.40) correct = correct + 1
        tested = tested + 1
        if (sp_line%lcount.eq.6) correct = correct + 1
        tested = tested + 1
        if (sp_line%rpstage.eq.3) correct = correct + 1
        tested = tested + 1
        if (all(sp_line%lindex(1:6).eq.(/1, 2, 3, 4, 5, 6/))) then
            correct = correct + 1
        endif
        tested = tested + 1
        if (all(sp_line%lvalue(1:6).eq.(/2.0_spdp, 4.0_spdp, 6.0_spdp,&
            8.0_spdp, 10.0_spdp, 12.0_spdp/))) then
            correct = correct + 1
        endif
        ! Trying to resize again
        call resize_sparse_line(sp_line,(/2, -2, 0/), 10, as, err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (as.eq.34) correct = correct + 1
        tested = tested + 1
        if (sp_line%lsize.eq.40) correct = correct + 1
        tested = tested + 1
        if (sp_line%lcount.eq.6) correct = correct + 1
        tested = tested + 1
        if (sp_line%rpstage.eq.3) correct = correct + 1
        tested = tested + 1
 
        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_resize_sparse_line

    subroutine test_resize_sparse_matrix()
        implicit none
        type(sparse_matrix)::sp_matrix
        integer(spip):: tested = 0, correct = 0, err_stat
        integer(spip):: nr = 10, nc = 20, is = 5
        integer(spip):: i

        write(*,"(a)",advance="no")"Testing resize sparse line:"
        tested = tested + 1
        call allocate_sparse_matrix(sp_matrix, nrows = nr,&
        ncols = nc, isize = is, mtype = "row", storage = "full",&
        resize_policy = (/2, 4, 0/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        do i=1, nr
            tested = tested + 1
            if (sp_matrix%line(i)%lsize.eq.is) correct = correct + 1
            tested = tested + 1
            if (sp_matrix%line(i)%length.eq.nc) correct = correct + 1
        enddo
        ! Test if invalid input is correctly identified
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/2,2,2,2/), stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
        ! Test case in which no lines are resized
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/2,2,2/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.is) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.is) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.is) correct = correct + 1
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/2,6,2/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.is) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.is) correct = correct + 1
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/6,6,6/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/6,10,6/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/6,11,6/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.(4*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/6,11,6/), stat = err_stat)
        if (err_stat.eq.0) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(1)%lsize.eq.(2*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(2)%lsize.eq.(4*is)) correct = correct + 1
        tested = tested + 1
        if (sp_matrix%line(3)%lsize.eq.(2*is)) correct = correct + 1
        ! Now let's test if it correctly identifies a fail in resize
        tested = tested + 1
        call resize_sparse_matrix(sp_matrix, indexes = (/1, 2, 3/),&
        necessary_size = (/6,21,6/), stat = err_stat)
        if (err_stat.ne.0) correct = correct + 1
 


        write(*,'(a,i2,a,i2,a)')" Passed [",correct,"/",tested,"]"
    end subroutine test_resize_sparse_matrix

    subroutine perform_all_dev_tests()
        implicit none
        call test_allocate_deallocate_line()
        call test_available_space()
        call test_pushing_to_line()
        call test_copy_line_terms()
        call test_quicksort()
        call test_assemble_line()
        call test_search_line()
        call test_sparse_line_to_array()
        call test_array_to_sparse_line()
        call test_allocate_deallocate_sparse_matrix()
        call test_resize_sparse_line()
        call test_resize_sparse_matrix()
    end subroutine perform_all_dev_tests

end module dev_tests

