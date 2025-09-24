    
    #include "definitions.h"

! --------------------------------------------------------------
!>        Solve Gauss Seidel in Blocks (Common)
! --------------------------------------------------------------
subroutine block_Point_Gauss_Seidel(variable, rhs, gamma)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start,j_end
    implicit none
    
    real(kind=prcs_var), intent(inout)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)     :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(in)     :: gamma            ! Coeff of U_p
    integer :: i, j, k

    k = 1
    !> Update variable at new step
    do j = j_start, j_end
        do i = i_start, i_end
            
            #IFDEF SQ_CELL
                variable(i,j,k) = (rhs(i,j,k) &
                                    - (variable(i,j-1,k)+variable(i,j+1,k)) &
                                    - (variable(i-1,j,k)+variable(i+1,j,k))  )/gamma
            #ELSE
                variable(i,j,k) = (rhs(i,j,k) &
                                    - (variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                                    - (variable(i-1,j,k)+variable(i+1,j,k))  )/gamma
            #ENDIF

        end do
    end do

end subroutine block_Point_Gauss_Seidel

subroutine block_Point_Gauss_Seidel_body(variable, rhs, gamma, i_range,flag)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start,j_end
    implicit none
    
    real(kind=prcs_var), intent(inout)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)     :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(in)     :: gamma            ! Coeff of U_p
    integer,             intent(in)     :: i_range(0:M+2, N+2, 1)    ! RHS of equation
    integer,             intent(in)     :: flag
    integer :: i, j, k, i1, l

    k = 1
    !> Update variable at new step
    do j = j_start, j_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k)

            
                #IFDEF SQ_CELL
                    variable(i,j,k) = (rhs(i,j,k) &
                                        - (variable(i,j-1,k)+variable(i,j+1,k)) &
                                        - (variable(i-1,j,k)+variable(i+1,j,k))  )/gamma
                #ELSE
                    variable(i,j,k) = (rhs(i,j,k) &
                                        - (variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                                        - (variable(i-1,j,k)+variable(i+1,j,k))  )/gamma
              #ENDIF

              !if(flag==1) then
               ! print*,"Cell(i,j):",i,j,"u(i,j):",variable(i,j,k)
                ! else if(flag==2) then 
                !print*,"Cell(i,j):",i,j,"v(i,j):",variable(i,j,k)   
                !else if(flag==3) then
                !print*,"Cell(i,j):",i,j,"p(i,j):",variable(i,j,k)    
               !end if 

            end do
        end do
    end do

end subroutine block_Point_Gauss_Seidel_body

! --------------------------------------------------------------
!>        Solve Weighted Jacobi in Blocks (Common)
! --------------------------------------------------------------
subroutine block_Weighted_Jacobi(variable, rhs, gamma, weight)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start,j_end
    implicit none
    
    real(kind=prcs_var), intent(inout)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)     :: rhs(M+2, N+2, 1)     ! RHS of equation
    real(kind=prcs_var), intent(in)     :: gamma    ! Coeff of U_p
    real(kind=prcs_var), intent(in)     :: weight   ! Weightage of updation


    integer, dimension(:) :: i(i_end - i_start + 1), j(j_end - j_start + 1), k(1)
    integer :: i_range, j_range, k_range
    real(kind=prcs_var) :: variable_new(M+2, N+2, 1)  ! Temp_array for variable_old

    ! allocate(i(i_end - i_start + 1))
    ! allocate(j(j_end - j_start + 1))
    ! allocate(k(1))
    i = [(i_range, i_range = i_start, i_end)]
    j = [(i_range, j_range = j_start, j_end)]
    k = [(k_range, k_range = 1, 1)]

    !> Update variable at intermediate step
    #IFDEF SQ_CELL
        variable_new(i,j,k) = (rhs(i,j,k) &
                        - (variable(i  ,j-1,k)+variable(i  ,j+1,k)) &
                        - (variable(i-1,j  ,k)+variable(i+1,j  ,k))  )/gamma
    #ELSE
        variable_new(i,j,k) = (rhs(i,j,k) &
                        - (variable(i  ,j-1,k)+variable(i  ,j+1,k)) * beta_sq &
                        - (variable(i-1,j  ,k)+variable(i+1,j  ,k))  )/gamma
    #ENDIF

    ! Update the variable based on the weight for new iteration
    variable(i,j,k) = (1.0_prcs_var-weight)*variable(i,j,k) + (weight)*variable_new(i,j,k)

end subroutine block_Weighted_Jacobi

! --------------------------------------------------------------
!>        Solve Multi-Grid in Blocks (Pressure)
! --------------------------------------------------------------
subroutine MG_restrict(tree, idx, mg_lvl_no)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module

    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: idx
    integer, intent(in)         :: mg_lvl_no

    integer :: i, j, k, i2, j2
    integer :: i_prolong_start, i_prolong_end, i_subtract, j_prolong_start, j_prolong_end, j_subtract
    integer :: child_idx, prev_mg_lvl, child_i
    real(kind=prcs_var) :: residual_1, residual_2, residual_3, residual_4
    integer :: tree_lvl_no

    tree_lvl_no = tree%blocks(idx)%lvl

    prev_mg_lvl = mg_lvl_no-1

    ! Blocks which are not refined need to inject values from PREVIOUS MG level of SAME block
    if (.not. tree%blocks(idx)%is_refined) then
        k = 1
        do j = j_start, j_end
            do i = i_start, i_end
                tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) = tree%blocks(idx)%multi_grid(prev_mg_lvl)%r(i,j,k)
                tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) * dx_sq_lvl(tree_lvl_no)
            end do
        end do

    ! Blocks which are refined need to restrict values from PREVIOUS MG level of CHILD blocks
    else

        do child_i = 1, 4

            if (child_i == 1) then
                i_prolong_start = i_start
                i_prolong_end   = i_mid
                i_subtract      = 1
                j_prolong_start = j_start
                j_prolong_end   = j_mid
                j_subtract      = 1
    
            else if (child_i == 2) then
                i_prolong_start = i_start
                i_prolong_end   = i_mid
                i_subtract      = 1
                j_prolong_start = j_mid+1
                j_prolong_end   = j_end
                j_subtract      = j_mid
                
            else if (child_i == 3) then
                i_prolong_start = i_mid+1
                i_prolong_end   = i_end
                i_subtract      = i_mid
                j_prolong_start = j_start
                j_prolong_end   = j_mid
                j_subtract      = 1
                
            else if (child_i == 4) then
                i_prolong_start = i_mid+1
                i_prolong_end   = i_end
                i_subtract      = i_mid
                j_prolong_start = j_mid+1
                j_prolong_end   = j_end
                j_subtract      = j_mid
                
            end if

            child_idx = tree%blocks(idx)%ptr_children(child_i)
    
            k = 1
            do j = j_prolong_start, j_prolong_end
                j2 = 2*(j-j_subtract)
                do i = i_prolong_start, i_prolong_end
                    i2 = 2*(i-i_subtract)
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) = &
                        (tree%blocks(child_idx)%multi_grid(prev_mg_lvl)%r(i2  ,j2  ,k) + &
                         tree%blocks(child_idx)%multi_grid(prev_mg_lvl)%r(i2+1,j2  ,k) + &
                         tree%blocks(child_idx)%multi_grid(prev_mg_lvl)%r(i2  ,j2+1,k) + &
                         tree%blocks(child_idx)%multi_grid(prev_mg_lvl)%r(i2+1,j2+1,k) ) * 0.25_prcs_var
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) = &
                        tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i,j,k) * dx_sq_lvl(tree_lvl_no)
                end do
            end do
            

        end do

    end if

end subroutine MG_restrict

subroutine MG_prolongate(tree, idx, prolong_max, mg_lvl_no)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module

    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: idx
    integer, intent(in)         :: mg_lvl_no
    real(kind=prcs_var), intent(inout) :: prolong_max

    integer :: i, j, k, i2, j2
    integer :: parent_idx, next_mg_lvl, own_lvl, max_lvl
    integer :: j_prolong_start, j_prolong_end, j_subtract
    integer :: i_prolong_start, i_prolong_end, i_subtract
    real(kind=prcs_var) :: p_prolong_temp
    real(kind=prcs_var) :: value_1, value_2, value_3, value_4, value_sum

    next_mg_lvl = mg_lvl_no+1
    own_lvl     = tree%blocks(idx)%lvl
    max_lvl     = tree%no_levels
    parent_idx  = tree%blocks(idx)%ptr_parent

    ! Blocks which were not refined need to inject values from NEXT MG level of SAME block
    if (own_lvl < (max_lvl-mg_lvl_no+1)) then
        k = 1
        do j = j_start, j_end
            do i = i_start, i_end
                value_sum = tree%blocks(idx)%multi_grid(next_mg_lvl)%p(i,j,k)
                tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i,j,k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no  )%p(i,j,k) + &
                    value_sum
                if (value_sum>prolong_max) prolong_max = value_sum
            end do
        end do

    ! Blocks which were refined need to prolongate values from NEXT MG level of PARENT block
    ! Ghost cells of finer grid is filled, (to be ignored)
    else
        
        if (idx == tree%blocks(parent_idx)%ptr_children(1)) then
            i_prolong_start = i_start
            i_prolong_end   = i_mid+1
            i_subtract      = 1
            j_prolong_start = j_start
            j_prolong_end   = j_mid+1
            j_subtract      = 1

        else if (idx == tree%blocks(parent_idx)%ptr_children(2)) then
            i_prolong_start = i_start
            i_prolong_end   = i_mid+1
            i_subtract      = 1
            j_prolong_start = j_mid+1
            j_prolong_end   = j_end+1
            j_subtract      = j_mid
            
        else if (idx == tree%blocks(parent_idx)%ptr_children(3)) then
            i_prolong_start = i_mid+1
            i_prolong_end   = i_end+1
            i_subtract      = i_mid
            j_prolong_start = j_start
            j_prolong_end   = j_mid+1
            j_subtract      = 1
            
        else if (idx == tree%blocks(parent_idx)%ptr_children(4)) then
            i_prolong_start = i_mid+1
            i_prolong_end   = i_end+1
            i_subtract      = i_mid
            j_prolong_start = j_mid+1
            j_prolong_end   = j_end+1
            j_subtract      = j_mid
            
        end if


        k = 1
        do j = j_prolong_start, j_prolong_end
            j2 = 2*(j-j_subtract)
            do i = i_prolong_start, i_prolong_end
                i2 = 2*(i-i_subtract)

                value_1 = tree%blocks(parent_idx)%multi_grid(next_mg_lvl)%p(i-1,j-1,k)
                value_2 = tree%blocks(parent_idx)%multi_grid(next_mg_lvl)%p(i-1,j  ,k)
                value_3 = tree%blocks(parent_idx)%multi_grid(next_mg_lvl)%p(i  ,j-1,k)
                value_4 = tree%blocks(parent_idx)%multi_grid(next_mg_lvl)%p(i  ,j  ,k)

                value_sum = ( &
                    (9.0_prcs_var) * value_1 + &
                    (3.0_prcs_var) * value_2 + &
                    (3.0_prcs_var) * value_3 + &
                    (1.0_prcs_var) * value_4 ) / 16.0_prcs_var
                tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2-1, j2-1, k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2-1, j2-1, k) + value_sum
                if (value_sum>prolong_max) prolong_max = value_sum

                value_sum = ( &
                    (3.0_prcs_var) * value_1 + &
                    (9.0_prcs_var) * value_2 + &
                    (1.0_prcs_var) * value_3 + &
                    (3.0_prcs_var) * value_4 ) / 16.0_prcs_var
                tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2-1, j2, k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2-1, j2, k) + value_sum
                if (value_sum>prolong_max) prolong_max = value_sum

                value_sum = ( &
                    (3.0_prcs_var) * value_1 + &
                    (1.0_prcs_var) * value_2 + &
                    (9.0_prcs_var) * value_3 + &
                    (3.0_prcs_var) * value_4 ) / 16.0_prcs_var
                tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2, j2-1, k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2, j2-1, k) + value_sum
                if (value_sum>prolong_max) prolong_max = value_sum

                value_sum = (&
                    (1.0_prcs_var) * value_1 + &
                    (3.0_prcs_var) * value_2 + &
                    (3.0_prcs_var) * value_3 + &
                    (9.0_prcs_var) * value_4 ) / 16.0_prcs_var
                tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2, j2, k) = &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i2, j2, k) + value_sum
                if (value_sum>prolong_max) prolong_max = value_sum

            end do
        end do

    end if
    

end subroutine MG_prolongate

subroutine MG_V_Cycle(tree, start_mg_lvl, prolong_max, residual)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module

    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: start_mg_lvl
    real(kind=prcs_var), intent(inout) :: residual
    real(kind=prcs_var), intent(inout) :: prolong_max

    integer :: max_lvl, mg_lvl_no, next_mg_lvl
    integer :: tree_lvl_no, last_tree_lvl, next_last_tree_lvl
    integer :: idx, id
    integer :: iter_mu
    character(len=2) :: string

    max_lvl = tree%no_levels



    ! Moving towards coarsest grid
    do mg_lvl_no = start_mg_lvl, max_lvl-1

        last_tree_lvl  = max_lvl - mg_lvl_no + 1


        ! RELAX Linear Equation on level mu_1 times
        do iter_mu = 1, mu_1
            do tree_lvl_no = last_tree_lvl, 1, -1
                do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                    idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                    if (MG_SOLVER == WEIGHTED_JACOBI) then
                        call block_Weighted_Jacobi( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p, &
                            weight)
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                    else if (MG_SOLVER == GAUSS_SEIDEL) then
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                        call block_Point_Gauss_Seidel( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p &
                        )
                    end if
                end do
                call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
            end do

        end do


        ! Calculate RESIDUAL at level
        residual = 0.0_prcs_var
        do tree_lvl_no = last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                call calculate_residue_array( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%r, &
                    dx_sq_lvl(tree_lvl_no), &
                    residual, &
                    gamma_p)
                
            end do
        end do

        ! Check convergence of residual
        residual = sqrt(residual)
        if (mg_lvl_no==1 .and. residual<=epsilon_p) return


        ! RESTRICT residual(r) to rhs(f) of next mg level
        ! INITIALISE variable(p) for next mg level
        next_mg_lvl         = mg_lvl_no + 1
        next_last_tree_lvl  = max_lvl - next_mg_lvl + 1
        do tree_lvl_no = next_last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(next_mg_lvl)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(next_mg_lvl)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call MG_restrict( &
                    tree, &
                    idx, &
                    next_mg_lvl)

                ! Initialise error variable (p) in next mg lvl to 0.0
                tree%blocks(idx)%multi_grid(next_mg_lvl)%p = 0.0_prcs_var

            end do
        end do


    end do


    ! call print_screen_real(tree%blocks(1)%multi_grid(1)%p, 'p-', 1, M+2, 1, N+2)
    
    ! Solving the coarsest grid (in Phase-1)
    mg_lvl_no   = max_lvl
    tree_lvl_no = 1
    do iter_mu = 1, mu_0
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
            idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
            if (MG_SOLVER == WEIGHTED_JACOBI) then
                call block_Weighted_Jacobi( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    gamma_p, &
                    weight)
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
            else if (MG_SOLVER == GAUSS_SEIDEL) then
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                call block_Point_Gauss_Seidel( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    gamma_p &
                )
            end if
            
        end do
        call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)

    end do


    ! call print_screen_real(tree%blocks(1)%multi_grid(1)%p, 'p-', 1, M+2, 1, N+2)

    !  ------------------------------------------------------------------------------------
    ! Calculate RESIDUAL at Coarsest level
    residual = 0.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
        idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
        call calculate_residue_array_body( &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%r, &
            dx_sq_lvl(tree_lvl_no), &
            residual, &
            gamma_p, &
            tree%blocks(idx)%i_range)
        
    end do

    ! Check convergence of residual
    residual = sqrt(residual)
    if (mg_lvl_no==1 .and. residual<=epsilon_p) return
    !  ------------------------------------------------------------------------------------

    ! Moving towards finest grid        
    do mg_lvl_no = max_lvl-1, start_mg_lvl, -1

        last_tree_lvl  = max_lvl - mg_lvl_no + 1

        ! PROLONGATE variable(p) to previous mg level
        prolong_max = 0.0_prcs_var
        do tree_lvl_no = last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call MG_prolongate( &
                    tree, &
                    idx, &
                    prolong_max, &
                    mg_lvl_no)
            end do
            call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
        end do


        ! RELAX Linear Equation on level mu_2 times
        do iter_mu = 1, mu_2
            do tree_lvl_no = last_tree_lvl, 1, -1
                do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                    idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                    if (MG_SOLVER == WEIGHTED_JACOBI) then
                        call block_Weighted_Jacobi( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p, &
                            weight)
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                    else if (MG_SOLVER == GAUSS_SEIDEL) then
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                        call block_Point_Gauss_Seidel( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p &
                        )
                    end if
                end do
                call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
            end do
        end do

    end do

end subroutine MG_V_Cycle

subroutine MG_V_Cycle_body(tree, start_mg_lvl, prolong_max, residual)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module

    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: start_mg_lvl
    real(kind=prcs_var), intent(inout) :: residual
    real(kind=prcs_var), intent(inout) :: prolong_max

    integer :: max_lvl, mg_lvl_no, next_mg_lvl
    integer :: tree_lvl_no, last_tree_lvl, next_last_tree_lvl
    integer :: idx, id, i1, i2
    integer :: iter_mu
    character(len=2) :: string
    logical :: found

    max_lvl = tree%no_levels

    ! Moving towards coarsest grid
    do mg_lvl_no = start_mg_lvl, max_lvl-1

        last_tree_lvl  = max_lvl - mg_lvl_no + 1

        ! if (mg_lvl_no == start_mg_lvl) then
        !     ! RELAX Linear Equation on level mu_1 times
        !     do iter_mu = 1, 50

        !         ! Update the body ghost cells
        !         do tree_lvl_no = tree%no_levels, 1, -1
        !             do id = 1, tree%levels(tree_lvl_no)%idx_body(0)
        !                 idx = tree%levels(tree_lvl_no)%idx_body(id)

        !                 found = .false.
        !                 do i1 = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
        !                     i2 = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(i1)
        !                     if(idx==i2) then
        !                         found = .true.
        !                         exit
        !                     end if
        !                 end do

        !                 if (found) then 
        !                     call interpolate_IP_variable(&
        !                         tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
        !                         tree%blocks(idx)%IP_interpol_cells, &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
        !                         tree%blocks(idx)%no_IB_ghost_cells(nbody) &
        !                     )
        !                     call update_IB_ghost_cells_variable(&
        !                         tree%blocks(idx)%IP_pos_val, &
        !                         tree%blocks(idx)%IB_ghost_cell, &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
        !                         tree%blocks(idx)%no_IB_ghost_cells, &
        !                         PRESSURE, &
        !                         tree%blocks(idx)%dl_GC_IP & 
        !                     )
        !                 end if

        !             end do
        !         end do


        !         do tree_lvl_no = last_tree_lvl, 1, -1
        !             do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
        !                 idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)

        !                 #IF (MG_SOLVER == WEIGHTED_JACOBI)
        !                     call block_Weighted_Jacobi( &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
        !                         gamma_p, &
        !                         weight)
        !                     call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)

        !                 #ELIF (MG_SOLVER == GAUSS_SEIDEL)
        !                     call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
        !                     call block_Point_Gauss_Seidel_body( &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
        !                         tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
        !                         gamma_p, &
        !                         tree%blocks(idx)%i_range &
        !                     )
        !                 #ENDIF
        !             end do
        !             call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
        !         end do

        !     end do
        ! else
            ! RELAX Linear Equation on level mu_1 times
            do iter_mu = 1, mu_1

                ! Update the body ghost cells
                do tree_lvl_no = tree%no_levels, 1, -1
                    do id = 1, tree%levels(tree_lvl_no)%idx_body(0)
                        idx = tree%levels(tree_lvl_no)%idx_body(id)

                        found = .false.
                        do i1 = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                            i2 = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(i1)
                            if(idx==i2) then
                                found = .true.
                                exit
                            end if
                        end do

                        if (found) then 
                            call interpolate_IP_variable(&
                                tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                                tree%blocks(idx)%IP_interpol_cells, &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                                tree%blocks(idx)%no_IB_ghost_cells(nbody) &
                            )
                            call update_IB_ghost_cells_variable(&
                                tree%blocks(idx)%IP_pos_val, &
                                tree%blocks(idx)%IB_ghost_cell, &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                                tree%blocks(idx)%no_IB_ghost_cells, &
                                PRESSURE, &
                                tree%blocks(idx)%dl_GC_IP & 
                            )
                        end if

                    end do
                end do


                do tree_lvl_no = last_tree_lvl, 1, -1
                    do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                        idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)

                        if (MG_SOLVER == WEIGHTED_JACOBI) then
                            call block_Weighted_Jacobi( &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                                gamma_p, &
                                weight)
                            call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)

                        else if (MG_SOLVER == GAUSS_SEIDEL) then
                            call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                            call block_Point_Gauss_Seidel_body( &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                                tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                                gamma_p, &
                                tree%blocks(idx)%i_range &
                            )
                        end if
                    end do
                    call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
                end do

            end do
        ! end if


        ! Calculate RESIDUAL at level
        residual = 0.0_prcs_var

        ! Update the body ghost cells
        do tree_lvl_no = tree%no_levels, 1, -1
            do id = 1, tree%levels(tree_lvl_no)%idx_body(0)
                idx = tree%levels(tree_lvl_no)%idx_body(id)
                
                found = .false.
                do i1 = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                    i2 = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(i1)
                    if(idx==i2) then
                        found = .true.
                        exit
                    end if
                end do

                if (found) then
                    call interpolate_IP_variable(&
                        tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                        tree%blocks(idx)%IP_interpol_cells, &
                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                        tree%blocks(idx)%no_IB_ghost_cells(nbody) &
                    )
                    call update_IB_ghost_cells_variable(&
                        tree%blocks(idx)%IP_pos_val, &
                        tree%blocks(idx)%IB_ghost_cell, &
                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                        tree%blocks(idx)%no_IB_ghost_cells, &
                        PRESSURE, &
                        tree%blocks(idx)%dl_GC_IP & 
                    )
                end if
            end do
        end do

        do tree_lvl_no = last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                call calculate_residue_array_body( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%r, &
                    dx_sq_lvl(tree_lvl_no), &
                    residual, &
                    gamma_p, &
                    tree%blocks(idx)%i_range)
                
            end do
        end do

        ! Check convergence of residual
        residual = sqrt(residual)
        if (mg_lvl_no==1 .and. residual<=epsilon_p) return


        ! RESTRICT residual(r) to rhs(f) of next mg level
        ! INITIALISE variable(p) for next mg level
        next_mg_lvl         = mg_lvl_no + 1
        next_last_tree_lvl  = max_lvl - next_mg_lvl + 1
        do tree_lvl_no = next_last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(next_mg_lvl)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(next_mg_lvl)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call MG_restrict( &
                    tree, &
                    idx, &
                    next_mg_lvl)

                ! Initialise error variable (p) in next mg lvl to 0.0
                tree%blocks(idx)%multi_grid(next_mg_lvl)%p = 0.0_prcs_var

            end do
        end do


    end do
    
    ! Solving the coarsest grid (in Phase-1)
    mg_lvl_no   = max_lvl
    tree_lvl_no = 1

    call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)

    do iter_mu = 1, mu_0
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
            idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)

            if (MG_SOLVER == WEIGHTED_JACOBI) then
                call block_Weighted_Jacobi( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    gamma_p, &
                    weight)
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)

            else if (MG_SOLVER == GAUSS_SEIDEL) then
                call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                call block_Point_Gauss_Seidel_body( &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                    gamma_p, &
                    tree%blocks(idx)%i_range &
                )

                ! call print_screen_real(tree%blocks(1)%multi_grid(1)%p, 'AA', 1, M+2*ngl, 1, N+2*ngl)
            end if
            
        end do
        call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
        
        ! print *, 'lvl: ', max_lvl, 'iter: ', iter_mu
        ! ! call print_screen_int(tree%blocks(1)%i_range, 'AA', 0, M+2*ngl, 1, N+2*ngl)
        ! call print_screen_real(tree%blocks(1)%multi_grid(1)%p, 'AA', 1, M+2*ngl, 1, N+2*ngl)

    end do


    !  ------------------------------------------------------------------------------------
    ! Calculate RESIDUAL at level
    residual = 0.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
        idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
        call calculate_residue_array_body( &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
            tree%blocks(idx)%multi_grid(mg_lvl_no)%r, &
            dx_sq_lvl(tree_lvl_no), &
            residual, &
            gamma_p, &
            tree%blocks(idx)%i_range)
        
    end do

    ! Check convergence of residual
    residual = sqrt(residual)
    if (mg_lvl_no==1 .and. residual<=epsilon_p) return

    ! !  ------------------------------------------------------------------------------------

    ! Moving towards finest grid        
    do mg_lvl_no = max_lvl-1, start_mg_lvl, -1

        last_tree_lvl  = max_lvl - mg_lvl_no + 1

        ! PROLONGATE variable(p) to previous mg level
        prolong_max = 0.0_prcs_var
        do tree_lvl_no = last_tree_lvl, 1, -1
            do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)
                call MG_prolongate( &
                    tree, &
                    idx, &
                    prolong_max, &
                    mg_lvl_no)
            end do
            call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
        end do


        ! RELAX Linear Equation on level mu_2 times
        do iter_mu = 1, mu_2

            ! Update the body ghost cells
            do tree_lvl_no = tree%no_levels, 1, -1
                do id = 1, tree%levels(tree_lvl_no)%idx_body(0)
                    idx = tree%levels(tree_lvl_no)%idx_body(id)
                    found = .false.
                    do i1 = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                        i2 = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(i1)
                        if(idx==i2) then
                            found = .true.
                            exit
                        end if
                    end do

                    if (found) then
                        call interpolate_IP_variable(&
                            tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                            tree%blocks(idx)%IP_interpol_cells, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%no_IB_ghost_cells(nbody) &
                        )
                        call update_IB_ghost_cells_variable(&
                            tree%blocks(idx)%IP_pos_val, &
                            tree%blocks(idx)%IB_ghost_cell, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%no_IB_ghost_cells, &
                            PRESSURE, &
                            tree%blocks(idx)%dl_GC_IP & 
                        )
                    end if
                end do
            end do


            do tree_lvl_no = last_tree_lvl, 1, -1
                do id = 1, size(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf)
                    idx = tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf(id)

                    if (MG_SOLVER == WEIGHTED_JACOBI) then
                        call block_Weighted_Jacobi( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p, &
                            weight)
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)

                    else if (MG_SOLVER == GAUSS_SEIDEL) then
                        call update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)
                        call block_Point_Gauss_Seidel( &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%p, &
                            tree%blocks(idx)%multi_grid(mg_lvl_no)%f, &
                            gamma_p &
                        )
                    end if
                end do
                call apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)
            end do
        end do

    end do

end subroutine MG_V_Cycle_body