
    #include "definitions.h"

! --------------------------------------------------------------
!         Solve Modified Momentum Equations
! --------------------------------------------------------------

!> Point Gauss Seidel Methods
subroutine solve_x_modified_momentum_equation_GS(tree, iter_u)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_u
    integer :: lvl, id, idx, iter_count,flag=1

    ! -- Calculate RHS of x modified momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 1 : Solve Modified Momentum Equation'
            print *, '--------------------------------------------------'
            print *
            print *, '>> Calculate RHS of x modified momentum equation'
        end if
    #ENDIF
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            call calculate_RHS_momentum( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%Cpx, &
                alpha(lvl), &
                dx_lvl(lvl), &
                dy_lvl(lvl), &
                dx_sq_lvl(lvl), &
                dy_sq_lvl(lvl), &
                coeff_u(lvl), &
                tree%blocks(idx)%i_range &
                )
        end do
    end do

    ! -- x-momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve x modified momentum equation'
            print *, '     (POINT GAUSS SEIDEL METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'
        end if
    #ENDIF

    ! ! Update ghost cell of U_VELOCITY around body (# not required as it is updated after time step completed)
    ! do lvl = tree%no_levels, 1, -1
    !     do id = 1, tree%levels(lvl)%idx_body_leaf(0)
    !         idx = tree%levels(lvl)%idx_body_leaf(id)

    !         call interpolate_IP_variable(&
    !             tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
    !             tree%blocks(idx)%IP_interpol_cells, &
    !             tree%blocks(idx)%u, &
    !             tree%blocks(idx)%no_IB_ghost_cells(nbody), &
    !             U_VELOCITY &
    !         )

    !         call update_IB_ghost_cells_variable(&
    !             tree%blocks(idx)%IP_pos_val, &
    !             tree%blocks(idx)%IB_ghost_cell, &
    !             tree%blocks(idx)%u, &
    !             tree%blocks(idx)%no_IB_ghost_cells, &
    !             U_VELOCITY, &
    !             tree%blocks(idx)%dl_GC_IP & 
    !         )

    !     end do
    ! end do


    residual_u = 1.0_prcs_var
    iter = 0
    do while(residual_u>epsilon_u .and. iter<max_iter)

        residual_u = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_u_vel(tree, idx)
                call block_Point_Gauss_Seidel_body( &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%rhs, &
                    gamma_u(lvl), &
                    tree%blocks(idx)%i_range, &
                    flag &
                    )
                call calculate_residual_body( &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%rhs, &
                    gamma_u(lvl), &
                    residual_u, &
                    tree%blocks(idx)%i_range &
                    )
            end do
        end do

        ! Update the Ghost cells with updated cell centre values
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_u_vel(tree, idx)
            end do
        end do

        call apply_boundary_condition_u_vel(tree)

        ! Update ghost cell of U_VELOCITY around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    U_VELOCITY &
                )

                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    U_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do

        residual_u = sqrt(residual_u)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_u
            end if
        #ENDIF
    
    end do


    iter_u = iter

    if (residual_u>epsilon_u) print '(A, F8.6, A, F10.8)', '      X - Mod Momentum Eq DID NOT CONVERGE at time: ',&
                                time,', residual: ', residual_u

end subroutine solve_x_modified_momentum_equation_GS

subroutine solve_y_modified_momentum_equation_GS(tree, iter_v)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_v
    integer :: lvl, id, idx, iter_count


    ! -- Calculate RHS of y modified momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *
            print *, '>> Calculate RHS of y modified momentum equation'
        end if
    #ENDIF

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            call calculate_RHS_momentum( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%Cpy, &
                alpha(lvl), &
                dx_lvl(lvl), &
                dy_lvl(lvl), &
                dx_sq_lvl(lvl), &
                dy_sq_lvl(lvl), &
                coeff_u(lvl), &
                tree%blocks(idx)%i_range &
            )
        end do
    end do
    
    ! -- y-momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve y modified momentum equation'
            print *, '     (POINT GAUSS SEIDEL METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'
        end if
    #ENDIF

    ! ! Update ghost cell of V_VELOCITY around body  (# not required as it is updated after time step completed)
    ! do lvl = tree%no_levels, 1, -1
    !     do id = 1, tree%levels(lvl)%idx_body_leaf(0)
    !         idx = tree%levels(lvl)%idx_body_leaf(id)

    !         call interpolate_IP_variable(&
    !             tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
    !             tree%blocks(idx)%IP_interpol_cells, &
    !             tree%blocks(idx)%v, &
    !             tree%blocks(idx)%no_IB_ghost_cells(nbody), &
    !             V_VELOCITY &
    !         )

    !         call update_IB_ghost_cells_variable(&
    !             tree%blocks(idx)%IP_pos_val, &
    !             tree%blocks(idx)%IB_ghost_cell, &
    !             tree%blocks(idx)%v, &
    !             tree%blocks(idx)%no_IB_ghost_cells, &
    !             V_VELOCITY, &
    !             tree%blocks(idx)%dl_GC_IP & 
    !         )

    !     end do
    ! end do

    residual_v = 1.0_prcs_var
    iter = 0
    do while(residual_v>epsilon_u .and. iter<max_iter)

        residual_v = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_v_vel(tree, idx)
                call block_Point_Gauss_Seidel_body( &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%rhs, &
                    gamma_u(lvl), &
                    tree%blocks(idx)%i_range &
                    )
                call calculate_residual_body( &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%rhs, &
                    gamma_u(lvl), &
                    residual_v, &
                    tree%blocks(idx)%i_range &
                    )
            end do
        end do

        ! Update the Ghost cells with updated cell centre values
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_v_vel(tree, idx)
            end do
        end do

        call apply_boundary_condition_v_vel(tree)

        ! Update ghost cell of V_VELOCITY around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    V_VELOCITY &
                )

                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    V_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do

        residual_v = sqrt(residual_v)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_v
            end if
        #ENDIF

    end do
    iter_v = iter

    if (residual_v>epsilon_u) print '(A, F8.6, A, F10.8)', '      Y - Mod Momentum Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_v

end subroutine solve_y_modified_momentum_equation_GS

subroutine solve_pressure_poisson_equation_GS(tree, iter_p)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_p
    integer :: lvl, id, idx, iter_count
    real(kind=prcs_var) :: rhs_max_val
    integer             :: rhs_max_loc(NDIM), max_idx


    ! -- Calculate RHS of pressure poisson equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 2 : Solve Pressure Poisson Equation'
            print *, '--------------------------------------------------'
            print * 
            print *, '>> Calculate RHS of pressure Poisson equation'
        end if
    #ENDIF    

    rhs_sum = 0.0_prcs_var
    rhs_max_val = 0.0_prcs_var
    max_idx = 0

    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve pressure Poisson equation'
            print *, '    (POINT GAUSS SEIDEL METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'
        end if
    #ENDIF

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            ! tree%blocks(idx)%rhs = 0.0_prcs_var ! Setting the source term as 0 to check the Poisson equation solver

            call calculate_rhs_pressure( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                coeff_p(lvl), &
                rhs_sum, &
                tree%blocks(idx)%i_range &
                )


        end do
    end do

    residual_p = 1.0_prcs_var
    iter = 0

    do while(residual_p>epsilon_p .and. iter<max_iter)

        residual_p = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_pressure(tree, idx, 1)

                call block_Point_Gauss_Seidel_body( &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%rhs, &
                    gamma_p, &
                    tree%blocks(idx)%i_range &
                    )

                call calculate_residual_body( &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%rhs, &
                    gamma_p, &
                    residual_p, &
                    tree%blocks(idx)%i_range &
                    )
            end do
        end do

        ! Update the Ghost cells with updated cell centre values
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_pressure(tree, idx)
            end do
        end do
        
        call apply_boundary_condition_pressure(tree)

        ! Update ghost cell of pressure around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    PRESSURE &
                )
                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    PRESSURE, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do
        
        residual_p = sqrt(residual_p)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_p
            end if
        #ENDIF

    end do
    iter_p = iter

    if (residual_p>epsilon_p) print '(A, F8.6, A, F10.8)', '      Pressure Poisson Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_p


end subroutine solve_pressure_poisson_equation_GS

! subroutine solve_pressure_poisson_equation_GS_singularity(tree, iter_p)

!     use precision_module, only: prcs_var
!     use data_type_module
!     use problem_module
!     use body_module
!     implicit none

!     type(t_tree), intent(inout) :: tree
!     integer, intent(inout)      :: iter_p
!     integer :: lvl, id, idx, iter_count
!     real(kind=prcs_var) :: rhs_max_val
!     integer             :: rhs_max_loc(NDIM), max_idx
!     integer             :: singular_cell_no, side
!     integer             :: i, j, k


!     ! -- Calculate RHS of pressure poisson equation
!     #IFDEF VERBOSE
!         if (mod(time_step+1,STEP_SKIP)==0) then
!             print *
!             print *, '--------------------------------------------------'
!             print *, '  STEP 2 : Solve Pressure Poisson Equation'
!             print *, '--------------------------------------------------'
!             print * 
!             print *, '>> Calculate RHS of pressure Poisson equation'
!         end if
!     #ENDIF    
    

!     rhs_sum = 0.0_prcs_var
!     rhs_max_val = 0.0_prcs_var
!     max_idx = 0

!     do lvl = tree%no_levels, 1, -1
!         do id = 1, no_leaf_idx_lvl(lvl)
!             idx = tree%levels(lvl)%idx_leaf(id)

!             ! tree%blocks(idx)%rhs = 0.0_prcs_var ! Setting the source term as 0 to check the Poisson equation solver

!             call calculate_rhs_pressure( &
!                 tree%blocks(idx)%rhs, &
!                 tree%blocks(idx)%Ue, &
!                 tree%blocks(idx)%Uw, &
!                 tree%blocks(idx)%Vn, &
!                 tree%blocks(idx)%Vs, &
!                 coeff_p(lvl), &
!                 rhs_sum, &
!                 tree%blocks(idx)%i_range &
!                 )

!             if (maxval(abs(tree%blocks(idx)%rhs)) > rhs_max_val) then
!                 rhs_max_val = maxval(abs(tree%blocks(idx)%rhs))
!                 rhs_max_loc = maxloc(abs(tree%blocks(idx)%rhs))
!                 max_idx = idx
!             end if

!         end do
!     end do

!     #IFDEF VERBOSE
!         if (mod(time_step+1,STEP_SKIP)==0) then
!             print *, '>> Solve pressure Poisson equation'
!             print *, '    (POINT GAUSS SEIDEL METHOD)'
!             print *
!             print '(T5, A, T20, A)', 'Iteration', 'Residual'
!             print *, '---------------------------------------'
!         end if
!     #ENDIF

!     residual_p = 1.0_prcs_var
!     iter = 0


!     ! print *, 'The sum of rhs of PPE: ', rhs_sum
!     ! print *, 'Max_val: ', rhs_max_val
!     ! print *, 'max_idx: ', max_idx
!     ! print *, 'Max_loc: ', rhs_max_loc

!     ! call print_screen_real(tree%blocks(max_idx)%rhs, 'rH', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%u, 'uc', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%Ue, 'Ue', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%Uw, 'Uw', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%Vn, 'Vn', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%Vs, 'Vs', 1, M+2, 1, N+2)
!     ! call print_screen_real(tree%blocks(max_idx)%p, 'Pr', 1, M+2, 1, N+2)

!     do while(residual_p>epsilon_p .and. iter<max_iter)

!         ! Update ghost cell of pressure around body
!         do lvl = tree%no_levels, 1, -1
!             do id = 1, tree%levels(lvl)%idx_body_leaf(0)
!                 idx = tree%levels(lvl)%idx_body_leaf(id)

!                 ! Update ghost cell of direction around body
!                 call interpolate_IP_variable(&
!                     tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
!                     tree%blocks(idx)%IP_interpol_cells, &
!                     tree%blocks(idx)%p, &
!                     tree%blocks(idx)%no_IB_ghost_cells(nbody), &
!                     PRESSURE &
!                 )
!                 call update_IB_ghost_cells_variable(&
!                     tree%blocks(idx)%IP_pos_val, &
!                     tree%blocks(idx)%IB_ghost_cell, &
!                     tree%blocks(idx)%p, &
!                     tree%blocks(idx)%no_IB_ghost_cells, &
!                     PRESSURE, &
!                     tree%blocks(idx)%dl_GC_IP & 
!                 )

!             end do
!         end do

!         ! Updating only the cells with common ghost cells
!         do singular_cell_no = 1, no_singularities

!             idx = body_singularity_cells(S_BLOCK, singular_cell_no)
!             i = body_singularity_cells(S_CELL_I, singular_cell_no)
!             j = body_singularity_cells(S_CELL_J, singular_cell_no)
!             k = body_singularity_cells(S_CELL_K, singular_cell_no)
!             side = body_singularity_cells(S_SIDE, singular_cell_no)

!             ! print *, 'side: ', side

!             #IFDEF SQ_CELL

!                 if (side == LEFT) then
!                     tree%blocks(idx)%p(i,j,k) =  (tree%blocks(idx)%rhs(i,j,k) - &
!                             ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) &
!                             +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i  ,j,k))) &
!                             )/gamma_p

!                 else if (side == RIGHT) then
!                     tree%blocks(idx)%p(i,j,k) =  (tree%blocks(idx)%rhs(i,j,k) - &
!                             ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) &
!                             +(tree%blocks(idx)%p(i  ,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k))) &
!                             )/gamma_p

!                 else if (side == BOTTOM) then
!                     tree%blocks(idx)%p(i,j,k) =  (tree%blocks(idx)%rhs(i,j,k) - &
!                             ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j  ,k)) &
!                             +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k))) &
!                             )/gamma_p

!                 else if (side == TOP) then
!                     tree%blocks(idx)%p(i,j,k) =  (tree%blocks(idx)%rhs(i,j,k) - &
!                             ((tree%blocks(idx)%p(i,j  ,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) &
!                             +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k))) &
!                             )/gamma_p

!                 end if

!             #ELSE 
!             ! Include RHS and modify as above
!                 if (side == LEFT) then
!                     tree%blocks(idx)%p(i,j,k) =  &
!                         ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) * beta_sq &
!                         +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i,j,k)) &
!                             +tree%blocks(idx)%p(i,j,k)*gamma_p)
                        
!                 else if (side == RIGHT) then
!                     tree%blocks(idx)%p(i,j,k) =  &
!                         ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) * beta_sq &
!                         +(tree%blocks(idx)%p(i,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k)) &
!                             +tree%blocks(idx)%p(i,j,k)*gamma_p)

!                 else if (side == BOTTOM) then
!                     tree%blocks(idx)%p(i,j,k) =  &
!                         ((tree%blocks(idx)%p(i,j-1,k) + &
!                             tree%blocks(idx)%p(i,j,k)) * beta_sq &
!                         +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k)) &
!                             +tree%blocks(idx)%p(i,j,k)*gamma_p)

!                 else if (side == TOP) then
!                     tree%blocks(idx)%p(i,j,k) =  &
!                         ((tree%blocks(idx)%p(i,j,k) + &
!                             tree%blocks(idx)%p(i,j+1,k)) * beta_sq &
!                         +(tree%blocks(idx)%p(i-1,j,k) + &
!                             tree%blocks(idx)%p(i+1,j,k)) &
!                             +tree%blocks(idx)%p(i,j,k)*gamma_p)

!                 end if

!             #ENDIF

!         end do


!         residual_p = 0.0_prcs_var
!         do lvl = tree%no_levels, 1, -1
!             do id = 1, no_leaf_idx_lvl(lvl)
!                 idx = tree%levels(lvl)%idx_leaf(id)

!                 call update_ghost_cell_pressure(tree, idx, 1)

!                 call block_Point_Gauss_Seidel_body( &
!                     tree%blocks(idx)%p, &
!                     tree%blocks(idx)%rhs, &
!                     gamma_p, &
!                     tree%blocks(idx)%i_range_singularity &
!                     )

!                 call calculate_residual_body( &
!                     tree%blocks(idx)%p, &
!                     tree%blocks(idx)%rhs, &
!                     gamma_p, &
!                     residual_p, &
!                     tree%blocks(idx)%i_range_singularity &
!                     )
!             end do
!         end do

!         ! Update the Ghost cells with updated cell centre values
!         do lvl = tree%no_levels, 1, -1
!             do id = 1, no_leaf_idx_lvl(lvl)
!                 idx = tree%levels(lvl)%idx_leaf(id)
!                 call update_ghost_cell_pressure(tree, idx)
!             end do
!         end do

        
!         call apply_boundary_condition_pressure(tree)


!         ! call print_screen_real(tree%blocks(max_idx)%p, 'Pr', 1, M+2, 1, N+2)

        
!         residual_p = sqrt(residual_p)

!         iter = iter + 1

!         #IFDEF VERBOSE
!             if (mod(time_step+1,STEP_SKIP)==0) then
!                 print '(T5, I6, T17, E14.4)', iter, residual_p
!             end if
!         #ENDIF

!     end do
!     iter_p = iter

!     if (residual_p>epsilon_p) print '(A, F8.6, A, F10.8)', '      Pressure Poisson Eq DID NOT CONVERGE at time: ',&
!                         time,', residual: ', residual_p


! end subroutine solve_pressure_poisson_equation_GS_singularity

subroutine solve_energy_equation_GS(tree, iter_T)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_T
    integer :: lvl, id, idx, iter_count


    ! -- Calculate RHS of energy equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 4 : Solve Energy Equation'
            print *, '--------------------------------------------------'
            print *
            print *, '>> Calculate RHS of energy equation'
        end if
    #ENDIF
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            call calculate_RHS_temperature( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%T,  &
                dx_lvl(lvl), &
                dy_lvl(lvl), &
                tree%blocks(idx)%CpT, &
                coeff_T(lvl) &
                )
        end do
    end do

    ! -- energy equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve energy equation'
            print *, '   (POINT GAUSS SEIDEL METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'
        end if
    #ENDIF

    residual_T = 1.0_prcs_var
    iter = 0
    do while(residual_T>epsilon_T .and. iter<max_iter)

        residual_T = 0.0_prcs_var        
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_temperature(tree, idx)
                call block_Point_Gauss_Seidel( &
                    tree%blocks(idx)%T, &
                    tree%blocks(idx)%rhs, &
                    gamma_T(lvl) &
                    )
                call calculate_residual( &
                    tree%blocks(idx)%T, &
                    tree%blocks(idx)%rhs, &
                    gamma_T(lvl), &
                    residual_T &
                    )
            end do
        end do

        ! Update the Ghost cells with updated cell centre values
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_temperature(tree, idx)
            end do
        end do

        call apply_boundary_condition_temperature(tree)

        residual_T = sqrt(residual_T)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_T
            end if
        #ENDIF


    end do
    iter_T = iter

    if (residual_T>epsilon_T) print '(A, F8.6, A, F10.8)', '      Energy Eq DID NOT CONVERGE at time: ',&
                                time,', residual: ', residual_T

end subroutine solve_energy_equation_GS


!> Conjugate Gradient Methods

subroutine solve_x_modified_momentum_equation_CG(tree, iter_u)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_u
    integer :: lvl, id, idx, iter_count
    integer :: i, j, k, i1, l
    real(kind=prcs_var) :: dAd, rTr, rTr_new
    real(kind=prcs_var) :: cg_alpha, cg_beta


    ! -- Calculate RHS of x modified momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 1 : Solve Modified Momentum Equation'
            print *, '--------------------------------------------------'
            print *
            print *, '>> Calculate RHS of x modified momentum equation'
        end if
    #ENDIF
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            call calculate_RHS_momentum( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%Cpx, &
                alpha(lvl), &
                dx_lvl(lvl), &
                dy_lvl(lvl), &
                dx_sq_lvl(lvl), &
                dy_sq_lvl(lvl), &
                coeff_u(lvl), &
                tree%blocks(idx)%i_range &
                )
        end do
    end do

    rTr      = 0.0_prcs_var
    cg_alpha = 0.0_prcs_var
    cg_beta  = 0.0_prcs_var


    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve x modified momentum equation'
            print *, '     (CONJUGATE GRADIENT METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'

        end if
    #ENDIF

    ! ! Update ghost cell of U_VELOCITY around body (# not required as value already updated after time-step)
    ! do lvl = tree%no_levels, 1, -1
    !     do id = 1, tree%levels(lvl)%idx_body_leaf(0)
    !         idx = tree%levels(lvl)%idx_body_leaf(id)

    !         call interpolate_IP_variable(&
    !             tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
    !             tree%blocks(idx)%IP_interpol_cells, &
    !             tree%blocks(idx)%u, &
    !             tree%blocks(idx)%no_IB_ghost_cells(nbody), &
    !             U_VELOCITY &
    !         )

    !         call update_IB_ghost_cells_variable(&
    !             tree%blocks(idx)%IP_pos_val, &
    !             tree%blocks(idx)%IB_ghost_cell, &
    !             tree%blocks(idx)%u, &
    !             tree%blocks(idx)%no_IB_ghost_cells, &
    !             U_VELOCITY, &
    !             tree%blocks(idx)%dl_GC_IP & 
    !         )

    !     end do
    ! end do

    ! Initialise r, d, and r.r
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        ! Calculate r(0)
                        #IFDEF SQ_CELL
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%u(i,j-1,k) + tree%blocks(idx)%u(i,j+1,k)) &
                                 +(tree%blocks(idx)%u(i-1,j,k) + tree%blocks(idx)%u(i+1,j,k)) &
                                  +tree%blocks(idx)%u(i,j,k)*gamma_u(lvl) )
                        #ELSE
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%u(i,j-1,k) + tree%blocks(idx)%u(i,j+1,k)) * beta_sq &
                                 +(tree%blocks(idx)%u(i-1,j,k) + tree%blocks(idx)%u(i+1,j,k)) &
                                  +tree%blocks(idx)%u(i,j,k)*gamma_u(lvl) )
                        #ENDIF

                        ! Calculate d(0)
                        tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k)

                        ! Calculate r.r
                        rTr = rTr + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                    end do

                end do
            end do
        end do
    end do

    ! Calculate dAd iniially to check if A is SPD matrix
    dAd = 0.0_prcs_var
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            call update_ghost_cell_direction(tree, idx)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        #IFDEF SQ_CELL
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                    +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                    + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                        #ELSE
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                    +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                    + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                        #ENDIF

                        dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                    end do

                end do
            end do

        end do
    end do
    
    if (dAd == 0.0_prcs_var) return

    residual_u = 1.0_prcs_var
    iter = 0
    do while(residual_u>epsilon_u .and. iter<max_iter)

        residual_u = 0.0_prcs_var

        call apply_boundary_condition_direction(tree)

        ! Update ghost cell of direction around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    U_VELOCITY &
                )

                call update_IB_ghost_cells_direction(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    U_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do

        ! Calculate dAd
        dAd = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_direction(tree, idx)

                k = 1
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            #IFDEF SQ_CELL
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                            #ELSE
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                            #ENDIF

                            dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                        end do

                    end do
                end do

            end do
        end do

        ! Calculate alpha
        cg_alpha = rTr/dAd

        ! Update U_VELOCITY, residual and rTr_new
        rTr_new = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                do  j = j_start, j_end
                    ! do i = i_start, i_end

                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%u(i,j,k) = tree%blocks(idx)%u(i,j,k) &
                                                        + cg_alpha*tree%blocks(idx)%d(i,j,k)

                            tree%blocks(idx)%r(i,j,k) = tree%blocks(idx)%r(i,j,k) - cg_alpha*tree%blocks(idx)%Ad(i,j,k)

                            rTr_new = rTr_new + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        ! Update ghost cell of U_VELOCITY between blocks
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_u_vel(tree, idx)

            end do
        end do

        ! Update ghost cell of U_VELOCITY around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    U_VELOCITY &
                )

                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    U_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )


            end do
        end do


        ! Applying Boundary Condition
        call apply_boundary_condition_u_vel(tree)


        ! Calculate beta
        cg_beta = rTr_new / rTr

        rTr = rTr_new

        ! Find new direction
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k) + cg_beta*tree%blocks(idx)%d(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        
        residual_u = sqrt(rTr_new)

        iter = iter + 1

        ! print *, 'Value after end of iteration', tree%blocks(1)%u(5,8,1)
        ! call print_screen_real(tree%blocks(1)%u, 'Uv', 1, M+2, 1, N+2)


        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_u
            end if
        #ENDIF


    end do

    ! print *, '-------------- MAX VAL and LOC --------------'
    ! print *,  maxval((tree%blocks(idx)%u(i_start:i_end, j_start:j_end, 1:1)))
    ! print *,  maxloc((tree%blocks(idx)%u(i_start:i_end, j_start:j_end, 1:1)))

    ! print *, '-------------- MIN VAL and LOC --------------'
    ! print *,  minval((tree%blocks(idx)%u(i_start:i_end, j_start:j_end, 1:1)))
    ! print *,  minloc((tree%blocks(idx)%u(i_start:i_end, j_start:j_end, 1:1)))

    iter_u = iter

    if (residual_u>epsilon_u) print '(A, F8.6, A, F16.6)', '      X - Mod Momentum Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_u


end subroutine solve_x_modified_momentum_equation_CG

subroutine solve_y_modified_momentum_equation_CG(tree, iter_v)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_v
    integer :: lvl, id, idx, iter_count
    integer :: i, j, k, i1, l
    real(kind=prcs_var) :: dAd, rTr, rTr_new
    real(kind=prcs_var) :: cg_alpha, cg_beta


    ! -- Calculate RHS of y modified momentum equation
    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '>> Calculate RHS of y modified momentum equation'
        end if
    #ENDIF
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            call calculate_RHS_momentum( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%Cpy, &
                alpha(lvl), &
                dx_lvl(lvl), &
                dy_lvl(lvl), &
                dx_sq_lvl(lvl), &
                dy_sq_lvl(lvl), &
                coeff_u(lvl), &
                tree%blocks(idx)%i_range &
                )
        end do
    end do

    rTr      = 0.0_prcs_var
    cg_alpha = 0.0_prcs_var
    cg_beta  = 0.0_prcs_var


    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve y modified momentum equation'
            print *, '     (CONJUGATE GRADIENT METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'

        end if
    #ENDIF

    ! ! Update ghost cell of V_VELOCITY around body (# not required, already updated after time-step)
    ! do lvl = tree%no_levels, 1, -1
    !     do id = 1, tree%levels(lvl)%idx_body_leaf(0)
    !         idx = tree%levels(lvl)%idx_body_leaf(id)

    !         call interpolate_IP_variable(&
    !             tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
    !             tree%blocks(idx)%IP_interpol_cells, &
    !             tree%blocks(idx)%v, &
    !             tree%blocks(idx)%no_IB_ghost_cells(nbody), &
    !             V_VELOCITY &
    !         )

    !         call update_IB_ghost_cells_variable(&
    !             tree%blocks(idx)%IP_pos_val, &
    !             tree%blocks(idx)%IB_ghost_cell, &
    !             tree%blocks(idx)%v, &
    !             tree%blocks(idx)%no_IB_ghost_cells, &
    !             V_VELOCITY, &
    !             tree%blocks(idx)%dl_GC_IP & 
    !         )

    !     end do
    ! end do

    ! Initialise r, d, and r.r
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        ! Calculate r(0)
                        #IFDEF SQ_CELL
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%v(i,j-1,k) + tree%blocks(idx)%v(i,j+1,k)) &
                                 +(tree%blocks(idx)%v(i-1,j,k) + tree%blocks(idx)%v(i+1,j,k)) &
                                  +tree%blocks(idx)%v(i,j,k)*gamma_u(lvl) )
                        #ELSE
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%v(i,j-1,k) + tree%blocks(idx)%v(i,j+1,k)) * beta_sq &
                                 +(tree%blocks(idx)%v(i-1,j,k) + tree%blocks(idx)%v(i+1,j,k)) &
                                  +tree%blocks(idx)%v(i,j,k)*gamma_u(lvl) )
                        #ENDIF

                        ! Calculate d(0)
                        tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k)

                        ! Calculate r.r
                        rTr = rTr + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                    end do

                end do
            end do
        end do
    end do

    ! Calculate dAd initially to check if A is SPD matrix
    dAd = 0.0_prcs_var
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            call update_ghost_cell_direction(tree, idx)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        #IFDEF SQ_CELL
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                    +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                    + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                        #ELSE
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                    +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                    + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                        #ENDIF

                        dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                    end do

                end do
            end do

        end do
    end do

    if (dAd == 0.0_prcs_var) return

    residual_v = 1.0_prcs_var
    iter = 0
    do while(residual_v>epsilon_u .and. iter<max_iter)

        residual_v = 0.0_prcs_var

        call apply_boundary_condition_direction(tree)

        ! Update ghost cell of direction around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    V_VELOCITY &
                )

                call update_IB_ghost_cells_direction(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    V_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do

        ! Calculate dAd
        dAd = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_direction(tree, idx)

                k = 1
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            #IFDEF SQ_CELL
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                            #ELSE
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_u(lvl))
                            #ENDIF

                            dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                        end do

                    end do
                end do

            end do
        end do

        ! Calculate alpha
        cg_alpha = rTr/dAd

        ! Update V_VELOCITY, residual and rTr_new
        rTr_new = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                do  j = j_start, j_end
                    ! do i = i_start, i_end

                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%v(i,j,k) = tree%blocks(idx)%v(i,j,k) &
                                                        + cg_alpha*tree%blocks(idx)%d(i,j,k)

                            tree%blocks(idx)%r(i,j,k) = tree%blocks(idx)%r(i,j,k) - cg_alpha*tree%blocks(idx)%Ad(i,j,k)

                            rTr_new = rTr_new + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        ! Update ghost cell of V_VELOCITY between blocks
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_v_vel(tree, idx)

            end do
        end do

        ! Update ghost cell of V_VELOCITY around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    V_VELOCITY &
                )

                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    V_VELOCITY, &
                    tree%blocks(idx)%dl_GC_IP & 
                )


            end do
        end do

        ! Applying Boundary Condition
        call apply_boundary_condition_v_vel(tree)

        ! Calculate beta
        cg_beta = rTr_new / rTr

        rTr = rTr_new

        ! Find new direction
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k) + cg_beta*tree%blocks(idx)%d(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        
        residual_v = sqrt(rTr_new)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_v
            end if
        #ENDIF


    end do
    iter_v = iter

    if (residual_v>epsilon_u) print '(A, F8.6, A, F16.6)', '      X - Mod Momentum Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_v


end subroutine solve_y_modified_momentum_equation_CG

subroutine solve_pressure_poisson_equation_CG(tree, iter_p)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(inout)      :: iter_p
    integer :: lvl, id, idx, iter_count
    real(kind=prcs_var) :: rhs_max_val
    integer             :: rhs_max_loc(NDIM), max_idx
    integer :: i, j, k, i1, l
    real(kind=prcs_var) :: dAd, rTr, rTr_new
    real(kind=prcs_var) :: cg_alpha, cg_beta


    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 2 : Solve Pressure Poisson Equation'
            print *, '--------------------------------------------------'
            print *
            print *, '>> Calculate RHS of pressure Poisson equation'
        end if
    #ENDIF

    ! -- Calculate RHS of pressure poisson equation and Initialise CG variables
    rhs_sum = 0.0_prcs_var
    ! rhs_max_val = 0.0_prcs_var
    ! max_idx = 0

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            ! tree%blocks(idx)%rhs = 0.0_prcs_var ! Setting the source term as 0 to check the Poisson equation solver

            call calculate_rhs_pressure( &
                tree%blocks(idx)%rhs, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                coeff_p(lvl), &
                rhs_sum, &
                tree%blocks(idx)%i_range &
                )

        end do
    end do

    ! print *, 'RHS Sum: ', rhs_sum


    rTr      = 0.0_prcs_var
    cg_alpha = 0.0_prcs_var
    cg_beta  = 0.0_prcs_var



    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve pressure Poisson equation'
            print *, '    (CONJUGATE GRADIENT METHOD)'
            print *
            print '(T5, A, T20, A)', 'Iteration', 'Residual'
            print *, '---------------------------------------'

        end if
    #ENDIF

    ! ! Update ghost cell of pressure around body (# not required, updated already after time-step)
    ! do lvl = tree%no_levels, 1, -1
    !     do id = 1, tree%levels(lvl)%idx_body_leaf(0)
    !         idx = tree%levels(lvl)%idx_body_leaf(id)

    !         call interpolate_IP_variable(&
    !             tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
    !             tree%blocks(idx)%IP_interpol_cells, &
    !             tree%blocks(idx)%p, &
    !             tree%blocks(idx)%no_IB_ghost_cells(nbody), &
    !             PRESSURE &
    !         )

    !         call update_IB_ghost_cells_variable(&
    !             tree%blocks(idx)%IP_pos_val, &
    !             tree%blocks(idx)%IB_ghost_cell, &
    !             tree%blocks(idx)%p, &
    !             tree%blocks(idx)%no_IB_ghost_cells, &
    !             PRESSURE, &
    !             tree%blocks(idx)%dl_GC_IP & 
    !         )

    !     end do
    ! end do

    ! Initialise r, d, and r.r
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            ! call update_ghost_cell_pressure(tree, idx, 1)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        ! Calculate r(0)
                        #IFDEF SQ_CELL
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%p(i,j-1,k) + tree%blocks(idx)%p(i,j+1,k)) &
                                 +(tree%blocks(idx)%p(i-1,j,k) + tree%blocks(idx)%p(i+1,j,k)) &
                                  +tree%blocks(idx)%p(i,j,k)*gamma_p )
                        #ELSE
                            tree%blocks(idx)%r(i,j,k) =  &
                                tree%blocks(idx)%rhs(i,j,k) - &
                                ( (tree%blocks(idx)%p(i,j-1,k) + tree%blocks(idx)%p(i,j+1,k)) * beta_sq &
                                 +(tree%blocks(idx)%p(i-1,j,k) + tree%blocks(idx)%p(i+1,j,k)) &
                                  +tree%blocks(idx)%p(i,j,k)*gamma_p )
                        #ENDIF

                        ! Calculate d(0)
                        tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k)

                        ! Calculate r.r
                        rTr = rTr + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                        ! ! Summing up the values of residue (or direction)
                        ! d_mean = d_mean + tree%blocks(idx)%d(i,j,k)

                    end do

                end do
            end do
        end do
    end do
    
    ! if (pressure_BC_neumann) then
    !     call reset_direction_non_SPD(tree)
    !     call reset_residual_non_SPD(tree)
    ! end if


    ! Calculate Ad
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            call update_ghost_cell_direction(tree, idx)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        #IFDEF SQ_CELL
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                   +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                   + tree%blocks(idx)%d(i,j,k)*gamma_p)
                        #ELSE
                            tree%blocks(idx)%Ad(i,j,k) =  ( &
                                    (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                   +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                   + tree%blocks(idx)%d(i,j,k)*gamma_p)
                        #ENDIF

                        ! dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                    end do

                end do
            end do

        end do
    end do

    ! if (pressure_BC_neumann) then
    !     call reset_Ad_non_SPD(tree)
    ! end if

    ! Calcaulte dAd
    dAd = 0.0_prcs_var
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                ! do i = i_start, i_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                        dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)

                    end do

                end do
            end do

        end do
    end do

    if (dAd == 0.0_prcs_var) return

    residual_p = 1.0_prcs_var
    iter = 0
    do while(residual_p>epsilon_p .and. iter<max_iter)

        residual_p = 0.0_prcs_var

        call apply_boundary_condition_direction(tree)

        ! Update ghost cell of direction around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    PRESSURE &
                )

                call update_IB_ghost_cells_direction(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%d, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    PRESSURE, &
                    tree%blocks(idx)%dl_GC_IP & 
                )

            end do
        end do

        ! Calculate Ad
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
    
                call update_ghost_cell_direction(tree, idx)
    
                k = 1
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
    
                            #IFDEF SQ_CELL
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_p)
                            #ELSE
                                tree%blocks(idx)%Ad(i,j,k) =  ( &
                                        (tree%blocks(idx)%d(i,j-1,k) + tree%blocks(idx)%d(i,j+1,k)) *beta_sq &
                                       +(tree%blocks(idx)%d(i-1,j,k) + tree%blocks(idx)%d(i+1,j,k)) &
                                       + tree%blocks(idx)%d(i,j,k)*gamma_p)
                            #ENDIF
    
                            ! dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)
    
                        end do
    
                    end do
                end do
    
            end do
        end do
    
        ! if (pressure_BC_neumann) then
        !     call reset_Ad_non_SPD(tree)
        ! end if
    
        ! Calcaulte dAd
        dAd = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
    
                k = 1
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
    
                            dAd = dAd + tree%blocks(idx)%d(i,j,k) * tree%blocks(idx)%Ad(i,j,k)
    
                        end do
    
                    end do
                end do
    
            end do
        end do

        ! Calculate alpha
        cg_alpha = rTr/dAd

        ! Update pressure, residual and rTr_new
        rTr_new = 0.0_prcs_var
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                do  j = j_start, j_end
                    ! do i = i_start, i_end

                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%p(i,j,k) = tree%blocks(idx)%p(i,j,k) &
                                                        + cg_alpha*tree%blocks(idx)%d(i,j,k)

                            tree%blocks(idx)%r(i,j,k) = tree%blocks(idx)%r(i,j,k) - cg_alpha*tree%blocks(idx)%Ad(i,j,k)

                            rTr_new = rTr_new + tree%blocks(idx)%r(i,j,k)*tree%blocks(idx)%r(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        ! Update ghost cell of pressure between blocks
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                call update_ghost_cell_pressure(tree, idx)

            end do
        end do

        ! Update ghost cell of pressure around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
                idx = tree%levels(lvl)%idx_body_leaf(id)

                call interpolate_IP_variable(&
                    tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                    tree%blocks(idx)%IP_interpol_cells, &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                    PRESSURE &
                )

                call update_IB_ghost_cells_variable(&
                    tree%blocks(idx)%IP_pos_val, &
                    tree%blocks(idx)%IB_ghost_cell, &
                    tree%blocks(idx)%p, &
                    tree%blocks(idx)%no_IB_ghost_cells, &
                    PRESSURE, &
                    tree%blocks(idx)%dl_GC_IP & 
                )


            end do
        end do


        ! Applying Boundary Condition
        call apply_boundary_condition_pressure(tree)

        ! Calculate beta
        cg_beta = rTr_new / rTr

        rTr = rTr_new

        ! Find new direction
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                do j = j_start, j_end
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%r(i,j,k) + cg_beta*tree%blocks(idx)%d(i,j,k)

                        end do
                    end do
                end do
            end do
        end do

        
        residual_p = sqrt(rTr_new)

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4)', iter, residual_p
            end if
        #ENDIF


    end do
    iter_p = iter


    if (residual_p>epsilon_p) print '(A, F8.6, A, F16.6)', '      Pressure Poisson Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_p


end subroutine solve_pressure_poisson_equation_CG


!> Multi Grid Methods
subroutine solve_pressure_poisson_equation_MG(tree, iter_p)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout)     :: tree
    integer, intent(inout)          :: iter_p
    integer :: lvl, id, idx, iter_count
    real(kind=prcs_var) :: prolong_max
    integer :: i, j, k

    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 2 : Solve Pressure Poisson Equation'
            print *, '--------------------------------------------------'
            print *
            print *, '>> Calculate RHS of pressure Poisson equation'
        end if
    #ENDIF
    
    ! -- Calculate RHS of pressure poisson equation
    rhs_sum = 0.0_prcs_var
    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            ! tree%blocks(idx)%multi_grid(1)%f = 0.0_prcs_var ! Setting the source term as 0 to check the Poisson equation solver

            call calculate_rhs_pressure( &
                tree%blocks(idx)%multi_grid(1)%f, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                coeff_p(lvl), &
                rhs_sum, &
                tree%blocks(idx)%i_range &
                )
        end do
    end do


    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *, '>> Solve pressure Poisson equation'
            print *, '      (MULTI GRID METHOD)'
            print *
            print '(T5, A, T20, A, T39, A)', 'Iteration', 'Residual', 'Max_Error'
            print *, '-------------------------------------------------------'
        end if
    #ENDIF


    residual_p = 1.0_prcs_var
    iter = 0
    do while(residual_p>epsilon_p .and. iter<max_iter)

        call MG_V_Cycle( &
            tree, & ! Tree
            1, &    ! start_mg_lv
            prolong_max, &
            residual_p) ! weight

        ! call MG_V_Cycle_body( &
        !     tree, & ! Tree
        !     1, &    ! start_mg_lv
        !     prolong_max, &
        !     residual_p) ! weight

        iter = iter + 1

        #IFDEF VERBOSE
            if (mod(time_step+1,STEP_SKIP)==0) then
                print '(T5, I6, T17, E14.4, T35, E14.4)', iter, residual_p, prolong_max
            end if
        #ENDIF

    end do
    iter_p = iter

    if (residual_p>epsilon_p) print '(A, F8.6, A, F16.6)', '      Pressure Poisson Eq DID NOT CONVERGE at time: ',&
                        time,', residual: ', residual_p


end subroutine solve_pressure_poisson_equation_MG
