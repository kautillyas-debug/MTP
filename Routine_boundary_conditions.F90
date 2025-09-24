
! --------------------------------------------------------------
!>        Apply Boundary Condition
! --------------------------------------------------------------
subroutine boundary_ghost_cell_calculate(BC_condition_flag, value_wall, gradient_wall, &
                                         array_ghost, array_interior, dn, direction_flag, arr_dim)

    #include "definitions.h"
    use precision_module, only: prcs_var

    implicit none

    integer, intent(in)             :: BC_condition_flag, &             !> Boundary Condition, 1-Dirichlet, 2-Neumann, 3-Mixed
                                       arr_dim                       !> Size of the dummy arrays
    real(kind=prcs_var), intent(in) :: value_wall, gradient_wall, &     !> Boundary value of the wall
                                       dn, &                            !> dx or dy normal to boundary
                                       direction_flag, &                !> -1 if left or bottom, +1 if top or right
                                       array_interior(arr_dim)          !> array next to ghost array inside the domain
    real(kind=prcs_var), intent(inout)  :: array_ghost(arr_dim)          !> boundary (ghost) array which is to be updated

    if (BC_condition_flag == DIRICHLET) then       !> Dirichlet condition
        array_ghost(:) = 2.0_prcs_var*value_wall - array_interior(:)

    else if(BC_condition_flag == NEUMANN) then   !> Neumann condition
        array_ghost(:) = dn*gradient_wall*direction_flag + array_interior(:)
    end if

end subroutine boundary_ghost_cell_calculate

!> ------------------- VELOCITIES -------------------
subroutine apply_boundary_condition_u_vel(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION

    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%left_blocks)
        idx = tree%left_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- u - Velocity ---
        call boundary_ghost_cell_calculate(u_BC_left, u_left, du_dx_left, &
                                      tree%blocks(idx)%u(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%u(i_start    ,:,1), &      ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
        
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%right_blocks)
        idx = tree%right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- u - Velocity ---
        call boundary_ghost_cell_calculate(u_BC_right, u_right, du_dx_right, &
                                      tree%blocks(idx)%u(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%u(i_end    ,:,1), &    ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%bottom_blocks)
        idx = tree%bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- u - Velocity ---
        call boundary_ghost_cell_calculate(u_BC_bottom, u_bottom, du_dy_bottom, &
                                      tree%blocks(idx)%u(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%u(:,j_start    ,1), &      ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%top_blocks)
        idx = tree%top_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- u - Velocity ---
        call boundary_ghost_cell_calculate(u_BC_top, u_top, du_dy_top, &
                                      tree%blocks(idx)%u(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%u(:,j_end    ,1), &    ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do

end subroutine apply_boundary_condition_u_vel

subroutine apply_boundary_condition_v_vel(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION

    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%left_blocks)
        idx = tree%left_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- v - Velocity ---
        call boundary_ghost_cell_calculate(v_BC_left, v_left, dv_dx_left, &
                                      tree%blocks(idx)%v(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%v(i_start    ,:,1), &      ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
        
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%right_blocks)
        idx = tree%right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- v - Velocity ---
        call boundary_ghost_cell_calculate(v_BC_right, v_right, dv_dx_right, &
                                      tree%blocks(idx)%v(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%v(i_end    ,:,1), &    ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%bottom_blocks)
        idx = tree%bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- v - Velocity ---
        call boundary_ghost_cell_calculate(v_BC_bottom, v_bottom, dv_dy_bottom, &
                                      tree%blocks(idx)%v(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%v(:,j_start    ,1), &      ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%top_blocks)
        idx = tree%top_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- v - Velocity ---
        call boundary_ghost_cell_calculate(v_BC_top, v_top, dv_dy_top, &
                                      tree%blocks(idx)%v(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%v(:,j_end    ,1), &    ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do

end subroutine apply_boundary_condition_v_vel

!> ------------------- PRESSURE -------------------
subroutine apply_boundary_condition_pressure(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION

    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%left_blocks)
        idx = tree%left_blocks(id)
        lvl = tree%blocks(idx)%lvl

        !>  --- Pressure ---
        if (PPE_SOLVER==MULTI_GRID) then
            call boundary_ghost_cell_calculate(p_BC_left, p_left, dp_dx_left, &
                                            tree%blocks(idx)%multi_grid(1)%p(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                            tree%blocks(idx)%multi_grid(1)%p(i_start    ,:,1), &      ! Interior array
                                            dx_lvl(lvl), direction_flag, N+2)
        else
            call boundary_ghost_cell_calculate(p_BC_left, p_left, dp_dx_left, &
                                            tree%blocks(idx)%p(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                            tree%blocks(idx)%p(i_start    ,:,1), &      ! Interior array
                                            dx_lvl(lvl), direction_flag, N+2)

        end if
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%right_blocks)
        idx = tree%right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- Pressure ---
        if (PPE_SOLVER==MULTI_GRID) then
            call boundary_ghost_cell_calculate(p_BC_right, p_right, dp_dx_right, &
                                            tree%blocks(idx)%multi_grid(1)%p(i_end+ngl,:,1), &      ! Boundary (ghost) array
                                            tree%blocks(idx)%multi_grid(1)%p(i_end    ,:,1), &      ! Interior array
                                            dx_lvl(lvl), direction_flag, N+2)
        else
            call boundary_ghost_cell_calculate(p_BC_right, p_right, dp_dx_right, &
                                            tree%blocks(idx)%p(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                            tree%blocks(idx)%p(i_end    ,:,1), &    ! Interior array
                                            dx_lvl(lvl), direction_flag, N+2)
        end if
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%bottom_blocks)
        idx = tree%bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- Pressure ---
        if (PPE_SOLVER==MULTI_GRID) then
            call boundary_ghost_cell_calculate(p_BC_bottom, p_bottom, dp_dy_bottom, &
                                            tree%blocks(idx)%multi_grid(1)%p(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                            tree%blocks(idx)%multi_grid(1)%p(:,j_start    ,1), &      ! Interior array
                                            dy_lvl(lvl), direction_flag, M+2)
        else
            call boundary_ghost_cell_calculate(p_BC_bottom, p_bottom, dp_dy_bottom, &
                                            tree%blocks(idx)%p(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                            tree%blocks(idx)%p(:,j_start    ,1), &      ! Interior array
                                            dy_lvl(lvl), direction_flag, M+2)
        end if
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%top_blocks)
        idx = tree%top_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- Pressure ---
        if (PPE_SOLVER==MULTI_GRID) then
            call boundary_ghost_cell_calculate(p_BC_top, p_top, dp_dy_top, &
                                                tree%blocks(idx)%multi_grid(1)%p(:,j_end+ngl,1), &      ! Boundary (ghost) array
                                                tree%blocks(idx)%multi_grid(1)%p(:,j_end    ,1), &      ! Interior array
                                                dy_lvl(lvl), direction_flag, M+2)
        else
            call boundary_ghost_cell_calculate(p_BC_top, p_top, dp_dy_top, &
                                                tree%blocks(idx)%p(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                                tree%blocks(idx)%p(:,j_end    ,1), &    ! Interior array
                                                dy_lvl(lvl), direction_flag, M+2)
        end if
    end do

end subroutine apply_boundary_condition_pressure

subroutine apply_boundary_condition_pressure_multi_grid(tree, mg_lvl_no)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: mg_lvl_no
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION IF mg_lvl_no is 1 (Finest grid)
    if (mg_lvl_no==1) then
        !> Left Boundary Blocks
        direction_flag = -1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_left_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_left_blocks(id)
            lvl = tree%blocks(idx)%lvl

            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_left, p_left, dp_dx_left, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_start    ,:,1), &      ! Interior array
                                        dx_lvl(lvl), direction_flag, N+2)
        end do
        
        !> Right Boundary Blocks
        direction_flag = 1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_right_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_right_blocks(id)
            lvl = tree%blocks(idx)%lvl
            
            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_right, p_right, dp_dx_right, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_end    ,:,1), &    ! Interior array
                                        dx_lvl(lvl), direction_flag, N+2)
        end do
        
        !> Bottom Boundary Blocks
        direction_flag = -1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(id)
            lvl = tree%blocks(idx)%lvl
            
            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_bottom, p_bottom, dp_dy_bottom, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_start    ,1), &      ! Interior array
                                        dy_lvl(lvl), direction_flag, M+2)
        end do
        
        !> Top Boundary Blocks
        direction_flag = 1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_top_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_top_blocks(id)
            lvl = tree%blocks(idx)%lvl

            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_top, p_top, dp_dy_top, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_end    ,1), &    ! Interior array
                                        dy_lvl(lvl), direction_flag, M+2)
        end do

    else
        !> Left Boundary Blocks
        direction_flag = -1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_left_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_left_blocks(id)
            lvl = tree%blocks(idx)%lvl

            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_left, 0.0_prcs_var, 0.0_prcs_var, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_start    ,:,1), &      ! Interior array
                                        dx_lvl(lvl), direction_flag, N+2)
        end do
        
        !> Right Boundary Blocks
        direction_flag = 1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_right_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_right_blocks(id)
            lvl = tree%blocks(idx)%lvl
            
            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_right, 0.0_prcs_var, 0.0_prcs_var, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(i_end    ,:,1), &    ! Interior array
                                        dx_lvl(lvl), direction_flag, N+2)
        end do
        
        !> Bottom Boundary Blocks
        direction_flag = -1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(id)
            lvl = tree%blocks(idx)%lvl
            
            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_bottom, 0.0_prcs_var, 0.0_prcs_var, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_start    ,1), &      ! Interior array
                                        dy_lvl(lvl), direction_flag, M+2)
        end do
        
        !> Top Boundary Blocks
        direction_flag = 1.0_prcs_var
        do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_top_blocks)
            idx = tree%mg_lvl(mg_lvl_no)%MG_top_blocks(id)
            lvl = tree%blocks(idx)%lvl

            !>  --- p - Pressure ---
            call boundary_ghost_cell_calculate(p_BC_top, 0.0_prcs_var, 0.0_prcs_var, &
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                        tree%blocks(idx)%multi_grid(mg_lvl_no)%p(:,j_end    ,1), &    ! Interior array
                                        dy_lvl(lvl), direction_flag, M+2)
        end do
    end if

end subroutine apply_boundary_condition_pressure_multi_grid

!> ------------------- DIRECTION -------------------
subroutine apply_boundary_condition_direction(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION

    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%left_blocks)
        idx = tree%left_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- pressure ---
        call boundary_ghost_cell_calculate(p_BC_left, 0.0_prcs_var, 0.0_prcs_var, &
                                      tree%blocks(idx)%d(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%d(i_start    ,:,1), &      ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
        
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%right_blocks)
        idx = tree%right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- pressure ---
        call boundary_ghost_cell_calculate(p_BC_right, 0.0_prcs_var, 0.0_prcs_var, &
                                      tree%blocks(idx)%d(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%d(i_end    ,:,1), &    ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%bottom_blocks)
        idx = tree%bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- pressure ---
        call boundary_ghost_cell_calculate(p_BC_bottom, 0.0_prcs_var, 0.0_prcs_var, &
                                      tree%blocks(idx)%d(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%d(:,j_start    ,1), &      ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%top_blocks)
        idx = tree%top_blocks(id)
        lvl = tree%blocks(idx)%lvl
        !>  --- pressure ---
        call boundary_ghost_cell_calculate(p_BC_top, 0.0_prcs_var, 0.0_prcs_var, &
                                      tree%blocks(idx)%d(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%d(:,j_end    ,1), &    ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do

end subroutine apply_boundary_condition_direction

!> ------------------- RHS -------------------
subroutine apply_boundary_condition_rhs(tree, mg_lvl_no)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer, intent(in)         :: mg_lvl_no
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION IF mg_lvl_no is 1 (Finest grid)
    
    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_left_blocks)
        idx = tree%mg_lvl(mg_lvl_no)%MG_left_blocks(id)
        lvl = tree%blocks(idx)%lvl

        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(p_BC_left, 0.0_prcs_var, 0.0_prcs_var, &
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i_start    ,:,1), &      ! Interior array
                                    dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_right_blocks)
        idx = tree%mg_lvl(mg_lvl_no)%MG_right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        
        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(p_BC_right, 0.0_prcs_var, 0.0_prcs_var, &
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(i_end    ,:,1), &    ! Interior array
                                    dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks)
        idx = tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        
        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(p_BC_bottom, 0.0_prcs_var, 0.0_prcs_var, &
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(:,j_start    ,1), &      ! Interior array
                                    dy_lvl(lvl), direction_flag, M+2)
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%mg_lvl(mg_lvl_no)%MG_top_blocks)
        idx = tree%mg_lvl(mg_lvl_no)%MG_top_blocks(id)
        lvl = tree%blocks(idx)%lvl

        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(p_BC_top, 0.0_prcs_var, 0.0_prcs_var, &
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                    tree%blocks(idx)%multi_grid(mg_lvl_no)%f(:,j_end    ,1), &    ! Interior array
                                    dy_lvl(lvl), direction_flag, M+2)
    end do

end subroutine apply_boundary_condition_rhs

!> ------------------- TEMPERATURE -------------------
subroutine apply_boundary_condition_temperature(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id
    real(kind=prcs_var) :: direction_flag   !> -1 if left or bottom, +1 if top or right


    !!   BOUNDARY CONDITION FOR DIRICHLET AND NEUMAN CONDITION

    !> Left Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%left_blocks)
        idx = tree%left_blocks(id)
        lvl = tree%blocks(idx)%lvl

        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(T_BC_left, T_left, dT_dx_left, &
                                      tree%blocks(idx)%T(i_start-ngl,:,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%T(i_start    ,:,1), &      ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
        
    end do
    
    !> Right Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%right_blocks)
        idx = tree%right_blocks(id)
        lvl = tree%blocks(idx)%lvl
        
        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(T_BC_right, T_right, dT_dx_right, &
                                      tree%blocks(idx)%T(i_end+ngl,:,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%T(i_end    ,:,1), &    ! Interior array
                                      dx_lvl(lvl), direction_flag, N+2)
    end do
    
    !> Bottom Boundary Blocks
    direction_flag = -1.0_prcs_var
    do id = 1, size(tree%bottom_blocks)
        idx = tree%bottom_blocks(id)
        lvl = tree%blocks(idx)%lvl
        
        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(T_BC_bottom, T_bottom, dT_dy_bottom, &
                                      tree%blocks(idx)%T(:,j_start-ngl,1), &      ! Boundary (ghost) array
                                      tree%blocks(idx)%T(:,j_start    ,1), &      ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do
    
    !> Top Boundary Blocks
    direction_flag = 1.0_prcs_var
    do id = 1, size(tree%top_blocks)
        idx = tree%top_blocks(id)
        lvl = tree%blocks(idx)%lvl

        !>  --- p - Pressure ---
        call boundary_ghost_cell_calculate(T_BC_top, T_top, dT_dy_top, &
                                      tree%blocks(idx)%T(:,j_end+ngl,1), &    ! Boundary (ghost) array
                                      tree%blocks(idx)%T(:,j_end    ,1), &    ! Interior array
                                      dy_lvl(lvl), direction_flag, M+2)
    end do

end subroutine apply_boundary_condition_temperature
