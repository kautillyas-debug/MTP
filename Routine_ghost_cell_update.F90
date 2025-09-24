! --------------------------------------------------------------
!         Ghost Cell Updation
! --------------------------------------------------------------
subroutine update_ghost_cell_u_vel(tree, idx)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(inout) :: idx
    ! integer, intent(in) :: iter
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    real(kind=prcs_var) :: sum_intrp

    ! ---------- Left values ----------

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%u(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%u(1,block_counter,1) = sum_intrp
        

    end do

    ! ---------- Right values ---------- 

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%u(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

            cell_counter = cell_counter + 1


        end do

        tree%blocks(idx)%u(M+2,block_counter,1) = sum_intrp

    end do


    ! ---------- Bottom values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%u(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%u(block_counter+1,1,1) = sum_intrp

    end do


    ! ---------- Top values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%u(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%u(block_counter+1,N+2,1) = sum_intrp

    end do


end subroutine update_ghost_cell_u_vel

subroutine update_ghost_cell_v_vel(tree, idx)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(inout) :: idx
    ! integer, intent(in) :: iter
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    real(kind=prcs_var) :: sum_intrp

    ! ---------- Left values ----------

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%v(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%v(1,block_counter,1) = sum_intrp

        

    end do

    ! ---------- Right values ---------- 

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%v(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

            cell_counter = cell_counter + 1


        end do

        tree%blocks(idx)%v(M+2,block_counter,1) = sum_intrp

    end do


    ! ---------- Bottom values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%v(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%v(block_counter+1,1,1) = sum_intrp

    end do


    ! ---------- Top values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%v(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%v(block_counter+1,N+2,1) = sum_intrp

    end do


end subroutine update_ghost_cell_v_vel

subroutine update_ghost_cell_pressure(tree, idx)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(inout) :: idx
    ! integer, intent(in) :: iter
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    real(kind=prcs_var) :: sum_intrp

    ! ---------- Left values ----------

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%p(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%p(1,block_counter,1) = sum_intrp

        

    end do

    ! ---------- Right values ---------- 

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%p(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

            cell_counter = cell_counter + 1


        end do

        tree%blocks(idx)%p(M+2,block_counter,1) = sum_intrp

    end do


    ! ---------- Bottom values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%p(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%p(block_counter+1,1,1) = sum_intrp

    end do


    ! ---------- Top values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%p(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%p(block_counter+1,N+2,1) = sum_intrp

    end do
end subroutine update_ghost_cell_pressure

subroutine update_ghost_cell_pressure_multi_grid(tree, idx, mg_lvl_no)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(in)     :: mg_lvl_no
    integer, intent(inout)  :: idx
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    integer :: max_tree_lvl
    real(kind=prcs_var) :: sum_intrp

    max_tree_lvl = tree%no_levels - mg_lvl_no + 1

    if (tree%blocks(idx)%lvl < max_tree_lvl) then
    ! ---------- Left values ----------

        cell_counter = 1
        do block_counter = 1, N+2

            sum_intrp = 0.0_prcs_var
            if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then
                cell_counter = cell_counter + 1
                cycle
            end if
            
            do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)
                block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
                cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
                cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)
                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)
                cell_counter = cell_counter + 1
            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(1,block_counter,1) = sum_intrp

        end do

    ! ---------- Right values ---------- 

        cell_counter = 1
        do block_counter = 1, N+2

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

                block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
                cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
                cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

                cell_counter = cell_counter + 1


            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(M+2,block_counter,1) = sum_intrp

        end do


    ! ---------- Bottom values ---------- 

        cell_counter = 1
        do block_counter = 1, M

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

                block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
                cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
                cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

                cell_counter = cell_counter + 1

            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(block_counter+1,1,1) = sum_intrp

        end do


    ! ---------- Top values ---------- 

        cell_counter = 1
        do block_counter = 1, M

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

                block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
                cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
                cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

                cell_counter = cell_counter + 1

            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(block_counter+1,N+2,1) = sum_intrp

        end do
    else !(when block at the max_tree_lvl of mg_lvl_no)
    ! ---------- Left values ----------

        
        cell_counter = 1
        do block_counter = 1, N+2

            sum_intrp = 0.0_prcs_var
            if(tree%blocks(idx)%MG_nb_block_intrp_left(cell_counter) == 0) then
                cell_counter = cell_counter + 1
                cycle
            end if
            do i = 1, tree%blocks(idx)%MG_no_nb_cells_left(block_counter)

                block_id = tree%blocks(idx)%MG_nb_block_intrp_left(cell_counter)
                cell_x   = tree%blocks(idx)%MG_nb_cells_intrp_left(1,cell_counter)
                cell_y   = tree%blocks(idx)%MG_nb_cells_intrp_left(2,cell_counter)
                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%MG_nb_cells_coeff_left(cell_counter)
                cell_counter = cell_counter + 1
            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(1,block_counter,1) = sum_intrp

        end do

    ! ---------- Right values ---------- 

        cell_counter = 1
        do block_counter = 1, N+2

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%MG_nb_block_intrp_right(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%MG_no_nb_cells_right(block_counter)

                block_id = tree%blocks(idx)%MG_nb_block_intrp_right(cell_counter)
                cell_x   = tree%blocks(idx)%MG_nb_cells_intrp_right(1,cell_counter)
                cell_y   = tree%blocks(idx)%MG_nb_cells_intrp_right(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%MG_nb_cells_coeff_right(cell_counter)

                cell_counter = cell_counter + 1


            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(M+2,block_counter,1) = sum_intrp

        end do


    ! ---------- Bottom values ---------- 

        cell_counter = 1
        do block_counter = 1, M

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%MG_nb_block_intrp_bottom(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%MG_no_nb_cells_bottom(block_counter)

                block_id = tree%blocks(idx)%MG_nb_block_intrp_bottom(cell_counter)
                cell_x   = tree%blocks(idx)%MG_nb_cells_intrp_bottom(1,cell_counter)
                cell_y   = tree%blocks(idx)%MG_nb_cells_intrp_bottom(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%MG_nb_cells_coeff_bottom(cell_counter)

                cell_counter = cell_counter + 1

            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(block_counter+1,1,1) = sum_intrp

        end do


    ! ---------- Top values ---------- 

        cell_counter = 1
        do block_counter = 1, M

            sum_intrp = 0.0_prcs_var

            if(tree%blocks(idx)%MG_nb_block_intrp_top(cell_counter) == 0) then

                cell_counter = cell_counter + 1
                cycle

            end if

            do i = 1, tree%blocks(idx)%MG_no_nb_cells_top(block_counter)

                block_id = tree%blocks(idx)%MG_nb_block_intrp_top(cell_counter)
                cell_x   = tree%blocks(idx)%MG_nb_cells_intrp_top(1,cell_counter)
                cell_y   = tree%blocks(idx)%MG_nb_cells_intrp_top(2,cell_counter)

                sum_intrp = sum_intrp &
                            + tree%blocks(block_id)%multi_grid(mg_lvl_no)%p(cell_x, cell_y, 1) &
                            * tree%blocks(idx)%MG_nb_cells_coeff_top(cell_counter)

                cell_counter = cell_counter + 1

            end do

            tree%blocks(idx)%multi_grid(mg_lvl_no)%p(block_counter+1,N+2,1) = sum_intrp

        end do

    end if

end subroutine update_ghost_cell_pressure_multi_grid

subroutine update_ghost_cell_temperature(tree, idx)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(inout) :: idx
    ! integer, intent(in) :: iter
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    real(kind=prcs_var) :: sum_intrp

    ! ---------- Left values ----------

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%T(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%T(1,block_counter,1) = sum_intrp
        

    end do

    ! ---------- Right values ---------- 

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%T(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

            cell_counter = cell_counter + 1


        end do

        tree%blocks(idx)%T(M+2,block_counter,1) = sum_intrp

    end do


    ! ---------- Bottom values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%T(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%T(block_counter+1,1,1) = sum_intrp

    end do


    ! ---------- Top values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%T(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%T(block_counter+1,N+2,1) = sum_intrp

    end do


end subroutine update_ghost_cell_temperature

subroutine update_ghost_cell_direction(tree, idx)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer, intent(inout) :: idx
    ! integer, intent(in) :: iter
    integer :: block_counter, cell_counter
    integer :: i, block_id, cell_x, cell_y
    real(kind=prcs_var) :: sum_intrp

    ! ---------- Left values ----------
    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_left(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_left(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_left(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_left(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_left(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%d(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_left(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%d(1,block_counter,1) = sum_intrp

        

    end do

    ! ---------- Right values ---------- 

    cell_counter = 1
    do block_counter = 1, N+2

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_right(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_right(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_right(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_right(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_right(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%d(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_right(cell_counter)

            cell_counter = cell_counter + 1


        end do

        tree%blocks(idx)%d(M+2,block_counter,1) = sum_intrp

    end do


    ! ---------- Bottom values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_bottom(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_bottom(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_bottom(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_bottom(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_bottom(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%d(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_bottom(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%d(block_counter+1,1,1) = sum_intrp

    end do


    ! ---------- Top values ---------- 

    cell_counter = 1
    do block_counter = 1, M

        sum_intrp = 0.0_prcs_var

        if(tree%blocks(idx)%nb_block_intrp_top(cell_counter) == 0) then

            cell_counter = cell_counter + 1
            cycle

        end if

        do i = 1, tree%blocks(idx)%no_nb_cells_top(block_counter)

            block_id = tree%blocks(idx)%nb_block_intrp_top(cell_counter)
            cell_x   = tree%blocks(idx)%nb_cells_intrp_top(1,cell_counter)
            cell_y   = tree%blocks(idx)%nb_cells_intrp_top(2,cell_counter)

            sum_intrp = sum_intrp &
                        + tree%blocks(block_id)%d(cell_x, cell_y, 1) * tree%blocks(idx)%nb_cells_coeff_top(cell_counter)

            cell_counter = cell_counter + 1

        end do

        tree%blocks(idx)%d(block_counter+1,N+2,1) = sum_intrp

    end do


end subroutine update_ghost_cell_direction
