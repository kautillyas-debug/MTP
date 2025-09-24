
    #include "definitions.h"

! --------------------------------------------------------------
!>        Update Intermediate Face Velocities
! --------------------------------------------------------------
subroutine calculate_intermediate_face_velocities_block(u, v, Ue, Uw, Vn, Vs, i_range)
    use precision_module, only: prcs_var
    use data_type_module, only: M, N
    implicit none

    real(kind=prcs_var), dimension(M+2,N+2,1), intent(inout)    ::  Ue, Uw, Vn, Vs
    real(kind=prcs_var), dimension(M+2,N+2,1), intent(in)       ::  u, v
    integer, dimension(0:M+2,N+2,1), intent(in)                 ::  i_range
    integer     ::  i, j, k, i1, l

    k = 1
    do j = 2, N+1
        ! do i = i_start, i_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k) 

                Ue(i,j,k) = 0.5_prcs_var*(u(i,j,k) + u(i+1,j,k))
                Uw(i,j,k) = 0.5_prcs_var*(u(i,j,k) + u(i-1,j,k))
                Vn(i,j,k) = 0.5_prcs_var*(v(i,j,k) + v(i,j+1,k))
                Vs(i,j,k) = 0.5_prcs_var*(v(i,j,k) + v(i,j-1,k))

            end do
        end do
    end do

end subroutine calculate_intermediate_face_velocities_block

subroutine calculate_intermediate_face_velocities(tree)
    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: lvl, id, idx, nb_idx, ch_id_1, ch_id_2
    integer :: i, j, k, i1, l, fine_cell
    integer :: ibody, ghost_no, VAR
    real(kind=prcs_var) :: Q_out

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            call calculate_intermediate_face_velocities_block( &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%Ue, &
                tree%blocks(idx)%Uw, &
                tree%blocks(idx)%Vn, &
                tree%blocks(idx)%Vs, &
                tree%blocks(idx)%i_range &
            )

            ! Enforcing Face Velocity of Coarse Block to be interpolated from fine blocks
                ! Left block is refined:
                nb_idx = tree%blocks(idx)%ptr_nb_blocks(2)
                if (nb_idx/=0) then
                    if (tree%blocks(nb_idx)%is_refined) then
                        ch_id_1 = tree%blocks(nb_idx)%ptr_children(3)
                        ch_id_2 = tree%blocks(nb_idx)%ptr_children(4)
                        do j = j_start, j_mid
                            fine_cell = 2*(j-1)
                            tree%blocks(idx)%Uw(i_start,j,1) = (tree%blocks(ch_id_1)%Ue(i_end,fine_cell  ,1) &
                                                    + tree%blocks(ch_id_1)%Ue(i_end,fine_cell+1,1))/2.0_prcs_var
                        end do
                        do j = j_mid+1, j_end
                            fine_cell = 2*(j-j_mid)
                            tree%blocks(idx)%Uw(i_start,j,1) = (tree%blocks(ch_id_2)%Ue(i_end,fine_cell  ,1) &
                                                    + tree%blocks(ch_id_2)%Ue(i_end,fine_cell+1,1))/2.0_prcs_var
                        end do
                    end if
                end if
                
                ! Right block is refined:
                nb_idx = tree%blocks(idx)%ptr_nb_blocks(6)
                if (nb_idx/=0) then
                    if (tree%blocks(nb_idx)%is_refined) then
                        ch_id_1 = tree%blocks(nb_idx)%ptr_children(1)
                        ch_id_2 = tree%blocks(nb_idx)%ptr_children(2)
                        do j = j_start, j_mid
                            fine_cell = 2*(j-1)
                            tree%blocks(idx)%Ue(i_end,j,1) = (tree%blocks(ch_id_1)%Uw(i_start,fine_cell  ,1) &
                                                    + tree%blocks(ch_id_1)%Uw(i_start,fine_cell+1,1))/2.0_prcs_var
                        end do
                        do j = j_mid+1, j_end
                            fine_cell = 2*(j-j_mid)
                            tree%blocks(idx)%Ue(i_end,j,1) = (tree%blocks(ch_id_2)%Uw(i_start,fine_cell  ,1) &
                                                    + tree%blocks(ch_id_2)%Uw(i_start,fine_cell+1,1))/2.0_prcs_var
                        end do
                    end if
                end if

                ! Bottom block is refined:
                nb_idx = tree%blocks(idx)%ptr_nb_blocks(8)
                if (nb_idx/=0) then
                    if (tree%blocks(nb_idx)%is_refined) then
                        ch_id_1 = tree%blocks(nb_idx)%ptr_children(2)
                        ch_id_2 = tree%blocks(nb_idx)%ptr_children(4)
                        do i = i_start, i_mid
                            fine_cell = 2*(i-1)
                            tree%blocks(idx)%Vs(i,j_start,1) = (tree%blocks(ch_id_1)%Vn(fine_cell  ,j_end,1) &
                                                    + tree%blocks(ch_id_1)%Vn(fine_cell+1,j_end,1))/2.0_prcs_var
                        end do
                        do i = i_mid+1, i_end
                            fine_cell = 2*(i-i_mid)
                            tree%blocks(idx)%Vs(i,j_start,1) = (tree%blocks(ch_id_2)%Vn(fine_cell  ,j_end,1) &
                                                    + tree%blocks(ch_id_2)%Vn(fine_cell+1,j_end,1))/2.0_prcs_var
                        end do
                    end if
                end if

                ! Top block is refined:
                nb_idx = tree%blocks(idx)%ptr_nb_blocks(4)
                if (idx==1) then
                end if
                if (nb_idx/=0) then
                    if (tree%blocks(nb_idx)%is_refined) then
                        ch_id_1 = tree%blocks(nb_idx)%ptr_children(1)
                        ch_id_2 = tree%blocks(nb_idx)%ptr_children(3)
                        do i = i_start, i_mid
                            fine_cell = 2*(i-1)
                            tree%blocks(idx)%Vn(i,j_end,1) = (tree%blocks(ch_id_1)%Vs(fine_cell  ,j_start,1) &
                                                    + tree%blocks(ch_id_1)%Vs(fine_cell+1,j_start,1))/2.0_prcs_var
                        end do
                        do i = i_mid+1, i_end
                            fine_cell = 2*(i-i_mid)
                            tree%blocks(idx)%Vn(i,j_end,1) = (tree%blocks(ch_id_2)%Vs(fine_cell  ,j_start,1) &
                                                    + tree%blocks(ch_id_2)%Vs(fine_cell+1,j_start,1))/2.0_prcs_var
                        end do
                    end if
                end if      
            !
            
        end do
    end do
    
    ! Enforcing 
    call enforce_no_slip_body_face(tree)

    ! checking the net mass flow into domain, and set the outlet face velocity Ue such that the net flow is zero
    if (pressure_BC_neumann .and. u_BC_right==NEUMANN) then
        
        Q_out = 0.0_prcs_var
        
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                k = 1
                do j = 2, N+1
                    ! do i = i_start, i_end
                    do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                        l = 2*i1 - 1
                        do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  

                            Q_out = Q_out + ((tree%blocks(idx)%Ue(i,j,k)-tree%blocks(idx)%Uw(i,j,k)) * dy_lvl(lvl) &
                                          +  (tree%blocks(idx)%Vn(i,j,k)-tree%blocks(idx)%Vs(i,j,k)) * dx_lvl(lvl))

                        end do
                    end do
                end do

            end do
        end do
        
        print *, '            Q_out before correction: ', Q_out

        ! Assuming the outlet is at the right boundary:
        do id = 1, size(tree%right_blocks)
            idx = tree%right_blocks(id)
            lvl = tree%blocks(idx)%lvl

            ! tree%blocks(idx)%Ue(i_end,:,1) = tree%blocks(idx)%Ue(i_end,:,1) - (Q_out / (dy_lvl(lvl)/min(1.0_prcs_var, domain(2)/domain(1))))
            tree%blocks(idx)%Ue(i_end,:,1) = tree%blocks(idx)%Ue(i_end,:,1) - (Q_out / domain(2))
            tree%blocks(idx)%u(i_end,:,1) = tree%blocks(idx)%Ue(i_end,:,1)
            
        end do
        
    end if


end subroutine calculate_intermediate_face_velocities


! --------------------------------------------------------------
!>        Update Cell-Centred and Face-Centred Velocities
! --------------------------------------------------------------
subroutine update_cell_face_velocities_block(u, v, Ue, Uw, Vn, Vs, p, dx, dy, i_range)
    use precision_module, only: prcs_var
    use data_type_module, only: M, N
    use problem_module, only: delta_time
    implicit none

    real(kind=prcs_var), dimension(M+2,N+2,1), intent(inout)    ::  u, v, Ue, Uw, Vn, Vs
    real(kind=prcs_var), dimension(M+2,N+2,1), intent(in)       ::  p
    integer, dimension(0:M+2,N+2,1), intent(in)                 ::  i_range
    real(kind=prcs_var), intent(in)       ::  dx, dy
    real(kind=prcs_var)                   ::  coeff_up, coeff_vp, coeff_U, coeff_V
    integer     ::  i, j, k, i1, l

    coeff_U = delta_time/dx
    coeff_V = delta_time/dy
    coeff_up = coeff_U/2.0_prcs_var
    coeff_vp = coeff_V/2.0_prcs_var


    k = 1
    do j = 2, N+1
        ! do i = i_start, i_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k)  

                u(i,j,k)  = u(i,j,k) - coeff_up*(p(i+1,j,k)-p(i-1,j,k))
                v(i,j,k)  = v(i,j,k) - coeff_vp*(p(i,j+1,k)-p(i,j-1,k))
                Ue(i,j,k) = Ue(i,j,k) - coeff_U*(p(i+1,j,k)-p(i,j,k))
                Uw(i,j,k) = Uw(i,j,k) - coeff_U*(p(i,j,k)-p(i-1,j,k))
                Vn(i,j,k) = Vn(i,j,k) - coeff_V*(p(i,j+1,k)-p(i,j,k))
                Vs(i,j,k) = Vs(i,j,k) - coeff_V*(p(i,j,k)-p(i,j-1,k))

            end do


        end do
    end do

    ! i = 3
    ! j = 4
    ! k = 1

    ! print *, coeff_vp
    ! print *, v(i,j,k), Vn(i,j,k), Vs(i,j,k)
    ! print *, p(i,j+1,k), p(i,j-1,k)

end subroutine update_cell_face_velocities_block

subroutine update_cell_face_velocities(tree)
    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: lvl, id, idx

    #IFDEF VERBOSE
        if (mod(time_step+1,STEP_SKIP)==0) then
            print *
            print *, '--------------------------------------------------'
            print *, '  STEP 3 : Update Face and Cell Centred Velocities'
            print *, '--------------------------------------------------'
        end if
    #ENDIF

    ! print *, 'i_range:::::'
    ! call print_screen_real(tree%blocks(1)%i_range, 'ir', 0, M+2, 1, N+2)

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            
            if (PPE_SOLVER == MULTI_GRID) then
                call update_cell_face_velocities_block( &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%Ue, &
                    tree%blocks(idx)%Uw, &
                    tree%blocks(idx)%Vn, &
                    tree%blocks(idx)%Vs, &
                    tree%blocks(idx)%multi_grid(1)%p, &
                    dx_lvl(lvl), &
                    dy_lvl(lvl), &
                    tree%blocks(idx)%i_range &
                    )
            
            else 
                call update_cell_face_velocities_block( &
                    tree%blocks(idx)%u, &
                    tree%blocks(idx)%v, &
                    tree%blocks(idx)%Ue, &
                    tree%blocks(idx)%Uw, &
                    tree%blocks(idx)%Vn, &
                    tree%blocks(idx)%Vs, &
                    tree%blocks(idx)%p, &
                    dx_lvl(lvl), &
                    dy_lvl(lvl), &
                    tree%blocks(idx)%i_range &
                    )
            end if

        end do
    end do

end subroutine update_cell_face_velocities