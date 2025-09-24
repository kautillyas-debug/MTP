#include "definitions.h"

! --------------------------------------------------------------
!>        Correcting the direction by subtracting it's mean
!>        Used when P_BC is neumann all around, and CG is now solving non SPD matrix
! --------------------------------------------------------------
subroutine reset_direction_non_SPD(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module, only: no_fluid_cells
    implicit none

    type(t_tree), intent(inout) :: tree
    integer                     :: lvl, id, idx
    integer                     :: i, j, k, i1, l
    real(kind=prcs_var)         :: d_mean

    d_mean = 0.0_prcs_var

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Summing up the values of residue (or direction)
                        d_mean = d_mean + tree%blocks(idx)%d(i,j,k)


                    end do

                end do
            end do


        end do
    end do


    ! Finding the mean of the values of d (or x) and subtracting the value from d (or x)
    d_mean = d_mean / no_fluid_cells

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Resetting mean of Ad to 0
                        tree%blocks(idx)%d(i,j,k) = tree%blocks(idx)%d(i,j,k) - d_mean


                    end do

                end do
            end do
        end do
    end do


end subroutine reset_direction_non_SPD


subroutine reset_residual_non_SPD(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module, only: no_fluid_cells
    implicit none

    type(t_tree), intent(inout) :: tree
    integer                     :: lvl, id, idx
    integer                     :: i, j, k, i1, l
    real(kind=prcs_var)         :: d_mean

    d_mean = 0.0_prcs_var

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Summing up the values of residue (or direction)
                        d_mean = d_mean + tree%blocks(idx)%r(i,j,k)


                    end do

                end do
            end do


        end do
    end do


    ! Finding the mean of the values of d (or x) and subtracting the value from d (or x)
    d_mean = d_mean / no_fluid_cells

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Resetting mean of r to 0
                        tree%blocks(idx)%r(i,j,k) = tree%blocks(idx)%r(i,j,k) - d_mean


                    end do

                end do
            end do
        end do
    end do


end subroutine reset_residual_non_SPD


subroutine reset_Ad_non_SPD(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module, only: no_fluid_cells
    implicit none

    type(t_tree), intent(inout) :: tree
    integer                     :: lvl, id, idx
    integer                     :: i, j, k, i1, l
    real(kind=prcs_var)         :: d_mean

    d_mean = 0.0_prcs_var

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Summing up the values of residue (or direction)
                        d_mean = d_mean + tree%blocks(idx)%Ad(i,j,k)


                    end do

                end do
            end do


        end do
    end do


    ! Finding the mean of the values of d (or x) and subtracting the value from d (or x)
    d_mean = d_mean / no_fluid_cells

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)
            k = 1
            do j = j_start, j_end
                do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
                    l = 2*i1 - 1
                    do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)  
                        
                        ! Resetting mean of Ad to 0
                        tree%blocks(idx)%Ad(i,j,k) = tree%blocks(idx)%Ad(i,j,k) - d_mean


                    end do

                end do
            end do
        end do
    end do


end subroutine reset_Ad_non_SPD


! --------------------------------------------------------------
!>        Enforcing the face vel at object to body_value
! --------------------------------------------------------------
subroutine enforce_no_slip_body_face(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module
    implicit none

    type(t_tree), intent(inout) :: tree
    integer                     :: lvl, id, idx
    integer                     :: i, j, k, VAR, ghost_no, ibody
    real(kind=prcs_var)         :: d_mean

    
    ! Enforcing Face Velocities at cells near object with the body boundary condition
    do lvl = 1, tree%no_levels
        ! do id = 1, tree%levels(lvl)%idx_body(0)
        !     idx = tree%levels(lvl)%idx_body(id)
        do id = 1, tree%levels(lvl)%idx_body_leaf(0)
            idx = tree%levels(lvl)%idx_body_leaf(id)

            do ibody = 1, nbody
            
                VAR = U_VELOCITY
                do ghost_no = tree%blocks(idx)%no_IB_ghost_cells(ibody-1)+1, tree%blocks(idx)%no_IB_ghost_cells(ibody)
                    i = tree%blocks(idx)%IB_ghost_cell(ghost_no, I_g)
                    j = tree%blocks(idx)%IB_ghost_cell(ghost_no, J_g)
                    k = tree%blocks(idx)%IB_ghost_cell(ghost_no, K_g)
                    tree%blocks(idx)%Ue(i-1,j,k) = body_data(ibody)%Boundary_value(VAR)
                    tree%blocks(idx)%Uw(i+1,j,k) = body_data(ibody)%Boundary_value(VAR)
                end do

                VAR = V_VELOCITY
                do ghost_no = tree%blocks(idx)%no_IB_ghost_cells(ibody-1)+1, tree%blocks(idx)%no_IB_ghost_cells(ibody)
                    i = tree%blocks(idx)%IB_ghost_cell(ghost_no, I_g)
                    j = tree%blocks(idx)%IB_ghost_cell(ghost_no, J_g)
                    k = tree%blocks(idx)%IB_ghost_cell(ghost_no, K_g)
                    tree%blocks(idx)%Vn(i,j-1,k) = body_data(ibody)%Boundary_value(VAR)
                    tree%blocks(idx)%Vs(i,j+1,k) = body_data(ibody)%Boundary_value(VAR)
                end do

            end do
                
        end do
    end do

end subroutine enforce_no_slip_body_face