! --------------------------------------------------------------
!         Update Ghost cells for IB arrays
! --------------------------------------------------------------

! !> \brief Update the ghost cells from the IB variables created
! !> \param[in] IP_pos_val    Array of Image Point values at specific positions
! !> \param[in] IB_ghost_cell   Array ghost cell indices of each body
! !> \param[in,out] variable  A 3D array of variables with dimensions (M, N, 1), representing the computational data.
! !> \param[in] no_IB_ghost_cells Number of immersed boundary ghost cells for each body used
! !> \param[in] VAR Integer value determining which variable is used u-1, v-2, p-3, T-4
! !> \param[in] dl_GC_IP The outward length from GC to IP
subroutine update_IB_ghost_cells_variable(&
        IP_pos_val, IB_ghost_cell, variable, no_IB_ghost_cells, &
        VAR, dl_GC_IP & 
    )

    #include "definitions.h"
    use precision_module, only: prcs_var
    use data_type_module, only: M, N, ngl
    use body_module

    implicit none
    integer,             intent(in)     :: no_IB_ghost_cells(0:nbody)
    real(kind=prcs_var), intent(in)     :: IP_pos_val(no_IB_ghost_cells(nbody),7)
    integer,             intent(in)     :: IB_ghost_cell(no_IB_ghost_cells(nbody), 4)
    real(kind=prcs_var), intent(inout)  :: variable(M+2*ngl, N+2*ngl, 1)
    integer,             intent(in)     :: VAR
    real(kind=prcs_var), intent(in)     :: dl_GC_IP(no_IB_ghost_cells(nbody))

    integer :: i, j, k
    integer :: ibody      ! iterate through nbody
    integer :: ghost_no   ! number of ghost cells in the block; iterating through ghost cells


    ! ! For each body in the domain
    do ibody = 1, nbody

        if (body_data(ibody)%Boundary_condition(VAR) == DIRICHLET) then
            do ghost_no = no_IB_ghost_cells(ibody-1)+1, no_IB_ghost_cells(ibody)

                i = IB_ghost_cell(ghost_no, I_g)
                j = IB_ghost_cell(ghost_no, J_g)
                k = IB_ghost_cell(ghost_no, K_g)

                ! print *, ibody, i, j, k, IP_pos_val(ghost_no,3+VAR)

                variable(i,j,k) = 2.0_prcs_var*body_data(ibody)%Boundary_value(VAR) - IP_pos_val(ghost_no,3+VAR)
                ! variable(i,j,k) = 0.0_prcs_var
                

            end do

        else if (body_data(ibody)%Boundary_condition(VAR) == NEUMANN) then
            do ghost_no = no_IB_ghost_cells(ibody-1)+1, no_IB_ghost_cells(ibody)

                i = IB_ghost_cell(ghost_no, I_g)
                j = IB_ghost_cell(ghost_no, J_g)
                k = IB_ghost_cell(ghost_no, K_g)

                variable(i,j,k) = - dl_GC_IP(ghost_no) * body_data(ibody)%Boundary_gradient(VAR)  &
                                  + IP_pos_val(ghost_no,3+VAR)
                ! variable(i,j,k) = 0.0_prcs_var


            end do

        end if

    end do
    
end subroutine update_IB_ghost_cells_variable

! !> \brief Update the ghost cells from the IB variables created
! !> \param[in] IP_pos_val    Array of Image Point values at specific positions
! !> \param[in] IB_ghost_cell   Array ghost cell indices of each body
! !> \param[in,out] direction  A 3D array of variables with dimensions (M, N, 1), representing the computational data.
! !> \param[in] no_IB_ghost_cells Number of immersed boundary ghost cells for each body used
! !> \param[in] VAR Integer value determining which variable is used u-1, v-2, p-3, T-4
! !> \param[in] dl_GC_IP The outward length from GC to IP
subroutine update_IB_ghost_cells_direction(&
        IP_pos_val, IB_ghost_cell, direction, no_IB_ghost_cells, &
        VAR, dl_GC_IP & 
    )

    #include "definitions.h"
    use precision_module, only: prcs_var
    use data_type_module, only: M, N, ngl
    use body_module

    implicit none
    integer,             intent(in)     :: no_IB_ghost_cells(0:nbody)
    real(kind=prcs_var), intent(in)     :: IP_pos_val(no_IB_ghost_cells(nbody),7)
    integer,             intent(in)     :: IB_ghost_cell(no_IB_ghost_cells(nbody), 4)
    real(kind=prcs_var), intent(inout)  :: direction(M+2*ngl, N+2*ngl, 1)
    integer,             intent(in)     :: VAR
    real(kind=prcs_var), intent(in)     :: dl_GC_IP(no_IB_ghost_cells(nbody))

    integer :: i, j, k
    integer :: ibody      ! iterate through nbody
    integer :: ghost_no   ! number of ghost cells in the block; iterating through ghost cells


    ! ! For each body in the domain
    do ibody = 1, nbody

        if (body_data(ibody)%Boundary_condition(VAR) == DIRICHLET) then
            do ghost_no = no_IB_ghost_cells(ibody-1)+1, no_IB_ghost_cells(ibody)

                i = IB_ghost_cell(ghost_no, I_g)
                j = IB_ghost_cell(ghost_no, J_g)
                k = IB_ghost_cell(ghost_no, K_g)

                ! print *, ibody, i, j, k, IP_pos_val(ghost_no,3+VAR)

                direction(i,j,k) = - IP_pos_val(ghost_no,3+VAR)
                ! direction(i,j,k) = 2.0_prcs_var*body_data(ibody)%Boundary_value(VAR) - IP_pos_val(ghost_no,3+VAR)
                ! variable(i,j,k) = 0.0_prcs_var
                

            end do

        else if (body_data(ibody)%Boundary_condition(VAR) == NEUMANN) then
            do ghost_no = no_IB_ghost_cells(ibody-1)+1, no_IB_ghost_cells(ibody)

                i = IB_ghost_cell(ghost_no, I_g)
                j = IB_ghost_cell(ghost_no, J_g)
                k = IB_ghost_cell(ghost_no, K_g)

                direction(i,j,k) = IP_pos_val(ghost_no,3+VAR)
                ! variable(i,j,k) = - dl_GC_IP(ghost_no) * body_data(ibody)%Boundary_gradient(VAR)  &
                !                   + IP_pos_val(ghost_no,3+VAR)
                ! variable(i,j,k) = 0.0_prcs_var


            end do

        end if

    end do
    
end subroutine update_IB_ghost_cells_direction