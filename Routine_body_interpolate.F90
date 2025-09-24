
    #include "definitions.h"

! --------------------------------------------------------------
!         Interpolation of variables at IP from neighbours
! --------------------------------------------------------------

! !> \todo (Need to interpolate from 4 points near the IP)
! !> \brief Interpolate the value of Pressure at IP from the neighbouring points
! !> \param[in] IP_interpol_cells Array containing the interpolation cells, with dimensions (3, 4, no_IB_ghost_cells).
! !> \param[in] variable          A 3D array of variables with dimensions (M, N, 1), representing the computational data.
! !> \param[in] no_IB_ghost_cells Number of immersed boundary ghost cells used in the interpolation process.
! !> \param[in,out] IP_pos_val    Array of interpolated values at specific positions, modified in-place, with size (no_IB_ghost_cells).
! !> \param[in] VAR Integer value determining which variable is used u-1, v-2, p-3, T-4
subroutine interpolate_IP_variable(IP_pos_val, IP_interpol_cells, variable, no_IB_ghost_cells, VAR)


    use precision_module, only: prcs_var
    use data_type_module, only: M, N, ngl, no_interp_pts_body

    implicit none
    integer,             intent(in)    :: no_IB_ghost_cells
    integer, intent(in)    :: IP_interpol_cells(3, no_interp_pts_body, no_IB_ghost_cells)
    real(kind=prcs_var), intent(in)    :: variable(M+2*ngl, N+2*ngl, 1)
    real(kind=prcs_var), intent(inout) :: IP_pos_val(no_IB_ghost_cells)
    integer,             intent(in)     :: VAR

    integer :: i, j, k
    integer :: ghost_no   ! number of ghost cells in the block; iterating through ghost cells
    integer :: pt_no      ! Iterate through the interpolation points
    real(kind=prcs_var) :: interpolation_sum ! sum the value of each point with it's coefficient


    ! Finding the image point of the blocks
    do ghost_no = 1, no_IB_ghost_cells

        ! Iterate  through all the points for interpolation
        interpolation_sum = 0.0_prcs_var

        if (VAR == PRESSURE) then
            do pt_no = 1, no_interp_pts_body
                i = IP_interpol_cells(I_g, pt_no, ghost_no) 
                j = IP_interpol_cells(J_g, pt_no, ghost_no) 
                k = IP_interpol_cells(K_g, pt_no, ghost_no)  

                ! print *, i, IP_interpol_cells(I_g, :, ghost_no)

                ! Currently Coefficients are considered as 1, would need to change later based on Vandermonde Matrix
                ! interpolation_sum = 1.0_prcs_var * variable(i,j,k)
                interpolation_sum = interpolation_sum + 0.5_prcs_var * variable(i,j,k)

            end do

        else

            do pt_no = 1, no_interp_pts_body
                i = IP_interpol_cells(I_g, 1, ghost_no) 
                j = IP_interpol_cells(J_g, 1, ghost_no) 
                k = IP_interpol_cells(K_g, 1, ghost_no)  

                ! print *, i, IP_interpol_cells(I_g, :, ghost_no)

                ! Currently Coefficients are considered as 1, would need to change later based on Vandermonde Matrix
                ! interpolation_sum = 1.0_prcs_var * variable(i,j,k)
                interpolation_sum = interpolation_sum + 0.5_prcs_var * variable(i,j,k)

            end do


        end if

        IP_pos_val(ghost_no) = interpolation_sum

    end do


end subroutine interpolate_IP_variable

