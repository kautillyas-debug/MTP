! --------------------------------------------------------------
!>        Calculate Residual of Variable L2(b-Ax)
! --------------------------------------------------------------
subroutine calculate_residual(variable, rhs, gamma, residual)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start, j_end
    implicit none
    
    real(kind=prcs_var), intent(in)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)     :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(in)     :: gamma            ! Coeff of U_p
    real(kind=prcs_var), intent(inout)  :: residual         ! L2 norm residual 
    integer :: i, j, k
    real(kind=prcs_var)     :: temp_residual

    k = 1
    !> Calculate residual
    temp_residual = 0.0_prcs_var
    do j = j_start, j_end
        do i = i_start, i_end

            #IFDEF SQ_CELL
                temp_residual = temp_residual + (rhs(i,j,k) &
                    - ((variable(i,j-1,k)+variable(i,j+1,k)) &
                    + (variable(i-1,j,k)+variable(i+1,j,k)) &
                    + gamma*variable(i,j,k)))**2
            #ELSE
                temp_residual = temp_residual + (rhs(i,j,k) &
                    - ((variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                    + (variable(i-1,j,k)+variable(i+1,j,k)) &
                    + gamma*variable(i,j,k)))**2
            #ENDIF

        end do
    end do

    residual = residual + temp_residual

end subroutine calculate_residual

subroutine calculate_residual_body(variable, rhs, gamma, residual, i_range)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start, j_end
    implicit none
    
    real(kind=prcs_var), intent(in)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)     :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(in)     :: gamma            ! Coeff of U_p
    real(kind=prcs_var), intent(inout)  :: residual         ! L2 norm residual 
    integer,             intent(in)     :: i_range(0:M+2, N+2, 1)    ! RHS of equation
    integer :: i, j, k, i1, l
    real(kind=prcs_var)     :: temp_residual

    k = 1
    !> Calculate residual
    temp_residual = 0.0_prcs_var
    do j = j_start, j_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k)

                #IFDEF SQ_CELL
                    temp_residual = temp_residual + (rhs(i,j,k) &
                        - ((variable(i,j-1,k)+variable(i,j+1,k)) &
                        + (variable(i-1,j,k)+variable(i+1,j,k)) &
                        + gamma*variable(i,j,k)))**2
                #ELSE
                    temp_residual = temp_residual + (rhs(i,j,k) &
                        - ((variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                        + (variable(i-1,j,k)+variable(i+1,j,k)) &
                        + gamma*variable(i,j,k)))**2
                #ENDIF

            end do

        end do
    end do

    residual = residual + temp_residual

end subroutine calculate_residual_body

subroutine calculate_residue_array(variable, rhs, residue, dx_sq, residual, gamma)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start, j_end
    implicit none
    
    real(kind=prcs_var), intent(in)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)  :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(out) :: residue(M+2, N+2, 1)    ! residue of equation
    real(kind=prcs_var), intent(in)  :: gamma            ! Coeff of U_p
    real(kind=prcs_var), intent(in)  :: dx_sq            ! dx_square
    real(kind=prcs_var), intent(inout)  :: residual            ! L2 norm of Residual
    integer :: i, j, k

    !> Calculate residual
    k = 1
    do j = j_start, j_end
        do i = i_start, i_end

            #IFDEF SQ_CELL
                residue(i,j,k) = (rhs(i,j,k) &
                    - ((variable(i,j-1,k)+variable(i,j+1,k)) &
                    + (variable(i-1,j,k)+variable(i+1,j,k)) &
                    + gamma*variable(i,j,k))) 
                residue(i,j,k) = residue(i,j,k) / dx_sq
                residual = residual + residue(i,j,k)*residue(i,j,k)
            #ELSE
                residue(i,j,k) = (rhs(i,j,k) &
                    - ((variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                    + (variable(i-1,j,k)+variable(i+1,j,k)) &
                    + gamma*variable(i,j,k)))
                residue(i,j,k) = residue(i,j,k) / dx_sq
                residual = residual + residue(i,j,k)*residue(i,j,k)
            #ENDIF

        end do
    end do

end subroutine calculate_residue_array

subroutine calculate_residue_array_body(variable, rhs, residue, dx_sq, residual, gamma, i_range)

    use precision_module, only: prcs_var
    use data_type_module, only: beta_sq, N, M, i_start, i_end, j_start, j_end
    implicit none
    
    real(kind=prcs_var), intent(in)  :: variable(M+2, N+2, 1)  ! Array to be updated
    real(kind=prcs_var), intent(in)  :: rhs(M+2, N+2, 1)    ! RHS of equation
    real(kind=prcs_var), intent(out) :: residue(M+2, N+2, 1)    ! residue of equation
    real(kind=prcs_var), intent(in)  :: gamma            ! Coeff of U_p
    real(kind=prcs_var), intent(in)  :: dx_sq            ! dx_square
    real(kind=prcs_var), intent(inout)  :: residual            ! L2 norm of Residual
    integer,             intent(in)     :: i_range(0:M+2, N+2, 1)    ! RHS of equation
    integer :: i, j, k, i1, l

    !> Calculate residual
    k = 1
    do j = j_start, j_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k)

                #IFDEF SQ_CELL
                    residue(i,j,k) = (rhs(i,j,k) &
                        - ((variable(i,j-1,k)+variable(i,j+1,k)) &
                        + (variable(i-1,j,k)+variable(i+1,j,k)) &
                        + gamma*variable(i,j,k))) 
                    residue(i,j,k) = residue(i,j,k) / dx_sq
                    residual = residual + residue(i,j,k)*residue(i,j,k)
                #ELSE
                    residue(i,j,k) = (rhs(i,j,k) &
                        - ((variable(i,j-1,k)+variable(i,j+1,k)) * beta_sq &
                        + (variable(i-1,j,k)+variable(i+1,j,k)) &
                        + gamma*variable(i,j,k)))
                    residue(i,j,k) = residue(i,j,k) / dx_sq
                    residual = residual + residue(i,j,k)*residue(i,j,k)
                #ENDIF
            end do
        end do
    end do

end subroutine calculate_residue_array_body
