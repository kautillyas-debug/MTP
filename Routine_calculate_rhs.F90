
! --------------------------------------------------------------
!>        Calculate RHS
! --------------------------------------------------------------
subroutine calculate_RHS_momentum(rhs, vel, Ue, Uw, Vn, Vs, Cp, alpha, dx, dy, dx_sq, dy_sq, coeff, i_range)
    
    use precision_module, only: prcs_var
    use data_type_module, only: N, M, i_start, i_end, j_start, j_end
    use problem_module,   only: Re
    implicit none
    
    real(kind=prcs_var), intent(inout), dimension(M+2, N+2, 1)  ::  rhs, Cp
    real(kind=prcs_var), intent(in), dimension(M+2, N+2, 1)     ::  vel, Ue, Uw, Vn, Vs 
    integer, dimension(0:M+2,N+2,1), intent(in)                 ::  i_range
    real(kind=prcs_var), intent(in) ::  alpha, dx, dy, dx_sq, dy_sq, coeff
    real(kind=prcs_var) ::  Cp_temp, Dp_temp
    integer         ::  i, j, k, i1, l

    k = 1
    do j = 2, N+1
        ! do i = i_start, i_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k) 

                ! Convective Term
                Cp_temp = ( ( (vel(i+1,j,k)+vel(i,j,k))*Ue(i,j,k) - (vel(i,j,k)+vel(i-1,j,k))*Uw(i,j,k) )/dx &
                          + ( (vel(i,j+1,k)+vel(i,j,k))*Vn(i,j,k) - (vel(i,j,k)+vel(i,j-1,k))*Vs(i,j,k) )/dy )*0.5_prcs_var

                ! Diffusive Term
                Dp_temp = ((vel(i+1,j,k)+vel(i-1,j,k)-2.0_prcs_var*vel(i,j,k))/dx_sq &
                         + (vel(i,j+1,k)+vel(i,j-1,k)-2.0_prcs_var*vel(i,j,k))/dy_sq) / Re

                ! Calculate the RHS
                rhs(i,j,k)  = coeff*vel(i,j,k) + alpha*(3.0_prcs_var*Cp_temp - Cp(i,j,k) - Dp_temp)

                Cp(i,j,k)   = Cp_temp

            end do
        end do
    end do

end subroutine calculate_RHS_momentum

subroutine calculate_RHS_pressure(rhs, Ue, Uw, Vn, Vs, coeff, sum_rhs, i_range)
    
    use precision_module, only: prcs_var
    use data_type_module, only: N, M, beta
    implicit none
    
    real(kind=prcs_var), intent(inout), dimension(M+2, N+2, 1)  ::  rhs
    real(kind=prcs_var), intent(in), dimension(M+2, N+2, 1)     ::  Ue, Uw, Vn, Vs 
    integer, dimension(0:M+2,N+2,1), intent(in)                 ::  i_range
    real(kind=prcs_var), intent(in)    ::  coeff
    real(kind=prcs_var), intent(inout) ::  sum_rhs
    integer :: i, j, k, i1, l

    rhs = 0.0_prcs_var

    k = 1
    do j = 2, N+1
        ! do i = i_start, i_end
        do i1 = 1, i_range(0,j,k)
            l = 2*i1 - 1
            do i = i_range(l,j,k), i_range(l+1,j,k)  

                #IFDEF SQ_CELL
                    rhs(i,j,k) = coeff*((Ue(i,j,k)-Uw(i,j,k)) + (Vn(i,j,k)-Vs(i,j,k)))
                #ELSE
                    rhs(i,j,k) = coeff*((Ue(i,j,k)-Uw(i,j,k)) + beta*(Vn(i,j,k)-Vs(i,j,k)))
                #ENDIF
            
            end do
        end do
    end do

    
    sum_rhs = sum_rhs + sum(rhs)

end subroutine calculate_RHS_pressure

subroutine calculate_RHS_temperature(rhs, Ue, Uw, Vn, Vs, T, dx, dy, CpT, coeff)
    
    use precision_module, only: prcs_var
    use data_type_module, only: N, M, i_start, i_end, j_start, j_end
    use problem_module,   only: delta_time
    implicit none
    
    real(kind=prcs_var), intent(inout), dimension(M+2, N+2, 1)  ::  rhs, CpT
    real(kind=prcs_var), intent(in), dimension(M+2, N+2, 1)     ::  T, Ue, Uw, Vn, Vs
    real(kind=prcs_var), intent(in) ::  coeff, dx, dy
    integer         ::  i, j, k

    k = 1
    do j = j_start, j_end
        do i = i_start, i_end

            ! Convective Term
            CpT(i,j,k) = ( ( (T(i+1,j,k)+T(i,j,k))*Ue(i,j,k) - (T(i,j,k)+T(i-1,j,k))*Uw(i,j,k) )/dx &
                        +  ( (T(i,j+1,k)+T(i,j,k))*Vn(i,j,k) - (T(i,j,k)+T(i,j-1,k))*Vs(i,j,k) )/dy )*0.5_prcs_var


            ! Calculate the RHS
            rhs(i,j,k)  = coeff * ( - T(i,j,k)/delta_time + CpT(i,j,k))


        end do
    end do

end subroutine calculate_RHS_temperature