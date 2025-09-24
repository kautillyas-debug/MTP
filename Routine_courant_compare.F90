! --------------------------------------------------------------
!>        Check the Courant Number of the flow
! --------------------------------------------------------------
subroutine check_courant_number(tree)
    use precision_module, only: prcs_var
    use data_type_module
    use problem_module, only: courant_u_max, courant_v_max, delta_time
    implicit none

    type(t_tree), intent(inout) :: tree
    integer                     :: lvl, id, idx
    real(kind=prcs_var)         :: courant_u, courant_v

    courant_u = 0.0_prcs_var
    courant_v = 0.0_prcs_var

    do lvl = tree%no_levels, 1, -1
        do id = 1, no_leaf_idx_lvl(lvl)
            idx = tree%levels(lvl)%idx_leaf(id)

            courant_u = maxval(abs(tree%blocks(idx)%u(i_start:i_end, j_start:j_end, 1)))*delta_time/dx_lvl(lvl)
            courant_v = maxval(abs(tree%blocks(idx)%v(i_start:i_end, j_start:j_end, 1)))*delta_time/dy_lvl(lvl)

            if (courant_u > courant_u_max) courant_u_max = courant_u
            if (courant_v > courant_v_max) courant_v_max = courant_v

        end do
    end do
    

end subroutine check_courant_number


! --------------------------------------------------------------
!>        Compare current output to actual value (from fine mesh)
! --------------------------------------------------------------
subroutine compare_with_actual(tree) 

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use iso_fortran_env
    implicit none

    type(t_tree), intent(inout) :: tree

    character(len=10)   :: line
    character(len=100)  :: temp_line
    character(len=1)    :: path_separator
    character(len=100)  :: line_2

    integer :: pos_I, pos_J
    integer :: idx, M_act, N_act

    real(kind=prcs_var), allocatable :: x_act(:), y_act(:), u_act(:), v_act(:), p_act(:), T_act(:), var_act(:)
    real(kind=prcs_var) :: dx_act, dy_act, X, Y
    real(kind=prcs_var) :: x_temp_low, x_temp_high, y_temp_low, y_temp_high
    real(kind=prcs_var) :: a1, a2, a3, a4
    real(kind=prcs_var) :: var1, var2, var3, var4, var_interp
    real(kind=prcs_var) :: error_act

    logical :: check_flag, done_flag
    
    integer :: no_leaf_blocks

    integer :: lvl, id, i, j, k

    call get_environment_variable("PATH", temp_line)
    i = index(temp_line, '/')
    if (i>0) then
        path_separator = '/'
    else
        path_separator = '\'
    end if

    ! Reading the Fine Mesh Data (Actual Values)
    open(0, file = 'Verify'//path_separator//'Output_vel_Re400_Eps-6_original.csv')

        read(0,*)   !   variables = x,y,u,v,p,T
        read(0,'(A)')   line_2
        pos_I = index(line_2, "I =") + 3
        pos_J = index(line_2, "J =") + 3
        read(line_2(pos_I:), '(I10)') M_act
        read(line_2(pos_J:), '(I10)') N_act

        allocate(x_act(M_act*N_act))
        allocate(y_act(M_act*N_act))
        allocate(u_act(M_act*N_act))
        allocate(v_act(M_act*N_act))
        allocate(p_act(M_act*N_act))
        allocate(T_act(M_act*N_act))
        allocate(var_act(M_act*N_act))

        dx_act = 1.0_prcs_var / (M_act-2)
        dy_act = 1.0_prcs_var / (N_act-2)
        
        do i = 1, M_act*N_act

            read(0,*) x_act(i), y_act(i), u_act(i), v_act(i), p_act(i), T_act(i) 

        end do

        
    close(0)
    
    error_act = 0.0_prcs_var
    no_leaf_blocks = 0
    check_flag = .False.
    done_flag  = .True.

    ! Setting which Variable to be compared with
    var_act = v_act

    do lvl = 1, tree%no_levels
        do id = 1, size(tree%levels(lvl)%idx_leaf)
            idx = tree%levels(lvl)%idx_leaf(id)
            no_leaf_blocks = no_leaf_blocks + 1

            do j = j_start, j_end
                do i = i_start, i_end

                    X = tree%blocks(idx)%position(1)/domain(1) - real(i_mid-i+0.5)/domain(1)/2**(lvl-1)
                    Y = tree%blocks(idx)%position(2)/domain(2) - real(j_mid-j+0.5)/domain(2)/2**(lvl-1)

                    x_temp_low  = (int(X/dx_act+0.5)-0.5)*dx_act
                    y_temp_low  = (int(Y/dy_act+0.5)-0.5)*dy_act
                    x_temp_high = (int(X/dx_act+0.5)-0.5)*dx_act + dx_act
                    y_temp_high = (int(Y/dy_act+0.5)-0.5)*dy_act + dy_act

                    ! Finding the pressure values at the four corners
                    do k = 1, M_act*N_act
                        if (x_act(k)==x_temp_low) then
                            if (y_act(k)==y_temp_low)  var1 = var_act(k)
                            if (y_act(k)==y_temp_high) var3 = var_act(k)
                        end if
                        if (x_act(k)==x_temp_high) then
                            if (y_act(k)==y_temp_low)  var2 = var_act(k)
                            if (y_act(k)==y_temp_high) var4 = var_act(k)
                        end if
                    end do

                    a1 = (x_temp_high - X)*(y_temp_high - Y)/dx_act/dy_act
                    a2 = (-x_temp_low + X)*(y_temp_high - Y)/dx_act/dy_act
                    a3 = (x_temp_high - X)*(-y_temp_low + Y)/dx_act/dy_act
                    a4 = (-x_temp_low + X)*(-y_temp_low + Y)/dx_act/dy_act

                    var_interp = a1*var1 + a2*var2 + a3*var3 + a4*var4

                    ! if (X>0.5 .and. Y>0.5) then
                    !     check_flag = .True.

                    !     if (check_flag .and. done_flag) then
                    !         print *, var1, var2, var3, var4
                    !         print *, var_interp, tree%blocks(idx)%u(i,j,1)
                    !     end if

                    !     done_flag = .False.

                    ! end if

                    ! ------------------- Checking with the Calculated Variable ----------------
                    error_act = error_act + (tree%blocks(idx)%v(i,j,1) - var_interp)**2

                end do
            end do
        end do
    end do


    error_act = sqrt(error_act / (M*N*no_leaf_blocks))

    print *, "***************** ERROR *************************"    
    print *, "The error in Variable is : ", error_act
    ! print *, "***************** ERROR *************************"    

end subroutine compare_with_actual

