
    #include "definitions.h"

! ------------- Print value of variable on screen ------------------
subroutine print_screen_real(variable, text, i_start, i_end, j_start, j_end)
    use precision_module, only: prcs_var
    use data_type_module, only: M, N, ngl

    integer, intent(in) :: i_start, i_end, j_start, j_end
    real(kind=prcs_var), intent(in), dimension(i_start:i_end,j_start:j_end,1) :: variable
    character(len=2), intent(in) :: text
    integer :: i, j, i_start_print, j_start_print, i_end_print, j_end_print
    character(len=20) :: line_char
    write(line_char, '(I0)') j_end - j_start
    line_char = '('//trim(adjustl(line_char))//'F13.9)'

    print *
    print *, 'Printing values of ', text

    i_start_print = i_start
    i_end_print   = i_end
    j_start_print = j_start
    j_end_print   = j_end

    ! i_start_print = 70
    ! i_end_print   = 90
    ! j_start_print = 35
    ! j_end_print   = 45
    
    do j = j_end_print, j_start_print, -1
        write(*, '(I2, A2)', advance="no") j, '|'
        do i = i_start_print, i_end_print
            write(*, line_char, advance="no") variable(i, j, 1)
        end do
        print *
    end do
    write(*, *) "  | _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ "
    write(*, '(A4)', advance="no") "   "
    do i = i_start_print, i_end_print
        write(*, line_char, advance="no") real(i, kind=prcs_var)
    end do

    print *
    print *
end subroutine print_screen_real

subroutine print_screen_int(variable, text, i_start, i_end, j_start, j_end)

    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer, intent(in), dimension(i_start:i_end,j_start:j_end,1) :: variable
    character(len=2), intent(in) :: text
    integer :: i, j, idx
    character(len=20) :: line_char
    write(line_char, '(I0)') j_end - j_start +1
    line_char = '('//trim(adjustl(line_char))//'I3)'

    print *
    print *, 'Printing values of ', text
    idx = 1
    
    do j = j_end, j_start, -1
        write(*, '(I2, A2)', advance="no") j, '|'
        do i = i_start, i_end
            write(*, line_char, advance="no") variable(i, j, 1)
        end do
        print *
    end do
    write(*, *) "  | _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ "
    write(*, '(A4)', advance="no") "   "
    do i = i_start, i_end
        write(*, line_char, advance="no") i
    end do

    print *
    print *
end subroutine print_screen_int

subroutine print_screen_int_1(variable, text, i_start, i_end, j_start, j_end)

    use data_type_module, only: M, N 

    integer, intent(in) :: i_start, i_end, j_start, j_end
    integer*1, intent(in), dimension(i_start:i_end,j_start:j_end,1) :: variable
    character(len=2), intent(in) :: text
    integer :: i, j, idx
    character(len=20) :: line_char
    write(line_char, '(I0)') j_end - j_start +1
    line_char = '('//trim(adjustl(line_char))//'I3)'

    print *
    print *, 'Printing values of ', text
    idx = 1
    
    do j = j_end, j_start, -1
        write(*, '(I2, A2)', advance="no") j, '|'
        do i = i_start, i_end
            write(*, line_char, advance="no") variable(i, j, 1)
        end do
        print *
    end do
    write(*, *) "  | _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ "
    write(*, '(A4)', advance="no") "   "
    do i = i_start, i_end
        write(*, line_char, advance="no") i
    end do

    print *
    print *
end subroutine print_screen_int_1


! --------------------------------------------------------------
!         Print Output Files 
! Techplot  : values are interpolated to coarsest cells of tree blocks
! Matplotlib: values are prints values in cell center (no interpolation)
! VTK       : values are prints values in nodes of cells at final blocks(interpolation performed)
! --------------------------------------------------------------
subroutine print_output_techplot(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none
    type(t_tree), intent(inout) :: tree
    ! real(kind=prcs_var), dimension(M+1,N+1) ::  u_node, v_node, vorticity
    real(kind=prcs_var) ::  x, y, skip, ox, oy
    integer :: i, j, lvl, id, idx ! k, 
    character(len=100)  :: path, filename
    character(len=1)    :: path_separator
    character(len=9)    :: epsilon_temp, Reynolds_temp
    character(len=:), allocatable :: epsilon, Reynolds

    call get_environment_variable("PATH", path)
    i = index(path, '/')
    if (i>0) then
        path_separator = '/'
    else
        path_separator = '\'
    end if

    write(Reynolds_temp, '(I4)') int(Re)
    i = len(trim(adjustl(Reynolds_temp)))
    allocate(character(len=i) :: Reynolds)
    Reynolds = trim(adjustl(Reynolds_temp))

    write(epsilon_temp, '(I6)') int((time_step-1))
    i = len(trim(adjustl(epsilon_temp)))
    allocate(character(len=i) :: epsilon)
    epsilon = trim(adjustl(epsilon_temp))

    filename = '.'//path_separator//'Verify'//path_separator//'Time_Step_Techplot'//&
        path_separator//'Output_Re'//Reynolds//'-Time-'//epsilon//'.plt'

    print *
    print *, 'Printing ', filename
    open(unit=4, file=filename, status='replace')
    

        write(4, '(A)') "variables = x,y,u,v,p,T"
        write(4,'(A, I6, A, I6)') 'Zone I =',M+1, ', J =',N+1

        do lvl = 1, tree%no_levels
            skip = 2**(lvl-1)
            do id = 1, size(tree%levels(lvl)%idx_leaf)
                idx = tree%levels(lvl)%idx_leaf(id)

                ! Write all cell centred values (including ghost cells)
                    ! ox = tree%blocks(idx)%position(1) - real(M+1)*dx_lvl(lvl)/2.0_prcs_var
                    ! oy = tree%blocks(idx)%position(2) - real(N+1)*dy_lvl(lvl)/2.0_prcs_var
                    ! do j = j_start-ngl, j_end+ngl
                    !     do i = i_start-ngl, i_end+ngl

                    !         x = ox + real(i-1)*dx_lvl(lvl)
                    !         y = oy + real(j-1)*dy_lvl(lvl)

                    !         if (PPE_SOLVER==MULTI_GRID) then
                    !             write(4, '(F, F, F, F, F, F)') &
                    !                 x, y, &
                    !                 tree%blocks(idx)%u(i,j,1), &        !> U-velocity
                    !                 tree%blocks(idx)%v(i,j,1), &        ! V-velocity
                    !                 tree%blocks(idx)%multi_grid(1)%p(i,j,1), &        ! pressure
                    !                 tree%blocks(idx)%t(i,j,1)                ! temperature
                    !         else 
                    !             write(4, '(F, F, F, F, F, F)') &
                    !                 x, y, &
                    !                 tree%blocks(idx)%u(i,j,1), &        !> U-velocity
                    !                 tree%blocks(idx)%v(i,j,1), &        ! V-velocity
                    !                 tree%blocks(idx)%p(i,j,1), &        ! pressure
                    !                 tree%blocks(idx)%t(i,j,1)           ! temperature

                    !         end if

                    !     end do
                    ! end do
                !

                ! Write values at the nodes
                    ox = tree%blocks(idx)%position(1) - real(M)*dx_lvl(lvl)/2.0_prcs_var
                    oy = tree%blocks(idx)%position(2) - real(N)*dy_lvl(lvl)/2.0_prcs_var
                    do j = j_start-ngl, j_end
                        do i = i_start-ngl, i_end

                            x = ox + real(i-1)*dx_lvl(lvl)
                            y = oy + real(j-1)*dy_lvl(lvl)

                            if (PPE_SOLVER==MULTI_GRID) then

                                write(4, '(F, F, F, F, F, F)') &
                                    x, y, &
                                        ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                        + tree%blocks(idx)%u(i+1,j,1) & 
                                        + tree%blocks(idx)%u(i,j+1,1) & 
                                        + tree%blocks(idx)%u(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                        + tree%blocks(idx)%v(i+1,j,1) & 
                                        + tree%blocks(idx)%v(i,j+1,1) & 
                                        + tree%blocks(idx)%v(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%multi_grid(1)%p(i,j,1) &        ! pressure
                                        + tree%blocks(idx)%multi_grid(1)%p(i+1,j,1) & 
                                        + tree%blocks(idx)%multi_grid(1)%p(i,j+1,1) & 
                                        + tree%blocks(idx)%multi_grid(1)%p(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                        + tree%blocks(idx)%T(i+1,j,1) & 
                                        + tree%blocks(idx)%T(i,j+1,1) & 
                                        + tree%blocks(idx)%T(i+1,j+1,1) &
                                        )/4.0

                            else 
                                write(4, '(F, F, F, F, F, F)') &
                                    x, y, &
                                        ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                        + tree%blocks(idx)%u(i+1,j,1) & 
                                        + tree%blocks(idx)%u(i,j+1,1) & 
                                        + tree%blocks(idx)%u(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                        + tree%blocks(idx)%v(i+1,j,1) & 
                                        + tree%blocks(idx)%v(i,j+1,1) & 
                                        + tree%blocks(idx)%v(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%p(i,j,1) &        ! pressure
                                        + tree%blocks(idx)%p(i+1,j,1) & 
                                        + tree%blocks(idx)%p(i,j+1,1) & 
                                        + tree%blocks(idx)%p(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                        + tree%blocks(idx)%T(i+1,j,1) & 
                                        + tree%blocks(idx)%T(i,j+1,1) & 
                                        + tree%blocks(idx)%T(i+1,j+1,1) &
                                        )/4.0
                            end if

                        end do
                    end do
                !

            end do
        end do

    close(4)


end subroutine print_output_techplot

subroutine print_output_matplotlib(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none
    type(t_tree), intent(inout) :: tree
    ! real(kind=prcs_var), dimension(M+1,N+1) ::  u_node, v_node, vorticity
    real(kind=prcs_var) ::  x, y, ox, oy, fluid
    integer :: i, j, lvl, id, idx, skip ! k, 
    character(len=100)  :: path, filename
    character(len=1)    :: path_separator
    character(len=8)    :: epsilon_temp, Reynolds_temp
    character(len=:), allocatable :: epsilon, Reynolds

    call get_environment_variable("PATH", path)
    i = index(path, '/')
    if (i>0) then
        path_separator = '/'
    else
        path_separator = '\'
    end if

    write(Reynolds_temp, '(I4)') int(Re)
    i = len(trim(adjustl(Reynolds_temp)))
    allocate(character(len=i) :: Reynolds)
    Reynolds = trim(adjustl(Reynolds_temp))

    write(epsilon_temp, '(I6)') int((time_step-1))
    i = len(trim(adjustl(epsilon_temp)))
    allocate(character(len=i) :: epsilon)
    epsilon = trim(adjustl(epsilon_temp))

    filename = '.'//path_separator//'Verify'//path_separator//'Time_Step_Techplot'//&
        path_separator//'Output_Re'//Reynolds//'-Time-'//epsilon//'_original.csv'

    print *
    print *, 'Printing ', filename
    open(unit=4, file=filename, status='replace')

        write(4, '(A)') "variables = x,y,u,v,p,T"
        write(4,'(A, I6, A, I6)') 'Zone I =',i_end+ngl, ', J =',j_end+ngl

        do lvl = 1, tree%no_levels
            skip = 2**(lvl-1)
            do id = 1, size(tree%levels(lvl)%idx_leaf)
                idx = tree%levels(lvl)%idx_leaf(id)

                ! Write all cell centred values (including ghost cells)
                    ! ox = tree%blocks(idx)%position(1) - real(M+1)*dx_lvl(lvl)/2.0_prcs_var
                    ! oy = tree%blocks(idx)%position(2) - real(N+1)*dy_lvl(lvl)/2.0_prcs_var
                    ! do j = j_start-ngl, j_end+ngl
                    !     do i = i_start-ngl, i_end+ngl

                    !         x = ox + real(i-1)*dx_lvl(lvl)
                    !         y = oy + real(j-1)*dy_lvl(lvl)

                    !         if (PPE_SOLVER==MULTI_GRID) then
                    !             write(4, '(F, F, F, F, F, F)') &
                    !                 x, y, &
                    !                 tree%blocks(idx)%u(i,j,1), &        !> U-velocity
                    !                 tree%blocks(idx)%v(i,j,1), &        ! V-velocity
                    !                 tree%blocks(idx)%multi_grid(1)%p(i,j,1), &        ! pressure
                    !                 tree%blocks(idx)%t(i,j,1)                ! temperature
                    !         else 
                    !             write(4, '(F, F, F, F, F, F)') &
                    !                 x, y, &
                    !                 tree%blocks(idx)%u(i,j,1), &        !> U-velocity
                    !                 tree%blocks(idx)%v(i,j,1), &        ! V-velocity
                    !                 tree%blocks(idx)%p(i,j,1), &        ! pressure
                    !                 tree%blocks(idx)%t(i,j,1)           ! temperature

                    !         end if

                    !     end do
                    ! end do
                !

                ! Write values at the nodes
                    ox = tree%blocks(idx)%position(1) - real(M)*dx_lvl(lvl)/2.0_prcs_var
                    oy = tree%blocks(idx)%position(2) - real(N)*dy_lvl(lvl)/2.0_prcs_var
                    do j = j_start-ngl, j_end
                        do i = i_start-ngl, i_end

                            x = ox + real(i-1)*dx_lvl(lvl)
                            y = oy + real(j-1)*dy_lvl(lvl)

                            if (PPE_SOLVER==MULTI_GRID) then

                                write(4, '(F, F, F, F, F, F)') &
                                    x, y, &
                                        ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                        + tree%blocks(idx)%u(i+1,j,1) & 
                                        + tree%blocks(idx)%u(i,j+1,1) & 
                                        + tree%blocks(idx)%u(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                        + tree%blocks(idx)%v(i+1,j,1) & 
                                        + tree%blocks(idx)%v(i,j+1,1) & 
                                        + tree%blocks(idx)%v(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%multi_grid(1)%p(i,j,1) &        ! pressure
                                        + tree%blocks(idx)%multi_grid(1)%p(i+1,j,1) & 
                                        + tree%blocks(idx)%multi_grid(1)%p(i,j+1,1) & 
                                        + tree%blocks(idx)%multi_grid(1)%p(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                        + tree%blocks(idx)%T(i+1,j,1) & 
                                        + tree%blocks(idx)%T(i,j+1,1) & 
                                        + tree%blocks(idx)%T(i+1,j+1,1) &
                                        )/4.0

                            else 
                                write(4, '(F, F, F, F, F, F)') &
                                    x, y, &
                                        ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                        + tree%blocks(idx)%u(i+1,j,1) & 
                                        + tree%blocks(idx)%u(i,j+1,1) & 
                                        + tree%blocks(idx)%u(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                        + tree%blocks(idx)%v(i+1,j,1) & 
                                        + tree%blocks(idx)%v(i,j+1,1) & 
                                        + tree%blocks(idx)%v(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%p(i,j,1) &        ! pressure
                                        + tree%blocks(idx)%p(i+1,j,1) & 
                                        + tree%blocks(idx)%p(i,j+1,1) & 
                                        + tree%blocks(idx)%p(i+1,j+1,1) &
                                        )/4.0, &
                                        ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                        + tree%blocks(idx)%T(i+1,j,1) & 
                                        + tree%blocks(idx)%T(i,j+1,1) & 
                                        + tree%blocks(idx)%T(i+1,j+1,1) &
                                        )/4.0
                            end if

                        end do
                    end do
                !

            end do
        end do

    close(4)


end subroutine print_output_matplotlib

subroutine print_output_vtk(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    implicit none
    type(t_tree), intent(inout) :: tree
    ! real(kind=prcs_var), dimension(M+1,N+1) ::  u_node, v_node, vorticity
    real(kind=prcs_var) ::  x, y, ox, oy
    integer :: i, j, lvl, id, idx, skip ! k, 
    character(len=100)  :: path, filename
    character(len=1)    :: path_separator
    character(len=9)    :: epsilon_temp, Reynolds_temp
    character(len=:), allocatable :: epsilon, Reynolds
    integer :: idx_count

    call get_environment_variable("PATH", path)
    i = index(path, '/')
    if (i>0) then
        path_separator = '/'
    else
        path_separator = '\'
    end if

    write(Reynolds_temp, '(I4)') int(Re)
    i = len(trim(adjustl(Reynolds_temp)))
    allocate(character(len=i) :: Reynolds)
    Reynolds = trim(adjustl(Reynolds_temp))

    ! write(epsilon_temp, '(I3)') int(log10(epsilon_p))
    ! i = len(trim(adjustl(epsilon_temp)))
    ! allocate(character(len=i) :: epsilon)
    ! epsilon = trim(adjustl(epsilon_temp))

    ! filename = '.'//path_separator//'Verify'//path_separator//&
    ! 'Output_vel_Re'//Reynolds//'_Eps'//epsilon//'_original_vtk.csv'

    write(epsilon_temp, '(I6)') int((time_step-1))
    i = len(trim(adjustl(epsilon_temp)))
    allocate(character(len=i) :: epsilon)
    epsilon = trim(adjustl(epsilon_temp))

    filename = '.'//path_separator//'Verify'//path_separator//'Time_Step_Techplot'//&
        path_separator//'Output_Re'//Reynolds//'-Time-'//epsilon//'_vtk.csv'


    idx_count = 0
    do lvl = 1, tree%no_levels
        idx_count = idx_count + size(tree%levels(lvl)%idx_leaf)
    end do


    print *
    print *, 'Printing ', filename
    open(unit=4, file=filename, status='replace')

        write(4, '(A)') "variables = i,j,u,v,p,T"
        write(4,'(A)') 'M, N, total_blocks'
        write(4,'(I6, I6)') M+1, N+1, idx_count
        write(4, *)

        do lvl = 1, tree%no_levels
            skip = 2**(lvl-1)
            do id = 1, size(tree%levels(lvl)%idx_leaf)
                idx = tree%levels(lvl)%idx_leaf(id)


                ! Origin point, bottom left corner of the bottom left cell within domain
                ox = tree%blocks(idx)%position(1) - real(i_mid-1)*cell_dim(1)/skip
                oy = tree%blocks(idx)%position(2) - real(j_mid-1)*cell_dim(2)/skip

                write(4, *) 'Block:', idx
                write(4, *) 'y-Bottom:', oy
                write(4, *) 'y-Top:', oy + cell_dim(2)*N/skip
                write(4, *) 'x-Left:', ox
                write(4, *) 'x-Right:', ox + cell_dim(1)*M/skip

                ! do j = j_start-ngl, j_end, skip
                !     do i = i_start-ngl, i_end, skip
                do j = j_start-ngl, j_end
                    do i = i_start-ngl, i_end
                        
                        ! x = ox + real(i)/M/skip
                        ! y = oy + real(j)/M/skip
                        if (PPE_SOLVER==MULTI_GRID) then
                            write(4, '(I, I, F, F, F, F)') &
                                i, j, &
                                    ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                    + tree%blocks(idx)%u(i+1,j,1) & 
                                    + tree%blocks(idx)%u(i,j+1,1) & 
                                    + tree%blocks(idx)%u(i+1,j+1,1) &
                                    )/4.0, &
                                    ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                    + tree%blocks(idx)%v(i+1,j,1) & 
                                    + tree%blocks(idx)%v(i,j+1,1) & 
                                    + tree%blocks(idx)%v(i+1,j+1,1) &
                                    )/4.0, &
                                    
                                    ( tree%blocks(idx)%multi_grid(1)%p(i,j,1) &        ! pressure
                                    + tree%blocks(idx)%multi_grid(1)%p(i+1,j,1) & 
                                    + tree%blocks(idx)%multi_grid(1)%p(i,j+1,1) & 
                                    + tree%blocks(idx)%multi_grid(1)%p(i+1,j+1,1) &
                                    )/4.0, &
                                    
                                    ( tree%blocks(idx)%p(i,j,1) &        ! pressure
                                    + tree%blocks(idx)%p(i+1,j,1) & 
                                    + tree%blocks(idx)%p(i,j+1,1) & 
                                    + tree%blocks(idx)%p(i+1,j+1,1) &
                                    )/4.0, &
                                    ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                    + tree%blocks(idx)%T(i+1,j,1) & 
                                    + tree%blocks(idx)%T(i,j+1,1) & 
                                    + tree%blocks(idx)%T(i+1,j+1,1) &
                                    )/4.0
                        else
                            write(4, '(I, I, F, F, F, F)') &
                                i, j, &
                                    ( tree%blocks(idx)%u(i,j,1) &        !> U-velocity
                                    + tree%blocks(idx)%u(i+1,j,1) & 
                                    + tree%blocks(idx)%u(i,j+1,1) & 
                                    + tree%blocks(idx)%u(i+1,j+1,1) &
                                    )/4.0, &
                                    ( tree%blocks(idx)%v(i,j,1) &        ! V-velocity
                                    + tree%blocks(idx)%v(i+1,j,1) & 
                                    + tree%blocks(idx)%v(i,j+1,1) & 
                                    + tree%blocks(idx)%v(i+1,j+1,1) &
                                    )/4.0, &
                                    
                                    ( tree%blocks(idx)%multi_grid(1)%p(i,j,1) &        ! pressure
                                    + tree%blocks(idx)%multi_grid(1)%p(i+1,j,1) & 
                                    + tree%blocks(idx)%multi_grid(1)%p(i,j+1,1) & 
                                    + tree%blocks(idx)%multi_grid(1)%p(i+1,j+1,1) &
                                    )/4.0, &
                                    
                                    ( tree%blocks(idx)%p(i,j,1) &        ! pressure
                                    + tree%blocks(idx)%p(i+1,j,1) & 
                                    + tree%blocks(idx)%p(i,j+1,1) & 
                                    + tree%blocks(idx)%p(i+1,j+1,1) &
                                    )/4.0, &
                                    ( tree%blocks(idx)%T(i,j,1) &        ! Temperature
                                    + tree%blocks(idx)%T(i+1,j,1) & 
                                    + tree%blocks(idx)%T(i,j+1,1) & 
                                    + tree%blocks(idx)%T(i+1,j+1,1) &
                                    )/4.0
                        end if
                    end do
                end do
            end do
        end do

    close(4)


end subroutine print_output_vtk

subroutine print_output_iblank(tree, idx)
    use precision_module, only: prcs_var
    use data_type_module
    use problem_module

    implicit none

    type(t_tree), intent(in) :: tree
    integer, intent(in) :: idx
    real(kind=prcs_var) :: x, y
    integer :: i, j, i1, l, k
    character(len=100) :: filename
    character(len=1) :: path_separator
    character(len=2) :: idx_text
    integer          :: i_range_temp(M+2, N+2)

    ! Set the path separator based on platform
    path_separator = '/'
    
    ! Format idx as 2-digit, zero-padded
    write(idx_text, '(I2.2)') idx

    ! Construct the filename
    filename = '.' // path_separator // 'Verify' // path_separator //'iblank_check' // path_separator // &
               'Output_Iblank_Idx-' // trim(idx_text) // '.csv'

    print *, "Printing to file:", trim(filename)


    ! calculating i_range_temp values
    i_range_temp = 0

    k = 1
    do j = 2, N+1
        ! do i = i_start, i_end
        do i1 = 1, tree%blocks(idx)%i_range(0,j,k)
            l = 2*i1 - 1
            do i = tree%blocks(idx)%i_range(l,j,k), tree%blocks(idx)%i_range(l+1,j,k)
                
                i_range_temp(i, j) = 1

            end do
        end do
    end do


    ! Open the file
    open(unit=4, file=filename, status='replace', action='write')

    ! Write the file header
    write(4, '(A)') "variables = x, y, i_blank, i_range"
    write(4, '(A, I6, A, I6)') 'Zone I =', M + 2, ', J =', N + 2

    ! Loop through indices and compute values
    do j = j_start - ngl, j_end + ngl
        do i = i_start - ngl, i_end + ngl
    ! do i = i_start - ngl, i_end + ngl
    !     do j = j_start - ngl, j_end + ngl
            x = tree%blocks(idx)%position(1) / domain(1) - real(i_mid - i) / M
            y = tree%blocks(idx)%position(2) / domain(2) - real(j_mid - j) / N
            write(4, '(F10.5, 1X, F10.5, 1X, I2.2, 1X, I2.2)') x, y, &
                    tree%blocks(idx)%i_blank(i, j, 1, 1), i_range_temp(i,j)
        end do
    end do

    ! Close the file
    close(4)
end subroutine print_output_iblank
