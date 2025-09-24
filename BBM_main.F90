
!        ____________ MAIN_PROGRAM ____________


program flow_solver_FSM

    #include "definitions.h"
    use iso_fortran_env
    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module

    implicit none

    type(t_tree) :: tree        !> Creating the tree which will store the data
    integer :: id, idx, lvl, no_time_steps
    integer :: iter_u, iter_v, iter_p, iter_T
    real :: cpu_start_time, cpu_end_time, cpu_prev_iter_time, cpu_iter_time

    integer :: i, j, k
    integer :: ibody

    integer :: values(8)
    character(len=20) :: formatted_date, formatted_time

    
    !> Read input files and initialise arrays
    call read_input_memory_allocate(tree)      ! Reading the input .dat file and allocate memory
    call array_initialise(tree)                ! Initialise the array to default values
    call init_multigrid_vars(tree)             ! Initialise multigrid variables and fill ghost cell values of sub-levels
    call init_body_vars(tree)                  ! Initialise the datastructure and variables for dealing with body in domain


    ! Get the current date and time
    call date_and_time(values=values)
    write(formatted_date, '(I2.2, "-", I2.2, "-", I4.4)') values(3), values(2), values(1)
    write(formatted_time, '(I2.2, ":", I2.2, ":", I2.2)') values(5), values(6), values(7)

    !> Print details of simulation:
        print *
        print *, '============================================================'
        print *
        print *, '--------------------- COMPILER INFORMATION ---------------------'
        write (*, '(A)')   '    Compiler version  :   '//trim(compiler_version())
        write (*, '(A)')   '    Compiler options  :   '//trim(compiler_options())
        print *, '----------------------------------------------------------------'
        print *
        print *, '--------------------- DATE-TIME INFORMATION --------------------'
        write(*, '(A)')   '    Date              :   '//trim(formatted_date)
        write(*, '(A)')   '    Time              :   '//trim(formatted_time)
        print *, '----------------------------------------------------------------'
        print *
        print *, '--------------------- FILE INFORMATION -------------------------'
        print *,           '   Mesh File         :   ', mesh_file
        print *,           '   Solver File       :   ', solver_file
        print *,           '   Body File         :   ', canonical_body_file
        print *, '----------------------------------------------------------------'
        print *    
        print *, '--------------------- PROBLEM INFORMATION ----------------------'
        print '(A, F8.0)', '    REYNOLDS NUMBER   :   ', Re
        print '(A, F8.3)', '    COURANT NUMBER    :   ', 1.0*delta_time/dx_lvl(tree%no_levels)  ! Assuming Umax = 1.0
        print '(A, F8.0)', '    1/dx at fine lvl  :   ', 1.0/dx_lvl(tree%no_levels)             ! Assuming Umax = 1.0
        print '(A, I8)',   '    Number of cells   :   ', tree%no_leaf_blocks * tree%n_cell(1) * tree%n_cell(2) 
        print '(A, I8)',   '    # solid cells     :   ', no_solid_cells
        print '(A, I8)',   '    # fluid cells     :   ', no_fluid_cells
        print *, '----------------------------------------------------------------'
    !

    !> Checking initial condition of CFL 
    if (max(u_left, u_right, u_top, u_bottom)*delta_time/dx_lvl(tree%no_levels) > 1.0) then
        print *
        print *, '-----------------------------------------------------------'
        print *, 'The Courant Number is above 1, Please reduce time step'
        print *, '    COURANT NUMBER       ::   ', 1.0*delta_time/dx_lvl(tree%no_levels)  ! Assuming Umax = 1.0
        print '(A, F10.6)', '    change delta_time to  ::    <', dx_lvl(tree%no_levels)/1.0_prcs_var  ! Assuming Umax = 1.0
        print *
        call EXIT()
    end if


    !> Initialise Time      
    time            = 0.0
    time_step       = 0.0
    no_time_steps   = int(max_time/delta_time)

    !> Initialise Iterations    
    iter_u          = 0
    iter_v          = 0
    iter_p          = 0
    iter_T          = 0

    !> Initialise Courant Numbers
    courant_u_max = 0.0_prcs_var
    courant_v_max = 0.0_prcs_var


    !> CPU's start time of execution
    call cpu_time(cpu_start_time)
    cpu_prev_iter_time = 0.0

    ! Applying Boundary Conditions
    call apply_boundary_condition_u_vel(tree)           ! Apply the boundary conditions on u - velocities
    call apply_boundary_condition_v_vel(tree)           ! Apply the boundary conditions on v - velocities
    call apply_boundary_condition_pressure(tree)        ! Apply the boundary conditions on pressure
    call apply_boundary_condition_temperature(tree)     ! Apply the boundary conditions on temperature
    call calculate_intermediate_face_velocities(tree)   ! Update face velocities with boundary conditions


    ! Update ghost cell of Variables around body
    do lvl = tree%no_levels, 1, -1
        do id = 1, tree%levels(lvl)%idx_body_leaf(0)
            idx = tree%levels(lvl)%idx_body_leaf(id)

            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                U_VELOCITY &
            )
            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                V_VELOCITY &
            )
            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%p, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                PRESSURE &
            )


            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                U_VELOCITY, &
                tree%blocks(idx)%dl_GC_IP & 
            )
            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                V_VELOCITY, &
                tree%blocks(idx)%dl_GC_IP & 
            )                
            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%p, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                PRESSURE, &
                tree%blocks(idx)%dl_GC_IP & 
            )

        end do
    end do


    !>                  -------------- FRACTIONAL STEP METHOD --------------

    ! do while(time < max_time)
    do while(time_step < no_time_steps)

        ! Printing the Heading for each Time Step
        #IFDEF DEBUG
            if (mod(time_step+1,STEP_SKIP)==0) then
                print *
                print *
                call cpu_time(cpu_end_time)
                cpu_iter_time = cpu_end_time - cpu_prev_iter_time
                cpu_iter_time = cpu_iter_time * (no_time_steps - time_step)
                write(formatted_time, '(I2.2, ":", I2.2, ":", I2.2)') &
                    int(cpu_iter_time / 3600.0), int(mod(cpu_iter_time, 3600.0) / 60.0), int(mod(cpu_iter_time, 60.0))
                cpu_prev_iter_time = cpu_end_time
                print *, "******************************************************************************"
                print '(A, I6, A, F10.6, A, F10.4, A, A)', '  Time-Step : ', time_step+1, '    Time : ', time+delta_time, &
                        '    CPU Time : ', cpu_end_time - cpu_start_time, '    Est. Time : ', trim(formatted_time)
                print *, "******************************************************************************"
            end if
        #ENDIF


        ! Calculate Courant Number (CFL)
        call check_courant_number(tree)
        #IFDEF DEBUG
            if (mod(time_step+1,STEP_SKIP)==0) then
                print *, '--------------------------------------------------------------------'
                print '(T24, A, F14.6)', 'Courant U max: ', courant_u_max
                print '(T24, A, F14.6)', 'Courant V max: ', courant_v_max
            end if
        #ENDIF



        !> ********** STEP 1 : Solve Modified Momentum Equations and Update Face Velocities **********

        !> X-Momentum Equation
        if (MME_SOLVER == GAUSS_SEIDEL) then
            call solve_x_modified_momentum_equation_GS(tree, iter_u)

        else if (MME_SOLVER == CONJUGATE_GRADIENT) then
            call solve_x_modified_momentum_equation_CG(tree, iter_u)
        
        end if
                    
        !> Y-Momentum Equation
        if (MME_SOLVER == GAUSS_SEIDEL) then
            call solve_y_modified_momentum_equation_GS(tree, iter_v)

        else if (MME_SOLVER == CONJUGATE_GRADIENT) then
            call solve_y_modified_momentum_equation_CG(tree, iter_v)
        
        end if

        !> Calculate Intermediate Face Velocities
        call calculate_intermediate_face_velocities(tree)
        !


        !> ********** STEP 2 : Solve Pressure Correction equation  **********
        if (PPE_SOLVER == GAUSS_SEIDEL) then
            call solve_pressure_poisson_equation_GS(tree, iter_p)   ! GAUSS SEIDEL
            ! call solve_pressure_poisson_equation_GS_singularity(tree, iter_p)   ! GAUSS SEIDEL

        else if (PPE_SOLVER == CONJUGATE_GRADIENT) then
            call solve_pressure_poisson_equation_CG(tree, iter_p)   ! CONJUGATE GRADIENT
            ! call solve_pressure_poisson_equation_CG_singularity(tree, iter_p)   ! CONJUGATE GRADIENT

        else if (PPE_SOLVER == MULTI_GRID) then
            call solve_pressure_poisson_equation_MG(tree, iter_p)   ! MULTI GRID
        
        end if



        !> ********** STEP 3 : Update cell and face centred velocities  **********
        call update_cell_face_velocities(tree)

        call enforce_no_slip_body_face(tree)

        !> Update boundary conditions after time-step calculations
        call apply_boundary_condition_u_vel(tree)
        call apply_boundary_condition_v_vel(tree)
        call apply_boundary_condition_pressure(tree, 1)
        call apply_boundary_condition_temperature(tree)

        ! Update the Ghost cells with updated cell centre values
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)
                call update_ghost_cell_u_vel(tree, idx)
                call update_ghost_cell_v_vel(tree, idx)
                call update_ghost_cell_pressure(tree, idx)
                call update_ghost_cell_temperature(tree, idx)
            end do
        end do

        ! Update ghost cell of Variables around body
        do lvl = tree%no_levels, 1, -1
            do id = 1, tree%levels(lvl)%idx_body_leaf(0)
            idx = tree%levels(lvl)%idx_body_leaf(id)

            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,U_VELOCITY+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                U_VELOCITY &
            )
            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,V_VELOCITY+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                V_VELOCITY &
            )
            call interpolate_IP_variable(&
                tree%blocks(idx)%IP_pos_val(:,PRESSURE+3), &
                tree%blocks(idx)%IP_interpol_cells, &
                tree%blocks(idx)%p, &
                tree%blocks(idx)%no_IB_ghost_cells(nbody), &
                PRESSURE &
            )


            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%u, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                U_VELOCITY, &
                tree%blocks(idx)%dl_GC_IP & 
            )
            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%v, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                V_VELOCITY, &
                tree%blocks(idx)%dl_GC_IP & 
            )                
            call update_IB_ghost_cells_variable(&
                tree%blocks(idx)%IP_pos_val, &
                tree%blocks(idx)%IB_ghost_cell, &
                tree%blocks(idx)%p, &
                tree%blocks(idx)%no_IB_ghost_cells, &
                PRESSURE, &
                tree%blocks(idx)%dl_GC_IP & 
            )

            end do
        end do



        !> ********** STEP 4 : Solve Energy Equations to calucalte temperatures **********
        ! call solve_energy_equation_GS(tree, iter_T)



        !> Incrementing time and time_step
        time = time + delta_time
        time_step = time_step + 1


        ! !> Printing output files at time step: 
        ! if (mod(time_step-1,1000)==0)  then
            ! call print_output_techplot(tree)       ! Techplot  : values are interpolated to coarsest cells of tree blocks
            ! call print_output_matplotlib(tree)     ! Matplotlib: values are prints values in cell center (no interpolation)
            ! call print_output_vtk(tree)            ! VTK       : values are prints values in nodes of cells at final blocks(interpolation performed)
        ! end if

        
        ! call print_screen_real(tree%blocks(1)%rhs, 'rH', 1, M+2, 1, N+2)
        ! call print_screen_real(tree%blocks(1)%u, 'uc', 1, M+2, 1, N+2)
        ! call print_screen_real(tree%blocks(1)%Ue, 'ue', 1, M+2, 1, N+2)
        ! call print_screen_real(tree%blocks(1)%Vs, 'vs', 1, M+2, 1, N+2)
        ! call print_screen_real(tree%blocks(1)%p, 'pr', 1, M+2, 1, N+2)


        ! Printing the Iterations and Divergence of Flow
        #IFDEF DEBUG
            if (mod(time_step+1,STEP_SKIP)==0) then
                
                ! Calculate RHS of pressure poisson rhs after time-step calculation
                rhs_sum_prev = rhs_sum
                rhs_sum = 0.0_prcs_var
                do lvl = tree%no_levels, 1, -1
                    do id = 1, no_leaf_idx_lvl(lvl)
                        idx = tree%levels(lvl)%idx_leaf(id)
                        call calculate_rhs_pressure(tree%blocks(idx)%rhs, &
                                                    tree%blocks(idx)%Ue, tree%blocks(idx)%Uw, &
                                                    tree%blocks(idx)%Vn, tree%blocks(idx)%Vs, &
                                                    coeff_p(lvl), rhs_sum, tree%blocks(idx)%i_range)                 
                    end do
                end do
           
                print *, '        Q_out *beta/dt (before PPE): ', rhs_sum_prev
                print *, '        Q_out *beta/dt (after  PPE): ', rhs_sum
                print *, '--------------------------------------------------------------------'
                print '(T6, A, T19, A, T38, A, T52, A)', 'Equation', 'Iteration', 'Residual', 'Solver'
                print *, '--------------------------------------------------------------------'
                print '(5X, A, T22, I6, T32, E14.4, T56, A)', 'x-mm', iter_u, residual_u, mme_solver_text
                print '(5X, A, T22, I6, T32, E14.4, T56, A)', 'y-mm', iter_v, residual_v, mme_solver_text
                print '(5X, A, T22, I6, T32, E14.4, T56, A)', 'pressure', iter_p, residual_p, ppe_solver_text
                print '(5X, A, T22, I6, T32, E14.4)', 'energy', iter_T, residual_T
                print *, '--------------------------------------------------------------------'
            end if

        #ENDIF


    end do

    ! IF sq body present at (.35, .35) to (.65, .65) these can print values at bottom left corner

    ! FOR 200x200 grid
    ! call print_screen_real(tree%blocks(1)%u, 'uc', 70, 80, 70, 80)
    ! call print_screen_real(tree%blocks(1)%v, 'vc', 70, 80, 70, 80)
    ! call print_screen_real(tree%blocks(1)%p, 'pr', 70, 80, 70, 80)

    ! ! For 100x100 grid
    ! call print_screen_real(tree%blocks(1)%u, 'uc', 35, 45, 35, 45)
    ! call print_screen_real(tree%blocks(1)%v, 'vc', 35, 45, 35, 45)
    ! call print_screen_real(tree%blocks(1)%p, 'pr', 35, 45, 35, 45)

    !> Store the end time of program and calculate time taken
    call cpu_time(cpu_end_time)
    print *
    print *, 'Simulation Ended ...'
    print '(A, F12.6, A)', 'Time Elapsed: ', cpu_end_time - cpu_start_time, ' Seconds'
    print *, '============================================================'
    print *


    !> Printing output files post simulation: 
    ! call print_output_techplot(tree)       ! Techplot  : values are interpolated to coarsest cells of tree blocks
    call print_output_matplotlib(tree)     ! Matplotlib: values are prints values in cell center (no interpolation)
    ! call print_output_vtk(tree)            ! VTK       : values are prints values in nodes of cells at final blocks(interpolation performed)

    ! !> ----- For Grid Convergence Study ------
    ! !> Compare with actual values
    ! call compare_with_actual(tree)

        
end program flow_solver_FSM