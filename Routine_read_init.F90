
    #include "definitions.h"

! --------------------------------------------------------------
!>        Read Input File and Allocate Size for Tree
! --------------------------------------------------------------
subroutine read_input_memory_allocate(tree) 

    use precision_module, only: prcs_var
    use data_type_module
    use problem_module
    use body_module
    use iso_fortran_env

    implicit none
    type(t_tree), intent(inout) :: tree
    integer :: i, j, temp_len, block_count
    integer :: id, idx, lvl
    integer :: ibody
    character(len=10)   :: line
    character(len=100)  :: temp_line
    character(len=1)    :: path_separator
    character*100 :: body_file_name
    real(kind=prcs_var) :: temp_num, temp_den

    call get_environment_variable("PATH", temp_line)
    i = index(temp_line, '/')
    path_separator = '\'
    if (i>0) path_separator = '/'

    open(0, file = '..'//path_separator//'Input_Files.txt')

        read(0,*)   !-----------------------------------------------------------------------
        read(0,*)   !		FSI Flow Solver (IIT-Guwahati)
        read(0,*)   !		Input files for the solver 
        read(0,*)   !-----------------------------------------------------------------------
        read(0,*)   !		        
        read(0,*)   ! Mesh Information (Block_Structured_Mesh.dat)
        read(0,'(A)')           temp_line
        temp_len = len_trim(temp_line)
        allocate(character(len=temp_len) :: mesh_file)
        mesh_file = temp_line(1:temp_len)
        read(0,*)   !

        read(0,*)   !  Solver Input (Input_Solver.txt)
        read(0,'(A)')           temp_line
        temp_len = len_trim(temp_line)
        allocate(character(len=temp_len) :: solver_file)
        solver_file = temp_line(1:temp_len)
        read(0,*)   !

        read(0,*)   !  Canonical Body Data (Canonical_Body_Data.txt)
        read(0,'(A)')           temp_line
        temp_len = len_trim(temp_line)
        allocate(character(len=temp_len) :: canonical_body_file)
        canonical_body_file = temp_line(1:temp_len)

    close(0)

    
    ! Reading the Mesh File
    open(1, file = '..'//path_separator//'Mesh_Generator'//path_separator//'Meshes'//path_separator//mesh_file)
    
        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !               BLOCK STRUCTURED MESH                                 
        read(1,*)    !
        read(1,*)    !   Date: xx-xx-xxxx		Time: xx:xx:xx
        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !Input Details 
        read(1,*)    !      NDIM:
        read(1,*)       NDIM

        ! Allocate size with NDIM
        allocate(domain(NDIM))
        allocate(cell_dim(NDIM))
        allocate(tree%n_cell(NDIM))

        read(1,*)    !      Domain :
        read(1,*)       domain
        read(1,*)    !      cell_dim :
        read(1,*)       cell_dim
        read(1,*)    !      n_cell: :
        read(1,*)       tree%n_cell
        read(1,*)    !      no_bound: :
        read(1,*)       tree%no_bound


        ngl = 1        ! Number of Ghost Layers used
        M = tree%n_cell(1)
        N = tree%n_cell(2)

        i_start = ngl + 1
        i_end   = ngl + M
        i_mid   = int(M/2+ngl)

        j_start = ngl + 1
        j_end   = ngl + N
        j_mid   = int(N/2+ngl)



        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !Tree Details 
        read(1,*)    !      no_levels:
        read(1,*)       tree%no_levels
        read(1,*)    !      no_blocks:
        read(1,*)       tree%no_blocks

        ! Allocate level and block
        allocate(tree%levels(tree%no_levels))
        allocate(tree%blocks(tree%no_blocks))
        allocate(alpha(tree%no_levels))
        allocate(dx_lvl(tree%no_levels))
        allocate(dy_lvl(tree%no_levels))
        allocate(dx_sq_lvl(tree%no_levels))
        allocate(dy_sq_lvl(tree%no_levels))
        allocate(no_leaf_idx_lvl(tree%no_levels))
        allocate(gamma_u(tree%no_levels))
        allocate(gamma_T(tree%no_levels))
        allocate(coeff_p(tree%no_levels))
        allocate(coeff_u(tree%no_levels))
        allocate(coeff_T(tree%no_levels))

        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !Level Details
        do i = 1, tree%no_levels
            read(1,*)   !   Level_i
            read(1,*)   !   idx_all
            read(1,*)   temp_len
            allocate(tree%levels(i)%idx_all(temp_len))
            read (1,*)  tree%levels(i)%idx_all
            read(1,*)   !   idx_leaf
            read(1,*)   no_leaf_idx_lvl(i)
            allocate(tree%levels(i)%idx_leaf(no_leaf_idx_lvl(i)))
            read (1,*)  tree%levels(i)%idx_leaf
        end do
        
        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !Boundary Blocks Details   
        read(1,*)    !            Left boundary:
        read(1,*)   temp_len
        allocate(tree%left_blocks(temp_len))
        read(1,*)  tree%left_blocks
        read(1,*)    !            Right boundary:
        read(1,*)   temp_len
        allocate(tree%right_blocks(temp_len))
        read(1,*)  tree%right_blocks
        read(1,*)    !            Bottom boundary:
        read(1,*)   temp_len
        allocate(tree%bottom_blocks(temp_len))
        read(1,*)  tree%bottom_blocks
        read(1,*)    !            Top boundary:
        read(1,*)   temp_len
        allocate(tree%top_blocks(temp_len))
        read(1,*)  tree%top_blocks

        read(1,*)    !---------------------------------------------------------------------- 
        read(1,*)    !Block Details

        block_count = 0
        do i = 1, tree%no_blocks

            read(1,*)    !      idx:
            read(1,*)       tree%blocks(i)%idx
            read(1,*)    !      lvl:
            read(1,*)       tree%blocks(i)%lvl
            read(1,*)    !      is_refined:

            read(1, '(A)') line
            if (trim(adjustl(line)) == "	True") then
                tree%blocks(i)%is_refined = .true.
            else
                tree%blocks(i)%is_refined = .false.
            end if

            if (.not. (tree%blocks(i)%is_refined)) then
                block_count = block_count + 1
            end if

            read(1,*)    !      position:
            allocate(tree%blocks(i)%position(NDIM))
            read(1,*)       tree%blocks(i)%position
            read(1,*)    !      ptr_parent:
            read(1,*)       tree%blocks(i)%ptr_parent
            read(1,*)    !      ptr_children:
            if(tree%blocks(i)%is_refined) then
                allocate(tree%blocks(i)%ptr_children(2**NDIM))
                read(1,*)   tree%blocks(i)%ptr_children
            else
                read(1,*)
            end if
            read(1,*)    !      ptr_nb_blocks:
            allocate(tree%blocks(i)%ptr_nb_blocks(8))       !! HARDCODED FOR 2D Case
            read(1,*)       tree%blocks(i)%ptr_nb_blocks


            ! ------------------- Left Values -------------------
            read(1,*)    !      no_nb_cells_left:
            allocate(tree%blocks(i)%no_nb_cells_left(N+2))
            read(1,*)          tree%blocks(i)%no_nb_cells_left
            read(1,*)    !      no_coeff_left:
            read(1,*)       tree%blocks(i)%no_coeff_left
            allocate(tree%blocks(i)%nb_block_intrp_left(tree%blocks(i)%no_coeff_left))
            allocate(tree%blocks(i)%nb_cells_intrp_left(NDIM, tree%blocks(i)%no_coeff_left))
            allocate(tree%blocks(i)%nb_cells_coeff_left(tree%blocks(i)%no_coeff_left))
            read(1,*)    !      nb_block_intrp_left:
            read(1,*)       tree%blocks(i)%nb_block_intrp_left
            read(1,*)    !      nb_cells_intrp_left_x_y:
            do j = 1, NDIM
                read(1,*)       tree%blocks(i)%nb_cells_intrp_left(j,:)
            end do
            read(1,*)    !      nb_cells_coeff_left:
            read(1,*)       tree%blocks(i)%nb_cells_coeff_left

            ! ------------------- Right Values -------------------
            read(1,*)    !      no_nb_cells_right:
            allocate(tree%blocks(i)%no_nb_cells_right(N+2))
            read(1,*)          tree%blocks(i)%no_nb_cells_right
            read(1,*)    !      no_coeff_right:
            read(1,*)       tree%blocks(i)%no_coeff_right
            allocate(tree%blocks(i)%nb_block_intrp_right(tree%blocks(i)%no_coeff_right))
            allocate(tree%blocks(i)%nb_cells_intrp_right(NDIM, tree%blocks(i)%no_coeff_right))
            allocate(tree%blocks(i)%nb_cells_coeff_right(tree%blocks(i)%no_coeff_right))
            read(1,*)    !      nb_block_intrp_right:
            read(1,*)       tree%blocks(i)%nb_block_intrp_right
            read(1,*)    !      nb_cells_intrp_right_x_y:
            do j = 1, NDIM
                read(1,*)       tree%blocks(i)%nb_cells_intrp_right(j,:)
            end do
            read(1,*)    !      nb_cells_coeff_right:
            read(1,*)       tree%blocks(i)%nb_cells_coeff_right
            
            ! ------------------- Bottom Values -------------------
            read(1,*)    !      no_nb_cells_bottom:
            allocate(tree%blocks(i)%no_nb_cells_bottom(M))
            read(1,*)          tree%blocks(i)%no_nb_cells_bottom
            read(1,*)    !      no_coeff_bottom:
            read(1,*)       tree%blocks(i)%no_coeff_bottom
            allocate(tree%blocks(i)%nb_block_intrp_bottom(tree%blocks(i)%no_coeff_bottom))
            allocate(tree%blocks(i)%nb_cells_intrp_bottom(NDIM, tree%blocks(i)%no_coeff_bottom))
            allocate(tree%blocks(i)%nb_cells_coeff_bottom(tree%blocks(i)%no_coeff_bottom))
            read(1,*)    !      nb_block_intrp_bottom:
            read(1,*)       tree%blocks(i)%nb_block_intrp_bottom
            read(1,*)    !      nb_cells_intrp_bottom_x_y:
            do j = 1, NDIM
                read(1,*)       tree%blocks(i)%nb_cells_intrp_bottom(j,:)
            end do
            read(1,*)    !      nb_cells_coeff_bottom:
            read(1,*)       tree%blocks(i)%nb_cells_coeff_bottom

            ! ------------------- Top Values -------------------
            read(1,*)    !      no_nb_cells_top:
            allocate(tree%blocks(i)%no_nb_cells_top(M))
            read(1,*)          tree%blocks(i)%no_nb_cells_top
            read(1,*)    !      no_coeff_top:
            read(1,*)       tree%blocks(i)%no_coeff_top
            allocate(tree%blocks(i)%nb_block_intrp_top(tree%blocks(i)%no_coeff_top))
            allocate(tree%blocks(i)%nb_cells_intrp_top(NDIM, tree%blocks(i)%no_coeff_top))
            allocate(tree%blocks(i)%nb_cells_coeff_top(tree%blocks(i)%no_coeff_top))
            read(1,*)    !      nb_block_intrp_top:
            read(1,*)       tree%blocks(i)%nb_block_intrp_top
            read(1,*)    !      nb_cells_intrp_top_x_y:
            do j = 1, NDIM
                read(1,*)       tree%blocks(i)%nb_cells_intrp_top(j,:)
            end do
            read(1,*)    !      nb_cells_coeff_top:
            read(1,*)       tree%blocks(i)%nb_cells_coeff_top
            read(1,*)    !  -------------------------------------
            
            ! print *,  '-----', i, '-------'
            ! print *, tree%blocks(i)%nb_cells_coeff_top

        end do
        read(1,*)    !---------------------------------------------------------------------- 

        tree%no_leaf_blocks = block_count

        ! Allocating the number of multi_grid variable as max levels
        do idx = 1, tree%no_blocks
            allocate(tree%blocks(idx)%multi_grid(tree%no_levels))
        end do

        ! Allocating variables only for leaf blocks
        do lvl = tree%no_levels, 1, -1
            do id = 1, no_leaf_idx_lvl(lvl)
                idx = tree%levels(lvl)%idx_leaf(id)

                ! #IF (PPE_SOLVER==MULTI_GRID)
                allocate(tree%blocks(idx)%multi_grid(1)%p(M+2*ngl, N+2*ngl, 1))   !> Pressure in MultiGrid Level 1
                allocate(tree%blocks(idx)%multi_grid(1)%f(M+2*ngl, N+2*ngl, 1))   !> RHS of Pressure in MultiGrid Level 1
                allocate(tree%blocks(idx)%multi_grid(1)%r(M+2*ngl, N+2*ngl, 1))   !> Residue of Pressure Eq in MultiGrid Level 1
                ! #ELSE
                allocate(tree%blocks(idx)%p(M+2*ngl, N+2*ngl, 1))   !> Pressure
                ! #ENDIF
                allocate(tree%blocks(idx)%u  (M+2*ngl, N+2*ngl, 1))   !> u velocity
                allocate(tree%blocks(idx)%v  (M+2*ngl, N+2*ngl, 1))   !> v velocity
                allocate(tree%blocks(idx)%T  (M+2*ngl, N+2*ngl, 1))   !> Temperature
                allocate(tree%blocks(idx)%Ue (M+2*ngl, N+2*ngl, 1))  !> Ue
                allocate(tree%blocks(idx)%Uw (M+2*ngl, N+2*ngl, 1))  !> Uw
                allocate(tree%blocks(idx)%Vn (M+2*ngl, N+2*ngl, 1))  !> Vn
                allocate(tree%blocks(idx)%Vs (M+2*ngl, N+2*ngl, 1))  !> Vs
                allocate(tree%blocks(idx)%Cpx(M+2*ngl, N+2*ngl, 1)) !> Cpx
                allocate(tree%blocks(idx)%Cpy(M+2*ngl, N+2*ngl, 1)) !> Cpy
                allocate(tree%blocks(idx)%CpT(M+2*ngl, N+2*ngl, 1)) !> CpT
                allocate(tree%blocks(idx)%rhs(M+2*ngl, N+2*ngl, 1)) !> RHS for modified momentum equation

                ! allocate(tree%blocks(idx)%residue(M+2*ngl, N+2*ngl, 1)) !> residue for pressure equation

                allocate(tree%blocks(idx)%d (M+2*ngl, N+2*ngl, 1)) !> For Conjugate Gradient
                allocate(tree%blocks(idx)%r (M+2*ngl, N+2*ngl, 1)) 
                allocate(tree%blocks(idx)%Ad(M+2*ngl, N+2*ngl, 1)) 


            end do
        end do

    close(1)

    ! Reading the Solver File
    open(2, file = '..'//path_separator//solver_file)

        read(2,*)    !  ----- Solver Input -----
        read(2,*)    !  
        read(2,*)    !  Physical Information:
        read(2,*)    ! ---------------------------------------- 
        read(2,*)    !		Re	(Reynolds Number)
        read(2,*)        Re
        read(2,*)    !
        read(2,*)    !		Pr	(Prandlt Number)
        read(2,*)        Pr
        read(2,*)    !
        read(2,*)    !	max_time	delta_time (Time Step)
        read(2,*)        max_time, delta_time
        read(2,*)    !
        read(2,*)    !	Boundary Information:
        read(2,*)    ! ---------------------------------------- 
        read(2,*)    ! 
        read(2,*)    ! ----- VELOCITIES ----- 
        read(2,*)    ! U_vel	left		right		bottom		top	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
        read(2,*)       u_BC_left, u_BC_right, u_BC_bottom, u_BC_top
        read(2,*)    ! 
        read(2,*)    ! U_left		U_right		U_bottom	U_top	( U - Velocity Values )
        read(2,*)       u_left, u_right, u_bottom, u_top
        read(2,*)    ! 
        read(2,*)    ! dU_dx_left	dU_dx_right	dU_dy_bottom	dU_dy_top	( U - Velocity gradient Values )
        read(2,*)       du_dx_left, du_dx_right, du_dy_bottom, du_dy_top
        read(2,*)    !
        read(2,*)    ! V_vel	left		right		bottom		top	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
        read(2,*)       v_BC_left, v_BC_right, v_BC_bottom, v_BC_top
        read(2,*)    ! 
        read(2,*)    ! V_left		V_right		V_bottom	V_top	( V - Velocity Values )
        read(2,*)       v_left, v_right, v_bottom, v_top
        read(2,*)    ! 
        read(2,*)    ! dV_dx_left	dV_dx_right	dV_dy_bottom	dV_dy_top	( V - Velocity gradient Values )
        read(2,*)       dv_dx_left, dv_dx_right, dv_dy_bottom, dv_dy_top
        read(2,*)    !
        read(2,*)    ! ----- PRESSURE ----- 
        read(2,*)    ! P	left		right		bottom		top	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
        read(2,*)       p_BC_left, p_BC_right, p_BC_bottom, p_BC_top
        read(2,*)    ! 
        read(2,*)    ! p_left		p_right		p_bottom	p_top	(Pressure Values)
        read(2,*)       p_left, p_right, p_bottom, p_top
        read(2,*)    ! 
        read(2,*)    ! dp_dx_left	dp_dx_right	dp_dy_bottom	dp_dx_top	(Pressure gradient Values)
        read(2,*)       dp_dx_left, dp_dx_right, dp_dy_bottom, dp_dy_top
        read(2,*)    !
        read(2,*)    ! ----- TEMPERATURE ----- 
        read(2,*)    !  T	left		right		bottom		top	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
        read(2,*)       T_BC_left, T_BC_right, T_BC_bottom, T_BC_top
        read(2,*)    ! 
        read(2,*)    !  T_left		T_right		T_bottom	T_top	(Temperature Values)
        read(2,*)       T_left,     T_right,    T_bottom,   T_top
        read(2,*)    ! 
        read(2,*)    !  dT_dx_left	dT_dx_right	 dT_dy_bottom  dT_dx_top	(Temperature gradient Values)
        read(2,*)       dT_dx_left, dT_dx_right, dT_dy_bottom, dT_dy_top
        read(2,*)    !
        read(2,*)    !	Convergence Information:
        read(2,*)    ! ---------------------------------------- 
        read(2,*)    ! 	epsilon_velocity	epsilon_pressure    epsilon_temperature	max_iter
        read(2,*)       epsilon_u,          epsilon_p,          epsilon_T,          max_iter
        read(2,*)    !
        read(2,*)    !	Multi Grid Information:
        read(2,*)    ! ---------------------------------------- 
        read(2,*)    ! 	mu_1		mu_2		mu_0		weight_nu	weight_den
        read(2,*)       mu_1,       mu_2,       mu_0,       temp_num,   temp_den
        read(2,*)    !
        read(2,*)    !	Solver Options: 1 GAUSS_SEIDEL, 2 WEIGHTED_JACOBI, 3 CONJUGATE_GRADIENT, 4 MULTI_GRID
        read(2,*)    ! ---------------------------------------- 
        read(2,*)    ! 	MME_SOLVER	PPE_SOLVER	MG_SOLVER
        read(2,*)       MME_SOLVER,	PPE_SOLVER,	MG_SOLVER

        if (MME_SOLVER==GAUSS_SEIDEL)       mme_solver_text = 'GS'
        if (MME_SOLVER==WEIGHTED_JACOBI)    mme_solver_text = 'WJ'
        if (MME_SOLVER==CONJUGATE_GRADIENT) mme_solver_text = 'CG'
        if (MME_SOLVER==MULTI_GRID)         mme_solver_text = 'MG'

        if (PPE_SOLVER==GAUSS_SEIDEL)       ppe_solver_text = 'GS'
        if (PPE_SOLVER==WEIGHTED_JACOBI)    ppe_solver_text = 'WJ'
        if (PPE_SOLVER==CONJUGATE_GRADIENT) ppe_solver_text = 'CG'
        if (PPE_SOLVER==MULTI_GRID)         ppe_solver_text = 'MG'

        pressure_BC_neumann = .False.
        if (p_BC_left==NEUMANN .and. p_BC_right==NEUMANN .and. p_BC_bottom==NEUMANN .and. p_BC_top==NEUMANN) then
            pressure_BC_neumann = .True.
        end if

        weight = temp_num / temp_den

        epsilon_u = 10**(epsilon_u)
        epsilon_p = 10**(epsilon_p)
        epsilon_T = 10**(epsilon_T)

        beta    = (cell_dim(1)/cell_dim(2))
        beta_sq = beta**2

        gamma_p = -2.0_prcs_var*(1.0_prcs_var+beta_sq)

        dx_lvl(1)    = cell_dim(1)
        dy_lvl(1)    = cell_dim(2)

        do i = 1, tree%no_levels

            dx_lvl(i)    = dx_lvl(1)/2.0_prcs_var**(i-1)
            dy_lvl(i)    = dy_lvl(1)/2.0_prcs_var**(i-1)
            
            dx_sq_lvl(i) = dx_lvl(i)**2
            dy_sq_lvl(i) = dy_lvl(i)**2 

            alpha(i)     = dx_sq_lvl(i) * Re

            gamma_u(i)   = - 2.0_prcs_var*alpha(i)/delta_time + gamma_p
            gamma_T(i)   = - Re*Pr*dx_sq_lvl(i)/delta_time + gamma_p

            coeff_p(i)   = dx_lvl(i)/delta_time
            coeff_u(i)   = - 2.0_prcs_var*alpha(i)/delta_time
            coeff_T(i)   = dx_sq_lvl(i) * Re * Pr

        end do

        epsilon_p_sq = epsilon_p**2
        epsilon_u_sq = epsilon_u**2
        epsilon_T_sq = epsilon_T**2

    close(2)

    ! Reading the Canonical File
    open(unit=3, file = '..'//path_separator//'Body_Data'//path_separator//canonical_body_file)

        read(3, *)    !----- Canonical Body Data -----
        read(3, *)    !
        read(3, *)    ! Number of Bodies:
        read(3, *)    !----------------------------------------
        read(3, *)    !nbody	nbody_solid	nbody_holes
        read(3, *) nbody, nbody_solid, nbody_holes
        read(3, *)    !
        read(3, *)    ! Dimension of body:
        read(3, *)    ! ----------------------------------------
        read(3, *)    ! body_dim
        read(3, *) body_dim
        read(3, *)    !
        read(3, *)    ! Point outside Body:
        read(3, *)    ! ----------------------------------------
        read(3, *)    ! pointOutsideBodyX, pointOutsideBodyY, pointOutsideBodyZ
        read(3, *) pointOutsideBodyX, pointOutsideBodyY, pointOutsideBodyZ
        read(3, *)    !
        read(3, *)    ! Thermal Data:
        read(3, *)    ! ----------------------------------------
        read(3, *)    ! ktherm, source
        read(3, *) ktherm, source

    close(3)

    ! Reading Individial Body Files
    allocate(body_data(nbody))
    do ibody = 1, nbody

        write(body_file_name,'("Structured_Body_",i1.1,".txt")'), ibody
        open(4, file = '..'//path_separator//'Body_Data'//path_separator//body_file_name)

            read(4,*)   !  ----- Body in Domain -----
            read(4,*)   !  
            read(4,*)   !  Rectangular Body:
            read(4,*)   ! ---------------------------------------- 
            read(4,*)   !		x_min	x_max	y_min	y_max
            read(4,*)       body_data(ibody)%body_x_min, body_data(ibody)%body_x_max, &
                            body_data(ibody)%body_y_min, body_data(ibody)%body_y_max
            read(4,*)   !  
            read(4,*)   !  Circular Body:
            read(4,*)   ! ---------------------------------------- 
            read(4,*)   !		center_x	center_y	radius
            read(4,*)       body_data(ibody)%body_center_x, body_data(ibody)%body_center_y, body_data(ibody)%body_radius
            read(4,*)   !  
            read(4,*)   !  Body_Selection:
            read(4,*)   ! ---------------------------------------- 
            read(4,*)   !		Body Selected: (0:Rectangular, 1:Circular)
            read(4,*)       body_data(ibody)%body_select
            read(4,*)   !
            read(4,*)   !   Boundary Information: (Velocity and Pressure)
            read(4,*)   !   ----------------------------------------
            read(4,*)   !   
            read(4,*)   !   ----- VELOCITIES -----
            read(4,*)   !   U_BC	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
            read(4,*)       body_data(ibody)%Boundary_condition(U_VELOCITY)
            read(4,*)   !
            read(4,*)   !   U_val( U - Velocity Values )
            read(4,*)       body_data(ibody)%Boundary_value(U_VELOCITY)
            read(4,*)   !
            read(4,*)   !   dU_dx	( U - Velocity gradient Values )
            read(4,*)       body_data(ibody)%Boundary_gradient(U_VELOCITY)
            read(4,*)   !
            read(4,*)   !   V_BC	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
            read(4,*)       body_data(ibody)%Boundary_condition(V_VELOCITY)
            read(4,*)   !
            read(4,*)   !   V_val( V - Velocity Values )
            read(4,*)       body_data(ibody)%Boundary_value(V_VELOCITY)
            read(4,*)   !
            read(4,*)   !   dV_dy	( V - Velocity gradient Values )
            read(4,*)       body_data(ibody)%Boundary_gradient(V_VELOCITY)
            read(4,*)   !
            read(4,*)   !   ----- PERESSURE -----
            read(4,*)   !   p_BC	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
            read(4,*)       body_data(ibody)%Boundary_condition(PRESSURE)
            read(4,*)   !
            read(4,*)   !   p_val (Pressure Values)
            read(4,*)       body_data(ibody)%Boundary_value(PRESSURE)
            read(4,*)   !
            read(4,*)   !   dp_dn	(Pressure gradient Values)
            read(4,*)       body_data(ibody)%Boundary_gradient(PRESSURE)
            read(4,*)   !
            read(4,*)   !   ----- TEMPERATURE -----
            read(4,*)   !   T_BC	(1-> Diritchlet; 2-> Neumann; 3-> Mixed)
            read(4,*)       body_data(ibody)%Boundary_condition(TEMPERATURE)
            read(4,*)   !
            read(4,*)   !   T_val (Temperature Values)
            read(4,*)       body_data(ibody)%Boundary_value(TEMPERATURE)
            read(4,*)   !
            read(4,*)   !   dT_dn	(Temperature gradient Values)
            read(4,*)       body_data(ibody)%Boundary_gradient(TEMPERATURE)

            ! Resetting x and y range if body selected is circular
            if (body_data(ibody)%body_select==1) then
                body_data(ibody)%body_x_min = body_data(ibody)%body_center_x - body_data(ibody)%body_radius
                body_data(ibody)%body_x_max = body_data(ibody)%body_center_x + body_data(ibody)%body_radius
                body_data(ibody)%body_y_min = body_data(ibody)%body_center_y - body_data(ibody)%body_radius
                body_data(ibody)%body_y_max = body_data(ibody)%body_center_y + body_data(ibody)%body_radius
            end if

        close(4)

    end do

end subroutine read_input_memory_allocate


! --------------------------------------------------------------
!>        Array Initialisation
! --------------------------------------------------------------
subroutine array_initialise(tree)

    use precision_module, only: prcs_var
    use data_type_module
    implicit none
    type(t_tree), intent(inout) :: tree
    integer :: idx, lvl, id

    ! Initialise values of all physical variables
    do lvl = tree%no_levels, 1, -1
        do id = 1, size(tree%levels(lvl)%idx_leaf)
            idx = tree%levels(lvl)%idx_leaf(id)

            ! tree%blocks(idx)%u(:,:,:) = 0.0_prcs_var      !> Initialise u-velocity
            tree%blocks(idx)%v(:,:,:) = 0.0_prcs_var      !> Initialise v-velocity
            tree%blocks(idx)%T(:,:,:) = 0.0_prcs_var      !> Initialise T
            tree%blocks(idx)%p(:,:,:) = 0.0_prcs_var      !> Initialise pressure
            tree%blocks(idx)%multi_grid(1)%p(:,:,:) = 0.0_prcs_var      !> Initialise pressure

            tree%blocks(idx)%u(:,:,:) = 0.0_prcs_var      !> Initialise u-velocity
            ! tree%blocks(idx)%u(:,:,:) = 0.5_prcs_var      !> Initialise u-velocity
            ! tree%blocks(idx)%v(:,:,:) = 0.5_prcs_var      !> Initialise v-velocity

            tree%blocks(idx)%Ue(:,:,:) = 0.0_prcs_var     !> Initialise Ue
            tree%blocks(idx)%Uw(:,:,:) = 0.0_prcs_var     !> Initialise Uw
            tree%blocks(idx)%Vn(:,:,:) = 0.0_prcs_var     !> Initialise Vn
            tree%blocks(idx)%Vs(:,:,:) = 0.0_prcs_var     !> Initialise Vs

            tree%blocks(idx)%Cpx(:,:,:) = 0.0_prcs_var    !> Initialise Cpx
            tree%blocks(idx)%Cpy(:,:,:) = 0.0_prcs_var    !> Initialise Cpy
            tree%blocks(idx)%CpT(:,:,:) = 0.0_prcs_var    !> Initialise CpT

        end do
    end do

end subroutine array_initialise


! --------------------------------------------------------------
!>        Multi-Grid Variable Initialisation
! --------------------------------------------------------------
subroutine init_multigrid_vars(tree)

    use precision_module, only: prcs_var
    use problem_module
    use data_type_module, only: N, M, t_tree, i_start, i_end, j_start,j_end, ngl, NDIM
    implicit none

    type(t_tree), intent(inout) :: tree
    integer :: i, j, k
    integer :: idx, nb_idx, parent_idx, prev_idx
    integer :: mg_lvl_no, own_mg_lvl, max_lvl
    integer :: max_tree_lvl, tree_lvl_no
    integer :: side
    integer :: count

    max_lvl = tree%no_levels

    ! Initialise variables for blocks
    do idx = 1, tree%no_blocks

        own_mg_lvl = max_lvl - tree%blocks(idx)%lvl + 1

        ! Allocating dimension of variables within multi_grid variables (only for levels where data is required)
        if (tree%blocks(idx)%is_refined) then
            allocate(tree%blocks(idx)%multi_grid(own_mg_lvl)%p(M+2*ngl, N+2*ngl, 1))
            allocate(tree%blocks(idx)%multi_grid(own_mg_lvl)%f(M+2*ngl, N+2*ngl, 1))
            allocate(tree%blocks(idx)%multi_grid(own_mg_lvl)%r(M+2*ngl, N+2*ngl, 1))
        else
            do mg_lvl_no = 2, own_mg_lvl
                allocate(tree%blocks(idx)%multi_grid(mg_lvl_no)%p(M+2*ngl, N+2*ngl, 1))
                allocate(tree%blocks(idx)%multi_grid(mg_lvl_no)%f(M+2*ngl, N+2*ngl, 1))
                allocate(tree%blocks(idx)%multi_grid(mg_lvl_no)%r(M+2*ngl, N+2*ngl, 1))
            end do
        end if


        ! Creating ghost cell interpolation data for when neighbouring finer blocks are restricted
        allocate(tree%blocks(idx)%MG_no_nb_cells_left(N+2))
        allocate(tree%blocks(idx)%MG_no_nb_cells_right(N+2))
        allocate(tree%blocks(idx)%MG_no_nb_cells_bottom(M))
        allocate(tree%blocks(idx)%MG_no_nb_cells_top(M))
        
        ! Find the lvl of the blocks in ptr_nb_blocks (stores info of nbs of same or lower/coarser levels only)
        allocate(tree%blocks(idx)%nb_lvl_high(2**NDIM))

        do side = 1, 2**NDIM
            nb_idx = tree%blocks(idx)%ptr_nb_blocks(2*side)
            tree%blocks(idx)%nb_lvl_high(side) = .False.
            if (nb_idx==0) then
                tree%blocks(idx)%nb_lvl_high(side) = .True.
            else if (tree%blocks(nb_idx)%lvl==tree%blocks(idx)%lvl) then
                tree%blocks(idx)%nb_lvl_high(side) = .True.
            end if
            ! print *, 'nb_idx: ', nb_idx, 'is_true? : ', tree%blocks(idx)%nb_lvl_high(side)
        end do

        ! Creating ghost cell interpolation changing any finer refinement to the same level
        ! if (own_mg_lvl==1) cycle

        if (tree%blocks(idx)%nb_lvl_high(LEFT)) then
            allocate(tree%blocks(idx)%MG_nb_block_intrp_left(N+2))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_left(NDIM, N+2))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_left(N+2))
            tree%blocks(idx)%MG_no_nb_cells_left(:)      = 1
            tree%blocks(idx)%MG_no_coeff_left            = N+2
            tree%blocks(idx)%MG_nb_block_intrp_left(:)   = tree%blocks(idx)%ptr_nb_blocks(2*LEFT)
            tree%blocks(idx)%MG_nb_cells_intrp_left(1,:) = i_end
            tree%blocks(idx)%MG_nb_cells_intrp_left(2,:) = [(j, j=j_start-ngl, j_end+ngl)]
            tree%blocks(idx)%MG_nb_cells_coeff_left(:)   = 1.0_prcs_var
        else
            allocate(tree%blocks(idx)%MG_nb_block_intrp_left(tree%blocks(idx)%no_coeff_left))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_left(NDIM, tree%blocks(idx)%no_coeff_left))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_left(tree%blocks(idx)%no_coeff_left))
            tree%blocks(idx)%MG_no_nb_cells_left(:)      = tree%blocks(idx)%no_nb_cells_left(:)
            tree%blocks(idx)%MG_no_coeff_left            = tree%blocks(idx)%no_coeff_left
            tree%blocks(idx)%MG_nb_block_intrp_left(:)   = tree%blocks(idx)%nb_block_intrp_left(:)
            tree%blocks(idx)%MG_nb_cells_intrp_left(1,:) = tree%blocks(idx)%nb_cells_intrp_left(1,:)
            tree%blocks(idx)%MG_nb_cells_intrp_left(2,:) = tree%blocks(idx)%nb_cells_intrp_left(2,:)
            tree%blocks(idx)%MG_nb_cells_coeff_left(:)   = tree%blocks(idx)%nb_cells_coeff_left(:)
        end if
    
        if (tree%blocks(idx)%nb_lvl_high(RIGHT)) then
            allocate(tree%blocks(idx)%MG_nb_block_intrp_right(N+2))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_right(NDIM, N+2))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_right(N+2))
            tree%blocks(idx)%MG_no_nb_cells_right(:)      = 1
            tree%blocks(idx)%MG_no_coeff_right            = N+2
            tree%blocks(idx)%MG_nb_block_intrp_right(:)   = tree%blocks(idx)%ptr_nb_blocks(2*RIGHT)
            tree%blocks(idx)%MG_nb_cells_intrp_right(1,:) = i_start
            tree%blocks(idx)%MG_nb_cells_intrp_right(2,:) = [(j, j=j_start-ngl, j_end+ngl)]
            tree%blocks(idx)%MG_nb_cells_coeff_right(:)   = 1.0_prcs_var
        else
            allocate(tree%blocks(idx)%MG_nb_block_intrp_right(tree%blocks(idx)%no_coeff_right))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_right(NDIM, tree%blocks(idx)%no_coeff_right))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_right(tree%blocks(idx)%no_coeff_right))
            tree%blocks(idx)%MG_no_nb_cells_right(:)      = tree%blocks(idx)%no_nb_cells_right(:)
            tree%blocks(idx)%MG_no_coeff_right            = tree%blocks(idx)%no_coeff_right
            tree%blocks(idx)%MG_nb_block_intrp_right(:)   = tree%blocks(idx)%nb_block_intrp_right(:)
            tree%blocks(idx)%MG_nb_cells_intrp_right(1,:) = tree%blocks(idx)%nb_cells_intrp_right(1,:)
            tree%blocks(idx)%MG_nb_cells_intrp_right(2,:) = tree%blocks(idx)%nb_cells_intrp_right(2,:)
            tree%blocks(idx)%MG_nb_cells_coeff_right(:)   = tree%blocks(idx)%nb_cells_coeff_right(:)
        end if
    
        if (tree%blocks(idx)%nb_lvl_high(BOTTOM)) then
            allocate(tree%blocks(idx)%MG_nb_block_intrp_bottom(M))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_bottom(NDIM, M))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_bottom(M))
            tree%blocks(idx)%MG_no_nb_cells_bottom(:)      = 1
            tree%blocks(idx)%MG_no_coeff_bottom            = M
            tree%blocks(idx)%MG_nb_block_intrp_bottom(:)   = tree%blocks(idx)%ptr_nb_blocks(2*BOTTOM)
            tree%blocks(idx)%MG_nb_cells_intrp_bottom(1,:) = [(i, i=i_start, i_end)]
            tree%blocks(idx)%MG_nb_cells_intrp_bottom(2,:) = j_end
            tree%blocks(idx)%MG_nb_cells_coeff_bottom(:)   = 1.0_prcs_var
        else
            allocate(tree%blocks(idx)%MG_nb_block_intrp_bottom(tree%blocks(idx)%no_coeff_bottom))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_bottom(NDIM, tree%blocks(idx)%no_coeff_bottom))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_bottom(tree%blocks(idx)%no_coeff_bottom))
            tree%blocks(idx)%MG_no_nb_cells_bottom(:)      = tree%blocks(idx)%no_nb_cells_bottom(:)
            tree%blocks(idx)%MG_no_coeff_bottom            = tree%blocks(idx)%no_coeff_bottom
            tree%blocks(idx)%MG_nb_block_intrp_bottom(:)   = tree%blocks(idx)%nb_block_intrp_bottom(:)
            tree%blocks(idx)%MG_nb_cells_intrp_bottom(1,:) = tree%blocks(idx)%nb_cells_intrp_bottom(1,:)
            tree%blocks(idx)%MG_nb_cells_intrp_bottom(2,:) = tree%blocks(idx)%nb_cells_intrp_bottom(2,:)
            tree%blocks(idx)%MG_nb_cells_coeff_bottom(:)   = tree%blocks(idx)%nb_cells_coeff_bottom(:)
        end if
    
        if (tree%blocks(idx)%nb_lvl_high(TOP)) then
            allocate(tree%blocks(idx)%MG_nb_block_intrp_top(M))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_top(NDIM, M))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_top(M))
            tree%blocks(idx)%MG_no_nb_cells_top(:)      = 1
            tree%blocks(idx)%MG_no_coeff_top            = M
            tree%blocks(idx)%MG_nb_block_intrp_top(:)   = tree%blocks(idx)%ptr_nb_blocks(2*TOP)
            tree%blocks(idx)%MG_nb_cells_intrp_top(1,:) = [(i, i=i_start, i_end)]
            tree%blocks(idx)%MG_nb_cells_intrp_top(2,:) = j_start
            tree%blocks(idx)%MG_nb_cells_coeff_top(:)   = 1.0_prcs_var
        else
            allocate(tree%blocks(idx)%MG_nb_block_intrp_top(tree%blocks(idx)%no_coeff_top))
            allocate(tree%blocks(idx)%MG_nb_cells_intrp_top(NDIM, tree%blocks(idx)%no_coeff_top))
            allocate(tree%blocks(idx)%MG_nb_cells_coeff_top(tree%blocks(idx)%no_coeff_top))
            tree%blocks(idx)%MG_no_nb_cells_top(:)      = tree%blocks(idx)%no_nb_cells_top(:)
            tree%blocks(idx)%MG_no_coeff_top            = tree%blocks(idx)%no_coeff_top
            tree%blocks(idx)%MG_nb_block_intrp_top(:)   = tree%blocks(idx)%nb_block_intrp_top(:)
            tree%blocks(idx)%MG_nb_cells_intrp_top(1,:) = tree%blocks(idx)%nb_cells_intrp_top(1,:)
            tree%blocks(idx)%MG_nb_cells_intrp_top(2,:) = tree%blocks(idx)%nb_cells_intrp_top(2,:)
            tree%blocks(idx)%MG_nb_cells_coeff_top(:)   = tree%blocks(idx)%nb_cells_coeff_top(:)
        end if
        
    end do
    
    ! Initialise variables for levels
    allocate(tree%mg_lvl(max_lvl))
    
    ! For each mg_level, we find the boundary blocks at the level required
    do mg_lvl_no = 1, max_lvl
        
        ! LEFT BOUNDARY BLOCKS
            parent_idx = tree%left_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            count = 1
            do i = 2, size(tree%left_blocks)
                parent_idx = tree%left_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                end if
            end do
            allocate(tree%mg_lvl(mg_lvl_no)%MG_left_blocks(count))
            parent_idx = tree%left_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            tree%mg_lvl(mg_lvl_no)%MG_left_blocks(1) = parent_idx
            count = 1
            do i = 2, size(tree%left_blocks)
                parent_idx = tree%left_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                    tree%mg_lvl(mg_lvl_no)%MG_left_blocks(count) = parent_idx
                end if
            end do
        !

        ! RIGHT BOUNDARY BLOCKS
            parent_idx = tree%right_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            count = 1
            do i = 2, size(tree%right_blocks)
                parent_idx = tree%right_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                end if
            end do
            allocate(tree%mg_lvl(mg_lvl_no)%MG_right_blocks(count))
            parent_idx = tree%right_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            tree%mg_lvl(mg_lvl_no)%MG_right_blocks(1) = parent_idx
            count = 1
            do i = 2, size(tree%right_blocks)
                parent_idx = tree%right_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                    tree%mg_lvl(mg_lvl_no)%MG_right_blocks(count) = parent_idx
                end if
            end do
        !

        ! BOTTOM BOUNDARY BLOCKS
            parent_idx = tree%bottom_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            count = 1
            do i = 2, size(tree%bottom_blocks)
                parent_idx = tree%bottom_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                end if
            end do
            allocate(tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(count))
            parent_idx = tree%bottom_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(1) = parent_idx
            count = 1
            do i = 2, size(tree%bottom_blocks)
                parent_idx = tree%bottom_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                    tree%mg_lvl(mg_lvl_no)%MG_bottom_blocks(count) = parent_idx
                end if
            end do
        !

        ! TOP BOUNDARY BLOCKS
            parent_idx = tree%top_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            count = 1
            do i = 2, size(tree%top_blocks)
                parent_idx = tree%top_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                end if
            end do
            allocate(tree%mg_lvl(mg_lvl_no)%MG_top_blocks(count))
            parent_idx = tree%top_blocks(1)
            do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                parent_idx = tree%blocks(parent_idx)%ptr_parent
            end do
            prev_idx = parent_idx
            tree%mg_lvl(mg_lvl_no)%MG_top_blocks(1) = parent_idx
            count = 1
            do i = 2, size(tree%top_blocks)
                parent_idx = tree%top_blocks(i)
                do while (tree%blocks(parent_idx)%lvl>(max_lvl-mg_lvl_no+1))
                    parent_idx = tree%blocks(parent_idx)%ptr_parent
                end do
                if (parent_idx==prev_idx) then 
                    cycle
                else
                    prev_idx = parent_idx
                    count = count + 1
                    tree%mg_lvl(mg_lvl_no)%MG_top_blocks(count) = parent_idx
                end if
            end do
        !

        ! Storing Index of blocks at each level while at mg_lvl_no
        max_tree_lvl = max_lvl-mg_lvl_no+1
        allocate(tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(max_tree_lvl))
        do tree_lvl_no = 1, max_tree_lvl
            if (tree_lvl_no==max_tree_lvl) then
                tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf = &
                    tree%levels(tree_lvl_no)%idx_all
            else 
                tree%mg_lvl(mg_lvl_no)%mg_lvl_idx(tree_lvl_no)%idx_leaf = &
                    tree%levels(tree_lvl_no)%idx_leaf
            end if
        end do

    end do


end subroutine init_multigrid_vars
