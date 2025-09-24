
    #include "definitions.h"
    
!        ____________ MODULES ____________

! --------------------------------------------------------------
!         Precision Module
! --------------------------------------------------------------
module precision_module

    implicit none
    integer, parameter    :: prcs_var   = selected_real_kind(p=8, r=10)     !> precision with precision of p digits, exponent range of r
    integer, parameter    :: prcs_cell  = selected_int_kind(2)              !> int wtih 2 bytes
    integer, parameter    :: prcs_block = selected_int_kind(4)              !> int wtih 4 bytes
    
end module precision_module


! --------------------------------------------------------------
!         Tree DataType Module
! --------------------------------------------------------------
module data_type_module

    use precision_module, only: prcs_var
    implicit none

    !> Input Details
    integer :: NDIM     !> Dimension of the problem
    integer :: M, N     !> Number of cells within the domain on x and y direction
    integer :: ngl      !> Number of ghost layers used
    integer :: i_start, i_end, i_mid, j_start, j_end, j_mid !> Start, End and Middle index for the arrays along i and j
    integer, allocatable                :: no_leaf_idx_lvl(:)         !> Number of leaf blocks in each level
    real(kind=prcs_var), allocatable    :: domain(:)               !> Domain size
    real(kind=prcs_var), allocatable    :: cell_dim(:)             !> Cell size
    real(kind=prcs_var), allocatable    :: dx_lvl(:), dy_lvl(:)    !> dx and dy at different levels
    real(kind=prcs_var), allocatable    :: dx_sq_lvl(:), dy_sq_lvl(:)    !> dx^2 and dy^2 at different levels
    real(kind=prcs_var)                 :: beta, beta_sq                 !> beta = (dx/dy)

    real(kind=prcs_var)                 :: rhs_sum, rhs_sum_prev

    !> Multigrid Data
    integer             :: mu_0, mu_1, mu_2
    real(kind=prcs_var) :: weight

    !> Body Data, number of interpolating points
    integer     :: no_interp_pts_body


    !> Type which stored the data of each block
    type t_block
        integer     :: idx              !> Index of the block
        integer     :: lvl              !> Level at which the block is at
        logical     :: is_refined       !> Refinement flag
        real, allocatable  :: position(:)   !> Position of center of block

        !> Pointers
        integer               :: ptr_parent       
        integer, allocatable  :: ptr_children(:)   
        integer, allocatable  :: ptr_nb_blocks(:)   

        !> Ghost Cell Interpolating data
        integer ::  no_coeff_left, &
                    no_coeff_right, &
                    no_coeff_bottom, &
                    no_coeff_top
        integer, allocatable :: no_nb_cells_left(:), &
                                no_nb_cells_right(:), &
                                no_nb_cells_bottom(:), &
                                no_nb_cells_top(:)
        integer, allocatable, dimension(:)  ::  nb_block_intrp_left, &
                                                nb_block_intrp_right, &
                                                nb_block_intrp_bottom, &
                                                nb_block_intrp_top
        integer, allocatable, dimension(:,:)    ::  nb_cells_intrp_left, &
                                                    nb_cells_intrp_right, &
                                                    nb_cells_intrp_bottom, &
                                                    nb_cells_intrp_top
        real(kind=prcs_var), allocatable, dimension(:)  ::  nb_cells_coeff_left, &
                                                            nb_cells_coeff_right, &
                                                            nb_cells_coeff_bottom, &
                                                            nb_cells_coeff_top
        
        !> Physical Quantities
        real(kind=prcs_var), allocatable   :: p(:,:,:)                                           !> Pressure
        real(kind=prcs_var), allocatable   :: u(:,:,:), v(:,:,:), &
                                              Ue(:,:,:), Uw(:,:,:), Vn(:,:,:), Vs(:,:,:), &      !> Velocties
                                              Cpx(:,:,:), Cpy(:,:,:), CpT(:,:,:)                 !> Convective terms
        real(kind=prcs_var), allocatable   :: T(:,:,:)                                           !> Temperature
        real(kind=prcs_var), allocatable, dimension(:,:,:)  :: rhs          !> RHS term of the equations A.x = RHS
            


        !> For Conjugate Gradient Method
        real(kind=prcs_var), allocatable   :: d(:,:,:)
        real(kind=prcs_var), allocatable   :: r(:,:,:)
        real(kind=prcs_var), allocatable   :: Ad(:,:,:)
        ! real(kind=prcs_var)   :: rTr 
        ! real(kind=prcs_var)   :: cg_alpha
        ! real(kind=prcs_var)   :: cg_beta

        !> For Multi-Grid Method
        type(t_multi_grid), allocatable :: multi_grid(:)
        logical, allocatable :: nb_lvl_high(:)
        ! Ghost cells for Multi Grid when finer nb blocks aren't used
        integer ::  MG_no_coeff_left, &
                    MG_no_coeff_right, &
                    MG_no_coeff_bottom, &
                    MG_no_coeff_top
        integer, allocatable :: MG_no_nb_cells_left(:), &
                                MG_no_nb_cells_right(:), &
                                MG_no_nb_cells_bottom(:), &
                                MG_no_nb_cells_top(:)
        integer, allocatable, dimension(:)  ::  MG_nb_block_intrp_left, &
                                                MG_nb_block_intrp_right, &
                                                MG_nb_block_intrp_bottom, &
                                                MG_nb_block_intrp_top
        integer, allocatable, dimension(:,:)    ::  MG_nb_cells_intrp_left, &
                                                    MG_nb_cells_intrp_right, &
                                                    MG_nb_cells_intrp_bottom, &
                                                    MG_nb_cells_intrp_top
        real(kind=prcs_var), allocatable, dimension(:)  ::  MG_nb_cells_coeff_left, &
                                                            MG_nb_cells_coeff_right, &
                                                            MG_nb_cells_coeff_bottom, &
                                                            MG_nb_cells_coeff_top

        !> To deal with body in the domain
        integer*1, allocatable  :: i_blank(:,:,:,:)    ! To store flag in cells within blocks which are solid (1) and fluid (0)
        integer*1, allocatable  :: i_blank_sum(:,:,:)  ! Common i_blank for all bodies
        integer  , allocatable  :: i_range(:,:,:)      ! To store the range at which i loop should run through
        
        ! integer*1, allocatable  :: i_blank_singularity(:,:,:)    ! i_blank modified to consider singularity cells as solid (1)
        ! integer  , allocatable  :: i_range_singularity(:,:,:)    ! i_range excluding the cells which have same ghost cells. (singularity condition)

        integer, allocatable    :: no_IB_ghost_cells(:)    ! Number of IB ghost cells in the block
        integer, allocatable    :: IB_ghost_cell(:,:)   ! The IB ghost cells in the block
        integer, allocatable    :: IP_interpol_cells(:,:,:) ! Cells from which the value of Image Point is interpolated
        real(kind=prcs_var), allocatable    :: IP_pos_val(:,:)  ! The position and value of the Image Point
        real(kind=prcs_var), allocatable    :: dl_GC_IP(:)    ! Outward lengths of GC to IP

    end type t_block

    !> Type list of ids of boxes (all, leaves, body in it)
    type t_level
        integer, allocatable    :: idx_all(:)     !> Index of the blocks at that level
        integer, allocatable    :: idx_leaf(:)    !> Index of the leaf blocks at that level
        integer, allocatable    :: idx_body(:)    !> Index of the blocks which have the body in it
        integer, allocatable    :: idx_body_leaf(:)    !> Index of the leaf blocks which have the body in it
    end type t_level

    !> Type that stores the data of the tree
    type t_tree
        integer, allocatable :: n_cell(:)   !> Number of cells per dimension in block
        integer :: no_levels                    
        integer :: no_blocks
        integer :: no_leaf_blocks
        integer :: no_bound
        integer, allocatable :: left_blocks(:), right_blocks(:), bottom_blocks(:), top_blocks(:)
        type(t_level), allocatable  :: levels(:)    !> Storing the ids of each level
        type(t_block), allocatable  :: blocks(:)    !> Storing all the blocks used in the tree

        !> Boundary blocks during each step of Multi_Grid methods
        type(t_multi_grid_lvl), allocatable :: mg_lvl(:)
    end type  t_tree

    !> Type storing the data for multigrid blocks
    type t_multi_grid
        real(kind=prcs_var), allocatable   :: p(:,:,:) ! in equation Ap = f
        real(kind=prcs_var), allocatable   :: f(:,:,:)
        real(kind=prcs_var), allocatable   :: r(:,:,:)  ! r = f - Ap
    end type t_multi_grid

    !> Type for data on each multigrid level
    type t_multi_grid_lvl
        integer, allocatable :: MG_left_blocks(:), MG_right_blocks(:), MG_bottom_blocks(:), MG_top_blocks(:)
        type(t_level), allocatable :: mg_lvl_idx(:)
    end type t_multi_grid_lvl

end module data_type_module


! --------------------------------------------------------------
!         Problem Inputs Module
! --------------------------------------------------------------
module problem_module

    use precision_module, only: prcs_var
    implicit none

    integer                :: MME_SOLVER, PPE_SOLVER, MG_SOLVER
    character(len=2)       :: mme_solver_text, ppe_solver_text

    real(kind=prcs_var)    :: Re                            !> Reynolds Number
    real(kind=prcs_var)    :: Pr                            !> Prandtl Number
    real(kind=prcs_var)    :: time, delta_time, max_time    !> Time Step

    integer                :: u_BC_left, u_BC_right, u_BC_bottom, u_BC_top         !> u_velocity Boundary Conditions
    real(kind=prcs_var)    :: u_left, u_right, u_top, u_bottom                     !> U Boundary velocity
    real(kind=prcs_var)    :: du_dx_left, du_dx_right, du_dy_top, du_dy_bottom     !> U Boundary velocity Gradient

    integer                :: v_BC_left, v_BC_right, v_BC_bottom, v_BC_top         !> v_velocity Boundary Conditions
    real(kind=prcs_var)    :: v_left, v_right, v_top, v_bottom                     !> v Boundary velocity
    real(kind=prcs_var)    :: dv_dx_left, dv_dx_right, dv_dy_top, dv_dy_bottom     !> v Boundary velocity Gradient

    integer                :: p_BC_left, p_BC_right, p_BC_bottom, p_BC_top         !> pressure Boundary Conditions
    real(kind=prcs_var)    :: p_left, p_right, p_top, p_bottom                     !> Boundary pressure
    real(kind=prcs_var)    :: dp_dx_left, dp_dx_right, dp_dy_top, dp_dy_bottom     !> Boundary pressure Gradient

    integer                :: T_BC_left, T_BC_right, T_BC_bottom, T_BC_top         !> Temperature Boundary Conditions
    real(kind=prcs_var)    :: T_left, T_right, T_top, T_bottom                     !> Boundary Temperature
    real(kind=prcs_var)    :: dT_dx_left, dT_dx_right, dT_dy_top, dT_dy_bottom     !> Boundary Temperature Gradient

    real(kind=prcs_var), allocatable    :: alpha(:)                     !> alpha value at each level (alpha = dx Re / dt)
    real(kind=prcs_var)                 :: gamma_p                      !> gamma_pressure value at each level    (gamma_p = -2(1+beta^2))
    real(kind=prcs_var), allocatable    :: gamma_u(:)                   !> gamma_velocity value at each level    (gamma_u = (-2*alpha/dt + gamma_p)
    real(kind=prcs_var), allocatable    :: gamma_T(:)                   !> gamma_temperature value at each level (gamma_T = (-Re*Pr*dx^2/dt + gamma_p)
    real(kind=prcs_var), allocatable    :: coeff_p(:)                   !> coeff for RHS in p (coeff_p = dx^2 / dt)
    real(kind=prcs_var), allocatable    :: coeff_u(:)                   !> coeff for RHS in v (coeff_u = -2*alpha/dt + gamma_p)
    real(kind=prcs_var), allocatable    :: coeff_T(:)                   !> coeff for RHS in T (coeff_T = dx^2*Re*Pr)
    
    real(kind=prcs_var)    :: epsilon_u, epsilon_p, epsilon_T,&         !> error cutoff value
                              epsilon_u_sq, epsilon_p_sq, epsilon_T_sq  !> error cutoff value squared        
    real(kind=prcs_var)    :: residual_u, residual_v, residual_p, residual_T     !> error values of velocities and pressure
                              
    logical :: pressure_BC_neumann
    integer :: no_fluid_cells
    integer :: no_solid_cells

    integer :: max_iter            !> max iter to run
    integer :: iter                !> iteratons run
    integer :: time_step           !> time step used for iteration

    real(kind=prcs_var)     :: courant_u_max, courant_v_max
    
    character(len=:), allocatable   :: mesh_file, solver_file, canonical_body_file

end module problem_module


! --------------------------------------------------------------
!         Body in Domain Inputs Module
! --------------------------------------------------------------
module body_module

    use precision_module, only: prcs_var

    integer :: nbody, nbody_solid, nbody_holes
    integer :: body_dim

    integer                          :: no_singularities
    integer            , allocatable :: body_singularity_cells(:,:)  !> (block, cell_x, cell_y, ghost_x, ghost_y, i_body), no_singularities
    real(kind=prcs_var), allocatable :: body_singularity_dl_GC_IP(:) !> no_singularities

    real(kind=prcs_var) :: pointOutsideBodyX, pointOutsideBodyY, pointOutsideBodyZ
    real(kind=prcs_var) :: ktherm, source

    type t_body
        !> Square / Rectangular Body
        real(kind=prcs_var) :: body_x_min
        real(kind=prcs_var) :: body_x_max
        real(kind=prcs_var) :: body_y_min
        real(kind=prcs_var) :: body_y_max

        !> Circular Body
        real(kind=prcs_var) :: body_center_x
        real(kind=prcs_var) :: body_center_y
        real(kind=prcs_var) :: body_radius

        !> Selected Body
        integer :: body_select

        !> Boundary condition of the bodies
        integer                :: Boundary_condition(4)         !> Boundary Conditions u, v, p, T
        real(kind=prcs_var)    :: Boundary_value(4)     !> Values
        real(kind=prcs_var)    :: Boundary_gradient(4)     !> Gradient

    end type t_body

    type(t_body), allocatable   :: body_data(:)

end module body_module