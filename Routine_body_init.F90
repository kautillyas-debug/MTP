
    #include "definitions.h"

! --------------------------------------------------------------
!         Initialise variables for body in domain
! --------------------------------------------------------------

!> \brief Initialise the datastructure and variables for dealing with body in domain
!> \param[in,out] tree: The tree structure used.
subroutine init_body_vars(tree)

    use data_type_module, only: t_tree
    use data_type_module, only: M, N

    implicit none
    type(t_tree), intent(inout) :: tree

    call find_idx_body(tree)        ! Find the indices which have the body in it
    call find_i_blank(tree)         ! Find the cells which have block in it
    call find_IB_ghost_cell(tree)   ! Find the IB ghost cells at each block
    call find_i_range(tree)         ! Find the range of i variable
    call find_Image_Point_IP(tree)      ! Find the range of i variable
    call find_IP_interpol_cells(tree)   ! Find the neighbouring cells which will be used to interpolate to IP

    ! call find_singluarity_cells(tree)   ! Find the cells (singularity cells) which have common ghost cells in the body
    ! call find_i_range_singularity(tree) ! Find the i_range for the loop without solid cells and singularity cells.

end subroutine init_body_vars


!> \brief Find the blocks which would have the boundary of the body in it
!> \param[in,out] tree: The tree structure used.
subroutine find_idx_body(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: block_count      ! Temp variable to count which blocks have objects
    integer :: leaf_block_count ! Temp variable to count which leaf blocks have objects
    integer :: no_blocks_lvl    ! Temp variable to find max number of blocks in a level
    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: ibody            ! Variable to iterate through nbody
    integer :: lvl              ! Variables to iterate through the levels present
    logical, allocatable    :: is_outside(:)    ! T/F to check if block is outside the bounding rectangle (body_x_min,... etc.)
    logical                 :: is_counted       ! T/F to check if a block is already counted for having a block in it
    real(kind=prcs_var), allocatable :: block_dim(:)    ! Dimensions of the block selected
    real(kind=prcs_var), allocatable :: idx_pos(:)      ! Position of centre of block idx

    !> \todo Read all the body data files here, and compile the idx_body for all.

    ! Allocataing variables based on dimension of problem
    allocate(is_outside(NDIM*2))
    allocate(block_dim(NDIM))
    allocate(idx_pos(NDIM))


    do lvl = 1, tree%no_levels

        block_count = 1
        leaf_block_count = 1
        no_blocks_lvl = size(tree%levels(lvl)%idx_all)
        block_dim(1)  = cell_dim(1)/(2.0_prcs_var**(lvl-1))*M
        block_dim(2)  = cell_dim(2)/(2.0_prcs_var**(lvl-1))*N


        ! Allocating idx_body as size of idx_all + 1 (1st index to record the number of blocks are selected)
        allocate(tree%levels(lvl)%idx_body(0:no_blocks_lvl))
        allocate(tree%levels(lvl)%idx_body_leaf(0:no_blocks_lvl))

        do id = 1, no_blocks_lvl
            idx = tree%levels(lvl)%idx_all(id)
            idx_pos(:) = tree%blocks(idx)%position
            is_counted = .False.

            do ibody = 1, nbody

                if (is_counted) cycle
                
                ! Checking if the block is outside the bounding rectangle (written for 2D domain)
                is_outside(1) = ((idx_pos(1) + block_dim(1)/2.0_prcs_var) < body_data(ibody)%body_x_min)
                is_outside(2) = ((idx_pos(1) - block_dim(1)/2.0_prcs_var) > body_data(ibody)%body_x_max)
                is_outside(3) = ((idx_pos(2) + block_dim(2)/2.0_prcs_var) < body_data(ibody)%body_y_min)
                is_outside(4) = ((idx_pos(2) - block_dim(2)/2.0_prcs_var) > body_data(ibody)%body_y_max)

                ! Finding if block is not outside the bounding rectangle
                if (.not. (is_outside(1) .or. is_outside(2) .or. is_outside(3) .or. is_outside(4))) then
                    tree%levels(lvl)%idx_body(block_count) = idx
                    block_count = block_count + 1
                    is_counted = .True.
                    if (.not. tree%blocks(idx)%is_refined) then
                        tree%levels(lvl)%idx_body_leaf(leaf_block_count) = idx
                        leaf_block_count = leaf_block_count + 1
                    end if

                    ! Allocating memory for i_blank and i_boundary only on blocks with body in it
                    allocate(tree%blocks(idx)%i_blank(M+2*ngl, N+2*ngl, 1, nbody))
                    allocate(tree%blocks(idx)%i_blank_sum(M+2*ngl, N+2*ngl, 1))
                    ! allocate(tree%blocks(idx)%i_blank_singularity(M+2*ngl, N+2*ngl, 1))

                    tree%blocks(idx)%i_blank             = FLUID_PRESENT
                    tree%blocks(idx)%i_blank_sum         = FLUID_PRESENT
                    ! tree%blocks(idx)%i_blank_singularity = FLUID_PRESENT

                end if

            end do

        end do

        ! Storing the number of idx_body variables to use
        tree%levels(lvl)%idx_body(0) = block_count - 1
        tree%levels(lvl)%idx_body_leaf(0) = leaf_block_count - 1

        ! print *
        ! print *, '---------------------------------------'
        ! print *, 'Level: ', lvl
        ! print *, '---------------------------------------'
        ! print *, tree%levels(lvl)%idx_body

    end do

end subroutine find_idx_body


!> \brief Find the cells within blocks (i_blank) which are solid (1) and fluid (0)
!> \param[in,out] tree: The tree structure used.
!> \todo : This if blocks should be changed when Geometry of Body is Input.
subroutine find_i_blank(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module
    use problem_module, only: no_solid_cells, no_fluid_cells

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: block_count      ! Variable to count which blocks have objects
    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present
    integer :: ibody            ! Variable to iterate through nbody
    real(kind=prcs_var), allocatable :: cell_dim_lvl(:)     ! Dimensions of the cells in level selected
    real(kind=prcs_var), allocatable :: cell_pos(:)         ! Position of the cell selected
    real(kind=prcs_var), allocatable :: cell_pos_start(:)   ! Position of the bottom left cell (ghost cell)

    logical, allocatable :: is_inside(:)    ! flag to check if the cell center is in the body

    integer             :: i, j, k
    real(kind=prcs_var) :: distance

    no_solid_cells = 0

    
    allocate(cell_dim_lvl(NDIM))
    allocate(cell_pos(NDIM))
    allocate(cell_pos_start(NDIM))
    allocate(is_inside(2*NDIM))

    do lvl = 1, tree%no_levels

        cell_dim_lvl(1)  = cell_dim(1)/(2.0_prcs_var**(lvl-1))
        cell_dim_lvl(2)  = cell_dim(2)/(2.0_prcs_var**(lvl-1))

        do id = 1, tree%levels(lvl)%idx_body(0)
            idx = tree%levels(lvl)%idx_body(id)


            cell_pos_start(1) = tree%blocks(idx)%position(1) - cell_dim_lvl(1)*real(M/2 + 0.5, kind=prcs_var)
            cell_pos_start(2) = tree%blocks(idx)%position(2) - cell_dim_lvl(2)*real(N/2 + 0.5, kind=prcs_var)

            k = 1
            do j = j_start-ngl, j_end+ngl
                do i = i_start-ngl, i_end+ngl

                    cell_pos(1) = cell_pos_start(1) + real((i-1)*cell_dim_lvl(1), kind=prcs_var)
                    cell_pos(2) = cell_pos_start(2) + real((j-1)*cell_dim_lvl(2), kind=prcs_var)

                    do ibody = 1, nbody

                        !> \todo : This if blocks should be changed when Geometry of Body is Input.
                        if (body_data(ibody)%body_select == RECTANGLE) then
                            is_inside(1) = (cell_pos(1) >= body_data(ibody)%body_x_min)
                            is_inside(2) = (cell_pos(1) <= body_data(ibody)%body_x_max)
                            is_inside(3) = (cell_pos(2) >= body_data(ibody)%body_y_min)
                            is_inside(4) = (cell_pos(2) <= body_data(ibody)%body_y_max)

                            if (is_inside(1) .and. is_inside(2) .and. is_inside(3) .and. is_inside(4)) then
                                tree%blocks(idx)%i_blank(i,j,k,ibody)       = BODY_PRESENT
                                tree%blocks(idx)%i_blank_sum(i,j,k)         = BODY_PRESENT
                                ! tree%blocks(idx)%i_blank_singularity(i,j,k) = BODY_PRESENT
                            end if

                        else if (body_data(ibody)%body_select == CIRCLE) then
                            distance = sqrt((cell_pos(1) - body_data(ibody)%body_center_x)**2 &
                                          + (cell_pos(2) - body_data(ibody)%body_center_y)**2)
                            if (distance <= body_data(ibody)%body_radius) then
                                tree%blocks(idx)%i_blank(i,j,k,ibody)       = BODY_PRESENT
                                tree%blocks(idx)%i_blank_sum(i,j,k)         = BODY_PRESENT
                                ! tree%blocks(idx)%i_blank_singularity(i,j,k) = BODY_PRESENT
                            end if

                        end if
                    end do

                    ! if (tree%blocks(idx)%i_blank(i,j,k,ibody)==BODY_PRESENT) print *, "Here it is present"

                end do
            end do


            no_solid_cells = no_solid_cells + sum(tree%blocks(idx)%i_blank_sum(i_start:i_end, j_start:j_end, 1))

        end do

    end do

    no_fluid_cells = tree%no_leaf_blocks*tree%n_cell(1)*tree%n_cell(2) - no_solid_cells

end subroutine find_i_blank


!> \brief Find the ghost cells in each block IB_ghost_cells; Allocated size for IB Variables
!> \param[in,out] tree: The tree structure used.
subroutine find_IB_ghost_cell(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module, only: nbody

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present

    integer :: i, j, k
    integer :: blank_sum, max_blank_sum     ! Variables to check the number of i_blank cells around
    integer :: no_IB_ghost_cells, ghost_no  ! number of ghost cells in the block; iterating through ghost cells
    integer :: ibody

    max_blank_sum = 2**NDIM

    !> \todo Change to 2**NDIM when using IBM method
    no_interp_pts_body = 2


    do lvl = 1, tree%no_levels
        do id = 1, tree%levels(lvl)%idx_body(0)
            idx = tree%levels(lvl)%idx_body(id)
            no_IB_ghost_cells = 0
            allocate(tree%blocks(idx)%no_IB_ghost_cells(0:nbody))
            tree%blocks(idx)%no_IB_ghost_cells(0) = 0

            ! Finding the number of IB_ghost cells within the block

            do ibody = 1, nbody
                k = 1
                do j = j_start, j_end
                    do i = i_start, i_end

                        if (tree%blocks(idx)%i_blank(i,j,k,ibody)==SOLID_PRESENT) then
                            blank_sum =   tree%blocks(idx)%i_blank(i+1,j,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i,j+1,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i-1,j,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i,j-1,k,ibody)

                            if (blank_sum<max_blank_sum) then
                                no_IB_ghost_cells = no_IB_ghost_cells + 1
                            end if

                        end if

                    end do
                end do
                tree%blocks(idx)%no_IB_ghost_cells(ibody) = no_IB_ghost_cells
            end do

            ! Allocating the IB variable based on the number of IB ghost cells
            
            allocate(tree%blocks(idx)%IB_ghost_cell(no_IB_ghost_cells, 4)) ! 3 for storing the i,j,k and ibody of the cell
            allocate(tree%blocks(idx)%IP_pos_val(no_IB_ghost_cells, 3+N_VARS)) ! 7 for storing the x,y,z and u,v,p,t of IP interpolated
            allocate(tree%blocks(idx)%IP_interpol_cells(3, no_interp_pts_body, no_IB_ghost_cells)) ! 3: i,j,k of no_interp_pts_body interpolating points
            allocate(tree%blocks(idx)%dl_GC_IP(no_IB_ghost_cells)) ! Outward lengths from GC to IP

            ! Finding the IB_ghost cells
            ghost_no = 1

            do ibody = 1, nbody
                k = 1
                do j = j_start, j_end
                    do i = i_start, i_end

                        if (tree%blocks(idx)%i_blank(i,j,k,ibody)==SOLID_PRESENT) then
                            blank_sum =   tree%blocks(idx)%i_blank(i+1,j,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i,j+1,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i-1,j,k,ibody) &
                                        + tree%blocks(idx)%i_blank(i,j-1,k,ibody)

                            if (blank_sum<max_blank_sum) then
                                tree%blocks(idx)%IB_ghost_cell(ghost_no, I_g) = i
                                tree%blocks(idx)%IB_ghost_cell(ghost_no, J_g) = j
                                tree%blocks(idx)%IB_ghost_cell(ghost_no, K_g) = k
                                tree%blocks(idx)%IB_ghost_cell(ghost_no, I_BODY) = ibody

                                ghost_no = ghost_no + 1

                            end if
                        end if
                    end do
                end do
            end do        

        end do
    end do

end subroutine find_IB_ghost_cell


!> \brief Find the cells which share the same ghost cells
!> \param[in,out] tree: The tree structure used.
! subroutine find_singluarity_cells(tree)

    !     use precision_module, only: prcs_var
    !     use data_type_module
    !     use body_module, only: nbody, no_singularities, body_singularity_cells, body_singularity_dl_GC_IP

    !     implicit none
    !     type(t_tree), intent(inout) :: tree

    !     integer :: id, idx          ! Variables to iterate through the blocks present
    !     integer :: lvl              ! Variable to iterate through the levels present

    !     integer :: i, j, k
    !     integer :: blank_sum, max_blank_sum     ! Variables to check the number of i_blank cells around
    !     integer :: singular_cell_no
    !     integer :: ibody

    !     max_blank_sum = 2**NDIM

    !     no_singularities = 0

    !     ! Finding the number of singularity cells
    !     do lvl = 1, tree%no_levels
    !         do id = 1, tree%levels(lvl)%idx_body(0)
    !             idx = tree%levels(lvl)%idx_body(id)
    !             k = 1
    !             do j = j_start, j_end
    !                 do i = i_start, i_end

    !                     if (tree%blocks(idx)%i_blank_singularity(i,j,k)==SOLID_PRESENT) then
    !                         blank_sum =   tree%blocks(idx)%i_blank_singularity(i+1,j,k) &
    !                                     + tree%blocks(idx)%i_blank_singularity(i,j+1,k) &
    !                                     + tree%blocks(idx)%i_blank_singularity(i-1,j,k) &
    !                                     + tree%blocks(idx)%i_blank_singularity(i,j-1,k)

    !                         if (blank_sum<max_blank_sum-1) then
    !                             no_singularities = no_singularities + (max_blank_sum-blank_sum)
    !                         end if

    !                     end if

    !                 end do
    !             end do
    !         end do
    !     end do


    !     ! Allocating the IB variable based on the no of singularities
    !     allocate(body_singularity_cells(S_VARS, no_singularities))   !> (block, cell_x, cell_y, ghost_x, ghost_y, i_body), no_singularities


    !     ! Updating i_blank_singularity and filling body_singularity_cells variables
    !     singular_cell_no = 1
    !     do lvl = 1, tree%no_levels
    !         do id = 1, tree%levels(lvl)%idx_body(0)
    !             idx = tree%levels(lvl)%idx_body(id)
    !             do ibody = 1, nbody
    !                 k = 1
    !                 do j = j_start, j_end
    !                     do i = i_start, i_end

    !                         if (tree%blocks(idx)%i_blank(i,j,k,ibody)==SOLID_PRESENT) then
    !                             blank_sum =   tree%blocks(idx)%i_blank(i+1,j,k,ibody) &
    !                                         + tree%blocks(idx)%i_blank(i,j+1,k,ibody) &
    !                                         + tree%blocks(idx)%i_blank(i-1,j,k,ibody) &
    !                                         + tree%blocks(idx)%i_blank(i,j-1,k,ibody)

    !                             if (blank_sum<max_blank_sum-1) then
                                    
    !                                 if (tree%blocks(idx)%i_blank(i-1,j,k,ibody) == FLUID_PRESENT) then
    !                                     body_singularity_cells(S_BLOCK, singular_cell_no) = idx
    !                                     body_singularity_cells(S_CELL_I, singular_cell_no) = i-1
    !                                     body_singularity_cells(S_CELL_J, singular_cell_no) = j
    !                                     body_singularity_cells(S_CELL_K, singular_cell_no) = k
    !                                     body_singularity_cells(S_SIDE, singular_cell_no) = LEFT
    !                                     tree%blocks(idx)%i_blank_singularity(i-1,j,k) = SOLID_PRESENT
    !                                     singular_cell_no = singular_cell_no + 1
    !                                 end if

    !                                 if (tree%blocks(idx)%i_blank(i+1,j,k,ibody) == FLUID_PRESENT) then
    !                                     body_singularity_cells(S_BLOCK, singular_cell_no) = idx
    !                                     body_singularity_cells(S_CELL_I, singular_cell_no) = i+1
    !                                     body_singularity_cells(S_CELL_J, singular_cell_no) = j
    !                                     body_singularity_cells(S_CELL_K, singular_cell_no) = k
    !                                     body_singularity_cells(S_SIDE, singular_cell_no) = RIGHT
    !                                     tree%blocks(idx)%i_blank_singularity(i+1,j,k) = SOLID_PRESENT
    !                                     singular_cell_no = singular_cell_no + 1
    !                                 end if


    !                                 if (tree%blocks(idx)%i_blank(i,j-1,k,ibody) == FLUID_PRESENT) then
    !                                     body_singularity_cells(S_BLOCK, singular_cell_no) = idx
    !                                     body_singularity_cells(S_CELL_I, singular_cell_no) = i
    !                                     body_singularity_cells(S_CELL_J, singular_cell_no) = j-1
    !                                     body_singularity_cells(S_CELL_K, singular_cell_no) = k
    !                                     body_singularity_cells(S_SIDE, singular_cell_no) = BOTTOM
    !                                     tree%blocks(idx)%i_blank_singularity(i,j-1,k) = SOLID_PRESENT
    !                                     singular_cell_no = singular_cell_no + 1
    !                                 end if


    !                                 if (tree%blocks(idx)%i_blank(i,j+1,k,ibody) == FLUID_PRESENT) then
    !                                     body_singularity_cells(S_BLOCK, singular_cell_no) = idx
    !                                     body_singularity_cells(S_CELL_I, singular_cell_no) = i
    !                                     body_singularity_cells(S_CELL_J, singular_cell_no) = j+1
    !                                     body_singularity_cells(S_CELL_K, singular_cell_no) = k
    !                                     body_singularity_cells(S_SIDE, singular_cell_no) = TOP
    !                                     tree%blocks(idx)%i_blank_singularity(i,j+1,k) = SOLID_PRESENT
    !                                     singular_cell_no = singular_cell_no + 1
    !                                 end if


    !                             end if
    !                         end if

    !                     end do

    !                 end do
    !             end do
    !         end do
    !     end do

    !     ! do singular_cell_no = 1, no_singularities

    !     !     print *, body_singularity_cells(:,singular_cell_no)

    !     ! end do

! end subroutine find_singluarity_cells


! !> \brief Find the values of i_range for the i loop to iterate only in the fluid domain with sum(i_range)
! !> \param[in,out] tree: The tree structure used.
subroutine find_i_range(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module, only: nbody

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present

    integer :: i, j, k, l, skip, i1

    ! Create default values for i_range
    do lvl = 1, tree%no_levels
        do id = 1, size(tree%levels(lvl)%idx_all)
            idx = tree%levels(lvl)%idx_all(id)

            allocate(tree%blocks(idx)%i_range(0:M+2, N+2, 1))

            k = 1
            do j = j_start-ngl, j_end+ngl
                
                tree%blocks(idx)%i_range(0,j,k) = 1
                tree%blocks(idx)%i_range(1,j,k) = i_start
                tree%blocks(idx)%i_range(2,j,k) = i_end

            end do


        end do
    end do


    ! Find the i_range for the blocks which have body in it
    do lvl = 1, tree%no_levels
        do id = 1, tree%levels(lvl)%idx_body(0)
            idx = tree%levels(lvl)%idx_body(id)
            tree%blocks(idx)%i_range = 0 ! resetting the number of ranges

            k = 1
            do j = j_start, j_end
                skip = 0
                do i = i_start, i_end
                    if (skip > 0) then
                        skip = skip-1
                        cycle
                    end if
                    
                    if (sum(tree%blocks(idx)%i_blank(i,j,k,:)) == FLUID_PRESENT) then
                        tree%blocks(idx)%i_range(0,j,k) = tree%blocks(idx)%i_range(0,j,k) + 1
                        l = 2*tree%blocks(idx)%i_range(0,j,k) - 1
                        tree%blocks(idx)%i_range(l,j,k) = i
                        tree%blocks(idx)%i_range(l+1,j,k) = i_end
                        skip = i_end - i
                        do i1 = i+1,i_end+1
                            if (sum(tree%blocks(idx)%i_blank(i1,j,k,:)) /= FLUID_PRESENT) then
                                tree%blocks(idx)%i_range(l+1,j,k) = i1-1
                                skip = tree%blocks(idx)%i_range(l+1,j,k) - i
                                Exit
                            end if
                        end do
                    end if

                end do
            end do

        end do
    end do


end subroutine find_i_range


! !> \brief Find the values of i_range for the i loop to iterate only in the fluid domain with i_range_sum
! !> \param[in,out] tree: The tree structure used.
subroutine find_i_range_1(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module, only: nbody

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present

    integer :: i, j, k, l, skip, i1

    ! Create default values for i_range
    do lvl = 1, tree%no_levels
        do id = 1, size(tree%levels(lvl)%idx_all)
            idx = tree%levels(lvl)%idx_all(id)

            allocate(tree%blocks(idx)%i_range(0:M+2, N+2, 1))

            k = 1
            do j = j_start-ngl, j_end+ngl
                
                tree%blocks(idx)%i_range(0,j,k) = 1
                tree%blocks(idx)%i_range(1,j,k) = i_start
                tree%blocks(idx)%i_range(2,j,k) = i_end

            end do


        end do
    end do


    ! Find the i_range for the blocks which have body in it
    do lvl = 1, tree%no_levels
        do id = 1, tree%levels(lvl)%idx_body(0)
            idx = tree%levels(lvl)%idx_body(id)
            tree%blocks(idx)%i_range = 0 ! resetting the number of ranges

            k = 1
            do j = j_start, j_end
                skip = 0
                do i = i_start, i_end
                    if (skip > 0) then
                        skip = skip-1
                        cycle
                    end if
                    
                    if (tree%blocks(idx)%i_blank_sum(i,j,k) == FLUID_PRESENT) then
                        tree%blocks(idx)%i_range(0,j,k) = tree%blocks(idx)%i_range(0,j,k) + 1
                        l = 2*tree%blocks(idx)%i_range(0,j,k) - 1
                        tree%blocks(idx)%i_range(l,j,k) = i
                        tree%blocks(idx)%i_range(l+1,j,k) = i_end
                        skip = i_end - i
                        do i1 = i+1,i_end+1
                            if (tree%blocks(idx)%i_blank_sum(i,j,k) == SOLID_PRESENT) then
                                tree%blocks(idx)%i_range(l+1,j,k) = i1-1
                                skip = tree%blocks(idx)%i_range(l+1,j,k) - i
                                Exit
                            end if
                        end do
                    end if

                end do
            end do

        end do
    end do


end subroutine find_i_range_1

! subroutine find_i_range_singularity(tree)

    !     use precision_module, only: prcs_var
    !     use data_type_module
    !     use body_module, only: nbody

    !     implicit none
    !     type(t_tree), intent(inout) :: tree

    !     integer :: id, idx          ! Variables to iterate through the blocks present
    !     integer :: lvl              ! Variable to iterate through the levels present

    !     integer :: i, j, k, l, skip, i1

    !     ! Create default values for i_range_singularity
    !     do lvl = 1, tree%no_levels
    !         do id = 1, size(tree%levels(lvl)%idx_all)
    !             idx = tree%levels(lvl)%idx_all(id)

    !             allocate(tree%blocks(idx)%i_range_singularity(0:M+2, N+2, 1))

    !             k = 1
    !             do j = j_start-ngl, j_end+ngl
                    
    !                 tree%blocks(idx)%i_range_singularity(0,j,k) = 1
    !                 tree%blocks(idx)%i_range_singularity(1,j,k) = i_start
    !                 tree%blocks(idx)%i_range_singularity(2,j,k) = i_end

    !             end do


    !         end do
    !     end do


    !     ! Find the i_range_singularity for the blocks which have body in it
    !     do lvl = 1, tree%no_levels
    !         do id = 1, tree%levels(lvl)%idx_body(0)
    !             idx = tree%levels(lvl)%idx_body(id)
    !             tree%blocks(idx)%i_range_singularity = 0 ! resetting the number of ranges

    !             k = 1
    !             do j = j_start, j_end
    !                 skip = 0
    !                 do i = i_start, i_end
    !                     if (skip > 0) then
    !                         skip = skip-1
    !                         cycle
    !                     end if
                        
    !                     if (tree%blocks(idx)%i_blank_singularity(i,j,k) == FLUID_PRESENT) then
    !                         tree%blocks(idx)%i_range_singularity(0,j,k) = tree%blocks(idx)%i_range_singularity(0,j,k) + 1
    !                         l = 2*tree%blocks(idx)%i_range_singularity(0,j,k) - 1
    !                         tree%blocks(idx)%i_range_singularity(l,j,k) = i
    !                         tree%blocks(idx)%i_range_singularity(l+1,j,k) = i_end
    !                         skip = i_end - i
    !                         do i1 = i+1,i_end+1
    !                             if (tree%blocks(idx)%i_blank_singularity(i1,j,k) == SOLID_PRESENT) then
    !                                 tree%blocks(idx)%i_range_singularity(l+1,j,k) = i1-1
    !                                 skip = tree%blocks(idx)%i_range_singularity(l+1,j,k) - i
    !                                 Exit
    !                             end if
    !                         end do
    !                     end if

    !                 end do
    !             end do

    !         end do
    !     end do


! end subroutine find_i_range_singularity


! !> \todo (Need to be changes later which uses the intercept and projection)
! !> \brief Find the Image Point of each ghost cell at every block 
! !> \param[in,out] tree: The tree structure used.
subroutine find_Image_Point_IP(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present

    integer :: i, j, k
    integer :: blank_sum, max_blank_sum     ! Variables to check the number of i_blank cells around
    integer :: no_IB_ghost_cells, ghost_no  ! number of ghost cells in the block; iterating through ghost cells
    ! integer :: ibody

        
    do lvl = 1, tree%no_levels
        do id = 1, tree%levels(lvl)%idx_body(0)
            idx = tree%levels(lvl)%idx_body(id)

            ! do ibody = 1, nbody
            no_IB_ghost_cells = tree%blocks(idx)%no_IB_ghost_cells(nbody)
            ! print *, 'no_IB_ghost_cells', no_IB_ghost_cells

            ! Finding the image point of the blocks
            do ghost_no = 1, no_IB_ghost_cells

                i = tree%blocks(idx)%IB_ghost_cell(ghost_no, I_g) 
                j = tree%blocks(idx)%IB_ghost_cell(ghost_no, J_g) 
                k = tree%blocks(idx)%IB_ghost_cell(ghost_no, K_g) 

                ! if (idx==39 .or. idx==17)   print *, idx, i, j, k

                ! Checking the IP in the order of #1:Left  #2:Right  #3:Bottom  #4:Top
                if (sum(tree%blocks(idx)%i_blank(i-1, j, k, :)) == FLUID_PRESENT) then
                    tree%blocks(idx)%IP_pos_val(ghost_no, I_g) = i-1   ! Might have to change to actual coordinates
                    tree%blocks(idx)%IP_pos_val(ghost_no, J_g) = j
                    tree%blocks(idx)%IP_pos_val(ghost_no, K_g) = k
                    tree%blocks(idx)%dl_GC_IP(ghost_no) = dx_lvl(lvl)
                    cycle

                else if (sum(tree%blocks(idx)%i_blank(i+1, j, k, :)) == FLUID_PRESENT) then
                    tree%blocks(idx)%IP_pos_val(ghost_no, I_g) = i+1
                    tree%blocks(idx)%IP_pos_val(ghost_no, J_g) = j
                    tree%blocks(idx)%IP_pos_val(ghost_no, K_g) = k
                    tree%blocks(idx)%dl_GC_IP(ghost_no) = dx_lvl(lvl)
                    cycle

                else if (sum(tree%blocks(idx)%i_blank(i, j-1, k, :)) == FLUID_PRESENT) then
                    tree%blocks(idx)%IP_pos_val(ghost_no, I_g) = i
                    tree%blocks(idx)%IP_pos_val(ghost_no, J_g) = j-1
                    tree%blocks(idx)%IP_pos_val(ghost_no, K_g) = k
                    tree%blocks(idx)%dl_GC_IP(ghost_no) = dy_lvl(lvl)
                    cycle

                else if (sum(tree%blocks(idx)%i_blank(i, j+1, k, :)) == FLUID_PRESENT) then
                    tree%blocks(idx)%IP_pos_val(ghost_no, I_g) = i
                    tree%blocks(idx)%IP_pos_val(ghost_no, J_g) = j+1
                    tree%blocks(idx)%IP_pos_val(ghost_no, K_g) = k
                    tree%blocks(idx)%dl_GC_IP(ghost_no) = dy_lvl(lvl)
                    cycle

                else
                    print *, "======================================================"
                    print *, "      IB_Ghost cells seem to be Selected Wrong        "
                    print *, "======================================================"

                end if


            end do
            ! end do
            
        end do
    end do

end subroutine find_Image_Point_IP


! !> \todo (Need to be changed later which will scan the surroundings and find the neighbors)
! !> \brief Find the neighbouring cells from which IP will be interpolated 
! !> \param[in,out] tree: The tree structure used.
subroutine find_IP_interpol_cells(tree)

    use precision_module, only: prcs_var
    use data_type_module
    use body_module

    implicit none
    type(t_tree), intent(inout) :: tree

    integer :: id, idx          ! Variables to iterate through the blocks present
    integer :: lvl              ! Variable to iterate through the levels present

    integer :: IP_i, IP_j, IP_k, GC_i, GC_j, GC_k
    integer :: no_IB_ghost_cells, ghost_no   ! number of ghost cells in the block; iterating through ghost cells
    integer :: pt_no           ! iterate through the points
    integer :: ibody
    integer :: blank_sum, max_blank_sum

    max_blank_sum = 2**NDIM

        
    do lvl = 1, tree%no_levels
        ! do id = 1, tree%levels(lvl)%idx_body(0)
        !     idx = tree%levels(lvl)%idx_body(id)
        do id = 1, tree%levels(lvl)%idx_body_leaf(0)
            idx = tree%levels(lvl)%idx_body_leaf(id)

            ! do ibody = 1, nbody
                no_IB_ghost_cells = tree%blocks(idx)%no_IB_ghost_cells(nbody)

                ! Finding the image point of the blocks
                do ghost_no = 1, no_IB_ghost_cells

                    GC_i = tree%blocks(idx)%IB_ghost_cell(ghost_no, I_g)
                    GC_j = tree%blocks(idx)%IB_ghost_cell(ghost_no, J_g)
                    GC_k = tree%blocks(idx)%IB_ghost_cell(ghost_no, K_g)

                    IP_i = tree%blocks(idx)%IP_pos_val(ghost_no, I_g) 
                    IP_j = tree%blocks(idx)%IP_pos_val(ghost_no, J_g) 
                    IP_k = tree%blocks(idx)%IP_pos_val(ghost_no, K_g)

                    blank_sum =   tree%blocks(idx)%i_blank_sum(GC_i+1,GC_j,GC_k) &
                                + tree%blocks(idx)%i_blank_sum(GC_i,GC_j+1,GC_k) &
                                + tree%blocks(idx)%i_blank_sum(GC_i-1,GC_j,GC_k) &
                                + tree%blocks(idx)%i_blank_sum(GC_i,GC_j-1,GC_k)

                    if (blank_sum==max_blank_sum-1) then
                        pt_no = 1
                        tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = IP_i
                        tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = IP_j
                        tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = IP_k

                        pt_no = 2
                        tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = IP_i
                        tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = IP_j
                        tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = IP_k

                    else
                        pt_no = 1

                        if (tree%blocks(idx)%i_blank_sum(GC_i-1,GC_j,GC_k) == FLUID_PRESENT) then
                            tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = GC_i-1
                            tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = GC_j
                            tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = GC_k
                            ! print *, idx, pt_no, GC_i-1, GC_j, GC_k
                            pt_no = pt_no + 1
                        end if

                        if (tree%blocks(idx)%i_blank_sum(GC_i+1,GC_j,GC_k) == FLUID_PRESENT) then
                            tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = GC_i+1
                            tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = GC_j
                            tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = GC_k
                            ! print *, idx, pt_no, GC_i+1, GC_j, GC_k
                            pt_no = pt_no + 1
                        end if

                        if (tree%blocks(idx)%i_blank_sum(GC_i,GC_j-1,GC_k) == FLUID_PRESENT) then
                            tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = GC_i
                            tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = GC_j-1
                            tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = GC_k
                            ! print *, idx, pt_no, GC_i, GC_j-1, GC_k
                            pt_no = pt_no + 1
                        end if

                        if (tree%blocks(idx)%i_blank_sum(GC_i,GC_j+1,GC_k) == FLUID_PRESENT) then
                            tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = GC_i
                            tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = GC_j+1
                            tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = GC_k
                            ! print *, idx, pt_no, GC_i, GC_j+1, GC_k
                            pt_no = pt_no + 1
                        end if

                        if (pt_no > 3) then
                            print *, '-----------------------------------------------------'
                            print *, '      MORE NUMBER OF INTERPOLANTS THAN EXPECTED      '
                            print *, '-----------------------------------------------------'
                        end if


                    end if

                

                end do
            ! end do
        end do
    end do


    ! idx = 1
    ! do ghost_no = 1, no_IB_ghost_cells

    !     do pt_no = 1, 2

    !         print *, tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no), &
    !                  tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no), &
    !                  tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no)
    !     end do

    !     print *, '----------------'
    

    ! end do


end subroutine find_IP_interpol_cells

! subroutine find_IP_interpol_cells(tree)

    !     use precision_module, only: prcs_var
    !     use data_type_module
    !     use body_module

    !     implicit none
    !     type(t_tree), intent(inout) :: tree

    !     integer :: id, idx          ! Variables to iterate through the blocks present
    !     integer :: lvl              ! Variable to iterate through the levels present

    !     integer :: i, j, k
    !     integer :: no_IB_ghost_cells, ghost_no   ! number of ghost cells in the block; iterating through ghost cells
    !     integer :: pt_no           ! iterate through the points
    !     integer :: ibody

            
    !     do lvl = 1, tree%no_levels
    !         do id = 1, tree%levels(lvl)%idx_body(0)
    !             idx = tree%levels(lvl)%idx_body(id)

    !             ! do ibody = 1, nbody
    !                 no_IB_ghost_cells = tree%blocks(idx)%no_IB_ghost_cells(nbody)

    !                 ! Finding the image point of the blocks
    !                 do ghost_no = 1, no_IB_ghost_cells

    !                     i = tree%blocks(idx)%IP_pos_val(ghost_no, I_g) 
    !                     j = tree%blocks(idx)%IP_pos_val(ghost_no, J_g) 
    !                     k = tree%blocks(idx)%IP_pos_val(ghost_no, K_g) 

    !                     ! (Create logic to find the neighbouring 2**NDIM cells which are close to i, j, k)
                        
    !                     ! Creating just one neighbour which lies on the IP
    !                     pt_no = 1
    !                     tree%blocks(idx)%IP_interpol_cells(I_g, pt_no, ghost_no) = i
    !                     tree%blocks(idx)%IP_interpol_cells(J_g, pt_no, ghost_no) = j
    !                     tree%blocks(idx)%IP_interpol_cells(K_g, pt_no, ghost_no) = k

    !                 end do
    !             ! end do
    !         end do
    !     end do

! end subroutine find_IP_interpol_cells
