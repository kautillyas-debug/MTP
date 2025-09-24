#============ COMPILER ============
FC = ifx

#============ COMPILER FLAG ============
CFLAG = -O3 ## Optimisation flag
CFLAG += -qmkl ## Intel Math Kernel Library
# CFLAG += -g ## Profiling Flag
# CFLAG += -qopenmp ## OpenMP Parallel Processing Flag
# CFLAG += -check bounds -traceback -g ## Error Debuggings
#------------------------------------------------------
CFLAG += -DDEBUG -DSTEP_SKIP=1 ## Debugging Flag
# CFLAG += -DVERBOSE ## Detailed Information
#------------------------------------------------------
CFLAG += -DSQ_CELL ## Sides of Cell are equal dx=dy
# CFLAG += -DNDIM=2 ## Dimension of Problem (Yet to be implemented)
#------------------------------------------------------

#============ TARGET ============
exe = binary

#============ OBJECTS ============
obj =   BBM_modules.o                   	\
		Routine_calculate_rhs.o				\
		Routine_iterative_solver.o			\
		Routine_boundary_conditions.o		\
		Routine_face_velocities.o			\
		Routine_print_output.o				\
		Routine_calculate_residual.o		\
		Routine_ghost_cell_update.o			\
		Routine_read_init.o					\
		Routine_SOLVE_EQS.o 				\
		Routine_courant_compare.o 			\
		Routine_body_init.o 				\
		Routine_body_interpolate.o			\
		Routine_body_ghost_cell_update.o	\
		Routine_corrections.o				\
		BBM_main.o				

#============ CREATING TARGET ============
${exe}: ${obj} makefile definitions.h
	${FC} ${CFLAG} -o ${exe} ${obj}

#============ CREATING OBJECTS ============
%.o: %.F90 makefile definitions.h
	${FC} ${CFLAG} -c  $< 

#============ CLEANING FILES ============
clean:
	rm -f *.o *.mod ${exe}