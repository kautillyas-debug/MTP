
#define RECTANGLE  0
#define CIRCLE     1

#define BODY_PRESENT   1
#define SOLID_PRESENT  1
#define FLUID_PRESENT  0

#define I_g    1
#define J_g    2
#define K_g    3
#define I_BODY 4

#define LEFT_BOTTOM -1
#define RIGHT_TOP    1

#define DIRICHLET 1
#define NEUMANN   2

#define N_VARS      5
#define U_VELOCITY  1
#define V_VELOCITY  2
#define PRESSURE    3
#define TEMPERATURE 4
#define DIRECTION   5

#define N_DIRS   4
#define LEFT    1
#define TOP     2
#define RIGHT   3
#define BOTTOM  4


#define GAUSS_SEIDEL       1
#define WEIGHTED_JACOBI    2
#define CONJUGATE_GRADIENT 3
#define MULTI_GRID         4

#define S_VARS   5
#define S_BLOCK  1
#define S_CELL_I 2
#define S_CELL_J 3
#define S_CELL_K 4
#define S_SIDE   5