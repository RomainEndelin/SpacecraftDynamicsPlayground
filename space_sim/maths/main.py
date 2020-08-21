from sympy import symbols
from sympy.physics.vector import ReferenceFrame

from space_sim.maths.helpers import (
    acceleration_algorithm,
    align_x_to_vector_dcm,
    motion_params,
    solve_R_T,
)

## Reference frames
N = ReferenceFrame("N")  # Inertial frame
A = ReferenceFrame("A")

## Reference variables
T = symbols("T", real=True, nonnegative=True)
R = symbols("R", real=True, nonnegative=True)
t = symbols("t", real=True, nonnegative=True)

## Acceleration
a_magnitude = symbols("|a|", nonnegative=True)
a_vector = acceleration_algorithm(A, a_magnitude, T, R)

## Initial velocity
v0_magnitude = symbols("|v_0|")
v0_vector = A.x * v0_magnitude

## Motion params
a_fn, v_fn, p_fn = motion_params(A, a_vector, v0_vector, t)

## Final position
distance = symbols("|d|")
pT_vector = A.x * distance

## Final velocity
vT_magnitude = symbols("|v_T|")
vT_vector = A.x * vT_magnitude


## Equation set
solved_R, solved_T = solve_R_T(A, p_fn, pT_vector, v_fn, vT_vector, t, R, T)

## Coordinates in N

p0x_n, p0y_n, p0z_n = symbols("p^N_{0x} p^N_{0y} p^N_{0z}")
pTx_n, pTy_n, pTz_n = symbols("p^N_{Tx} p^N_{Ty} p^N_{Tz}")
p0_vector_N = N.x * p0x_n + N.y * p0y_n + N.z * p0z_n
pT_vector_N = N.x * pTx_n + N.y * pTy_n + N.z * pTz_n

## Orienting A relative to N

direction_vector = pT_vector_N - p0_vector_N
dcm = align_x_to_vector_dcm(N, direction_vector)

A.orient(N, "DCM", dcm)
