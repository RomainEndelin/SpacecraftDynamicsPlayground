from sympy import Matrix, Max, Min, Piecewise, eye, nonlinsolve, solveset, symbols
from sympy.physics.units import gravitational_constant as G
from sympy.physics.vector import ReferenceFrame, dynamicsymbols, get_motion_params

## Reference frames
N = ReferenceFrame("N")  # Inertial frame
A = ReferenceFrame("A")

## Reference variables
T = symbols("T", real=True, nonnegative=True)
R = symbols("R", real=True, nonnegative=True)
t = symbols("t", real=True, nonnegative=True)

## Acceleration
# F stands for "force", I need to find a better name
compensate_y, compensate_z = symbols("F_y F_z")
a_magnitude = symbols("|a|", nonnegative=True)

a_vector = (
    (
        A.x * Piecewise((-1, T * R - dynamicsymbols._t < 0), (1, True))
        + A.y * compensate_y
        + A.z * compensate_z
    ).normalize()
    * a_magnitude
    * G
)

## Initial velocity
v0x_a, v0y_a, v0z_a = symbols("v^A_{0x} v^A_{0y} v^A_{0z}")
v0_vector = A.x * v0x_a + A.y * v0y_a + A.z * v0z_a


## Motion params
def replace_max_0_RT(f):
    # Given that T >= 0 and R in [0, 1]:
    #  We can say that max(0, R * T) = R * T
    #  However sympy is not able to infer this, so we replace manually
    return f.replace(Max, lambda a, b: (b if a == 0 else a))


original_a_fn, original_v_fn, original_p_fn = get_motion_params(
    A, acceleration=a_vector, velocity=v0_vector
)

a_fn = original_a_fn.subs(dynamicsymbols._t, t)
v_fn = original_v_fn.subs(dynamicsymbols._t, t).simplify()
p_fn = original_p_fn.subs(dynamicsymbols._t, t).applyfunc(replace_max_0_RT).simplify()

## Final position
distance = symbols("|d|")
pT_vector = A.x * distance

## Final velocity
vTx_a, vTy_a, vTz_a = symbols("v^A_{Tx} v^A_{Ty} v^A_{Tz}")
vT_vector = A.x * vTx_a + A.y * vTy_a + A.z * vTz_a


## Equation set
def replace_min_T(f):
    # Given that T >= 0 and R in [0, 1]:
    #  We can say that min(t, R * T) = R * T
    #  However sympy is not able to infer this, so we replace manually
    return f.replace(Min(T, R * T), R * T)


# TODO: Re-introduce v0_y and vT_y
v_fn_T = v_fn.subs(t, T).applyfunc(
    replace_min_T
)  # .subs({v0y_a: 0, v0z_a: 0, vTy_a: 0, vTz_a: 0})
p_fn_T = p_fn.subs(t, T).applyfunc(
    replace_min_T
)  # .subs({v0y_a: 0, v0z_a: 0, vTy_a: 0, vTz_a: 0})

eq_matrix = (
    vT_vector.to_matrix(A).col_join(pT_vector.to_matrix(A))
    - v_fn_T.to_matrix(A).col_join(p_fn_T.to_matrix(A))
).simplify()

# nonlinsolve(
#     [eq_matrix[1] - eq_matrix[4], eq_matrix[2] - eq_matrix[5]],
#     [compensate_y, compensate_z]
# )
# from sympy import S
# solved_Fy = solveset(eq_matrix[1] - eq_matrix[4], compensate_y, domain=S.Reals)
# solved_Fy
# tmp = eq_matrix.subs({compensate_y: solved_Fy})
# solveset(tmp[2] - tmp[5], compensate_z)

## Solving Fy and Fz
result_compensate = nonlinsolve(eq_matrix, [compensate_y, compensate_z])

# TODO: Turn it into a singleton
solved_Fy, solved_Fz = list(result_compensate)[0]

## Solving R

result_R = solveset(
    eq_matrix[0].subs({compensate_y: solved_Fy, compensate_z: solved_Fz}), R
)

solved_R = list(result_R)[0]

## Solving T

intermediate_matrix = eq_matrix.subs(
    {compensate_y: solved_Fy, compensate_z: solved_Fz, R: solved_R}
).simplify()

# result_T = nonlinsolve(intermediate_matrix, [T])
result_T = solveset(intermediate_matrix[3], T)
# solved_T = list(result_T)[0]

# FIXME: We expect intermediate_matrix[1] and intermediate_matrix[2] to solve depending of T

## Coordinates in N

p0x_n, p0y_n, p0z_n = symbols("p^N_{0x} p^N_{0y} p^N_{0z}")
pTx_n, pTy_n, pTz_n = symbols("p^N_{Tx} p^N_{Ty} p^N_{Tz}")
p0_vector_N = N.x * p0x_n + N.y * p0y_n + N.z * p0z_n
pT_vector_N = N.x * pTx_n + N.y * pTy_n + N.z * pTz_n

v0x_n, v0y_n, v0z_n = symbols("v^N_{0x} v^N_{0y} v^N_{0z}")
vTx_n, vTy_n, vTz_n = symbols("v^N_{Tx} v^N_{Ty} v^N_{Tz}")
v0_vector_N = N.x * v0x_n + N.y * v0y_n + N.z * v0z_n
vT_vector_N = N.x * vTx_n + N.y * vTy_n + N.z * vTz_n

## Orienting A relative to N

direction_vector = (pT_vector_N - p0_vector_N).normalize()

Nx = N.x.to_matrix(N)
NA1 = direction_vector.to_matrix(N)

v = Nx.cross(NA1)
c = Nx.dot(NA1)  # cosine of angle
K = Matrix(3, 3, [0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0])
rotation_matrix = eye(3) + K + (K * K) * (1 - c) / (1 - c ** 2)

NA2 = (rotation_matrix * N.y.to_matrix(N)).simplify()

dcm = Matrix(3, 3, [*NA1, *NA2, *NA1.cross(NA2)])

A.orient(N, "DCM", dcm.transpose())
# TODO: Why do I need to transpose?
