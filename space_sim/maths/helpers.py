from functools import partial

from sympy import Matrix, Max, Min, Piecewise, eye, solveset
from sympy.physics.vector import dynamicsymbols, get_motion_params

from space_sim.maths.constants import G

# TODO: Break it down. For now, this file is just here to help breaking down
# the code into functions


## Orienting A relative to N
def make_rotation_matrix(origin_vector, destination_vector):
    v = origin_vector.cross(destination_vector)
    c = origin_vector.dot(destination_vector)  # cosine of angle
    K = Matrix(3, 3, [0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0])
    return eye(3) + K + (K * K) * (1 - c) / (1 - c ** 2)


def align_x_to_vector_dcm(frame, vector):
    unit_vector = vector.normalize()
    Nx = frame.x.to_matrix(frame)
    NA1 = unit_vector.to_matrix(frame)

    rotation_matrix = make_rotation_matrix(Nx, NA1)
    NA2 = (rotation_matrix * frame.y.to_matrix(frame)).simplify()

    dcm = Matrix(3, 3, [*NA1, *NA2, *NA1.cross(NA2)])
    # TODO: Why do I need to transpose?
    return dcm.transpose()


## Motion
def acceleration_algorithm(frame, G_multiplier, T, R):
    return (
        frame.x
        * Piecewise((-1, T * R - dynamicsymbols._t < 0), (1, True))
        * G_multiplier
        * G
    )


def motion_params(frame, a, v0, t):
    original_a_fn, original_v_fn, original_p_fn = get_motion_params(
        frame, acceleration=a, velocity=v0
    )

    a_fn = original_a_fn.subs(dynamicsymbols._t, t)
    v_fn = original_v_fn.subs(dynamicsymbols._t, t).simplify()
    p_fn = (
        original_p_fn.subs(dynamicsymbols._t, t).applyfunc(replace_max_0_RT).simplify()
    )
    return a_fn, v_fn, p_fn


def solve_R_T(frame, p_t, p_T, v_t, v_T, t, R, T):
    computed_v_T = v_t.subs(t, T).applyfunc(partial(replace_min_T, R=R, T=T))
    computed_p_T = p_t.subs(t, T).applyfunc(partial(replace_min_T, R=R, T=T))

    eq_matrix = (
        v_T.to_matrix(frame).col_join(p_T.to_matrix(frame))
        - computed_v_T.to_matrix(frame).col_join(computed_p_T.to_matrix(frame))
    ).simplify()

    ## Solving R

    result_R = solveset(eq_matrix[0], R)
    solved_R = list(result_R)[0]

    ## Solving T

    intermediate_matrix = eq_matrix.subs({R: solved_R}).simplify()

    solved_T1, solved_T2 = solveset(intermediate_matrix[3], T)
    solved_T = Max(solved_T1, solved_T2)
    return solved_R, solved_T


## Utilities
def replace_max_0_RT(f):
    # Given that T >= 0 and R in [0, 1]:
    #  We can say that max(0, R * T) = R * T
    #  However sympy is not able to infer this, so we replace manually
    return f.replace(Max, lambda a, b: (b if a == 0 else a))


def replace_min_T(f, R, T):
    # Given that T >= 0 and R in [0, 1]:
    #  We can say that min(t, R * T) = R * T
    #  However sympy is not able to infer this, so we replace manually
    return f.replace(Min(T, R * T), R * T)
