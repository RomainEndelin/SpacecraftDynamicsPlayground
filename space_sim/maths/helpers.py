from functools import partial

from sympy import Matrix, Max, Min, Piecewise, eye, solveset
from sympy.physics.vector import dynamicsymbols, get_motion_params

from space_sim.maths.constants import G

# TODO: Break it down. For now, this file is just here to help breaking down
# the code into functions


def make_vector(frame, x, y, z):
    return frame.x * x + frame.y * y + frame.z * z


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
def acceleration_algorithm(frame, T, Rx, Ry, Rz, Mx, My, Mz):
    return (
        frame.x
        * Piecewise(
            (-1, T * Rx - dynamicsymbols._t < 0),
            (1, T * Rx - dynamicsymbols._t > 0),
            (0, True),
        )
        * Mx
        * G
        + frame.y
        * Piecewise(
            (-1, T * Ry - dynamicsymbols._t < 0),
            (1, T * Ry - dynamicsymbols._t > 0),
            (0, True),
        )
        * My
        * G
        + frame.z
        * Piecewise(
            (-1, T * Rz - dynamicsymbols._t < 0),
            (1, T * Rz - dynamicsymbols._t > 0),
            (0, True),
        )
        * Mz
        * G
    )


def acceleration_algorithm2(frame, M, T, Rx, Ry, Rz, Mx, My, Mz):
    return (
        (
            frame.x * Piecewise((-1, T * Rx - dynamicsymbols._t < 0), (1, True)) * Mx
            + frame.y * Piecewise((-1, T * Ry - dynamicsymbols._t < 0), (1, True)) * My
            + frame.z * Piecewise((-1, T * Rz - dynamicsymbols._t < 0), (1, True)) * Mz
        ).normalize()
        * M
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


def eval_in_T(vector, t, T, Rx, Ry, Rz):
    return (
        vector.subs(t, T)
        .applyfunc(partial(replace_min_T, R=Rx, T=T))
        .applyfunc(partial(replace_min_T, R=Ry, T=T))
        .applyfunc(partial(replace_min_T, R=Rz, T=T))
    )


def solve_RT_1d(axis, computed_pT, computed_vT, expected_pT, expected_vT, R, T):
    result_R = solveset((computed_vT - expected_vT).dot(axis), R)
    solved_R = list(result_R)[0]

    result_T = solveset((computed_pT - expected_pT).dot(axis).subs({R: solved_R}), T)
    # TODO: Why is it two solutions?
    solved_T = list(result_T)[1]
    return solved_R, solved_T


def make_equation_matrix(
    frame, computed_pT, computed_vT, computed_M, expected_pT, expected_vT, expected_M
):
    return (
        expected_vT.to_matrix(frame)
        .col_join(expected_pT.to_matrix(frame))
        .col_join(expected_M)
        - computed_vT.to_matrix(frame)
        .col_join(computed_pT.to_matrix(frame))
        .col_join(computed_M)
    ).simplify()


def make_equation_matrix2(frame, computed_pT, computed_vT, expected_pT, expected_vT):
    return (
        expected_vT.to_matrix(frame).col_join(expected_pT.to_matrix(frame))
        - computed_vT.to_matrix(frame).col_join(computed_pT.to_matrix(frame))
    ).simplify()


def solve_system(frame, p_t, p_T, v_t, v_T, M, t, unknowns):
    computed_v_T = eval_in_T(
        v_t, t, unknowns["T"], unknowns["Rx"], unknowns["Ry"], unknowns["Rz"]
    )
    computed_p_T = eval_in_T(
        p_t, t, unknowns["T"], unknowns["Rx"], unknowns["Ry"], unknowns["Rz"]
    )
    computed_M = (
        unknowns["Mx"] * frame.x + unknowns["My"] * frame.y + unknowns["My"] * frame.z
    ).magnitude()

    eq_matrix = make_equation_matrix(
        frame, computed_p_T, computed_v_T, computed_M, p_T, v_T, M
    )

    ## Solving R

    solved_Rx = list(solveset(eq_matrix[0], unknowns["Rx"]))[0]
    solved_Ry = list(solveset(eq_matrix[1], unknowns["Ry"]))[0]
    solved_Rz = list(solveset(eq_matrix[2], unknowns["Rz"]))[0]

    ## Solving M

    intermediate_matrix = eq_matrix.subs(
        {
            unknowns["Rx"]: solved_Rx,
            unknowns["Ry"]: solved_Ry,
            unknowns["Rz"]: solved_Rz,
        }
    ).simplify()

    solved_Mx = list(solveset(intermediate_matrix[0], unknowns["Mx"]))[0]
    solved_My = list(solveset(intermediate_matrix[1], unknowns["My"]))[0]
    solved_Mz = list(solveset(intermediate_matrix[2], unknowns["Mz"]))[0]

    ## Solving T

    intermediate_matrix2 = intermediate_matrix.subs(
        {
            unknowns["Mx"]: solved_Mx,
            unknowns["My"]: solved_My,
            unknowns["Mz"]: solved_Mz,
        }
    ).simplify()

    solved_T1, solved_T2 = solveset(intermediate_matrix2[3], unknowns["T"])
    solved_T = Max(solved_T1, solved_T2)
    return {
        "Rx": solved_Rx,
        "Ry": solved_Ry,
        "Rz": solved_Rz,
        "Mx": solved_Mx,
        "My": solved_My,
        "Mz": solved_Mz,
        "T": solved_T,
    }


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
