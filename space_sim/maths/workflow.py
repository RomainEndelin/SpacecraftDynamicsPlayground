from sympy import symbols
from sympy.physics.vector import ReferenceFrame

from space_sim.maths.helpers import (
    acceleration_algorithm,
    align_x_to_vector_dcm,
    make_vector,
    motion_params,
    solve_R_T,
)


def create_body_frame(inertial_frame, p0, pT):
    body_frame = ReferenceFrame("A")
    dcm = align_x_to_vector_dcm(inertial_frame, pT - p0)
    body_frame.orient(inertial_frame, "DCM", dcm)

    return body_frame


def initialize_unknowns():
    T = symbols("T", real=True, nonnegative=True)
    R = symbols("R", real=True, nonnegative=True)
    return T, R


def initialize_motion_vectors(body_frame, T, R):
    a_magnitude = symbols("|a|", nonnegative=True)
    a_vector = acceleration_algorithm(body_frame, a_magnitude, T, R)
    v0_magnitude, vT_magnitude = symbols("|v_0| |v_T|")
    v0_vector = body_frame.x * v0_magnitude
    vT_vector = body_frame.x * vT_magnitude
    return (
        a_vector,
        v0_vector,
        vT_vector,
        {
            "a_magnitude": a_magnitude,
            "v0_magnitude": v0_magnitude,
            "vT_magnitude": vT_magnitude,
        },
    )


def body_motion(body_frame, pT_vector, t):
    T, R = initialize_unknowns()

    a_vector, v0_vector, vT_vector, variables = initialize_motion_vectors(
        body_frame, T, R
    )

    a_fn, v_fn, p_fn = motion_params(body_frame, a_vector, v0_vector, t)

    solved_R, solved_T = solve_R_T(
        body_frame, p_fn, pT_vector, v_fn, vT_vector, t, R, T
    )

    return [
        solved_T,
        a_fn.subs({R: solved_R, T: solved_T}),
        v_fn.subs({R: solved_R, T: solved_T}),
        p_fn.subs({R: solved_R, T: solved_T}),
        variables,
    ]


def initialize_system():
    p0x, p0y, p0z = symbols("p^N_{0x} p^N_{0y} p^N_{0z}")
    pTx, pTy, pTz = symbols("p^N_{Tx} p^N_{Ty} p^N_{Tz}")

    N = ReferenceFrame("N")  # Inertial frame

    p0 = make_vector(N, p0x, p0y, p0z)
    pT = make_vector(N, pTx, pTy, pTz)

    return (
        N,
        p0,
        pT,
        {"p0x": p0x, "p0y": p0y, "p0z": p0z, "pTx": pTx, "pTy": pTy, "pTz": pTz},
    )


def main(t):
    N, p0, pT, variables = initialize_system()

    A = create_body_frame(N, p0, pT)

    direction = (pT - p0).express(A)
    T, a_fn, v_fn, p_fn, motion_variables = body_motion(A, direction, t,)

    return (
        T,
        a_fn,
        v_fn,
        p_fn + p0,
        {"N": N, "A": A},
        {**variables, **motion_variables},
    )
