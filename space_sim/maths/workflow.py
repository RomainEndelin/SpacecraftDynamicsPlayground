from sympy import symbols
from sympy.physics.vector import ReferenceFrame

from space_sim.maths.helpers import (
    acceleration_algorithm,
    align_x_to_vector_dcm,
    make_vector,
    motion_params,
    solve_system,
)


def create_body_frame(inertial_frame, p0, pT):
    body_frame = ReferenceFrame("A")
    dcm = align_x_to_vector_dcm(inertial_frame, pT - p0)
    body_frame.orient(inertial_frame, "DCM", dcm)

    return body_frame


def initialize_unknowns():
    T = symbols("T", real=True, nonnegative=True)
    Rx, Ry, Rz = symbols("R_x R_y R_z", real=True, nonnegative=True)
    Mx, My, Mz = symbols("M_x M_y M_z", real=True, nonnegative=True)
    return {"T": T, "Rx": Rx, "Ry": Ry, "Rz": Rz, "Mx": Mx, "My": My, "Mz": Mz}


def body_motion(body_frame, pT_vector, v0_vector, vT_vector, t):
    unknowns = initialize_unknowns()

    M = symbols("M", nonnegative=True)
    a_vector = acceleration_algorithm(body_frame, **unknowns)

    a_fn, v_fn, p_fn = motion_params(body_frame, a_vector, v0_vector, t)

    solved_unknowns = solve_system(
        body_frame, p_fn, pT_vector, v_fn, vT_vector, M, t, unknowns
    )
    solve_map = {unknowns[key]: solved_unknowns[key] for key in unknowns.keys()}

    return [
        solved_unknowns["T"],
        a_fn.subs(solve_map),
        v_fn.subs(solve_map),
        p_fn.subs(solve_map),
        {"M": M},
    ]


def initialize_system():
    p0x, p0y, p0z = symbols("p^N_{0x} p^N_{0y} p^N_{0z}")
    pTx, pTy, pTz = symbols("p^N_{Tx} p^N_{Ty} p^N_{Tz}")
    v0x, v0y, v0z = symbols("v^N_{0x} v^N_{0y} v^N_{0z}")
    vTx, vTy, vTz = symbols("v^N_{Tx} v^N_{Ty} v^N_{Tz}")

    N = ReferenceFrame("N")  # Inertial frame

    p0 = make_vector(N, p0x, p0y, p0z)
    pT = make_vector(N, pTx, pTy, pTz)

    v0 = make_vector(N, v0x, v0y, v0z)
    vT = make_vector(N, vTx, vTy, vTz)

    return (
        N,
        p0,
        pT,
        v0,
        vT,
        {"p0x": p0x, "p0y": p0y, "p0z": p0z, "pTx": pTx, "pTy": pTy, "pTz": pTz},
    )


def main(t):
    N, p0, pT, v0, vT, variables = initialize_system()

    A = create_body_frame(N, p0, pT)

    direction = (pT - p0).express(A)
    v0_A = v0.express(A)
    vT_A = vT.express(A)
    T, a_fn, v_fn, p_fn, motion_variables = body_motion(A, direction, v0_A, vT_A, t)

    return (
        T,
        a_fn,
        v_fn,
        p_fn + p0,
        {"N": N, "A": A},
        {**variables, **motion_variables},
    )
