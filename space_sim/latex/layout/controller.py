import jinja2
from sympy import latex, symbols

# import os
# import templates
from space_sim.latex.latex_to_pdf import compile_tex
from space_sim.maths.helpers import (
    eval_in_T,
    make_equation_matrix,
    motion_params,
    solve_R_T,
)
from space_sim.maths.workflow import (
    create_body_frame,
    initialize_motion_vectors,
    initialize_system,
    initialize_unknowns,
)

env = jinja2.Environment(loader=jinja2.PackageLoader("space_sim.latex"))
env.globals["latex"] = latex

template = env.get_template("layout.tex.j2")

t = symbols("t", real=True, nonnegative=True)
N, p0, pT, _ = initialize_system()
A = create_body_frame(N, p0, pT)
direction = (pT - p0).express(A)

T, R = initialize_unknowns()
a_vector, v0_vector, vT_vector, variables = initialize_motion_vectors(A, T, R)
a_fn, v_fn, p_fn = motion_params(A, a_vector, v0_vector, t)
distance = symbols("|d|")

solved_R, solved_T = solve_R_T(A, p_fn, direction, v_fn, vT_vector, t, R, T)
# These variables are intermediate variables from solve_R_T
computed_v_T = eval_in_T(v_fn, t, T, R)
computed_p_T = eval_in_T(p_fn, t, T, R)
eq_matrix = make_equation_matrix(A, computed_p_T, computed_v_T, direction, vT_vector)


# latex_jinja_env = jinja2.Environment(
#   block_start_string='\BLOCK{',
#   block_end_string='}',
#   variable_start_string='\VAR{',
#   variable_end_string='}',
#   comment_start_string='\#{',
#   comment_end_string='}',
#   line_statement_prefix='%%',
#   line_comment_prefix='%#',
#   trim_blocks=True,
#   autoescape=False,
#   loader =  PackageLoader('space_sim.latex.layout')
#   # loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(templates.__file__), 'latex'))
# )

# def print_template(template):
#   print(latex_jinja_env.loader.get_source(latex_jinja_env, template)[0])

rendered_tex = template.render(
    N=N,
    A=A,
    a_magnitude=variables["a_magnitude"],
    compensate_y=0,
    compensate_z=0,
    a_vector=a_vector,
    v0x_a=variables["v0_magnitude"],
    v0y_a=0,
    v0z_a=0,
    v0_vector=v0_vector,
    a_fn=a_fn,
    v_fn=v_fn,
    p_fn=p_fn,
    T=solved_T,
    R=R,
    distance=distance,
    pT_vector=pT.express(A),
    vTx_a=variables["vT_magnitude"],
    vTy_a=0,
    vTz_a=0,
    vT_vector=vT_vector,
    v_fn_T=computed_v_T,
    p_fn_T=computed_p_T,
    eq_matrix=eq_matrix,
    result_compensate=0,
    result_R=solved_R,
    v0_vector_N=v0_vector.express(A),
    vT_vector_N=vT_vector.express(A),
    p0_vector_N=p0,
    pT_vector_N=pT,
    direction_vector=direction,
)

compile_tex(rendered_tex, ".")
