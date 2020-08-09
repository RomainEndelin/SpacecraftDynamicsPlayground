from sympy import latex
import jinja2
# import os
# import templates
from space_sim.latex.latex_to_pdf import compile_tex

env = jinja2.Environment(loader=jinja2.PackageLoader("space_sim.latex"))
env.globals['latex'] = latex

template = env.get_template('layout.tex.j2')

from space_sim.maths.main import (
  N, A,
  a_magnitude, compensate_y, compensate_z, a_vector,
  v0x_a, v0y_a, v0z_a, v0_vector,
  a_fn, v_fn, p_fn,
  T, R,
  distance, pT_vector,
  vTx_a, vTy_a, vTz_a, vT_vector,
  v_fn_T, p_fn_T,
  eq_matrix, result_compensate, result_R,
  v0_vector_N, vT_vector_N, p0_vector_N, pT_vector_N, direction_vector
)

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
  N=N, A=A, a_magnitude=a_magnitude, compensate_y=compensate_y,
  compensate_z=compensate_z, a_vector=a_vector, v0x_a=v0x_a, v0y_a=v0y_a,
  v0z_a=v0z_a, v0_vector=v0_vector, a_fn=a_fn, v_fn=v_fn, p_fn=p_fn, T=T, R=R,
  distance=distance, pT_vector=pT_vector, vTx_a=vTx_a, vTy_a=vTy_a,
  vTz_a=vTz_a, vT_vector=vT_vector, v_fn_T=v_fn_T, p_fn_T=p_fn_T,
  eq_matrix=eq_matrix, result_compensate=result_compensate, result_R=result_R,
  v0_vector_N=v0_vector_N, vT_vector_N=vT_vector_N, p0_vector_N=p0_vector_N,
  pT_vector_N=pT_vector_N, direction_vector=direction_vector
)

compile_tex(rendered_tex, '.')