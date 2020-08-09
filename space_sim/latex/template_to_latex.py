import os
import jinja2
import templates

latex_jinja_env = jinja2.Environment(
  block_start_string='\BLOCK{',
  block_end_string='}',
  variable_start_string='\VAR{',
  variable_end_string='}',
  comment_start_string='\#{',
  comment_end_string='}',
  line_statement_prefix='%%',
  line_comment_prefix='%#',
  trim_blocks=True,
  autoescape=False,
  loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(templates.__file__), 'latex'))
)

def print_template(template):
  print(latex_jinja_env.loader.get_source(latex_jinja_env, template)[0])