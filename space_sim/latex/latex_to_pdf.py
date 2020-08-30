import os
import shutil
from subprocess import Popen


def compile_tex(rendered_tex, out_pdf_path):
    latex_build_dir = os.path.dirname(os.path.join("builds", "latex"))
    in_tmp_path = os.path.join(latex_build_dir, "assets", "rendered.tex")
    with open(in_tmp_path, "w") as outfile:
        outfile.write(rendered_tex)
    out_build_path = os.path.join(latex_build_dir, "build")
    out_build_file = os.path.join(out_build_path, "rendered.pdf")
    p = Popen(["pdflatex", "-output-directory", out_build_file, in_tmp_path])
    p.communicate()
    shutil.copy(out_build_file, out_pdf_path)
