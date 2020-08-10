import os
import shutil
import tempfile
from subprocess import Popen


def compile_tex(rendered_tex, out_pdf_path):
    tmp_dir = tempfile.mkdtemp()
    in_tmp_path = os.path.join(tmp_dir, "rendered.tex")
    with open(in_tmp_path, "w") as outfile:
        outfile.write(rendered_tex)
    out_tmp_path = os.path.join(tmp_dir, "rendered.pdf")
    p = Popen(["pdflatex", "-output-directory", tmp_dir, in_tmp_path])
    p.communicate()
    shutil.copy(out_tmp_path, out_pdf_path)
    shutil.rmtree(tmp_dir)
