import subprocess
import tempfile
import papermill as pm
import os.path as path

def _exec_notebook(path):
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=1000",
                "--ExecutePreprocessor.kernel_name=python3",
                "--output", fout.name, path]
        subprocess.check_call(args)


def _exec_papermill(input_nb, args):
    output_nb = path.join('test',input_nb)
    pm.execute_notebook(input_nb, output_nb, parameters=args)


def test():
    print('Testing Jupyter notebooks...')
    _exec_notebook('Fig01_Diagram.ipynb')
    #_exec_notebook('Fig02_Phonon-Pulse-Template.ipynb') #Comment out incomplete notebook.
    _exec_notebook('Fig03_Calibration-Lines.ipynb')
    _exec_notebook('Fig04_Detector-Resolution-Model.ipynb')
    _exec_notebook('Fig05_dN-Distribution-Fit.ipynb')
    _exec_notebook('Fig06_Estimated-Cut-Efficiencies.ipynb')
    _exec_notebook('Fig07_Trigger-Write-Efficiencies.ipynb')
    _exec_notebook('Fig08_Energy-Yield.ipynb')
    _exec_notebook('Fig09_Measured-Spectra.ipynb')
    _exec_notebook('Fig10_Distribution-Yield-Curves-and-Components')
    _exec_notebook('Fig11_Example-Corner-Plot.ipynb')
    _exec_notebook('Fig12_Example-Parameter-MCMC-Chains.ipynb')
    _exec_notebook('Fig13_Range-BestFit.ipynb')
