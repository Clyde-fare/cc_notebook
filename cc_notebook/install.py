# Adapted from mdtraj: github.com/mdtraj/mdtraj

import os
import warnings
from IPython.display import display, Javascript
from IPython.html.nbextensions import install_nbextension
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from pkg_resources import resource_filename

__all__ = ['enable_notebook']

activate_nb = Javascript('''
require(['nbextensions/cc_notebook'])
require(['nbextensions/jsmol/JSmol.min.nojq'])
''')

def enable_notebook(verbose=0):
    """Enable IPython notebook widgets to be displayed.

    This function should be called before using cc_notebook.
    """
    libs = ['cc_notebook.js', 
            'jsmol']
    fns = [resource_filename('cc_notebook', os.path.join('static', f)) for f in libs]
    install_nbextension(fns, verbose=verbose, overwrite=False)
    display(activate_nb)
