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


def cc_notebook_init(verbose=0):
    """Initilises cc_notebook

    Checks that .cc_notebook.ini exists on the home directory, if not places an empty one there

    Installs javascript dependencies"""

    default_cc_ini = ( "[pbs]",
                       "user=username",
                       "server=servername",
                       "[gaussian]",
                       "gauss_host=user@server",
                       "gauss_home=remote_homedirectory",
                       "gauss_scratch=remote_scratchdirectory",
                       "[ase]",
                       "ase_home=local_homedirectory",
                       "ase_scratch=local_scratchdirectory",
                       "[smart logs]",
                       "click=open -e",
                       "shift_click=prog1",
                       "cntrl_click=prog2",
                       "cntrl_shift_click=prog3" )

    cc_ini_fl_pth = os.path.expanduser('~/.cc_notebook.ini')

    if not os.path.isfile(cc_ini_fl_pth):
        with open(cc_ini_fl_pth, 'w') as cc_ini:
            cc_ini.write('\n'.join(default_cc_ini))

        print('\n{f} constructed please fill in details\n'.format(f=cc_ini_fl_pth))

    else:
        print('\n{f} already exists not overwriting\n'.format(f=cc_ini_fl_pth))

def enable_notebook(verbose=0):
    """Enables notebook widgets to be displayed.

    This function should be called before using cc_notebook.
    """
    cc_ini_fl_pth = os.path.expanduser('~/.cc_notebook.ini')
    if not os.path.isfile(cc_ini_fl_pth):
        raise(RuntimeError('cc_notebook not initialised please run cc_notebook_init from the terminal'))

    libs = ['cc_notebook.js',
            'jsmol']
    fns = [resource_filename('cc_notebook', os.path.join('static', f)) for f in libs]
    install_nbextension([fns[0]], verbose=verbose, overwrite=True)
    
    # jsmol is large so we won't overwrite
    install_nbextension(fns[1:], verbose=verbose, overwrite=False)
    display(activate_nb)

if __name__ == '__main__':
    cc_notebook_init()