#if ~/cc_notebook.ini does not exist importing will cause a runtime error
# the initialisation exit point that generate ~/.cc_notebook.ini lives inside install.py
# in order for it to execute we can't break here.
#TODO throw a more specific exception than RuntimeError when there is no .ini file though
#TODO make sure additional calls to functions/methods that require ~/.cc_notebook.ini all throw exceptions

try:
    from cc_notebook_utils import *
except RuntimeError:
    pass