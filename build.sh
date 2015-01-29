#!/bin/bash

if [ ! -f ~/.cc_notebook.ini ]; then
   cp .cc_notebook.ini ~
fi
$PYTHON setup.py install

# Add more build steps here, if they are necessary.
# See
# https://github.com/ContinuumIO/conda/blob/master/conda/builder/README.txt

