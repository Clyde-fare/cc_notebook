#!/usr/bin/env python

import os
from setuptools import setup

#get static files we want to install
# (includes multed nested folders within jsmol so rather than list them manually we traverse the directory)
my_data = [(x[0], [x[0]+'/'+f for f in x[2]]) for x in os.walk('cc_notebook/static')]

setup(name='cc_notebook',
      version='0.1',
      packages=['cc_notebook'],
      data_files = my_data,

      entry_points = {
        'console_scripts': [
            'cc_notebook_init = cc_notebook.install:cc_notebook_init',
        ],
        }
      )

