"""
====
besl
====

.. moduleauthor :: Brian Svoboda <svobodb@email.arizona.edu>

"""

import os
from ConfigParser import ConfigParser


# Config files
config_name = 'besl.cfg'
local_filen = os.path.expanduser(os.path.join('~', '.' + config_name))
root_filen = os.path.abspath(os.path.join(
    os.getcwd(), os.path.pardir, config_name))
# Parse config
cfg = ConfigParser()
cfg.read([root_filen, local_filen])


__all__ = ['bplot', 'catalog', 'coord', 'image', 'math', 'misc', 'util']

__author__ = ['Brian Svoboda, svobodb@email.arizona.edu']
__version__ = '0.0.2'
