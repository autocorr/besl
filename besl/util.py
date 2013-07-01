"""
=========
Utilities
=========

General routines and utility functions

"""

import pandas as _pd
import catalog

def bgps_import_check(bgps, exten='all'):
    if len(bgps) == 0:
        catalog.read_bgps(exten=exten)
        return bgps
    elif isinstance(bgps, _pd.core.frame.DataFrame):
        return bgps
    else:
        raise Exception("Not a pandas DataFrame: {0}".format(type(bgps)))


