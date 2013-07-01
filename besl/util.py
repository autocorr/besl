"""
=========
Utilities
=========

General routines and utility functions

"""

import pandas as _pd
from catalog import read_bgps

def bgps_import_check(bgps, exten='all'):
    if len(bgps) == 0:
        # If empty return new BGPS catalog
        return read_bgps(exten=exten)
    elif isinstance(bgps, _pd.core.frame.DataFrame):
        # Return non-empty pandas DataFrame
        return bgps
    else:
        raise Exception("Not a pandas DataFrame: {0}".format(type(bgps)))


