"""
=========
Utilities
=========

General routines and utility functions

"""

import pandas as pd
from matplotlib import pyplot as plt
from tempfile import NamedTemporaryFile
from .catalog import read_bgps


def savefig(outname, close=False, dpi=300):
    """
    Save a figure to a PDF/EPS/PNG (300 dpi default).

    Parameters
    ----------
    outname : string
        The filename head to save to
    close : bool
        Close the matplotlib figure
    dpi : number
        DPI to save to
    """
    plt.savefig(outname + '.pdf')
    plt.savefig(outname + '.eps')
    plt.savefig(outname + '.png', dpi=dpi)
    if close:
        plt.close()


def bgps_import_check(bgps, exten='all'):
    if len(bgps) == 0:
        # If empty return new BGPS catalog
        return read_bgps(exten=exten)
    elif isinstance(bgps, pd.core.frame.DataFrame):
        # Return non-empty pandas DataFrame
        return bgps
    else:
        raise Exception("Not a pandas DataFrame: {0}".format(type(bgps)))


def convert_endian(df, keep_index=False):
    """
    Convert a DataFrame to the system Endian by writing to a temporary file
    and reading it back again.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to convert

    Returns
    -------
    df : pandas.DataFrame
        DataFrame converted to system endian
    """
    with NamedTemporaryFile() as temp:
        df.to_csv(temp.name, index=keep_index)
        df = pd.read_csv(temp.name)
    return df


