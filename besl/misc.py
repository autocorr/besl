"""
Authors: Brian Svoboda
Date: 2012, 12, 10

Miscellaneous routines.
"""

def logit():
    """
    Log IPython session to log file tagged by date and time:
    ipython_log_YY-MM-DD_HH:MM.py.
    """
    import os
    import time
    from IPython.core.interactiveshell import InteractiveShell
    cwd = os.getcwd()
    t = time.localtime()
    log_filen = '.ipython_log_{year}-{mon}-{day}_{hour}:{minu}.py'.format(
            year=str(t.tm_year)[2:],
            mon=t.tm_mon,
            day=t.tm_mday,
            hour=t.tm_hour,
            minu=t.tm_min)
    get_ipython().magic(u'logstart {}'.format(log_filen))
    return

