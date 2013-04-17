"""
Routines for writing proposals.
"""

import pandas as pd
from besl.coord import dec2sexstr

def create_gbt_prop_out(file_name='bgps_master1.cat', out_name='source_list'):
    """
    Make target list for online NRAO GBT Proposal Submission Tool.
    """
    df = pd.read_csv(file_name)
    obs = df[(df.nh3f == 2) & (df.hcof != 0) & (df.hco_tmb > 0.3) &
        (df.robit == 0) & (df.msx == 0) & (df.ego == 0) & ((df.glon_peak
        < 65) | (df.glon_peak > 295))]
    obs = obs[['name', 'ra_x', 'dec_x', 'hco_v']]
    obs = obs.sort(columns='ra_x')
    obs['Group Names'] = ''
    obs['Coordinate System'] = 'Equatorial'
    obs['Epoch'] = 'J2000'
    obs['Longitude'] = obs.ra_x.apply(dec2sexstr, sfigs=2, hd='h')
    obs['Latitude'] = obs.dec_x.apply(dec2sexstr, sfigs=2, hd='d')
    obs['Reference Frame'] = 'LSRK'
    obs['Convention'] = 'Radio'
    obs['Calibrator'] = 'False'
    obs['Place Holder'] = ''
    obs = obs.rename(columns={'name': 'Source Name', 'hco_v': 'Velocity'})
    obs['Velocity'] = obs['Velocity'].apply('{:0.1f}'.format)
    obs = obs[['Source Name',
                 'Group Names',
                 'Coordinate System',
                 'Epoch',
                 'Longitude',
                 'Latitude',
                 'Reference Frame',
                 'Convention',
                 'Velocity',
                 'Calibrator',
                 'Place Holder']]
    obs.to_csv(out_name + '.cat', sep=';', index=False)
    print 'File printed to ' + out_name + '.cat'
    print obs.shape
    return obs

def create_hht_prop_out(file_name='bgps_master1.cat', out_name='source_list'):
    """
    Make target list for LaTeX ARO SMT proposal.
    """
    df = pd.read_csv(file_name)
    # starless
    #obs = df[(df.hcof != 0) & (df.nnhf != 0) & (df.h2of == 0) & (df.nnh_tmb >
    #    0.3) & (df.robit == 0) & (df.msx == 0) & (df.ego == 0) & ((df.glon_peak
    #    < 65) | (df.glon_peak > 295))]
    # all nnh strong
    obs = df[(df.hcof != 0) & (df.hcof != -999) & (df.nnhf != 0) & (df.h2of !=
        2) & (df.nnh_tmb > 0.5) & ((df.glon_peak < 65) | (df.glon_peak > 295))]
    # egos
    #obs = df[(df.hcof != 0) & (df.ego != 0) & (df.ego != -999) &
    #    ((df.glon_peak < 65) | (df.glon_peak > 295))]
    obs = obs[['name', 'ra_x', 'dec_x', 'hco_v']]
    obs = obs.sort(columns='ra_x')
    obs['Sun Avoidance'] = 'none'
    obs['RA'] = obs.ra_x.apply(dec2sexstr, hd='h')
    obs['Dec'] = obs.dec_x.apply(dec2sexstr, hd='d', lead_psign=True)
    obs = obs.rename(columns={'name': 'Source Name', 'hco_v': 'Velocity'})
    obs['Velocity'] = obs['Velocity'].apply('{:0.1f}'.format)
    obs = obs[['Source Name',
               'RA',
               'Dec',
               'Velocity',
               'Sun Avoidance']]
    obs.to_csv(out_name + '.cat', sep=',', index=False)
    print 'File printed to ' + out_name + '.cat'
    print obs.shape
    return obs
