"""
=================
Proposal Routines
=================

Routines for writing proposals.

"""

import pandas as pd
from .catalog import read_cat
from .coord import (dec2sexstr, pd_gal2eq, pd_eq2gal, sep_coords)


def create_gbt_prop_out(obs, out_name='source_list', ra_col='ra',
        dec_col='dec', velo_col='mol_hco_v'):
    """
    Make target list for online NRAO GBT Proposal Submission Tool.
    """
    #obs = df[(df.nh3f == 2) & (df.hcof != 0) & (df.hco_tmb > 0.3) &
    #    (df.robit == 0) & (df.msx == 0) & (df.ego == 0) & ((df.glon_peak
    #    < 65) | (df.glon_peak > 295))]
    obs = obs[['name', ra_col, dec_col, velo_col]]
    obs = obs.sort(columns=ra_col)
    obs['Group Names'] = ''
    obs['Coordinate System'] = 'Equatorial'
    obs['Epoch'] = 'J2000'
    obs['Longitude'] = obs[ra_col].apply(dec2sexstr, sfigs=2, hd='h')
    obs['Latitude'] = obs[dec_col].apply(dec2sexstr, sfigs=2, hd='d')
    obs['Reference Frame'] = 'LSRK'
    obs['Convention'] = 'Radio'
    obs['Calibrator'] = 'False'
    obs['Place Holder'] = ''
    obs = obs.rename(columns={'name': 'Source Name', velo_col: 'Velocity'})
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


def create_hht_prop_out(obs, out_name='source_list', ra_col='ra',
        dec_col='dec', velo_col='mol_hco_v'):
    """
    Make target list for LaTeX ARO SMT proposal.
    """
    obs = obs[['name', ra_col, dec_col, velo_col]]
    obs = obs.sort(columns=ra_col)
    obs['Sun Avoidance'] = 'none'
    obs['RA'] = obs[ra_col].apply(dec2sexstr, hd='h')
    obs['Dec'] = obs[dec_col].apply(dec2sexstr, hd='d', lead_psign=True)
    obs = obs.rename(columns={'name': 'Source Name', velo_col: 'Velocity'})
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


class HhtTargetList(object):
    vref = 'LSR'
    vsys = 'RAD'
    ctyp = 'J2000.0'
    out_line = '{ra} {dec} {ctyp} {name:<9} {velo:>+7.2f} {vref} {vsys}'
    nearest_off = True
    offs = read_cat('hht_offs')

    def __init__(self, df, coord_cols=['glon_peak', 'glat_peak'],
        velo_col='mol_mol_vlsr', name_col='v210cnum', coord_type='gal'):
        assert coord_type in ['gal', 'eq']
        self.df = df
        self.coord_cols = coord_cols
        self.velo_col = velo_col
        self.name_col = name_col
        self.coord_type = coord_type
        self._convert_coords()
        self._make_text()

    def _convert_coords(self):
        if self.coord_type == 'gal':
            self.df = pd_gal2eq(self.df, self.coord_cols,
                                new_labels=['_ra', '_dec'])
            self.coord_cols = ['_ra', '_dec']
        self.df[self.coord_cols[0] + '_str'] = \
            self.df[self.coord_cols[0]].apply(lambda x: dec2sexstr(x, hd='h'))
        self.df[self.coord_cols[1] + '_str'] = \
            self.df[self.coord_cols[1]].apply(lambda x: dec2sexstr(x, hd='d',
                lead_psign=True))

    def _make_text(self):
        lines = []
        for ii, row in self.df.iterrows():
            lines.append(self.out_line.format(
                         ra=row[self.coord_cols[0] + '_str'],
                         dec=row[self.coord_cols[1] + '_str'],
                         ctyp=self.ctyp,
                         name=row[self.name_col],
                         velo=row[self.velo_col],
                         vref=self.vref,
                         vsys=self.vsys))
            if self.nearest_off:
                coord = (row[self.coord_cols[0]], row[self.coord_cols[1]])
                ii = self._get_nearest_off(coord)
                lines.append(self.out_line.format(
                             ra=self.offs.iloc[ii]['ra_str'],
                             dec=self.offs.iloc[ii]['dec_str'],
                             ctyp=self.ctyp,
                             name=self.offs.iloc[ii]['name'],
                             velo=row[self.velo_col],
                             vref=self.vref,
                             vsys=self.vsys))
        self.text = '\n'.join(lines)

    def _get_nearest_off(self, coord):
        inearest = sep_coords(coord, self.offs[['ra', 'dec']].values).argmin()
        return inearest

    def write(self, out_filen='target_list'):
        out_filen += '.cat'
        self.out_filen = out_filen
        with open(out_filen, 'w') as handle:
            handle.write(self.text)
        print '-- File written to: {0}'.format(out_filen)


