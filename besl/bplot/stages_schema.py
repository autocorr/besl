"""
=============================
Evolutionary Stages Schematic
=============================

Plot a schematic of the evolutionary stage flag categories.

"""

from __future__ import division

from matplotlib import rc
rc('font', family='serif', size=10)
rc('text', usetex=True)
import daft


x = 6
y = x * 1.6
x_off = -1
y_off = -0.7
pgm = daft.PGM([x, y], origin=[x_off, y_off])

# BGPS Plate
bg_y = 4 / 5 * y
box_dx, box_dy = 2 / 3 * x, y / 14
bg_label = r'{\LARGE ${\rm BGPS} \ \ \lambda = 1.1 \ {\rm mm}$}'
pgm.add_plate(daft.Plate([0, bg_y, box_dx, box_dy], label=bg_label,
    label_offset=[44, 12.5]))
pgm.add_node(daft.Node('bgps_b', '', x / 3, bg_y - 0.1, scale=0))

# Bolocat
bo_y = 3.2 / 5 * y
bo_label = r'{\LARGE $\tt Bolocat \large \ \ {\rm Extracts \ \ Clumps}$}'
pgm.add_plate(daft.Plate([0, bo_y, box_dx, box_dy], label=bo_label,
    label_offset=[20, 12.5]))
pgm.add_node(daft.Node('bolocat_t', '', x / 3, bo_y + box_dy + 0.1, scale=0))
pgm.add_node(daft.Node('bolocat_b', '', x / 3, bo_y - 0.1, scale=0))
pgm.add_edge('bgps_b', 'bolocat_t')

# HCO+
hc_y = 2.4 / 5 * y
hc_label = r'{\LARGE ${\rm HCO}^+ \ J=3-2 \large \ \ {\rm Detection}$}'
pgm.add_plate(daft.Plate([0, hc_y, box_dx, box_dy], label=hc_label,
    label_offset=[19, 13]))
pgm.add_node(daft.Node('hco_t', '', x / 3, hc_y + box_dy + 0.1, scale=0))
pgm.add_edge('bolocat_b', 'hco_t')

# Column Density
cd_y = 2 / 5 * y
hc_label = r'{\LARGE $N({\rm H}_2) > 10^{23} \ [{\rm cm}^{-2}]$}'
pgm.add_plate(daft.Plate([0, cd_y, box_dx, box_dy], label=hc_label,
    label_offset=[41, 12]))
pgm.add_node(daft.Node('cold_b', '', x / 3, cd_y - 0.1, scale=0))

# Dense/Evo Plate
dn_label = r'{\LARGE ${\rm Dense \ Clumps}$}'
pgm.add_plate(daft.Plate([0.77 * x_off, 0.4 * y + y_off, 0.9 * x, 0.575 * y],
    label=dn_label))
fg_label = r'{\LARGE ${\rm Evolutionary \ Flags \ \ 10^{\circ} < \ell < 60^{\circ}}$}'
pgm.add_plate(daft.Plate([0.77 * x_off, 0.75 * y_off, 0.9 * x, 0.375 * y],
    label=fg_label))

# Flag groups
wide = 1.5
f_x, f_y = 0.16 * x, 2.4
pgm.add_node(daft.Node('hg70', r'HG70', 0 * f_x, f_y, scale=wide))
pgm.add_node(daft.Node('ir', r'IR', 1 * f_x, f_y, scale=wide))
pgm.add_node(daft.Node('h2o', r'H2O', 2 * f_x, f_y, scale=wide))
pgm.add_node(daft.Node('ch3oh', r'CH3OH', 3 * f_x, f_y, scale=wide))
pgm.add_node(daft.Node('uchii', r'UCHII', 4 * f_x, f_y, scale=wide))
pgm.add_edge('hg70', 'ir', directed=False,
    plot_params={'linestyle': 'dotted'})
pgm.add_edge('ir', 'h2o', directed=False,
    plot_params={'linestyle': 'dotted'})
pgm.add_edge('h2o', 'ch3oh', directed=False,
    plot_params={'linestyle': 'dotted'})
pgm.add_edge('ch3oh', 'uchii', directed=False,
    plot_params={'linestyle': 'dotted'})

# Auxillary flag groups
pgm.add_node(daft.Node('akari', r'AKARI', 4 * f_x, f_y - 1.05 * f_x, scale=wide,
    observed=True))
pgm.add_node(daft.Node('mipsgal', r'MIPS', 4 * f_x, f_y - 2.1 * f_x, scale=wide,
    observed=True))
pgm.add_edge('akari', 'mipsgal', directed=False,
    plot_params={'linestyle': 'dotted'})

# Starless
sl_dx = 0.25
sl_label = r'\LARGE ${\rm Starless \ Candidates}$'
pgm.add_plate(daft.Plate([0.6 * x_off, 0.05, 0.60 * x, 0.175 * y],
    label=sl_label))
pgm.add_node(daft.Node('irdark', 'IR \n \ Quiesc.', sl_dx + 0 * f_x, 0.107 * y, scale=wide))
pgm.add_node(daft.Node('knownd', 'Known \n \ Dist.', sl_dx + 1 * f_x, 0.107 * y, scale=wide))
pgm.add_node(daft.Node('sl_tk', 'NH$_3$ \n T$_{\\rm K}$', sl_dx + 2 * f_x, 0.107 * y,
    offset=[1, -4.5], scale=wide))
pgm.add_edge('irdark', 'knownd', directed=False,
    plot_params={'linestyle': 'dotted'})
pgm.add_edge('knownd', 'sl_tk', directed=False,
    plot_params={'linestyle': 'dotted'})

# Save Figure
pgm.render()
pgm.figure.savefig('stages_schema.pdf')
pgm.figure.savefig('stages_schema.png', dpi=150)
