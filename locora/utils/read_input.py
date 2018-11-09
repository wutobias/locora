import os
from collections import OrderedDict

def read_infile(file_path):

    option_dict    = OrderedDict()
    optional_opts  = OrderedDict()

    optional_opts['start']       = ['0']
    optional_opts['stop']        = ['-1']
    optional_opts['start']       = ['0']
    optional_opts['unitcell']    = ['unitcell.dat']
    optional_opts['pop']         = ['pop.dat']
    optional_opts['center']      = ['center.dat']
    optional_opts['origin']      = ['origin.dat']
    optional_opts['O_idxs']      = ['O_idxs.dat']
    optional_opts['theta']       = ['theta.dat']
    optional_opts['phi']         = ['phi.dat']
    optional_opts['psi']         = ['psi.dat']
    optional_opts['xx1_wat']     = ['xx1_wat.dat']
    optional_opts['xx2_wat']     = ['xx2_wat.dat']
    optional_opts['yy_wat']      = ['yy_wat.dat']
    optional_opts['zz_wat']      = ['zz_wat.dat']
    optional_opts['O_frac']      = ['O_frac.dat']
    optional_opts['H1_frac']     = ['H1_frac.dat']
    optional_opts['H2_frac']     = ['H2_frac.dat']
    optional_opts['frames']      = 'frames.dat'
    optional_opts['dims']        = ['10', '10', '10']
    optional_opts['xx_ref']      = ['None']
    optional_opts['zz_ref']      = ['None']
    optional_opts['water']       = ['water']
    optional_opts['center_sele'] = ['None']

    required_opts = ['trajin',
                     'parm',
                     'xx',
                     'zz']

    if not os.path.exists(file_path):
        raise IOError("File %s not found." %file_path)

    with open(file_path, 'r') as file:
        for line in file:
            l = line.rstrip().lstrip().split()
            
            if len(l) == 0:
                continue
            if l[0].startswith('#'):
                continue

            if len(l) == 1:
                raise IOError("Option %s not understood." %l[0])
            else:
                option_dict[l[0]] = list()
                for val in l[1:]:
                    option_dict[l[0]].append(val)

    for key in required_opts:
        if key not in option_dict.keys():
            raise IOError("Keyword %s not found." %key)

    for key, value in optional_opts.items():
        if key not in option_dict.keys():
            option_dict[key] = value

    option_dict['dims'][0] = float(option_dict['dims'][0])
    option_dict['dims'][1] = float(option_dict['dims'][1])
    option_dict['dims'][2] = float(option_dict['dims'][2])

    option_dict['start'][0] = int(option_dict['start'][0])
    option_dict['stop'][0]  = int(option_dict['stop'][0])

    sele_opts = ['xx',
                 'zz', 
                 'center_sele', 
                 'xx_ref', 
                 'zz_ref', 
                 'water']

    for sele_opt in sele_opts:
        for i in range(len(option_dict[sele_opt])-1):
            option_dict[sele_opt][i] += " "

    return option_dict

example_input = """
### Parameter and trajectory data
trajin my_traj.nc
parm   my_parameter.prmtop
start  0
stop   -1

### Unitcell definitions
### (selection masks)
xx resid 205 and (name CZ  or name CG)
zz resid 205 and (name CD1 or name CD2 or name CE1 or name CE2)
center_sele resid 205 and (name CD1 or name CD2 or name CE1 or name CE2)
dims 10 10 10
xx_ref None
zz_ref resid 192 and (name CG or name CD1 or name NE1 or name CE2 or name CZ2 or name CH2 or name CZ3 or name CE3 or name CD2)
water  water

### Output data
unitcell   unitcell.dat
pop        pop.dat
center     center.dat
origin     origin.dat
O_idxs     O_idxs.dat
theta      theta.dat
phi        phi.dat
psi        psi.dat
xx1_wat    xx1_wat.dat
xx2_wat    xx2_wat.dat
yy_wat     yy_wat.dat
zz_wat     zz_wat.dat
frames     frames.dat
O_frac     O_frac.dat
H1_frac    H1_frac.dat
H2_frac    H2_frac.dat
"""
