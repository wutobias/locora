import numpy as np
from collections import OrderedDict
import mdtraj as md

from locora.grid_solvent import solvent_field
from locora.grid_solvent.solvent_field import solvent_field

from locora.utils.io_operations import write_files

def grid_solvent(option_dict, _start, _stop):

    _results  = OrderedDict()

    traj_path  = option_dict['trajin'][0]
    parm_path  = option_dict['parm'][0]
    zz_sele    = "".join(option_dict['zz'])
    xx_sele    = "".join(option_dict['xx'])
    cntr_sele  = "".join(option_dict['center_sele'])
    dims       = np.array(option_dict['dims'])

    zz_ref_sele = "".join(option_dict['zz_ref'])
    xx_ref_sele = "".join(option_dict['xx_ref'])
    water_str   = "".join(option_dict['water'])

    image       = not option_dict['noimage']

    verbose     = option_dict['verbose']

    if verbose:
        print "Preparing selection masks..."
    no_zz_ref   = False
    zz_ref_crds = None
    if zz_ref_sele.rstrip().lstrip().startswith("None"):
        no_zz_ref = True
    no_xx_ref   = False
    xx_ref_crds = None
    if xx_ref_sele.rstrip().lstrip().startswith("None"):
        no_xx_ref = True

    topo      = md.load_topology(parm_path)

    ### Z axis
    zz_indxs    = topo.select(zz_sele)
    if not no_zz_ref:
        zz_ref_indxs = topo.select(zz_ref_sele)
    ### X axis
    xx_indxs    = topo.select(xx_sele)
    if not no_xx_ref:
        xx_ref_indxs = topo.select(xx_ref_sele)
    ### Center
    if cntr_sele.rstrip().lstrip().startswith("None"):
        center_idxs = np.unique(np.concatenate((xx_indxs, zz_indxs)))
        solvent     = topo.select("%s and not (%s or %s)" %(water_str,zz_sele,xx_sele))
        solvent_O   = topo.select("(%s and name O) and not (%s or %s)" %(water_str,zz_sele,xx_sele))
    else:
        center_idxs = topo.select(cntr_sele)
        solvent     = topo.select("%s and not (%s or %s or %s)" %(water_str,zz_sele,xx_sele,cntr_sele))
        solvent_O   = topo.select("(%s and name O) and not (%s or %s or %s)" %(water_str,zz_sele,xx_sele,cntr_sele))

    sites      = solvent_O[2]-solvent_O[1]
    solv_field = solvent_field(solvent_O.shape[0],
                               sites,
                               dims)

    frame_range = range(_start,_stop)
    N_frames    = _stop-_start

    uc_data      = np.zeros((N_frames*3,3), dtype=np.float)
    pop_data     = np.zeros(N_frames, dtype=np.int)
    center_data  = np.zeros((N_frames,3), dtype=np.float)
    origin_data  = np.zeros((N_frames,3), dtype=np.float)
    O_idxs_data  = None
    theta_data   = None
    phi_data     = None
    psi_data     = None
    xx1_wat_data = None
    xx2_wat_data = None
    yy_wat_data  = None
    zz_wat_data  = None
    O_frac_data  = None
    H1_frac_data = None
    H2_frac_data = None
    frame_data   = None

    if image:
        uc = np.eye(3,3)
    else:
        uc = None

    zz_crds      = np.zeros((zz_indxs.shape[0],3),    dtype=np.float)
    xx_crds      = np.zeros((xx_indxs.shape[0],3),    dtype=np.float)
    solv_crds    = np.zeros((solvent.shape[0],3),     dtype=np.float)
    cntr_crds    = np.zeros((center_idxs.shape[0],3), dtype=np.float)

    if verbose:
        print "Start processing trajectory..."
    with md.open(traj_path) as md_traj:
        started_fill = False
        for i in range(N_frames):
            frame_i = frame_range[i]
            if verbose and i%100==0:
                print "Frame %d..." %frame_i
            if traj_path.endswith(".pdb"):
                frame = md.load_frame(traj_path, index=frame_i, top=topo)
            else:
                md_traj.seek(frame_i)
                frame = md_traj.read_as_traj(topo, n_frames=1, stride=1)

            if not no_zz_ref:
                zz_ref_crds = np.mean(frame.xyz[0][zz_ref_indxs]*10., axis=0)
            if not no_xx_ref:
                xx_ref_crds = np.mean(frame.xyz[0][xx_ref_indxs]*10., axis=0)

            zz_crds   = frame.xyz[0][zz_indxs]*10.
            xx_crds   = frame.xyz[0][xx_indxs]*10.
            solv_crds = frame.xyz[0][solvent]*10.
            cntr_crds = frame.xyz[0][center_idxs]*10.

            if image:
                uc[:] = frame.unitcell_vectors[0]*10.

            solv_field.set_axis(xx_crds, zz_crds, xx_ref_crds, zz_ref_crds)
            solv_field.set_center(cntr_crds.mean(axis=0))
            solv_field.update_field(solv_crds, uc)

            j = i*3
            uc_data[j:j+3,:] = solv_field.get_nice_frac2real()/solv_field.delta
            pop_data[i]      = solv_field.N_inside
            center_data[i]   = solv_field.center
            origin_data[i]   = solv_field.origin

            if solv_field.N_inside > 0:

                if not started_fill:
                    O_idxs_data = np.copy(solv_field.inside_idxs)
                    theta_data  = np.copy(solv_field.theta)
                    phi_data    = np.copy(solv_field.phi)
                    psi_data    = np.copy(solv_field.psi)

                    xx1_wat_data = np.copy(solv_field.xx1_wat)
                    xx2_wat_data = np.copy(solv_field.xx2_wat)
                    yy_wat_data  = np.copy(solv_field.yy_wat)
                    zz_wat_data  = np.copy(solv_field.zz_wat)

                    O_frac_data  = np.copy(solv_field.O_crds_frac)
                    H1_frac_data = np.copy(solv_field.H1_crds_frac)
                    H2_frac_data = np.copy(solv_field.H2_crds_frac)

                    frame_data   = np.empty(solv_field.N_inside, dtype=np.int)
                    frame_data.fill(frame_i)

                    started_fill = True

                else:
                    O_idxs_data      = np.concatenate((O_idxs_data, solv_field.inside_idxs))
                    theta_data       = np.concatenate((theta_data, solv_field.theta))
                    phi_data         = np.concatenate((phi_data, solv_field.phi))
                    psi_data         = np.concatenate((psi_data, solv_field.psi))

                    xx1_wat_data     = np.concatenate((xx1_wat_data, solv_field.xx1_wat))
                    xx2_wat_data     = np.concatenate((xx2_wat_data, solv_field.xx2_wat))
                    yy_wat_data      = np.concatenate((yy_wat_data, solv_field.yy_wat))
                    zz_wat_data      = np.concatenate((zz_wat_data, solv_field.zz_wat))

                    O_frac_data  = np.concatenate((O_frac_data, solv_field.O_crds_frac))
                    H1_frac_data = np.concatenate((H1_frac_data, solv_field.H1_crds_frac))
                    H2_frac_data = np.concatenate((H2_frac_data, solv_field.H2_crds_frac))

                    __frame_data = np.empty(solv_field.N_inside, dtype=np.int)
                    __frame_data.fill(frame_i)
                    frame_data   = np.concatenate((frame_data, __frame_data))

                if verbose and i%100==0:
                    write_files(XYZ=solv_field.O_crds_frac, Format='PDB', Filename='frac_frame%d.pdb' %frame_i)
                    write_files(XYZ=solv_field.O_crds, Format='PDB', Filename='real_frame%d.pdb' %frame_i)
                    

    _results['uc_data']      = uc_data
    _results['pop_data']     = pop_data
    _results['center_data']  = center_data
    _results['origin_data']  = origin_data
    _results['O_idxs_data']  = O_idxs_data
    _results['theta_data']   = theta_data
    _results['phi_data']     = phi_data
    _results['psi_data']     = psi_data
    _results['xx1_wat_data'] = xx1_wat_data
    _results['xx2_wat_data'] = xx2_wat_data
    _results['yy_wat_data']  = yy_wat_data
    _results['zz_wat_data']  = zz_wat_data
    _results['O_frac_data']  = O_frac_data
    _results['H1_frac_data'] = H1_frac_data
    _results['H2_frac_data'] = H2_frac_data
    _results['frame_data']   = frame_data

    return _results

def grid_solvent_out(option_dict, results):

    with open(option_dict['unitcell'][0], 'wb') as f:
        np.savetxt(f, results['uc_data'])

    with open(option_dict['pop'][0], 'wb') as f:
        np.savetxt(f, results['pop_data'])

    with open(option_dict['center'][0], 'wb') as f:
        np.savetxt(f, results['center_data'])

    with open(option_dict['origin'][0], 'wb') as f:
        np.savetxt(f, results['origin_data'])

    with open(option_dict['O_idxs'][0], 'wb') as f:
        np.savetxt(f, results['O_idxs_data'])

    with open(option_dict['theta'][0], 'wb') as f:
        np.savetxt(f, results['theta_data'])

    with open(option_dict['phi'][0], 'wb') as f:
        np.savetxt(f, results['phi_data'])

    with open(option_dict['psi'][0], 'wb') as f:
        np.savetxt(f, results['psi_data'])

    with open(option_dict['xx1_wat'][0], 'wb') as f:
        np.savetxt(f, results['xx1_wat_data'])

    with open(option_dict['xx2_wat'][0], 'wb') as f:
        np.savetxt(f, results['xx2_wat_data'])

    with open(option_dict['yy_wat'][0], 'wb') as f:
        np.savetxt(f, results['yy_wat_data'])

    with open(option_dict['zz_wat'][0], 'wb') as f:
        np.savetxt(f, results['zz_wat_data'])

    with open(option_dict['O_frac'][0], 'wb') as f:
        np.savetxt(f, results['O_frac_data'])

    with open(option_dict['H1_frac'][0], 'wb') as f:
        np.savetxt(f, results['H1_frac_data'])

    with open(option_dict['H2_frac'][0], 'wb') as f:
        np.savetxt(f, results['H2_frac_data'])

    with open(option_dict['frames'][0], 'wb') as f:
        np.savetxt(f, results['frame_data'])
