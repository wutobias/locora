import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from collections import OrderedDict

from locora.process_data.utils import fit_lifetime
from locora.process_data.utils import int_lifetime
from locora.process_data._utils_ext import lifetime_distr_ext
from locora.process_data._utils_ext import lifetime_distr_majmin_ext
from locora.process_data.utils import _make_1dhist, _make_2dhist

from locora.utils.misc import are_you_numpy
from locora.utils._read_write_ext import read_table_ext

from locora.process_data.stationary_block_bootstrap import resample

class solvent_temporal_analysis(object):

    def __init__(self, option_dict, start, stop):

        self.cutoff      = float(option_dict['cutoff'])
        self.dims        = np.array(option_dict['dims'], dtype=float)
        self.N_resample  = int(option_dict['bootstrap'])
        self.window_size = int(option_dict['window'])
        self.timestep    = float(option_dict['timestep'])
        self.nonexp      = option_dict['nonexp']
        self.plot        = option_dict['plot']
        self.transient   = 0
        self.aniso       = option_dict['anisotropic']
        self.minstep     = int(option_dict['minstep'])

        self.start       = int(start)
        self.stop        = int(stop)
        self.traj_length = self.stop-self.start

        self.N_windows   = int(self.traj_length)/int(self.window_size)

        self.do_bootstrap = False

        if self.N_resample>0:
            self.N_windows    = self.N_resample
            self.do_bootstrap = True
            np.random.seed(np.random.randint(99999999))

        self.prefix  = str(option_dict['prefix'])

        self.verbose = option_dict['verbose']

        self.check()

        #__table_load = np.loadtxt
        __table_load = read_table_ext

        if self.verbose:
            print "Loading %s ..." %option_dict['O_idxs'][0]
        self.O_idxs_data  = __table_load(option_dict['O_idxs'][0]).astype(np.int)
        if self.verbose:
            print "Loading %s ..." %option_dict['xx1_wat'][0]
        self.xx1_wat_data = __table_load(option_dict['xx1_wat'][0]).astype(np.float)
        if self.verbose:
            print "Loading %s ..." %option_dict['xx2_wat'][0]
        self.xx2_wat_data = __table_load(option_dict['xx2_wat'][0]).astype(np.float)
        if self.verbose:
            print "Loading %s ..." %option_dict['yy_wat'][0]
        self.yy_wat_data  = __table_load(option_dict['yy_wat'][0]).astype(np.float)
        if self.verbose:
            print "Loading %s ..." %option_dict['zz_wat'][0]
        self.zz_wat_data  = __table_load(option_dict['zz_wat'][0]).astype(np.float)
        if self.verbose:
            print "Loading %s ..." %option_dict['O_frac'][0]
        self.O_frac_data  = __table_load(option_dict['O_frac'][0]).astype(np.float)
        if self.verbose:
            print "Loading %s ..." %option_dict['frames'][0]
        self.frame_data   = __table_load(option_dict['frames'][0]).astype(np.int)

        self.transient = 0
        self.plane     = [[], []]
        self.ortho     = -1
        self.sign      = "+-"

        self.selection = []

        self.P_order   = 1

        self.tau_trans = np.zeros(self.N_windows, dtype=float)
        self.a_trans   = np.zeros(self.N_windows, dtype=float)
        self.b_trans   = np.zeros(self.N_windows, dtype=float)
        self.c_trans   = np.zeros(self.N_windows, dtype=float)
        if self.nonexp:
            self.d_trans = np.zeros(self.N_windows, dtype=float)

        self.tau_orient = np.zeros((3, self.N_windows), dtype=float)
        self.a_orient   = np.zeros((3, self.N_windows), dtype=float)
        self.b_orient   = np.zeros((3, self.N_windows), dtype=float)
        self.c_orient   = np.zeros((3, self.N_windows), dtype=float)
        if self.nonexp:
            self.d_orient = np.zeros((3, self.N_windows), dtype=float)

        self.header_trans  = "# t* "
        self.header_orient = "# t* "

        self.scope_trans_avg  = ""
        self.scope_trans_std  = ""
        self.scope_orient_avg = ""
        self.scope_orient_std = ""

        ### Check out water molecules that are within cutoff of center
        center_dist = np.linalg.norm(self.O_frac_data-self.dims*0.5, axis=1)
        self.center = np.where(center_dist<self.cutoff)[0]
        del center_dist


    def check(self):

        if self.verbose:
            print "Sanity checking..."
        if self.cutoff>self.dims[0]:
            raise ValueError("Cutoff is set to %s, but must be <%s." (self.cutoff[0], self.dims[0]))
        if self.cutoff>self.dims[1]:
            raise ValueError("Cutoff is set to %s, but must be <%s." (self.cutoff[1], self.dims[1]))
        if self.cutoff>self.dims[2]:
            raise ValueError("Cutoff is set to %s, but must be <%s." (self.cutoff[2], self.dims[2]))
        if self.N_resample<0:
            raise ValueError("bootstrap must be set to >=0. bootstrap=0 means don't perform bootstrap.")
        if self.window_size<0:
            raise ValueError("window must be set to >=0. window_size=0 means perform calculations \
                              over whole trajectory and not over seperate windows.")
        if self.traj_length<self.window_size:
            raise ValueError("window cannot be longer than length of trajectory %d." %self.traj_length)


    def set_order(self, P_order):

        self.P_order = P_order

    def set_transient(self, transient):

        if transient<0:
            raise ValueError("Transient must be >=0.")
        self.transient = int(transient)


    def set_plane(self, planes):

        self.ortho = -1
        for i in range(3):
            if i not in planes:
                self.ortho = i
                break

        self.plane = [[], []]
        for plane in planes:
            if plane>2:
                raise ValueError("No element in axis must be greater than 2.")
            self.plane[0].append(plane)
            if plane==0:
                self.plane[1].append("X")
            if plane==1:
                self.plane[1].append("Y")
            if plane==2:
                self.plane[1].append("Z")

        if self.ortho == -1:
            raise ValueError("Could not assign correct plane.")


    def set_sign(self, sign):

        if sign not in ["+", "-"]:
            raise ValueError("Sign must be one of [\"+\", \"-\"]")
        self.sign = sign


    def select(self):

        if not self.aniso:
            self.selection = self.center

        elif self.sign=="+":
            self.selection = np.where(self.O_frac_data[self.center,self.ortho]>=self.dims[self.ortho]*0.5)[0]
            self.selection = self.center[self.selection]

        elif self.sign=="-":
            self.selection = np.where(self.O_frac_data[self.center,self.ortho]<self.dims[self.ortho]*0.5)[0]
            self.selection = self.center[self.selection]

        else:
            raise Exception("Could not make anisotropic selection.")

        if self.verbose:
            print "Found %d solvent molecules in" %self.selection.shape[0],
            if self.aniso:
                print "anisotropic",
            else:
                print "isotropic",
            print "selection ..."


    def plot_timeaverage(self):

        if self.aniso:
            ax=self.plane[0]
            ax_name=self.plane[1]
            title="%s%s %s" %(ax_name[0],ax_name[1],self.sign)
            filename="histogramm2d_%s%s_%s.png" %(ax_name[0],ax_name[1],self.sign)
            try:
                _make_2dhist(self.O_frac_data[self.selection][:,ax],
                            [self.dims[ax[0]], self.dims[ax[1]]],
                            ax_name,
                            title,
                            filename,
                            self.prefix)
            except:
                print "Could not make 2d histogram for %s_%s" %(self.prefix, filename)

        else:
            axis=[[0,1],[0,2],[1,2]]
            axis_names=[["X","Y"],["X","Z"],["Y","Z"]]
            for ax, ax_name in zip(axis, axis_names):
                title="%s%s" %(ax_name[0], ax_name[1])
                filename="histogramm2d_%s%s.png" %(ax_name[0],ax_name[1])
                try:
                    _make_2dhist(self.O_frac_data[self.selection][:,ax],
                                [self.dims[ax[0]], self.dims[ax[1]]],
                                ax_name,
                                title,
                                filename,
                                self.prefix)
                except:
                    print "Could not make 2d histogram for %s_%s" %(self.prefix, filename)


    def out_temporal_model(self):

        self.header_trans  = ""
        self.header_orient = ""

        self.scope_trans_avg  = ""
        self.scope_trans_std  = ""
        self.scope_orient_avg = ""
        self.scope_orient_std = ""

        plot_trans=OrderedDict()
        plot_orient=OrderedDict()

        valids            = np.invert(np.isnan(self.tau_trans))

        plot_trans['tau'] = self.tau_trans[valids]
        plot_trans['a']   = self.a_trans[valids]
        plot_trans['b']   = self.b_trans[valids]
        plot_trans['c']   = self.c_trans[valids]
        if self.nonexp:
            plot_trans['d'] = self.d_trans[valids]

        plot_orient['tau-x1'] = self.tau_orient[0,valids]
        plot_orient['a-x1']   = self.a_orient[0,valids]
        plot_orient['b-x1']   = self.b_orient[0,valids]
        plot_orient['c-x1']   = self.c_orient[0,valids]
        if self.nonexp:
            plot_orient['d-x1'] = self.d_orient[0,valids]

        plot_orient['tau-y'] = self.tau_orient[1,valids]
        plot_orient['a-y']   = self.a_orient[1,valids]
        plot_orient['b-y']   = self.b_orient[1,valids]
        plot_orient['c-y']   = self.c_orient[1,valids]
        if self.nonexp:
            plot_orient['d-y'] = self.d_orient[1,valids]

        plot_orient['tau-z'] = self.tau_orient[2,valids]
        plot_orient['a-z']   = self.a_orient[2,valids]
        plot_orient['b-z']   = self.b_orient[2,valids]
        plot_orient['c-z']   = self.c_orient[2,valids]
        if self.nonexp:
            plot_orient['d-z'] = self.d_orient[2,valids]

        for name, data in plot_trans.items():
            
            title="trans "+name+" t*=%d" %self.transient
            filename="trans_"+name+"_t=%d" %self.transient

            self.header_trans += name

            if self.aniso:
                ax=self.plane[0]
                ax_name=self.plane[1]
                title+=" %s%s %s " %(ax_name[0],ax_name[1],self.sign)
                filename+="_%s%s_%s.png" %(ax_name[0],ax_name[1],self.sign)
                self.header_trans+="_%s%s_%s" %(ax_name[0],ax_name[1],self.sign)
            else:
                filename+=".png"

            if self.plot:
                try:
                    _make_1dhist(data, title, filename, self.prefix)
                except:
                    print "Could not make 1d histogram for %s_%s" %(self.prefix, filename)

            self.scope_trans_avg  += "%6.3f " %data.mean()
            self.scope_trans_std  += "%6.3f " %data.std()

            self.header_trans += " "

        for name, data in plot_orient.items():
            
            title="orient "+name+" t*=%d" %self.transient
            filename="orient_"+name+"_t=%d" %self.transient

            self.header_orient += name

            if self.aniso:
                ax=self.plane[0]
                ax_name=self.plane[1]
                title+=" %s%s %s" %(ax_name[0],ax_name[1],self.sign)
                filename+="_%s%s_%s.png" %(ax_name[0],ax_name[1],self.sign)
                self.header_orient+="_%s%s_%s" %(ax_name[0],ax_name[1],self.sign)
            else:
                filename+=".png"

            if self.plot:
                try:
                    _make_1dhist(data, title, filename, self.prefix)
                except:
                    print "Could not make 1d histogram for %s_%s" %(self.prefix, filename)

            self.scope_orient_avg  += "%6.3f " %data.mean()
            self.scope_orient_std  += "%6.3f " %data.std()

            self.header_orient     += " "


    def _process(self):

        if self.verbose:
            ax_name=self.plane[1]
            print "Processing with",
            print "t*=%d" %self.transient,
            if self.aniso:
                print "plane=%s%s" %(ax_name[0],ax_name[1]),
                print "sign=%s" %self.sign,
            print "..."

        if self.do_bootstrap:
            _startlist = np.random.randint(self.start,
                                            self.stop-self.window_size, 
                                            self.N_windows)

        _start = 0
        _stop  = 0
        for i in range(self.N_windows):
            if self.do_bootstrap:
                _start = _startlist[i]
            else:
                _start = _stop
            _stop  = _start+self.window_size
            if self.verbose:
                print "Processing frame %d to %d ..." %(_start, _stop)
            window_idxs   = np.where((self.frame_data[self.selection]>=_start) *\
                                     (self.frame_data[self.selection]<_stop))[0]
            if self.verbose:
                print "Found %d water molecules in time interval ..." %window_idxs.shape[0]
            if window_idxs.shape[0]<self.minstep:
                if self.verbose:
                    print "Skipped lifetime calculation for this window due to low occupancy."
                self.tau_trans[i] = np.nan
                self.a_trans[i]   = np.nan
                self.b_trans[i]   = np.nan
                self.c_trans[i]   = np.nan
                if self.nonexp:
                    self.d_trans[i] = np.nan

                self.tau_orient[:, i] = np.nan
                self.a_orient[:, i]   = np.nan
                self.b_orient[:, i]   = np.nan
                self.c_orient[:, i]   = np.nan
                if self.nonexp:
                    self.d_orient[:, i] = np.nan
                continue

            _frame_data   = self.frame_data[self.selection][window_idxs]
            _O_idxs_data  = self.O_idxs_data[self.selection][window_idxs]
            _xx1_wat_data = self.xx1_wat_data[self.selection][window_idxs]
            _xx2_wat_data = self.xx2_wat_data[self.selection][window_idxs]
            _yy_wat_data  = self.yy_wat_data[self.selection][window_idxs]
            _zz_wat_data  = self.zz_wat_data[self.selection][window_idxs]

            ### |--------------------| ###
            ### | Translational part | ###
            ### |--------------------| ###
            corr_distr  = lifetime_distr_ext(_O_idxs_data,
                                             _frame_data,
                                             self.transient,
                                             0).astype(float)
            corr_distr /= corr_distr[0]

            self.tau_trans[i] = int_lifetime(corr_distr, self.timestep)

            popt, pcov  = fit_lifetime(corr_distr, 
                                       dt=self.timestep, 
                                       non_exponential=self.nonexp)

            if popt[0] == np.nan or popt.shape[0]==1:
                self.a_trans[i] = np.nan
                self.b_trans[i] = np.nan
                self.c_trans[i] = np.nan
                if self.nonexp:
                    self.d_trans[i] = np.nan

            else:
                self.a_trans[i] = popt[0]
                self.b_trans[i] = popt[1]
                self.c_trans[i] = popt[2]
                if self.nonexp:
                    self.d_trans[i] = popt[3]


            ### |--------------------| ###
            ### | Orientational part | ###
            ### |--------------------| ###
            for _wat_data_i, _wat_data in enumerate([_xx1_wat_data, _yy_wat_data, _zz_wat_data]):

                corr_distr  = lifetime_distr_majmin_ext(_O_idxs_data,
                                                        _frame_data,
                                                        _wat_data,
                                                        self.transient, 
                                                        self.P_order,
                                                        0).astype(float)
                corr_distr /= corr_distr[0]

                self.tau_orient[_wat_data_i, i] = int_lifetime(corr_distr, self.timestep)

                popt, pcov  = fit_lifetime(corr_distr,
                                           dt=self.timestep,
                                           non_exponential=self.nonexp)

                if popt[0] == np.nan or popt.shape[0]==1:
                    self.a_orient[_wat_data_i, i] = np.nan
                    self.b_orient[_wat_data_i, i] = np.nan
                    self.c_orient[_wat_data_i, i] = np.nan
                    if self.nonexp:
                        self.d_orient[_wat_data_i, i] = np.nan

                else:
                    self.a_orient[_wat_data_i, i] = popt[0]
                    self.b_orient[_wat_data_i, i] = popt[1]
                    self.c_orient[_wat_data_i, i] = popt[2]
                    if self.nonexp:
                        self.d_orient[_wat_data_i, i] = popt[3]

        self.out_temporal_model()