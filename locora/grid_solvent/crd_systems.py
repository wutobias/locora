import numpy as np

from locora.utils.misc import are_you_numpy
from locora.grid_solvent.spatial import field

from locora.grid_solvent._crd_systems_ext import axis_paral
from locora.grid_solvent._crd_systems_ext import axis_ortho

#from locora.grid_solvent.spatial import py_axis_ortho
#from locora.grid_solvent.spatial import py_axis_paral


class internal_rectangular(field):

    """
    Class for creating an internal rectangular coordinate system
    """

    def __init__(self, Dims, Verbose=False):

        ### Note: We should check the correct call of super
        ###	      in python3.X
        self._dims = Dims
        center = np.zeros(3, dtype=np.float_)
        delta  = np.ones(3, dtype=np.float_)
        super(internal_rectangular, self).__init__(Bins=Dims,
                                                   Center=center,
                                                   Delta=delta)

        self.xx     = np.zeros(3, dtype=np.float_)
        self.yy     = np.zeros(3, dtype=np.float_)
        self.zz     = np.zeros(3, dtype=np.float_)

        self.verbose = Verbose

    def set_axis(self, x_set, z_set, x_ref, z_ref):

        """
        x_set: will be transformed to parallel coordinates

        z_set: will be transformed to orthogonal coordinates
        """

        x_use_ext = 1
        z_use_ext = 1

        if type(x_ref) == type(None):
            x_use_ext = 0
            x_ref     = np.zeros(3, dtype=np.float)
        if type(z_ref) == type(None):
            z_use_ext = 0
            z_ref     = np.zeros(3, dtype=np.float)

        axis_paral(x_set, self.xx, x_use_ext, x_ref)
        axis_ortho(z_set, self.zz, z_use_ext, z_ref)
        #py_axis_paral(x_set, self.xx)
        #py_axis_ortho(z_set, self.zz)
        self.yy  = np.cross(self.xx, self.zz)
        self.yy /= np.linalg.norm(self.yy)

        self.frac2real      = np.zeros((3,3), dtype=np.float)
        self.frac2real[:,0] = np.copy(self.xx)
        self.frac2real[:,1] = np.copy(self.yy)
        self.frac2real[:,2] = np.copy(self.zz)

        self.real2frac       = np.linalg.inv(self.frac2real)

        if self.verbose:
            print "frac2real matrix:"
            print self.frac2real[0,:]
            print self.frac2real[1,:]
            print self.frac2real[2,:]
            print "real2frac matrix:"
            print self.real2frac[0,:]
            print self.real2frac[1,:]
            print self.real2frac[2,:]

    def set_center(self, center):

        self.center = center
        self.origin = np.zeros(3)
        self.origin = self.center - self.get_real(self.bins/2)