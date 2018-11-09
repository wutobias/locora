import numpy as np

from locora.grid_solvent.spatial import field
from locora.grid_solvent.spatial import set_euler as _set_euler
from locora.grid_solvent.spatial import set_quaternion as _set_quaternion
from locora.grid_solvent.crd_systems import internal_rectangular

class solvent_field(internal_rectangular):

    """
    Class for calculating and manipulating boxes of solvent molecules.
    """

    def __init__(self, N_solvent, N_sites, Dims, Verbose=False):

        super(solvent_field, self).__init__(Dims)

        self.N_solvent   = N_solvent
        self.N_sites     = N_sites
        self.crds_shape  = (N_solvent*N_sites, 3)
        self.crds        = np.zeros(self.crds_shape, dtype=np.float)
        self.inside_idxs = np.array([], dtype=np.int) ### Array containing oxygen water idxs starting
                                                      ### with indx 0 as the first occuring water molecule
                                                      ### in the selection.

        self.uc          = np.eye(3,3)                    ### uc: Simulation box vectors as matrix presentation
        self.mic_idx     = np.zeros(N_solvent, dtype=int) ### mic_idx: current simulation box idx for each water 
                                                          ###          according minimum-image-convetion

        self.N_inside    = 0

        ### These coordinates and vectors are being calculated in frac space
        self.xx1_wat     = 0              ### O-H1 vector
        self.xx2_wat     = 0              ### O-H2 vector
        self.yy_wat      = 0              ### cross product zz_wat and xx1_wat
        self.zz_wat      = 0              ### O-H1/O-H2 orthogonal vector

        self.O_crds_frac  = 0
        self.H1_crds_frac = 0
        self.H2_crds_frac = 0

        ### These coordinates are being calculated in real space
        self.O_crds      = 0
        self.H1_crds     = 0
        self.H2_crds     = 0

        self.q_dist      = 0
        self.theta       = 0
        self.phi         = 0
        self.psi         = 0

        self.a = [-1.,0.,1.]
        self.b = [-1.,0.,1.]
        self.c = [-1.,0.,1.]

        self.verbose     = Verbose

        if self.verbose:
            print "N_solvent:", self.N_solvent
            print "N_sites:  ", self.N_sites

    
    def update_field(self, crds, uc=None):

        if crds.shape != self.crds_shape:
            raise Warning("New crds array has shape %s. Original crds array had shape %s. Both must have equal shape." %(crds.shape, self.crds_shape))

        self.crds[:] = crds

        image=False
        if type(uc)!=None:
            self.uc[:] = uc
            image = True

        self._set_inside_crds(image)
        self._set_watcrds()
        self._set_euler()
        self._set_quaternion()


    def _set_inside_crds(self, image):

        if image:

            ###FIXME: Port this routine to a C extension.

            uc_inv    = np.linalg.inv(self.uc)
            crds      = self.crds[::self.N_sites]
            _crds     = np.zeros_like(crds)
            crds_inv  = crds.dot(uc_inv)
            crds_inv  = crds_inv - np.floor(crds_inv)
            _crds_inv = np.zeros_like(crds_inv)

            solvent_frac   = np.zeros_like(crds_inv)
            _solvent_frac  = np.zeros_like(crds_inv)
            nearest_frac2  = np.zeros(self.N_solvent, dtype=float)
            _nearest_frac2 = np.zeros(self.N_solvent, dtype=float)

            center_inv  = self.center.dot(uc_inv)
            center_inv  = center_inv - np.floor(center_inv)
            self.center = center_inv.dot(self.uc)

            self.origin = np.zeros(3)
            self.origin = self.center - self.get_real(self.bins/2)

            i = 0
            for a_i in self.a:
                for b_i in self.b:
                    for c_i in self.c:
                        _crds_inv[:]     = crds_inv + np.array([a_i, b_i, c_i])
                        _crds[:]         = _crds_inv.dot(self.uc)
                        _solvent_frac[:] = self.get_frac(_crds)

                        if i==0:
                            nearest_frac2[:]  = np.power(_solvent_frac-self.bins*0.5, 2).sum(axis=1)
                            solvent_frac[:]   = np.copy(_solvent_frac)
                            self.mic_idx[:]   = 0

                        else:
                            _nearest_frac2[:]     = np.power(_solvent_frac-self.bins*0.5, 2).sum(axis=1)
                            update                = np.where(_nearest_frac2<nearest_frac2)[0]
                            nearest_frac2[update] = np.copy(_nearest_frac2[update])
                            solvent_frac[update]  = _solvent_frac[update]
                            self.mic_idx[update]  = i

                        i += 1

            valids = np.where( (solvent_frac[:,0] >= 0.) * (solvent_frac[:,0] < self.bins[0]) * \
                               (solvent_frac[:,1] >= 0.) * (solvent_frac[:,1] < self.bins[1]) * \
                               (solvent_frac[:,2] >= 0.) * (solvent_frac[:,2] < self.bins[2]) )[0]

            self.N_inside    = valids.shape[0]
            self.inside_idxs = valids*self.N_sites

            if self.verbose:
                print "Found %d water molecules inside grid." %self.N_inside

            O_crds_inv  = self.crds[self.inside_idxs].dot(uc_inv)
            H1_crds_inv = self.crds[self.inside_idxs+1].dot(uc_inv)
            H2_crds_inv = self.crds[self.inside_idxs+2].dot(uc_inv)

            O_crds_inv  = O_crds_inv  - np.floor(O_crds_inv)
            H1_crds_inv = H1_crds_inv - np.floor(H1_crds_inv)
            H2_crds_inv = H2_crds_inv - np.floor(H2_crds_inv)

            _mic_idx    = self.mic_idx[valids]

            self.O_crds_frac  = np.zeros_like(O_crds_inv)
            self.H1_crds_frac = np.zeros_like(H1_crds_inv)
            self.H2_crds_frac = np.zeros_like(H2_crds_inv)

            ### These are the reals space coordinates
            self.O_crds      = np.zeros_like(O_crds_inv)
            self.H1_crds     = np.zeros_like(H1_crds_inv)
            self.H2_crds     = np.zeros_like(H2_crds_inv)

            i=0
            for a_i in self.a:
                for b_i in self.b:
                    for c_i in self.c:
                        sele = np.where(_mic_idx==i)[0]
                        cell = np.array([a_i, b_i, c_i])
                        if sele.shape[0]>0:
                            crds_mic_inv           = O_crds_inv[sele] + cell
                            crds_mic               = crds_mic_inv.dot(self.uc)
                            solvent_frac           = self.get_frac(crds_mic)
                            self.O_crds_frac[sele] = np.copy(solvent_frac)
                            self.O_crds[sele]      = np.copy(crds_mic)

                            crds_mic_inv            = H1_crds_inv[sele] + cell
                            crds_mic                = crds_mic_inv.dot(self.uc)
                            solvent_frac            = self.get_frac(crds_mic)
                            self.H1_crds_frac[sele] = np.copy(solvent_frac)
                            self.H1_crds[sele]      = np.copy(crds_mic)

                            crds_mic_inv            = H2_crds_inv[sele] + cell
                            crds_mic                = crds_mic_inv.dot(self.uc)
                            solvent_frac            = self.get_frac(crds_mic)
                            self.H2_crds_frac[sele] = np.copy(solvent_frac)
                            self.H2_crds[sele]      = np.copy(crds_mic)

                        i += 1

        else:

            crds         = self.crds[::self.N_sites]
            solvent_frac = self.get_frac(crds)

            valids = np.where( (solvent_frac[:,0] >= 0.) * (solvent_frac[:,0] < self.bins[0]) * \
                               (solvent_frac[:,1] >= 0.) * (solvent_frac[:,1] < self.bins[1]) * \
                               (solvent_frac[:,2] >= 0.) * (solvent_frac[:,2] < self.bins[2]) )[0]

            self.N_inside    = valids.shape[0]
            self.inside_idxs = valids*self.N_sites

            if self.verbose:
                print "Found %d water molecules inside grid." %self.N_inside

            ### These coordinates are being calculated in real space
            self.O_crds      = self.crds[self.inside_idxs]
            self.H1_crds     = self.crds[self.inside_idxs+1]
            self.H2_crds     = self.crds[self.inside_idxs+2]

            self.O_crds_frac  = self.get_frac(self.O_crds)
            self.H1_crds_frac = self.get_frac(self.H1_crds)
            self.H2_crds_frac = self.get_frac(self.H2_crds)


    def _set_watcrds(self):
    
        self.xx1_wat = self.H1_crds_frac - self.O_crds_frac
        xx1_norm     = np.linalg.norm(self.xx1_wat, axis=-1)
        self.xx1_wat = np.einsum('ij,i->ij', self.xx1_wat, 1./xx1_norm)

        self.xx2_wat = self.H2_crds_frac - self.O_crds_frac
        xx2_norm     = np.linalg.norm(self.xx2_wat, axis=-1)
        self.xx2_wat = np.einsum('ij,i->ij', self.xx2_wat, 1./xx2_norm)
    
        self.zz_wat  = np.cross(self.xx1_wat, self.xx2_wat)
        zz_norm      = np.linalg.norm(self.zz_wat, axis=-1)
        self.zz_wat  = np.einsum('ij,i->ij', self.zz_wat, 1./zz_norm)
    
        self.yy_wat  = np.cross(self.xx1_wat, self.zz_wat)
        yy_norm      = np.linalg.norm(self.yy_wat, axis=-1)
        self.yy_wat  = np.einsum('ij,i->ij', self.yy_wat, 1./yy_norm)
    

    def _set_euler(self):

        self.theta = np.zeros(self.N_inside, dtype=np.float)
        self.phi   = np.zeros(self.N_inside, dtype=np.float)
        self.psi   = np.zeros(self.N_inside, dtype=np.float)

        _set_euler(self.O_crds,
                   self.H1_crds,
                   self.H2_crds,
                   self.xx,
                   self.yy,
                   self.zz,
                   self.theta,
                   self.phi,
                   self.psi)


    def _set_quaternion(self):
    
        self.q_dist = np.zeros((self.N_inside, 4), dtype=np.float)

        _set_quaternion(self.q_dist,
                        self.theta, 
                        self.phi, 
                        self.psi)

