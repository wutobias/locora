import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

def generate_timeseries(N, k):

    """
    Generate binary time series with exponential decay.
    Can be used for a hypothetical reaction like A -> B, where N(t)=1
    if in state A and N(t)=0 if in state B.
    
    k: decay constant (1/T)
    N: Total Number of attempts
    
    Idea:
    -----
    At every time t, the number of particles leaving state A within t and d+dt equals
    -k*N(t):

    dN(t)/dt = -k*N(t)
    
    The decay probability depends only on the time interval [t,t+dt]. The probability 
    for a particle to undergo the transformation A --> B therefore is equal to
    
    P(dt) = k*dt

    """
    s = list()
    k = float(k)
    for n in range(N):
        _new = 1
        while _new>0:
            r = np.random.random()
            if r < k:
                _new = 0
                #print r, k
            else:
                _new = 1
            s.append(_new)
    return np.array(s)


def generate_timeseries_2(N,k):
    
    """
    Like generate_timeseries, but returns array of fixed length N. Also, the random variable can switch back
    from 0 to 1 (0->1) with probability 1-k.
    """

    s     = np.zeros(N, dtype=np.int)
    state = 1
    s[0]  = 1
    for i in range(1,N):
        r = np.random.random()
        if r<k:
            state = 0
        else:
            state = 1
        s[i] = state
    return s


def lifetime_distr(series):
    
    """
    Method 1:
    Identifies lifetimes for continous series of data.
    A continous series of data ends, as soon as it drops to zero.
    """

    _sum         = 1.
    _lifetime    = np.zeros(series.shape)
    _i           = 0
    series_work  = np.concatenate(([0], series))
    while _sum > 0.:
        series_d      = np.diff(series_work)
        series_start  = np.where(series_d>0)[0]
        _sum          = series_start.shape[0]
        _lifetime[_i] = _sum
        _i           += 1

        series_work[series_start+1] = 0        
    return _lifetime


def lifetime_distr_2(series):
    
    """
    Method 2:
    Identifies lifetimes for continous series of data.
    A continous series of data ends, as soon as it drops to zero.
    """
    
    series_d     = np.concatenate(([0], series))
    series_d     = np.diff(series_d)
    series_start = np.where(series_d>0.)[0]
    series_stop  = np.where(series_d<0.)[0]
    if series_stop.shape[0] < series_start.shape[0]:
        series_stop = np.concatenate((series_stop, [series.shape[0]]))
    lt_values = series_stop - series_start
    lt_series = np.zeros(series.shape)
    for lt_i in range(series.shape[0]):
        lt_series[lt_i] = np.where(lt_values>lt_i)[0].shape[0]
        if lt_series[lt_i] == 0.:
            break    
    return lt_series


def int_lifetime(X, dt=1.):
    
    """
    Integrate lifetime distribution from 0 to the upper boundary.
    """
    
    #L      = X.shape[0]
    #t      = np.arange(0,L*dt, dt)
    
    #return np.trapz(X, t)
    return np.trapz(X, dx=dt)


def func_3(x, a, b, c):
    
    """
    Analytical form of a first order lifetime distribution.
    """
    
    return a * np.exp(-(x/b)) + c


def func_4(x, a, b, c, d):
    
    """
    Analytical form of a (pseudo-)first order lifetime distribution.
    The parameter c corrects for non-exponential behaviour.
    """
    
    return a * np.exp(-(x/b)**d) + c


def fit_lifetime(X, dt=2., non_exponential=False):
    
    """
    Method for fitting data to a first order 
    lifetime distribution function.
    """

    L      = X.shape[0]
    t      = np.arange(0, L, dtype=np.float)*dt

    valids = ~np.isnan(X)

    try:

        if non_exponential:
            popt, pcov = curve_fit(func_4, t[valids], X[valids],
                                    bounds=(0.,np.inf), 
                                    method='trf')
        else:
            popt, pcov = curve_fit(func_3, t[valids], X[valids],
                                    bounds=(0.,np.inf), 
                                    method='trf')

    except:

        popt = np.array([np.nan])
        pcov = np.array([np.nan])

    return popt, pcov


def get_stats(X, dt=2.):
    
    """
    Retrieve minimum, maximum and average value of a
    distrubtion of autocorrelation coefficients.
    """
    
    L      = X.shape[0]
    t      = np.arange(0, L, dtype=np.float)*dt

    valids = np.where(X > 0.)
    _min   = np.min(t[valids])
    _max   = np.max(t[valids])
    _avg   = np.sum(X*t)/np.sum(X)
    
    return _min, _max, _avg


class hist2d_it(object):
    
    """
    Generate 2d histograms.
    """
    
    def __init__(self, data, bins, xlim=None, ylim=None):
    
        self.data = data
        self.bins = bins
        
        self.hist   = None
        self.xedges = None
        self.yedges = None
        
        self.xlim   = xlim
        self.ylim   = ylim
        
        self.dx     = None
        self.dy     = None
        
        if self.xlim==None:

            self.xlim = np.min(self.data[:,0]), np.max(self.data[:,0])
            
        if self.ylim==None:
            
            self.ylim = np.min(self.data[:,1]), np.max(self.data[:,1])
            
    def hist_it(self):
    
        self.hist, self.xedges, self.yedges = np.histogram2d(self.data[:,0],
                                                             self.data[:,1],
                                                             bins=self.bins,
                                                             range=[self.xlim, self.ylim],
                                                             normed=True)

        self.hist = self.hist.T
        
        self.dx   = self.xedges[1] - self.xedges[0]
        self.dy   = self.yedges[1] - self.yedges[0]


class hist1d_it(object):

    """
    Generate 1d histograms
    """

    def __init__(self, x, xlim=None, num_bins=200):

        if xlim == None:

            self.xmin, self.xmax = x.min(), x.max()

        else:

            self.xmin, self.xmax = xlim[0], xlim[1]

        self.x           = x
        self.num_bins    = num_bins
        self.kernel      = stats.gaussian_kde(self.x, bw_method="scott")
        self.xx          = np.linspace(self.xmin, self.xmax, num_bins)
        self.f           = self.kernel.evaluate(self.xx)


    def plot(self, xlab="X0", ylab="Normalized Density", ylim=None, title="Normalized Density Plot", name="output.png", color=None, kde=True, norm=True, weights=None, vertical=None, return_plt=False, label=None):

        fig       = plt.figure()
        ax        = fig.gca()

        if not kde:

            if label != None:

                n, bins, patches = plt.hist(self.x, self.xx, color=color, weights=weights, density=norm, facecolor='green', alpha=0.5, label=label)

            else:

                n, bins, patches = plt.hist(self.x, self.xx, color=color, weights=weights, density=norm, facecolor='green', alpha=0.5)

            ax.set_xlim(self.xmin, self.xmax)

        else:

            if label != None:

                ax.plot(self.xx, self.f, color=color, label=label)

            else:

                ax.plot(self.xx, self.f, color=color)

            ax.set_xlim(self.xmin, self.xmax)
            ax.set_ylim(0.0, np.max(self.f)+0.0001)
      
        ax.set_xlabel(r"%s"%xlab)
        ax.set_ylabel(r"%s"%ylab)

        if ylim != None:

            ax.set_ylim(ylim[0], ylim[1])

        if vertical != None:

            ax.axvline(x = vertical, linewidth=2, linestyle="dashed", color="r")

        # Tweak spacing to prevent clipping of ylabel
        fig.subplots_adjust(left=0.15)
        #ax.title(title)

        if return_plt:

            return fig, ax

        fig.savefig(name, dpi=1000)
        plt.clf()


def _make_1dhist(data, title, filename, prefix, return_plt=False):

    hist_object = hist1d_it(x=data, num_bins=25)
    plt.plot(hist_object.xx,
            hist_object.f)
    plt.hist(hist_object.x,
            hist_object.xx,
            normed=True,
            facecolor='green',
            alpha=0.5)
    plt.title(title)
    if not return_plt:
        name = filename
        if prefix!="":
            name = prefix+"_"+name
        plt.savefig(name, dpi=1000)
        del hist_object
        plt.close('all')


def _make_2dhist(data, dims, ax_name, title, filename, prefix, return_plt=False):

    hist_object = hist2d_it(data,
                            bins=[int(dims[0]*5),int(dims[1]*5)], 
                            xlim=[0.,dims[0]], 
                            ylim=[0.,dims[1]])
    hist_object.hist_it()

    fig1       = plt.figure()
    ax1        = fig1.gca()
    plt.imshow(hist_object.hist,
               interpolation="nearest",
               extent=[hist_object.xedges[0],\
                       hist_object.xedges[-1],\
                       hist_object.yedges[0],\
                       hist_object.yedges[-1]],
               origin='low',
               cmap=cm.jet)
    cfset = ax1.contourf(hist_object.xedges[:-1], 
                         hist_object.yedges[:-1], 
                         hist_object.hist, 
                         levels=np.linspace(0.,1.,101),
                         cmap=cm.jet)
    cbar  = plt.colorbar(cfset, label="Relative Frequency")
    plt.xlabel(r"$%s [\AA]$" %ax_name[0])
    plt.ylabel(r"$%s [\AA]$" %ax_name[1])
    plt.title(title)
    if not return_plt:
        name = filename
        if prefix!="":
            name = prefix+"_"+name
        plt.savefig(name, dpi=1000)

        del hist_object
        plt.close('all')

    else:
        return fig1, ax1