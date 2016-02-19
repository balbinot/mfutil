# -*- coding: utf-8 -*-

import numpy as np

class pwpl():
    """A piece-wise powerlaw mass function f(x) ~ x^-a

    Args:
        a (list): the powerlaw exponents for each piece
        mlim (list): the mass ranges for each exponent

    NOT:
       default parameters to Kroupa IMF

    Returns:
        Object

    """

    def __init__(self, a=[1.3, 2.3, 2.30], mlim=[0.08, 0.5, 2.0, 120.0], normalize=True):

        self._a = a
        self._mlim = mlim

        norm = np.zeros(len(a))
        aa = np.zeros(len(a))
        area = np.zeros(len(a))

        ## Assure piecewise continuity
        aa[0] = 1.0
        if len(a) > 1.0:
            for i in range(len(a))[1:]:
                aa[i] = aa[i-1] * np.power(mlim[i], a[i] - a[i-1])

        # Loop through pieces
        for i in range(len(a)):
            if a[i] != 1.0:
               norm[i] = (np.power(mlim[i+1], 1.0 - a[i]) -
                          np.power(mlim[i], 1.0 - a[i]))/(1.0 - a[i])
            else:
                np.log(mlim[i+1]) - np.log(mlim[i])

            area[i] = aa[i] * norm[i]

        #get area for each piece
        # aa[0] /= np.sum(area)
        for i in range(len(aa)):
            aa[i] /= np.sum(area)
        self._aa = aa

        #now normalize the PDF to unity/normalization
        for i in range(len(a)):
            area[i] *= aa[0]*normalize

        # this is all needed for sampling
        self._norm = 1.0/norm # scale
        self._area = area # piece information


    def sample(self, n):
        """Samples n times from the mass function.

        Args:
            n (int): number of samples

        Return:
            array (numpy.ndarray): array containing the mass values sampled.

        """

        area = self._area
        norm = self._norm
        a = self._a
        mlim = self._mlim

        x = np.random.rand(n)
        mass = np.array([])

        if len(a) == 1.0:
            i = 0
            j = (x <= area[i])
            nseg = len(x[j])
            if nseg > 0.0:
                m = self._getmass(np.random.rand(nseg), a[i], norm[i], mlim[i],
                                  mlim[i+1])
                mass = np.r_[mass, m]
        else:
            for i in range(len(a)):
                j = (x <= np.sum(area[0:i+1])) & (x > np.sum(area[0:i]))
                nseg = len(x[j])
                if nseg > 0.0:
                    m = self._getmass(np.random.rand(nseg), a[i], norm[i],
                                      mlim[i], mlim[i+1])
                    mass = np.r_[mass, m]

        return mass

    def _getmass(self, x, slope, norm, xmin, xmax):

        mass = (1.0 - slope)*x/norm + xmin**(1.0 - slope)
        Z = 1.0/(1.0-slope)
        mass = mass**Z
        return mass

    def eval(self, m):

        if ~isinstance(m, (np.ndarray, np.generic)):
            m = np.array(m)

        a = self._a
        mlim = self._mlim
        aa = self._aa

        j = (m <= mlim[-1])&(m >= mlim[0])
        m = m[j]

        F = np.array([])
        X = np.array([])
        for i in range(len(a)):
            j = (m <= mlim[i+1])&(m >= mlim[i])
            x = m[j]
            f = np.power(x, -a[i])*(aa[i])
            X = np.r_[X, x]
            F = np.r_[F, f]

        return X, F

if __name__ == "__main__":
    pass
