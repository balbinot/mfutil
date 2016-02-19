# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import quad

class TaperedExp():
    """An tapered exponential mass function f(x) ~ x^a * (1 - Exp(-(x/mc)**b)) """

    def __init__(self, a=-2.3, b=1.3, mc=0.33, mmin=0.08, mmax=120.0, normalize=True):
        """
        NAME:
            __init__
        PURPOSE:
            Initialize a tapered exponential DF f(x) ~ x^a * (1 - Exp(-(x/mc)**b))
        INPUT:
            normalize -- if True, normalize the DF to 1 (one).-
                         If number, normalize to number.
            a -- float (default = -2.0)
            b -- float (default = 2.5)
            mc -- float (default = 0.15)
            mmin, mmax -- mass range (large range equals low performance)
        OUTPUT:
        HISTORY:
            2011-03-30 - Written - Balbinot (UFRGS)
        """
        self._a = a
        self._b = b
        self._mc = mc
        self.mmin = mmin
        self.mmax = mmax
        self.normalize = normalize

        if normalize==True:
            self._norm = self._normalize(1)
        else:
            self._norm = self._normalize(self.normalize)
        return None

    def eval(self, x):
        """
        NAME:
            eval
        PURPOSE:
            Evaluate the normalized DF at x
        INPUT:
            x -- number
        OUTPUT:
            DF(x) -- number-
        HISTORY:
            2011-03-30 - Written - Balbinot (UFRGS)
        """

        return self._norm*((x**self._a)*(1.0 - np.exp(-(x/self._mc)**self._b)))

    def _integrand(self, x):
        """ Whatenever needed for _normalize """

        return ((x**self._a)*(1 - np.exp(-(x/self._mc)**self._b)))

    def _normalize(self, N):
        """
        NAME
            _normalize
        PURPOSE:
            Normalize the DF
        INPUT:
        OUTPUT:
            number
        HISTORY:
            2011-03-30 - Written - Balbinot (UFRGS)
        """
        norm = N
        norm /= quad(self._integrand, self.mmin, self.mmax)[0]
        return norm
