#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pwpl import pwpl
from tapered_exp import TaperedExp

__all__ = ['pwpl', 'TaperedExp']

if __name__=='__main__':
    import matplotlib.pyplot as p
    import numpy as np

    mf = TaperedExp(normalize=1)
    mf2 = pwpl(normalize=1)

    x = np.arange(0.08,120,0.01)
    y = mf.eval(x)
    _, y2 = mf2.eval(x)

#    p.plot(x,y*x*np.log(10), 'k-')
#    p.plot(x,y2*x*np.log(10), 'r-')
    tm = 1.6e4
    p.plot(x,tm*y, 'k-')
    p.plot(x,tm*y2, 'r-')

    from scipy.integrate import quad

    tmp = lambda x: x*tm*mf.eval(x)
    tmp2 = lambda x: tm*mf.eval(x)
    mass = quad(tmp, 0.08, 120)[0]
    N = quad(tmp2, 0.08, 120)[0]
    print "Mass %lf; Number %lf" % (mass, N)

    p.loglog()
    p.show()

