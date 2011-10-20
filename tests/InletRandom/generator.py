#!/usr/bin/env python

class RandInlet(object):
  r"""Generate random variables from the PDF

  .. math:: f_x = \frac{X-a}{N} * e^{-X^2*b^2}

  where

  .. math:: N = b^2 e^{-\frac{a^2}{2 b^2}} + a \sqrt{\frac{\pi b^2}{2}}
                \left(1+\mathrm{erf}\left(\frac{a}{\sqrt{2 b^2}}\right)\right)

  using numerical inversion of the CDF with Newton-Raphson.

  Example
  -------
  >>> import matplotlib.pyplot as plt
  >>> r = RandOutlet(10, 0.5)
  >>> a = [r() for i in range(100000)]
  >>> plt.hist(a, bins=100, normed=True

  Note
  -----
  Using iterative methods is generally slow and inaccurate. It would be better
  to use an efficient tabulation algorithm, such as the one described by
  Hoermann and Leydold (2003) in http://dx.doi.org/10.1145/945511.945517. This
  algorithm is implemented e.g. in UNU.RAN (http://statmath.wu.ac.at/unuran).

  """
  def __init__(self, a, b):
    import numpy as np
    from scipy.special import erf
    self._a = a
    self._b = b
    # helpers
    abSqrtPi = a*b*np.sqrt(np.pi)
    erfab = erf(a*b)
    expa2b2 = np.exp(-a**2*b**2)
    denom = expa2b2+abSqrtPi*(1+erfab)
    b22 = 2*b**2
    # objective function (CDF - LHS)
    self._F = lambda x, U: (expa2b2-np.exp(-b**2*(x-a)**2)+abSqrtPi*(erfab+erf(b*(x-a))))/denom-U
    # derivative of objective function (PDF)
    b22denom = b22/denom
    self._f = lambda x, U: b22denom*np.exp(-b**2*(x-a)**2)*x
    # location of the peak value
    self._m = 0.5*(a*b+np.sqrt(2.+a**2*b**2))/b
    self._nCall = 0
    self._nIterTotal = 0
    self._nIterMax = 0
    self._nIterMin = 100000

  def _newton(self, func, x0, fprime, args=(), tol=1.48e-8, maxiter=50):
    self._nCall += 1
    p0 = 1.0 * x0
    for iter in range(maxiter):
      myargs = (p0,) + args
      fval = func(*myargs)
      fder = fprime(*myargs)
      if fder == 0:
        print "derivative was zero."
        return p0
      p = p0 - fval/fder
      if abs(p - p0) < tol:
        self._nIterTotal += iter+1
        self._nIterMax = max(self._nIterMax, iter+1)
        self._nIterMin = min(self._nIterMin, iter+1)
        return p
      p0 = p

  def __call__(self):
    import numpy as np
    import numpy.random as nr
    from scipy.optimize import newton
    U = nr.rand()
    # perturbation away from peak position
    d = -np.sign(self._F(self._m, U))*1e-12
    return self._newton(self._F, self._m+d, self._f, (U,))
