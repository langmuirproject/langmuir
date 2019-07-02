"""
Copyright 2018
    Sigvald Marholm <marholm@marebakken.com>
    Diako Darian <diakod@math.uio.no>

This file is part of langmuir.

langmuir is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

langmuir is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with langmuir.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division
from langmuir.analytical import *
from langmuir.geometry import *
from langmuir.species import *
from langmuir.misc import *
from scipy.interpolate import griddata
from scipy.constants import value as constants
from copy import deepcopy
import numpy as np
import os

def finite_length_current_density(geometry, species, V=None, eta=None,
                                  z=None, zeta=None, normalize=None):

    if isinstance(species, list):
        if normalize == True:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        i = 0
        for s in species:
            i += finite_length_current_density(geometry, s, V, eta, z, zeta)
        return i

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    k = constants('Boltzmann constant')

    if kappa != float('inf') or alpha != 0:
        logger.error("Finite length effect data only available for Maxwellian")

    if V is not None:
        eta_isarray = isinstance(V, (np.ndarray, list, tuple))
        V = make_array(V)
        eta = q*V/(k*T) # V must be array (not list) to allow this
    else:
        eta_isarray = isinstance(eta, (np.ndarray, list, tuple))
        eta = make_array(eta)

    eta = eta[:, None] # Make eta rows

    if not isinstance(geometry, Cylinder):
        raise ValueError('Geometry not supported: {}'.format(geometry))

    if zeta is None:
        zeta = 0.5*geometry.l/species.debye

    if z is not None:
        zeta_isarray = isinstance(z, (np.ndarray, list, tuple))
        z = make_array(z)
        zeta = z/species.debye # z must be array (not list) to allow this
    else:
        zeta_isarray = isinstance(zeta, (np.ndarray, list, tuple))
        zeta = make_array(zeta)

    lambd_p = geometry.l/species.debye      # Normalized probe length
    lambd_l = geometry.lguard/species.debye # Normalized left guard length
    lambd_r = geometry.rguard/species.debye # Normalized right guard length
    lambd_t = lambd_l + lambd_p + lambd_r   # Normalized total length

    # fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'cache.npz')
    # file = np.load(fname)
    # lambds = file['ls']
    # etas = file['etas']
    # As = file['popts'][:,0]
    # alphas = file['popts'][:,1]
    # gammas = file['popts'][:,4]
    # Cs = file['popts'][:,5]

    # A     = griddata((lambds, etas), As    , (lambd_t, eta))
    # alpha = griddata((lambds, etas), alphas, (lambd_t, eta))
    # gamma = griddata((lambds, etas), gammas, (lambd_t, eta))
    # C     = griddata((lambds, etas), Cs    , (lambd_t, eta))

    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params.npz')
    file = np.load(fname)
    lambds = file['lambds']
    etas = file['etas']
    Cs = file['Cs']
    As = file['As']
    alphas = file['alphas']
    deltas = file['deltas']

    lambd_coeff = min(lambd_t, max(lambds))

    print("ETA: {}".format(eta.shape))

    C     = griddata((lambds, etas), Cs    , (lambd_coeff, eta))
    A     = griddata((lambds, etas), As    , (lambd_coeff, eta))
    alpha = griddata((lambds, etas), alphas, (lambd_coeff, eta))
    delta = griddata((lambds, etas), deltas, (lambd_coeff, eta))

    print("ETA: {}, C: {}".format(eta.shape, C.shape))
    print(C)

    # def f(z):
    #     return A*np.exp(-alpha*z)*(z**gamma)

    # def g(z, l):
    #     return C+f(z)+f(l-z)

    if normalize=='th': # i0 = i_th => i = i_th * g
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        i0 = OML_current(geonorm, species, eta=eta, normalize=True)
    elif normalize=='OML': # i0 = 1 => i = g
        i0 = 1
    else: # i0 = i_OML => i = i_OML * g
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        i0 = OML_current(geonorm, species, eta=eta)

    # if lambd_l==float('inf'):
    #     i = i0*additive_model_noleft(zeta, lambd_p+lambd_r, C, A, alpha, delta)
    # else:
    #     # i = i0*g(zeta, lambd)
    #     # i = i0*additive_model(zeta, lambd, A, alpha, 0, 1, gamma, C)
    #     i = i0*additive_model(lambd_l+zeta, lambd_t, C, A, alpha, delta)

    # return i
    # return i[0] if len(i) == 1 else i

    def h(zeta):
        """
        Implements the function

            h(zeta) = A*(zeta-delta+alpha**(-1))*np.exp(-alpha*zeta),

        with h(inf)=0. Supports matrix operations.
        """
        res = np.zeros((len(alpha), len(zeta)))
        ind = np.where(zeta!=np.inf)[0] # res=0 where zeta=inf
        zeta = zeta[ind]
        res[:,ind] = A*(zeta-delta+alpha**(-1))*np.exp(-alpha*zeta)
        return res

    # def g(zeta_p, lambd_l, lambd_p, lambd_r):
    #     """
    #     Implements the function

    #         g(zeta) = C*(1+h(zeta)+h(lambda-zeta)).

    #     with lambda = lambda_l + lambda_p + lambda_r where the terms
    #     represent the left guard, the probe, and the right guard,
    #     respectively, and zeta_p = zeta + lambda_l is the position
    #     on the probe, and zeta is the position from the leftmost point.
    #     This decomposition allows the left/right edge-function h to be
    #     written independetly on lambda_r/lambda_l, which would cause
    #     loss of precision as the guard tends to infinity.
    #     """
    #     return C*(1+h(lambd_l+zeta_p)+h(lambd_p+lambd_r-zeta_p))

    # i = i0*g(zeta, lambd_l, lambd_p, lambd_r)

    g = C*(1+h(lambd_l+zeta)+h(lambd_p+lambd_r-zeta))
    i = i0*g

    # i = i0*additive_model(lambd_l+zeta, lambd_t, C, A, alpha, delta)

    if zeta_isarray:
        if eta_isarray:
            return i
        else:
            return i.ravel()
    else:
        if eta_isarray:
            return i.ravel()
        else:
            return i[0][0]

def finite_length_current(geometry, species,
                          V=None, eta=None, normalize=None):

    if isinstance(species, list):
        if normalize not in (False, None):
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        I = 0
        for s in species:
            I += finite_length_current(geometry, s, V, eta)
        return I

    q, m, n, T = species.q, species.m, species.n, species.T
    kappa, alpha = species.kappa, species.alpha
    k = constants('Boltzmann constant')

    if kappa != float('inf') or alpha != 0:
        logger.error("Finite length effect data only available for Maxwellian")

    # if V is not None:
    #     eta = q*V/(k*T)

    if V is not None:
        V = make_array(V)
        eta = q*V/(k*T)
    else:
        eta = make_array(eta)

    if not isinstance(geometry, Cylinder):
        raise ValueError('Geometry not supported: {}'.format(geometry))

    lambd_p = geometry.l/species.debye      # Normalized probe length
    lambd_l = geometry.lguard/species.debye # Normalized left guard length
    lambd_r = geometry.rguard/species.debye # Normalized right guard length
    lambd_t = lambd_l + lambd_p + lambd_r   # Normalized total length

    # fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'cache.npz')
    # file = np.load(fname)
    # lambds = file['ls']
    # etas = file['etas']
    # As = file['popts'][:,0]
    # alphas = file['popts'][:,1]
    # gammas = file['popts'][:,4]
    # Cs = file['popts'][:,5]

    # A     = griddata((lambds, etas), As    , (lambd_t, eta))
    # alpha = griddata((lambds, etas), alphas, (lambd_t, eta))
    # gamma = griddata((lambds, etas), gammas, (lambd_t, eta))
    # C     = griddata((lambds, etas), Cs    , (lambd_t, eta))

    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params.npz')
    file = np.load(fname)
    lambds = file['lambds']
    etas = file['etas']
    Cs = file['Cs']
    As = file['As']
    alphas = file['alphas']
    deltas = file['deltas']

    lambd_coeff = min(lambd_t, max(lambds))

    C     = griddata((lambds, etas), Cs    , (lambd_coeff, eta))
    A     = griddata((lambds, etas), As    , (lambd_coeff, eta))
    alpha = griddata((lambds, etas), alphas, (lambd_coeff, eta))
    delta = griddata((lambds, etas), deltas, (lambd_coeff, eta))

    if normalize=='th': # I = actual current / I_th
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        I0 = OML_current(geonorm, species, eta=eta, normalize=True)/geometry.l
        # Notice the difference between i_th [A/m] and I_th [A]
        # geonorm.l = 1 makes it per unit length, i.e. [A/m].
        # Hence OML_current returns division by i_th and not I_th.
        # To correct this we need to divide by geometry.l
    elif normalize=='OML': # I = integral of g
        I0 = 1/geometry.l
    else: # I0 = I_OML => I = I_OML * integral of g = actual current
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        I0 = OML_current(geonorm, species, eta=eta)

    # I = I0*species.debye*(int_additive_model(lambd_l+lambd_p, lambd_t, A, alpha, 0, 1, gamma, C)
    #                      -int_additive_model(lambd_l   , lambd_t, A, alpha, 0, 1, gamma, C))

    # if lambd_l==float('inf'):
    #     I = I0*species.debye*(int_additive_model_noleft(lambd_p, lambd_p+lambd_r, C, A, alpha, delta)
    #                          -int_additive_model_noleft(0      , lambd_p+lambd_r, C, A, alpha, delta))
    # else:
    #     I = I0*species.debye*(int_additive_model(lambd_l+lambd_p, lambd_t, C, A, alpha, delta)
    #                          -int_additive_model(lambd_l        , lambd_t, C, A, alpha, delta))

    # I = I0*species.debye*(int_additive_model(lambd_l+lambd_p, lambd_t, C, A, alpha, delta)
    #                      -int_additive_model(lambd_l        , lambd_t, C, A, alpha, delta))

    def H(zeta):
        if zeta==float('inf'): return np.zeros_like(alpha)
        return A*(alpha*(delta-zeta)-2)*np.exp(-alpha*zeta)/alpha**2

    int_g = C*(lambd_p+H(lambd_p+lambd_l)+H(lambd_p+lambd_r)-H(lambd_l)-H(lambd_r))
    I = I0*species.debye*int_g

    # return I
    return I[0] if len(I) == 1 else I

# def Gamma(a, x):
#     return special.gammaincc(a, x)*special.gamma(a)

# def h(zeta, alpha, gamma):
#     return np.exp(-alpha*zeta)*(zeta**gamma)

# # Indefinite integral of h
# def H(zeta, alpha, gamma):
#     if zeta==0: zeta=np.finfo(float).eps
#     return -(zeta**gamma)*((alpha*zeta)**(-gamma))*Gamma(1+gamma,alpha*zeta)/alpha

# def additive_model(zeta, lambd, A, alpha, B, beta, gamma, C):
#     return C + A*h(zeta, alpha, gamma) + B*h(zeta, beta, gamma) \
#              + A*h(lambd-zeta, alpha, gamma) + B*h(lambd-zeta, beta, gamma)

# def int_additive_model(zeta, lambd, A, alpha, B, beta, gamma, C):
#     return C*zeta + A*H(zeta, alpha, gamma) + B*H(zeta, beta, gamma) \
#                   - A*H(lambd-zeta, alpha, gamma) - B*H(lambd-zeta, beta, gamma)

def h(zeta, alpha, delta):
    res = np.zeros((len(alpha), len(zeta)))
    ind = np.where(zeta!=np.inf)[0] # res=0 where zeta=inf
    zeta = zeta[ind]
    res[:,ind] = (zeta-delta+alpha**(-1))*np.exp(-alpha*zeta)
    return res

# def H(zeta, alpha, delta):
#     res = np.zeros_like(alpha)
#     ind = np.where(zeta!=np.inf)[0]

def H(zeta, alpha, delta):
    if zeta==float('inf'): return np.zeros_like(alpha)
    return (alpha*(delta-zeta)-2)*np.exp(-alpha*zeta)/alpha**2

def additive_model(zeta, lambd, C, A, alpha, delta):
    return C * (1 + A*h(zeta, alpha, delta) + A*h(lambd-zeta, alpha, delta))

# def additive_model_noleft(zeta, lambd, C, A, alpha, delta):
#     return C * (1 + A*h(lambd-zeta, alpha, delta))

def int_additive_model(zeta, lambd, C, A, alpha, delta):
    print("zeta={}, labmda={}".format(zeta, lambd))
    return C * (zeta + A*H(zeta, alpha, delta) - A*H(lambd-zeta, alpha, delta))

# def int_additive_model_noleft(zeta, lambd, C, A, alpha, delta):
#     return C * (zeta - A*H(lambd-zeta, alpha, delta))
