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
from scipy.interpolate import griddata, RectBivariateSpline, interp1d
from scipy.constants import value as constants
from copy import deepcopy
import numpy as np
import os

def finite_length_current_density(geometry, species, V=None, eta=None,
                                  z=None, zeta=None, normalization=None):
    """
    Current collected per unit length on a cylindrical probe according to the
    Marholm-Marchand finite-length model. Assumes small radius compared to the
    Debye length.

    Parameters
    ----------
    geometry: Cylinder
        Probe geometry

    species: Species or array-like of Species
        Species constituting the background plasma

    V: float or array-like of floats
        Probe voltage(s) in [V]. Overrides eta.

    eta: float or array-like of floats
        Probe voltage(s) normalized by k*T/q, where q and T are the species'
        charge and temperature and k is Boltzmann's constant.

    z: float or array-like of floats
        Position(s) along the probe in [m]. Overrides zeta.

    zeta: float or array-like of floats
        Position(s) along the probe normalized by the Debye length.

    normalization: 'th', 'thmax', 'oml', None
        Wether to normalize the output current per unit length by,
        respectively, the thermal current per unit length, the Maxwellian
        thermal current per unit length, the OML current per unit length, or
        not at all, i.e., current in [A/m].

    Returns
    -------
    float if voltage and position are floats. 1D array of floats corresponding
    to voltage or position if one of them is array-like. 2D array of floats if
    voltage and position are both array-like, one row per voltage.
    """

    if isinstance(species, list):
        if normalization is not None:
            logger.error('Cannot normalize current to more than one species')
            return None
        if eta is not None:
            logger.error('Cannot normalize voltage to more than one species')
            return None
        if zeta is not None:
            logger.error('Cannot normalize position to more than one species')
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
        eta = -q*V/(k*T) # V must be array (not list) to allow this
    else:
        eta_isarray = isinstance(eta, (np.ndarray, list, tuple))
        eta = make_array(eta)

    # eta = eta[:, None] # Make eta rows

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

    C, A, alpha, delta = get_lerped_coeffs(lambd_t, eta)

    # Make these column vectors, i.e. one eta per row
    C = C[:, None]
    A = A[:, None]
    alpha = alpha[:, None]
    delta = delta[:, None]
    eta = eta[:, None]

    if normalization is None: # i0 = i_OML => i = i_OML * g
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        i0 = OML_current(geonorm, species, eta=eta)
    elif normalization.lower() in ['th', 'thmax']: # i0 = i_th => i = i_th * g
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        i0 = OML_current(geonorm, species, eta=eta, normalization='th')
    elif normalization.lower()=='oml': # i0 = 1 => i = g
        i0 = 1
    else:
        raise ValueError('Normalization not supported: {}'.format(normalization))

    def h(zeta):
        res = np.zeros((len(alpha), len(zeta)))
        ind = np.where(zeta!=np.inf)[0] # res=0 where zeta=inf
        zeta = zeta[ind]
        res[:,ind] = A*(zeta-delta+alpha**(-1))*np.exp(-alpha*zeta)
        return res

    g = C*(1+h(lambd_l+zeta)+h(lambd_p+lambd_r-zeta))
    i = i0*g

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
                          V=None, eta=None, normalization=None,
                          interpolate='I'):
    """
    Current collected by a cylindrical probe according to the Marholm-Marchand
    finite-length model. Assumes small radius compared to the Debye length.

    Parameters
    ----------
    geometry: Cylinder
        Probe geometry

    species: Species or array-like of Species
        Species constituting the background plasma

    V: float or array-like of floats
        Probe voltage(s) in [V]. Overrides eta.

    eta: float or array-like of floats
        Probe voltage(s) normalized by k*T/q, where q and T are the species'
        charge and temperature and k is Boltzmann's constant.

    normalization: 'th', 'thmax', 'oml', None
        Whether to normalize the output current by, respectively, the thermal
        current, the Maxwellian thermal current, the OML current, or not at
        all, i.e., current in [A/m].

    interpolate: 'I', 'g'
        Whether to interpolate the coefficients of the profile function g and
        then integrate g to get the current (faster), or if g is integrated it
        its present grid to get a grid of currents which can be interpolated
        from. This makes the interpolation linear in I, and avoids the
        irregular behaviour sometimes experienced for shorter probes due to
        irregularities in the coefficients.

    Returns
    -------
    float if voltage is float. array of floats corresponding to voltage if
    voltage is array-like.
    """

    if isinstance(species, list):
        if normalization is not None:
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

    if V is not None:
        eta_isarray = isinstance(V, (np.ndarray, list, tuple))
        V = make_array(V)
        eta = -q*V/(k*T)
    else:
        eta_isarray = isinstance(eta, (np.ndarray, list, tuple))
        eta = make_array(eta)

    if(np.any(eta>100.)):
        logger.warning('Finite-length theory yields erroneous results for eta>100')

    if not isinstance(geometry, Cylinder):
        raise ValueError('Geometry not supported: {}'.format(geometry))

    lambd_p = geometry.l/species.debye      # Normalized probe length
    lambd_l = geometry.lguard/species.debye # Normalized left guard length
    lambd_r = geometry.rguard/species.debye # Normalized right guard length
    lambd_t = lambd_l + lambd_p + lambd_r   # Normalized total length

    if normalization is None: # I0 = I_OML => I = I_OML * integral of g = actual current
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        I0 = OML_current(geonorm, species, eta=eta)
    elif normalization.lower() in ['th', 'thmax']: # I = actual current / I_th
        geonorm = deepcopy(geometry)
        geonorm.l = 1
        I0 = OML_current(geonorm, species, eta=eta, normalization='th')/geometry.l
        # Notice the difference between i_th [A/m] and I_th [A]
        # geonorm.l = 1 makes it per unit length, i.e. [A/m].
        # Hence OML_current returns division by i_th and not I_th.
        # To correct this we need to divide by geometry.l
    elif normalization.lower()=='oml': # I = integral of g
        I0 = (1/geometry.l)*np.ones_like(eta)
    else:
        raise ValueError('Normalization not supported: {}'.format(normalization))

    if interpolate.lower() == 'g':
        C, A, alpha, delta = get_lerped_coeffs(lambd_t, eta)

        def H(zeta):
            if zeta==float('inf'): return np.zeros_like(alpha)
            return A*(alpha*(delta-zeta)-2)*np.exp(-alpha*zeta)/alpha**2

        int_g = C*(lambd_p+H(lambd_p+lambd_l)+H(lambd_p+lambd_r)-H(lambd_l)-H(lambd_r))
        I = I0*species.debye*int_g

    elif interpolate.lower() in ['i', 'i2']:

        fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params_structured.npz')
        file = np.load(fname)

        lambds = file['lambds']
        etas = file['etas']

        ind = np.where(lambds<lambd_t)[0][-1]

        if ind==len(lambds)-1: # extrapolate from highest lambda

            lambds = lambds[ind]
            As     = file['As'][ind]
            Cs     = file['Cs'][ind]
            alphas = file['alphas'][ind]
            deltas = file['deltas'][ind]

            def H(zeta):
                res = np.zeros_like(alphas)
                if(zeta==float('inf')):
                    return np.zeros_like(alphas)
                else:
                    return As*(alphas*(deltas-zeta)-2)*np.exp(-alphas*zeta)/alphas**2

            int_gs = Cs*(lambd_p+H(lambd_p+lambd_l)+H(lambd_p+lambd_r) \
                                -H(lambd_l)-H(lambd_r))

        else: # interpolate between lambdas

            As       = file['As'][ind:ind+2]
            Cs       = file['Cs'][ind:ind+2]
            alphas   = file['alphas'][ind:ind+2]
            deltas   = file['deltas'][ind:ind+2]

            # Stretch whole probe, including guards
            lambd_ts = lambds[ind:ind+2]
            lambd_ls = lambd_l*lambd_ts/lambd_t
            lambd_ps = lambd_p*lambd_ts/lambd_t
            lambd_rs = lambd_r*lambd_ts/lambd_t

            # Stretch only probe segment, constant guards
            # lambd_ts = lambds[ind:ind+2]
            # lambd_ps = lambd_ts-lambd_l-lambd_r
            # lambd_ls = np.array([lambd_l, lambd_l])
            # lambd_rs = np.array([lambd_r, lambd_r])

            # if lambd_ps[0]<0:
            #     # Probe segment can't be shorter than zero.
            #     # Let it be zero, and (rightly) let it collect zero current.
            #     lambd_ps[0] = 0
            #     lambd_ts[0] = lambd_l+lambd_r
            #     Cs[0] = 0
            #     print("WARNING")

            def H(zetas):
                return As*(alphas*(deltas-zetas[:,None])-2) \
                         *np.exp(-alphas*zetas[:,None])/alphas**2

            int_gs = Cs*( lambd_ps[:,None]     \
                         +H(lambd_ps+lambd_ls) \
                         +H(lambd_ps+lambd_rs) \
                         -H(lambd_ls)          \
                         -H(lambd_rs))

            weight = (lambd_ps[1]-lambd_p)/(lambd_ps[1]-lambd_ps[0])
            int_gs = weight*int_gs[0] + (1-weight)*int_gs[1]

        attracted = np.where(eta>=0.)[0]
        repelled  = np.where(eta< 0.)[0]
        over      = np.where(eta>100.)[0]

        eta[over] = 100

        I = I0*species.debye*np.ones_like(eta)
        func = interp1d(etas, int_gs)

        I[attracted] *= func(eta[attracted])
        I[repelled]  *= lambd_p

    else:
        raise ValueError("interpolate must be either 'g' or 'I'")

    return I if eta_isarray else I[0]

# def get_lerped_coeffs(lambd, eta):
#     """
#     Fetches and interpolates Marholm-Marchand coefficients.

#     Parameters
#     ----------
#     lambd: float
#         Normalized probe length (lambda) to get coefficients for

#     eta: float
#         Normalized probe voltage (eta) to get coefficients for

#     Returns
#     -------
#     4-tuple of coefficients (C, A, alpha, delta)
#     """

#     fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params.npz')
#     file = np.load(fname)

#     lambds = file['lambds']
#     etas = file['etas']
#     Cs = file['Cs']
#     As = file['As']
#     alphas = file['alphas']
#     deltas = file['deltas']

#     # Extrapolate for larger lambda by using the largest available
#     lambd_coeff = min(lambd, max(lambds))

#     # Tabulated values contains datapoints as described in paper
#     C     = griddata((lambds, etas), Cs    , (lambd_coeff, eta))
#     A     = griddata((lambds, etas), As    , (lambd_coeff, eta))
#     alpha = griddata((lambds, etas), alphas, (lambd_coeff, eta))
#     delta = griddata((lambds, etas), deltas, (lambd_coeff, eta))

#     # Make repelled species identical to OML through these coefficients
#     ind = np.where(eta>=0)[0]
#     C[ind] = 1
#     A[ind] = 0
#     alpha[ind] = 1
#     delta[ind] = 1

#     return C, A, alpha, delta

class lerper(RectBivariateSpline):
    def __init__(self, coeffname, degree=1):

        assert coeffname in ['C', 'A', 'alpha', 'delta']
        coeffname += 's'

        fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params_structured.npz')
        file = np.load(fname)

        lambds = file['lambds']
        etas = file['etas']
        coeff = file[coeffname]

        super(lerper, self).__init__(lambds, etas, coeff, kx=degree, ky=degree)

def get_max_lambd():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),'params_structured.npz')
    file = np.load(fname)
    lambds = file['lambds']
    return max(lambds)

lerp_C = lerper('C')
lerp_A = lerper('A')
lerp_alpha = lerper('alpha')
lerp_delta = lerper('delta')
max_lambd = get_max_lambd()

def get_lerped_coeffs_new(lambd, eta):

    # Extrapolate for larger lambda by using the largest available
    lambd_coeff = min(lambd, max_lambd)

    # Make repelled species identical to OML through these coefficients
    C     = np.ones_like(eta)
    A     = np.zeros_like(eta)
    alpha = np.ones_like(eta)
    delta = np.ones_like(eta)

    # Tabulated values contains datapoints as described in paper
    ind = np.where(eta>0.0)[0]
    C[ind] = lerp_C(lambd_coeff, eta[ind], grid=False)
    A[ind] = lerp_A(lambd_coeff, eta[ind], grid=False)
    alpha[ind] = lerp_alpha(lambd_coeff, eta[ind], grid=False)
    delta[ind] = lerp_delta(lambd_coeff, eta[ind], grid=False)

    return C, A, alpha, delta

get_lerped_coeffs=get_lerped_coeffs_new
