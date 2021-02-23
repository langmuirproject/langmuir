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
from langmuir.finite_length import *
from scipy.optimize import root_scalar

def generate_synthetic_data(geometry,
                            V,
                            model=finite_length_current,
                            V0=None,
                            alt_range=(100,500),
                            noise=1e-5,
                           ):

    fname = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'iri.npz')

    iri_data = np.load(fname)

    res = {}
    ind, = np.where(np.logical_and(iri_data['alt']>=alt_range[0],
                                   iri_data['alt']<=alt_range[1]))
    res['alt'] = iri_data['alt'][ind]
    res['Te'] = iri_data['Te'][ind]
    res['Ti'] = iri_data['Ti'][ind]
    res['ne'] = iri_data['ne'][ind]
    res['nO+'] = iri_data['nO'][ind]
    res['nO2+'] = iri_data['nO2'][ind]
    res['nNO+'] = iri_data['nNO'][ind]
    res['nH+'] = iri_data['nH'][ind]

    V = make_array(V)

    I = np.zeros((len(res['alt']), len(V)))
    V0_ = np.ones(len(res['alt']))
    if V0 is not None:
        V0_ *= V0

    def cap(x):
        if np.isnan(x):
            return 1 # FIXME: Finite-length model does not support zero density
        else:
            return x

    for i in range(len(res['alt'])):

        Te = res['Te'][i]
        Ti = res['Ti'][i]

        plasma = [Electron(n=cap(res['ne'][i]), T=Te)
                 ,Oxygen(n=cap(res['nO+'][i]), T=Ti)
                 ,Oxygen(n=cap(res['nO2+'][i]), T=Ti, Z=2)
                 ,Species(n=cap(res['nNO+'][i]), T=Ti, amu=30)
                 ,Hydrogen(n=cap(res['nH+'][i]), T=Ti)
                 ]

        if V0 is None:
            # This is just to have some varying number for the floating potential.
            # It doesn't have to be realistic, just within a resonable value
            residual_geometry = Sphere(r=10*debye(plasma))
            def residual(V):
                return OML_current(residual_geometry, plasma, V)
            sol = root_scalar(residual, x0=-10, x1=0)
            V0_[i] = sol.root

        try:
            I[i] = model(geometry, plasma, V=V0_[i]+V)
        except:
            I[i] = float('nan')

    I += noise*np.sqrt(np.abs(I))*np.random.randn(*I.shape)

    res['I'] = I
    res['V0'] = V0_

    return res

