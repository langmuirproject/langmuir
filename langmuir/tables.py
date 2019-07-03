import numpy as np
from scipy.interpolate import griddata
from langmuir.species import *
from langmuir.geometry import *
from langmuir.analytical import *

# https://stackoverflow.com/questions/11144513/numpy-cartesian-product-of-x-and-y-array-points-into-single-array-of-2d-points
def cartesian_product(*arrays):
    """
    Takes the cartesian product of multiple numpy arrays. Example::

        >>> a = np.array([100,200])
        >>> b = np.array([10,20])
        >>> c = np.array([1,2,3])
        >>> cartesian_product(a, b, c)
        array([[100,  10,   1],
               [100,  10,   2],
               [100,  10,   3],
               [100,  20,   1],
               [100,  20,   2],
               [100,  20,   3],
               [200,  10,   1],
               [200,  10,   2],
               [200,  10,   3],
               [200,  20,   1],
               [200,  20,   2],
               [200,  20,   3]])
    """
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def get_coords_and_value(table, *indices):
    """
    Return the coordinates and value corresponding to a set of indices in a
    table. Example::

        >>> t = get_table('laframboise sphere')
        >>> x = get_coords_and_value(t, 4, 3)
        >>> np.allclose(x, [1, -0.6, 1.595])
        True
    """
    coords = [a[i] for a,i in zip(table['axes'], indices)]
    coords.append(table['values'][indices])
    return coords

def get_table(name, provide_points=True):
    """
    Returns tabulated normalized attracted-species currents for finite radii
    spherical or cylindrical probes. If you do not need the raw tables,
    consider using the higher level function finite_radius_current() instead.

    The following tables exists:

        - laframboise sphere
        - laframboise cylinder
        - darian-marholm uncomplete sphere
        - darian-marholm uncomplete cylinder
        - darian-marholm sphere
        - darian-marholm cylinder
        - laframboise-darian-marholm sphere
        - laframboise-darian-marholm cylinder

    The Laframboise tables are tables 5c and 6c in Laframboise's thesis.
    These only cover Maxwellian velocity distributions. The Darian-Marholm
    tables cover Kappa-Cairns (and subtype) distributions, although it is
    not as accurate and with as wide input domain for Maxwellian as Laframboise.
    The uncomplete version of Darian-Marholm tables are the ones presented in
    the Darian-Marholm paper, whereas the completed version is the same but
    with analytical values inserted for zero radius (OML theory) which are
    not covered by the uncomplete version, and analytical thermal currents for
    zero potential since these are not as accurate as the other currents in
    the uncomplete Darian-Marholm paper. The Laframboise-Darian-Marholm tables
    are the same as the complete Darian-Marholm tables except that the
    Maxwellian values are replaced by those of Laframboise, which are more
    accurate. These composed tables therefore offer the greatest accuracy and
    range for any case. Since the Laframboise case has a larger grid, the
    grid is no longer regular for the Laframboise-Darian-Marholm tables.

    The returned dictionary ``table`` has the following keys::

        - ``table['axes']`` is a tuple of lists, each list containing the
          grid values along that axis. The axes are, in this order, 1/kappa,
          alpha, R, eta, where kappa and alpha are the spectral indices of
          the Kappa-Cairns distribution, R is the probe length in terms of
          Debye lengths, and eta is the normalized voltage qV/kT. For
          the Laframboise tables the first to axes do not exist.

        - ``table['values']`` is a 4D array (2D for Laframboise) of values.
          ``table['values'][i][j][k][l] is the value corresponding to
          table['axes'][0][i], table['axes'][1][j], and so forth. For
          non-regular grids (the Laframboise-Darian-Marholm tables) it is
          flattened to 1D.

        - ``table['points']`` is a flatted array of 4-tuples (2-tuples for
          Laframboise tables), of 1/kappa, alpha, R, eta corresponding to
          the values in a flattened array. For non-regular grids,
          ``table['values']`` are already this flattened array, but in any
          case, ``table['values'].reshape(-1)`` is always correct. This
          is only provided when the grid is non-regular or when
          ``provide_points==True``.
    """

    name = name.lower()
    tol = 1e-6 # For float comparisons

    if name == 'laframboise sphere':

        # Table 5c in Laframboise's thesis

        Rs =  [0, 0.2, 0.3, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 50, 100]

        etas =[0.0  , 0.1  , 0.3  , 0.6  , 1.0  , 1.5  , 2.0  , 3.0  , 5.0  , 7.5  , 10.0  , 15.0  , 20.0  , 25.0]
        Is = [[1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 26.000],   # R=0
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 25.763],   # R=0.2
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.500, 3.000, 4.000, 6.000, 8.500, 11.000, 16.000, 21.000, 25.462],   # R=0.3
              [1.000, 1.100, 1.300, 1.600, 2.000, 2.493, 2.987, 3.970, 5.917, 8.324, 10.704, 15.403, 20.031, 24.607],   # R=0.5
              [1.000, 1.0999,1.299, 1.595, 1.987, 2.469, 2.945, 3.878, 5.687, 7.871,  9.990, 14.085, 18.041, 21.895],   # R=1
              [1.000, 1.0999,1.299, 1.584, 1.955, 2.399, 2.824, 3.632, 5.126, 6.847,  8.460, 11.482, 14.314, 17.018],   # R=2
              [1.000, 1.0999,1.293, 1.572, 1.922, 2.329, 2.709, 3.406, 4.640, 6.007,  7.258,  9.542, 11.636, 13.603],   # R=3
              [1.000, 1.099, 1.288, 1.552, 1.869, 2.219, 2.529, 3.086, 3.957, 4.887,  5.710,  7.167,  8.473,  9.676],   # R=5
              [1.000, 1.099, 1.288, 1.552, 1.869, 2.219, 2.529, 3.086, 3.957, 4.094,  4.658,  5.645,  6.518,  7.318],   # R=7.5
              [1.000, 1.098, 1.280, 1.518, 1.783, 2.050, 2.226, 2.609, 3.119, 3.620,  4.050,  4.796,  5.453,  6.053],   # R=10
              [1.000, 1.098, 1.280, 1.518, 1.783, 2.050, 2.226, 2.609, 3.119, 3.620,  4.050,  4.796,  4.318,  4.719],   # R=15
              [1.000, 1.097, 1.269, 1.481, 1.694, 1.887, 2.030, 2.235, 2.516, 2.779,  3.002,  3.383,  3.716,  4.018],   # R=20
              [1.000, 1.095, 1.255, 1.433, 1.592, 1.719, 1.803, 1.910, 2.037, 2.148,  2.241,  2.397,  2.532,  2.658],   # R=50
              [1.000, 1.094, 1.245, 1.402, 1.534, 1.632, 1.694, 1.762, 1.833, 1.891,  1.938,  2.022,  2.097,  2.166]]   # R=100

        Rs = np.array(Rs)
        etas = -np.array(etas)
        Is = np.array(Is)

        table = {'axes': (Rs, etas), 'values': Is}

    elif name == 'laframboise cylinder':

        # Table 6c in Laframboise's thesis

        Rs =  [0, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 20, 30, 40, 50, 100]

        etas =  [0.0  , 0.1   , 0.3   , 0.6   , 1.0   , 1.5   , 2.0   , 3.0   , 5.0   , 7.5   , 10.0  , 15.0  , 20.0  , 25.0]
        Is = [[1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.7388, 4.5114, 5.1695, 5.7526],   # R=0
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.7388, 4.5114, 5.1695, 5.7525],   # R=1
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.2417, 2.7555, 3.2846, 3.735 , 4.493 , 5.141 , 5.711 ],   # R=1.5
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.247 , 2.750 , 3.266 , 3.703 , 4.439 , 5.060 , 5.607 ],   # R=2
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.7551, 1.9320, 2.237 , 2.731 , 3.227 , 3.645 , 4.342 , 4.936 , 5.462 ],   # R=2.5
              [1.000, 1.0804, 1.2101, 1.3721, 1.5560, 1.754 , 1.928 , 2.226 , 2.701 , 3.174 , 3.567 , 4.235 , 4.789 , 5.291 ],   # R=3
              [1.000, 1.0804, 1.2101, 1.3721, 1.554 , 1.747 , 1.913 , 2.192 , 2.626 , 3.050 , 3.402 , 3.990 , 4.489 , 4.926 ],   # R=4
              [1.000, 1.0804, 1.2100, 1.371 , 1.549 , 1.735 , 1.893 , 2.151 , 2.544 , 2.920 , 3.231 , 3.749 , 4.183 , 4.565 ],   # R=5
              [1.000, 1.0803, 1.208 , 1.362 , 1.523 , 1.677 , 1.798 , 1.98  , 2.22  , 2.442 , 2.622 , 2.919 , 3.166 , 3.384 ],   # R=10
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.319 , 2.455 , 2.576 ],   # R=20
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.082 , 2.177 , 2.262 ],   # R=30
              [1.000, 1.0803, 1.205 , 1.348 , 1.486 , 1.605 , 1.689 , 1.801 , 1.940 , 2.060 , 2.157 , 2.082 , 2.025 , 2.092 ],   # R=40
              [1.000, 1.0803, 1.198 , 1.327 , 1.439 , 1.523 , 1.576 , 1.638 , 1.703 , 1.756 , 1.798 , 1.868 , 1.929 , 1.983 ],   # R=50
              [1.000, 1.0803, 1.194 , 1.314 , 1.409 , 1.478 , 1.518 , 1.561 , 1.599 , 1.628 , 1.650 , 1.686 , 1.719 , 1.747 ]]   # R=100

        Rs = np.array(Rs)
        etas = -np.array(etas)
        Is = np.array(Is)

        table = {'axes': (Rs, etas), 'values': Is}

    elif name == 'darian-marholm uncomplete sphere':

        Rs = [0.2, 1.0, 2.0, 3.0, 5.0, 10.0]
        etas = [0, 1, 2, 3, 5, 10, 15, 20, 25]
        kappa_recips = [0, 1.0/6, 1.0/4]
        alphas = [0, 0.2]
        Is = [[

            # kappa=inf
            # alpha=0.0

            [[ 0.965,  0.972,  0.982,  0.983, 1.008, 0.992],
             [ 2.019,  1.950,  1.958,  1.916, 1.851, 1.749],
             [ 3.021,  2.918,  2.813,  2.692, 2.492, 2.244],
             [ 4.042,  3.828,  3.618,  3.440, 3.059, 2.583],
             [ 5.964,  5.667,  5.109,  4.658, 3.949, 3.110],
             [10.961,  9.956,  8.454,  7.283, 5.703, 4.086],
             [15.937, 13.982, 11.432,  9.720, 7.121, 4.870],
             [21.155, 17.886, 14.380, 11.743, 8.605, 5.507],
             [26.083, 21.704, 17.040, 13.822, 9.660, 6.076]]

            ,

            # kappa=inf
            # alpha=0.2

            [[ 1.461,  1.398,  1.423,  1.416, 1.453, 1.445],
             [ 2.150,  2.062,  2.069,  2.063, 2.073, 2.034],
             [ 2.759,  2.709,  2.690,  2.657, 2.614, 2.520],
             [ 3.418,  3.362,  3.325,  3.233, 3.106, 2.883],
             [ 4.678,  4.600,  4.443,  4.259, 3.920, 3.425],
             [ 7.919,  7.605,  7.015,  6.490, 5.439, 4.316],
             [11.318, 10.467,  9.434,  8.355, 6.852, 4.941],
             [14.467, 13.403, 11.637, 10.058, 7.780, 5.577],
             [17.535, 16.176, 13.781, 11.663, 8.924, 6.106]]

            ],[

            # kappa=6.0
            # alpha=0.0

            [[ 0.960,  0.946,  0.952,  0.965,  0.968, 0.969],
             [ 2.027,  2.016,  1.982,  1.927,  1.829, 1.744],
             [ 3.147,  3.043,  2.910,  2.758,  2.492, 2.238],
             [ 4.189,  4.032,  3.745,  3.457,  3.061, 2.551],
             [ 6.381,  5.965,  5.348,  4.773,  3.974, 3.132],
             [11.858, 10.506,  8.798,  7.541,  5.863, 4.175],
             [17.293, 14.840, 11.990,  9.935,  7.297, 4.989],
             [22.569, 18.968, 14.941, 12.084,  8.837, 5.642],
             [28.061, 22.968, 17.595, 14.172, 10.146, 6.414]]

            ,

            # kappa=6.0
            # alpha=0.2

            [[ 1.757,  1.769,  1.787,  1.778, 1.821, 1.794],
             [ 2.372,  2.359,  2.369,  2.333, 2.330, 2.281],
             [ 2.990,  2.917,  2.918,  2.864, 2.770, 2.718],
             [ 3.613,  3.486,  3.435,  3.390, 3.261, 3.094],
             [ 4.687,  4.568,  4.479,  4.288, 4.051, 3.701],
             [ 7.659,  7.297,  6.750,  6.346, 5.647, 4.773],
             [10.475,  9.871,  8.996,  8.255, 6.940, 5.541],
             [13.511, 12.433, 11.031,  9.881, 8.139, 6.301],
             [16.399, 14.965, 12.972, 11.532, 9.218, 6.796]]

            ],[

            # kappa=4.0
            # alpha=0.0

            [[ 0.921,  0.922,  0.930,  0.934,  0.940, 0.940],
             [ 2.067,  2.050,  1.990,  1.951,  1.823, 1.729],
             [ 3.288,  3.127,  2.958,  2.783,  2.502, 2.189],
             [ 4.435,  4.165,  3.840,  3.509,  3.069, 2.556],
             [ 6.726,  6.163,  5.441,  4.877,  3.963, 3.149],
             [12.408, 10.891,  9.050,  7.675,  5.989, 4.199],
             [18.143, 15.408, 12.237, 10.206,  7.628, 5.098],
             [23.866, 19.660, 15.297, 12.484,  8.982, 5.769],
             [29.487, 23.781, 18.182, 14.588, 10.406, 6.515]]

            ,

            # kappa=4.0
            # alpha=0.2

            [[ 2.399,  2.322,  2.344,  2.328, 2.397, 2.352],
             [ 2.867,  2.818,  2.846,  2.825, 2.804, 2.766],
             [ 3.330,  3.329,  3.308,  3.273, 3.254, 3.157],
             [ 3.858,  3.816,  3.753,  3.704, 3.555, 3.497],
             [ 4.903,  4.787,  4.649,  4.535, 4.255, 4.054],
             [ 7.420,  7.115,  6.593,  6.352, 5.727, 5.146],
             [10.150,  9.375,  8.592,  8.009, 7.026, 6.038],
             [12.619, 11.533, 10.400,  9.499, 8.223, 6.782],
             [15.165, 13.686, 12.074, 10.986, 9.245, 7.427]]

            ]]

        Rs = np.array(Rs)
        etas = -np.array(etas)
        kappa_recips = np.array(kappa_recips)
        alphas = np.array(alphas)
        Is = np.array(Is)
        Is = np.transpose(Is, (0,1,3,2))

        table = {'axes': (kappa_recips, alphas, Rs, etas), 'values': Is}

    elif name == 'darian-marholm uncomplete cylinder':

        Rs = [1.0, 2.0, 3.0, 5.0, 10.0]
        etas = [0, 1, 2, 3, 5, 10, 15, 20, 25]
        kappa_recips = [0, 1.0/6, 1.0/4]
        alphas = [0, 0.2]
        Is = [[

            # kappa=inf
            # alpha=0.0

            [[0.974, 0.962, 0.956, 0.957, 0.940],
             [1.553, 1.543, 1.538, 1.545, 1.483],
             [1.945, 1.937, 1.920, 1.895, 1.788],
             [2.269, 2.257, 2.231, 2.176, 1.957],
             [2.800, 2.717, 2.729, 2.566, 2.226],
             [3.764, 3.682, 3.581, 3.254, 2.661],
             [4.562, 4.374, 4.207, 3.807, 2.922],
             [5.163, 4.978, 4.831, 4.251, 3.132],
             [5.688, 5.492, 5.207, 4.520, 3.325]]

            ,

            # kappa=inf
            # alpha=0.2

            [[1.496, 1.488, 1.482, 1.483, 1.471],
             [1.944, 1.931, 1.914, 1.923, 1.896],
             [2.263, 2.275, 2.290, 2.278, 2.248],
             [2.575, 2.551, 2.595, 2.555, 2.509],
             [3.044, 3.031, 3.060, 2.991, 2.812],
             [3.989, 3.987, 3.994, 3.765, 3.346],
             [4.750, 4.713, 4.734, 4.261, 3.585],
             [5.367, 5.346, 5.308, 4.759, 3.748],
             [5.937, 5.897, 5.684, 5.134, 4.083]]

            ],[

            # kappa=6.0
            # alpha=0.0

            [[0.946, 0.935, 0.934, 0.944, 0.931],
             [1.540, 1.540, 1.509, 1.541, 1.486],
             [1.937, 1.938, 1.889, 1.849, 1.748],
             [2.257, 2.241, 2.213, 2.145, 1.958],
             [2.790, 2.758, 2.723, 2.568, 2.131],
             [3.770, 3.691, 3.552, 3.239, 2.640],
             [4.460, 4.338, 4.098, 3.813, 2.932],
             [5.115, 4.975, 4.815, 4.276, 3.232],
             [5.692, 5.350, 5.213, 4.635, 3.525]]

            ,

            # kappa=6.0
            # alpha=0.2

            [[1.827, 1.804, 1.797, 1.825, 1.797],
             [2.209, 2.190, 2.165, 2.185, 2.154],
             [2.541, 2.504, 2.503, 2.503, 2.430],
             [2.791, 2.768, 2.751, 2.791, 2.730],
             [3.242, 3.251, 3.203, 3.208, 3.082],
             [4.153, 4.132, 4.058, 4.033, 3.682],
             [4.880, 4.844, 4.772, 4.540, 4.175],
             [5.483, 5.440, 5.332, 5.202, 4.327],
             [6.050, 5.961, 5.875, 5.523, 4.753]]

            ],[

            # kappa=4.0
            # alpha=0.0

            [[0.934, 0.927, 0.925, 0.920, 0.911],
             [1.538, 1.521, 1.510, 1.514, 1.433],
             [1.938, 1.910, 1.915, 1.865, 1.756],
             [2.248, 2.222, 2.164, 2.139, 1.958],
             [2.773, 2.755, 2.674, 2.584, 2.188],
             [3.736, 3.657, 3.551, 3.235, 2.742],
             [4.518, 4.348, 4.295, 3.788, 2.971],
             [5.145, 4.890, 4.709, 4.202, 3.326],
             [5.645, 5.356, 5.347, 4.630, 3.547]]

            ,

            # kappa=4.0
            # alpha=0.2

            [[2.219, 2.217, 2.198, 2.198, 2.159],
             [2.559, 2.530, 2.547, 2.551, 2.495],
             [2.849, 2.809, 2.810, 2.841, 2.800],
             [3.113, 3.059, 3.043, 3.074, 2.986],
             [3.567, 3.509, 3.465, 3.490, 3.351],
             [4.379, 4.347, 4.341, 4.264, 3.948],
             [5.082, 5.033, 4.976, 4.866, 4.454],
             [5.692, 5.614, 5.530, 5.359, 4.894],
             [6.219, 6.185, 6.006, 5.848, 5.174]]

            ]]

        Rs = np.array(Rs)
        etas = -np.array(etas)
        kappa_recips = np.array(kappa_recips)
        alphas = np.array(alphas)
        Is = np.array(Is)
        Is = np.transpose(Is, (0,1,3,2))

        table = {'axes': (kappa_recips, alphas, Rs, etas), 'values': Is}

    elif name == 'darian-marholm cylinder':

        table = get_table('darian-marholm uncomplete cylinder', False)
        kappa_recips, alphas, Rs, etas = table['axes']

        # Add R=0 column
        Rs = np.concatenate([[0], Rs])
        vals = np.zeros((len(kappa_recips), len(alphas), len(Rs), len(etas)))
        vals[:,:,1:,:] = table['values']

        # Insert analytical results for R=0 and eta=0
        for i, kappa_recip in enumerate(kappa_recips):
            kappa = 1/kappa_recip if kappa_recip != 0 else float('inf')
            for j, alpha in enumerate(alphas):
                sp = Species(n=1e11, T=1000, kappa=kappa, alpha=alpha)
                geo = Cylinder(1.0, 1)
                vals[i,j,0,:] = OML_current(geo, sp, eta=etas, normalization='thmax')
                for k, R in enumerate(Rs):
                    geo = Cylinder(R*sp.debye, 1)
                    vals[i,j,k,0] = thermal_current(geo, sp, normalization='thmax')

        table['values'] = vals
        table['axes'] = (kappa_recips, alphas, Rs, etas)

    elif name == 'darian-marholm sphere':

        table = get_table('darian-marholm uncomplete sphere', False)
        kappa_recips, alphas, Rs, etas = table['axes']

        # Add R=0 column
        Rs = np.concatenate([[0], Rs])
        vals = np.zeros((len(kappa_recips), len(alphas), len(Rs), len(etas)))
        vals[:,:,1:,:] = table['values']

        # Insert analytical results for R=0 and eta=0
        for i, kappa_recip in enumerate(kappa_recips):
            kappa = 1/kappa_recip if kappa_recip != 0 else float('inf')
            for j, alpha in enumerate(alphas):
                sp = Species(n=1e11, T=1000, kappa=kappa, alpha=alpha)
                geo = Sphere(1.0)
                vals[i,j,0,:] = OML_current(geo, sp, eta=etas, normalization='thmax')
                for k, R in enumerate(Rs):
                    geo = Sphere(R*sp.debye)
                    vals[i,j,k,0] = thermal_current(geo, sp, normalization='thmax')

        table['values'] = vals
        table['axes'] = (kappa_recips, alphas, Rs, etas)

    elif name == 'laframboise-darian-marholm sphere':

        # Get darian-marholm table
        table = get_table('darian-marholm sphere')

        # Grid will no longer be regular.
        table['values'] = table['values'].reshape(-1)
        del table['axes']

        # Remove points with kappa=inf and alpha=0
        keep = map(lambda x: abs(x[0])>tol or abs(x[1])>tol, table['points'])
        keep = np.where(list(keep))
        table['points'] = table['points'][keep]
        table['values'] = table['values'][keep]

        # Add all values from laframboise 
        t = get_table('laframboise sphere')
        points = np.zeros((len(t['points']),4))
        points[:,2:] = t['points']
        values = t['values'].reshape(-1)
        table['points'] = np.concatenate((table['points'], points))
        table['values'] = np.concatenate((table['values'], values))

    elif name == 'laframboise-darian-marholm cylinder':

        # Get darian-marholm table
        table = get_table('darian-marholm cylinder')

        # Grid will no longer be regular.
        table['values'] = table['values'].reshape(-1)
        del table['axes']

        # Remove points with kappa=inf and alpha=0
        keep = map(lambda x: abs(x[0])>tol or abs(x[1])>tol, table['points'])
        keep = np.where(list(keep))
        table['points'] = table['points'][keep]
        table['values'] = table['values'][keep]

        # Add all values from laframboise 
        t = get_table('laframboise cylinder')
        points = np.zeros((len(t['points']),4))
        points[:,2:] = t['points']
        values = t['values'].reshape(-1)
        table['points'] = np.concatenate((table['points'], points))
        table['values'] = np.concatenate((table['values'], values))

    else:

        raise ValueError('Table {} does not exist'.format(name))

    if provide_points and 'points' not in table:
        table['points'] = cartesian_product(*table['axes'])

    return table
