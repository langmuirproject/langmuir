import numpy as np
from scipy.interpolate import griddata

def cartesian_product(*arrays):
    return np.dstack(np.meshgrid(*arrays)).reshape(-1, len(arrays))

def get_table(name, provide_points=True):

    name = name.lower()

    if name == 'laframboise sphere':

        # Table 5c in Laframboise's thesis

        Rs =  [0, 0.2, 0.3, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 50, 100]

        Vs =  [0.0  , 0.1  , 0.3  , 0.6  , 1.0  , 1.5  , 2.0  , 3.0  , 5.0  , 7.5  , 10.0  , 15.0  , 20.0  , 25.0]
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
        Vs = np.array(Vs)
        Is = np.array(Is).T.reshape(-1)

        table = {'axes': (Rs, Vs), 'values': Is}

    elif name == 'laframboise cylinder':

        # Table 6c in Laframboise's thesis

        Rs =  [0, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 20, 30, 40, 50, 100]

        Vs =  [0.0  , 0.1   , 0.3   , 0.6   , 1.0   , 1.5   , 2.0   , 3.0   , 5.0   , 7.5   , 10.0  , 15.0  , 20.0  , 25.0]
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
        Vs = np.array(Vs)
        Is = np.array(Is).T.reshape(-1)

        table = {'axes': (Rs, Vs), 'values': Is}

    elif name == 'darian-marholm sphere':

        raise NotImplementedError
        # Results from paper. Table with four axes: R, V (eta), kappa, alpha.

    elif name == 'darian-marholm cylinder':

        raise NotImplementedError
        # Results from paper. Table with four axes: R, V (eta), kappa, alpha.

    elif name == 'laframboise-darian-marholm sphere':

        raise NotImplementedError
        # Same as darian-marholm but with the values for the pure Maxwellian
        # distribution removed and replaced by those of Laframboise, for the
        # best accuracy. Since his table contains more data points the values
        # will not be on a uniform grid in 4D space, and the dictionary can
        # not have the 'axes' entry. Provide the 'points' entry instead.
        # Do not duplicate values, but obtain them from the other tables.

    elif name == 'laframboise-darian-marholm cylinder':

        raise NotImplementedError
        # Same as darian-marholm but with the values for the pure Maxwellian
        # distribution removed and replaced by those of Laframboise, for the
        # best accuracy. Since his table contains more data points the values
        # will not be on a uniform grid in 4D space, and the dictionary can
        # not have the 'axes' entry. Provide the 'points' entry instead.
        # Do not duplicate values, but obtain them from the other tables.

    else:

        raise ValueError('Table {} does not exist'.format(name))

    if provide_points and 'points' not in table:
        table['points'] = cartesian_product(*table['axes'])

    return table
