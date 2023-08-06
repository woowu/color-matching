#!/usr/bin/python3
from scipy.optimize import *
import math
import random

# number of lambda
#
l_num = len(range(400, 720, 20))

# number of fibers
#
c_num = len(range(3))

# r object
r_objects = [
    [
        0.507045, 0.532145, 0.546225, 0.557051,
        0.493464, 0.430796, 0.370752, 0.359875,
        0.365956, 0.513965, 0.716593, 0.837977,
        0.888602, 0.900323, 0.906528, 0.913298,
    ],
    [
        0.681044, 0.724153, 0.749374, 0.752377,
        0.733554, 0.699726, 0.663874, 0.652066,
        0.642280, 0.689994, 0.720488, 0.725711,
        0.732905, 0.771153, 0.838096, 0.859082,
    ],
    [
        0.629420, 0.668642, 0.690550, 0.693808,
        0.654061, 0.603730, 0.553318, 0.539833,
        0.534764, 0.616401, 0.680087, 0.695650,
        0.705559, 0.747784, 0.823197, 0.846680,
    ],
    [
        0.598105, 0.636542, 0.647724, 0.634758,
        0.638956, 0.610761, 0.589352, 0.556888,
        0.500971, 0.457357,	0.439055, 0.430547,
        0.437718, 0.503085, 0.646008, 0.692619, 
    ],
    [
        0.641290, 0.680826, 0.695806, 0.685813,
        0.690465, 0.665507, 0.646113, 0.616788,
        0.565565, 0.524583, 0.507192, 0.499046,
        0.506062, 0.568136, 0.697296, 0.738285,
    ],
    [
        0.507637, 0.534544, 0.550889, 0.560809, 0.488216, 0.420199, 0.357506, 0.345772, 0.351452, 0.499036, 0.703689, 0.825218, 0.874593, 0.889770, 0.901908, 0.910100,
    ],
    [
        0.376815, 0.379640, 0.371046, 0.387641, 0.401820, 0.412266, 0.412390, 0.425963, 0.439564, 0.576515, 0.736379, 0.810591, 0.836221, 0.858021, 0.886589, 0.898088,
    ],
    [
        0.611623, 0.636383, 0.644250, 0.654583, 0.661286, 0.657130, 0.645990, 0.645381, 0.639611, 0.688578, 0.719943, 0.725500, 0.732772, 0.771024, 0.837919, 0.858908,
    ],
    [
        0.616192, 0.654238, 0.667361, 0.660191, 0.663057, 0.638999, 0.618215, 0.593015, 0.547792, 0.519088, 0.506045, 0.498777, 0.505973, 0.568058, 0.697160, 0.738137,
    ],
    [
        0.457151, 0.471454, 0.465966, 0.472114, 0.495952, 0.506516, 0.517486, 0.510463, 0.471248, 0.447707, 0.436959, 0.430048, 0.437543, 0.502925, 0.645742, 0.692336,
    ],
]


# RGB stimulation matrix:
#

rgb_stimulation = [
    [0.136, 0.014, 0.613],
    [1.644, 0.172, 7.820],
    [3.463, 0.560, 17.755],
    [3.065, 1.300, 17.697],
    [0.803, 2.530, 7.703],
    [0.036, 4.337, 2.056],
    [1.062, 6.870, 0.548],
    [3.385, 8.644, 0.123],
    [6.069, 8.583, 0.000],
    [8.361, 7.163, 0.000],
    [8.707, 5.100, 0.000],
    [6.463, 3.004, 0.000],
    [3.109, 1.295, 0.000],
    [1.053, 0.416, 0.000],
    [0.275, 0.107, 0.000],
    [0.059, 0.023, 0.000],
]

# ks for base fiber, each element for a lambda
#
ks_base = [
    0.039492, 0.025906, 0.017964, 0.015092,
    0.011439, 0.009515, 0.007961, 0.006947,
    0.006284, 0.005889, 0.005238, 0.004948,
    0.004626, 0.004247, 0.004100, 0.003617
];

# from a single r to a k/s value, it's a function ks = m*r + a,
# but there a three type of r for R,G,B respectively.
#
ks_fun_coef = [
    [
        [71.3592, 0.6580],
        [66.6862, 0.2800],
        [63.8362, 0.1607],
        [62.0659, 0.1797],
        [107.1160, 0.2298],
        [166.7703, 0.2925],
        [246.7703, 0.3590],
        [266.7703, 0.3347],
        [256.7703, 0.3793],
        [107.3490, 0.0464],
        [23.3392, 0.0115],
        [4.9626, 0.0090],
        [0.9681, 0.0068],
        [0.4930, 0.0052],
        [0.1675, 0.0067],
        [0.0675, 0.0064],
    ],
    [
        [159.1008, 0.4882],
        [169.1008, 0.4985],
        [189.1008, 0.4889],
        [170.1008, 0.4227],
        [129.1008, 0.3267],
        [84.0990, 0.2121],
        [39.0990, 0.1126],
        [14.9213, 0.0478],
        [5.9179, 0.0240],
        [2.4191, 0.0098],
        [0.6910, 0.0047],
        [0.2210, 0.0024],
        [0.0994, 0.0020],
        [0.0594, 0.0018],
        [0.0394, 0.0017],
        [0.0324, 0.0014],
    ],
    [
        [37.0638, 0.0492],
        [29.2607, 0.0196],
        [28.2607, 0.0177],
        [36.2607, 0.0160],
        [39.2607, 0.0181],
        [53.8878, 0.0347],
        [67.8878, 0.0393],
        [88.0412, 0.0337],
        [128.0412, 0.0361],
        [167.5274, 0.0568],
        [187.5274, 0.0641],
        [197.1128, 0.0688],
        [189.5274, 0.0598],
        [128.0412, 0.0293],
        [49.2607, 0.0177],
        [34.2607, 0.0046],
    ],
]

obj_r = [ random.random() for _ in range(l_num) ];

# c is a 3-d fiber vector
# Returns a vector of ks, each for a lambda
#
def c_to_ks(c):
    ks_vector = []
    for l in range(l_num):
        ks_sum = 0;
        for c_index in range(len(c)):
            ks_fun = ks_fun_coef[c_index][l]
            ks = ks_fun[0] * c[c_index] + ks_fun[1]
            #print(' check ks value for:', l, c_index, ':', ks_fun, c[c_index], ks)
            ks_sum += ks * c[c_index]
        ks_sum += ks_base[l]
        ks_vector.append(ks_sum)
    #print('ks vector:', ks_vector)
    return ks_vector

# for a given ks vector, calculate the corresponding r
# vector according the one-constant K-M model.
#
def ks_to_r(ks):
    r = []
    for e in ks: 
        # the quadratic formular is:
        #   ks = (1 - r)^2 / 2r
        # hence solving it for r is:
        #   r = (1 + ks) +- sqrt( ks(2 + ks) )
        # i like to choose the smaller r as the root since it's looks reasonable.
        #
        d = e*(2 + e)        # d must >= 0 because ks is non-negative
        d = math.sqrt(d)
        #print('solve for r when ks = ' + str(e) + ': (a) ' + str(1 + e - d))
        r.append(1 + e - d)
    return r

def c_to_r(c):
    return ks_to_r(c_to_ks(c))

# r is a reflectance vector, each for a lambda.
#
def r_to_xyz(r):
    dl = 20
    k = .1

    xyz = []
    for color_index in range(3):
        sum = 0
        for l in range(l_num):  # 16 lambda's
            stimulation = rgb_stimulation[l][color_index]
            sum += stimulation * r[l] * dl
        xyz.append(sum * k)
    #print('xyz:', xyz)
    return xyz

def xyz_to_lab(xyz):
    idea_white = [94.83, 100, 107.38]
    option = 1
    index = 0
    for e in xyz:
        if e / idea_white[index] <= 0.008856:
            option = 2
            break
        ++index

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    x0 = idea_white[0]
    y0 = idea_white[1]
    z0 = idea_white[2]
    if option == 1:
        return [
            116 * math.pow(y/y0, 1/3) - 16,
            500 * (math.pow(x/x0, 1/3) - math.pow(y/y0, 1/3)),
            200 * (math.pow(y/y0, 1/3) - math.pow(z/z0, 1/3)),
            ]
    else:
        return [
            093.3 * y/y0,
            3893.5 * (x/x0 - y/y0),
            1557.4 * (y/y0 - z/z0),
            ]

def delta_e(lab1, lab2):
    sum = 0;
    for i in range(len(lab1)):
        e1 = lab1[i]
        e2 = lab2[i]
        sum += (e1 - e2)**2
    return math.sqrt(sum)

def main():
    for i in range(10):
        r_object = r_objects[i]

        def obj(c, params=None):
            lab = xyz_to_lab(r_to_xyz(ks_to_r(c_to_ks(c))))
            lab_object = xyz_to_lab(r_to_xyz(r_object))
            #print('lab_object:', lab_object)
            #print('lab_var:', lab)
            return delta_e(lab, lab_object)

        result = minimize(obj, method='Nelder-Mead', x0=[0.1e-2, 0.1e-2, 0.11e-2], bounds=((0, 5e-2*1.2), (0, 5e-2*1.2), (0, 5e-2*1.2)), tol=1e-10)
        print(result)

if __name__ == '__main__':
    main()

