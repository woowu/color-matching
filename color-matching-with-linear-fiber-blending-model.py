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
r_object = [
    0.507045, 0.532145, 0.546225, 0.557051,
    0.493464, 0.430796, 0.370752, 0.359875,
    0.365956, 0.513965, 0.716593, 0.837977,
    0.888602, 0.900323, 0.906528, 0.913298,
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

# The RGB stimulation constans over light spectrum
# lambda: 400-700 nm, Dlamda = 20nm
#
rgb_spectrum = [];

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
        [1.7860, 0.8317],
        [1.9043, 1.4825],
        [1.9883, 1.6976],
        [2.0520, 1.6397],
        [1.1888, 1.7357],
        [0.7632, 1.7869],
        [0.5156, 1.8259],
        [0.4769, 1.8516],
        [0.4955, 1.8232],
        [1.1861, 1.9543],
        [5.4442, 1.9518],
        [25.6587, 1.7787], 
        [131.2240, 1.1279], 
        [257.6183, 0.6877], 
        [757.9442, -3.0313], 
        [1794.6916, -9.2860], 
    ],
    [
        [0.7993, 1.6222],
        [0.7520, 1.6376],
        [0.6724, 1.6838],
        [0.7476, 1.6964],
        [0.9851, 1.6904],
        [1.5126, 1.6907],
        [3.2541, 1.6445],
        [8.5218, 1.6054],
        [21.4947, 1.4948],
        [51.5542, 1.5561],
        [184.0496, 1.1489],
        [574.1510, 0.6345],
        [1278.8192, -0.5633],
        [2139.7924, -1.7802],
        [3224.9009, -3.4094],
        [3927.2819, -3.6809],
    ],
    [
        [3.4242, 1.8489],
        [4.3507, 1.9243],
        [4.5049, 1.9297],
        [3.5100, 1.9541],
        [3.2415, 1.9519],
        [2.3654, 1.9244],
        [1.8767, 1.9338],
        [1.4432, 1.9656],
        [0.9926, 1.9779],
        [0.7591, 1.9691],
        [0.6781, 1.9687],
        [0.6452, 1.9678],
        [0.6710, 1.9721],
        [0.9926, 1.9846],
        [2.5829, 1.9653],
        [3.7151, 1.9928],
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
            #print(l, c_index, ':', ks_fun, c[c_index], ks)
            ks_sum += ks
        ks_sum += ks_base[l]
        ks_vector.append(ks_sum)
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
        print('solve for r when ks = ' + str(e) + ': (a) ' + str(1 + e - d))
        print('solve for r when ks = ' + str(e) + ': (b) ' + str(1 + e + d))
        if d <= 1 + e:      # i guess the smaller root should be real reflectance
            r.append(1 + e - d)
        else:
            r.append(1 + e + d)
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
            500 * math.pow(x/x0, 1/3) - math.pow(y/y0, 1/3),
            200 * math.pow(y/y0, 1/3) - math.pow(z/z0, 1/3),
            ]
    else:
        return [
            093.3 * math.exp(y/y0, 1/3),
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

def obj(c, params=None):
    lab = xyz_to_lab(r_to_xyz(ks_to_r(c_to_ks(c))))
    lab_object = xyz_to_lab(r_to_xyz(r_object))
    return delta_e(lab, lab_object)

def main():
    print(obj([.05e-2, .05e-2, .05e-2]))

if __name__ == '__main__':
    main()

#option = minimize(obj, method='COBYLA', x0=[0.001, 0.001, 0.001], bounds=(0, 0.05))
#print(option)


