
import numpy as np
from astropy.coordinates import Distance


prsc_lindegren = {
    # Weighted average estimate using Bailer-Jones distances for the
    # selected members of this cluster: 4000 +- 1300 pc
    'vdBH85': (2.769, 4.153, 5.536),
    'vdBH73': (4.475, 4.919, 5.287),
    'vdBH92': (2.331, 2.434, 2.514),
    'RUP85': (4.447, 4.640, 4.831),
    'vdBH106': (4.372, 4.771, 5.157),
    'TR13': (4.439, 4.584, 4.711),
    'NGC4349': (1.980, 1.996, 2.035),
    'RUP162': (4.185, 4.365, 4.553),
    'TR12': (3.505, 3.631, 3.758),
    'vdBH87': (2.200, 2.261, 2.323),
}

prsc_shuangjing = {
    # 'vdBH85': (5.341, 5.645, 5.954),
    'vdBH73': (3.708, 4.054, 4.372),
    'vdBH92': (2.099, 2.173, 2.248),
    'RUP85': (3.695, 3.830, 3.972),
    'vdBH106': (3.740, 4.057, 4.345),
    'TR13': (3.660, 3.751, 3.847),
    'NGC4349': (1.807, 1.828, 1.863),
    'RUP162': (3.542, 3.662, 3.806),
    'TR12': (3.017, 3.109, 3.207),
    'vdBH87': (1.996, 2.048, 2.095),
}

prsc_schonrich = {
    # 'vdBH85': (5.982, 6.283, 6.619),
    'vdBH73': (4.147, 4.460, 4.769),
    'vdBH92': (2.216, 2.282, 2.351),
    'RUP85': (4.013, 4.164, 4.319),
    'vdBH106': (3.995, 4.311, 4.663),
    'TR13': (3.983, 4.100, 4.207),
    'NGC4349': (1.863, 1.897, 1.945),
    'RUP162': (3.798, 3.936, 4.101),
    'TR12': (3.223, 3.314, 3.424),
    'vdBH87': (2.077, 2.133, 2.185),
}

prsc_nooffset = {
    # 'vdBH85': (7.798, 8.228, 8.660),
    'vdBH73': (5.041, 5.482, 5.931),
    'vdBH92': (2.499, 2.606, 2.722),
    'RUP85': (5.063, 5.298, 5.539),
    'vdBH106': (5.034, 5.407, 5.804),
    'TR13': (5.100, 5.250, 5.421),
    'NGC4349': (2.105, 2.115, 2.154),
    'RUP162': (4.788, 4.973, 5.179),
    'TR12': (3.941, 4.084, 4.221),
    'vdBH87': (2.356, 2.428, 2.495),
}

dm_asteca = {
    'vdBH85': (13.317, 0.12085, 9.872),
    'vdBH106': (13.437, 0.36267, 9.482),
    'TR12': (12.721, 0.08997, 8.822),
    'RUP85': (13.556, 0.14048, 8.444),
    'vdBH73': (13.498, 0.26366, 8.893),
    'vdBH87': (11.437, 0.06779, 8.307),
    'TR13': (13.411, 0.15, 8.044),
    'NGC4349': (11.708, 0.090487, 8.708),
    'vdBH92': (12.066, 0.092136, 7.276),
    'RUP162': (13.23, 0.10038, 8.924),
}


for cl, dm in dm_asteca.items():
    d_ast = 10 ** ((dm[0] + 5.) / 5.)
    e_ast = round((dm[1] * d_ast * .2 * np.log(10)) / 1000., 2)
    d_ast = round(d_ast / 1000., 2)

    d_lin = round(prsc_lindegren[cl][1], 2)
    e_lin = round(.5 * (prsc_lindegren[cl][2] - prsc_lindegren[cl][0]), 2)
    if cl == 'vdBH85':
        d_noo, e_noo, d_sch, e_sch, d_xu, e_xu = [np.nan] * 6
    else:
        d_noo = round(prsc_nooffset[cl][1], 2)
        e_noo = round(.5 * (prsc_nooffset[cl][2] - prsc_nooffset[cl][0]), 2)
        d_sch = round(prsc_schonrich[cl][1], 2)
        e_sch = round(.5 * (prsc_schonrich[cl][2] - prsc_schonrich[cl][0]), 2)
        d_xu = round(prsc_shuangjing[cl][1], 2)
        e_xu = round(.5 * (prsc_shuangjing[cl][2] - prsc_shuangjing[cl][0]), 2)

    print(r"{}  ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$\\".format(
        cl, d_ast, e_ast, d_noo, e_noo, d_lin, e_lin, d_sch, e_sch, d_xu, e_xu))
