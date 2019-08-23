
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import matplotlib.cm as cm


def main(
    dm_asteca, prsc_nooffset, prsc_lindegren, prsc_schonrich,
        prsc_shuangjing):
    """
    """
    parsec_data, plx_data = {}, {}
    for cl, asteca_data in dm_asteca.items():

        # ASteCA distance modulus to parsecs
        dm, e_dm, log_age = asteca_data[:3]
        plx_asteca = 1. / (10 ** ((dm + 5.) / 5.))
        pc50_asteca = 1. / plx_asteca
        e_pc_asteca = .2 * np.log(10.) * pc50_asteca * e_dm
        # e_plx_asteca = 1. / e_pc_asteca**2
        pc16_asteca = pc50_asteca - e_pc_asteca
        pc84_asteca = pc50_asteca + e_pc_asteca

        # No offset parallax distances
        plx16_nooffset, plx50_nooffset, plx84_nooffset =\
            1. / (1000. * np.array(prsc_nooffset[cl]))
        # e_plx_nooffset = .5 * (plx16_nooffset + plx84_nooffset)

        # Order: ASteCA, no offset
        plx_data[cl] = np.array([plx_asteca, plx50_nooffset, log_age])

        # Lindegren parsec distances
        pc16_lindegren, pc50_lindegren, pc84_lindegren = 1000. * np.array(
            prsc_lindegren[cl])

        # Schonrich parsec distances
        pc16_schonrich, pc50_schonrich, pc84_schonrich = 1000. * np.array(
            prsc_schonrich[cl])

        # Shuangjing parsec distances
        pc16_shuangjing, pc50_shuangjing, pc84_shuangjing = 1000. * np.array(
            prsc_shuangjing[cl])

        # No offset parsec distances
        pc16_nooffset, pc50_nooffset, pc84_nooffset = 1000. * np.array(
            prsc_nooffset[cl])

        # Order: ASteCA, Lindegren, Shuangjing, no offset
        parsec_data[cl] = np.array([
            [pc16_asteca, pc50_asteca, pc84_asteca, log_age],
            [pc16_lindegren, pc50_lindegren, pc84_lindegren],
            [pc16_schonrich, pc50_schonrich, pc84_schonrich],
            [pc16_shuangjing, pc50_shuangjing, pc84_shuangjing],
            [pc16_nooffset, pc50_nooffset, pc84_nooffset]])

    # Make final plot
    print("Plotting")
    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(10, 10)
    # plt.style.use('seaborn-darkgrid')
    xymin, xymax = 1050., 6690.

    gsi = ((0, 2, 2, 4), (0, 2, 4, 6), (0, 2, 6, 8))
    txt = (
        'Lindegren et al. (0.029 mas)', 'Schönrich et al. (0.054 mas)',
        'Xu et al. (0.075 mas)')
    labels = ((True, False), (False, False), (False, True))
    for j in range(1, 4):
        dpcPlot(
            gs, gsi[j - 1], xymin, xymax, txt[j - 1], parsec_data,
            labels[j - 1], j)

    fig.tight_layout()
    plt.savefig('dist_comparision1.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close("all")

    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(10, 10)
    # plt.style.use('seaborn-darkgrid')

    dpcPlot(
        gs, (0, 2, 0, 2), xymin, xymax, 'No bias correction', parsec_data,
        (True, True), 4)

    ax = plt.subplot(gs[2:4, 0:2])
    # ax.set_title(r"No bias correction", fontsize=10)
    for cl, plxs_age in plx_data.items():
        # plxs = (ASteCA, Bayes)
        x, y, _ = plxs_age * 1000.
        if cl == 'vdBH87':
            xo = -.035
        elif cl == 'NGC4349':
            xo = -.045
        else:
            xo = 0.
        cl = cl + '*' if cl in ('TR12', 'RUP162', 'vdBH106') else cl
        ax.annotate(cl, (y + xo, x - y + .005), fontsize=9)
    x, y, lage = np.array(list(plx_data.values())).T
    x, y = x * 1000., y * 1000.
    deltas = x - y
    plt.scatter(y, x - y, zorder=4, c=lage, marker='s', s=50)
    plt.axhline(np.mean(deltas), c='r', ls='--', lw=.7, zorder=1,
                label=r"$\Delta_{{mean}}$   = {:.3f} [mas]".format(
                    np.mean(deltas)))
    plt.axhline(np.median(deltas), c='blue', ls='--', lw=.7, zorder=1,
                label=r"$\Delta_{{median}}$ = {:.3f} [mas]".format(
                    np.median(deltas)))
    plt.grid(ls=':', lw=.5, c='grey')
    # plt.xlabel(r'$Plx_{ASteCA}$ [mas]', fontsize=12)
    plt.xlabel(r'$Plx_{Bayes}$ [mas]', fontsize=12)
    plt.ylabel(r'$\Delta\;Plx_{[ASteCA-Bayes]}$ [mas]', fontsize=12)
    plt.legend(loc=0, fontsize=12)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=9)

    fig.tight_layout()
    plt.savefig('dist_comparision2.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close("all")


def dpcPlot(gs, gsi, xymin, xymax, txt, parsec_data, labels, j):
    ax = plt.subplot(gs[gsi[0]:gsi[1], gsi[2]:gsi[3]])
    deltas, xyall = [], []
    for cl, dists in parsec_data.items():
        exl, x, exh, lage = dists[0]
        exl, exh = x - exl, exh - x
        eyl, y, eyh = dists[j]
        eyl, eyh = y - eyl, eyh - y
        deltas.append(x - y)
        xyall.append([x, y, lage])
        ax.errorbar(
            x, y, xerr=np.array([[exl, exh]]).T,
            yerr=np.array([[eyl, eyh]]).T, fmt='.', elinewidth=.85, ms=4,
            ecolor='b', label=None)
        if cl == 'vdBH87':
            xo, yo = -300., -250.
        elif cl == 'vdBH106':
            xo, yo = -900., -200.
        else:
            xo, yo = (50., 50.)
        cl = cl + '*' if cl in ('TR12', 'RUP162', 'vdBH106') else cl
        ax.annotate(cl, (x + xo, y + yo), fontsize=9)
    x, y, lage = np.array(xyall).T
    plt.scatter(x, y, marker='s', s=20, c=lage, zorder=6)  # , cmap='seismic')
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=9)
    ax.set_title(
        r"{}, $\Delta_{{mean}}=${:.0f} [pc]".format(
            txt, np.mean(deltas)), fontsize=10)  # , x=.1, y=.95)
    plt.plot((xymin, xymax), (xymin, xymax), ls='--', c='r', lw=.5)
    plt.grid(ls=':', lw=.5, c='grey')
    plt.xlabel(r'$d_{ASteCA}$ [pc]', fontsize=12)
    if labels[0]:
        plt.ylabel(r'$d_{Bayes}$ [pc]', fontsize=12)
    else:
        ax.set_yticklabels([])
    if labels[1]:
        cbar.set_label(r"$\log(age)$", fontsize=10, labelpad=4)
    plt.xlim(xymin, xymax)
    plt.ylim(xymin, xymax)


if __name__ == '__main__':
    prsc_lindegren = {
        'TR12': (3.446, 3.567, 3.710),
        'RUP162': (3.608, 3.779, 3.952),
        'vdBH92': (2.213, 2.322, 2.421),
        'vdBH106': (4.150, 4.505, 4.902),
        'vdBH87': (2.151, 2.213, 2.286),
        'TR13': (4.240, 4.385, 4.515),
        'vdBH73': (4.394, 4.719, 5.081),
        'NGC4349': (2.065, 2.137, 2.241),
        'RUP85': (4.446, 4.642, 4.815),
        'vdBH85': (5.288, 5.593, 5.919),
        #
        'LYNGA15': (2.763, 2.924, 3.096),
        'LODEN565': (2.477, 2.680, 2.864),
        'RUP87': (1.515, 1.764, 2.681),
        'vdBH91': (2.092, 2.266, 2.435),
        'NGC4230': (2.593, 2.866, 3.144),
        'RUP88': (3.854, 4.349, 4.826)
    }

    prsc_shuangjing = {
        'TR12': (2.960, 3.056, 3.148),
        'RUP162': (3.053, 3.201, 3.350),
        'vdBH92': (1.997, 2.097, 2.189),
        'vdBH106': (3.451, 3.793, 4.153),
        'vdBH87': (1.949, 2.006, 2.065),
        'TR13': (3.565, 3.672, 3.792),
        'vdBH73': (3.673, 3.976, 4.312),
        'NGC4349': (1.871, 1.930, 2.030),
        'RUP85': (3.705, 3.855, 4.003),
        'vdBH85': (4.453, 4.728, 4.999),
        #
        'NGC4230': (2.271, 2.551, 2.879),
        'LYNGA15': (2.424, 2.545, 2.683),
        'LODEN565': (2.170, 2.382, 2.526),
        'RUP88': (3.337, 3.764, 4.248),
        'RUP87': (1.423, 1.746, 2.744),
        'vdBH91': (1.896, 2.043, 2.196)
    }

    prsc_schonrich = {
        'TR12': (3.159, 3.264, 3.376),
        'RUP162': (3.258, 3.408, 3.580),
        'vdBH92': (2.095, 2.202, 2.284),
        'vdBH106': (3.865, 4.212, 4.545),
        'vdBH87': (2.036, 2.097, 2.163),
        'TR13': (3.834, 3.950, 4.070),
        'vdBH73': (3.994, 4.317, 4.620),
        'NGC4349': (1.953, 2.012, 2.094),
        'RUP85': (3.996, 4.163, 4.330),
        'vdBH85': (4.812, 5.099, 5.345),
        #
        'NGC4230': (2.327, 2.625, 2.913),
        'LYNGA15': (2.563, 2.702, 2.863),
        'LODEN565': (2.298, 2.500, 2.652),
        'RUP88': (3.550, 3.998, 4.463),
        'RUP87': (1.452, 1.773, 2.878),
        'vdBH91': (1.970, 2.121, 2.312)
    }

    prsc_nooffset = {
        'TR12': (3.838, 3.972, 4.120),
        'RUP162': (3.961, 4.184, 4.403),
        'vdBH92': (2.353, 2.479, 2.618),
        'vdBH106': (4.633, 5.060, 5.547),
        'vdBH87': (2.313, 2.374, 2.444),
        'TR13': (4.837, 5.019, 5.193),
        'vdBH73': (4.895, 5.317, 5.759),
        'NGC4349': (2.208, 2.294, 2.450),
        'RUP85': (5.048, 5.303, 5.544),
        'vdBH85': (5.792, 6.149, 6.525),
        #
        'LYNGA15': (3.014, 3.183, 3.365),
        'LODEN565': (2.625, 2.873, 3.096),
        'RUP88': (4.277, 4.847, 5.437),
        'RUP87': (2.809, 3.741, 4.418),
        'vdBH91': (2.202, 2.404, 2.627),
        'NGC4230': (2.739, 3.093, 3.459)
    }

    dm_asteca = {
        'vdBH87': (11.272, 0.059715, 8.505, 0.2164, 0.58, 0.0112, 954, 252.02),
        'vdBH106': (13.997, 0.21338, 9.987, 0.03759, 0.296, 0.0254, 1845, 349.8),
        'vdBH92': (11.462, 0.12332, 7.92, 0.2157, 0.655, 0.033, 399, 74.563),
        'RUP162': (12.251, 0.12524, 7.485, 0.1558, 0.474, 0.0242, 523, 92.22),
        'TR12': (12.409, 0.11919, 8.903, 0.09872, 0.288, 0.0275, 718, 97.),
        'vdBH85': (13.09, 0.12321, 9.869, 0.05458, 0.284, 0.0273, 2767, 270.65),
        'RUP85': (13.479, 0.17536, 8.615, 0.06377, 1.06, 0.0417, 2727, 649.47),
        'vdBH73': (12.106, 0.29436, 9.783, 0.08208, 0.623, 0.0368, 978, 805.69),
        'NGC4349': (12.194, 0.13254, 8.598, 0.1348, 0.401, 0.0251, 1106, 239.63),
        'TR13': (13.003, 0.124, 7.427, 0.177, 0.554, 0.014, 505, 148),
        #
        # 'LYNGA15': (11.035, 0.081641, 9.938, 0.06299, 0.271, 0.0228, 491, 92.),
        # 'RUP88': (14.746, 0.37555, 9.541, 0.1185, 0.331, 0.056, 1047, 293.66),
        # 'NGC4230': (15.176, 0.25646, 9.631, 0.07867, 0.21, 0.039, 1360, 485.47),
        # 'LODEN565': (12.595, 0.13702, 8.1, 0.2375, 0.891, 0.0149, 273, 69.569),
        # 'RUP87': (13.522, 0.31632, 9.142, 0.06687, 0.228, 0.0387, 191, 57.954),
        # 'vdBH91': (13.399, 0.245, 8.631, 0.2927, 0.746, 0.0904, 473, 126.66)
    }

    main(
        dm_asteca, prsc_nooffset, prsc_lindegren, prsc_schonrich,
        prsc_shuangjing)
