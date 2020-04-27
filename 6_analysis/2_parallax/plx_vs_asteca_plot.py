
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import matplotlib.cm as cm


def main(
    dm_asteca, prsc_nooffset, prsc_lindegren, prsc_schonrich,
        prsc_shuangjing, offset_flag=None, all_flag=False):
    """
    """
    parsec_data, plx_data = {}, {}
    for cl, asteca_data in dm_asteca.items():
        print(cl)

        # ASteCA distance modulus to parsecs
        dm, e_dm, log_age = asteca_data

        # Possible BV correction: BV_G = BV + 0.0153
        # Ag = cg * Av = cg * 3.1 * E_BV
        # E_BV' = E_BV + 0.0153
        #
        # Ag = cg * 3.1 * E_BV; cg = 0.829
        # Using the corrected E_BV':
        # Ag' = cg * 3.1 * E_BV' = cg * 3.1 * E_BV + (cg * 3.1 * 0.0153) -->
        # Ag' = Ag + 0.03932
        # The observed magnitude is created as:
        # m_obs = M_int + Ax + dist_mod
        # using the corrected Ag':
        # G_obs = G_int + Ag + 0.03932 + dist_mod
        # moving the correction into the distance modulus:
        # dist_mod' = 0.03932 + dist_mod

        if offset_flag is None:
            offset = 0.
        elif offset_flag == 'indiv':
            # Individual Gaia offsets
            mu_offset = {
                'vdBH85': 0.004, 'vdBH106': -0.027, 'TR12': -0.022,
                'RUP85': 0.028, 'vdBH73': 0.014, 'vdBH87': 0.041,
                'TR13': 0.044, 'NGC4349': 0.02, 'vdBH92': 0.034,
                'RUP162': 0.035}
            offset = mu_offset[cl] * 3.1 * 0.829
        elif offset_flag == 'mean':
            # General mean offset
            offset = 0.03932

        # Apply offset
        dm = dm + offset

        plx_asteca = 1. / (10 ** ((dm + 5.) / 5.))
        pc50_asteca = 1. / plx_asteca
        e_pc_asteca = .2 * np.log(10.) * pc50_asteca * e_dm
        # e_plx_asteca = 1. / e_pc_asteca**2
        pc16_asteca = pc50_asteca - e_pc_asteca
        pc84_asteca = pc50_asteca + e_pc_asteca

        # Lindegren parsec distances
        pc16_lindegren, pc50_lindegren, pc84_lindegren = 1000. * np.array(
            prsc_lindegren[cl])

        if cl != 'vdBH85':
            # No offset parallax distances in [mas]
            plx16_nooffset, plx50_nooffset, plx84_nooffset =\
                1. / (1000. * np.array(prsc_nooffset[cl]))
            # e_plx_nooffset = .5 * (plx16_nooffset + plx84_nooffset)

            # Order: ASteCA, no offset
            plx_data[cl] = np.array([
                plx_asteca, plx50_nooffset, plx16_nooffset, plx84_nooffset,
                log_age])

            # Schonrich parsec distances
            pc16_schonrich, pc50_schonrich, pc84_schonrich = 1000. * np.array(
                prsc_schonrich[cl])

            # Shuangjing parsec distances
            pc16_shuangjing, pc50_shuangjing, pc84_shuangjing =\
                1000. * np.array(prsc_shuangjing[cl])

            # No offset parsec distances
            pc16_nooffset, pc50_nooffset, pc84_nooffset = 1000. * np.array(
                prsc_nooffset[cl])
        else:
            plx_data[cl] = np.array([np.nan for _ in range(5)])
            pc16_schonrich, pc50_schonrich, pc84_schonrich, pc16_shuangjing,\
                pc50_shuangjing, pc84_shuangjing, pc16_nooffset,\
                pc50_nooffset, pc84_nooffset = [np.nan for _ in range(9)]

        # Order: ASteCA, Lindegren, Shuangjing, no offset
        parsec_data[cl] = np.array([
            [pc16_asteca, pc50_asteca, pc84_asteca, log_age],
            [pc16_lindegren, pc50_lindegren, pc84_lindegren],
            [pc16_schonrich, pc50_schonrich, pc84_schonrich],
            [pc16_shuangjing, pc50_shuangjing, pc84_shuangjing],
            [pc16_nooffset, pc50_nooffset, pc84_nooffset]])

    asteca_lindegren = np.array(
        [(_[0][1], _[1][1]) for cl, _ in parsec_data.items()
         if cl != 'vdBH85']).T
    print("(ASteCA - Lindegren) mean difference without vdBH85: {:.0f}".format(
        np.mean(asteca_lindegren[0] - asteca_lindegren[1])))

    # Make final plot
    print("Plotting")
    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(10, 10)
    # plt.style.use('seaborn-darkgrid')
    xymin, xymax = 1550., 5995.

    if all_flag:
        # Lindegren + Schönrich + Xu plot
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
        plt.savefig('dist_comparision_all.png', dpi=150, bbox_inches='tight')
        plt.clf()
        plt.close("all")
        return

    # fig = plt.figure(figsize=(25, 25))
    # gs = gridspec.GridSpec(20, 10)
    fig = plt.figure(figsize=(35, 25))
    gs = gridspec.GridSpec(10, 12)
    # plt.style.use('seaborn-darkgrid')

    dpcPlot(
        gs, (0, 2, 0, 2), xymin, xymax, 'No bias correction', parsec_data,
        (True, True), 4)

    ax = plt.subplot(gs[0:2, 2:4])
    ax.tick_params(axis='both', which='major', labelsize=8)
    for cl, plxs_d in plx_data.items():
        # plxs = (ASteCA, Bayes)
        x, y, _, _, _ = plxs_d * 1000.
        xo, yo = -0.01, 0.0035
        if cl == 'RUP85':
            xo, yo = -0.02, -.007
        elif cl == 'NGC4349':
            xo = -0.035
        elif cl == 'vdBH106':
            xo, yo = 0.015, -.001
        elif cl == 'vdBH73':
            xo, yo = -0.028, -.0065
        elif cl == 'TR13':
            xo, yo = 0.01, -.003
        cl = cl + '*' if cl in ('RUP162', 'vdBH106') else cl
        ax.annotate(cl, (y + xo, x - y + yo), fontsize=9)

    x, y, eyl, eyh, lage = np.array(list(plx_data.values())).T * 1000.
    lage = lage / 1000.
    # Reverse 16th and 84th percentiles
    eyl, eyh = eyl - y, y - eyh
    deltas = x - y
    plt.scatter(y, deltas, zorder=4, c=lage, marker='s', s=50)

    _16th, _84th = np.nanpercentile(deltas, 16), np.nanpercentile(deltas, 84)
    plt.axhline(_16th, c='blue', ls='--', lw=.7, zorder=1)
    plt.axhline(_84th, c='blue', ls='--', lw=.7, zorder=1)
    plt.axhline(
        np.nanmean(deltas), c='r', ls='--', lw=.7, zorder=1,
        label=r"$\Delta_{{mean}}$={:.3f}$_{{{:.3f}}}^{{{:.3f}}}$ [mas]".format(
            np.nanmean(deltas), _16th, _84th))
    # plt.axhline(
    #     np.nanmedian(deltas), c='g', ls='--', lw=.7, zorder=1,
    #     label=r"$\Delta_{{median}}$={:.3f} [mas]".format(
    #         np.nanmedian(deltas)))
    plt.errorbar(
        y, deltas, xerr=np.array([eyl, eyh]),
        fmt='.', elinewidth=.85, ms=4, ecolor='b', label=None)
    # plt.grid(ls=':', lw=.5, c='grey')
    # plt.xlabel(r'$Plx_{ASteCA}$ [mas]', fontsize=12)
    plt.xlabel(r'$Plx_{Bayes}$ [mas]', fontsize=12)
    plt.ylabel(r'$\Delta\;Plx_{[ASteCA-Bayes]}$ [mas]', fontsize=12)
    plt.legend(loc=2, fontsize=10)
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_ticks([7.5, 8., 8.5, 9., 9.5])

    # Lindegren plot
    dpcPlot(
        gs, (0, 2, 4, 6), xymin, xymax, 'Lindegren et al. (0.029 mas)',
        parsec_data, (True, False), 1)

    fig.tight_layout()
    plt.savefig('dist_comparision.png', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close("all")


def dpcPlot(gs, gsi, xymin, xymax, txt, parsec_data, labels, j):
    ax = plt.subplot(gs[gsi[0]:gsi[1], gsi[2]:gsi[3]])
    ax.tick_params(axis='both', which='major', labelsize=8)
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
            ecolor='b', label=None, zorder=1)
        xo, yo = 50., 50.
        if cl == 'vdBH87':
            xo, yo = -80., 120.
        elif cl == 'NGC4349':
            xo, yo = 120., 0.
        elif cl == 'RUP162':
            xo, yo = -650., -130.
        elif cl == 'TR13':
            xo, yo = -400., 40.
        elif cl == 'RUP85':
            xo, yo = 100., -170.
        elif cl == 'vdBH106':
            xo, yo = -700., 50.
        elif cl == 'vdBH73':
            xo, yo = 40., 50.
        cl = cl + '*' if cl in ('RUP162', 'vdBH106') else cl
        ax.annotate(cl, (x + xo, y + yo), fontsize=9, zorder=4)
    x, y, lage = np.array(xyall).T
    plt.scatter(x, y, marker='s', s=20, c=lage, zorder=6)  # , cmap='seismic')
    cbar = plt.colorbar(pad=.01, fraction=.02, aspect=50)
    cbar.ax.tick_params(labelsize=9)
    ax.set_title(
        r"{}, $\Delta_{{mean}}=${:.0f} [pc]".format(
            txt, np.nanmean(deltas)), fontsize=10)  # , x=.1, y=.95)
    plt.plot((xymin, xymax), (xymin, xymax), ls='--', c='r', lw=.5)
    # plt.grid(ls=':', lw=.5, c='grey')
    plt.xlabel(r'$d_{ASteCA}$ [pc]', fontsize=12)
    if labels[0]:
        plt.ylabel(r'$d_{Bayes}$ [pc]', fontsize=12)
    else:
        ax.set_yticklabels([])
    if labels[1]:
        cbar.set_label(r"$\log(age)$", fontsize=10, labelpad=4)
    cbar.set_ticks([7.5, 8., 8.5, 9., 9.5])
    plt.xlim(xymin, xymax)
    plt.ylim(xymin, xymax)


if __name__ == '__main__':
    prsc_lindegren = {
        # Weighted average estimate using Bailer-Jones distances for the
        # selected members of this cluster: 4000 +- 1300 pc
        'vdBH85': (2.769, 4.153, 5.536),
        'vdBH73': (4.475, 4.919, 5.287),
        'vdBH92': (2.331, 2.434, 2.514),
        'RUP85': (4.447, 4.640, 4.831),
        'vdBH106': (4.372, 4.771, 5.157),
        'TR13': (4.439, 4.584, 4.711),
        'NGC4349': (1.920, 1.920, 1.940),  # Actually uniform prior
        'RUP162': (4.185, 4.365, 4.553),
        'TR12': (3.505, 3.631, 3.758),
        'vdBH87': (2.200, 2.261, 2.323),
        #
        'LYNGA15': (3.113, 3.312, 3.559),
        'LODEN565': (2.748, 2.915, 3.105),
        'NGC4230': (2.359, 2.694, 3.089),
        'RUP88': (4.593, 5.090, 5.574),
        'RUP87': (1.482, 1.601, 1.810),
        'vdBH91': (2.739, 2.890, 3.045)
    }

    prsc_shuangjing = {
        # 'vdBH85': (5.341, 5.645, 5.954),
        'vdBH73': (3.708, 4.054, 4.372),
        'vdBH92': (2.099, 2.173, 2.248),
        'RUP85': (3.695, 3.830, 3.972),
        'vdBH106': (3.740, 4.057, 4.345),
        'TR13': (3.660, 3.751, 3.847),
        'NGC4349': (1.752, 1.762, 1.773),  # Actually uniform prior
        'RUP162': (3.542, 3.662, 3.806),
        'TR12': (3.017, 3.109, 3.207),
        'vdBH87': (1.996, 2.048, 2.095),
        #
        'NGC4230': (2.095, 2.407, 2.750),
        'LYNGA15': (2.701, 2.874, 3.058),
        'vdBH91': (2.430, 2.547, 2.688),
        'LODEN565': (2.415, 2.562, 2.704),
        'RUP88': (3.778, 4.153, 4.598),
        'RUP87': (1.359, 1.478, 1.643)
    }

    prsc_schonrich = {
        # 'vdBH85': (5.982, 6.283, 6.619),
        'vdBH73': (4.147, 4.460, 4.769),
        'vdBH92': (2.216, 2.282, 2.351),
        'RUP85': (4.013, 4.164, 4.319),
        'vdBH106': (3.995, 4.311, 4.663),
        'TR13': (3.983, 4.100, 4.207),
        'NGC4349': (1.812, 1.832, 1.842),  # Actually uniform prior
        'RUP162': (3.798, 3.936, 4.101),
        'TR12': (3.223, 3.314, 3.424),
        'vdBH87': (2.077, 2.133, 2.185),
        #
        'LYNGA15': (2.847, 3.045, 3.258),
        'RUP88': (4.080, 4.551, 5.055),
        'RUP87': (1.437, 1.560, 1.778),
        'LODEN565': (2.533, 2.691, 2.876),
        'vdBH91': (2.561, 2.707, 2.830),
        'NGC4230': (2.198, 2.540, 2.926)
    }

    prsc_nooffset = {
        # MP>0. and a Gaussian prior
        'vdBH85': (7.798, 8.228, 8.660),
        'vdBH73': (5.041, 5.482, 5.931),
        'vdBH92': (2.499, 2.606, 2.722),
        'RUP85': (5.160, 5.393, 5.619),  # Actually uniform prior
        'vdBH106': (5.034, 5.407, 5.804),
        'TR13': (5.100, 5.250, 5.421),
        'NGC4349': (2.010, 2.040, 2.060),  # Actually uniform prior
        'RUP162': (4.788, 4.973, 5.179),
        'TR12': (3.941, 4.084, 4.221),
        'vdBH87': (2.352, 2.423, 2.487),  # Actually uniform prior
        #
        'LYNGA15': (3.442, 3.681, 3.898),
        'LODEN565': (3.015, 3.248, 3.500),
        'RUP88': (5.169, 5.722, 6.211),
        'RUP87': (5.468, 6.188, 6.975),  # Actually uniform prior
        'vdBH91': (2.987, 3.157, 3.333),
        'NGC4230': (2.593, 2.968, 3.384)
    }

    dm_asteca = {
        'vdBH85': (13.317, 0.12085, 9.872),
        'vdBH106': (13.437, 0.36267, 9.482),
        'TR12': (12.721, 0.08997, 8.822),
        'RUP85': (13.405, 0.11791, 8.253),
        'vdBH73': (13.498, 0.26366, 8.893),
        'vdBH87': (11.588, 0.098, 8.399),
        'TR13': (13.411, 0.15, 8.044),
        'NGC4349': (11.375, 0.11, 8.470),
        'vdBH92': (12.066, 0.092136, 7.276),
        'RUP162': (13.23, 0.10038, 8.924),
        #
        # 'LYNGA15': (11.741, 0.15512, 9.787),
        # 'RUP88': (13.704, 0.52557, 9.741),
        # 'NGC4230': (13.167, 0.58768, 9.904),
        # 'LODEN565': (13.247, 0.251, 7.019),
        # 'RUP87': (12.96, 0.26, 9.49),
        # 'vdBH91': (11.03, 0.17141, 8.916)
    }


if __name__ == '__main__':
    offset_flag = None # 'mean' #'indiv'
    all_flag = False
    main(
        dm_asteca, prsc_nooffset, prsc_lindegren, prsc_schonrich,
        prsc_shuangjing, offset_flag, all_flag)
