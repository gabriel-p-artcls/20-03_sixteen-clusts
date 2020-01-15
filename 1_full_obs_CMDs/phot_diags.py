
from astropy.io import ascii
import os
from os import getcwd
from os.path import join, realpath, dirname
import numpy as np
from scipy.spatial import distance
# from scipy import stats
import warnings
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import matplotlib.offsetbox as offsetbox


def main():
    """
    Generate CMDs of full observed frames.
    """

    zams = ascii.read("zams.txt")
    z_UB, z_BV = zams['(U-B)o'], zams['(B-V)o']

    names_dict = {
        'ngc4230_match.dat': 'NGC4230', 'rup85_match.dat': 'RUP85',
        'bh106_match.dat': 'vdBH106', 'trumpler13_match.dat': 'TR13',
        'rup88_match.dat': 'RUP88', 'rup87_match.dat': 'RUP87',
        'ngc4349_match.dat': 'NGC4349', 'rup162_match.dat': 'RUP162',
        'bh91_match.dat': 'vdBH91', 'trumpler12_match.dat': 'TR12',
        'lynga15_match.dat': 'LYNGA15', 'bh87_match.dat': 'vdBH87',
        'bh85_match.dat': 'vdBH85', 'loden565_match.dat': 'LODEN565',
        'bh92_match.dat': 'vdBH92', 'bh73_match.dat': 'vdBH73'}

    rootdir = '/'.join(
        realpath(join(getcwd(), dirname(__file__))).split('/')[:-1]) +\
        '/3_new_match/2_cross_match_gaiadr2'
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            filepath = subdir + os.sep + file

            if filepath.endswith((".dat")):
                data = ascii.read(filepath)
                x, y, G, BV, UB, VI = data['x'], data['y'], data['Gmag'],\
                    data['BV'], data['UB'], data['VI']
                name = names_dict[file]
                makePlots(z_UB, z_BV, file, name, x, y, G, BV, UB, VI)


def makePlots(z_UB, z_BV, file, name, x, y, G, BV, UB, VI):
    print(name)

    cents = {
        'bh73': (2085, 2056, 305), 'bh85': (1048, 1425, 450),
        'bh87': (1089, 1283, 700), 'bh91': (1100, 1300, 500),
        'bh92': (931, 1284, 380), 'bh106': (1100, 1100, 500),
        'loden565': (1400, 1140, 400), 'lynga15': (2520, 1690, 600),
        'ngc4349': (1915, 2213, 600), 'ngc4230': (970, 970, 400),
        'rup85': (1090, 1320, 440), 'rup87': (1020, 1420, 300),
        'rup88': (1250, 1050, 300), 'rup162': (1300, 1500, 600),
        'trumpler12': (1000, 1150, 400), 'trumpler13': (1170, 1280, 500)
    }
    file = file.replace('_match.dat', '').replace('_match.dat', '')

    fig = plt.figure(figsize=(30, 25))
    gs = gridspec.GridSpec(10, 12)

    ax = plt.subplot(gs[2:4, 0:2])
    mask = (VI < 50.)

    coords = np.array([x[mask], y[mask]]).T
    dists = distance.cdist(
        np.array([[cents[file][0], cents[file][1]]]), coords)
    clmsk = dists[0] <= cents[file][2]

    Gf, VIf = G[mask], VI[mask]
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'mag', VIf, Gf)
    plt.scatter(
        VIf[~clmsk], Gf[~clmsk], c='grey', s=5, alpha=.5, zorder=1,
        label=r"$N_{{field}}={}$".format(len(VIf[~clmsk])))
    plt.scatter(
        VIf[clmsk], Gf[clmsk], c='k', s=15, zorder=3, lw=0.2, edgecolor='w',
        label=r"$N_{{clust}}={}$".format(len(VIf[clmsk])))
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$(V-I)$', fontsize=10)
    plt.ylabel('$G$', fontsize=10)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    plt.legend(loc='upper right')

    # ax = plt.subplot(gs[1])
    ax = plt.subplot(gs[2:4, 2:4])
    mask = (BV < 50.) & (UB < 50.)

    coords = np.array([x[mask], y[mask]]).T
    dists = distance.cdist(
        np.array([[cents[file][0], cents[file][1]]]), coords)
    clmsk = dists[0] <= cents[file][2]

    BVf, UBf = BV[mask], UB[mask]
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'col', BVf, UBf)
    plt.plot(z_BV, z_UB, c='r', ls='--', zorder=4)
    # plt.scatter(BVf, UBf, c='k', s=5, zorder=3)

    plt.scatter(
        BVf[~clmsk], UBf[~clmsk], c='grey', s=5, alpha=.5, zorder=1,
        label=r"$N_{{field}}={}$".format(len(BVf[~clmsk])))
    plt.scatter(
        BVf[clmsk], UBf[clmsk], c='k', s=15, zorder=3, lw=0.2, edgecolor='w',
        label=r"$N_{{clust}}={}$".format(len(BVf[clmsk])))

    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$(B-V)$', fontsize=10)
    plt.ylabel('$(U-B)$', fontsize=10)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    plt.legend(loc='upper right')
    # Add text box.
    text = '{}'.format(name)
    ob = offsetbox.AnchoredText(text, pad=0.2, loc=9, prop=dict(size=12))
    ob.patch.set(alpha=0.7)
    ax.add_artist(ob)

    ax = plt.subplot(gs[2:4, 4:6])
    mask = (BV < 50.)

    coords = np.array([x[mask], y[mask]]).T
    dists = distance.cdist(
        np.array([[cents[file][0], cents[file][1]]]), coords)
    clmsk = dists[0] <= cents[file][2]

    Gf, BVf = G[mask], BV[mask]
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'mag', BVf, Gf)
    plt.scatter(
        BVf[~clmsk], Gf[~clmsk], c='grey', s=5, alpha=.5, zorder=1,
        label=r"$N_{{field}}={}$".format(len(BVf[~clmsk])))
    plt.scatter(
        BVf[clmsk], Gf[clmsk], c='k', s=15, zorder=3, lw=0.2, edgecolor='w',
        label=r"$N_{{clust}}={}$".format(len(BVf[clmsk])))
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$(B-V)$', fontsize=10)
    plt.ylabel('$G$', fontsize=10)
    # Set minor ticks
    ax.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    plt.legend(loc='upper right')

    fig.tight_layout()
    fig.savefig('obs_' + name + '.png', dpi=150, bbox_inches='tight')
    plt.close()


def diag_limits(yaxis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    x_median, x_std = np.nanmedian(phot_x), np.nanstd(phot_x)
    x_min_cmd, x_max_cmd = x_median - 4.5 * x_std, x_median + 4.5 * x_std

    # Use stars within the x limits defined. This prevents stars far away
    # from the x median from affecting the limit in y.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xmsk = (phot_x > x_min_cmd) & (phot_x < x_max_cmd)

    phot_y_msk = np.array(phot_y)[xmsk]
    y_median, y_std = np.nanmedian(phot_y_msk), np.nanstd(phot_y_msk)

    # y limits.
    if yaxis == 'mag':
        y_min_cmd = y_median + 1.25 * y_std + .75
        # If photometric axis y is a magnitude, make sure the brightest star
        # is always plotted.
        y_max_cmd = np.nanmin(phot_y_msk) - 1.
    else:
        y_max_cmd, y_min_cmd = y_median - 4.5 * y_std, y_median + 4.5 * y_std

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


# def diag_limits(yaxis, phot_x, phot_y):
#     '''
#     Define plot limits for *all* photometric diagrams.
#     '''
#     x_v, y_v = kde_limits(phot_x, phot_y)

#     # Define diagram limits.
#     x_min_cmd, x_max_cmd = min(x_v) - 1., max(x_v) + 1.
#     y_min_cmd = max(y_v) + 1.
#     # If photometric axis y is a magnitude, make sure the brightest star
#     # is always plotted.
#     if yaxis == 'mag':
#         y_max_cmd = min(phot_y) - 1.
#     else:
#         y_max_cmd = min(y_v) - 1.

#     return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


# def kde_limits(phot_x, phot_y):
#     '''
#     Return photometric diagram limits taken from a 2D KDE.
#     '''

#     xmin, xmax = min(phot_x), max(phot_x)
#     ymin, ymax = min(phot_y), max(phot_y)
#     # Stack photometric data.
#     values = np.vstack([phot_x, phot_y])
#     # Obtain Gaussian KDE.
#     kernel = stats.gaussian_kde(values)
#     # Grid density (number of points).
#     gd = 25
#     gd_c = complex(0, gd)
#     # Define x,y grid.
#     x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
#     positions = np.vstack([x.ravel(), y.ravel()])
#     # Evaluate kernel in grid positions.
#     k_pos = kernel(positions)

#     # Generate 30 contour lines.
#     plt.figure()
#     cs = plt.contour(x, y, np.reshape(k_pos, x.shape), 30)
#     plt.close()
#     # Extract (x,y) points delimiting each line.
#     x_v, y_v = np.asarray([]), np.asarray([])
#     # Only use the outer curve.
#     col = cs.collections[0]
#     # If more than one region is defined by this curve (ie: the main sequence
#     # region plus a RC region or some other detached region), obtain x,y from
#     # all of them.
#     for lin in col.get_paths():
#         x_v = np.append(x_v, lin.vertices[:, 0])
#         y_v = np.append(y_v, lin.vertices[:, 1])

#     return x_v, y_v


if __name__ == '__main__':
    main()
