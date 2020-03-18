
from os import getcwd
from os.path import join, realpath, dirname
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.table import Table
import numpy as np


def main():
    """
    Read  image files, extract the proper plots, and generate a single image
    for each cluster.
    """

    names_dict = {
        'rup85': 'RUP85', 'bh106': 'vdBH106', 'rup88': 'RUP88',
        'rup87': 'RUP87', 'bh91': 'vdBH91', 'trumpler12': 'TR12',
        'lynga15': 'LYNGA15', 'bh85': 'vdBH85', 'loden565': 'LODEN565',
        'bh92': 'vdBH92', 'bh73': 'vdBH73', 'ngc4349': 'NGC4349',
        'ngc4230': 'NGC4230', 'trumpler13': 'TR13', 'rup162': 'RUP162',
        'bh87': 'vdBH87'}

    # # Print the per-cluster V, B, BV mean differences
    # name = '16clusts_print'

    # # Plot the combined transformations for the observed clusters
    # name = '16clusts'

    # # Plot the combined transformations for Carrasco's Landolt set
    # name = 'carrasco'

    # Plot the combined transformations for APASS fields
    name = 'APASS'

    if name == 'carrasco':
        data_dict = loadCarrasco()

    elif name == 'APASS':
        data_dict = loadAPASS(names_dict)

    elif name == '16clusts':
        data_dict = myData(names_dict)

    elif name == '16clusts_print':
        for fname, name in names_dict.items():
            data_dict = myData({fname: name})
            transf_dict = transf(data_dict)
            printDeltas(name, data_dict, transf_dict)
        print("Finished")
        return

    transf_dict = transf(data_dict)
    makePlots(name, data_dict, transf_dict)
    print("Finished")


def myData(names_dict):
    """
    Load data on the 16 clusters
    """
    r_path = realpath(join(getcwd(), dirname(__file__)))
    path = '/'.join(r_path.split('/')[:-1]) +\
        '/3_new_match/2_cross_match_gaiadr2/'

    data_dict = {
        'Gm': [], 'BPRPm': [], 'Um': [], 'Bm': [], 'Vm': [], 'Im': [],
        'UBm': [], 'BVm': [], 'VIm': []}
    for fname, name in names_dict.items():
        fname = join(path, fname + '_match.dat')

        Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm = readData(fname, 'mydata')
        data_dict['Gm'] += Gm
        data_dict['BPRPm'] += BPRPm
        data_dict['Um'] += Um
        data_dict['Bm'] += Bm
        data_dict['Vm'] += Vm
        data_dict['Im'] += Im
        data_dict['UBm'] += UBm
        data_dict['BVm'] += BVm
        data_dict['VIm'] += VIm

        # print(fname.split('/')[-1], len(Gm))

    data = {}
    for k, v in data_dict.items():
        data[k] = np.array(v)

    return data


def loadAPASS(names_dict):
    """
    Load APASS data on a section of the 16 clusters
    """
    path = realpath(join(getcwd(), dirname(__file__), 'APASS_GAIA'))
    path = path.replace('in_repo/', '')

    data_dict = {
        'Gm': [], 'BPRPm': [], 'Um': [], 'Bm': [], 'Vm': [], 'Im': [],
        'UBm': [], 'BVm': [], 'VIm': []}
    for fname, name in names_dict.items():
        fname = join(path, fname + '_match_apass_match.dat')

        Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm = readData(fname, 'APASS')
        data_dict['Gm'] += Gm
        data_dict['BPRPm'] += BPRPm
        data_dict['Um'] += Um
        data_dict['Bm'] += Bm
        data_dict['Vm'] += Vm
        data_dict['Im'] += Im
        data_dict['UBm'] += UBm
        data_dict['BVm'] += BVm
        data_dict['VIm'] += VIm

        print(fname.split('/')[-1], len(Gm))

    data = {}
    for k, v in data_dict.items():
        data[k] = np.array(v)

    return data


def loadCarrasco():
    """
    Load Carrasco's list of cross-matched Landolt standards.
    """
    Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm = readData(
        'DR2G13xLandolt09-13.dat', 'carrasco')

    data_dict = {
        'Gm': Gm, 'BPRPm': BPRPm, 'Um': Um, 'Bm': Bm, 'Vm': Vm, 'Im': Im,
        'UBm': UBm, 'BVm': BVm, 'VIm': VIm}

    return data_dict


def readData(fname, dataID):
    """
    """

    t = Table.read(
        fname, format='ascii',
        fill_values=[('', '0'), ('NA', '0'), ('INDEF', '0')])

    if dataID == 'mydata':
        Gm, eGm, BPRPm, e_BP, e_RP, eEBPRP, G_BP, G_RP =\
            t['Gmag'], t['e_Gmag'], t['BP-RP'],\
            t['e_BPmag'], t['e_RPmag'], t['E_BR_RP_'], t['BPmag'], t['RPmag']
        UBm, BVm, VIm, Vm = t['UB'], t['BV'], t['VI'], t['V']
        eV, eUB, eBV, eVI = t['Verr'], t['UBerr'], t['BVerr'], t['VIerr']

        # Mask on our photometry
        msk_o = (Vm < 50.) & (UBm < 50.) & (BVm < 50.) & (VIm < 50.) &\
            (eV < .05) & (eUB < .05) & (eBV < .05) & (eVI < .05)
        # General mask
        msk_G, msk_eG, msk_eEBPRP = Gm < 13., eGm < 0.01,\
            eEBPRP < 1.5 + 0.03 * (G_BP - G_RP)**2
        msk_BPRP = (e_BP < 0.01) & (e_RP < 0.01)
        msk = msk_o & msk_G & msk_eG & msk_eEBPRP & msk_BPRP

        # Masked photom
        Gm, BPRPm, UBm, BVm, VIm, Vm = Gm[msk], BPRPm[msk], UBm[msk],\
            BVm[msk], VIm[msk], Vm[msk]
        Bm = BVm + Vm
        Um = UBm + Bm
        Im = Vm - VIm

        arrs = [list(_) for _ in (Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm)]
        return arrs

    if dataID == 'APASS':
        # Gaia data
        Gm, eGm, BPRPm, e_BP, e_RP, eEBPRP, G_BP, G_RP =\
            t['Gmag'], t['e_Gmag'], t['BP-RP'],\
            t['e_BPmag'], t['e_RPmag'], t['E_BR_RP_'], t['BPmag'], t['RPmag']
        # APASS data
        Vm, Bm = t['Johnson_V (V)'], t['Johnson_B (B)']
        # APASS does not contain this data in the Johnson system, so load SLOAN
        Um, Im = t['Sloan_u (SU)'], t['Sloan_i (SI)']

        # General GAIA mask
        msk_G, msk_eG, msk_eEBPRP = Gm < 13., eGm < 0.01,\
            eEBPRP < 1.5 + 0.03 * (G_BP - G_RP)**2
        msk_BPRP = (e_BP < 0.01) & (e_RP < 0.01)
        msknan = ~np.isnan(Gm.filled(np.nan))
        msk = msk_G & msk_eG & msk_eEBPRP & msk_BPRP & msknan

        # Masked photom
        Gm, BPRPm, Um, Bm, Im, Vm = Gm[msk], BPRPm[msk], Um[msk],\
            Bm[msk], Im[msk], Vm[msk]

        UBm, BVm, VIm = Um - Bm, Bm - Vm, Vm - Im

        arrs = [list(_) for _ in (Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm)]
        return arrs

    elif dataID == 'carrasco':
        G, eG, V, UB, BV, VR, RI, BP, e_BP, RP, e_RP =\
            t['phot_g_mean_mag'], t['phot_g_mean_mag_error'], t['Vmag'],\
            t['U-B'], t['B-V'], t['V-R'], t['R-I'], t['phot_bp_mean_mag'],\
            t['phot_bp_mean_mag_error'], t['phot_rp_mean_mag'],\
            t['phot_rp_mean_mag_error']

        B = BV + V
        VI = VR + RI
        U, Imag = UB + B, V - VI
        BPRP = BP - RP
        Rmag = RI + Imag

        # General mask
        msk = (G < 13.) & (eG < 0.01) & (e_BP < 0.01) & (e_RP < 0.01)
        Gm, Um, Bm, Vm, Im, Rm, BPRPm, UBm, BVm, VIm = [
            _[msk] for _ in (G, U, B, V, Imag, Rmag, BPRP, UB, BV, VI)]

        # from scipy.optimize import least_squares
        # def fun(coeffs, yd):
        #     """
        #     Obtain transformation coefficients for the U,B,I magnitude as:
        #     G - X = f(BP-RP)
        #     """
        #     v = (Gm - yd) - np.polyval(coeffs, BPRPm)
        #     return v

        # poly_order = 4
        # UBI_coeffs = []
        # for yd in (Bm,):  # Um, Vm, Rm, Im
        #     res = least_squares(fun, [.1] * poly_order, args=([yd]))
        #     UBI_coeffs.append(np.array(res.x))

        #     # Parameters errors. See:
        #     # https://stackoverflow.com/a/21844726/1391441
        #     # https://stackoverflow.com/a/14857441/1391441

        #     # Modified Jacobian
        #     J = res.jac
        #     # Reduced covariance matrix
        #     red_cov = np.linalg.inv(J.T.dot(J))
        #     # RMS of the residuals
        #     RMS = (res.fun**2).mean()
        #     # Covariance Matrix
        #     cov = red_cov * RMS
        #     # Standard deviation of the parameters
        #     p_std = np.sqrt(np.diagonal(cov))

        #     print(res.x)
        #     print(p_std)
        #     print(np.sqrt((res.fun**2).mean()))

        return [np.array(_) for _ in (
            Gm, BPRPm, Um, Bm, Vm, Im, UBm, BVm, VIm)]


def transf(data, Ninterp=500):
    """
    Carrasco Gaia DR2 transformations:
    https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/
    chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
    """
    Gm, BPRPm = data['Gm'], data['BPRPm']

    x = np.linspace(np.nanmin(BPRPm), np.nanmax(BPRPm), Ninterp)
    transf_dict = {'x': x}

    # G-U vs BP-RP
    U_coeffs = [-0.34971232, 2.26474081, -4.66713317, 2.69990085, -1.45412643,
                0.12850817]
    # [0.30457846, -1.14012079, 0.64847919, -1.5527081, 0.34591795]
    # [0.21240278, -1.01742603, -1.14593491, 0.48686377]
    U_poly = np.polyval(U_coeffs, x)
    U_Gaia = Gm - np.polyval(U_coeffs, BPRPm)
    transf_dict.update({'U_poly': U_poly, 'U_Gaia': U_Gaia})

    # G-B vs BP-RP
    B_coeffs = [0.06719282, -0.42375974, -0.63802333, 0.00282631]
    B_poly = np.polyval(B_coeffs, x)
    B_Gaia = Gm - np.polyval(B_coeffs, BPRPm)
    transf_dict.update({'B_poly': B_poly, 'B_Gaia': B_Gaia})

    # G-V vs BP-RP
    V_coeffs = [-0.1732, -0.006860, -0.01760]
    V_poly = np.polyval(V_coeffs, x)
    V_Gaia = Gm - np.polyval(V_coeffs, BPRPm)
    transf_dict.update({'V_poly': V_poly, 'V_Gaia': V_Gaia})

    # G-I vs BP-RP
    I_coeffs = [-0.09631, 0.7419, 0.02085]
    I_poly = np.polyval(I_coeffs, x)
    I_Gaia = Gm - np.polyval(I_coeffs, BPRPm)
    transf_dict.update({'I_poly': I_poly, 'I_Gaia': I_Gaia})

    return transf_dict


def printDeltas(name, data, transf):
    """
    """
    Bm, Vm, BV = data['Bm'], data['Vm'], data['BVm']

    minmax = .8

    delta_mag = []
    for mag, mag_name in ((Vm, 'V'), (Bm, 'B')):
        delta_M = transf[mag_name + '_Gaia'] - mag
        msk = (-minmax < delta_M) & (delta_M < minmax)
        delta_mag.append(delta_M[msk])

    BV_gaia = transf['B_Gaia'] - transf['V_Gaia']
    delta_col = BV_gaia - BV
    msk = (-minmax < delta_col) & (delta_col < minmax)
    delta_col = delta_col[msk]

    # V, B, BV
    print(
        r"{} & {:.3f}$\pm${:.3f} & {:.3f}$\pm${:.3f} & {:.3f}$\pm${:.3f} & {}".format(
            name, np.mean(delta_mag[0]), np.std(delta_mag[0]),
            np.mean(delta_mag[1]), np.std(delta_mag[1]),
            np.mean(delta_col), np.std(delta_col), len(delta_mag[0])))
        # len(delta_mag[1]), len(delta_col)))


def makePlots(name, data, transf):
    """
    """
    Gm, BPRPm, Um, Bm, Vm, Im, UB, BV, VI = data['Gm'], data['BPRPm'],\
        data['Um'], data['Bm'], data['Vm'], data['Im'], data['UBm'],\
        data['BVm'], data['VIm']

    plt.style.use('seaborn-darkgrid')
    plt.set_cmap('viridis')
    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(4, 4)

    def polyPlot(gs_ax, mag, mag_name):
        plt.subplot(gs[gs_ax])
        plt.xlabel(r"$BP-RP$", fontsize=12)
        plt.ylabel(r"$G-{}$".format(mag_name), fontsize=12)
        plt.plot(transf['x'], transf[mag_name + '_poly'], c='k', zorder=5)
        plt.scatter(BPRPm, Gm - mag, label="N={}".format(len(Gm)), c=Vm)
        plt.xlim(max(-.4, min(BPRPm) - .05), min(2.8, max(BPRPm) + .05))

    # # G-X vs BP-RP
    # polyPlot(0, Um, 'U')
    # polyPlot(1, Bm, 'B')
    # polyPlot(2, Vm, 'V')
    # polyPlot(3, Im, 'I')

    minmax = .8

    def magDiffs(gs_ax, mag, mag_name, minmax=minmax, BPRPm=BPRPm):
        plt.subplot(gs[gs_ax])
        plt.xlabel(r"${}$".format(mag_name), fontsize=12)  # _{{Gaia}}
        plt.ylabel(
            r"${}_{{Gaia}} - {}$".format(mag_name, mag_name), fontsize=12)
        delta_M = transf[mag_name + '_Gaia'] - mag

        msk = (-minmax < delta_M) & (delta_M < minmax)
        plt.title(
            "N={}, mask=({}, {})".format(msk.sum(), -minmax, minmax),
            fontsize=12)

        delta_M, BPRPm, mag = delta_M[msk], BPRPm[msk], mag[msk]
        plt.axhline(
            np.mean(delta_M), c='red', ls='--', lw=1.5, zorder=1,
            label=r"$\Delta {}_{{mean}}=${:.2f}$\pm${:.2f}".format(
                mag_name, np.mean(delta_M), np.std(delta_M)))
        plt.axhline(
            y=np.nanmedian(delta_M), ls=':', c='k',
            label="Median = {:.2f}".format(np.nanmedian(delta_M)))
        plt.scatter(mag, delta_M, s=7, c=BPRPm)
        plt.legend(fontsize=12)

    # Delta plots for magnitudes
    # magDiffs(4, Um, 'U')
    magDiffs(1, Bm, 'B')
    magDiffs(0, Vm, 'V')
    # magDiffs(7, Im, 'I')
    # magDiffs(0, Vm, 'V')
    # magDiffs(1, Bm, 'B')

    def colDiffs(gs_ax, col_data, col, minmax=minmax, Vm=Vm, BPRPm=BPRPm):
        ax = plt.subplot(gs[gs_ax])
        col_gaia = transf[col[0] + '_Gaia'] - transf[col[1] + '_Gaia']
        delta_col = col_gaia - col_data

        msk = (-minmax < delta_col) & (delta_col < minmax)
        plt.title(
            "N={}, mask=({}, {})".format(msk.sum(), -minmax, minmax),
            fontsize=12)
        plt.xlabel(r"$V}$", fontsize=12)
        plt.ylabel(r"${}_{{Gaia}} - {}$".format(col, col), fontsize=12)
        plt.axhline(
            np.mean(delta_col[msk]), c='red', ls='--', lw=1.5, zorder=1,
            label=r"$\Delta {}_{{mean}}=${:.2f}$\pm${:.2f}".format(
                col, np.mean(delta_col[msk]), np.std(delta_col[msk])))
        plt.axhline(
            y=np.nanmedian(delta_col[msk]), ls=':', c='k',
            label="Median = {:.2f}".format(np.nanmedian(delta_col[msk])))
        im = plt.scatter(Vm[msk], delta_col[msk], s=7, c=BPRPm[msk])
        plt.legend(fontsize=12)

        if col == 'BV':
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='2%', pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            cbar.ax.set_ylabel('(BP-RP)', fontsize=10)
            cbar.ax.tick_params(labelsize=8)

    # Delta plots for colors
    # colDiffs(8, UB, 'UB')
    colDiffs(2, BV, 'BV')
    # colDiffs(10, VI, 'VI')
    # colDiffs(2, BV, 'BV')

    fig.tight_layout()
    fig.savefig(
        join('gaia_transf_' + name + '.png'), dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
