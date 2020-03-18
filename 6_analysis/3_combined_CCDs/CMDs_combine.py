
from os import getcwd, walk
from os.path import join, realpath, dirname
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


names_dict = {
    # 'ngc4230_match_D2.png': 'NGC4230', 'rup85_match_D2.png': 'RUP85',
    # 'bh106_match_D2.png': 'vdBH106', 'trumpler13_match_D2.png': 'TR13',
    # 'rup88_match_D2.png': 'RUP88', 'rup87_match_D2.png': 'RUP87',
    # 'rup162_match_D2.png': 'RUP162',
    # 'bh91_match_D2.png': 'vdBH91', 'trumpler12_match_D2.png': 'TR12',
    # 'lynga15_match_D2.png': 'LYNGA15', 'bh87_match_D2.png': 'vdBH87',
    # 'bh85_match_D2.png': 'vdBH85', 'loden565_match_D2.png': 'LODEN565',
    # 'bh92_match_D2.png': 'vdBH92', 'bh73_match_D2.png': 'vdBH73',
    'ngc4349_match_D2.png': 'NGC4349'}

final_params = {
    'bh85_match': (0.014, 0.002, 9.87, 0.05, 0.30, 0.03, 13.32, 0.12, 2200, 500),
    'bh106_match': (0.012, 0.003, 9.48, 0.12, 0.30, 0.04, 13.44, 0.36, 500, 200),
    'trumpler12_match': (0.009, 0.002, 8.82, 0.06, 0.31, 0.03, 12.72, 0.09, 700, 100),
    'lynga15_match': (0.027, 0.007, 9.79, 0.08, 0.17, 0.05, 11.74, 0.16, 400, 130),
    'rup85_match': (0.027, 0.004, 8.44, 0.13, 1.07, 0.04, 13.56, 0.14, 3500, 520),
    'bh73_match': (0.019, 0.004, 8.89, 0.05, 1.06, 0.04, 13.50, 0.26, 2600, 950),
    'bh87_match': (0.021, 0.003, 8.31, 0.13, 0.56, 0.02, 11.44, 0.07, 1100, 120),
    'rup87_match': (0.015, 0.005, 9.70, 0.20, 0.09, 0.05, 11.64, 0.89, 110, 180),
    'trumpler13_match': (0.007, 0.004, 8.04, 0.07, 0.56, 0.02, 13.41, 0.15, 800, 250),
    'ngc4230_match': (0.008, 0.006, 9.90, 0.25, 0.12, 0.05, 13.17, 0.59, 180, 370),
    'rup88_match': (0.017, 0.006, 9.74, 0.25, 0.15, 0.06, 13.70, 0.53, 320, 450),
    'ngc4349_match': (0.01143, 0.00438, 8.47, 0.13, 0.41, 0.05, 11.38, 0.11, 2000, 60),
    'loden565_match': (0.021, 0.005, 7.02, 0.18, 0.85, 0.05, 13.25, 0.25, 330, 70),
    'bh92_match': (0.009, 0.004, 7.28, 0.33, 0.65, 0.03, 12.07, 0.09, 440, 80),
    'rup162_match': (0.009, 0.002, 8.92, 0.09, 0.54, 0.03, 13.23, 0.10, 1200, 200),
    'bh91_match': (0.009, 0.002, 8.92, 0.27, 0.20, 0.05, 11.03, 0.17, 200, 80),
}


def main():
    """
    """
    r_path = realpath(join(getcwd(), dirname(__file__)))

    # Path to UB-BV CCD
    path = '5_ASteCA/D1_D2/1st_run/'
    p_UBBV = join('/'.join(r_path.split('/')[:-2]), path)

    # Path to G-BV-VI CCD
    # path = '5_ASteCA/D1_D2/3rd_run/'
    path = '5_ASteCA/D1_D2/5th_run/'
    p1 = join('/'.join(r_path.split('/')[:-2]), path)

    for dirpath, dirnames, filenames in walk(p1):
        for fn in [f for f in filenames if f.endswith("_D2.png")]:
            print(fn)

            # Read entire image file
            im = plt.imread(join(p_UBBV, dirpath.split('/')[-1], fn))
            # Cut out UB-BV CCD
            xmin, ymin, w, h = 15, 1605, 717, 740  # 1570
            UB_BV = im[ymin:ymin + h, xmin:xmin + w]
            # plt.imshow(UB_BV);plt.show()
            cluster = [UB_BV]

            # Read entire image file
            im = plt.imread(join(dirpath, fn))

            xmin, ymin, w, h = 18, 90, 717, 750  # 55
            V_VI = im[ymin:ymin + h, xmin:xmin + w]

            ymin = 845  # 810
            V_BV = im[ymin:ymin + h, xmin:xmin + w]

            cluster += [V_BV, V_VI]

            makePlots(names_dict[fn], cluster, final_params[fn[:-7]])


def makePlots(name, data, params):
    fig = plt.figure(figsize=(10.5, 10.5))
    gs = gridspec.GridSpec(4, 4)

    ax = plt.subplot(gs[0])
    ax.imshow(data[-1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    ax = plt.subplot(gs[1])
    ax.imshow(data[0])
    ax.set_title(name, fontsize="5", x=.55, y=.9,
                 bbox=dict(facecolor='white', alpha=0.75))
    txt = (
        r"$z={:.4f}\pm{:.3f}$" + '\n' + r"$\log(age)={:.2f}\pm{:.2f}$" + '\n' +
        r"$E_{{(B-V)}}={:.2f}\pm{:.2f}$" + '\n' + r"$\mu_0={:.2f}\pm{:.2f}$" +
        '\n' + r"$M(M_{{\odot}})={:.0f}\pm{:.0f}$").format(*params)
    ax.text(
        0.19, 0.15, txt, transform=ax.transAxes, fontsize="5.5",
        bbox=dict(facecolor='white', alpha=0.85))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    ax = plt.subplot(gs[2])
    ax.imshow(data[-2])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    fig.tight_layout()
    fig.savefig('cmds_' + name + '.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
