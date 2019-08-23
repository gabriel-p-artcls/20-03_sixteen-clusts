
from os import getcwd, walk
from os.path import join, realpath, dirname
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


names_dict = {
    'ngc4230_match_D2.png': 'NGC4230', 'rup85_match_D2.png': 'RUP85',
    'bh106_match_D2.png': 'vdBH106', 'trumpler13_match_D2.png': 'TR13',
    'rup88_match_D2.png': 'RUP88', 'rup87_match_D2.png': 'RUP87',
    'ngc4349_match_D2.png': 'NGC4349', 'rup162_match_D2.png': 'RUP162',
    'bh91_match_D2.png': 'vdBH91', 'trumpler12_match_D2.png': 'TR12',
    'lynga15_match_D2.png': 'LYNGA15', 'bh87_match_D2.png': 'vdBH87',
    'bh85_match_D2.png': 'vdBH85', 'loden565_match_D2.png': 'LODEN565',
    'bh92_match_D2.png': 'vdBH92', 'bh73_match_D2.png': 'vdBH73'}

final_params = {
    'vdBH87': (0.024994, 0.0042559, 8.505, 0.2164, 0.58, 0.0112, 11.272, 0.059715, 954, 252.02),
    'vdBH91': (0.02118, 0.0047991, 8.631, 0.2927, 0.746, 0.0904, 13.399, 0.245, 473, 126.66),
    'vdBH106': (0.0088567, 0.0040885, 9.987, 0.03759, 0.296, 0.0254, 13.997, 0.21338, 1845, 349.8),
    'RUP88': (0.013038, 0.0051975, 9.541, 0.1185, 0.331, 0.056, 14.746, 0.37555, 1047, 293.66),
    'NGC4230': (0.007701, 0.0050347, 9.631, 0.07867, 0.21, 0.039, 15.176, 0.25646, 1360, 485.47),
    'LODEN565': (0.027462, 0.0066517, 8.1, 0.2375, 0.891, 0.0149, 12.595, 0.13702, 273, 69.569),
    'RUP87': (0.024552, 0.0051166, 9.142, 0.06687, 0.228, 0.0387, 13.522, 0.31632, 191, 57.954),
    'vdBH92': (0.0064851, 0.0029284, 7.92, 0.2157, 0.655, 0.033, 11.462, 0.12332, 399, 74.563),
    'RUP162': (0.0027856, 0.0015016, 7.485, 0.1558, 0.474, 0.0242, 12.251, 0.12524, 523, 92.22),
    'TR12': (0.006047, 0.0013128, 8.903, 0.09872, 0.288, 0.0275, 12.409, 0.11919, 718, 97.164),
    'vdBH85': (0.015203, 0.0023697, 9.869, 0.05458, 0.284, 0.0273, 13.09, 0.12321, 2767, 270.65),
    'RUP85': (0.025893, 0.006639, 8.615, 0.06377, 1.06, 0.0417, 13.479, 0.17536, 2727, 649.47),
    'LYNGA15': (0.0059159, 0.00086094, 9.938, 0.06299, 0.271, 0.0228, 11.035, 0.081641, 491, 92.207),
    'vdBH73': (0.020808, 0.004578, 9.783, 0.08208, 0.623, 0.0368, 12.106, 0.29436, 978, 805.69),
    'NGC4349': (0.027674, 0.0058084, 8.598, 0.1348, 0.401, 0.0251, 12.194, 0.13254, 1106, 239.63),
    'TR13': (0.00726, 0.00335, 7.427, 0.177, 0.554, 0.014, 13.003, 0.124, 505, 148)}


def main():
    """
    """
    r_path = realpath(join(getcwd(), dirname(__file__)))

    # Path to UB-BV CCD
    path = 'UB_BV/'
    p_UBBV = join(r_path, path)

    # Path to G-BV-VI CCD
    path = '5_ASteCA/D_fundam_params/3rd_run/'
    p1 = join('/'.join(r_path.split('/')[:-2]), path)
    for dirpath, dirnames, filenames in walk(p1):
        for fn in [f for f in filenames if f.endswith("_D2.png")]:
            print(names_dict[fn])

            # Read entire image file
            im = plt.imread(join(p_UBBV, fn))
            # Cut out UB-BV CCD
            xmin, ymin, w, h = 15, 1570, 717, 740
            UB_BV = im[ymin:ymin + h, xmin:xmin + w]
            # plt.imshow(UB_BV);plt.show()
            cluster = [UB_BV]

            # Read entire image file
            im = plt.imread(join(dirpath, fn))

            xmin, ymin, w, h = 18, 55, 717, 750
            V_VI = im[ymin:ymin + h, xmin:xmin + w]

            ymin = 810
            V_BV = im[ymin:ymin + h, xmin:xmin + w]

            cluster += [V_BV, V_VI]

            makePlots(names_dict[fn], cluster, final_params[names_dict[fn]])


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
        r"$z={:.4f}\pm{:.5f}$" + '\n' + r"$\log(age)={:.2f}\pm{:.3f}$" + '\n' +
        r"$E_{{(B-V)}}={:.2f}\pm{:.2f}$" + '\n' + r"$\mu_0={:.1f}\pm{:.2f}$" +
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
