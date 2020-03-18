
from os import getcwd
from os.path import join, realpath, dirname
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image
from astropy import units as u
from astropy.coordinates import SkyCoord


"""
Convert PDF to PNG:

convert -density 300 bh73.PDF -quality 90 bh73.png

Then do one (or two or three) 'Autocrop' with Kolourpaint.
"""


def makePlot(N, _img):
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(4, 3)
    for i, g in enumerate(gs):
        ax = plt.subplot(g)
        try:
            ax.imshow(_img[i])
        except IndexError:
            pass

        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        # ax.set_aspect('auto')

    fig.tight_layout()
    fig.savefig('frames_' + str(N) + '.png', dpi=300, bbox_inches='tight')
    plt.close()


r_path = realpath(join(getcwd(), dirname(__file__))) + '/'

lon_lat = [
    ('n4230', 298.025, 7.445), ('rup87', 279.372, 4.883),
    ('bh106', 286.048, 4.7), ('bh85', 276.914, 4.544),
    ('loden565', 297.65, 1.71), ('bh73', 273.634, 0.951),
    ('n4349', 299.719, 0.83), ('bh92', 282.984, 0.438),
    ('rup85', 280.15, 0.16), ('bh87', 280.719, 0.059),
    ('lynga15', 295.053, -0.672), ('bh91', 284.03, -1.6),
    ('tr13', 285.515, -2.353), ('rup162', 289.638, -2.545),
    ('tr12', 283.828, -3.698), ('rup88', 286.661, -5.186)]

set1 = [[], []]
for cl, lon, lat in lon_lat:
    c = SkyCoord(lon * u.degree, lat * u.degree, frame='galactic')
    # print(c.fk5.ra.hms, c.fk5.ra, c.fk5.dec)
    # h_dec = c.fk5.ra.hms.h + (c.fk5.ra.hms.m + c.fk5.ra.hms.s / 60.) / 60.
    set1[0].append(cl)
    set1[1].append(c.l.deg)

# Order by lon
set1 = [x for _, x in sorted(zip(set1[1], set1[0]))]
# print(set1)

_img = []
# fix_w, fix_h = 1904, 1846
# for im in set1[:6]:
for im in set1:
    # _img.append(mimage.imread(r_path + im + '.png'))
    print(r_path + im + '.png')
    im_opn = Image.open(r_path + im + '.png')
    # print(im, im_opn.getbbox())
    # w, h = im.size
    # dx, dy = w - fix_w, h - fix_h
    # print(-dx, -dy, w, h)
    # im2 = im.crop((-dx, -dy, w, h))
    im2 = im_opn.crop(im_opn.getbbox())
    # print(im2.size)
    _img.append(im2)

makePlot(0, _img)

_img = []
for im in set1[12:]:
    print(r_path + im + '.png')
    im_opn = Image.open(r_path + im + '.png')
    # print(im, im_opn.getbbox())
    im2 = im_opn.crop(im_opn.getbbox())
    # print(im2.size)
    _img.append(im2)

makePlot(1, _img)

# _img = []
# for im in set1[12:]:
#     im_opn = Image.open(r_path + im + '.png')
#     print(im, im_opn.getbbox())
#     im2 = im_opn.crop(im_opn.getbbox())
#     print(im2.size)
#     _img.append(im2)

# makePlot(2, _img)
