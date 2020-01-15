
from astropy.io import ascii
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec


#
path = "/media/gabriel/rest/Dropbox_nosync/Papers/2019/11-sixteen_clusts/" +\
    "in_repo/3_new_match/2_cross_match_gaiadr2/"

cluster = 'rup162'

fra, fde = .05, .1
bw_kde = 0.05
nbins1 = 350

print(cluster)
fname = path + cluster + '_match.dat'
data = ascii.read(fname)
arr = np.array([data['pmRA'], data['pmDE']])

xmedian, ymedian = np.nanmedian(data['pmRA']), np.nanmedian(data['pmDE'])
xstddev, ystddev = np.nanstd(data['pmRA']), np.nanstd(data['pmDE'])
xmin, xmax = xmedian - xstddev, xmedian + xstddev
ymin, ymax = ymedian - ystddev, ymedian + ystddev

# Perform the kernel density estimate
xx, yy = np.mgrid[xmin:xmax:70j, ymin:ymax:70j]
positions = np.vstack([xx.ravel(), yy.ravel()])
values = np.vstack([data['pmRA'], data['pmDE']])
kernel = gaussian_kde(values, bw_method=bw_kde)
k_pos = kernel(positions)
kde = np.reshape(k_pos.T, xx.shape)
kde_max = positions.T[np.argmax(k_pos)]

print(kde_max)
print(xstddev, ystddev)

msk1 = (data['pmRA'] > kde_max[0] - fra * xstddev) &\
    (data['pmRA'] < kde_max[0] + fra * xstddev)
msk2 = (data['pmDE'] > kde_max[1] - fde * ystddev) &\
    (data['pmDE'] < kde_max[1] + fde * ystddev)

#
fig = plt.figure(figsize=(30, 25))
gs = gridspec.GridSpec(10, 12)

ax = plt.subplot(gs[0:2, 0:2])
ax.minorticks_on()
ext_range = [xmin, xmax, ymin, ymax]
im = plt.imshow(
    np.rot90(kde), cmap=plt.get_cmap('RdYlBu_r'), extent=ext_range)
plt.contour(xx, yy, kde, colors='k', linewidths=.5)
# (bottom, left), width, height
width, height = 2. * (fra * xstddev), 2. * (fde * ystddev)
rect = patches.Rectangle((
    np.min(data['pmRA'][msk1]), np.min(data['pmDE'][msk2])),
    width, height, linewidth=1, edgecolor='g', facecolor='none')
ax.add_patch(rect)
ax.set_aspect(aspect='auto')
plt.xlabel('pmRA')
plt.ylabel('pmDE')

ax = plt.subplot(gs[0:2, 2:4])
ax.minorticks_on()
plt.scatter(
    data['pmRA'], data['pmDE'], s=4, lw=.1, edgecolor='k')
plt.scatter(*kde_max, marker='x', s=25, c='r', lw=.5, zorder=5)
# (bottom, left), width, height
width, height = 2. * (fra * xstddev), 2. * (fde * ystddev)
rect = patches.Rectangle((
    np.min(data['pmRA'][msk1]), np.min(data['pmDE'][msk2])),
    width, height, linewidth=1, edgecolor='r', facecolor='none')
ax.add_patch(rect)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel('pmRA')
plt.ylabel('pmDE')

ax = plt.subplot(gs[0:1, 4:6])
ax.minorticks_on()
plt.hist(data['pmRA'], nbins1)
plt.axvline(
    kde_max[0], c='k',
    label=r'${:.3f}_{{{:.3f}}}^{{{:.3f}}}$'.format(
        kde_max[0], kde_max[0] - fra * xstddev, kde_max[0] + fra * xstddev))
plt.axvline(kde_max[0] - fra * xstddev, c='r')
plt.axvline(kde_max[0] + fra * xstddev, c='r')
plt.xlim(kde_max[0] - xstddev, kde_max[0] + xstddev)
plt.xlabel('pmRA')
plt.legend(fontsize=10)

ax = plt.subplot(gs[1:2, 4:6])
ax.minorticks_on()
plt.hist(data['pmDE'], nbins1)
plt.axvline(
    kde_max[1], c='k',
    label=r'${:.3f}_{{{:.3f}}}^{{{:.3f}}}$'.format(
        kde_max[1], kde_max[1] - fde * ystddev, kde_max[1] + fde * ystddev))
plt.axvline(kde_max[1] - fde * ystddev, c='r')
plt.axvline(kde_max[1] + fde * ystddev, c='r')
plt.xlim(kde_max[1] - ystddev, kde_max[1] + ystddev)
plt.xlabel('pmDE')
plt.legend(fontsize=10)

msk = msk1 & msk2

ax = plt.subplot(gs[2:4, 0:2])
ax.minorticks_on()
plt.scatter(
    data['RA_ICRS'], data['DE_ICRS'], 10, lw=.2, edgecolor='w', zorder=0,
    label='N={}'.format(len(data['RA_ICRS'])))
plt.scatter(
    data['RA_ICRS'][msk], data['DE_ICRS'][msk], 10, c='r', zorder=5,
    label='N={}'.format(len(data['DE_ICRS'][msk])))
plt.legend()
plt.xlabel('RA')
plt.ylabel('DE')
ax.set_aspect(aspect='auto')
plt.gca().invert_xaxis()

ax = plt.subplot(gs[2:4, 2:4])
ax.minorticks_on()
ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
plt.scatter(
    data['BV'], data['Gmag'], 10, c='b', edgecolor='w', lw=.2,
    zorder=0)
plt.scatter(
    data['BV'][msk], data['Gmag'][msk], 10, edgecolor='w', c='r', lw=.2,
    zorder=5)
plt.gca().invert_yaxis()
med, std = np.nanmedian(data['BV']), np.nanstd(data['BV'])
plt.xlim(med - 5. * std, med + 5. * std)
plt.xlabel('BV')
plt.ylabel('G')

ax = plt.subplot(gs[2:4, 4:6])
ax.minorticks_on()
ax.grid(b=True, which='major', color='gray', linestyle='--', lw=.5)
plt.scatter(data['Plx'], data['Gmag'], marker='^', c='grey', s=5)
plt.scatter(
    data['Plx'][msk], data['Gmag'][msk], marker='o', c='r', s=10)
med = np.average(data['Plx'][msk], weights=1. / data['e_Plx'][msk])
std = np.nanstd(data['Plx'][msk])
t1 = r'$Plx_{{wa}}$={:.3f}$\pm${:.3f}'.format(med, std)
d_pc = 1000. / med
t2 = r"$\mu={:.2f}$".format(-5 + 5 * np.log10(d_pc))
txt = t1 + '\n' + t2 + '\nd={:.0f} pc'.format(d_pc)
plt.axvline(med, c='k', label=txt)
plt.axvline(med - std, c='b', ls=':')
plt.axvline(med + std, c='b', ls=':')

med, std = np.nanmedian(data['Plx']), np.nanstd(data['Plx'])
plt.xlim(med - std, med + std)
plt.xlabel('Plx')
plt.ylabel('gray')
plt.legend()
plt.gca().invert_yaxis()

fig.tight_layout()
plt.savefig(cluster + '_explore.png', dpi=200, bbox_inches='tight')
