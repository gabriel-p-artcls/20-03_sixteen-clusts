
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


converters = {'ID': [ascii.convert_numpy(np.str)]}

cl = 'bh85'
# Read MP data
mp_data = ascii.read(cl + '_match_memb.dat', converters=converters)
# Read matched data.
BJ_dist_pc = ascii.read(cl + '_match_match.dat', converters=converters)

BJ_lst = list(BJ_dist_pc['ID'])

d_pc, e_d_pc, e_plx, Gmag = [], [], [], []
for star in mp_data:
    if star['sel'] == 1:
        try:
            idx = BJ_lst.index(star['ID'])
            # mps.append(star['MP'])
            d_pc.append(BJ_dist_pc['rest'][idx])
            # Error as average of upper and lower uncertainty estimates
            e_d_pc.append(
                .5 * (BJ_dist_pc['B_rest'][idx] - BJ_dist_pc['b_rest'][idx]))
            e_plx.append(BJ_dist_pc['e_Plx'][idx])
            Gmag.append(BJ_dist_pc['Gmag'][idx])
        except:
            pass

plt.scatter(d_pc, e_d_pc)
plt.show()
print(np.mean(d_pc), np.std(d_pc))
wa = np.average(d_pc, weights=e_d_pc)
# https://stackoverflow.com/a/2415343/1391441
wa_stddev = np.sqrt(np.average((d_pc - wa)**2, weights=e_d_pc))
print(wa - wa_stddev, wa, wa + wa_stddev)
# Another approach for the STDDEV with slightly different result
print(np.sqrt(np.cov(d_pc, aweights=e_d_pc)))

print("Mean Plx error for selected stars:", np.nanmean(e_plx))
# 75% of estimated members have e_plx > 0.1 mas
print("25% percentile ePlx for selected stars:", np.percentile(e_plx, 25))
# 61% of stars have Gmag>18
print("39% percentile Gmag for selected stars:", np.percentile(Gmag, 39))
