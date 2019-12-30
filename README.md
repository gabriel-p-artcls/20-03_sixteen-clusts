
# Sixteen clusters analysis

Contents of folders outside of repo:

* `1_old_analysis`: contains the results of the 2017-2019 original analysis.
  * `original_phot`: originally received photometry for the remaining 15 clusters
* `2_bh73_phot`: photometric analysis for this cluster.

Contents of folders inside of repo:

* `0_obs_clusters`: images and code to generate the three figures of the
  original observed V frames.
* `2_bh73_phot`: `11_short_exp` contains the final photometry used for this
  cluster. More details in the `README.md` file in the `2_bh73_phot` folder
  not synced with the repo.
* `3_new_match`: results of the new Gaia DR2 match obtained using my scripts
  `astrometry` and `catalog_match`.
* `4_transformations`: Carrasco transformations for all the clusters.
* `5_ASteCa`: full analysis with ASteCA
* `6_analysis`: parallax analysis and some scripts to generate plots


## New analysis (October 2019)

Re-processed all the clusters with the latest version of ASteCA (0.2.7). This fixes the issue with outliers in PMs producing weird values in membership probabilities, and greatly improves the isochrone averaging and the IMF sampling.

#### Membership analysis

The membership analysis was re-done with the now cross-matched data, using the combination Plx + PMs (no photometry). The results are stored in the
`C1_mps_analysis/` folder. All max uncertainties were fixed to:
eG_max=0.01, eVI_max=0.1, eBV_max=0.1, eUB_max=0.2.

#### Fundamental parameters analysis

##### First run

We use GBVU photometry to restrict the E_BV parameter (particularly through the U filter). Results are stored in the `1st_run/` folder.

The clusters vdBH73 and RUP88 do not have enough stars in their photometric diagrams for this analysis to be of use. For these clusters, the second run does not use the information gathered here.

##### Second run run

The second run uses the restricted range in E_BV found above (stored in the
`E_BV_range_2nd_run.ods` file), fixes the metallicity to (approx) solar value
(0.015), and fits the remaining parameters. Results are stored in the
`2nd_run/` folder.

##### Third run

Ranges for all parameters except binarity (fixed to 0.3) and metallicity (set to the entire range), are restricted as shown in the
`3rd_input_param_ranges.ods` file.

##### Fourth run

Run again for vdBH85 and TR13.

* vdBH85: shows a **very** large discrepancy in its distance estimate between the photometric (ASteCA) distance of ~4600 pc (13.32 mag) and the Gaia DR2 distance of ~8200 pc (14.57 mag). I attempted a new run, but the issue is with the parallaxes, not with the dm so I discarded it.

* TR13: the previous run chose a very young age (log(age)=7), which results in a rather large distance (13.954 mag). Since this cluster has been studied before, we run it here again now restricting the log(age) to 8-9. The resulting best fit choses a log(age)~8, and a distance of 13.41 mag, much closer to Gaia DR2's estimate.

These results are merged into the 3rd run folder.







## New analysis (June 2019)

All clusters were re-analyzed with ASteCA. Results are stored in the
`/3_June_2019/` folder.

#### Membership analysis

The membership analysis was re-done with the now cross-matched data, using the combination G + UBVI + Plx + PMs. The results are stored in the
`C1_mps_analysis/` folder.

The vdBH73 cluster differs from the rest in the sense that the max uncertainty cut for this cluster was done at e_BV=0.2, e_UB=0.2, instead of 0.1. This is because this cluster has larger uncertainties due to the short exposure photometry used.

#### Fundamental parameters analysis

##### First run

We use GBVU photometry to restrict the E_BV parameter (particularly through the U filter). Results are stored in the `1st_run/` folder.

The clusters vdBH73, vdBH85, vdBH106, NGC4230, and RUP88 do not have enough stars in their photometric diagrams for this analysis to be of use. For these clusters, the second run does not use the information gathered here.

Special cases:

* **RUP162**: the max E_BV was restricted to `0.2<=E_BV<=0.6` (range found through the UB-BV diagram's sequence) and `log(age)<8` to avoid non-reasonable solutions.

##### Second run

The second run uses the restricted range in E_BV found above (stored in the
`2nd_E_BV_ranges.ods` file), fixes the metallicity to (approx) solar value
(0.0155), and fits the remaining parameters. Results are stored in the
`2nd_run/` folder.

Special cases:

* **vdBH73, vdBH85, vdBH106, NGC4230, and RUP88**: the Schlafly and Finkbeiner (2011) `E_BV` maximum values are used.
* **RUP162**: the age was again restricted to  `log(age)<8`, due to the large contamination present in the diagram. This age range means the synthetic sequence in the VI vs BV CCD is reasonably covered.

##### Third run

Ranges for all parameters except binarity (fixed to 0.3) and metallicity (set to the entire range), are restricted as shown in the
`3rd_input_param_ranges.ods` file.




## New analysis (May 2019)

Solving the mentioned issues in the previous analysis.

### vdBH73

To solve the issue with vdBH73, I used the photometry obtained with the short exposure times (since the problems arise from mixing short and long exposures taken at very different times, and the standard frames are for the short exposures). This results in a shallower photometry with larger uncertainties, but at least this way the cluster can be matched with GaiA DR2 with reasonable results. The files used to generate the final data to be analyzed by ASteCA are stored in the `2_bh73_phot/11_short_exp/` folder.

### Gaia DR2 match

Initially I wrote a `clean.py` script (in the folder `0_original_phot/gaiadr2_match/`) to remove those low-mass stars with bad matches. The percentage of stars removed by this process with a maximum magnitude difference cut of 1 mag are:

| Cluster  | % lost |
| -------- | ------:|
| rup85    |  27.2  |
| tr13     |  26.5  |
| rup162   |  23.5  |
| tr12     |  22.9  |
| loden565 |  17.4  |
| rup88    |  17.2  |
| bh91     |  15.5  |
| rup87    |  11.7  |
| bh106    |  11.2  |
| lynga15  |  11.2  |
| ngc4230  |  10.8  |
| bh85     |   9.9  |
| bh87     |   9.8  |
| bh92     |   9.1  |
| n4349    |   8.0  |

I realized that I should try to re-match with Gaia (the original matching came with the photometry, probably done by Edgard or Giovanni) to see if I could avoid losing that many stars.

The match with Gaia DR2 was re-done from scratch using my scripts `astrometry` and `catalog_match`. I used a maximum match radius of 2-4 arcsec, depending on the cluster. The results are much better than with the (ra,dec) coordinates that came with the photometry. The percentage of stars lost (not matched) are:

| Cluster  | max d | max V | % lost |
| -------- |:-----:|:-----:| ------:|
| rup162   | 2.5'' |   1   |    7   |
| lynga15  |  4''  |   1   |    6   |
| tr12     |  4''  |   1   |    6   |
| tr13     |  2''  |   1   |    4   |
| ngc4230  |  3''  |   1   |    4   |
| bh92     |  3''  |   1   |    4   |
| rup85    |  2''  |   1   |    2   |
| rup87    |  2''  |   1   |    2   |
| rup88    |  2''  |   1   |    2   |
| bh106    |  2''  |   1   |    2   |
| bh85     |  2''  |   1   |    2   |
| bh87     |  2''  |   1   |    2   |
| loden565 |  2''  |   1   |    2   |
| bh91     |  2''  |   1   |    1   |
| n4349    |  2''  |   1   |    1   |

The results are stored in the `3_new_match/` folder.

### Transformations with Gaia DR2

Using the data with the new cross-match, I re-did the transformations now using a cut on 0.05 over our photometry (eV, eBV, eVI). The results are a bit better than before, now all clusters have median differences with the Carrasco transformations below 0.06.

The images are stored in the `4_transformations/` folder.



## Old analysis (2017 - 2019)

These clusters were initially analyzed in 2017/2018. The results were used to generate a first draft of the article in March 2019, that was criticized by André.

The results, scripts, and outputs from this first run are stored in the
(not synced) `1_old_analysis/` folder.

There were two main problems with this first analysis:

1. The photometry for vdBH73 shows a very poor match with Gaia DR2.
2. Most of the clusters show bad matches with Gaia DR2 for a group of low mass stars. This throws off the uncertainty analysis.
3. The Carrasco transformations showed poor values for some clusters.
4. The kinematic data from Gaia DR2 (Plx+PMs) was not used in the membership analysis.







