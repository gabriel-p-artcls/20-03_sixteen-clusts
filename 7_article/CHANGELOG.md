
# Article

*Sixteen overlooked open clusters*, Perren et al. (2019)

## Versions

#### v1.3-AA - 2020/04/03

* Version uploaded to SAGA


#### v1.2-AA - 2020/04/02

* Version with language editing and new figures


#### v1.1-AA - 2020/03/27

* Version sent to A&A
* Slightly modified to follow US convention
* Removed comments from tex


#### v1.0-arxiv - 2020/03/26

* Version sent to the arXiv


#### v0.9 - 2020/03/18

* Second set of corrections by the referee.


#### v0.8 - 2020/02/09

* Corrections by the referee.


#### v0.7 - 2019/11/17

* Corrections by André
 
- [x] Fix significant figures in distances in Table 6 and figures 7, 11, 15,..
- [x] Title: "Sixteen overlooked open clusters in the fourth Galactic quadrant: A combined analysis of UBVI photometry and Gaia DR2 with ASteCA"
- [x] Mention outlier detection in Section 4.2 about the membership algorithm
- [x] Add Andre's corrections

* Version uploaded to A&A


#### v0.6 - 2019/10/18

- Included corrections by André 02/10 (in sixteen_clusts_v3_AM.pdf)

##### André corrections

0. Abstract

- [x] Una de las cosas importantes era darle un giro al abstract para hacer los objetivos más claros. Sugiero algo como este:
(hay un texto entre "**" que es contenido nuevo y que tenemos que discutir)

Aims.
This paper has three main objectives: (1) To search for tracers of the Perseus arm and of the Galactic warp in the poorly known sector, at 270  < l < 300 in the fourth quadrant of the Milky Way. (2) To determine the distances and ages of sixteen faint and almost unknown catalogued open clusters, for use as probes of the structures under investigation. (3) To address previously reported systematics in Gaia DR2 parallaxes by comparing the cluster distances derived from photometry with those derived from parallaxes.

Methods. A deep UBVI photometric survey of the 16 cluster fields was performed. The photometric observations were reduced and analyzed in an automatic way using the ASteCA package, resulting in an homogeneous set of distances, reddening, masses, ages and metallicities. Distances obtained from the photometric analysis were compared to those obtained from a Bayesian analysis of Gaia DR2 parallaxes.

Results. Seven out of the sixteen catalogued clusters are confirmed as true entities, while three are considered as probable ones. Three of these ten objects are young and therefore potential arm tracers. The existence of clusters in the remining fields is discarded, being their patterns due to random fluctuations in the stellar distribution. The comparison of distances from photometry and parallaxes data revealed a variable level of disagreement.

Conclusions.
We have detected two objects that appear related to the Carina feature. No trace of Perseus has been found. (** Beyond 4 kpc, we find one young object, Trumpler 13 appearing to trace the warp at ~200 pc below the plane, and two older clusters more than 300 pc above the plane (~500 if we consider the plane to be warped down). **explain these two or raise the question.***)

Various zero point corrections for Gaia DR2 parallax data from recent literature were considered in our comparison between photometric and parralax based distances. We find that results tend to improve with these corrections. Our photometric distances suggest a correction of ~ 0.07 mas (to be added to the parrallaxes). The correction may have a more intricate distance dependency, but addressing that level of detail will require larger sample.

- [x] In several cases, the (asteca) analysis appears to be dominated by faint stars G> 17.5, where there is already considerable spread in the CMD.  It is worrisome that the high probability stars are almost always at the fainter not-so-reliable end of the CMD. This seems to happen when there is a big concentration of stars in the lower CMD.
> This has been addressed by the new membership processing.

- [x] It seems that the memberships in Version 1 of the paper were better (illusion of colour scale?) Membership colour scale should be kept the same.
> The corrected membership processing gives much better results. The color scale is there to make the more probable members pop-out. If we normalize this to a 0-1 scale for all clusters, most of the detail is lost.

- [x] BH73 - Based on faint stars. Parameters do not seem very reliable (2.64 kpc). Compare with plx based isochrone?
> The new processing fixes this discrepancy.

- [x] NGC 4349 - Based on faint stars. Seems that distance erros are very large because of that (2.75 kpc), although general fit in CMD seems OK (but with low probability members).
> The new processing results on a much clearer clean sequence and much better parallax-photometric distance fit.

- [x] BH87 - Good probabilities but fit looks awkward
> Much better probabilities and fit in the new results.

- [x] TR12 -  weird. Holes in CMD
> Fixed by the new membership processing.

- [x] Check BH106 - d_asteca too large? 10^10 yr?
> Both parameters are now smaller. The distance is a better match to the parallax distance, and the age is now ~3^9 yr. Still, this is a dubious cluster as is RUP162.

- [x] RUP162 - High probability stars everywhere in the CMD. Problems with probability assignment. Or is it the colour scale?
* This is probably because there is no cluster here. The membership algorithm works under the assumption that a cluster exists in the region. If there is no cluster, the algorithm still tries to assign large probabilities to stars in the cluster region that don't "look like" stars in the surrounding field region. The case of this cluster/region is complicated because there does seem to be an overdensity and its photometric distribution does look somewhat different to that of the field regions, but there also does not seem to be a clear cluster sequence. There is also no clear overdensity in either the parallax or the proper motions. The AD test also supports the hypothesis that there is no cluster here. In any case, I suggest keeping this as a "dubious" cluster.

3.1. Photometric reduction process
- [x] "Not very good.. What is the quoted error (possible systematics) for the Carrasco transformations?"
> I think it actually is very good. The Carrasco transformations mention a sigma of 0.0467 for the G-V vs V-I diagram, and 0.06285 for the G-V vs B-V. Our mean differences for those transformations are 0.0244 and 0.0341.

4. Photometric data analysis process: Gaia data and the ASteCA code
- [x] Is it over distances or over parallaxes? Just to make sure, the final distance for the cluster is the inverse of the global parallax (and not the average of the individual stellar distances), correct?
> The model uses the individual star's parallax values (which are known) and *marginalizes* the individual star distances (which are unknown).Ie: we integrates those parameters out of the problem because we do not care about the individual distances to the stars (only about the single distance to the cluster).
The final distance for the cluster is neither the inverse of the global parallax nor the average of the individual stellar distances, but the estimated 'r_c' parameter in the model that can be seen in Eq 20 here: https://github.com/agabrown/astrometry-inference-tutorials/blob/master/single-source/tutorial/resources/cluster_inference.pdf

7. Discussion of results and concluding remarks

- [x] circular argument? The models have been calibrated using z_solar = 0.0153, so re-finding it just means that most clusters have solar metallicity, not that the value is 0.0153. Forv that we would need spectroscopy.
> There appears to be some confusion here. We are *not* stating that we can estimate the solar metallicity value (which is indeed 0.0152 for the PARSEC isochrones), but that the fitting process tells us that our set of clusters has on average solar metallicity, and that this is a reasonable results given the region being explored.
Just to clarify: the models are only fitted using a fixed z value in the first two runs, were we look for suitable ranges for the remaining parameters. The third and final run fits the 'z' parameter using a large range, as a completely free parameter.
Is this clearer now? Perhaps André would like to suggest a re-wording of this paragraph?
- [x] what else can we say? Plot also other studied clusters and see them ours in context? It’s OK to be incremental here. Comment not only spiral structure, but general structure: warp (done below), but also clusters high above plane.
> This paragraph was modified and this is no longer one of our aims.
- [x] Perhaps a comment here: what range of distances?
> This sentence was removed as now the agreement is very good overall.

#### Changes in clusters' parameters

##### vdBH73
* the distance is now larger (2.6 kpc --> 5 kpc) and so is the excess (0.6 --> 1)
* Gaia parallax distance now coincides with photometric
* Cluster is now classified as 'intermediate age' instead of old (6x10^9 yrs --> 7.8x10^8 yrs)

#### NGC4349
* distance changed sligthly (2.7 kpc --> 2.2 kpc), now better agreement with Gaia
* age is a bit larger (4x10^8 --> 5.1x10^8 yrs)

#### RUP85
* MPs changed so I removed the sentence: "Intermediate probability members are found around the giant branch but we suggest that care should be taken with this fact."

#### vdBH85
* MPs changed so I removed the sentence: "Also at $G=14$ mag and almost located at the center of the CMDs, a couple of stars with MPs around 0.8 could be classified as blue stragglers associated to van den Bergh-Hagen 85." as those stars are no longer there

#### vdBH87
* MPs changed so I removed the sentence: "The CCD no the other hand does not show a clear cluster sequence." as it now does

#### TR12
This cluster improved markedly with the new MPs. Somewhat large changes to the text were required to clearly assert this clusters existence.
* I removed the sentence: "We draw the attention to the fact that only 142 stars were detected inside the potential cluster zone."
* I removed the sentence: "There exists a clear difficulty to separate field stars from cluster members, reflected in the noisy distribution visible in the CMDs."
* I removed the sentence: "in case TR 12 is a real cluster"
* I removed the paragraph: "The noisy profile that does not allow establishing a clear cluster extension and the poor main sequence requires us to be cautious. More deep photometry ($(U-B)$ particularly) is needed to achieve a precise result as for the true nature of this object."

#### vdBH91
* I removed the sentence: "It is interesting that most stars with lower probabilities outline a fictitious sequence below $G=16$ mag in these diagrams"

#### vdBH 106
* No longer the oldest cluster in the sample, vdBH85 is now.

#### RUP 162
The memberships assigned to this cluster changed considerably, so the text had to be adapted. I believe there is no cluster here, but we can keep it as "dubious" anyway.

#### Loden 565
* No clear sequence is present now. Changed the text to reflect this.


#### v0.5 - 2019/07/30

- Included corrections by André 26/07 (in sixteen_clusts_v2_AM.pdf)

##### André corrections

1. Intro
- [x] "to uncover its fundamental parameters" -> removed

3. Photometric observations
- [x] "copy-paste from another article..." -> modified the paragraph

3.1. Photometric reduction process
- [x] "have no optical photometric studies" -> "have no dedicated photometric studies"
- [x] "not so good for bluer (brighter?) stars" -> the comparison between
Gaia G mag and our transformed G mag is good. On average it is -0.01986,
with min=-0.0532, max=0.0255. Not worse for brighter stars than for low mass stars. Added average median differences to the text, and G_gaia vs G_Transf data.

4.1. Gaia data
- [x] "reference" -> the reference is Luri et al. (2018), this is an accompanying tutorial as stated
- [x] "This is a fundamental problem!" -> The part about the negative values was removed, and text was added to better explain this process.
- [x] "in ASteCA?" -> added "in our analysis". The current version of ASteCA includes a more sophisticated version of this test, but the analysis in the article is previous to the one in ASteCA.

4.2. The way ASteCA works
- [x] "Specify" -> more data given
- [x] "e.g." -> the referenced article is not an example, this is where the actual likelihood used by ASteCA is described
- [x] "Not clear." -> re-structured the paragraph
- [x] "e.g." -> added

5. Cluster-by-cluster discussion on structural and intrinsic parameters provided by ASteCA
- [x] "Isn't this what the A-D test is supposed to solve?" -> The A-D test is an independent method to estimate the probability that the cluster region and the field regions come from the same parent distribution. As with any statistical test, it is not definitive. What ASteCA does is basically the minimization of a function (the likelihood) through a mathematical optimizer
(the genetic algorithm). Such a set up will *always* return a "best fit", and it is up to the researcher to combine and asses all the evidence (i.e.: combine the results of the A-D test, the membership analysis, and the fundamental parameters analysis) into a single conclusion for each cluster.
- [x] "There are so few true clusters that they should be spelled out." -> mentioned in the text
- [x] "anomalous, electronic" -> removed

5.1. van den Bergh-Hagen 73
- [x] "~4% is not huge.." -> the *difference* is 2.7 Kpc, that's almost a 50% difference in the distance estimates
- [x] "not good.." (Fig 6) -> it's not the same figure
- [x] "Plx ASteCA is clearly off..." (Fig 7) -> as mentioned in the article, the photometric data for vdBH73 shows only ~1 mag of its main sequence. This means the photometric analysis is necessarily hard and error-prone. Still, we do not want to restrain the distance range in the photometric distance analysis using the information given by the parallax distance, as this would bias the photometric distance result, and then it would not be a fair comparison between photometric and parallax distances.

5.2. NGC 4349
- [x] "all high prob in the lower MS?" (Fig 10) -> This cluster shows an overdensity in its low-mass star region, compared to the field regions. The algorithm thus expectedly identifies stars in this region as very probable cluster members.

6. Gaia parallax distances analysis
- [x] "Briefly remind the reasons evoked for these offsets" -> added comment
- [x] "vdBH 73 and van den Bergh-Hagen 106" -> corrected
- [x] "a more extended sample" -> analysis of a more extended sample"

7. Discussion of results and concluding remarks
- [x] Added Valle (2005) missing reference


##### Not clear about the correction

5.1. van den Bergh-Hagen 73
- [ ] "significant digits" -> more, less?

##### Highlighted but no comment given:

1. Intro
* "define, strongly, reasonably, truly, very large"
3.1. Photometric reduction process
* "astrometry"
4. Photometric data analysis process: Gaia data and the ASteCA code
* "the optimal"
4.2. The way ASteCA works
* "may return incorrect values"
* "Hence, it should be regarded as a lower limit on the actual initial mass value"
5. Cluster-by-cluster discussion on structural and intrinsic parameters provided by ASteCA
* "others are clusters which are faint"
5.2. NGC 4349
* "slightly"
5.3. Ruprecht 87
* "poor and boring"


#### v0.4 - 2019/07/13

- Included modifications by Ruben (in sixteen_clusts_rev.pdf)
- Added Carrasco Gaia transformations
- Shortened abstract
- Lots of other minor changes


#### v0.3 - 2019/07/05

At the request of André, the analysis was re-done using Gaia DR2 data for the membership analysis, and PARSEC v1.2s isochrones.

Sent to Rubén for corrections.

##### Things that were revised/modified in this version

André corrections:

0. Abstract
- [x] GAIA --> Gaia (entire article)
- [x] Changed to proper A&A formating, added some keywords

1. Introduction
- [x] Added Moitinho (2010), Moitinho et al .(2006)
- [x] Added Bossini et al. (2019), Soubiran et al. 2018
- [ ] "Firstly, we wanted to study Galactic structure in a poorly known Galactic sector" --> should this be our **primary** purpose? We say nothing about the Galactic structure.
- [ ] "a bit more specific" 
- [ ] "Not clear on what this statement is based on"
- [ ] "This is a bit surprising"
- [x] "to make our analysis more robust..." --> sentence was modified
- [x] "On the other side, for the most distant clusters ..." --> sentence was removed

3. Photometric observations
- [ ] Table 3 --> shouldn't this table include the names of the clusters that correspond to each observing run?

3.1. Photometric reduction process
- [x] "cross-correlated" --> cross-matched
- [ ] "ooops… this is really bad photometry.." --> We can either: 1. keep APASS section and add Gaia DR2 transformations data, 2. Remove APASS entirely and replace with the Gaia DR2 transformations data.

4. Photometric data analysis process: the ASteCA code. The use of GAIA data.
- [x] Modified the name of the section and moved the 'Gaia data' sub-section before the 'The way ASteCA works' sub-section for clarity.

4.1. GAIA data
- [x] "This is one of the major problems with the paper." --> The sub-section was modified to accommodate this changes.
- [x] "Anderson- Darling quick explanation" --> added more info
- [x] "We followed the advice given in Lindegren..." --> edited this paragraph to include one more recent article where an intermediate offset is proposed

4.2. The way ASteCA works
- [x] "why not use Gaia proper motions at this stage?", ASteCA is not equipped to do this analysis (yet). It would be interesting, but it would force us to re-do (again) the entire processing of the data (after I code this into ASteCA, which would also take some time) It's a good idea, but I don't think it is worth the trouble right now.
- [x] "Not really needed" --> removed sentence.
- [x] "hundredths" --> hundreds
- [x] "by user request" --> defined by the user"
- [x] "ASteCA does not fit isochrones to observed stellar clusters." --> "ASteCA does not fit isochrones to cluster sequences in photometric diagrams"
- [x] "synthetic color-magnitude" --> "synthetic clusters in color-magnitude"
- [x] modified the description of the Bayesian membership algorithm to reflect that we also use Plx and PMs data
- [x] updated isochrone set used Marigo --> PARSEC
- [x] updated the section with more details about the analysis.

5. Cluster-by-cluster discussion on structural and intrinsic parameters provided by ASteCA
- [x] "by" --> with
- [x] "clean --> field decontaminated"
- [ ] Revise all the clusters' descriptions since the new parameter values are now different from the previous values

5.1. van den Bergh-Hagen 73
- [x] "free-absorption" --> absorption-free
- [x] "2.911" --> 2.9
- [x] "2511x10^6" --> 2.5x10^9
- [ ] "few stars above V=15 mag" --> not true anymore?
- [ ] "stars with negative (U-B) values appear strongly affected by variable reddening" --> not true anymore?
- [x] "subtending 3 magnitudes" --> subtending 1.5 magnitudes
- [x] "and a giant branch with stars up to $V = 15.5$ mag" --> and a faint giant branch with stars up to $G=15$ mag
- [x] "clean CCD in the middle panel shows the locus occupied by stars in vdBH 73, completely unnoticed in the equivalent CCD of Fig. 3" --> clean
CCD in the middle panel shows only a handful of stars with no clear locus defined"
- [x] "2.9 kpc from the Sun" --> 2.6 kpc from the Sun
- [x] "difference in distance reaching up 1.5 kpc" --> 2.7 kpc

5.2. NGC 4349
- [x] "appears in the three diagrams a narrow cluster" --> appears in the three diagrams a somewhat narrow cluster
- [x] "is remarkably close to the photometric distance" --> is lower than the photometric distance (..) and gets lower as larger offsets are applied
- [ ] "High probability values for stars inside the overdensity confirm the
true nature of this object since the over density and the density profile
are followed by a very well defined and extended photometric counterpart." <-- this sentence is hard to follow

6. Gaia parallax distances analysis
- [x] Added this separate section

7. Discussion of results and concluding remarks
- [x] "determined the true nature" --> analyzed the fields
- [x] "in precise $UBVI$ photometry carried out" --> on precise UBVI photometry analyzed
- [x] "semi-automatic" --> automatic (I don't think fixing center and radii is enough to call the procedure "semi")
- [x] "Perhaps it’s distance dependent?" --> there does not seem to be enough evidence to either support or discard this
- [x] "Again, it might depend on distance or some other parameter" --> it might, but we can not conclude anything with certainty
- [x] "circular argument since the isochrones are based on that?" --> no they are not, the metallicity is fitted in the final run. The new values for the mean and stddev metallicity are updated to: 0.01537, 0.009
- [x] "drastically" --> removed
- [ ] "but not among the older ones" <-- reference?
- [x] "The origin of these differences may be due to offsets in parallaxes in the regions where they are immersed" <-- removed
- [x] Table 5 --> parameters updated. Instead of using the Lindegren distances in d_Bayes, I put here the parallax distances estimated with no added bias
- [x] The cluster van den Bergh-Hagen 106 is said to be dubious in Appendix I, but it was not marked with an asterisk in Table 5, so I added the asterisk to the table entry
- [ ] This entire section was heavily edited, so it needs to be carefully checked

* Changes to figures
- [x] Fig 5, 9, 13, ...; updated CMD plots for the 16 clusters
- [x] Fig 6, 10, 14, ...; updated plx no offset analysis for the 16 clusters
- [x] Fig 15; update with new distances from ASteCA
- [x] Fig 16; update with new distances from ASteCA
- [ ] Fig 17; update with new distances from ASteCA
- [ ] Fig 17; add missing reference

* Appendix A: Ruprecht 85
- [x] updated parameter values
* Appendix B: vdBH 85
- [x] "Three magnitudes above the cluster turn-off.." --> changed this sentence since the new MPs show a different scenario.
* Appendix C: vdBH 87
- [x] updated parameter values and text
* Appendix D: Trumpler 12
- [x] updated parameter values andtext
* Appendix E: vdBH 91
- [x] updated parameter values and text
* Appendix F: Ruprecht 88
- [x] updated parameter values and text
* Appendix G: vdBH 92
- [x] updated parameter values and text
* Appendix H: Trumpler 13
- [x] updated parameter values andtext
* Appendix I: vdBH 106
- [x] updated parameter values and text
* Appendix J: Ruprecht 162
- [x] updated parameter values and text, mainly the AD test changed
* Appendix K: Lynga 15
- [x] updated parameter values and text
- [ ] Reference missing
* Appendix L: Loden 565
- [x] updated parameter values and text
* Appendix L: NGC 4230
- [x] updated parameter values and text. Added Tadross reference


#### v0.2 - 2019/04/03

* Corrections by Ruben added. Still a couple of corrections in the appendix are missing.
* Sent to Giovanni, Ruben, and André 2nd round of corrections.


#### v0.1 - 2019/03/19

* Initial version of the manuscript, written using A&A's LaTeX template.
* Sent to Ruben for 1st round of corrections.
