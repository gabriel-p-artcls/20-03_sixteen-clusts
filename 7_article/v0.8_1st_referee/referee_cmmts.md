--------------------
Referee Report

This paper analyses the photometric data of 16 open cluster candidates in the direction of the Perseus and Carina Milky Way arms. The colour-magnitude diagrams and colour-colour diagrams are analysed combined with proper motions and parallaxes from Gaia to determine the physical parameters of the clusters. There is good data, a good software, and a good plan, nevertheless the analysis needs to be done again and the writing of the text would benefit from a thorough internal consistency check, more formal and objective language.

Please find below some major and minor points:


MAJOR POINTS:

Section 3.1, par.4
"Table 4 shows the mean...":
Your point is clear here, however the mean differences may mask some colour dependency on the magnitude, for example. A more accurate comparison would be to use APASS stars as standard stars and then compare the coefficients with those obtained from your own calibration. If you do not detect any systematics, only adding a sentence saying so would suffice.


Section 4:
"we should see an increase in the star density (not always true)": if this is a limitation of the method, I wonder how sure you can be to confirm that a candidate is indeed not a cluster. Please clarify.



Fig.2
- "Right": There is a clear trend in this plot. The differences are very small, maybe of the order of your photometric uncertainties, but this needs to become clear. Otherwise you would need to check for a colour-dependence on the magnitude for example. This seems a detail, but in this paper you discuss systematics on distances obtained with Gaia and ground-based photometry, therefore I recommend verify this trend.



Section 4.2, par.9:
- "The only information...": previously, when the Bayesian approach with 7 dimensions was mentioned, it gave the impression that parameters were set free (presumably with some prior, limiting intervals etc), although not explicitly said. Now, the description says that some parameters are kept fixed in order to determine others. Please discuss the impact on the final parameters, uncertainties and degeneracy if all 7 dimensions are kept free during the fit. If this was done in previous papers, then just a summary and references is fine.

Section 5.0, par.2:
- footnote 12: vdBH85, RUP85 seems to have a round-ish shape and a reasonable radial density profile. I understand this statement for some clusters, but it is defintely not a general case. Please remove it and leave the discussion to the case by case subsections.
- "show as well the isochrone curves": do you generate isochrones for the exact age and metallicity final values? Or you choose a pre-selected isochrone with approximate values? If the latter, then the age and metallicity differences between the final results and the displayed isochrone will affect the visualisation. Please clarify.


Fig.5:
- "estimate the field star properties": Stars outside the larger black square have densities ranging from less than 20 to about 70 stars/arcmin2. Your final estimate for the field density was 16 stars/arcmin2. Please explain how you found this very low value. It does not seem to represent the field density. Also: I find it at least suspicious that the borders of the image have the lowest background. It happens in all cases, for large and compact clusters. I wonder if this is vignetting from the observations or some other artifact not considered. Please verify this and comment on the paper.
- "star number per square min": right panel shows the RDP from 80 stars/arcmin2 down to 16 within 5arcmin. But the left panels shows that the large square has about 11x11 arcmin2 and a minimum density of 35 stars/arcmin2. It does not make sense. Please explain and correct.


Fig.9:
- again, as in Fig.5, left and right panels are not consistent with each other. The smaller square on the left panel has 7.5x7.5 arcmin size and the lowest density of about 30 stars/arcmin. In the right panel the density drops from 40 to 20 in less than 2 arcmin. Please correct in all figures and revise all analysis.


Fig.10:
- well-defined CMD and CCDs. However, there is a 93% probable star at Gmag ~10.3mag and some main-sequence turnoff stars that are not fitted by the isochrone. I wonder if a slightly younger isochrone would be a better fit. Please revise this analysis.


Section 6:
- this whole analysis assumes that the photometric distances are accurate. However, as pointed out before in the analysis of the three selected sample clusters, the extinction and isochrone fitting could be improved. Therefore it seems that the photometric distances may (or may not) suffer from biases in the analysis. As a consequence, this analysis could not be used as a reference to define the parallax offset to Gaia results. Although the values agree with previous studies, this is not an argument to prove the quality of the photometric analysis. I strongly recommend revise the analysis with the specifics detailed before.
































MINOR POINTS:


Abstract, results:
"young and therefor potential arm tracers": this is highlighted as a good result. Nevertheless, it would be good to connect the results with the aims presented at the beginning of the abstract and hightlight other relevant results.

Introduction, par.3:
- "combining UBVI" --> "combining ground-based UBVI"
- "space based astrometry" --> "space-based astrometry"


Introduction, par.4:
- "absolute coordinates" --> "equatorial coordinates" (change it everywhere it appears in the paper)

Section 2:
- please explain how these clusters were selected to be analysed here. There are many clusters in this quadrant, do these 16 clusters represent the quadrant somehow, or are they just a subsample?


Section 3:
par.1:
- YALO: define acronym

par4:
- "All observations were made under air mass ranges from 1.08 to 1.092 for U, 1.072 for B and 1.077 to 1.164 for V and I."
- "are indicated using the usual vdBH abbreviation." --> "are indicated by vdBH." --> This information could be moved and completed to Tables 2 and 3 for all clusters and observations.


Table 1:
- According to the A&A style, the position of the legend should be on the top of tables, instead of below as in the case of figures. Please check all tables.
- "Note: van den Bergh-Hagen clusters": although it is well-known, please write the full reference, presumably van den Bergh & Hagen (1975, vdBH). Same for the other acronym definitions. Same for the other clusters and acronyms.

Table 2-3:
Tables 2 and 3 could be merged. There is repeated information about exposure times that could be summarized. YALO was also at CTIO. Suggestion: one (or more) line per cluster following Table 1, with columns: Date, airmass, DIMM seeing (or image quality), filters, exposure times.

Section 3.1, par.2
- after the equation the sentece starting with "where" should be preceded by "\noindent" latex command
- sentence starting with "In each case detector..." should be a new paragraph, as the topic is WCS calibration and astrometry. Please provide more details on this step: was there any distortion correction? What is the precision of the match? This information is useful to endorse your results based on the match between your data and that from Gaia.

Section 3.1, par.3
- "All-Sky Survey2)." --> "All-Sky Survey2), that has a magnitude limit of XXX, enough to identify the XXX [upper RGB?] stars in our sample clusters."
- "are for the most part very faint." --> "are mostly very faint."



Section 3.1, par.5
"our photometry" --> "our photometry with that from Gaia DR2"
"display the median differences": Why table 4 shows the mean and table 5 the median? Please clarify your choices.

Section 3.1, par.6
"Coordinates should not be read from images.": Please clarify why not.



Fig.2:
- Figures read well, but the colour bar between left and central panels may add some confusion. I suggest to move the colour bar to the right of the middle panel instead to make it clear that tha y-axis in both panels are exactly the same.



Section 4.1, par.2
- "following the advice": please add a few words to mention the advice given.
- "the model for the cluster": please specify which type of model your are talking about: 3D spatial distribution? colour-magnitude diagram...?
- the sentence "Our model marginalizes...decontaminated cluster region" is a bit confusing. Maybe because the definition of the model is missing above. Please consider rephrasing after defining what is the model.
- "maximum likelihood estimate": please clarify how this is done and what is compared to obtain the likelihood for the prior.
- "shown in a Plx vs G plot for each cluster in Sect 5." --> "be shown in Sect 5."
- "We also show": please add a reference to where this is shown. Presumably on Section 5? If so, please leave the details to Section 5 with the proper reference to the respective figures.


Section 4.1, par.3
"are indicated": please add reference to where these results are indicated.



Section 4.2:
- I agree it is important to provide a brief summary on the main points of ASteCA and refer the interested reader to the original papers for more details. Please highlight here any improvements or adaptation on the original code/paper to include more dimensions (e.g. N=7, rather than 4) on the analysis from Gaia (parallax, proper motions).


Section 4.2, par.3:
- "square rings": It is not clear why a square was chosen instead of a circle. This could be the reason why you have outliers, they may come from the corners of the squares. Even though the details of ASteCA are in previous papers, this particular choice needs a sentence to justify.
- ", i.e.:" --> ", i.e.,"


Section 4.2, par.7:
- isochrones shown: please add reference to where this is shown.


Section 4.2, par.8:
- "the effects of star loss at large magnitudes": do you mean photometric incompleteness or actual low-mass runaway stars from the clusters? Presumably the former; please rephrase for clarification.



Section 4.2, par.10:
- "was" --> "is"
- "finally compared to around ~2x10^7 synthetic clusters.": please clarify what this step means. Is a stochastic variation of star positions on the parameters space to derive uncertainties? Are you actually varying parameters from the initial isochrones...?



Section 5:
par 1, 4 are OK, but par 2, 3 are confusing to read before presenting the actual results. My suggestion is to dillute this detailed description (par 2,3) into Section 5.1 in the example of the first cluster, to be more objective.


Section 5, par 2:
- "four figures": please add explicit reference to which figures you are referring to.
- " "clean" color-color ": please remove quotes and be more speciifc on what you mean by clean CCD.
- "color-color diagram" --> "CCD"
- "In these three panels": it is very confusing to imagine all these figures with this level of details without actually looking at the figures. Please add references to which specific figures you are describing here.



Section 5, par 3:
- "is shown": please specify where this is shown.
- "fourth figure": please describe the so-called "fourth figure" in the same paragraph (either the previous or a new one, but avoid splitting the description).



Section 5, par 4:
- "of two extreme types of cluster" --> "of three extreme types of cluster"


Section 5.1, par 1:
- the sentence "However, the inspection ... presence of a cluster there." is misleading. From fig.3 it is clear that there a significantly larger fraction of field stars. If the CMD and CCD are not decontaminated in any sense, naturally you will not see the tiny cluster there. Specially a relatively older cluster. Please rephrase.
- "larger magnitudes the CMDs strongly widen surely due to the presence of increasing visual absorption.": visual absorption does change with the magnitude of the stars... The wider colour distribution for fainter stars is due to photometric errors. Please correct the argument in this paragraph or simply erase it.
- "Even some blue stars": this sentence seems to pass the message of a surprise that blue stars are specially affected by reddening. I do not see why this is a surprise. Consider rephrase the sentence.


Fig.4:
- I do not see the point of this figure. The final fit is done with Gaia filters as presented in Fig.6. Showing the polluted CMD and CCD of all stars together without any further info and different filters does not add much to the paper (that already has too many figures). My suggestion is to show CCDs and CMDs that may contain the field stars as background dots and the cluster member stars as foreground points. And display only the combination of filters and colours used to fit the parameters discussed in the paper.


Fig.5:
- "decimal format": Please be consistent: use the same units/format for RA, DEC in all your figures. Fig.3 displays RA in decimal hours, fig. 5 displays RA in decimal degrees, for example. Update all figures.
- "color scale": please also add label in the plot.
- "Right panel": what is the magnitude limit of the stars used to draw the RDP? This is important information as it directly affects the fitted parameters.
- "vertical black line": do you account for photometric incompleteness? It makes a lot of difference in the final structural parameters. Please comment about it in the paper.


Section 5.1, par 3:
- "0.38": this is a very specific number. Usually people use a cut in 50%. Please clarify how this number was defined.
- "well detached": please be consistent: this cluster was one of three examples selected because it was "poorly defined". Please correct and maybe reconsider the selection of this cluster as the "poorly defined" example.
- "low chance to be confused with the stellar field.": please be consistent: this cluster was listed on Sect.5 as "easy to confuse with the background"


Section 5.1, item b:
- please provide uncertainties for distance. Fig.6 shows error bars, it is unclear why they are neglected here.


Section 5.1, second last par:
- please provide uncertainties for Gaia distance.
- "Anderson-Darling test": please discuss the actual results from the test, for example, say explicitly the p-value in the text.


Section 5.1, last par:
- please provide uncertainties for the age. Fig.6 shows error bars, it is unclear why they are neglected here.


Fig.6:
- Please be consistent: on fig.5 you write B-V and there Bmag-Vmag. The same for other filters. Correct in all figures.
- "(B - V) vs (U - B)": the isochrone doe not fit well the points. Please verify the process. A wrong extinction directly affects the age, distance, metallicity determinations.
- "color bar": add label to the color bar "memb. prob." or something to identify its meaning.
- "Insert" --> "Inset"


Fig.7:
- "see text": The text that refers to this Fig.7 is Sect.5.1, and it does not explain these fits. Please define them briefly here in the figure caption to ease the reading of the figure. Further discussions on the implications of the different fits can be left to the text.


Section 5.2, par.1:
- "Fig.8": same comment as in sect.5.1 for fig.4. NGC4349 is younger and therefore the upper main sequence stands-out, but it does not mean much in Fig.8 if there is no information on membership probability, or at least distance from the cluster center. Consider remove or merge fig.8 with fig.10. similarly for all clusters as suggested before.



Section 5.2, par.2:
- "(say above 0.8": avoid informal language writting.


Section 5.2, item a:
- "one concludes that": if this conclusion is true, then you would see the cluster well detached from the background field stars in the CCD, but it does not seem to be the case. It could be the window used by S&F2011, or the fit. Please verify.


Section 5.2, par. 3:
- "marginally lower": please write down the Gaia distance explicitly. Avoid subjective comparisons.


Section 5.2, par. 5:
- "around ~0.51x10^9 years old": please add uncertainties
- "reporting log(t) = 8:32 equivalent to ~ 0.21x10^9 yrs." --> "reporting log(t) = 8:32 equivalent to ~ 0.21x10^9 yrs."


Fig.11:
- the fits seem to be at lower Plx with respect to the distribution of redder points. Please check.


Section 5.3, par. 1:
"boring": avoid subjective adjectives.
"respective color color maganitude resemble": correct and rephrase


Section 5.3, par. 2:
- "obviously seen in Fig.13": ... and strong hints in figs. 14, 15 by eye looking at the redder points with higher membership probability. I wonder if this is an actual cluster and the code was not able to find it for some unknown reason for me. Please check again the analysis.
- "we decide to focus the attention in the region encircled by a green line.": it is unclear what it was done here. ASteCA did not look around the circle and squares defined in fig.13? Please clarify.
- "significant portion of those may also belong to the surrounding field.": true. Also true: a significant portion belongs to the cluster. Please rethink the analysis selecting the stars with a larger membership probability.


Fig.12:
- same comments as for Fig.4 and all similar others.


Fig.14:
- the isochrones are significantly off the redder point distribution. Does ASteCA consider the membership probability as weight to the fit? Please check.


Fig.15:
- the fitted Plx are significantly off the well defined sequence of redder points. Based on figs.3, 14, 15, it could be a low-mass older cluster, but the photometry is not deep enough to see it. Please check the analysis again.



Sect.6, par.2:
- "~2195 pc vs. 2115 pc": please add uncertainties, the same presented as error bars in fig.16.
- "means that its sequence is not only": presumably main sequence?


Sect.6, par.4:
- the sentence "If we add to the parallax data ... worsens, as shown in Fig. 17." is redundant. Fig.16 and the obtained offset was 0.025. Therefore the previous paragraph together with Fig.16 are clear enough, no need to make another figure. Consider removing fig.17.
- "In the case of vdBH85 we use the individual distance estimates obtained in Bailer-Jones et al. (2018).": it is unclear why the analysis of Bailer-Jones was used. They basically perform a similar analysis as it is done here, but applying an offset of 0.029mas to the Gaia parallaxes. Here you perform this analysis and find an offset of 0.025mas. You may only reduce this paragraph to a sentence saying that using a similar analysis with a similar offset, Bailer-Jones found similar distances as you.


Sect.6, par.5:
- "The mean...": attach this paragraph to the end of the paragraph starting with "A number of..."


Sect.6, par.6:
- the sentence "This analysis points to ... given by Lindegren et al." is redundant. It was already stated when first discussing fig.16. Consider removing.



Sect. 7, par1:
- "inconvenient": presumably inconvenience?



Sect. 7, par3:
- "billions" --> "billion years"
- "case of vdBH 85": if true, this result should be highlighted. It is very interesting finding to have an open cluster this old.


Fig.17:
- this figure does not add much to the discussion. Fig.16 plus the discussion in the text is enough.


Fig.18:
- draw also the Perseus arm as it was mentioned in Section 3, for reference.
