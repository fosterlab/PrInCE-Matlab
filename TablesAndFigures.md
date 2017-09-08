# Tables

__Interactions__ (Output/Data/Interactions/)

* Final_Interactions_list_50_precision.csv. Pairwise protein interactions predicted from all replicates and conditions. The "Interaction score (avg.)" column is very useful, since it lets you create interaction lists of any precision. For example, to generate an interaction list with 75% precision, take all interactions with an interaction score greater than 0.75. For 90% precision, take all interactions with a score greater than 0.90. "Interaction probability" is the output of the Naive Bayes classifier. A true positive interaction has "1" in the "Interaction in CORUM" column; a false positive interaction has "0" in the "Interaction in CORUM" column and "1" in the "Both proteins in CORUM" column.

* Final_Interactions_list_condition1only_50_precision.csv. A subset of Final_Interactions_list_50_precision.csv, with only the interactions found in condition1. Each condition will have its own file.

* Precision_across_channels.csv. Summary of interactions found in any condition ("at least 1 condition") or more conditions (e.g. "at least 2 condition"). There should be fewer interactions found in more conditions, but with higher precision.

* Summary_Interactions_50pc_replicate.csv. Summary of the true positive and false positive interactions in Final_Interactions_list_50_precision.csv.


__Complexes__ (Output/Data/Complexes/)

* Complexes_All_interactions.csv. Complexes predicted from interactions in Final_Interactions_list_50_precision.csv, and their best matching known (CORUM) complexes. Predicted complexes with less than two members in common with a known complex are "Novel". "Complex density" is the fraction of possible intra-complex connections listed in Final_Interactions_list_50_precision.csv. "N overlap" is the number of members shared between predicted and known complexes; "CORUM coverage" is "N overlap" divided by the size of the "Best CORUM match".

* More_corum_matches_All_interactions.csv. Some predicted complexes share members with more than one known (CORUM) complex. The best match, as determined by the number of members in common, is listed here along with the second best match. The best match is the same as "Best CORUM match" in "Complexes_All_interactions.csv".

* Summary_complexes.csv. Summary of predicted complexes, including the number, average size, density, matching ratio and geometric accuracy. The latter two measure how closely the predicted complexes resemble known complexes.


__FoldChanges__

* Final_2foldchange_list_between_condition1andcondition2. Changes in protein amount between condition1 and condition2, expressed in log2 units. A fold change ("Fold change (raw data)") is given for each fitted Gaussian, which is assumed to be a complex. For example, a single protein that participates in two complexes can increase for the first complex and decrease for the second; this protein would be listed twice, with the elution fraction of each complex shwon in the "Center" column. If there was not enough data to calculate a fold change (1 or fewer data points within three fractions of "Center"), the fold change is "unquantifiable". Fold changes >1 are an "increase", fold changes <-1 are a "decrease", and all others are "no change".

* Gaussian_list.csv. List of all Gaussians fitted to co-fractionation profiles, along with the Gaussian parameters ("Center", "Width", and "Height"). "adjrsquared" is the R^2 between all Gaussians fit to a co-fractionation profile and the profile.

* Summary_gaussian_trend_analysis_protein_replicate.csv. Summary of fold changes that are a consistent increase, decrease, or do not change change, where "consistent" means all Gaussians (i.e. complexes) for a protein have the same behaviour (increase, decrease, or no change).

* Summary_gaussians_detected_between_replicates.csv. Summary of all Gaussians fit to co-fractionation profiles, and whether an increase (log2(fold change) >1), decrease (log2(fold change) <-1), or no change was detected.




# Figures

__Interactions__

* Final_Precision_50.png. A variation of a P-R curve, with precision (aka "interaction score") on the y-axis, and number of interactions on the x-axis.

* Number_interactions_per_channel.png. Top: The number of true positive, false positive, and novel interactions in at least 1, 2, etc conditions. Percentages above bars show precision. For example, x-axis=2 shows the number of interactions predicted in at least 2 conditions. Bottom: Number of interactions in each condition. For example, x-axis=2 shows the number of interactions (true positive, false positive, novel) predicted in condition 2.

* BestInteractions. This folder has plots of co-fractionation profiles for the top 200 interactions. For now, all replicates are plotted in the same figure (this is likely to change!).


__Complexes__

* Hairball2_predicted_vs_corum.png and .svg. Force-directed graphs for all predicted interactions, along with the best matching known (e.g. CORUM) complexes corresponding to each predicted complex. Orange, novel complex members, i.e. members of the predicted complex that are not in the best matching known complex; purple, predicted complex members that are also members of the best matching known complex; black, members of the known complex that are not in the predicted complex. For example, a graph with 1 orange, 2 purple, and 2 black dots is a three-protein predicted complex (1 orange + 2 purple) that partially overlaps with a known four-protein complex (2 purple + 2 black).

* Predicted_vs_corum_All_interactions_Complex_1.png and .svg. Force-directed graphs for the first predicted complex. The complex number (e.g. "1") corresponds to the "ID" column in Complexes_All_interactions.csv. For example, the figure of predicted complex 14 ("ID"=4) is Predicted_vs_corum_All_interactions_Complex_14. The same colour scheme as Hairball2_predicted_vs_corum is used.


__FoldChanges__

* Fold_changes_Gaussians.png. Log2 fold changes, measured using raw data within 2 fractions of each fitted Gaussian.

* Histogram_changes_per_fraction.png. Fold changes detected in each fraction, separated into increasing (log2(fold change)>1), decreasing (<-1), and no change. Fold changes are assigned to a fraction using the location of the fitted Gaussians.

* IndividualProteins. This folder includes a plot of each protein's co-fractionation curves, one plot per replicate. For example, 3_P62753_foldchange.png shows P62753 in replicate 3. Top: raw data. Middle: Gaussians fitted to raw data. Bottom: log2-fold changes, calculated using raw data within 2 fraction of the center of each Gaussian. Gaussians within 2 fractions are merged, and a fold change is calculated using data within 2 fractions of their average center. No fold change is calculated (grey bar) if there is less than 2 data points around the Gaussian center. x-axis numbers correspond to 
