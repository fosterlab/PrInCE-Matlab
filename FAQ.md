# FAQ

__Where are my interactions?__

All interactions are in __Final_Interactions_XX_precision.csv__ ("XX" is the chosen precision of the interaction list, set in experimental_design.rtf). This is in the folder __Output/Data/Interactions/__. This folder also has other interaction tables, such as Final_Interactions_list_condition1only_50_precision.csv, which is the interactions only predicted in condition1. Interaction figures are in the folder Output/Figures/Interactions - particularly interesting figures are precision-vs-number-of-interactions (Final_precision_XX.png) and elution curves for the 200 best interactions (in the BestInteractions folder).

__Where are my complexes?__

Complexes built from all predicted interactions are listed in __Complexes_All_interactions.csv__. This is in the __Output/Data/Complexes/__ folder. Figures of complexes are in the Output/Figures/Complexes/ folder.

__Where are my fold changes (between conditions)?__

Fold changes are listed in __Final_2foldchange_list_between_condition1andcondition2.csv__ in folder __Output/Data/FoldChanges/__. (If there are more than two conditions, "condition1andcondition2" might be e.g. "condition1andcondition3".) Figures for every individual protein are in Output/Figures/FoldChanges/IndividualProteins. If there are no FoldChange tables or figures, it could be because you only have one condition. If there are no figures in the IndividualProteins folder, you probably have "Plots of fold changes (FoldChanges): No" in experimental_design.rtf.

__What are all these output files?__

See [TablesAndFigures](TablesAndFigures.md).
