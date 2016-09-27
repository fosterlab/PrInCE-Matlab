# PCP-SILAC
Bioinformatics pipeline for analyzing PCP-SILAC and other co-elution experiments.

* Predict protein-protein interactions
* Predict protein complexes
* Calculate protein abundance differences

Original versions of the analysis were written by Anders Kristensen<sup>1</sup> and Nichollas Scott<sup>2</sup>.


# How to use

### 1. Download pipeline.

In a browser, go to to [https://github.com/fosterlab/PCP-SILAC](https://github.com/fosterlab/PCP-SILAC). Click the green "Clone or Download" button in the top right. Then click "Download ZIP".

![Download pipeline from github](/ReadmeFigures/01download.jpg?raw=true)

### 2. Format your data.

#### Data files
This pipeline is designed to work on multiple csv files. Each data file contains data from a single biological condition and one or more replicates. For example, data from a PCP-SILAC experiment with M/L and H/L ratios would be formatted into two csv files. **Important: Ensure that files are "saved as csv" in whatever program you use, e.g. Excel.**

Columns 1 and 2 are *protein ID* and *replicate number*, respectively. *protein ID* must be a string that matches to the reference database file (see below). *replicate number* must be an integer. Subsequent columns are protein amounts (e.g. isotopologue ratio for PCP-SILAC experiments) from fractionated data. For example, an experiment with 50 fractions would have files with 52 columns.

The first row in each file is reserved for the header. Each subsequent row is a chromatogram from a protein (or protein group) from a replicate and condition. An example data file for an experiment with two replicates and four protein IDs is below. Note: Real data should have many more proteins and fractions!

![Format your data files like this](/ReadmeFigures/examplefile1.jpg?raw=true)

#### Reference database of known complexes
This pipeline needs a reference database of known interactions, e.g. CORUM. This reference database must be a csv file in the same format as CORUM's *allComplexes.csv* file (downloadable [here](http://mips.helmholtz-muenchen.de/genre/proj/corum/)). That is, the reference database file must:
* have a header
* be semicolon-separated
* each reference complex must be on a separate line in the fourth column

Note: protein IDs in data files must match a subset of protein IDs in the reference database.


### 3. Organize your files.

Create a folder for storing pipeline output (e.g. "PCPanalysis/"). Unpack the zipped folder downloaded in step 1. From this unzipped directory, copy the Code/ folder and pcpsilac.m into the PCPanalysis/ folder. In PCPanalysis/, create a folder called "Input". Place all data files (one for each condition) and the reference database file in the PCPanalysis/Input/ folder. The resulting file structure should look like this:

  * PCPanalysis/
    * Code/
    * Input/
      * Reference database file (e.g *allComplexes.csv*)
      * Data file 1 (e.g *condition1.csv*)
      * Data file 2 (e.g *condition2.csv*)
      * Data file 3 (e.g *condition1.csv*)
      * ...
    * pcpsilac.m



### 4. Enter the details of your experiment in pcpsilac.m.


### 5. Run analysis.
Open Matlab. In Matlab, *cd* to the PCPanalysis/ folder. In the command line type

```
pcpsilac.m
```


## References

1. Kristensen AR, Gpsoner J, Foster LJ. A high-throughput approach for measuring temporal changes in the interactome. Nat Methods. 2012 Sep;9(9):907-9. doi: 10.1038/nmeth.2131.
2. Scott NE, Brown LM, Kristensen AR, Foster LJ. Development of a computational framework for the analysis of protein correlation profiling and spatial proteomics experiments. J Proteomics. 2015 Apr 6;118:112-29. doi: 10.1016/j.jprot.2014.10.024.
