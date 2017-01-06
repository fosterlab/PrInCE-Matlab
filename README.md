# PRInCE

Pipeline for Recovering Interactomes via Co-Elution. Bioinformatics pipeline for analyzing PCP-SILAC and other co-elution experiments.

* Predict protein-protein interactions
* Predict protein complexes
* Calculate protein abundance differences

Original versions of the analysis were written by Anders Kristensen<sup>1</sup> and Nichollas Scott<sup>2</sup>.


# How to use

1. Format your data.
2. Download pipeline.
3. Put your data in the input folder..
4. Enter your experiment details in prince.m.
5. Run prince.m.

### 1. Format your data.

#### Data files
You will be making one csv file for each biological condition. Each row, except for the header, is co-fractionation data from a single protein and replicate, e.g. a chromatogram. Each file is formatted like this:

* Column 1: Protein ID (must match with reference)
* Column 2: Replicate number (integer)
* Columns 3-end: Co-fractionated protein amount, e.g. isotopologue ratio for PCP-SILAC

Here's a simple dataset with one condition, two replicates, and five fractions:

![Format your data files like this](/ReadmeFigures/examplefile1.jpg?raw=true)

Important: Ensure that files are "saved as csv", e.g. Excel --> "Save as" --> "Save as csv"

#### Reference database of known complexes
This pipeline needs a reference database of known protein complexes, e.g. CORUM. For now, this reference database must be a file in the same format as CORUM's *allComplexes.csv* file (downloadable [here](http://mips.helmholtz-muenchen.de/genre/proj/corum/)), and we recommend this for mammalian datasets. In case a custom reference must be made, it must follow this format:

* must be a csv file
* must have a header
* reference complexes are in the fourth column, and each member is semicolon-separated

Note: protein IDs in data files must match a subset of protein IDs in the reference database.


### 2. Download pipeline.

In a browser, go to to [https://github.com/fosterlab/PRInCE](https://github.com/fosterlab/PRInCE). Click the green "Clone or Download" button in the top right. Then click "Download ZIP".

![Download pipeline from github](/ReadmeFigures/01download.jpg?raw=true)


### 3. Unzip the pipeline, and put your data in the Input folder.

In a working directory, unzip the downloaded zip file (PrInCE-master.zip). This creates four folders and a few files inside the working directory, like this:

  * a working directory/
    * Code/
    * Input/
      * Reference database file (e.g *allComplexes.csv*)
      * Data file 1 (e.g *condition1.csv*)
      * Data file 2 (e.g *condition2.csv*)
      * ...
    * Output/
    * prince.m

Place all data files (there will be one for each biological condition) and the reference database file (typically allcomplexes.csv from CORUM) in the Input/ folder.

### 4. Enter the details of your experiment in prince.m.


### 5. Run analysis.
Open Matlab. In Matlab, *cd* to the folder that contains pcpsilac.m, i.e. the working directory you made in 3. In the command line type

```
prince.m
```

As the pipeline runs, output figures and tables will be deposited in the Output/ folder.


## References

1. Kristensen AR, Gpsoner J, Foster LJ. A high-throughput approach for measuring temporal changes in the interactome. Nat Methods. 2012 Sep;9(9):907-9. doi: 10.1038/nmeth.2131.
2. Scott NE, Brown LM, Kristensen AR, Foster LJ. Development of a computational framework for the analysis of protein correlation profiling and spatial proteomics experiments. J Proteomics. 2015 Apr 6;118:112-29. doi: 10.1016/j.jprot.2014.10.024.
