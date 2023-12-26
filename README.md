# ORCA

Python scripts that predict and plot the location of the origin of replication (*oriC*) of circular bacterial genomes based on Z-curve, GC-skew, *dnaA*-box, and gene location analyses. This README will not explain all of ORCA's methods. All functions, main or helper, are labelled with docs-strings and type-hinting for ease-of-use. Most functions in ORCA can be used separately as well. For example, `calculate_disparities_curves` can be used outside of ORCA, too. As well as all functions for fetching files from NCBI, or plotter functions, etc. Please, see their provided docs-strings for more info.

## ORCA class

This script predicts the origin of replication for circular bacterial DNA. It makes use of a combination of [Z-curve](https://en.wikipedia.org/wiki/Z_curve) and [GC-skew](https://en.wikipedia.org/wiki/GC_skew), *dnaA*-box, and gene location analyses. You can load the required FASTA files yourself, or simply provide an accession and NCBI-account email and ORCA will fetch them. The docs-string of the function shows more information on what is needed to use ORCA. See the code for the function `example_use` for more information. **NOTE:** the model provided in this repository has been compressed due to file size limitations. This example also uses the [Joblib](https://joblib.readthedocs.io/en/stable/) package for unpickling the model. This is not necessary, and Python's standard pickle package can be used.

Please, make sure to load the proper file formats into the functions, otherwise ORCA will not work. A lot of invalid parameters will throw warnings or errors, but it is not improbable that a few were missed. More of ORCA's functionality and input methods are discussed in the paper.

### Example use of ORCA

```python
>>> import joblib
>>> from ORCA import ORCA
>>> email = "example@email.com"
>>> model = joblib.load("data/output/machine_learning/24k_set_model.pkl")
>>> orca = ORCA.from_accession("NC_000913.3", email=email, model=model)
>>> orca.find_oriCs(show_info=True, show_plot=False)
Accession    : NC_000913
predictions  :  0.99412
Z_scores     :      1.0
G_scores     :      0.5
D_scores     :      1.0
oriC_middles :  3927479
The best-scoring potential oriC was found at: 3927479 bp.
```

In the case of *Escherichia coli* (*E. coli*) K-12, only one potential origin was found and according to the model, this could also be classified as a true origin. Any prediction value >= 0.5, means that the model believes there is a 50 % or more that the corresponding origin is a true origin. It is possible that multiple candidate origins have a probability of being a true origin that is larger than 50 %. Then simply use the origin corresponding to the highest probability, as this was the origin that the model deemed most likely to be correct.

This repository also includes a pickled `SeqRecord` of the *E. coli* K-12 chromosome. If one wanted to use that, only a different constructor would have to be used.

```python
>>> import joblib
>>> from ORCA import ORCA
>>> model = joblib.load("data/output/machine_learning/24k_set_model.pkl")
>>> orca = ORCA.from_pkl("data/input/NC_000913_3.pkl", model=model)
>>> orca.find_oriCs(show_info=False, show_plot=False)
```

## ORCA's parameters

The parameters listed below are parameters that can be used as they are or can be fine-tuned for specific use cases. The standard parameters reflect ORCA's performance as shown in the application note. Retraining of the Random Forest Classifier is recommended if any parameters are changed. Otherwise, the same performance as is the paper can not be guaranteed.

- `dnaa_boxes`: If None, will use the consensus DnaA-box: `TTAT(A|C|G|T)CACA` (See Section: [DnaA](#dnaa)). Else, provide a list of 9 base strings. Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
- `max_mismatches`: Maximum allowed mismatches allowed in a dnaa_box for it still to be read as such. Recommended max: 2. ORCA uses 0 for use with the consensus DnaA-box.
- `genes_of_interest`: List of gene names to consider as 'oriC-proximal' and use for helping estimate the location of the oriC. This parameter is case insensitive.
- `max_point_spread`: Maximum distance between points in a group can have when looking for connected groups across the disparity curves. Default is 5 % of the total chromosome length.
- `windows`: The windows around peaks of skew curves to consider. Defaults are 1, 3, and 5 % of the total chromosome length. ORCA checks each of the given windows.
- `model`: A fitted scikit-learn classifier. Recommended to use the one provided on in this repository.

## Libraries and Versions used

- [Python 3.11.1](https://www.python.org/)
- [SciPy 1.10.0](https://scipy.org/)
- [Pandas 1.5.3](https://pandas.pydata.org/)
- [NumPy 1.24.2](https://numpy.org/)
- [Biopython 1.80](https://biopython.org/)
- [Scikit-learn 1.2.1](https://scikit-learn.org/)
- [Matplotlib 3.6.3](https://matplotlib.org/)

## Peak class

The `Peak` class. This class is used in handling potential *oriC*s. The `oriCs`-attribute of an `ORCA` object consists of a list of `Peak` objects. Each `Peak` represents a potential *oriC* and has attributes for its Z-, G-, and D-score as well as the confidence a potential machine learning model has.

## BioFile

Script with useful functions for I/O. These include functions for downloading, parsing, and saving files.

## Plotter

There are 3 general functions in this file which can be used to plot any generic 1D-np.array. To use these functions, make sure to have [matplotlib](https://matplotlib.org/) installed

- `plot_Z_curve_3D`: Makes a 3D-plot of the Z-curve.
- `plot_curves`: Can plot a maximum of four curves in a single 2D-plot. This one is useful for plotting single or multiple Z-curve or GC-skew components against each other.
- `plot_skew`: Does the same as `plot_curves`, except only takes one array.

## Machine_learning

This folder contains a lot of data used in analysing, training, and testing ORCA and its Random Forest Classifier. The DoriC data can be downloaded from: <http://tubic.tju.edu.cn/doric/public/index.php> on May 16th, 2023. Each script comes with docs-string explaining their functionality.

## data

### input

This folder contains all the input data. This includes the data from DoriC 12.0. It also includes the experimental_dataset CSV. This is a collection of oriCs that have been experimentally verified. It also includes sources for each chromosome. This dataset was made in order to check the quality of DoriC.

This folder also includes a pickled `SeqRecord` object of the *E. coli* K-12 chromosome. This file was used for quick testing and demonstration.

### output

This folder contains output from the various performance test we ran on ORCA. The machine_learning sub-folder contains the provided Random Forest Classifier used in the application note.

## DnaA

Useful articles about DnaA(-boxes): [[1]](https://doi.org/10.1101/cshperspect.a012922), [[2]](https://doi.org/10.1093/nar/gkr832), [[3]](https://doi.org/10.3389/fmicb.2018.00319), [[4]](https://doi.org/10.1046/j.1365-2958.1996.6481362.x).

We use the first entry in the table below and allow for 0 mismatches in this sequence.

|                                Sequence                          |                    Paper                       | Year | Notes                            |
| ---------------------------------------------------------------- | ---------------------------------------------- | ---- | -------------------------------- |
| TTAT(A\|C\|G\|T)CACA                                             | [[5]](https://doi.org/10.1093/dnares/dsm017), [[6]](https://doi.org/10.1093/bib/bbn031) | 2007/8 | In (roughly) all eubacteria, not just *B. subtilis* or *E. coli* |
| TTATCCACA, TTTTCCACA, TTATCAACA, TCATTCACA, TTATACACA, TTATCCAAA | [[7]](https://doi.org/10.1093/emboj/16.21.6574)      | 1997 | Affinity study                   |
| (T\|C)(T\|C)(A\|T\|C)T(A\|C)C(A\|G)(A\|C\|T)(A\|C)               | [[8]](https://doi.org/10.1007/BF00273584)            | 1991 | Only in E. coli K12. Do not use. |
| TGTG(G\|T)ATAAC                                                  | [[9]](https://doi.org/10.1016/0022-2836(85)90299-2)  | 1985 | Matsui-box                       |
| TTAT(A\|C)CA(A\|C)A                                              | [[10]](https://doi.org/10.1016/0092-8674(84)90284-8) | 1984 | The first consensus              |
