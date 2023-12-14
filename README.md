# ORCA
Python scripts that predict and plot the location of the origin of replication (*oriC*) of circular bacterial genomes based on Z-curve, GC-skew, *dnaA*-box, and gene location analyses. This README will not explain all of ORCA's methods. All functions, main or helper, are labelled with docs-strings and type-hinting for ease-of-use. Most functions in ORCA can be used seperately as well. For example, `calculate_disparities_curves` can be used outside of ORCA, too. As well as all functions for fetching files from NCBI, or plotter functions, etc. Please, see their provided docs-strings for more info.

## ORCA class
This script predicts the origin of replication for circular bacterial DNA. It makes use of a combination of [Z-curve](https://en.wikipedia.org/wiki/Z_curve) and [GC-skew](https://en.wikipedia.org/wiki/GC_skew), *dnaA*-box, and gene location analyses. You can load the required FASTA files yourself, or simply provide an accession and NCBI-account email and ORCA will fetch them. The docs-string of the function shows more information on what is needed to use ORCA. See the code fot the function `example_use` for more information.

Please, make sure to load the proper file formats into the functions, otherwise ORCA will not work. A lot of invalid parameters will throw warnings or errors, but it is not improbable that a few were missed. More of ORCA's functionality and input methods are discussed in the paper.

**Example use of ORCA**
```python
>>> import joblib, ORCA
>>> email = 'example@email.com'
>>> model = joblib.load('Machine_Learning/75_train_model.pkl')
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

In the case of *Escherichia coli* K-12, only one potential origin was found and according to the model, this could also be classified as a true origin. Any prediction value >= 0.5, means that the model believes there is a 50 % or more that the corresponding origin is a true origin. It is possible that multiple candidate origins have a probability of being a true origin that is larger than 50 %. Then simply use the origin corresponding to the highest probability, as this was the origin that the model deemed most likely to be correct.

**Libraries and Versions used:**
- [Python 3.11.1](https://www.python.org/)
  - [SciPy 1.10.0](https://scipy.org/)
  - [Pandas 1.5.3](https://pandas.pydata.org/)
  - [NumPy 1.24.2](https://numpy.org/)
  - [Biopython 1.80](https://biopython.org/)
  - [Scikit-learn 1.2.1](https://scikit-learn.org/)
  - [Matplotlib 3.6.3](https://matplotlib.org/)


### Peak class
The `Peak` class. This class is used in handling potential *oriC*s. The `oriCs`-attribute of an `ORCA` object consists of a list of `Peak` objects. Each `Peak` represents a potential *oriC* and has attributes for its Z-, G-, and D-score as well as the confidence a potential machine learing model has.


### BioFile
Script with a lot of useful functions used for I/O.


### Plotter
There are 3 general functions in this file which can be used to plot any generic 1D-np.array. To use these functions, make sure to have [matplotlib](https://matplotlib.org/) installed
- `plot_Z_curve_3D`: Makes a 3D-plot of the Z-curve.
- `plot_curves`: Can plot a maximum of four curves in a single 2D-plot. This one is useful for plotting single or multiple Z-curve or GC-skew components agaist each other.
- `plot_skew`: Does the same as `plot_curves`, except only takes one array.


### Machine_learning
This folder contains a lot of data used in analysing, training, and testing ORCA and its Random Forest Classifier. The DoriC data can be downloaded from: http://tubic.tju.edu.cn/doric/public/index.php. Each sceipt comes with docs-string explaining their functionality.


### data

#### input

This folder contains all the input data. This includes the data from DoriC 12.0. It also includes the experimental_dataset CSV. This is a collection of oriCs that have been experimentally verified. It also includes sources for each chromosome. This dataset was made in order to check the quality of DoriC.

This folder also includes a pickled `SeqRecord` object of the *E. coli* K-12 chromosome. This file was used for quick testing and demostration.

#### output

This folder contains output from the various performance test we ran on ORCA.


### DoriC_processing
- `data_prep_doric.py`: Processes DoriC CSV into something more usable.
- `download_files.py`: This script was used to download the necessary FASTA-files for ORCA to use when predicting all accessions in DoriC.
