# ORCA

Python scripts that predict and plot the location of the origin of replication (*oriC*) of circular bacterial genomes based on Z-curve, GC-skew, *dnaA*-box, and gene location analyses. This README will not explain all of ORCA's methods. All functions, main or helper, are labelled with docs-strings and type-hinting for ease of use. Most functions in ORCA can be used separately as well. For example, `calculate_disparity_curves` can be used outside of ORCA, too. As well as all functions for fetching files from NCBI, or plotter functions, etc. Please, see their provided docs-strings for more info.

## Installing ORCApy

To install ORCApy, simply use:

```sh
pip install orcapy
```

Or download it directly from the PyPi website [here](https://pypi.org/project/orcapy/). The installation and download only include the files present in the `src` folder.

To download the provided Random Forest Classifier, simply use:

```sh
wget https://github.com/ZoyavanMeel/ORCA/blob/main/data/output/machine_learning/ORCA_RFC_model_1_4_0.pkl.gz
```

Or use a similar method, like `curl`. You can also download it manually by going [here](https://github.com/ZoyavanMeel/ORCA/blob/main/data/output/machine_learning/ORCA_RFC_model_1_4_0.pkl.gz), clicking on the three horizontal dots in the top right and selecting 'Download'. See the section [Data Output](#output) for more provided models.

## ORCA class

This script predicts the origin of replication for circular bacterial DNA. It makes use of a combination of [Z-curve](https://en.wikipedia.org/wiki/Z_curve) and [GC-skew](https://en.wikipedia.org/wiki/GC_skew), *dnaA*-box, and gene location analyses. You can load the required FASTA files yourself, or simply provide an accession and NCBI-account email and ORCA will fetch them. The docs-string of the function shows more information on what is needed to use ORCA. See the [Example Use](#example-use-of-orca) section for more information. **NOTE:** the model provided in this repository has been compressed due to file size limitations. This example also uses the [Joblib](https://joblib.readthedocs.io/en/stable/) package for unpickling the model. This is not necessary, and Python's standard pickle and gzip package can be used as well.

### Input

Please, make sure to load the proper file formats into the functions, otherwise ORCA will not work. A lot of invalid parameters will throw warnings or errors, but it is not improbable that a few were missed. More of ORCA's functionality and input methods are discussed in the paper.

There are four valid input formats for ORCA. Provided one has an internet connection, ORCA can function using only an accession number (`from_accession`). When using only an accession, ORCA will fetch the corresponding GenBank file from NCBI using Biopython's Entrez module. Please adhere to NCBI's rules on making use of their API. This module functions as a Python interface for NCBI's Entrez Direct. The version number of the accession can be specified. If this is omitted, the most recent version is automatically chosen. This method is overall the slowest and takes about one minute to download and analyse for a chromosome of 5Mbp. One can also provide the GenBank file themselves (`from_gbk`). This file will be processed the same way, but does not require an internet connection.

Another option is to provide Biopython SeqRecord objects (`from_pkl`). ORCA uses Biopython's SeqIO module for parsing the GenBank file. This means that input from GenBank will go through this module anyway and uses SeqRecords. Providing a SeqRecord that was made by parsing a Genbank file using this same module is therefore also possible. Pickled SeqRecord objects take up more memory, but they can be parsed quicker. This is a consideration that can be taken into account when processing large genomic datasets.

Lastly, one can also simply provide a DNA-sequence string (`from_string`). This same input method allows for indicating the location of the indicator genes. These gene locations are not inferred from the provided sequence and will have to be provided. This is the fastest input format, but takes the most work from the user. All input parameters are assumed to be in the correct specified format and will not be processed using Biopython. Further explanation on each of these formats can also be found in the documentation of the code.

### Example use of ORCA

With the RFC, we provide a general use-case *oriC* prediction tool as outlined in the application note. However, as shown in the code documentation, there are many parameters that can be fine-tuned and possibly improved. One of the reasons of ORCA being open-source, is to provide not just a transparent *oriC*-prediction tool, but also a building block for further research.

It is possible to tune all parameters and train models for highly accurate prediction of the *oriC*s of bacterial species of interest. This could be done by incorporating more indicator genes for that species, or changing any number of other parameters. It could also be possible to adapt ORCA for the use in the prediction of *oriC*s in linear prokaryotic chromosomes, or archaea. it is our hope that ORCA can help research into any or all of these avenues with its adaptability or provide a lightweight easy-to-use tool with good out-of-the-box performance.

```python
>>> import joblib
>>> from orcapy import ORCA
>>> email = "example@email.com"
>>> model = joblib.load("path/to/model.pkl.gz")
>>> orca = ORCA.from_accession("NC_000913.3", email=email, model=model)
>>> orca.find_oriCs(show_info=True, plot_path=None)
Accession    : NC_000913
predictions  :  0.99412
Z_scores     :      1.0
G_scores     :      0.5
D_scores     :      1.0
oriC_middles :  3927479
The best-scoring potential oriC was found at: 3927479 bp.
```

If you do not wish to use joblib, you can also open the model using:

```python
import gzip, pickle
with gzip.open("path/to/model.pkl.gz", "rb") as fh:
    model = pickle.load(fh)
```

In the case of *Escherichia coli* (*E. coli*) K-12, only one potential origin was found and according to the model, this could also be classified as a true origin. Any prediction value >= 0.5, means that the model believes there is a 50 % or more that the corresponding origin is a true origin. It is possible that multiple candidate origins have a probability of being a true origin that is larger than 50 %. Then simply use the origin corresponding to the highest probability, as this was the origin that the model deemed most likely to be correct.

This repository also includes a pickled `SeqRecord` of the *E. coli* K-12 chromosome ([here](https://github.com/ZoyavanMeel/ORCA/blob/main/data/input/NC_000913_3.pkl), download like the RFC model). If one wanted to use that, only a different constructor would have to be used.

```python
>>> from orcapy import ORCA
>>> orca = ORCA.from_pkl("data/input/NC_000913_3.pkl", model=model)
>>> orca.find_oriCs(show_info=False, plot_path="path/to/plot.png")
```

There are four constructors in total: `from_accession`, `from_gbk`, `from_pkl`, `from_string`. Each comes with extensive documentation of how to use them. Once the `ORCA` object has been instantiated, call `find_oriCs` to analyse the sequence.

## ORCA's parameters

The parameters listed below are parameters that can be used as they are or can be fine-tuned for specific use cases. The standard parameters reflect ORCA's performance as shown in the application note. Retraining of the Random Forest Classifier is recommended if any parameters are changed. Otherwise, the same performance as is the paper can not be guaranteed.

- `dnaa_boxes`: If None, will use the consensus DnaA-box: `TTAT(A|C|G|T)CACA` (see Section [DnaA](#dnaa)). Else, provide a list of 9 base strings. Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
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
- [Scikit-learn 1.4.1](https://scikit-learn.org/)
- [Matplotlib 3.6.3](https://matplotlib.org/)

## Peak class

The `Peak` class. This class is used in handling potential *oriC*s. The `oriCs`-attribute of an `ORCA` object consists of a list of `Peak` objects. Each `Peak` represents a potential *oriC* and has attributes for its Z-, G-, and D-score as well as the confidence a potential machine learning model has.

## BioFile

Script with useful functions for I/O. These include functions for downloading, parsing, and saving files.

- `save_gbk`: Fetches and saves the Genbank file of the given accession in the given path.
- `save_pkl`: Fetches and saves the given accession in the given path as a pickled `SeqRecord` generated from parsing a Genbank file with Biopython.
- `fetch_file`: Fetches the given accession in the provided file format from NCBI. Will return a `TextIOWrapper` to be used.
- `parse_SeqRecord`: Extracts the locations of the provided genes of interest into a dictionary.

To use these function call:

```python
from orcapy import BioFile
```

## Plotter

There are 3 general functions in this file which can be used to plot any generic 1D-np.array. To use these functions, make sure to have [matplotlib](https://matplotlib.org/) installed.

- `plot_Z_curve_3D`: Makes a 3D-plot of the Z-curve.
- `plot_curves`: Can plot a maximum of four curves in a single 2D-plot. This one is useful for plotting single or multiple Z-curve or GC-skew components against each other.
- `plot_skew`: Does the same as `plot_curves`, except only takes one array.

To use these functions call:

```python
from orcapy import Plotter
```

## Machine_learning

This folder contains a lot of data used in analysing, training, and testing ORCA and its Random Forest Classifier. The DoriC data can be downloaded from: <http://tubic.tju.edu.cn/doric/public/index.php> on May 16th, 2023. Each script comes with docs-string explaining their functionality.

## data

### input

This folder contains all the input data. This includes the data from DoriC 12.0. It also includes the `experimental_dataset` CSV. This is a collection of oriCs that have been experimentally verified. It also includes sources for each chromosome. This dataset was made in order to check the quality of DoriC.

This folder also includes a pickled `SeqRecord` object of the *E. coli* K-12 chromosome. This file was used for quick testing and demonstration.

### output

This folder contains output from the various performance test we ran on ORCA. The machine_learning sub-folder contains the provided Random Forest Classifiers. Use the gzip and pickle from the standard Python library to load the models or, alternatively, use [joblib](https://joblib.readthedocs.io/en/stable/). Please be mindful of your Scikit-learn version, as models are not always backwards-compatible.

- `ORCA_RFC_model_1_4_0.pkl.gz`: This model has been trained on the full DoriC 12.0 dataset. 1_4_0 refers to the scikit-learn version that was used to train the model (1.4.0). This repository also contains a fully trained model on scikit-learn 1.2.1.
- `ORCA_RFC_model_70_percent.pkl.gz`: This model has been trained on roughly 70 % of the DoriC 12.0 dataset. First all the experimentally verified origins were subtracted and the remaining dataset was split 70:30 stratified. This model is has been included for posterity. Use the other model for any analyses. This model was trained on scikit-learn version 1.2.1.

## DnaA

The table below shows a small overview of common *dnaA*-boxes and relevant papers associated with them. This table is by no means complete, but can serve as a starting point for further research. We also include some more papers that were useful in researching the effects of the *dnaA* protein and its binding sites: [[1]](https://doi.org/10.1101/cshperspect.a012922), [[2]](https://doi.org/10.1093/nar/gkr832), [[3]](https://doi.org/10.3389/fmicb.2018.00319), [[4]](https://doi.org/10.1046/j.1365-2958.1996.6481362.x).

We use the first entry in the table below and allow for 0 mismatches in this sequence.

|                                Sequence                          |                    Paper                       | Year | Notes                            |
| ---------------------------------------------------------------- | ---------------------------------------------- | ---- | -------------------------------- |
| TTAT(A\|C\|G\|T)CACA                                             | [[5]](https://doi.org/10.1093/dnares/dsm017), [[6]](https://doi.org/10.1093/bib/bbn031) | 2007/8 | In (roughly) all eubacteria, not just *B. subtilis* or *E. coli* |
| TTATCCACA, TTTTCCACA, TTATCAACA, TCATTCACA, TTATACACA, TTATCCAAA | [[7]](https://doi.org/10.1093/emboj/16.21.6574)      | 1997 | Affinity study                   |
| (T\|C)(T\|C)(A\|T\|C)T(A\|C)C(A\|G)(A\|C\|T)(A\|C)               | [[8]](https://doi.org/10.1007/BF00273584)            | 1991 | Only in E. coli K12. Do not use. |
| TGTG(G\|T)ATAAC                                                  | [[9]](https://doi.org/10.1016/0022-2836(85)90299-2)  | 1985 | Matsui-box                       |
| TTAT(A\|C)CA(A\|C)A                                              | [[10]](https://doi.org/10.1016/0092-8674(84)90284-8) | 1984 | The first consensus              |
