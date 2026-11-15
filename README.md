# BarcodeAmpliconPython

Python tools for barcode amplicon analysis, specifically designed for processing 7mer barcodes. This package provides functionality for barcode reading, thresholding, and plotting.

## Installation

To install this package, ensure you have Python 3.9 or higher. Then, run:

```bash
pip install .
```

### Dependencies

- numpy
- pandas
- matplotlib

### Install in a conda environment

To set up a conda environment with all dependencies for this package, run from the project directory:

```bash
mamba env create -f environment.yml
```

This creates an environment named `barcode_amplicon` with Python 3.11 and the required dependencies (numpy, pandas, matplotlib, pigz) from conda-forge. The Slurm scripts in `chunked_pooled_library_scripts/` activate this environment by name.

After creating the environment, activate it with:

```bash
mamba activate barcode_amplicon
```

Then, install the package itself from the project directory. An editable install (`-e`) means a `git pull` updates the installed tools without reinstalling:

```bash
pip install --no-deps -e .
```


## Usage

After installation, you can use the command-line tools:

- `barcodeAmplicon-barcode-reader`: Read barcodes from FASTQ files
- `barcodeAmplicon-threshold-derivative`: Find thresholds using derivative method
- `barcodeAmplicon-threshold-regression`: Find thresholds using regression method
- `barcodeAmplicon-plot-threshold`: Plot threshold results

For example:

```bash
barcodeAmplicon-barcode-reader --input sample.fastq --output results.tsv
```

You can also import the package in Python:

```python
import barcodes7mer
```

## License

This project is licensed under the terms in the LICENSE file.
