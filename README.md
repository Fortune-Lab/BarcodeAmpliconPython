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

To set up a conda environment with all dependencies for this package, run the following command:

```bash
conda create -n barcode_amplicon_env python=3.9 numpy pandas matplotlib -c conda-forge -y
```

This creates a new environment named `barcode_amplicon_env` with Python 3.9 (meeting the minimum requirement) and installs the required dependencies (numpy, pandas, matplotlib) from the conda-forge channel.

After creating the environment, activate it with:

```bash
conda activate barcode_amplicon_env
```

Then, install the package itself from the project directory:

```bash
pip install .
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
