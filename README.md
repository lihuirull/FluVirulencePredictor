# fluvp: Influenza Marker Extraction, Annotation, and Virulence Prediction Tool

`fluvp` is a versatile command-line utility designed for the extraction of influenza markers, annotation of sequence data using DIAMOND BLAST, and prediction of virulence levels. It processes individual files or entire directories and utilizes a pre-trained model for virulence prediction.

## Quick Start

```bash
conda create -n fluvp-env python=3.6
conda activate fluvp-env
git clone https://github.com/lihuirull/FluVirulencePredictor.git
cd fluvp
pip install -r requirements.txt
fluvp anno -i your_data.fasta -o output_directory
```

## Installation

`fluvp` can be installed using the following method:

### From Source

First, clone the repository and install Python dependencies:

```shell
git clone https://github.com/lihuirull/FluVirulencePredictor.git
cd FluHostPredictor
pip install .
```

This will also attempt to install the `diamond` dependency if it is available from the `bioconda` channel. If not available, you will need to install `diamond` manually.

Please note that installing with `pip` may not install non-Python dependencies such as `diamond`, which you will need to install separately.

## Usage

`fluvp` includes three subcommands: `anno`, `extract`, and `pred`.

### Annotate (`anno`)

Annotate a FASTA file or all FASTA files in a directory using DIAMOND BLAST against a flu database.

**Example:**

```bash
fluvp anno -i path/to/input.fasta -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input FASTA file or directory containing FASTA files (required).
- `-o, --output_directory`: Directory to save the output files (default: current directory).
- `-p, --prefix`: Prefix for the output filenames (default: none).
- `-e, --evalue`: E-value threshold for DIAMOND BLAST hits (default: 1e-5).
- `-u, --update_file`: If set, updates the FASTA file with annotations (flag).
- `-t, --threads`: Number of threads for DIAMOND BLAST (default: 10).

### Extract (`extract`)

Extract and process protein annotations from annotated FASTA files.

**Example:**

```bash
fluvp extract -i path/to/input.fasta -a path/to/annotations.csv -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input FASTA file or directory containing FASTA files (required).
- `-a, --anno_path`: Input annotation CSV file or directory containing annotation CSV files (required).
- `-o, --output_directory`: Directory to save the output files (default: current directory).
- `-p, --prefix`: Prefix for the output filenames (default: none).

### Predict (`pred`)

Predict the virulence level of new data using a trained model.

**Example:**

```bash
fluvp pred -i path/to/marker_data.csv -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input CSV file with marker data or directory containing such files (required).
- `-m, --model_path`: Path to the trained model file.
- `-th, --threshold`: Probability threshold for model prediction (default: 0.5).
- `-o, --output_directory`: Directory to save the prediction results (default: current directory).
- `-p, --prefix`: Prefix for the output filenames of the predictions (default: none).

## Dependencies

Python dependencies are listed in the `requirements.txt` file and can be installed using `pip`. However, `fluvp` also requires the `diamond` tool, which is not a Python package and needs to be installed separately. Instructions for installing `diamond` can be found at its official [documentation](https://github.com/bbuchfink/diamond/wiki).

For a comprehensive setup, including all dependencies, please follow the installation instructions provided in the sections above.

## Output

The `fluvp` tool generates the following outputs based on the subcommand executed:

- For `anno`: Annotated FASTA files or a directory of annotated FASTA files with added annotations based on DIAMOND BLAST results.
- For `extract`: A CSV file with extracted markers and annotations formatted as: `Strain ID,Virulence Markers,Number of Virulence Markers,Protein Type`
- For `pred`: This parameter outputs not only the predicted class labels for the new marker data based on the trained model but also the prediction probabilities. The current study employs a very small probability threshold.

## Notes

Ensure that the DIAMOND database and the trained model are compatible with the sequences you wish to annotate or predict. It is essential to use the same versions of the software and model that were used during the initial training and annotation for consistency and accuracy.

For further assistance or to report issues, please visit the `fluvp` GitHub repository issues section.