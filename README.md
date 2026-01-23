# Ratatoskr

<p align="center">
  <img src="./Icon.png" width="200" alt="Ratatoskr logo">
</p>

Retrieve 16S rRNA sequences, genome sequences, and metadata for taxonomic type strains

## Quick Start

The easiest way to install Ratatoskr is via conda:

```
*doesn't work yet*
conda install -c conda-forge -c bioconda ratatoskr
```

Run ratatoskr as follows:

```
ratatoskr run --input <name of taxa that should be retrieved> \
              --output <path to generate output to> \
```

## Installation

Before using Ratatoskr, please ensure you have an NCBI API key, and an account with DSMZ to use their LPSN and BacDive APIs. You'll need these to access certain functionalities within Ratatoskr. If you don't have an NCBI API key or DSMZ account yet, you can obtain these by following these steps:

### NCBI API Key

1. Visit the [NCBI website](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us).
2. Follow the instructions provided to create an NCBI account or log in if you already have one.
3. Generate a new API key or copy an existing one if available.

### DSMZ account

1. Visit the [LPSN API website](https://api.lpsn.dsmz.de/) and click on the link to 'login or register'
2. On the new page, click on the 'Register' button at the bottom of the prompt box
3. Follow the registration steps

The NCBI API key and DSMZ email/password will be required to use Ratatoskr. A prompt requesting these information will appear during the Ratatoskr run (Note: Ratatoskr does not store this information; it's only used to interface with the NCBI/LPSN/BacDive APIs).

### Bioconda

```
*doesn't work yet*
mamba create -n ratatoskr -c conda-forge -c bioconda ratatoskr
```

### Source

If installing from source, Ratatoskr requires the following dependencies to also be installed within the users $PATH:

#### Required dependencies

- [ncbi datasets](https://github.com/ncbi/datasets)

To install from source:

```
git clone "git@github.com:Fabian-Bastiaanssen/Ratatoskr.git";
cd Ratatoskr;
pip install -e .
```


## Usage

### Help

To access the help menu use the `-h` option:

```
Usage: ratatoskr.py [OPTIONS] COMMAND [ARGS]...

Options:
  --version  -v  Show the version and exit.
  --help     -h  Show this message and exit.

Commands:
  run       Run the ratatoskr pipeline.
```

### Running Ratatoskr

To run Ratatoskr use the `run` command:

```
Usage: ratatoskr.py run [OPTIONS]

    Run the ratatoskr pipeline.

Options:
  --input           -i  TEXT                                Name of the taxa of interest. [required]
  --output_path     -o  PATH                                Specify the desired path for creation of the output folder [required]
  --level           -l  [domain|kingdom|phylum|class|       Specify the taxanomic level to check for your taxon. [default: auto]
                          order|family|genus|species|auto]  
  --threads         -t  INTEGER RANGE [x>=1]                Number of threads to be used. [default: 1]
  --force           -f                                      Force overwrite of output directory.
  --dev_mode        -d                                      Run in development mode.
  --skip_download   -s                                      Skip sequence download steps.
  --version         -v                                      Show the version and exit.
  --help            -h                                      Show this message and exit.       
```

### Outputs

Ratatoskr produces ouputs in three main folders:

- `taxonomy/`:
  - `taxonomic_metadata.tsv`: TSV containing full taxonomy, plus additional taxonomy-associated information
- `sequences/`:
  - `sequence_metadata.tsv`: TSV containing 16S rRNA accessions, genome accessions, and genome completeness information. 
  - `16S/`:
    - `16S.fasta`: Single file containing all downloaded 16S rRNA sequences in FASTA format
  - `genomes/`:
    - `fastas/`: Folder containing all downloaded genomes in FASTA format
    - `genbank/`: Folder containing all downloaded genomes in GenBank format
- `characteristics/`
  - `API_results/`: Folder containing TSVs describing the results for any API tests associated with the taxon of interest (one file per API test)
  - `fatty_acid_profile_metadata.tsv`: TSV describing any information on fatty acid profiles for the taxon of interest
  - `general_characteristics`: TSV describing any general characteristics for the taxon of interest (e.g. Gram stain, motility, cell shape/width/length, Spore formation, and Oxygen_tolerance)
  - `metabolite_utilisation_metadata.tsv`: TSV describing any metabolite utilisation information associated with the taxon of interest

### Example

To retrieve sequences and metadata for taxonomic type strains belonging to the Lachnospiraceae: 

```
ratatoskr run --input Lachnospiraceae --output output/Lachnospiraceae
```

## Citation

If you publish results from Ratatoskr please cite the following:

https://github.com/Fabian-Bastiaanssen/Ratatoskr
