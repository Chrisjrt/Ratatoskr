# Ratatoskr

<p align="center">
  <img src="./Icon.png" width="300" alt="Ratatoskr logo">
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

```
git clone "git@github.com:Fabian-Bastiaanssen/Ratatoskr.git";
cd Ratatoskr;
pip install -e .
```

If installing from source, Ratatoskr requires the following dependencies to also be installed within the users $PATH:

#### Required dependencies

- [ncbi datasets](https://github.com/ncbi/datasets)

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
  --version         -v                                      Show the version and exit.
  --help            -h                                      Show this message and exit.       
```

### Outputs

Ratatoskr produces ouputs in three folders:

- `16S`:
  - This folder will contain a single file `16S.fasta` that contains all 16S rRNA sequences found for the taxonomic types.
- `genomes`:
  - This will contain two subfolders `fastas` and `genbank` that will contain all genomes found for the taxonomic types in FASTA and GenBank format, respectively.
  - Each file name genome file will be named after the associated GenBank accession e.g. `GCA_000012305.1.fna` or `GCA_000012305.1.gbff`
- `metadata`:
  - This will contain three files: `sequence_metadata.tsv`, `taxonomic_metadata.tsv`, and `phenotypic_metadata.tsv`
  - `sequence_metadata.tsv` contains information on whether any genome/16S rRNA sequence was found for a type requested, and if so the associated accesssions.
  - `taxonomic_metadata.tsv` contains the full taxonomic information on the types retrieved, including full taxonomic lineage, references for their provenance and list inclusions, authority, and binomial synonyms. 
  - `phenotypic_metadata.tsv` contains all phenotypic data from BacDive, including any results from API tests, colony and cellular morphology information, and any culture information that was available. 

### Example

To retrieve sequences and metadata for taxonomic type strains belonging to the Lachnospiraceae: 

```
ratatoskr run --input Lachnospiraceae --output output/Lachnospiraceae
```

## Citation

If you publish results from Ratatoskr please cite the following:

https://github.com/Fabian-Bastiaanssen/Ratatoskr
