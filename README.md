# FLAMS: Find Lysine Acylations & other Modification Sites

A bioinformatics tool to analyze the conservation of post-translational modifications (PTMs), by means of a position-based search against the Compendium of Protein Lysine Modifications (CPLM database) v.4 and the experimental PTM sites in dbPTM and SCOP3P. FLAMS is available as command-line tool and as a [web service](https://www.biw.kuleuven.be/m2s/cmpg/research/CSB/tools/flams/).

# Table of contents

1.  [Introduction](#introduction)
2.  [System requirements](#system-requirements)
    1.  [General dependencies](#general-dependencies)
    2.  [Third-party dependencies](#third-party-dependencies)
3.  [Installation](#installation)
4.  [Usage](#usage)
    1. [Example use case](#example-use-case)
5.  [Output](#output)
6.  [Supported PTMs](#supported-ptms)
    1. [Supported PTM databases](#supported-ptm-databases)
    2. [Supported PTM types](#supported-ptm-types)
    3. [Local CPLM, dbPTM and SCOP3P install](#local-cplm-dbptm-and-scop3p-install)
7.  [Contact](#contact)
8.  [References](#references)
9.  [License](#license)

## Introduction

FLAMS is a bioinformatics tool to analyze the conservation of post-translational modifications, by means of a position-based search against CPLM, dbPTM* and SCOP3P*. FLAMS can be used (i) to quickly verify whether modifications in a specific protein have been reported before, (ii) to assess whether findings in one species might translate to other species, and (iii) to systematically assess the novelty and conservation of reported  modification sites.

The tool takes as input a protein (identifier or sequence) and the position of an amino acid. This repository contains the command-line tool `FLAMS`, which obtains an overview of the previously reported post-translational modifications matching your query, by using the following scripts:

* *input.py*: processing the user-provided input
* *cplmv4.py*, *dbptm.py*, *scop3p.py*, *uniprot.py* and *setup.py*: downloading and preparing the modification-specific databases
* *run_blast.py*: searching your query against the databases of proteins with post-translational modifications
* *display.py*: formatting the list of conserved post-translational modifications to a tab delimited output file
* *utils.py*: dealing with OS-dependent directory systems

FLAMS is also available as a web service at https://www.biw.kuleuven.be/m2s/cmpg/research/CSB/tools/flams/ .

*: only data with direct experimental support from these databases is incorporated in FLAMS.

## System requirements

Linux 64-bit, Windows and Mac OS supported.

### General dependencies

* Python3 (>=3.10, <3.12)

### Third-party dependencies

* [BLAST+ (>=2.13, tested until 2.16)](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)

## Installation

The recommended installation for Mac OS and Linux is through conda:

`conda install -c conda-forge -c bioconda flams`

It is also possible to install FLAMS through pip (recommended installation for Windows):

`pip install flams`

Please note that the pip install requires users to have BLAST+ installed locally and available in PATH. For more information on how to install BLAST+ on Windows, click [here](https://www.ncbi.nlm.nih.gov/books/NBK52637/) .

## Usage

Run the tool:

`FLAMS [-h] (--in inputFilePath | --id UniProtID | --batch batchFilePath) [-p position] [--range errorRange] [-o outputFilePath] [-d dataDir] [-t threadsBLAST] [-e evalueBLAST]
            [-m modification [modification ...]] `

Required argument:
* one of:
  * `inputFilePath`, used with `--in`, is the path to a .fasta file with the protein you wish to query against. (has to contain only 1 protein)
  * `UniProtID`, used with `--id`, is the UniProt ID of the protein you wish to query against.
  * `batchFilePath`, used with `--batch`, is the path to a tab seperated file for batch runs. The file should contain 1 entry per line, with UniProt IDs in the 1st column, and positions in the 2nd column.

Required argument when running FLAMS with --id/--in:
* `position` is the position of a (modified) residue in the protein, which you want to query against.

Optional arguments:
* `errorRange` is an number of positions before and after `position` to also search for modifications. [default: 0]
* `outputFilePath` is the path to where the result will be saved (in a .tsv file format). [default: out.tsv] If FLAMS is run with --batch, the specified -o/--output is used as preposition, followed by '\_$UniProtID\_$position.tsv'. [default: '']
* `dataDir` is the path to the directory where intermediate files (the UniProt sequence files) are stored. [default: $PWD/data]
* `threadsBLAST` is a BLAST parameter, allows you to speed up the search by multithreading. [default: 1]
* `evalueBLAST` is a BLAST parameter, allows you to filter out low quality BLAST hits. [default: 0.01]
* `modification` is a space-separated list of modifications (all lower case) to search for at the given position. Possible values are any (combinations) of the CPLM, dbPTM, and SCOP3P modifications. We also provide aggregated combinations for each amino acid ($AA-All), and the CPLM combinations. For a full list of all supported PTMs, and how they are named, see the [Supported PTM types](#supported-ptm-types) section of the README. In general, PTMs are written all lowercase, and spaces within a PTM name are replaced by underscores. [default: K-All]

### Example use case

We provide two example use cases for FLAMS:

With the following command, you search whether the TatA (UniProt ID: A0A916NWA0) acetylation on K66 in *Dehalococcoide mccartyi* strain CBDB1, as described by [Greiner-Haas (2021)](https://doi.org/10.3390/microorganisms9020365), had been previously detected.

`FLAMS --in A0A916NWA0.fa -p 66 -m acetylation -o tatA.tsv`

With the following command, you search whether the *Mycobabcterium smegmatis*' FadD2 (UniProt ID: A0QQ22) K537 is known to carry any modifications of the 'acylations' category, similar to what was reported by [Xu (2020)](https://doi.org/10.1128/mSystems.00424-19).

`FLAMS --id A0QQ22 -p 537 -m CPLM-Acylations -o FadD2.tsv`

You can find the example input and output data in the folder `test_data`. The output data is organized in folders reflecting the FLAMS version used to generate it, as the output can vary depending on the exact FLAMS version (due to FLAMS database updates).

For more example use cases, see the Supplementary information of the paper.

## Output

The output file is a .tsv containing one row per modification that matched the query, i.e., a modification aligning (within the user-specified range) to the query position, in a protein similar to the query protein. In case of batch jobs (ran with --batch), one output file per query (= a single line in the batch job file) will be generated.

The output file contains 16 columns:

* UniProt ID: UniProt identifier of the matched protein
* Protein name: protein name of the matched protein
* Modification: the type of modification found in the matched protein
* $AA location: the location of this matched modification in the matched protein
* $AA window: the local sequence containing the conserved modification (window of five amino acids before and after°)
* Species: the textual description of the species of the matched protein
* BLAST E Value: E value of BLASTp search of the matched protein against your query protein
* BLAST identity: % identity of BLASTp search of the matched protein against your query protein
* BLAST coverage: % coverage of BLASTp search of the matched protein against your query protein
* CPLM ID: CPLM ID of matched protein modification (if found in CPLM, otherwise empty)
* CPLM evidence code: CPLM evidence code of matched protein modification. Can be Exp(erimental), Dat(abase) or both. (if found in CPLM, otherwise empty)
* CPLM evidence links: CPLM evidence link of matched protein modification. Can be PubMed ID (for Exp.), or a database identifier (for Dat) or both. (if found in CPLM, otherwise empty)
* dbPTM evidence code: dbPTM evidence code of matched protein modification (if found in dbPTM, otherwise empty)
* dbPTM evidence links: dbPTM evidence link of matched protein modification. Refers to PubMed IDs. (if found in dbPTM, otherwise empty)
* SCOP3P evidence code: SCOP3P evidence code of matched protein modification. Can be Exp(erimental), Comb(ined) or both. PTMs derived from PRIDE have been assigned to Exp(erimental). (if found in SCOP3P, otherwise empty)
* SCOP3P evidence links: SCOP3P evidence link of matched protein modification. Refers to PubMed IDs. (if found in SCOP3P, otherwise empty)

°: window can be smaller than the [-5;+5] window if the sequence alignment ends sooner, which can happen for modified sites near the start/end of the protein

## Supported PTMs

### Supported PTM databases

FLAMS updates its search databases regularly. To get an overview of the supported databases, see the table below.

|FLAMS version|CPLM version|dbPTM version|SCOP3P version|database available for download|UniProt release|
|:----|:----|:----|:----|:----|:----|
|v1.1.7|v4 (Feb '25 update)|2026_January|2026_01|[yes](https://doi.org/10.5281/zenodo.18478351)|2026_01|
|v1.1.6|v4 (Feb '25 update)|2025_July|(not included)|[yes](https://doi.org/10.5281/zenodo.16737546)|2025_03|
|v1.1.5|v4|2025_January|(not included)|[yes](https://doi.org/10.5281/zenodo.14616210)|2024_06|
|v1.1.4|v4|2024_April|(not included)|[yes](https://doi.org/10.5281/zenodo.10958721)|2024_02|
|v1.1.0-3|v4|2023_November|(not included)|[yes](https://doi.org/10.5281/zenodo.10171879)|2023_05|
|v1.0|v4|(not included)|(not included)|[yes](https://cplm.biocuckoo.cn/Download.php)|NA|

Please note that only part of dbPTM and SCOP3P is integrated into FLAMS. For dbPTM, FLAMS only fetches the PTM sites with experimental evidence, as found [here](https://biomics.lab.nycu.edu.tw/dbPTM/download.php). For SCOP3P, FLAMS filters out any PTM entries that are derived from UniProt entries that have only sequence similarity-based evidence. As neither dbPTM, nor SCOP3P store complete protein sequences, these are fetched during database creation based on UniProt identifiers reported by dbPTM/SCOP3P and the UniProt release available at the time of database creation. As a consequence, FLAMS database updates can change the content of the PTM databases, beyond the simple addition of new SCOP3P, dbPTM and/or CPLM entries, reflecting changes in UniProt. The most common UniProt changes affecting FLAMS databases are removed UniProt entries (leading to the removal of PTM entries on the affected protein in our database) and sequence updates. We are aware of this issue, impacting the completeness and interpretation of FLAMS' results, and will consider solutions in future FLAMS releases.

Instructions on how to download the CPLM, dbPTM and SCOP3P database yourself are in [section 'Local CPLM, dbPTM and SCOP3P install'](#local-cplm-dbptm-and-scop3p-install). This is not recommended, as it takes multiple hours to generate some databases.

### Supported PTM types

FLAMS allows searches for all PTM types included in CPLM, and for those with experimental evidence in dbPTM and SCOP3P. An overview of the PTM types, how to call them in FLAMS, how they are called in CPLM and/or dbPTM and/or SCOP3P, and on which amino acid they can be found is given in the table below. This table can also be found as a tab seperated file named FLAMS_supported_ptms_v11.txt .

|FLAMS PTM name|CPLM name|dbPTM name|SCOP3P name|A (Ala)|C (Cys)|D (Asp)|E (Glu)|F (Phe)|G (Gly)|H (His)|I (Ile)|K (Lys)|L (Leu)|M (Met)|N (Asn)|P (Pro)|Q (Gln)|R (Arg)|S (Ser)|T (Thr)|V (Val)|W (Trp)|Y (Tyr)|CPLM-Acylations|CPLM-Ubs|CPLM-Others|CPLM-All|
|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
|acetylation|Acetylation|Acetylation| |X|X|X|X| |X| | |X| |X| |X| |X|X|X|X| |X|X| | |X|
|adp-ribosylation| |ADP-ribosylation| | |X|X|X| |X|X| |X| | |X| | |X|X| | | |X| | | | |
|amidation| |Amidation| |X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| | | | |
|ampylation| |AMPylation| | | | | | | | | | | | | | | | |X|X| | |X| | | | |
|benzoylation|Benzoylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|beta-hydroxybutyrylation|β-Hydroxybutyrylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|biotinylation|Biotinylation|Biotinylation| | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|blocked_amino_end| |Blocked amino end| |X|X|X|X| |X|X|X| |X|X|X|X|X|X|X|X|X| | | | | | |
|butyrylation|Butyrylation|Butyrylation| | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|carbamidation| |Carbamidation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|carboxyethylation|Carboxyethylation|Carboxyethylation| | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|carboxylation|Carboxylation|Carboxylation| | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|carboxymethylation|Carboxymethylation| | | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|cholesterol_ester| |Cholesterol ester| |  | | | | |X| | | | | | | | | | | | | | | | | | |
|citrullination| |Citrullination| | | | | | | | |  | | | | | | |X| | | | | | | | | |
|crotonylation|Crotonylation|Crotonylation| | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|c-linked_glycosylation| |C-linked Glycosylation| | | | | | | | | | | | | | | | | | | |X| | | | | |
|deamidation| |Deamidation| | | | | | | | | | | | |X| |X| | | | | | | | | | |
|deamination| |Deamination| | | | | | | | | |X| | | | | | | | | | | | | | | |
|decanoylation| |Decanoylation| | | | | | | | | | | | | | | | |X|X| | | | | | | |
|decarboxylation| |Decarboxylation| | | |X| | | | | | | | | | | | | |X| | | | | | | |
|dephosphorylation| |Dephosphorylation| | | | | | | | | | | | | | | | |X|X| | |X| | | | |
|dietylphosphorylation|Dietylphosphorylation| | | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|disulfide_bond| |Disulfide bond| | |X| | | | | | | | | | | | | | | | | | | | | | |
|d-glucuronylation| |D-glucuronoylation| | | | | | |X| | | | | | | | | | | | | | | | | | |
|farnesylation| |Farnesylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|formation_of_an_isopeptide_bond| |Formation of an isopeptide bond| | | | |X| | | | | | | | | |X| | | | | | | | | | |
|formylation|Formylation|Formylation| | | | | | |X| | |X| |X| | | | | | | | | |X| | |X|
|gamma-carboxyglutamic_acid| |Gamma-carboxyglutamic acid| | | | |X| | | | | | | | | | | | | | | | | | | | |
|geranylgeranylation| |Geranylgeranylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|glutarylation|Glutarylation|Glutarylation| | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|glutathionylation| |Glutathionylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|glycation|Glycation| | | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|gpi-anchor| |GPI-anchor| |X|X|X| | |X| | | | | |X| | | |X|X| | | | | | | |
|hmgylation|HMGylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|hydroxyceramide_ester| |Hydroxyceramide ester| | | | | | | | | | | | | | |X| | | | | | | | | | |
|hydroxylation|Hydroxylation|Hydroxylation| | |X|X|X|X| |X|X|X|X| |X|X| |X|X|X|X|X|X| | |X|X|
|iodination| |Iodination| | | | | | | | | | | | | | | | | | | | |X| | | | |
|lactoylation| |Lactoylation| | | | | | | | | |X| | | | | | | | | | | | | | | |
|lactylation|Lactylation|Lactylation| | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|lipoylation|Lipoylation|Lipoylation| | | | | | | | | |X| | | | | | | | | | | | | |X|X|
|malonylation|Malonylation|Malonylation| | | | | | | | | |X| |  | | | | | | | | | |X| | |X|
|methylation|Methylation|Methylation| |X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X|X| |X| | |X|X|
|mgcylation|MGcylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|mgylation|MGylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|myristoylation| |Myristoylation| | |X| | | |X| | |X| | | | | | | | | | | | | | | |
|neddylation|Neddylation|Neddylation| | | | | | | | | |X| | | | | | | | | | | | |X| |X|
|nitration| |Nitration| | | | | | | | | | | | | | | | | | | | |X| | | | |
|n-carbamoylation| |N-carbamoylation| |X| | | | | | | | | | | | | | | | | | |  | | | | |
|n-linked_glycosylation| |N-linked Glycosylation| | | |X| | | | |X|X| | |X| | |X|X|X|X|X| | | | | |
|n-palmitoylation| |N-palmitoylation| | |X| | | |X| | |X| | | | | | | | | | | | | | | |
|octanoylation| |Octanoylation| | | | | | | | | | | | | | | | |X|X| | | | | | | |
|oxidation| |Oxidation| | |X| | | | | | | |X|X| | | | |X| | |X| | | | | |
|o-linked_glycosylation| |O-linked Glycosylation| | | | | | | | | |X| | | |X| | |X|X| | |X| | | | |
|o-palmitoleoylation| |O-palmitoleoylation| | | | | | | | | | | | | | | | |X| | | | | | | | |
|o-palmitoylation| |O-palmitoylation| | | | | | | | | | | | | | | | |X|X| | | | | | | |
|phosphatidylethanolamine_amidation| |Phosphatidylethanolamine amidation| | | | | | |X| | | | | | | | | | | | | | | | | | |
|phosphoglycerylation|Phosphoglycerylation| | |  | | | | | | | |X| | |  | | | | | | | | | | |X|X|
|phosphorylation| |Phosphorylation|Phosphorylation|X|X|X|X|X|X|X|X|X|X| |X|X|X|X|X|X|X|X|X| | | | |
|propionylation|Propionylation|Propionylation| | | | | | | | | |X| | | | | | | | | | | |X| | |X|
|pupylation|Pupylation| | | | | | | | | | |X| | | | | | | | | | | | |X| |X|
|pyrrolidone_carboxylic_acid| |Pyrrolidone carboxylic acid| | | | |X| | | | | | | | | |X| | | | | | | | | | |
|pyrrolylation| |Pyrrolylation| | |X| | | | | | |  | | | | | | | | | | | | | | | |
|pyruvate| |Pyruvate| | |X| | | | | | | | | | | | | |X| | | | | | | | |
|serotonylation| |Serotonylation| | | | |  | | | | | | | | | |X| | | | | | | | | | |
|stearoylation| |Stearoylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|succinylation|Succinylation|Succinylation| | |X| | | | | | |X| | | | | | | | | |X| |X| | |X|
|sulfation| |Sulfation| | |X| | | | | | | | | | | | | |X|X| |  |X| | | | |
|sulfhydration| |Sulfhydration| |  |X| | | | | | | | | | | | | | | | | | | | | | |
|sulfoxidation| |Sulfoxidation| | | | | | | | | | | |X| | | | | | | | | | | | | |
|sumoylation|Sumoylation|Sumoylation| | | | | | | | | |X| | | | | | | | | | | | |X| |X|
|s-archaeol| |S-archaeol| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-carbamoylation| |S-carbamoylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-cyanation| |S-Cyanation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-cysteinylation| |S-cysteinylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-diacylglycerol| |S-diacylglycerol| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-linked_glycosylation| |S-linked Glycosylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-nitrosylation| |S-nitrosylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|s-palmitoylation| |S-palmitoylation| | |X| | | | | | | | | | | | | | | | | | | | | | |
|thiocarboxylation| |Thiocarboxylation| | | | | | |X| | | | | | | | | | | | | | | | | | |
|ubiquitination|Ubiquitination|Ubiquitination| | |X| | | | | | |X| | | | | |X|X| | | | | |X| |X|
|umpylation| |UMPylation| | | | | | | | | | | | | | | | |X|X| | |X| | | | |
|2-hydroxyisobutyrylation|2-Hydroxyisobutyrylation| | | | | | | | | | |X| | | | | | | | | | | |X| | |X|

### Local CPLM, dbPTM and SCOP3P install

It is possible to install the CPLM, dbPTM and SCOP3P databases directly, instead of using the pre-generated databases that are hosted on Zenodo. This is however **not** recommended as the download takes several hours for larger databases, such as phosphorylation, ubiquitination and acetylation.

However, if desired, follow these instructions to modify the scripts:

0. Make sure you are working in an environment with the correct dependencies:

  * on Linux/MacOS: create a FLAMS conda environment, get all dependencies by installing FLAMS as specified for Linux/MacOS.
  * on Windows: create a FLAMS conda environment, get all dependencies by installing FLAMS as specified for Windows. Make sure BLAST+ is correctly installed.

1. Download the latest FLAMS version from GitHub.

2. Adapt the scripts:

  * on a fresh install (= never ran FLAMS before, so no FLAMS databases yet):
    - go to `src/flams/databases/setup.py`
    - comment out lines 515-520 (function `_generate_blastdb_if_not_up_to_date` - the try/except _get_fasta_from_zenodo)
    - uncomment line 524 (function `_generate_blastdb_if_not_up_to_date` - the _get_fasta_for_blast)

  * on a FLAMS version with previously generated BLAST databases:
    - go to `src/flams/databases/setup.py`
    - change the version numbers of the databases you wish to update on lines 76-472.  E.g.:

    `"2-hydroxyisobutyrylation": ModificationType(
        "2-hydroxyisobutyrylation", 1.0, [ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation")],
        ["K"]
      ),`

    becomes

    `"2-hydroxyisobutyrylation": ModificationType(
        "2-hydroxyisobutyrylation", 2.0, [ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation")],
        ["K"]
      ),`

    - comment out lines 515-520 (function `_generate_blastdb_if_not_up_to_date` - the try/except _get_fasta_from_zenodo)
    - uncomment line 524 (function `_generate_blastdb_if_not_up_to_date` - the _get_fasta_for_blast)

3. Install your adapted FLAMS version locally:

    `python -m pip install ./PathToLocalFLAMS`

## Contact

Computational Systems Biology, KU Leuven.

## References

If you use FLAMS in your work, please cite:

Longin, H. *et al* (2024) "FLAMS: Find Lysine Acylations and other Modification Sites." Bioinformatics. 40(1):btae005.

In addition, FLAMS relies on third-party software & databases:

Altschul, S.F. *et al* (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.

Chung, C.-R. *et al* (2025) "dbPTM 2025 update: comprehensive integration of PTMs and proteomic data for advanced insights into cancer research." Nucleic Acids Research. 53(D1):D377–D386.

Ramasamy, P. *et al.* (2020) "Scop3P: A Comprehensive Resource of Human Phosphosites within Their Full Context." J. Proteome Res. 19(8):3478-3486.

Zhang, W. *et al* (2021) "CPLM 4.0: an updated database with rich annotations for protein lysine modifications." Nucleic Acids Research. 44(5):243–250.

## License

FLAMS is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
