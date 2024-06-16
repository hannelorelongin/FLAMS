#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen
"""

import argparse
import logging
import os
import re
import requests
import sys
from pathlib import Path
from typing import Tuple

from Bio import SeqIO

from .databases import setup as db_setup

""" setup
This script deals with parsing the input and checking the validity of all provided arguments.
"""

logging.basicConfig(
        level = logging.INFO,
        format = '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt = '%d/%m %H:%M:%S',
        force=True
    )

# Functions called in flams.py (top-level)
def parse_args(sys_args) -> argparse.Namespace:
    """
    This function parses all provided arguments and runs all protein-independent argument checks.
    If no issues are found, it returns the arguments. Else, it throws an error.

    Parameters
    ----------
    sys_args:
        Arguments passed to FLAMS

    """
    parser = create_args_parser()
    args = parser.parse_args(sys_args)

    validate_modifications(args.modification)
    check_out_file(args.output)
    check_pos_provided(args)

    return args

def validate_input_fasta(args, fasta, pos) -> Path:
    """
    This function checks whether all FASTA file related arguments pass the checks.
    If all checks are passed, it returns the path to the user-specified FASTA file.

    Parameters
    ----------
    args:
        Arguments passed to flams
    fasta: Path
        Path to FASTA file containing query protein.
    pos: int
        User provided position in the query protein

    """
    check_fasta(fasta)
    check_pos_in_range(fasta, pos)
    check_amino_acid(fasta, pos, get_intersection(args.modification))
    logging.info(f"Based on your provided information, we are screening for the following modifications: {args.modification} (on amino acid {get_position_in_seq(pos, fasta)}).")
    return fasta


def validate_modifications(modifications):
    """
    This function checks whether the user provided modification arguments pass the checks.
    If not, it returns an error.

    Parameters
    ----------
    modifications: List[modification]
        List of modifications provided by the user

    """
    if modifications:
        # deal with umbrella terms for modifications
        untangle_modifications(modifications)
        # remove duplicate modifications
        modifications = list(set(modifications))
        # check if modification is known to FLAMS
        check_modification_exist(modifications)
        # check if modifications actually share a carrier amino acid
        check_modifications_share_carrier(modifications)
        logging.info(f"The set of amino acids capable of carrying all your specified modifications {modifications} was determined as {get_intersection(modifications)}.")


def validate_input_uniprot(args, uniprot, pos) -> Path:
    """
    This function checks whether all UniProt-related arguments pass the checks.
    If all checks are passed, it returns the path to the downloaded fasta file corresponding to this UniProt ID.

    Parameters
    ----------
    args:
        Arguments passed to flams
    uniprot: str
        UniProt ID for protein containing the modification
    pos: int
        User provided position in the query protein

    """
    check_data_dir(args)
    protein_file = get_uniprot_file(args, uniprot)
    validate_input_fasta(args, protein_file, pos)

    return protein_file

def validate_input_batch(batch) -> Path:
    """
    This function checks whether all batch file-related arguments pass the checks.
    If all checks are passed, it returns the path to the batch file.

    Parameters
    ----------
    batch: Path
        Path to batch file containing queries as tab seperated protein- (as UniProt ID) and position-information, with one query per line.

    """
    check_batch_file(batch)
    check_batch_content(batch)
    return batch

# Function to create args parser
def create_args_parser():
    """
    This function creates an argument parser.

    """
    parser = argparse.ArgumentParser(
        description="Find Lysine Acylations & other Modification Sites."
    )

    # query proteins
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--in",
        dest="input",
        type=Path,
        help="Path to input .fasta file.",
        metavar="inputFilePath"
    )
    group.add_argument(
        "--id",
        type=str,
        help="UniProt ID of input protein.",
        metavar="UniProtID"
    )
    group.add_argument(
        "--batch",
        type=Path,
        help="Path to tab seperated input file for batch processing (1st column UniProt ID, 2nd column position). One query (UniProtID + position) per line." ,
        metavar="batchFilePath"
    )

    # position
    parser.add_argument(
        "-p",
        "--pos",
        type=int,
        help="Position in input protein that will be searched for conserved modifications.",
        metavar="position"
    )

    parser.add_argument(
        "--range",
        type=int,
        default=0,
        help="Allowed error range for position. [default: 0]",
        metavar="errorRange"
    )

    # File arguments
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("out.tsv"),
        help="Path to output .tsv file. [default: out.tsv] If FLAMS is run with --batch, the specified -o/--output is used as preposition, followed by '_$UniProtID_$position.tsv'. [default: '']",
        metavar="outputFilePath"
    )

    parser.add_argument(
        "-d",
        "--data_dir",
        type=Path,
        default=Path(os.getcwd()) / "data",
        help="Path to directory where intermediate files should be " +
        "saved. [default: $PWD/data]",
        metavar="dataDir"
    )

    # BLAST settings
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads to run BLAST with. [default: 1]",
        metavar="threadsBLAST"
    )


    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0.01,
        help="Desired E-value of BLAST run. [default: 0.01]",
        metavar="evalueBLAST"
    )

    parser.add_argument(
        "-m",
        "--modification",
        nargs="+",
        default=["K-All"],
        help="Space-seperated list of modifications (all lower case) to search for at the given position. Possible values are any (combinations) of the CPLM and dbPTM modifications. " +
        "We also provide aggregated combinations for each amino acid (AA-All), and the CPLM combinations. For a full list of all supported PTMs, and how they are named, see the " +
        "Supported PTM types section of the README. In general, PTMs are written all lowercase, and spaces within a PTM name are replaced by underscores. [default: K-All]",
        metavar="modification"
    )

    return parser


# Check functions (order alphabetically)
def check_amino_acid(fasta, pos, intersection):
    """
    This function checks whether the user provided position actually points to an amino acid that can carry all specified modifications.
    If not, it returns an error.

    Parameters
    ----------
    fasta: Path
        Path to FASTA file containing query protein
    pos: int
        User provided position in the query protein
    intersection: set
        Set of amino acids that can carry all user-specified modifications.

    """
    try:
        if not is_position_carrier_of_mod(pos, fasta, intersection):
            logging.error(
                f"Position {pos} does not point to an amino acid carrying your specified modifications: {_get_position_display_str(pos, fasta)} "
                )
            logging.error(f"Please provide a position that corresponds to one of the following amino acids {intersection}.")
            sys.exit()
    except IndexError as e:
        logging.error(f"{e}. Please provide an amino acid position smaller than the size of your protein.")
        sys.exit()

def check_batch_content(batch):
    """
    This function checks whether the provided batch file contents follow the format of 'UniProtID \tab position'.
    It checks whether the first column contains strings matching the UniProt ID REGEX (as obtained from UniProt), and the second column an integer.
    If not, it throws an error.

    Parameters
    ----------
    batch: Path
        Path to batch file containing queries as tab seperated protein- (as UniProt ID) and position-information, with one query per line.

    """
    batchFileContent = open(batch, 'r')
    entries = batchFileContent.readlines()
    line = 1
    for entry in entries:
        uniprotToTest = entry.split('\t')[0]
        posToTest = entry.split('\t')[1].strip()
        if re.fullmatch("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",uniprotToTest) is None:
            logging.error(f"Provided UniProt ID {uniprotToTest} on line {line} of batch file does not fulfill the UniProt ID formatting requirements. Please provide a correct UniProt ID. Exiting FLAMS...")
            sys.exit()
        try:
            int(posToTest)
        except ValueError:
            logging.error(f"Provided position {posToTest} on line {line} of batch file is not an integer. Please provide an integer. Exiting FLAMS...")
            sys.exit()
        line += 1


def check_batch_file(batch):
    """
    This function checks whether the provided batch file exists and is not a directory.
    If not, it throws an error.

    Parameters
    ----------
    batch: Path
        Path to batch file containing queries as tab seperated protein- (as UniProt ID) and position-information, with one query per line.

    """
    if not batch.exists():
        logging.error(f"Input file {batch} does not exist. Please provide the correct path to the batch input file. Exiting FLAMS...")
        sys.exit()
    if batch.is_dir():
        logging.error(f"Provided batch file {batch} is a directory name, not a file name. Please provide a batch filename instead. Exiting FLAMS...")
        sys.exit()


def check_data_dir(args):
    """
    This function checks if the data directory exists when a UniProtID is provided.
    If not, it throws an error.

    Parameters
    ----------
    args:
        Arguments passed to FLAMS

    """
    exit_message = "Please make sure the provided path is correct. Exiting FLAMS..."

    if args.data_dir.exists() and not args.data_dir.is_dir():
        logging.error("Provided data directory is not a " +
                    f"directory: {args.data_dir}. {exit_message}")
        sys.exit()

    if not args.data_dir.is_dir():
        if args.data_dir.parent.is_dir():
            os.mkdir(args.data_dir)
            logging.info(f"Data directory created: {args.data_dir}")
        else:
            logging.error("Provided path is not an existing " +
                        f"path: {args.data_dir}. {exit_message}")
            sys.exit()


def check_fasta(fasta):
    """
    This function checks whether the provided input FASTA file exists and is a valid FASTA file.
    If not, it throws an error.

    Parameters
    ----------
    fasta: Path
        Path to FASTA file containing query protein.

    """
    if not fasta.exists():
        logging.error(f"Input file {fasta} does not exist. Please provide the correct path to the input file or use --id (instead of --in) when specifying a UniProt ID. Exiting FLAMS...")
        sys.exit()

    if not is_valid_fasta_file(fasta):
        logging.error(f"Input file {fasta} is not a valid FASTA file. Exiting FLAMS...")
        sys.exit()


def check_modification_exist(modifications):
    """
    This function checks whether the modifications are part of the list of modifications known to FLAMS.
    If not, it throws an error.

    Parameters
    ----------
    modifications: List[modification]
        List of modifications provided by the user

    """
    for i in modifications:
        if i not in db_setup.MODIFICATIONS:
            logging.error(f"Invalid modification type {i}. Please choose a modification from the list specified in the docs. Exiting FLAMS... ")
            sys.exit()


def check_modifications_share_carrier(modifications):
    """
    This function checks whether the combination of modifications share a common carrier amino acid.
    If not, it throws an error.

    Parameters
    ----------
    modifications: List[modification]
        List of modifications provided by the user

    """
    if len(get_intersection(modifications)) == 0:
        logging.error(f'This combination of modifications ({modifications}) does not share an amino acid that carries them.')
        logging.error('Please only provide combinations of modifications that can occur on the same amino acid. Exiting FLAMS...')
        sys.exit()


def check_out_file(output_file):
    """
    This function checks whether the provided output file is not a directory name.
    If the output file is a directory, it throws an error.

    Parameters
    ----------
    output_file:
        Output file name

    """
    if output_file and output_file.is_dir():
        logging.error(f"Provided output: {output_file} is a directory name, not a file name. Please provide an output filename instead. Exiting FLAMS...")
        sys.exit()


def check_pos_in_range(fasta, pos):
    """
    This function checks whether the user provided position is actually part of the protein.
    If not, it returns an error.

    Parameters
    ----------
    fasta: Path
        Path to FASTA file containing query protein
    pos: int
        User provided position in the query protein

    """

    if not is_within_range(pos, fasta):
        logging.error(
            f"Please provide a lysine position smaller than the size " +
            f"of your protein ({_get_length_protein(fasta)})."
            )
        sys.exit()


def check_pos_provided(args):
    """
    This function checks whether a position argument is provided when FLAMS is being called with --id/--in, and absent when called with --batch.
    If the position argument is lacking with --in/--id, it throws an error.
    If the position argument is provided with --pos, it throws a warning.

    Parameters
    ----------
    args:
        Arguments passed to FLAMS

    """
    if args.input and args.pos is None:
        logging.error("Use of parameter --in requires specification of position parameter -p/--pos. Please provide parameter -p/--pos. Exiting FLAMS...")
        sys.exit()
    if args.id and args.pos is None:
        logging.error("Use of parameter --id requires specification of position parameter -p/--pos. Please provide parameter -p/--pos. Exiting FLAMS...")
        sys.exit()
    if args.batch and args.pos:
        logging.warning("Use of parameter --batch overwrites the parameter -p/--pos. Position is being read from batch file.")


# Is ... functions (ordered alphabetically)
def is_position_carrier_of_mod(pos: int, fasta: Path, intersection) -> bool:
    """
    This function assess whether the user provided position actually points to an amino acid that is part of the set of amino acids that can carry aall user specified modifications.
    It returns True if the checks pass, False otherwise.

    Parameters
    ----------
    pos: int
        User provided position in the query protein
    fasta: Path
        Path to FASTA file containing query protein
    intersection: set
        Set of amino acids that can carry all user-specified modifications.

    """
    # user provides position in 1-based indexing system
    return get_position_in_seq(pos, fasta) in intersection


def is_position_lysine(pos: int, fasta: Path) -> bool:
    """
    This function assess whether the user provided position actually points to a lysine in the query protein.
    It returns True if the checks pass, False otherwise.

    Parameters
    ----------
    pos: int
        User provided position in the query protein
    fasta: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    position_idx = pos - 1
    input_seq = SeqIO.read(fasta, "fasta").seq
    return input_seq[position_idx] == "K"


def is_valid_fasta_file(fasta: Path) -> bool:
    """
    This function checks whether the provided input FASTA file is a valid FASTA file.
    It returns True if the checks pass, False otherwise.

    Parameters
    ----------
    fasta: Path
        Path to FASTA file containing query protein.

    """
    try:
        SeqIO.read(fasta, "fasta")
        return True
    except Exception:
        return False


def is_within_range(pos: int, fasta: Path) -> bool:
    """
    This function assess whether the user provided position is actually part of the query protein.
    It returns True if the checks pass, False otherwise.

    Parameters
    ----------
    pos: int
        User provided position in the query protein
    fasta: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    position_idx = pos - 1
    length = _get_length_protein(fasta)
    return position_idx < length


# Get functions (ordered alphabetically)
def get_intersection(modifications):
    """
    This function returns the intersection of carrier amino acids for a list of modifications.

    Parameters
    ----------
    modifications: List[modification]
        List of modifications provided by the user

    """
    allowed_amino_acids = []
    for i in modifications:
        allowed_amino_acids.append(db_setup.MODIFICATIONS.get(i).aas)
    intersect_allowed = set.intersection(*[set(x) for x in allowed_amino_acids])
    return intersect_allowed


def _get_length_protein(fasta: Path) -> int:
    """
    This function returns the length of the protein stored in the FASTA file.

    Parameters
    ----------
    fasta: Path
        Path to FASTA file containing query protein

    """
    prot_seq = SeqIO.read(fasta, "fasta").seq
    return len(prot_seq)


def get_position_in_seq(pos: int, fasta: Path) -> str:
    """
    This function returns the amino acid at position pos of the protein stored in the FASTA file.

    Parameters
    ----------
    pos: int
        User provided position in the query protein
    fasta: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    position_idx = pos - 1
    input_seq = SeqIO.read(fasta, "fasta").seq
    return input_seq[position_idx]


def get_uniprot_file(args, uniprot) -> Path:
    """
    This function retrieves the protein input file, by either:
    - returning the path to the already downloaded FASTA file corresponding to the user provided UniProt ID (this happens if FLAMS is run multiple times for the same protein)
    - downloading the FASTA file from UniProt, based on the user provided UniProt ID, then returning the path to the downloaded protein fasta

    Parameters
    ----------
    args:
        Arguments passed to flams
    uniprot: str
        UniProt ID for protein containing the modification

    """
    filename = args.data_dir / f"{uniprot}.fasta.tmp"
    if os.path.isfile(args.data_dir / f"{uniprot}.fasta.tmp"):
        logging.info(f"Found FASTA file for UniProt ID {uniprot} at {filename}")
        return filename

    try:
        return retrieve_protein_from_uniprot(args, uniprot)
    except requests.HTTPError:
        logging.error("Non-existing UniProt ID. Please provide a valid UniProt ID. Exiting FLAMS...")
        sys.exit()


def _get_position_display_str(pos: int, fasta: Path) -> str:
    """
    This function returns a fragment of the sequence around a chosen position.

    Parameters
    ----------
    pos: int
        User provided position in the query protein
    fasta: Path
        Path to FASTA file containing query protein

    """
    # user provides position in 1-based indexing system
    pos_idx = pos - 1
    seq = SeqIO.read(fasta, "fasta").seq
    lower = max(0, pos_idx - 3)
    upper = min(len(seq), pos_idx + 3)
    prefix = "..." if lower > 0 else ""
    sufix = "..." if upper < len(seq) - 1 else ""
    pos_idx = len(prefix) + (pos_idx - lower)
    seq_row = f"{prefix}{seq[lower:upper]}{sufix}"
    pointer_row = " " * pos_idx + "^"
    return "".join(["\n", seq_row, "\n", pointer_row])


# Retrieve function
def retrieve_protein_from_uniprot(args,uniprot) -> Path:
    """
    This function downloads the FASTA file from UniProt, based on the provided UniProt ID.
    It returns the path to the downloaded protein FASTA file.

    Parameters
    ----------
    args:
        Arguments passed to flams
    uniprot: str
        UniProt ID for protein containing the modification

    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
    logging.info(f"Retrieving FASTA file for Uniprot ID {uniprot} at {url}")
    r = requests.get(url)

    r.raise_for_status()

    filename = args.data_dir / f"{uniprot}.fasta.tmp"
    with filename.open("w+") as f:
        f.write(r.text)

    logging.info(f"Stored FASTA file for UniProt ID {uniprot} at {filename}")
    return filename


# Other Functions
def untangle_modifications(modifications):
    """
    This function transforms the amino-acid & CPLM aggregate modification options Ubs, Acylations, Others, and All to their respective collection of modifications.

    Parameters
    ----------
    modifications: List[modification]
        List of modifications provided by the user

    """
    if (('CPLM-Ubs' in modifications) | ('CPLM-Acylations' in modifications) | ('CPLM-Others' in modifications) | ('CPLM-All' in modifications)
        | ('A-All' in modifications) | ('C-All' in modifications) | ('D-All' in modifications) | ('E-All' in modifications) | ('F-All' in modifications)
        | ('G-All' in modifications) | ('H-All' in modifications) | ('I-All' in modifications) | ('K-All' in modifications) | ('L-All' in modifications)
        | ('M-All' in modifications) | ('N-All' in modifications) | ('P-All' in modifications) | ('Q-All' in modifications) | ('R-All' in modifications)
        | ('S-All' in modifications) | ('T-All' in modifications) | ('V-All' in modifications) | ('W-All' in modifications) | ('Y-All' in modifications)) :
        if 'CPLM-Ubs' in modifications:
            modifications.remove('CPLM-Ubs')
            modifications.extend(['ubiquitination','sumoylation','pupylation','neddylation'])
        if 'CPLM-Acylations' in modifications:
            modifications.remove('CPLM-Acylations')
            modifications.extend(['lactylation','acetylation','succinylation','crotonylation','malonylation',
            'beta-hydroxybutyrylation','benzoylation','propionylation','2-hydroxyisobutyrylation','formylation',
            'hmgylation','mgcylation','mgylation','glutarylation','butyrylation'])
        if 'CPLM-Others' in modifications:
            modifications.remove('CPLM-Others')
            modifications.extend(['methylation','hydroxylation','phosphoglycerylation','biotinylation','lipoylation',
            'dietylphosphorylation','glycation','carboxymethylation','carboxyethylation','carboxylation'])
        if 'CPLM-All' in modifications:
            modifications.remove('CPLM-All')
            modifications.extend(['ubiquitination','sumoylation','pupylation','neddylation',
            'lactylation','acetylation','succinylation','crotonylation','malonylation',
            'beta-hydroxybutyrylation','benzoylation','propionylation','2-hydroxyisobutyrylation','formylation',
            'hmgylation','mgcylation','mgylation','glutarylation','butyrylation',
            'methylation','hydroxylation','phosphoglycerylation','biotinylation','lipoylation',
            'dietylphosphorylation','glycation','carboxymethylation','carboxyethylation','carboxylation'])
        if 'A-All' in modifications:
            modifications.remove('A-All')
            modifications.extend(['phosphorylation','acetylation','gpi-anchor','amidation','blocked_amino_end',
            'methylation','n-carbamoylation'])
        if 'C-All' in modifications:
            modifications.remove('C-All')
            modifications.extend(['phosphorylation','acetylation','pyruvate','gpi-anchor',
            'adp-ribosylation','amidation','hydroxylation','blocked_amino_end',
            'oxidation','methylation','sulfation','ubiquitination','carbamidation',
            'farnesylation','geranylgeranylation','glutathionylation','myristoylation',
            'n-palmitoylation','pyrrolylation','s-archaeol','s-carbamoylation',
            's-cyanation','s-cysteinylation','s-diacylglycerol','s-linked_glycosylation',
            's-nitrosylation','s-palmitoylation','stearoylation',
            'succinylation','sulfhydration','disulfide_bond'])
        if 'D-All' in modifications:
            modifications.remove('D-All')
            modifications.extend(['phosphorylation','acetylation','gpi-anchor',
            'adp-ribosylation','amidation','hydroxylation','blocked_amino_end',
            'methylation','n-linked_glycosylation','decarboxylation'])
        if 'E-All' in modifications:
            modifications.remove('E-All')
            modifications.extend(['phosphorylation','acetylation','adp-ribosylation',
            'amidation','hydroxylation','blocked_amino_end','methylation',
            'formation_of_an_isopeptide_bond','gamma-carboxyglutamic_acid',
            'pyrrolidone_carboxylic_acid'])
        if 'F-All' in modifications:
            modifications.remove('F-All')
            modifications.extend(['phosphorylation','amidation','hydroxylation',
            'methylation'])
        if 'G-All' in modifications:
            modifications.remove('G-All')
            modifications.extend(['phosphorylation','d-glucuronoylation',
            'cholesterol_ester','phosphatidylethanolamine_amidation',
            'thiocarboxylation','blocked_amino_end','formylation',
            'adp-ribosylation','gpi-anchor','myristoylation','n-palmitoylation',
            'amidation','methylation','acetylation'])
        if 'H-All' in modifications:
            modifications.remove('H-All')
            modifications.extend(['phosphorylation','blocked_amino_end',
            'adp-ribosylation','amidation','hydroxylation','methylation'])
        if 'I-All' in modifications:
            modifications.remove('I-All')
            modifications.extend(['phosphorylation','blocked_amino_end',
            'amidation','hydroxylation','methylation','n-linked_glycosylation'])
        if 'K-All' in modifications:
            modifications.remove('K-All')
            modifications.extend(['phosphorylation', 'deamination','formylation',
            'ubiquitination','sumoylation','pupylation','adp-ribosylation','lactoylation',
            'neddylation','lactylation','acetylation','succinylation','myristoylation',
            'crotonylation','malonylation','beta-hydroxybutyrylation','benzoylation',
            'propionylation','2-hydroxyisobutyrylation','formylation','n-palmitoylation',
            'hmgylation','mgcylation','mgylation','glutarylation','butyrylation',
            'methylation','hydroxylation','phosphoglycerylation','biotinylation',
            'lipoylation','dietylphosphorylation','glycation','carboxymethylation',
            'carboxyethylation','carboxylation','amidation','o-linked_glycosylation',
            'n-linked_glycosylation'])
        if 'L-All' in modifications:
            modifications.remove('L-All')
            modifications.extend(['phosphorylation','blocked_amino_end','oxidation',
            'amidation','hydroxylation','methylation'])
        if 'M-All' in modifications:
            modifications.remove('M-All')
            modifications.extend(['blocked_amino_end','formylation','oxidation',
            'amidation','sulfoxidation','methylation','acetylation'])
        if 'N-All' in modifications:
            modifications.remove('N-All')
            modifications.extend(['phosphorylation','blocked_amino_end',
            'deamidation','adp-ribosylation','gpi-anchor','amidation',
            'hydroxylation','methylation','n-linked_glycosylation'])
        if 'P-All' in modifications:
            modifications.remove('P-All')
            modifications.extend(['phosphorylation','blocked_amino_end','amidation',
            'hydroxylation','o-linked_glycosylation','methylation','acetylation'])
        if 'Q-All' in modifications:
            modifications.remove('Q-All')
            modifications.extend(['phosphorylation','hydroxyceramide_ester','serotonylation',
            'blocked_amino_end','formation_of_an_isopeptide_bond','deamidation',
            'pyrrolidone_carboxylic_acid','amidation','methylation'])
        if 'R-All' in modifications:
            modifications.remove('R-All')
            modifications.extend(['phosphorylation','blocked_amino_end','citrullination',
            'adp-ribosylation','amidation','hydroxylation','methylation',
            'n-linked_glycosylation','acetylation','ubiquitination'])
        if 'S-All' in modifications:
            modifications.remove('S-All')
            modifications.extend(['phosphorylation','decanoylation','octanoylation',
            'o-palmitoylation','umpylation','ampylation','blocked_amino_end',
            'o-palmitoleoylation','adp-ribosylation','gpi-anchor','sulfation',
            'oxidation','pyruvate','amidation','hydroxylation','o-linked_glycosylation',
            'methylation','n-linked_glycosylation','acetylation','ubiquitination',
            'dephosphorylation'])
        if 'T-All' in modifications:
            modifications.remove('T-All')
            modifications.extend(['phosphorylation','decarboxylation','decanoylation',
            'octanoylation','o-palmitoylation','umpylation','ampylation',
            'blocked_amino_end','gpi-anchor','sulfation',
            'amidation','hydroxylation','o-linked_glycosylation',
            'methylation','n-linked_glycosylation','acetylation', 'dephosphorylation'])
        if 'V-All' in modifications:
            modifications.remove('V-All')
            modifications.extend(['phosphorylation','blocked_amino_end','amidation',
            'hydroxylation','methylation','n-linked_glycosylation','acetylation'])
        if 'W-All' in modifications:
            modifications.remove('W-All')
            modifications.extend(['phosphorylation','c-linked_glycosylation','oxidation',
            'amidation','hydroxylation','succinylation','n-linked_glycosylation'])
        if 'Y-All' in modifications:
            modifications.remove('Y-All')
            modifications.extend(['phosphorylation','umpylation','iodination',
            'ampylation','adp-ribosylation','sulfation','nitration','amidation',
            'hydroxylation','o-linked_glycosylation','methylation','acetylation',
            'dephosphorylation'])
