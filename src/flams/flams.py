#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen, Retro212
"""

import logging
import os
import shutil
import sys
from pathlib import Path

from .databases.setup import update_db_for_modifications
from .display import create_output
from .input import parse_args, validate_input_fasta, validate_input_uniprot, validate_input_batch, get_position_in_seq, _get_length_protein
from .run_blast import run_blast

""" FLAMS
FLAMS is a tool that serves to find previously identified amino acid modification sites,
by enabling a position-based search of the CPLM database (CPLM v.4, Zhang, W. et al. Nucleic Acids Research. 2021, 44 (5): 243–250.).
FLAMS can be used
(i) to quickly verify whether modifications in a specific protein have been reported before,
(ii) to assess whether findings in one species might translate to other species, and
(iii) to systematically assess the novelty and conservation of all reported amino acid modification sites.

The tool takes as input a protein (ID/sequence) and the position of an amino acid.
"""

def is_available(program):
    """
    This function verifies the installation of third-party dependencies and prints out the result to users.

    Parameters
    ----------
    program: program
        Third-party dependency program.

    """
    if shutil.which(program) is not None:
        logging.info("Checking third-party dependencies. Installation of " + program + " : OK.")
    else:
        logging.error("Checking third-party dependencies. Installation of " + program + " failed verification: it is not available on the path... exiting FLAMS.")
        sys.exit()

def flams_blast(args, fasta, amino_acid, pos, out, pwd):
    """
    This function carries out FLAMS' main functionalities: running BLAST and formatting and outputting the results.

    Parameters
    ----------
    args:
        Arguments passed to flams
    fasta: Path
        Path to FASTA file containing query protein
    amino_acid: str
        Amino acid containing the post-translational modification under investigation
    pos: int
        User provided position in the query protein
    uniprot: str
        UniProt ID for protein containing the modification
    out: str
        Output file name
    pwd: Path
        Path of current working directory

    """
    # BLAST run
    result = run_blast(
        input=fasta,
        amino_acid_x=amino_acid,
        modifications=args.modification,
        x_pos=pos,
        x_range=args.range,
        num_threads=args.num_threads,
        evalue=args.evalue
        )
    # change working directory back (BLAST changes this)
    os.chdir(pwd)
    # format and generate output file with results
    create_output(output_filename=out, amino_acid_x=amino_acid, blast_records=result, len=_get_length_protein(fasta))
    # let user know where results can be found
    logging.info("Succesfully ran FLAMS! You can find your results in file: " + str(out) + " at " + str(pwd))


def main():
    """ Main function of FLAMS
    """
# set logger
    logging.basicConfig(
        level = logging.INFO,
        format = '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt = '%d/%m %H:%M:%S',
        force=True
    )

# check BLAST installation, local modification databases, validate arguments (non-protein specific)
    args = parse_args(sys.argv[1:])
    is_available('blastp')
    update_db_for_modifications(args.modification)

# get path from which FLAMS is run (need to save as BLAST will change this)
    pathCWD = Path().absolute()

# protein specific argument checks & run of BLAST
    # individual protein mode
    if not args.batch:
        # set output_file and protein_file (after checks)
        output_file = args.output
        if args.input:
            protein_file = validate_input_fasta(args, args.input, args.pos)
        if args.id:
            protein_file = validate_input_uniprot(args, args.id, args.pos)
        amino_acid_x = get_position_in_seq(args.pos, protein_file)
        # run BLAST, return results
        flams_blast(args, protein_file, amino_acid_x, args.pos, output_file, pathCWD)

    # batch mode
    if args.batch:
        # set outputfile
        preceedingOut = f"{args.output}_"
        if preceedingOut == "out.tsv_":
            preceedingOut = ""
        # read batch file (after checks)
        batchFile = open(validate_input_batch(args.batch), 'r')
        entries = batchFile.readlines()
        for entry in entries:
            # set protein_file (after checks), output_file, position
            uniprot = entry.split('\t')[0]
            pos = int(entry.split('\t')[1].strip())
            protein_file = validate_input_uniprot(args, uniprot, pos)
            output_file = f"{preceedingOut}{uniprot}_{pos}.tsv"
            amino_acid_x = get_position_in_seq(pos, protein_file)
            # run BLAST, return results
            flams_blast(args, protein_file, amino_acid_x, pos, output_file, pathCWD)


if __name__ == "__main__":
    main()
