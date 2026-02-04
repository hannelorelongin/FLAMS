#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen
"""

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Record import Blast, Alignment

from .databases import setup as db_setup
from .utils import get_data_dir

""" run_blast
This script contains all functions necessary to search through the proteins stored in the FLAMS database, with BLAST,
and to retrieve those that contain conserved amino acid X modifications.
"""

def run_blast(input, amino_acid_x, modifications, x_pos, x_range = 0, evalue = 0.01, num_threads = 1, **kwargs,):
    """
    This function runs the BLAST search and the following filter steps for each modification.
    Ultimately, it only returns conserved (within the specified range) protein modifications for similar proteins.
    It flattens the results to an array.

    Parameters
    ----------
    input: fasta
        Sequence file of query protein
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    modifications: str
        Space-seperated list of modifications (which are keys to any of the ModificationType's stored in the MODIFICATIONS dictionary)
    x_pos: int
        Position of amino acid X in query that is under investigation for conservation
    x_range: int (default: 0)
        Error margin for conservation of x_pos
    evalue: float (default: 0.01)
        BLAST parameter, e-value for BLAST run
    num_threads: int (default: 1)
        BLAST parameter, number of threads that can be used by BLAST

    """
    results = []
    input = input.absolute()
    for m in modifications:
        result = _run_blast(input, amino_acid_x, m, x_pos, x_range, evalue, num_threads)
        for r in result:
            results.append(r)
    return results


@dataclass
class ModificationHeader:
    """
    This dataclass consists of the different components contained in the header of each modification entry, and a function to parse it.

    Parameters
    ----------
    uniprot_id: str
        UniProt ID for protein containing the modification
    position: int
        Position at which the modification was detected
    length: int
        Length of protein in which the modification was detected
    protein_name: str
        Protein name of protein containing the modification
    modification: str
        Post-translational modification found at $position in protein with $uniprot_id
    database: str
        Indicator of database from which PTM was fetched (either CPLM or dbPTM)
    species: str
        Species that encodes the protein containing the modification
    db_id: str
        If database is CPLM: CPLM ID for each modification. If database is dbPTM or SCOP3P, this just states "dbPTM"/"SCOP3P".
    evidence_code: str
        Evidence code for each modification (either Exp., Dat. or a combination thereof for CPLM; Exp. for dbPTM; Exp., Comb. for SCOP3P)
    evidence_link: str
        Evidence link for each modification (PubMed ID for Exp./Comb., database code for Dat.)

    """
    uniprot_id: str
    position: int
    length: int
    database: str
    protein_name: str
    modification: str
    species: str
    db_id: str
    evidence_code: str
    evidence_link: str

    @staticmethod
    def parse(title: str) -> "ModificationHeader":

        regex = (
            r"(?P<uniprot_id>\S+)\|"
            r"(?P<position>\d+)\|"
            r"(?P<length>\d+)\|"
            r"(?P<database>\S+)"
            r" (?P<protein_name>\S+)\|(?P<modification>\S+)\|(?P<species>\S+) \[(?P<db_id>\S+)\|(?P<evidence_code>\S+)\|(?P<evidence_link>.+)\]"
        )
        vars = re.match(regex, title).groupdict()
        vars["position"] = int(vars["position"])
        vars["length"] = int(vars["length"])
        return ModificationHeader(**vars)


def _run_blast(input, amino_acid_x, modification, x_pos, x_range, evalue, num_threads = 1):
    """
    This function runs the BLAST search and the following filter steps for 1 modification.
    Ultimately, it only returns conserved (within the specified range) protein modifications for similar proteins.

    Parameters
    ----------
    input: fasta
        Sequence file of query protein
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    modification: str
        Modification for which you search (which is the key to any of the ModificationType's stored in the MODIFICATIONS dictionary)
    x_pos: int
        Position of amino acid X in query that is under investigation for conservation
    x_range: int
        Error margin for conservation of x_pos
    evalue: float
        BLAST parameter, e-value for BLAST run
    num_threads: int (default: 1)
        BLAST parameter, number of threads that can be used by BLAST

    """
    # Get BLASTDB name for selected modification + get a temporary path for output
    BLASTDB = db_setup.get_blastdb_name_for_modification(modification)
    BLAST_OUT = "temp.xml"

    # Adjust working directory conditions
    os.chdir(get_data_dir())

    logging.info(f"Running BLAST search for {input} against local {modification} BLAST database.")
    # Run BLAST
    blast_exec = NcbiblastpCommandline(
        query = input,
        db = BLASTDB,
        evalue = evalue,
        outfmt = 5,
        out = BLAST_OUT,
        num_threads = num_threads,
    )
    blast_exec()

    with open(BLAST_OUT) as handle:
        blast_records = list(NCBIXML.parse(handle))

    logging.info(f"Filtering results of BLAST search for {input} against local {modification} BLAST database.")
    return [_filter_blast(i, amino_acid_x, x_pos, x_range, evalue) for i in blast_records]


def _filter_blast(blast_record, amino_acid_x, x_pos, x_range, evalue) -> Blast:
    """
    This function filters the BLAST results.
    First, it filters out any BLAST results where the alignment does not contain:
    - the queried x_pos (in the protein query)
    - the modification position (in the aligned protein)
    Then, it filters out results where the queried modified amino acid X does not align with the modified amino acid X in the aligned protein.
    Ultimately, it only returns a BLAST record containing conserved (within range) protein modifications for similar proteins.

    Parameters
    ----------
    blast_record: Blast
        Blast record containing all similar proteins to the queried one, that are in the specific modification database
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    x_pos: int
        Position of amino acid X in query that is under investigation for conservation
    x_range: int
        Error margin for conservation of x_pos
    evalue: float
        BLAST parameter, e-value for BLAST run

    """
    # Create new Blast Record where we append filtered matches.
    filtered = Blast()

    for a in blast_record.alignments:
        # Parse FASTA title where post-translational modification info is stored
        mod = ModificationHeader.parse(a.title)

        # Append matching High Scoring partners here, which will then be added to the 'filtered' BLAST frame
        filter1_hsps = [] ## Filter1: filters out all hsps which do not contain the modification (both in query and hit)
        filter2_hsps = [] ## Filter2: filters out hsps that do not contain CONSERVED modification

        for hsp in a.hsps:
            if hsp.expect < evalue and _is_modHit_in_alignment(hsp, mod.position) and _is_modQuery_in_alignment(hsp, x_pos):
                # WEE! we have a match.
                filter1_hsps.append(hsp)

        for hsp in filter1_hsps:
            # To assess whether a hsp contains a conserved modification, we need to
            # (1) find the location of the query modification in the aligned query
            if hsp.query.find("-") == -1:
            # (2) find out if the aligned position (+- range) in the hit is amino acid X
                if len(_findXs_in_alignedHit(hsp, amino_acid_x, x_pos, x_range)) != 0:
            # (3) if this aligned position is amino acid X, was this the amino acid X carrying the modification
                    _add_conservedModX_to_listConsHsp(hsp, amino_acid_x, x_pos, x_range, mod, filter2_hsps)
            # (1) find the location of the query modification in the aligned query
            elif (hsp.query_start + hsp.query.find("-") + 1) > x_pos:
            # (2) find out if the aligned position (+- range) in the hit is amino acid X
                if len(_findXs_in_alignedHit(hsp, amino_acid_x, x_pos, x_range)) != 0:
            # (3) if this aligned position is amino acid X, was this the amino acid X carrying the modification
                    _add_conservedModX_to_listConsHsp(hsp, amino_acid_x, x_pos, x_range, mod, filter2_hsps)
            # (1) find the location of the query modification in the aligned query
            else:
            #    should adapt amino acid X position here to match number of gaps before
                countGapBefore = hsp.query[0:x_pos+1].count("-")
                newSeq = hsp.query[0:x_pos+1].replace("-", "") + hsp.query[x_pos+1:len(hsp.query)]
                while newSeq[0:x_pos+1].find('-') != -1:
                    newSeq = newSeq[0:x_pos+1].replace("-", "") + newSeq[x_pos+1:len(newSeq)]
                    countGapBefore += 1
            # (2) find out if the aligned position (+- range) in the hit is amino acid X
                if len(_findXs_in_alignedHit(hsp, amino_acid_x, x_pos + countGapBefore, x_range))  != 0:
            # (3) if this aligned position is amino acid X, was this the amino acid X carrying the modification
                    _add_conservedModX_to_listConsHsp(hsp, amino_acid_x, x_pos + countGapBefore, x_range, mod, filter2_hsps)

        # If some HSPS matched, let's append that to the filtered BLAST frame for future processing.
        if filter2_hsps:
            new_alignment = Alignment()
            new_alignment.title = a.title
            new_alignment.hsps = filter2_hsps
            filtered.alignments.append(new_alignment)

    return filtered

def _is_modHit_in_alignment(hsp, mod_pos) -> bool:
    """
    This function asserts that the aligned hit does contain its modification in the aligned portion of the protein.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    mod_pos: int
        Position of amino acid X in the aligned protein that is known to be modified

    """
    return hsp.sbjct_start <= mod_pos <= hsp.sbjct_end

def _is_modQuery_in_alignment(hsp, query_pos) -> bool:
    """
    This function asserts that the aligned portion of the query protein contains the modification being queried.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    query_pos: int
        Position of amino acid X in query that is under investigation for conservation

    """
    return hsp.query_start <= query_pos <= hsp.query_end

def _findXs_in_alignedHit(hsp, amino_acid_x, x_pos, x_range):
    """
    This function finds the relative positions of specified amino acid X in the neighbourhood of the position of the residue aligned to the amino acid being queried.
    It returns a list of relative positions, all within the x_range.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    x_pos: int
        Position of amino acid X in query that is under investigation for conservation
    x_range: int
        Error margin for conservation of x_pos

    """
    rangeX = []
    for i in range(-x_range, x_range + 1):
        ##need to check that we do not try to access an index out of range for this subject
        if (x_pos - hsp.query_start + i <= len(hsp.sbjct) - 1) and (x_pos - hsp.query_start + i >= 0):
            if hsp.sbjct[x_pos - hsp.query_start + i] == amino_acid_x:
                rangeX.append(i)
    return rangeX

def _add_conservedModX_to_listConsHsp(hsp, amino_acid_x, x_pos, x_range, modification, listHsp):
    """
    This function adds the hsps of modification database entries with conserved modified amino acid X to a list, namely the listHsp.

    Parameters
    ----------
    hsp: hsp
        High Scoring partners, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    x_pos: int
        Position of amino acid X in query that is under investigation for conservation
    x_range: int
        Error margin for conservation of x_pos
    modification: ModificationHeader
        Modification for which you search
    listHsp: list
        List that will be used to append hsps of modification database entries with conserved modified amino acid X to

    """
    for i in _findXs_in_alignedHit(hsp, amino_acid_x, x_pos, x_range):
        indexXhit = x_pos - hsp.query_start + i
        numGapUntilX = hsp.sbjct[0:x_pos - hsp.query_start + i].count("-")
        coordXOriginalSubject = indexXhit - numGapUntilX + hsp.sbjct_start
        if modification.position == coordXOriginalSubject:
            listHsp.append(hsp)
