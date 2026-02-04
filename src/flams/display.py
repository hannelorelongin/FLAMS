#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, Retro212
"""

import csv
import logging
import os
import pandas as pd

""" setup
This script deals with returning the results of FLAMS to the user in a .tsv file.
"""

def create_output(output_filename, amino_acid_x, blast_records, len):
    """
    This function creates the final output .tsv file containing all conserved modification sites, based on a specific FLAMS run.
    It first creates the temporary file with duplicates, before creating a deduplicated output file.
    Then, it deletes the temporary file with duplicates.

    Parameters
    ----------
    output_filename: str
        Output file name
    blast_records: array
        Array containing BLAST records that met search criteria of FLAMS run.
    len: int
        Length of query protein (to calculate query coverage)

    """
    output_pre_dedupl = f"{output_filename}.tmp"
    _display_result(output_pre_dedupl, amino_acid_x, blast_records, len)
    _deduplicate_output(amino_acid_x, output_pre_dedupl, output_filename)
    os.remove(output_pre_dedupl)

def _display_result(output_filename, amino_acid_x, blast_records, len):
    """
    This function creates a .tsv file containing all conserved modification sites, based on a specific FLAMS run.

    Parameters
    ----------
    output_filename: str
        Output file name
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    blast_records: array
        Array containing BLAST records that met search criteria of FLAMS run.
    len: int
        Length of query protein (to calculate query coverage)

    """
    logging.info(f"Writing .tsv output file with all conserved {amino_acid_x} modifications.")
    with open(output_filename, "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter = "\t")
        tsv_writer.writerow(
            [
                "Uniprot ID",
                "Protein name",
                "Modification",
                f"{amino_acid_x} location",
                f"{amino_acid_x} window",
                "Species",
                "BLAST E-value",
                "BLAST identity",
                "BLAST coverage",
                "CPLM ID",
                "CPLM evidence code",
                "CPLM evidence links",
                "dbPTM evidence code",
                "dbPTM evidence links",
                "SCOP3P evidence code",
                "SCOP3P evidence links"
            ]
        )
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Parsing header of format UniProtID|proteinName|xPosition|db modificationType|speciesNoSpaces [db_ID|evidenceCode|evidenceLink]
                    headerSplitSpace = (alignment.title).split()  # split up header into list, seperated by space
                    generalDescr1 = headerSplitSpace[0]
                    generalDescr2 = headerSplitSpace[1]
                    dbDescr = headerSplitSpace[2]
                    # Split generalDescr1
                    uniprot_id = generalDescr1.split("|")[0]
                    x_location = int(generalDescr1.split("|")[1])
                    protein_length = int(generalDescr1.split("|")[2])
                    # Split generalDescr2
                    protein_name = generalDescr2.split("|")[0].replace("__"," ")
                    modification_type = generalDescr2.split("|")[1].replace("__"," ")
                    species = generalDescr2.split("|")[2].replace("__"," ")
                    # Split dbDescr
                    if generalDescr1.split("|")[3] == "CPLM":
                        cplm_id = dbDescr.split("|")[0][1:]
                        cplm_evidenceCode = dbDescr.split("|")[1]
                        cplm_evidenceLink = dbDescr.split("|")[2][:-1]
                        dbptm_evidenceCode = "NA"
                        dbptm_evidenceLink = "NA"
                        scop3p_evidenceCode = "NA"
                        scop3p_evidenceLink = "NA"
                    elif generalDescr1.split("|")[3] == "dbPTM":
                        cplm_id = "NA"
                        cplm_evidenceCode = "NA"
                        cplm_evidenceLink = "NA"
                        dbptm_evidenceCode = dbDescr.split("|")[1]
                        dbptm_evidenceLink = dbDescr.split("|")[2][:-1]
                        scop3p_evidenceCode = "NA"
                        scop3p_evidenceLink = "NA"
                    elif generalDescr1.split("|")[3] == "SCOP3P":
                        cplm_id = "NA"
                        cplm_evidenceCode = "NA"
                        cplm_evidenceLink = "NA"
                        dbptm_evidenceCode = "NA"
                        dbptm_evidenceLink = "NA"
                        scop3p_evidenceCode = dbDescr.split("|")[1]
                        scop3p_evidenceLink = dbDescr.split("|")[2][:-1]

                    # BLAST properties
                    eval = hsp.expect
                    percentageIdentity = "{:.2%}".format(hsp.identities/hsp.align_length)
                    percentageCoverage = "{:.2%}".format((hsp.query_end-hsp.query_start+1)/len)

                    # Write output
                    tsv_writer.writerow(
                        [
                            uniprot_id,
                            protein_name,
                            modification_type,
                            x_location,
                            _getSequenceWindow(hsp, x_location),
                            species,
                            eval,
                            percentageIdentity,
                            percentageCoverage,
                            cplm_id,
                            cplm_evidenceCode,
                            cplm_evidenceLink,
                            dbptm_evidenceCode,
                            dbptm_evidenceLink,
                            scop3p_evidenceCode,
                            scop3p_evidenceLink
                        ]
                    )


def _deduplicate_output(amino_acid_x, output_pre_dedupl, output_filename):
    """
    This function creates a .tsv file containing all conserved modification sites, based on a specific FLAMS run, without duplicates.
    In this case, duplicates refer to identical hits, found once in CPLM and once in dbPTM.
    The created .tsv file merges the info from the two rows with the identical hits.

    Parameters
    ----------
    amino_acid_x: str
        Amino acid containing the post-translational modification under investigation
    output_pre_dedupl: str
        Output file name of the file before deduplication
    output_filename: str
        Output file name

    """
    df = pd.read_table(output_pre_dedupl, dtype = {"dbPTM evidence links": str, "CPLM evidence links": str,  "SCOP3P evidence links": str})
    df["Protein name"] = df["Protein name"].fillna("not found")
    df2 = df.groupby(["Uniprot ID", "Modification", f"{amino_acid_x} location", f"{amino_acid_x} window", "Species",
                        "BLAST E-value", "BLAST identity", "BLAST coverage"]).agg({"Protein name" : "first",
                        "CPLM ID" : "first", "CPLM evidence code" : "first", "CPLM evidence links" : "first",
                        "dbPTM evidence code" : "first", "dbPTM evidence links" : "first",
                        "SCOP3P evidence code" : "first", "SCOP3P evidence links": "first"}).reset_index()
    df2 = df2.reindex(["Uniprot ID", "Protein name", "Modification", f"{amino_acid_x} location", f"{amino_acid_x} window", "Species",
    "BLAST E-value", "BLAST identity", "BLAST coverage",
    "CPLM ID", "CPLM evidence code", "CPLM evidence links",
    "dbPTM evidence code", "dbPTM evidence links",
    "SCOP3P evidence code", "SCOP3P evidence links"], axis = 1)
    df2.to_csv(output_filename, sep = "\t", index = False)


def _getSequenceWindow(hsp, x_location):
    """
    This function generates the sequence window around the modified amino acid X.
    If the modified amino acid X is not near the end (neither in the query nor in the aligned sequence),
        it simply returns the window containing the 5 amino acids before and after the modified amino acid X.
    However, if the modified amino acid X is near either the start or the end of the aligned sequence, the sequence window can only contain part of this window,
        and this function makes sure this limit is respected.

    Parameters
    ----------
    hsp: hsp
        High Scoring partner, contains information on the alignment between the query protein and one of the aligned entries of the modification database
    x_location: int
        Position of amino acid X in the aligned protein that is known to be modified

    """
    sequence = hsp.sbjct.replace("-", "")
    protSize = len(sequence)
    modPos = x_location - hsp.sbjct_start
    xWindowMax = modPos + 6
    xWindowMin = modPos - 5
    if modPos + 6 > protSize:
        xWindowMax = protSize
    if modPos - 6 < 0:
        xWindowMin = 0
    windowString = (str(xWindowMin+hsp.sbjct_start) + "-" + sequence[xWindowMin:xWindowMax] + "-" + str(xWindowMax+hsp.sbjct_start-1))
    return windowString
