#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hannelorelongin, kasgel
"""

import csv
import logging
import requests
from io import BytesIO, StringIO
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


""" cplmv4
This script downloads the different contents of the CPLM database, and transforms them into a fasta format.
Script developed to work with CPLM database version 4.
"""

URL = "http://cplm.biocuckoo.cn/Download/{0}.zip"

def get_fasta(descriptor, location):
    """
    This function downloads the entries of the CPLM database for a specific modification (according to $descriptor),
    and saves it as a fasta format in $location.

    Parameters
    ----------
    descriptor: str
        Description of a specific modification
    location: str
        Output file

    """
    # HTTP request with stream. This way, we get the size of the file first and can begin downloading it in chunks.
    req = requests.get(URL.format(descriptor), stream = True)

    # Raise an exception if HTTP request failed.
    req.raise_for_status()

    size_in_mb = int(req.headers.get("content-length")) / 1048576

    logging.info(f"Downloading CPMLv4 {descriptor} Database, please wait. Size: {size_in_mb:.1f} MB")

    with ZipFile(BytesIO(req.content)) as myzip:
        # Extract the single txt file and return as UTF-8 string
        plm = myzip.read(myzip.namelist()[0]).decode("UTF-8")
        # SeqIO can not write greek letter beta, immediately change this
        if descriptor == "β-Hydroxybutyrylation":
            plm = plm.replace("β", "beta")

    # actually create fasta file from PTM records
    with open(location, "a", encoding = "UTF-8") as out:
        SeqIO.write(_convert_plm_to_fasta(plm), out, "fasta")

    logging.info(f"Converted and stored CPMLv4 {descriptor} Database entries as FASTA entries for the local {descriptor} BLAST database format.")


def _convert_plm_to_fasta(plm):
    """
    This function converts the string containing all entries of the CPLM database for a specific modification to a fasta format.
    It stores relevant data on the entries in the sequence records.

    Parameters
    ----------
    plm: str
        Content of the text file detailing all entries of the CPLM database for a specific modification

    """
    recs = []
    reader = csv.reader(StringIO(plm), delimiter = "\t")
    # create FLAMS records based on the CPLM information
    for row in reader:
        # properties that can be directly extracted
        seq = Seq(row[6])
        length = len(row[6])
        # properties for which we fix issue with spaces by casting them to underscores
        proteinNoSpaces = f"{row[4]}".replace(" ", "__")
        speciesNoSpaces = f"{row[5]}".replace(" ", "__")
        # actually create record
        id = f"{row[1]}|{row[2]}|{length}|CPLM"
        rec = SeqRecord(
            seq,
            id = id,
            description = f"{proteinNoSpaces}|{row[3]}|{speciesNoSpaces} [{row[0]}|{row[7]}|{row[8]}]",
        )
        recs.append(rec)
    return recs
