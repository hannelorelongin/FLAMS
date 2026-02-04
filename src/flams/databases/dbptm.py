#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hannelorelongin
"""

import csv
import logging
import re
import requests
import os
from io import BytesIO, StringIO 
from pathlib import Path
from zipfile import ZipFile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .uniprot import fetch_seqRecord, fetch_UniProt


""" dbptm
This script downloads the different contents of the dbPTM database, and transforms them into a fasta format.
Script developed to work with dbPTM database versions 2023-2026.
"""

# This URL refers to the dbPTM databases before the 2025 release
# URL = "https://awi.cuhk.edu.cn/dbPTM/download/experiment/{0}.zip"
# The following URL refers to the dbPTM databases since the 2025 update
URL = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/{0}.zip"

def get_fasta(descriptor, location):
    """
    This function downloads the entries of the dbPTM database for a specific modification (according to $descriptor),
    and saves it as a fasta format in $location.

    Parameters
    ----------
    descriptor: str
        Description of a specific modification
    location: str
        Output file

    """
    # Accomodate old TSL v1.2 protocol by adding cipher
    requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS += ':RSA+AESGCM'
    try:
        requests.packages.urllib3.contrib.pyopenssl.DEFAULT_SSL_CIPHER_LIST += ':RSA+AESGCM'
    except AttributeError:
    # no pyopenssl support used / needed / available
        pass

    # HTTP request with stream. This way, we get the size of the file first and can begin downloading it in chunks.
    req = requests.get(URL.format(descriptor), stream = True)

    # Raise an exception if HTTP request failed.
    req.raise_for_status()

    size_in_mb = int(req.headers.get("content-length")) / 1048576

    logging.info(f"Downloading dbPTM {descriptor} Database, please wait. Size: {size_in_mb:.1f} MB")

    with ZipFile(BytesIO(req.content)) as myzip:
        # Extract the single txt file and return as UTF-8 string
        pm = myzip.read(myzip.namelist()[0]).decode("UTF-8")

    # actually create fasta file from PTM records
    with open(location, "a", encoding = "UTF-8") as out:
        SeqIO.write(_convert_dbptm_to_fasta(pm), out, "fasta")

    logging.info(f"Converted and stored dbPTM {descriptor} Database entries as FASTA entries for the local {descriptor} BLAST database format.")


def _convert_dbptm_to_fasta(pm):
    """
    This function converts the string containing all entries of the dbPTM database for a specific modification to a fasta format.
    It stores relevant data on the entries in the sequence records.

    Parameters
    ----------
    pm: str
        Content of the text file detailing all entries of the dbPTM database for a specific modification

    """
    recs = []
    reader = csv.reader(StringIO(pm), delimiter = "\t")
    # create FLAMS records based on the dbPTM information
    for row in reader:
        # properties that can be directly extracted
        uniprot = row[1].strip()
        # based on UniProt identifier, fetch the relevant information form UniProt
        try:
            record = fetch_seqRecord(fetch_UniProt(uniprot))
        # deal with with UniProt FASTA file issues
        except requests.HTTPError:
            # happens for isoforms
            logging.warning(f"UniProt FASTA file for protein with UniProt ID {uniprot} does not exist on UniProt. Hence, this protein was not added to the database.")
            continue
        except ValueError:
            # happens when UniProt ID is obsolete
            logging.warning(f"UniProt FASTA file for protein with UniProt ID {uniprot} is empty. Hence, this protein was not added to the database.")
            continue
        # properties that can be extracted from UniProt
        seq = record.seq
        length = len(seq)
        species = re.search(r"OS=(.*?) OX=", record.description).group(1)
        # deal with FASTA files missing GN= field (no ORF name or gene name defined)
        try:
            protein = re.search(r"GN=(.*?) PE=", record.description).group(1)
        except AttributeError:
            protein = "NA"
            logging.warning(f"UniProt FASTA file for protein with UniProt ID {row[1]} does not contain a field GN=... Hence, no name could be assigned to the protein.")
        # properties for which we fix issue with spaces by casting them to underscores
        proteinNoSpaces = protein.replace(" ", "__")
        speciesNoSpaces = species.replace(" ", "__")
        modificationNoSpaces =  row[3].replace(" ", "__")
        # actually create record
        id = f"{uniprot}|{row[2]}|{length}|dbPTM"
        rec = SeqRecord(
            seq,
            id = id,
            description = f"{proteinNoSpaces}|{modificationNoSpaces}|{speciesNoSpaces} [dbPTM|Exp.|{row[4]}]",
        )
        recs.append(rec)
        logging.info("added rec.")

    # delete temporary files fetched from UniProt
    list_tmp_files = os.listdir(Path().absolute())
    for item in list_tmp_files:
        if item.endswith(".fasta.tmp"):
            os.remove(item)

    return recs



