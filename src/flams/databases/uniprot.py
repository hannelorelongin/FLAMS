#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hannelorelongin
"""

import os
import requests
from pathlib import Path

from Bio import SeqIO

""" uniprot
This script contains the functions to fetch (missing) information from UniProt to create the FLAMS databases.
"""

def fetch_UniProt(uniprot):
    """
    This function gets the FASTA file from UniProt, based on the provided UniProt ID.
    If it is already downloaded, it just provides the filename. If not, it will download the file.
    It returns the filename of the FASTA file.

    Parameters
    ----------
    uniprot: str
        UniProt ID for protein containing the modification

    """
    filename = Path(f"{uniprot}.fasta.tmp")
    if os.path.isfile(filename):
        return filename
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
        sess = requests.Session()
        adapter = requests.adapters.HTTPAdapter(max_retries = 20)
        sess.mount("http://", adapter)
        r = sess.get(url)
        r.raise_for_status()
        filename = Path(f"{uniprot}.fasta.tmp")
        with filename.open("w+") as f:
            f.write(r.text)
        return filename


def fetch_seqRecord(filename):
    """
    This function reads the FASTA file from UniProt, based on the provided filename.
    It returns the SeqRecord object containing all information from this FASTA file.

    Parameters
    ----------
    filename: str
        Filename of UniProt FASTA file for protein containing the modification

    """
    record = SeqIO.read(filename, "fasta")
    return record
