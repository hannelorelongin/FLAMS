#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hannelorelongin
"""

import logging
import os
import re
import requests
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .uniprot import fetch_seqRecord, fetch_UniProt

""" scop3p
This script downloads the contents of the SCOP3P database, and transforms them into a fasta format.
Script developed to work with SCOP3 database versions from 2026 on.
"""

URL = "https://iomics.ugent.be/scop3p/api/get-all-modifications"

def get_fasta(descriptor, location):
    """
    This function downloads the entries of the SCOP3P database,
    and saves it as a fasta format in $location.

    Parameters
    ----------
    location: str
        Output file

    """
    # use shared SCOP3P modification file if it exists in the working directory
    if os.path.exists("Scop3P_Modification_2026_01.txt"):
        logging.info("Parsing SCOP3P Database file, please wait.")
        mod_list = _modlist_from_file("Scop3P_Modification_2026_01.txt")
    # if not, use API call to fetch data
    else:
        # HTTP request with stream. This way, we get the size of the file first and can begin downloading it in chunks.
        req = requests.get(URL, stream = True)
        # Raise an exception if HTTP request failed.
        req.raise_for_status()
        logging.info("Downloading SCOP3P Database, please wait.")
        # extract list of all modifications
        mod_list = req.json().get("modifications")

    # filter the modifications to retain only those based on experimental evidence
    exp_mod_list = _filter_modifications(mod_list)

    # actually create fasta file from PTM records
    with open(location, "a", encoding = "UTF-8") as out:
        SeqIO.write(_convert_pm_list_to_fasta(exp_mod_list), out, "fasta")

    logging.info("Converted and stored SCOP3P Database entries as FASTA entries for the local phosphorylation BLAST database format.")

def _modlist_from_file(scop3p_file):
    """
    This function creates a pm_list from a SCOP3P modification file stored locally, 
    in such a way that it replicates the necessary properties of protein modification data fetched from the SCOP3P API.

    Parameters
    ----------
    scop3p_file: Path
        Path to the SCOP3P modification file

    """
    # dictionaries to help parse the SCOP3P file
    # important note: confusingly, the fields source and reference are swapped in the SCOP3P file compared to the API
    selected_keys = {0: "uniprotId", 1: "position", 3: "reference", 4: "evidence", 5: "source"}
    defaults = {"name": "phosphorylation"}
    diff_api = {"": None, "Combinatorial": "Combined"}

    # empty list to store modifications
    pm_list = [] 

    # actually parsing the provided modification file
    with open(scop3p_file, "r", encoding = "utf-8") as in_file: 
        # skip first line containing headers
        next(in_file) 
        # go through contents
        for line in in_file: 
            parts = line.strip("\n").split("\t")
            pm = {}
            for col_index, key in selected_keys.items(): 
                # making sure we deal with subtle differences in API vs file content
                if parts[col_index] in diff_api.keys():
                    value = diff_api.get(parts[col_index])
                else: 
                    value = parts[col_index]
                pm[key] = value 
            # Add default fields 
            for k, v in defaults.items(): 
                pm.setdefault(k, v) 
            pm_list.append(pm)

    return pm_list


def _filter_modifications(pm_list):
    """
    This function filters the list containing all entries of the SCOP3P database, to retain only those with experimental evidence.
    Specifically, this means storing (i) those that are derived from the PRIDE database, 
    (ii) those that are derived from UniProt with experimental evidence.
    It returns a filtered list, containing only SCOP3P entries with experimental evidence.

    Parameters
    ----------
    pm_list: list
        Content of the SCOP3P database, as fetched from URL

    """
    filtered_pm_list = list()

    for pm in pm_list:
        # extract source of PTM data
        source = pm.get("reference")
        # keep those from PRIDE
        if "PRIDE" in source:
            filtered_pm_list.append(pm)
        # for those listing UniProt    
        elif "UP" in source:
            # extract evidence
            evidence = pm.get("evidence")
            # keep if experimental evidence (evidence code is "Experimental" or "Combined")
            if (evidence == "Experimental") or (evidence == "Combined"):
                filtered_pm_list.append(pm)
    
    return filtered_pm_list

def _convert_pm_list_to_fasta(pm_list):
    """
    This function converts the list containing all entries of the SCOP3P database to a fasta format.
    It stores relevant data on the entries in the sequence records.

    Parameters
    ----------
    pm_list: list of dictionaries
        Content of the SCOP3P database, as fetched from URL and filtered

    """
    recs = []
    # create FLAMS records based on the SCOP3P information
    for pm in pm_list:
        # properties that can be directly extracted
        uniprot = pm.get("uniprotId")
        pubmed = pm.get("source")
        # conform the pubmed list to adhere to same format as dbPTM/CPLM records
        if pubmed != None:
            pubmed = pm.get("source").replace("PubMed:", "").replace(",", ";")
        evidence = pm.get("evidence")
        # deal with the fact that SCOP3P uses None value for evidence if record comes from PRIDE
        if evidence == None:
            evidence = "Experimental"
        # conform the pubmed list to adhere to same format as dbPTM/CPLM records
        evidence = evidence.replace("Experimental", "Exp.").replace("Combined", "Comb.")
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
            logging.warning(f"UniProt FASTA file for protein with UniProt ID {uniprot} does not contain a field GN=... Hence, no name could be assigned to the protein.")
        # properties for which we fix issue with spaces by casting them to underscores
        proteinNoSpaces = protein.replace(" ", "__")
        speciesNoSpaces = species.replace(" ", "__")
        modificationNoSpaces =  pm.get("name").capitalize().replace(" ","__")
        # actually create record
        id = f"{uniprot}|{pm.get('position')}|{length}|SCOP3P"
        rec = SeqRecord(
            seq,
            id = id,
            description = f"{proteinNoSpaces}|{modificationNoSpaces}|{speciesNoSpaces} [SCOP3P|{evidence}|{pubmed}]",
        )
        recs.append(rec)
        logging.info("added rec.")

    # delete temporary files fetched from UniProt
    list_tmp_files = os.listdir(Path().absolute())
    for item in list_tmp_files:
        if item.endswith(".fasta.tmp"):
            os.remove(item)
    return recs
