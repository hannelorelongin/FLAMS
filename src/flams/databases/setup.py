#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: annkamsk, hannelorelongin, kasgel, MaartenLangen
"""

import subprocess
import logging
import requests
import sys
from io import BytesIO
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List
from zipfile import ZipFile

from Bio import SeqIO

from . import cplmv4, dbptm, scop3p
from ..utils import get_data_dir

""" setup
This script contains our modification database, and all associated functions to generate and maintain BLAST databases for each modification type.
"""

@dataclass
class ModificationDatabase:
    """
    This dataclass contains a list of tuples, containing the information necessary to get the FASTA file of the modification.

    Parameters
    ----------
    module: Any
        Refers to database module, necessary to retrieve the fasta files on the modification described by Descriptor
    descriptor:
        Label to identify the modification

    """
    module: Any
    descriptor: str


@dataclass
class ModificationType:
    """
    This dataclass contains the different types of modifications, which can be identied by type and database version, and contains their modification databases.

    Parameters
    ----------
    type: str
        Label to identify modification
    version: float
        Label to identify CPLM database version
    dbs: List[ModificationDatabase]
        List of modification database

    """
    type: str
    version: float
    dbs: List[ModificationDatabase]
    aas: List[str]


# Here we store a dict of Zenodo URLs that can be queried for.
version_urls = {
    1.0: "https://zenodo.org/records/10143464/files/{0}-{1}.zip?download=1",
    1.1: "https://zenodo.org/records/10171879/files/{0}-{1}.zip?download=1",
    1.2: "https://zenodo.org/records/10958721/files/{0}-{1}.fasta.zip?download=1",
    1.3: "https://zenodo.org/records/14616210/files/{0}-{1}.zip?download=1",
    1.4: "https://zenodo.org/records/16737546/files/{0}-{1}.zip?download=1",
    1.5: "https://zenodo.org/records/18478351/files/{0}-{1}.zip?download=1"
    }

# Here we store a dict of modifications that can be queried for.
    # sorted alphabetically
MODIFICATIONS = {
    "acetylation": ModificationType(
        "acetylation", 1.5,
        [ModificationDatabase(cplmv4, "Acetylation"), ModificationDatabase(dbptm,"Acetylation")],
        ["A","C","D","E","G","K","M","P","R","S","T","V","Y"]
    ),
    "adp-ribosylation": ModificationType(
        "adp-ribosylation", 1.5,
        [ModificationDatabase(dbptm, "ADP-ribosylation")],
        ["C","D","E","G","H","K","N","R","S","Y"]
    ),
    "amidation": ModificationType(
        "amidation", 1.5,
        [ModificationDatabase(dbptm, "Amidation")],
        ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    ),
    "ampylation" : ModificationType(
        "ampylation", 1.5,
        [ModificationDatabase(dbptm, "AMPylation")],
        ["S","T","Y"]
    ),
    "benzoylation": ModificationType(
        "benzoylation", 1.5, [ModificationDatabase(cplmv4, "Benzoylation")],
        ["K"]
    ),
    "beta-hydroxybutyrylation": ModificationType(
        "beta-hydroxybutyrylation", 1.5, [ModificationDatabase(cplmv4, "β-Hydroxybutyrylation")],
        ["K"]
    ),
    "biotinylation": ModificationType(
        "biotinylation", 1.5,
        [ModificationDatabase(cplmv4, "Biotinylation"), ModificationDatabase(dbptm, "Biotinylation")],
        ["K"]
    ),
    "blocked_amino_end": ModificationType(
        "blocked_amino_end", 1.5,
        [ModificationDatabase(dbptm, "Blocked amino end")],
        ["A","C","D","E","G","H","I","L","M","N","P","Q","R","S","T","V"]
    ),
    "butyrylation": ModificationType(
        "butyrylation", 1.5,
        [ModificationDatabase(cplmv4, "Butyrylation"), ModificationDatabase(dbptm, "Butyrylation")],
        ["K"]
    ),
    "carbamidation": ModificationType(
        "carbamidation", 1.5,
        [ModificationDatabase(dbptm, "Carbamidation")],
        ["C"]
    ),
    "carboxyethylation": ModificationType(
        "carboxyethylation", 1.5,
        [ModificationDatabase(cplmv4, "Carboxyethylation"), ModificationDatabase(dbptm, "Carboxyethylation")],
        ["K"]
    ),
    "carboxylation": ModificationType(
        "carboxylation", 1.5,
        [ModificationDatabase(cplmv4, "Carboxylation"), ModificationDatabase(dbptm, "Carboxylation")],
        ["K"]
    ),
    "carboxymethylation": ModificationType(
        "carboxymethylation", 1.5, [ModificationDatabase(cplmv4, "Carboxymethylation")],
        ["K"]
    ),
    "cholesterol_ester": ModificationType(
        "cholesterol_ester", 1.5,
        [ModificationDatabase(dbptm, "Cholesterol ester")],
        ["G"]
    ),
    "citrullination": ModificationType(
        "citrullination", 1.5,
        [ModificationDatabase(dbptm, "Citrullination")],
        ["R"]
    ),
    "crotonylation": ModificationType(
        "crotonylation", 1.5,
        [ModificationDatabase(cplmv4, "Crotonylation"), ModificationDatabase(dbptm, "Crotonylation")],
        ["K"]
    ),
    "c-linked_glycosylation": ModificationType(
        "c-linked_glycosylation", 1.5,
        [ModificationDatabase(dbptm, "C-linked Glycosylation")],
        ["W"]
    ),
    "deamidation": ModificationType(
        "deamidation", 1.5,
        [ModificationDatabase(dbptm, "Deamidation")],
        ["N","Q"]
    ),
    "deamination": ModificationType(
        "deamination", 1.5,
        [ModificationDatabase(dbptm, "Deamination")],
        ["K"]
    ),
    "decanoylation": ModificationType(
        "decanoylation", 1.5,
        [ModificationDatabase(dbptm, "Decanoylation")],
        ["S","T"]
    ),
    "decarboxylation": ModificationType(
        "decarboxylation", 1.5,
        [ModificationDatabase(dbptm, "Decarboxylation")],
        ["D","T"]
    ),
    "dephosphorylation": ModificationType(
        "dephosphorylation", 1.5,
        [ModificationDatabase(dbptm, "Dephosphorylation")],
        ["S","T","Y"]
    ),
    "dietylphosphorylation": ModificationType(
        "dietylphosphorylation", 1.5, [ModificationDatabase(cplmv4, "Dietylphosphorylation")],
        ["K"] 
    ),
    "disulfide_bond": ModificationType(
        "disulfide_bond", 1.5,
        [ModificationDatabase(dbptm, "Disulfide bond")],
        ["C"]
    ),
    "d-glucuronoylation": ModificationType(
        "d-glucuronoylation", 1.5,
        [ModificationDatabase(dbptm, "D-glucuronoylation")],
        ["G"]
    ),
    "farnesylation": ModificationType(
        "farnesylation", 1.5,
        [ModificationDatabase(dbptm, "Farnesylation")],
        ["C"]
    ),
    "formation_of_an_isopeptide_bond": ModificationType(
        "formation_of_an_isopeptide_bond", 1.5,
        [ModificationDatabase(dbptm, "Formation of an isopeptide bond")],
        ["E","Q"]
    ),
    "formylation": ModificationType(
        "formylation", 1.5,
        [ModificationDatabase(cplmv4, "Formylation"), ModificationDatabase(dbptm, "Formylation")],
        ["G","K","M"]
    ),
    "gamma-carboxyglutamic_acid": ModificationType(
        "gamma-carboxyglutamic_acid", 1.5,
        [ModificationDatabase(dbptm, "Gamma-carboxyglutamic acid")],
        ["E"]
    ),
    "geranylgeranylation": ModificationType(
        "geranylgeranylation", 1.5,
        [ModificationDatabase(dbptm, "Geranylgeranylation")],
        ["C"]
    ),
    "glutarylation": ModificationType(
        "glutarylation", 1.5,
        [ModificationDatabase(cplmv4, "Glutarylation"), ModificationDatabase(dbptm, "Glutarylation")],
        ["K"]
    ),
    "glutathionylation": ModificationType(
        "glutathionylation", 1.5,
        [ModificationDatabase(dbptm, "Glutathionylation")],
        ["C"]
    ),
    "glycation": ModificationType(
        "glycation", 1.5, [ModificationDatabase(cplmv4, "Glycation")],
        ["K"]
    ),
    "gpi-anchor": ModificationType(
        "gpi-anchor", 1.5, [ModificationDatabase(dbptm, "GPI-anchor")],
        ["A","C","D","G","N","S","T"]
    ),
    "hmgylation": ModificationType(
        "hmgylation", 1.5, [ModificationDatabase(cplmv4, "HMGylation")],
        ["K"] 
    ),
    "hydroxyceramide_ester": ModificationType(
        "hydroxyceramide_ester", 1.5, [ModificationDatabase(dbptm, "Hydroxyceramide ester")],
        ["Q"]
    ),
    "hydroxylation": ModificationType(
        "hydroxylation", 1.5,
        [ModificationDatabase(cplmv4, "Hydroxylation"), ModificationDatabase(dbptm, "Hydroxylation")],
        ["C","D","E","F","H","I","K","L","N","P","R","S","T","V","W","Y"]
    ),
    "iodination": ModificationType(
        "iodination", 1.5,
        [ModificationDatabase(dbptm, "Iodination")],
        ["Y"]
    ),
    "lactoylation": ModificationType(
        "lactoylation", 1.5,
        [ModificationDatabase(dbptm, "Lactoylation")],
        ["K"]
    ),
    "lactylation": ModificationType(
        "lactylation", 1.5,
        [ModificationDatabase(cplmv4, "Lactylation"), ModificationDatabase(dbptm, "Lactylation")],
        ["K"]
    ),
    "lipoylation": ModificationType(
        "lipoylation", 1.5,
        [ModificationDatabase(cplmv4, "Lipoylation"), ModificationDatabase(dbptm, "Lipoylation")],
        ["K"]
    ),
    "malonylation": ModificationType(
        "malonylation", 1.5,
        [ModificationDatabase(cplmv4, "Malonylation"), ModificationDatabase(dbptm, "Malonylation")],
        ["K"]
    ),
    "methylation": ModificationType(
        "methylation", 1.5,
        [ModificationDatabase(cplmv4, "Methylation"), ModificationDatabase(dbptm, "Methylation")],
        ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","Y"]
    ),
    "mgcylation": ModificationType(
        "mgcylation", 1.5, [ModificationDatabase(cplmv4, "MGcylation")],
        ["K"]
    ),
    "mgylation": ModificationType(
        "mgylation", 1.5, [ModificationDatabase(cplmv4, "MGylation")],
        ["K"]
    ),
    "myristoylation": ModificationType(
        "myristoylation", 1.5, [ModificationDatabase(dbptm, "Myristoylation")],
        ["C","G","K"]
    ),
    "neddylation": ModificationType(
        "neddylation", 1.5,
        [ModificationDatabase(cplmv4, "Neddylation"), ModificationDatabase(dbptm, "Neddylation")],
        ["K"]
    ),
    "nitration": ModificationType(
        "nitration", 1.5, [ModificationDatabase(dbptm, "Nitration")],
        ["Y"]
    ),
    "n-carbamoylation": ModificationType(
        "n-carbamoylation", 1.5, [ModificationDatabase(dbptm, "N-carbamoylation")],
        ["A"]
    ),
    "n-linked_glycosylation": ModificationType(
        "n-linked_glycosylation", 1.5,
        [ModificationDatabase(dbptm, "N-linked Glycosylation")],
        ["D","I","K","N","R","S","T","V","W"]
    ),
    "n-palmitoylation": ModificationType(
        "n-palmitoylation", 1.5, [ModificationDatabase(dbptm, "N-palmitoylation")],
        ["C","G","K"]
    ),
    "octanoylation": ModificationType(
        "octanoylation", 1.5,
        [ModificationDatabase(dbptm, "Octanoylation")],
        ["S","T"]
    ),
    "oxidation": ModificationType(
        "oxidation", 1.5,
        [ModificationDatabase(dbptm, "Oxidation")],
        ["C","L","M","S","W"]
    ),
    "o-linked_glycosylation": ModificationType(
        "o-linked_glycosylation", 1.5,
        [ModificationDatabase(dbptm, "O-linked Glycosylation")],
        ["K","P","S","T","Y"]
    ),
    "o-palmitoleoylation": ModificationType(
        "o-palmitoleoylation", 1.5,
        [ModificationDatabase(dbptm, "O-palmitoleoylation")],
        ["S"]
    ),
    "o-palmitoylation": ModificationType(
        "o-palmitoylation", 1.5,
        [ModificationDatabase(dbptm, "O-palmitoylation")],
        ["S","T"]
    ),
    "phosphatidylethanolamine_amidation": ModificationType(
        "phosphatidylethanolamine_amidation", 1.5, [ModificationDatabase(dbptm, "Phosphatidylethanolamine amidation")],
        ["G"]
    ),
    "phosphoglycerylation": ModificationType(
        "phosphoglycerylation", 1.5,
        [ModificationDatabase(cplmv4, "Phosphoglycerylation")],
        ["K"]
    ),
    "phosphorylation": ModificationType(
        "phosphorylation", 1.5,
        [ModificationDatabase(dbptm, "Phosphorylation"), ModificationDatabase(scop3p, "Phosphorylation")],
        ["A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","W","Y"]
    ),
    "propionylation": ModificationType(
        "propionylation", 1.5,
        [ModificationDatabase(cplmv4, "Propionylation"), ModificationDatabase(dbptm, "Propionylation")],
        ["K"]
    ),
    "pupylation": ModificationType(
        "pupylation", 1.5, [ModificationDatabase(cplmv4, "Pupylation")],
        ["K"]
    ),
    "pyrrolidone_carboxylic_acid": ModificationType(
        "pyrrolidone_carboxylic_acid", 1.5, [ModificationDatabase(dbptm, "Pyrrolidone carboxylic acid")],
        ["E","Q"]
    ),
    "pyrrolylation": ModificationType(
        "pyrrolylation", 1.5, [ModificationDatabase(dbptm, "Pyrrolylation")],
        ["C"]
    ),
    "pyruvate": ModificationType(
        "pyruvate", 1.5, [ModificationDatabase(dbptm, "Pyruvate")],
        ["C","S"]
    ),
    "serotonylation" : ModificationType(
        "serotonylation", 1.5,
        [ModificationDatabase(dbptm, "Serotonylation")],
        ["Q"]
    ),
    "stearoylation" : ModificationType(
        "stearoylation", 1.5,
        [ModificationDatabase(dbptm, "Stearoylation")],
        ["C"]
    ),
    "succinylation": ModificationType(
        "succinylation", 1.5,
        [ModificationDatabase(cplmv4, "Succinylation"), ModificationDatabase(dbptm, "Succinylation")],
        ["C","K","W"]
    ),
    "sulfation" : ModificationType(
        "sulfation", 1.5,
        [ModificationDatabase(dbptm, "Sulfation")],
        ["C","S","T","Y"]
    ),
    "sulfhydration" : ModificationType(
        "sulfhydration", 1.5,
        [ModificationDatabase(dbptm, "Sulfhydration")],
        ["C"]
    ),
    "sulfoxidation" : ModificationType(
        "sulfoxidation", 1.5,
        [ModificationDatabase(dbptm, "Sulfoxidation")],
        ["M"]
    ),
    "sumoylation": ModificationType(
        "sumoylation", 1.5,
        [ModificationDatabase(cplmv4, "Sumoylation"), ModificationDatabase(dbptm, "Sumoylation")],
        ["K"]
    ),
    "s-archaeol" : ModificationType(
        "s-archaeol", 1.5,
        [ModificationDatabase(dbptm, "S-archaeol")],
        ["C"]
    ),
    "s-carbamoylation" : ModificationType(
        "s-carbamoylation", 1.5,
        [ModificationDatabase(dbptm, "S-carbamoylation")],
        ["C"]
    ),
    "s-cyanation" : ModificationType(
        "s-cyanation", 1.5,
        [ModificationDatabase(dbptm, "S-Cyanation")],
        ["C"]
    ),
    "s-cysteinylation" : ModificationType(
        "s-cysteinylation", 1.5,
        [ModificationDatabase(dbptm, "S-cysteinylation")],
        ["C"]
    ),
    "s-diacylglycerol" : ModificationType(
        "s-diacylglycerol", 1.5,
        [ModificationDatabase(dbptm, "S-diacylglycerol")],
        ["C"]
    ),
    "s-linked_glycosylation": ModificationType(
        "s-linked_glycosylation", 1.5,
        [ModificationDatabase(dbptm, "S-linked Glycosylation")],
        ["C"]
    ),
    "s-nitrosylation" : ModificationType(
        "s-nitrosylation", 1.5,
        [ModificationDatabase(dbptm, "S-nitrosylation")],
        ["C"]
    ),
    "s-palmitoylation" : ModificationType(
        "s-palmitoylation", 1.5,
        [ModificationDatabase(dbptm, "S-palmitoylation")],
        ["C"]
    ),
    "thiocarboxylation" : ModificationType(
        "thiocarboxylation", 1.5,
        [ModificationDatabase(dbptm, "Thiocarboxylation")],
        ["G"]
    ),
    "ubiquitination": ModificationType(
        "ubiquitination", 1.5,
        [ModificationDatabase(cplmv4, "Ubiquitination"), ModificationDatabase(dbptm, "Ubiquitination")],
        ["C","K","R","S"]
    ),
    "umpylation" : ModificationType(
        "umpylation", 1.5,
        [ModificationDatabase(dbptm, "UMPylation")],
        ["S","T","Y"]
    ),
    "2-hydroxyisobutyrylation": ModificationType(
        "2-hydroxyisobutyrylation", 1.5, [ModificationDatabase(cplmv4, "2-Hydroxyisobutyrylation")],
        ["K"]
    ),
    }


def update_db_for_modifications(list_of_mods_to_check: List[str]):
    """
    This function updates the local BLASTDB for a given modification.

    Parameters
    ----------
    list_of_mods_to_check: List[str]
        List of modifications (which are keys to any of the ModificationType's stored in the MODIFICATIONS dictionary),
        for which the database should be updated.

    """
    for m in list_of_mods_to_check:
        _generate_blastdb_if_not_up_to_date(MODIFICATIONS[m])


def _generate_blastdb_if_not_up_to_date(modification: ModificationType):
    """
    This function generates a new local BLASTDB when a newer database version for a specific modification is available.

    Parameters
    ----------
    modification: ModificationType
        ModificationType for which a BLASTDB will be generated

    """
    data_dir = get_data_dir()

    BLASTDB_PATH = get_blastdb_name_for_modification(
        modification.type, modification.version
    )

    # If an up-to-date BLASTDB for the given modification already exists, do nothing.
    if Path(f"{data_dir}/{BLASTDB_PATH}.pdb").exists():
        return

    # If no up-to-date BLASTDB exists, check whether a FASTA file exists with the modification data
    fasta_location_mod = f"{data_dir}/{modification.type}-{modification.version}.fasta"

    if not Path(fasta_location_mod).exists():
        # If this does not exist, create a FASTA file with the modification data
        try:
            _get_fasta_from_zenodo(modification, fasta_location_mod)
        except requests.HTTPError:
            logging.error(f"Could not fetch FLAMS {modification.type} Database {modification.version} from Zenodo. Please try again later. Exiting FLAMS...")
            logging.error("If needed, you can also download the databases from CPLM and dbPTM yourself, see the docs. (not recommended: very slow!)")
            sys.exit()

# FOR OWN DATABASE DOWNLOAD - on a fresh install:
# Comment out try/except above (lines 515-520), uncomment the following line of code
#        _get_fasta_for_blast(modification, data_dir, fasta_location_mod)

    # Generate local BLASTDB from FASTA in fasta_location_mod
    _generate_blastdb(data_dir, modification)


def _generate_blastdb(data_dir, modification: ModificationType):
    """
    This function generates a local BLASTDB for a given modification.

    Parameters
    ----------
    data_dir: directory
        Platform-specific directory that stores app data. The local BLAST database will be stored here.
    modification: ModificationType
        ModificationType for which a local BLAST database will be generated

    """
    try:
        # We presume that the FASTA is stored in a file {modification.type}.fasta inside the data_dir.
        # We will write the local BLASTDB to out_path

        out_db_name = get_blastdb_name_for_modification(
            modification.type, modification.version
        )
        
        subprocess.run(["makeblastdb", "-in", f"{modification.type}-{modification.version}.fasta", "-dbtype", "prot", 
                        "-input_type", "fasta", "-parse_seqids", "-out", str(out_db_name)],
                        cwd = str(data_dir), check = True)

    except FileNotFoundError as e:
        raise FileNotFoundError(
            "You need to install BLAST and include it in system PATH."
        ) from e


def get_blastdb_name_for_modification(modification: str, version = None):
    """
    This function gets the name of the local BLASTDB for a given modification.

    Parameters
    ----------
    modification: str
        Description of a specific modification,
        must be the key to any of the ModificationType's stored in the MODIFICATIONS dictionary.
    version: float
        Database version

    """
    # If version was not specified, get the current
    if not version:
        version = MODIFICATIONS[modification].version

    return f"{modification}-{version}"


def _get_fasta_from_zenodo(modification, fasta_location_mod):
    """
    This function downloads the .zip file containing the fasta file with information from all databases containing information on the specified modification.
    It will create and store a .fasta file containing all entries of all databases related to this to modification file in the fasta_location.

    Parameters
    ----------
    modification: ModificationType
        ModificationType for which a fasta file will be generated
    fasta_location_mod:
        Output .fasta file containing all modification entries of all databases for the specified modification

    """
    URL = version_urls.get(modification.version)
    req = requests.get(URL.format(modification.type, modification.version), stream = True)
    # Raise an exception if HTTP request failed.
    req.raise_for_status()
    size_in_mb = int(req.headers.get("content-length")) / 1048576
    logging.info(f"Downloading FLAMS {modification.type} Database {modification.version}, please wait. Size: {size_in_mb:.1f} MB")

    with ZipFile(BytesIO(req.content)) as myzip:
        # Extract the single txt file and return as UTF-8 string
        ptm = myzip.read(myzip.namelist()[0]).decode("UTF-8")

    filename = Path(fasta_location_mod)
    with filename.open("w+", encoding = "UTF-8") as f:
        f.write(ptm)


#############################################
# Optional: downloading your own databases #
#############################################

def _get_fasta_for_blast(modification: ModificationType, data_dir, fasta_location_mod):
    """
    This function creates a fasta file combining the information from all databases containing information on the specified modification.
    It will create and store a fasta file containing all entries of all databases related to this to modification file in the fasta_location.

    Parameters
    ----------
    modification: ModificationType
        ModificationType for which a fasta file will be generated
    data_dir: directory
        Platform-specific directory that stores app data. The local BLAST database will be stored here.
    fasta_location_mod:
        Output .fasta file containing all modification entries of all databases for the specified modification

    """
    # create fasta
    for db in modification.dbs:
        module_name = db.module.__name__.split(".")[-1]
        fasta_location_per_database = f"{data_dir}/{modification.type}-{modification.version}-{module_name}.fasta"
        file_before_deduplication = f"{fasta_location_mod.removesuffix('.fasta')}-before_deduplication.fasta"
        if not Path(fasta_location_per_database).exists():
            _get_fasta_from_dbs(db, fasta_location_per_database)
        fileMod = open(file_before_deduplication, "a")
        fileDb = open(fasta_location_per_database, "r")
        fileMod.write(fileDb.read())
        fileMod.close()
        fileDb.close()

    # remove duplicate entries, that will create errors when creating BLAST DB
        # should only happen in rare cases, dbPTM and CPLM try to avoid duplicates, but they do happen
    seen = set()
    records = []
    for record in SeqIO.parse(file_before_deduplication, "fasta"):
        if record.id not in seen:
            seen.add(record.id)
            records.append(record)
        else:
            logging.warning(f"Your database contained a duplicate entry for {record.id}, only the first one was added to your BLAST database.")
    SeqIO.write(records, fasta_location_mod, "fasta")


def _get_fasta_from_dbs(db, fasta_location):
    """
    This function calls on the get_fasta function from a specific module for a modification.
    It will create and store a fasta file containing all entries in the specified database related to this to modification file in the fasta_location.

    Parameters
    ----------
    db:
        PTM database for which all entries will be converted to a fasta file (either dbptm, cplmv4 or scop3p)
    fasta_location:
        Output .fasta file containing all modification entries of the specified database for the specified modification

    """
    db.module.get_fasta(db.descriptor, fasta_location)
    
