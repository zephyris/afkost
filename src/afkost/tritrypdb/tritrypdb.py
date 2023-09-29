import os
import requests
from functools import cached_property

class TriTrypDB:
    def __init__(self, version: str = None):
        """
        Initialises the tritrypdb class for fasta data acces, including cache directory for data storage.

        Named arguments:
        version -- tritrypdb version (default None, which gives latest version)
        """
        self.cache_path = "_tritrypdb_cache"
        if not os.path.isdir(self.cache_path):
            os.mkdir(self.cache_path)
        self._version = None
    
    @cached_property
    def version(self, _version = None):
        """
        Interprets the tritrypdb version to use from `self._version`
        """
        if self._version is None:
            # get latst release version
            url = "https://tritrypdb.org/common/downloads/Current_Release/Build_number"
            try:
                r = requests.get(url)
                return r.text.strip()
            except requests.exceptions.RequestException as e:
                raise SystemExit(e)
        else:
            return str(self._version)

    def fetch_fasta(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: species/strain name, as found on tritrypdb
        """
        if not os.path.isfile(species + ".fasta"):
            # construct url
            url = "https://tritrypdb.org/common/downloads/release-" + self.version + "/" + species + "/fasta/data/TriTrypDB-" + self.version + "_" + species + "_AnnotatedProteins.fasta"
            # download and save
            try:
                r = requests.get(url)
                if r.status_code == 404:
                    raise ValueError("Data not found for species: " + species + ", version: " + self.version + " at tritrypdb.org")
                with open(self.fasta_path(species), "w") as fasta_file:
                    fasta_file.write(r.text)
            except requests.exceptions.RequestException as e:
                raise SystemExit(e)
    
    def fasta_path(self, species: str):
        """
        Path to a fasta file for a given species.

        Required arguments:
        species: species/strain name, as found on tritrypdb
        """
        self.fetch_fasta(species)
        return os.path.join(self.cache_path, species + "." + self.version +".fasta")
    
    def sequences(self, species: str):
        """
        Returns sequences for a tritrypdb species as a Fasta instance

        Required arguments:
        species: species/strain name, as found on tritrypdb

        Returns:
        `Fasta` instance containing the sequences for that species
        """
        from afkost import Fasta
        self.fetch_fasta(species)
        fasta = Fasta(os.path.join(self.cache_path, species + "." + self.version +".fasta"))
        return fasta
