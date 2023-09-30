import os
import requests
import gzip
from functools import cached_property

class UniProt:
    def __init__(self):
        """
        Initialises the UniProt class, including cache directory for data storage.
        """
        self.cache_path = "_uniprot_cache"
        if not os.path.isdir(self.cache_path):
            os.mkdir(self.cache_path)
        # version is always `None`, as random access to older versions is not possible
        self._version = None

    @cached_property
    def version(self):
        """
        Returns the version, or fetches the latest version if `self._version` is `None`
        """
        if self._version is not None:
            return self._version
        elif not os.path.isfile(os.path.join(self.cache_path, "_proteomes.txt")):
            # if _proteomes.txt does not exist then no interim index has been retrieved
            self.fetch_index()
            return self.version
        else:
            # interpret version from _proteomes.txt and return
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "r") as index_file:
                lines = index_file.read().splitlines()
                lines = [x for x in lines if x]
                for i in range(len(lines)):
                    if len(lines[i]) > len("Release"):
                        if lines[i][0:len("Release")] == "Release":
                            return lines[i].split()[1][:-1]

    def fetch_index(self):
        """
        Downloads and saves to disk the uniprot readme which includes an index of species
        """
        if self._version is None:
            url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"
        else:
            url = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-" + self._version + "/knowledgebase/reference_proteomes/README"
        try:
            r = requests.get(url)
            # write to temporary path
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "w") as index_file:
                index_file.write(r.text)
                index_file.close()
            # rewrite to a new path with version (after writing to a tempory path, as self.version may need to read `_proteomes.txt`)
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "r") as index_file:
                with open(os.path.join(self. cache_path, "_proteomes." + self.version + ".txt"), "w") as new_index_file:
                    new_index_file.write(index_file.read())
            # remove `_proteomes.txt`
            os.remove(os.path.join(self.cache_path, "_proteomes.txt"))
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)

    @cached_property
    def proteome_index(self):
        self.fetch_index()
        with open(os.path.join(self.cache_path, "_proteomes." + self.version + ".txt"), "r") as index_file:
            lines = index_file.read().splitlines()
            lines = [x for x in lines if x]

            proteomes = {}
            i = 0
            #Find header line of table
            while lines[i] != "Proteome_ID\tTax_ID\tOSCODE\tSUPERREGNUM	#(1)\t#(2)\t#(3)\tSpecies Name":
                i += 1
            i += 1
            #Until next section
            while lines[i][0:2] == "UP" and i < len(lines):
                line = lines[i].split("\t")
                result = {
                    "proteome_id": line[0],
                    "tax_id": line[1],
                    "oscode": line[2],
                    "supergenum": line[3],
                    "main_entries": int(line[4]),
                    "main_entries": int(line[5]),
                    "gene2acc_entries": int(line[6]),
                    "name": line[7]
                }
                proteomes[result["proteome_id"]] = result
                i += 1
            # return index
            return proteomes

    def fetch_fasta(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        if not os.path.join(self.cache_path, species + "." + self.version + ".gzip"):
            # construct url
            if self._version is None:
                url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s/%s_%s.fasta.gz" % (self.proteome_index[species]["supergenum"].capitalize(), species, species, self.proteome_index[species]["tax_id"])
            else:
                url = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-" + self._version + "/knowledgebase/reference_proteomes/%s/%s/%s_%s.fasta.gz" % (self.proteome_index[species]["supergenum"].capitalize(), species, species, self.proteome_index[species]["tax_id"])
            # download and save
            try:
                # download as binary file, gzip compressed
                r = requests.get(url, stream=True)
                with open(os.path.join(self.cache_path, species + "." + self.version + ".gzip"), "wb") as gzip_file:
                    for chunk in r.iter_content(chunk_size=1024): 
                        if chunk:
                            gzip_file.write(chunk)
                # decompress to plain text
                with gzip.open(os.path.join(self.cache_path, species + "." + self.version + ".gzip"), mode="rt") as gzip_file:
                    with open(os.path.join(self.cache_path, species + "." + self.version + ".fasta"), mode="w") as fasta_file:
                        fasta_file.write(gzip_file.read())
            except requests.exceptions.RequestException as e:
                raise SystemExit(e)

    def gzip_path(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        self.fetch_fasta(species)
        return os.path.join(self.cache_path, species + "." + self.version + ".gzip")

    def fasta_path(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        self.fetch_fasta(species)
        return os.path.join(self.cache_path, species + "." + self.version + ".fasta")

    def sequences(self, species: str):
        """
        Returns sequences for a uniprot species as a Fasta instance

        Required arguments:
        species: species/strain name, as found on uniprot

        Returns:
        `Fasta` instance containing the sequences for that species
        """
        from afkost import Fasta
        self.fetch_fasta(species)
        fasta = Fasta(os.path.join(self.cache_path, species + "." + self.version +".fasta"))
        return fasta

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
