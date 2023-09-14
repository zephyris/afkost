import os
import requests
import gzip
from functools import cached_property

class UniProt:
    def __init__(self, version: str = None):
        """
        Initialises the UniProt class, including cache directory for data storage.

        Named arguments:
        version -- use this version of uniprot (default latest)
        """
        self.cache_path = "_uniprot"
        if not os.path.isdir(self.cache_path):
            os.mkdir(self.cache_path)
        self._version = version

    @cached_property
    def version(self):
        """
        Returns the version, or fetches the latest version if `self._version` is `None`
        """
        if self._version is not None:
            return self._version
        else:
            self.fetch_index()
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "r") as index_file:
                lines = index_file.read().splitlines()
                lines = [x for x in lines if x]
                for i in range(len(lines)):
                    if len(lines[i]) > "Release":
                        if lines[i:len("Release")] == "Release":
                            return lines[i].split()[1][:-1]

    @cached_property
    def proteome_index(self):
        self.fetch_index()
        with open(os.path.join(self.cache_path, "_proteomes.txt"), "r") as index_file:
            lines = index_file.read().splitlines()
            lines = [x for x in lines if x]

            proteomes = {}
            i = 0
            #Find header line of table
            while lines[i] != "Proteome_ID\tTax_ID\tOSCODE\tSUPERREGNUM	#(1)\t#(2)\t#(3)\tSpecies Name":
                i += 1
            i += 1
            #Until next section
            while lines[i][0] != "-":
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

    def fetch_index(self):
        """
        Downloads and saves to disk the uniprot readme which includes an index of species
        """
        if self._version is None:
            url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"
        else:
            url = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-/" + self._version
        try:
            r = requests.get(url)
            # write to temporary path
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "w") as index_file:
                index_file.write(r.text)
                index_file.close()
            # rewrite to a new path with version (self.version may need to read `_proteomes.txt`)
            with open(os.path.join(self.cache_path, "_proteomes.txt"), "r") as index_file:
                with open(os.path.join(self. cache_path, "_proteomes." + self.version + ".txt")) as new_index_file:
                    new_index_file.write(index_file.read())
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)

    def fetch_fasta(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        if not os.path.isfile(species + ".gzip"):
            # construct url
            if self._version is None:
                url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s/%s_%s.fasta.gz" % (self.proteome_index[species]["supergenum"].capitalize(), species, species, self.proteome_index[species]["tax_id"])
            else:
                url = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-" + self._version + "/knowledgebase/reference_proteomes/%s/%s/%s_%s.fasta.gz" % (self.proteome_index[species]["supergenum"].capitalize(), species, species, self.proteome_index[species]["tax_id"])
            # download and save
            try:
                # download as binary file, gzip compressed
                r = requests.get(url, stream=True)
                with open(self.gzip_path(species), "wb") as gzip_file:
                    for chunk in r.iter_content(chunk_size=1024): 
                        if chunk:
                            gzip_file.write(chunk)
                # decompress to plain text
                with gzip.open(self.gzip_path(species), mode="rt") as gzip_file:
                    with open(self.fasta_path(species), mode="w") as fasta_file:
                        fasta_file.write(gzip_file.read())
            except requests.exceptions.RequestException as e:
                raise SystemExit(e)

    def gzip_path(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        return os.path.join(self.cache_path, species + "." + self.version + ".gzip")

    def fasta_path(self, species: str):
        """
        Downloads and saves to disk a protein sequence fasta file for the specified species.

        Required arguments:
        species: Species/strain ID, as listed in uniprot
        """
        return os.path.join(self.cache_path, species + "." + self.version + ".fasta")

    def sequences(self, species: str):
        """
        Returns sequences for a uniprot species as a _Fasta instance

        Required arguments:
        species: species/strain name, as found on uniprot

        Returns:
        `_Fasta` instance containing the sequences for that species
        """
        from afkost import _Fasta
        self.fetch_fasta(species)
        fasta = _Fasta(os.path.join(self.cache_path, species + "." + self.version +".fasta"))
        return fasta
