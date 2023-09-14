import requests
import os

tritrypdb_species = {
    "tbrucei": "TbruceiTREU927",
    "lmexicana": "LmexicanaMHOMGT2001U1103"
}

for species in tritrypdb_species:
    tritrypdb_version = str(65)
    tritrypdb_url = "https://tritrypdb.org/common/downloads/release-" + tritrypdb_version + "/" + tritrypdb_species[species] + "/fasta/data/TriTrypDB-" + tritrypdb_version + "_" + tritrypdb_species[species] + "_AnnotatedProteins.fasta"
    if not os.path.isfile(species + ".fasta"):
        r = requests.get(tritrypdb_url)
        with open(species + ".fasta", "w") as fasta_file:
            fasta_file.write(r.text)

from afkost import _Fasta
tbrucei = _Fasta("tbrucei.fasta")
lmexicana = _Fasta("lmexicana.fasta")

from afkost import Sequence
esb1 = Sequence(tbrucei.sequences["Tb927.10.3800:mRNA-p1"])
esb1.plot(save_path="esb1.tb")
rbp6 = Sequence(tbrucei.sequences["Tb927.3.2930:mRNA-p1"])
rbp6.plot(save_path="rbp6.tb")
kharon1 = Sequence(tbrucei.sequences["Tb927.10.8940:mRNA-p1"])
kharon1.plot(save_path="kharon1.tb")
basalin = Sequence(tbrucei.sequences["Tb927.7.3130:mRNA-p1"])
basalin.plot(save_path="basalin.tb")
vex1 = Sequence(tbrucei.sequences["Tb927.11.16920:mRNA-p1"])
vex1.plot(save_path="vex1.tb")
marp = Sequence(tbrucei.sequences["Tb927.10.10360:mRNA-p1"])
marp.plot(save_path="marp.tb")
gapdh = Sequence(tbrucei.sequences["Tb927.6.4280:mRNA-p1"])
gapdh.plot(save_path="gapdh.tb")

from afkost import KmerMatrix
matrix = KmerMatrix()
matrix_tbrucei = matrix.composition_from_fasta("tbrucei.fasta")
gapdh.composition_dissimilarity_stats(matrix_tbrucei)
gapdh.kmer_outlier_stats(matrix_tbrucei)

esb1 = Sequence(lmexicana.sequences["LmxM.30.0660.1-p1"])
esb1.plot(save_path="esb1.lm")
rbp6 = Sequence(lmexicana.sequences["LmxM.08_29.2830.1-p1"])
rbp6.plot(save_path="rbp6.lm")
kharon1 = Sequence(lmexicana.sequences["LmxM.36.5850.1-p1"])
kharon1.plot(save_path="kharon1.lm")
basalin = Sequence(lmexicana.sequences["LmxM.22.1070.1-p1"])
basalin.plot(save_path="basalin.lm")
