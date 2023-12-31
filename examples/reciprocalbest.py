from afkost import blast
from afkost.databases import TriTrypDB

#blast = blast.Blast(search_tool="blast")
blast = blast.Blast()
tritrypdb = TriTrypDB()

ref = "TbruceiTREU927"
species = [
    "LmexicanaMHOMGT2001U1103",
    "TcruziDm28c2018"
]
for specie in species:
    tritrypdb.fetch_fasta(specie)

ref_seqs = tritrypdb.sequences(ref)
for specie in species:
    for gene in ref_seqs.sequences:
        result = blast.reciprocal_search(tritrypdb.fasta_path(ref), tritrypdb.fasta_path(specie), gene, ref_seqs.sequences[gene])
        print(result)
