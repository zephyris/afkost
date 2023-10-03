import subprocess
import os
from afkost import Fasta

class Blast:
    def __init__(self, search_tool: str = "diamond", sequence_type: str = "protein"):
        """
        Initialise a sequence search instance

        Named arguments:
        search_tool -- the search tool to use, either `"blast"` or `"diamond"`
        sequence_type -- the sequence type, either `"nucleotide"` or `"protein"`
        """
        # check search_tool value is valid
        search_tools = ["diamond", "blast"]
        self.search_tool = search_tool
        if self.search_tool not in search_tools:
            raise ValueError("`search_tool` must be " + " or ".join(search_tools))
        # check sequence_type value is valid
        sequence_types = ["protein", "nucleotide"]
        self.sequence_type = sequence_type
        if self.sequence_type not in sequence_types:
            raise ValueError("`search_tool` must be " + " or ".join(sequence_types))
        self.minimum_evalue = "0.00001"

    def _check_install(self, verbose: bool = False):
        """
        Check if the necessary programs are installed, and thow an error if not.
        ncbi-blast+ is required for `"blast"` search_type
        diamond is required for `"diamond"` search_type

        Named arguments:
        verbose -- print a verbose output, boolean, default `False`
        """
        command = {
            "diamond": ["diamond", "--help"], 
            "blast": ["makeblastdb", "-help"]
        }
        proc = subprocess.Popen(command[self.search_tool], stderr = subprocess.PIPE, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
        result = proc.communicate()[0].decode("utf-8").splitlines()
        if proc.returncode != 0:
            print(result)
            assert AssertionError(command[self.search_tool] + " did not run correctly, is it installed and available in the system path?")
        if verbose:
            print("'" + " ".join(command[self.search_tool]) + "' ran successfully, " + self.search_tool + " appears to be installed correctly")

    def search(self, query_sequence: str, fasta_path: str):
        """
        Do a blast or diamond sequence search.

        Required arguments:
        query_sequence -- query as a string
        fasta_path -- fasta file to use as a database to search against

        Returns a list of hits, in the form [subject sequence name, evalue, subject sequence hit], and the additional list entry [full subject sequence hit] if using diamond
        """
        # TODO: Error handling from failure to run blasts
        # check install
        self._check_install()
        # setup query name, currently unused
        query_name = None
        if query_name is None:
            query_name = "none"
        # select search program name
        search_program = {
            "protein": "blastp",
            "nucleotide": "blastn"
        }
        if self.search_tool == "diamond":
            # make the database
            if not os.path.isfile(fasta_path + ".dmnd"):
                command = ["diamond", "makedb", fasta_path, "--quiet"]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
                #if proc.returncode != 0:
                #    print(results)
                #    assert AssertionError("diamond database generation failed")
            # do the search
            search_command = ["-d", fasta_path, "--quiet", "-e", self.minimum_evalue, "-q", "-", "--outfmt", "6", "sseqid", "evalue", "sseq", "full_sseq"]
            command = ["diamond", search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            proc.stdin.write((">%s\n%s" % (query_name, query_sequence)).encode('utf-8'))
            results = proc.communicate()[0].decode("utf-8")
            #if proc.returncode != 0:
            #    print(results)
            #    assert AssertionError("diamond search failed")
        if self.search_tool == "blast":
            # make the database
            if not os.path.isfile(fasta_path + ".phr"):
                dbtype = {
                    "protein": "prot",
                    "nucleotide": "nucl"
                }
                command = ["makeblastdb", "-dbtype", dbtype[self.sequence_type], "-in", fasta_path]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
                #if proc.returncode != 0:
                #    print(results)
                #    assert AssertionError("blast database generation failed")
            # do the search
            search_command = ["-db", fasta_path, "-evalue", self.minimum_evalue, "-query", "-", "-outfmt", "6 sseqid evalue sseq"]
            command = [search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            proc.stdin.write((">%s\n%s" % (query_name, query_sequence)).encode('utf-8'))
            results = proc.communicate()[0].decode("utf-8")
            #if proc.returncode != 0:
            #    print(results)
            #    assert AssertionError("diamond search failed")
        results = results.splitlines()
        results = [x.split("\t") for x in results]
        return results

    def reciprocal_search(self, origin_fasta_path: str, subject_fasta_path: str, query_name: str = None, query_sequence: str = None, verbose: bool = False):
        """
        Do a reciprocal best blast or diamond sequence search. Either `query_name` or `query_sequence` must be met.
        Validity of `query_name` and `query_sequence` are NOT checked, ie. it is trusted that, if set, they exactly match an entry in the fasta file at `origin_fasta_path`.
        Providing `query_name` and `query_sequence` is fasta, as it prevents additional lookup from the fasta file at `origin_fasta_path`.
        Using diamond is faster than blast as it returns full subject sequence which can be used in the reciprocal search.

        Required arguments:
        origin_fasta_path -- fasta file to use as the originating database
        subject_fasta_path -- fasta file to use as the target database

        Named arguments:
        query_sequence -- query as a string
        query_name -- query name
        verbose -- print additional information for when a reciprocal blast was not found
    
        Returns None if there is no reciprocal best hit, list in the form [query_name, evalue_forward, subject_name, evalue_reverse] if there is.
        """
        # check query
        if query_name is None and query_sequence is None:
            raise ValueError("Either `query_name` or `query_sequence` must be given as an argument")
        # if necessary, look up query sequence
        if query_sequence is None:
            query_fasta = Fasta(origin_fasta_path)
            if query_name not in query_fasta.sequences:
                raise ValueError(query_name + " name not found in the fasta file at " + origin_fasta_path)
            query_sequence = query_fasta.sequences[query_name]
        # if necessary, look up query name
        if query_name is None:
            query_fasta = Fasta(origin_fasta_path)
            for name in query_fasta.sequences:
                if query_fasta.sequences[name] == query_sequence:
                    query_name = name
            if query_name is None:
                raise ValueError(query_sequence + "sequence not found in the fasta file at " + origin_fasta_path)
        # do the forward search, if no hits then return None
        forward_result = self.search(query_sequence, subject_fasta_path)
        if len(forward_result) == 0:
            if verbose: print("no forward hits")
            return None
        # do the reverse search, if no hits then return None
        forward_result = forward_result[0]
        if len(forward_result) < 4:
            # if index 3 not in result (blast only) then look up full subject sequence ID
            subject_fasta = Fasta(subject_fasta_path)
            forward_result.append(subject_fasta.sequences[forward_result[0]])
        reverse_result = self.search(forward_result[3], origin_fasta_path)
        if len(reverse_result) == 0:
            if verbose: print("no reverse hits")
            return None
        # check for reciprocality
        reverse_result = reverse_result[0]
        if reverse_result[0] != query_name:
            if verbose: print("no reciprocal match")
            return None
        else:
            if verbose: print("reciprocal match")
            return [query_name, forward_result[1], forward_result[0], forward_result[1]]

    def reciprocal_genomewide(self, fasta1_path: str, fasta2_path: str):
        """
        List all reciprocal best blast or diamond sequence search, pairwise for two entire fasta files.

        Required arguments:
        fasta1_path -- fasta file to use as the originating database
        fasta2_path -- fasta file to use as the target database

        Returns a list of reciprocal best blast hits, each entry in the form [query_name, forward_evalue, forward_pident, subject_name, reverse_evalue, reverse_pident]
        """
        # TODO: Error handling from failure to run blasts
        # check install
        self._check_install()
        # setup query name, currently unused
        query_name = None
        if query_name is None:
            query_name = "none"
        # select search program name
        search_program = {
            "protein": "blastp",
            "nucleotide": "blastn"
        }
        if self.search_tool == "diamond":
            # make the databases
            if not os.path.isfile(fasta1_path + ".dmnd"):
                command = ["diamond", "makedb", fasta1_path, "--quiet"]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
            if not os.path.isfile(fasta2_path + ".dmnd"):
                command = ["diamond", "makedb", fasta2_path, "--quiet"]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
            # do the searches
            # forward
            search_command = ["-d", fasta1_path, "--quiet", "-e", self.minimum_evalue, "-q", fasta2_path, "--outfmt", "6", "qseqid", "sseqid", "evalue", "pident"]
            command = ["diamond", search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            results1 = proc.communicate()[0].decode("utf-8")
            # reverse
            search_command = ["-d", fasta2_path, "--quiet", "-e", self.minimum_evalue, "-q", fasta1_path, "--outfmt", "6", "qseqid", "sseqid", "evalue", "pident"]
            command = ["diamond", search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            results2 = proc.communicate()[0].decode("utf-8")
        if self.search_tool == "blast":
            # make the databases
            dbtype = {
                "protein": "prot",
                "nucleotide": "nucl"
            }
            if not os.path.isfile(fasta1_path + ".phr"):
                command = ["makeblastdb", "-dbtype", dbtype[self.sequence_type], "-in", fasta1_path]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
            if not os.path.isfile(fasta2_path + ".phr"):
                command = ["makeblastdb", "-dbtype", dbtype[self.sequence_type], "-in", fasta2_path]
                proc = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
                results = proc.communicate()[0].decode("utf-8")
            # do the search
            # forward
            search_command = ["-db", fasta1_path, "-evalue", self.minimum_evalue, "-query", fasta2_path, "-outfmt", "6 qseqid sseqid evalue pident"]
            command = [search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            results1 = proc.communicate()[0].decode("utf-8")
            # reverse
            search_command = ["-db", fasta2_path, "-evalue", self.minimum_evalue, "-query", fasta1_path, "-outfmt", "6 qseqid sseqid evalue pident"]
            command = [search_program[self.sequence_type]] + search_command
            proc = subprocess.Popen(command, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
            results1 = proc.communicate()[0].decode("utf-8")
        results1 = [x.split("\t") for x in results1.splitlines()]
        results2 = [x.split("\t") for x in results2.splitlines()]
        # make dict of result1 (forward search) best hits, indexed by qseqid
        results1best = {}
        for i in range(len(results1)):
            # hits are in e value order, so earliest line is best hit
            if results1[i][0] not in results1best:
                results1best[results1[i][0]] = results1[i]
        # find reciprocal best hits
        reciprocalbest = []
        results2best = {}
        for i in range(len(results2)):
            # check for and record best hit
            if results2[i][0] not in results2best:
                results2best[results2[i][0]] = results2[i]
                # if it's a best hit, check if reverse search gave the same pair
                # if results2 (reverse search) sseqid is in forward search index
                if results2[i][1] in results1best:
                    # if reverse search qseqid is the same as the best hit from the forward search
                    if results1best[results2[i][1]][1] == results2[i][0]:
                        reciprocalbest.append([results2[i][0], results1best[results2[i][1]][2], results1best[results2[i][1]][3], results2[i][1], results2[i][2], results2[i][3]])
        return reciprocalbest
