# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

# Load the necessary modules
import csv, json, subprocess
from Bio import Entrez, SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO

def ncbi2dict(file):
    
    dictionary = {}
    with open(file, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter = ",")
        for row in csvreader:
            organism = ("_".join(row["Organism Name"].split(" "))).replace(".","")
            replicons, replicons_fv = row["Replicons"].split(";"), []
            for replicon in replicons:
                replicons_fv.append(replicon.split(":")[-1].split(".")[0])
            dictionary[" ".join(replicons_fv)] = organism
        return dictionary
    
def gb2faa(input_file, output_file):
    
    '''
    Download fasta_cds_aa using a list of tsv-formated GenBank accession numbers
    '''
    db_ids = input_file.split(" ")
    # save fasta_cds_aa of each GenBank into a single FASTA file
    out_handle = open(output_file, "w")
    for nucleotide_id in db_ids:
        try:
            print ("Downloading "+nucleotide_id+"...")
            handle = Entrez.efetch(db="nucleotide", id = nucleotide_id, rettype="fasta_cds_aa")
            out_handle.write(handle.read())
            handle.close()
        except:
            pass
    out_handle.close()


def blastdb(input_file, output_file):
    
    '''
    Construct the blast_db using the bash script makeblastdb.sh
    '''
    
    bashCommand = "bash makeblastdb.sh "+input_file+" "+output_file
    subprocess.call(bashCommand.split())
 
    
def blast_search(query, database, cutoff, nhits, coverage):
    
    """ Local BLASTp search to detect orthologues.
        Receives a query protein accession, an e-value cut off a coverage cut off
        and the maximum number of hits to be retrieved.

        Returns a list containing the protein accessions for the BLASTP hits.
    """
    
    #blastp against local database
    blastp_cline = NcbiblastpCommandline(query=query, db=database,\
                                        evalue=cutoff, outfmt='"6 std qcovs"',\
                                        num_alignments = nhits)
    stdout, sterr = blastp_cline()
    
    blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",\
                                   fields=['std','qcovs'])
        
    #return the hits passing retrieved thresholds
    validated_TF = []
    for blast_record in blast_records:
        for alignment in blast_record:
            if alignment.query_coverage >= coverage:
                acc = str(alignment.id_all[0])
                validated_TF.append(acc)
    
    validated_TF = list(set(validated_TF))
    return validated_TF


############################################################################

# Load the configuration file
with open("test_input/1_test_input_rBLAST.json") as json_conf : 
    conf = json.load(json_conf)

# tell to who I am
Entrez.email = conf["email"]
Entrez.api_key = conf["api_key"]

# copy the makeblastdb.sh script to the current dir
bashCommand = "cp "+conf["makeblastdb_file"]+" ."
subprocess.call(bashCommand.split())

PBP5_db = conf["blastp_db"]
PBP5_queries = {}
for record in SeqIO.parse(conf["query_seqs"], "fasta"):
    PBP5_queries[record.id] = str(record.seq)

dictionary = ncbi2dict(conf["query_genomes"])

out_handle_F = open(conf["output_file"], "w")
for key in dictionary.keys():
    # get faa file
    gb2faa(key, "TEMP_ncbi.fasta")
    # & construct the BLAST_DB
    blastdb("TEMP_ncbi.fasta", "TEMP_ncbi")
    
    subject_db = {}
    for record in SeqIO.parse("TEMP_ncbi.fasta", "fasta"):
        subject_db[record.id] = str(record.seq)
    
    # rBLAST
    for seq in PBP5_queries:
        out_handle_seq = open("query_test.faa", "w")
        out_handle_seq.write(">"+seq+"\n")
        out_handle_seq.write(PBP5_queries[seq]+"\n")
        out_handle_seq.close()
        
        blast1_list = blast_search("query_test.faa", "TEMP_ncbi", 1e-20,100,75)
        
        blast2_list = []
        try:
            for blast1 in blast1_list:
                
                out_handle_seq = open("query_test2.faa", "w")
                out_handle_seq.write(">"+blast1+"\n")
                out_handle_seq.write(subject_db[blast1]+"\n")
                out_handle_seq.close()
                
                blast2 = blast_search("query_test2.faa", PBP5_db, 1e-20,1,75)

                if blast2[0] == seq:
                    blast2_list.append(blast1)
        
            blast2_list = list(set(blast2_list))
            print (dictionary[key]+"\t"+key+"\t"+seq+"\t"+";".join(blast2_list))
            out_handle_F.write(dictionary[key]+"\t"+key+"\t"+seq+"\t"+";".join(blast2_list)+"\n")
            
        except:
            pass            
out_handle_F.close()
