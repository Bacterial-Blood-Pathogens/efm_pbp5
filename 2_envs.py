# -*- coding: utf-8 -*-

### miquelsanchezosuna
### Implemented for Python3

# Load the necessary modules
import csv
import os
import subprocess
import time
from Bio import  Entrez, SeqIO


def split_accession(acc):
    
    """ Splits a genbank/refseq record into prefix (e.g. NZ_MDWE or MDWE) and
        numeric component (e.g. 0100092), and returns them as a list.
    """

    #go through the record until we get the first numeric digit
    #store the prefix and use its length to get the numeric region
    prefix = ''
    for c in acc:
        if c.isdigit():
            break
        else:
            prefix=prefix + c

    suffix = acc[len(prefix):len(acc)]

    return [prefix, suffix]


def genome_record_retrieval(ortholog_acc, nucleotide, sleepy):
    
    """ Takes a protein accession as an input. Retrieves its IPG record.

        The idea here is to obtain, prioritarily, data from complete genome 
        records if they exist, from RefSeq (AC_ and NC_ accessions) . If no 
        RefSeq is available, then select complete genome records from GenBank 
        (AE, CP, CY accessions). Otherwise, select contigs or WGS scaffolds from 
        RefSeq (NT_, NW_, NZ_). If that fails, get contigs or WGS scaffolds from 
        GenBank (AAAA-AZZZ). Only when nothing else is available, select direct 
        GenBank submissions (U, AF, AY, DQ).

        Prioritizes each type of accession and returns the best record.

        Priority indices range from 7 (best, for a complete RefSeq record) and 6
        (complete GenBank record), to 5 (for complete RefSeq WGS) and all the 
        way to 3 (undetermined GenBank records).
        
        It returns a composite record with the nucleotide accession number for
        the "best" coding region, and the position and orientation of the CDS 
        within that accession, as well as the prioritization score obtained.
    """

    #Download IPG record for the specific ortholog
    records = None
    while records == None:
        try:
            records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
                                        rettype='ipg', retmode='xml'))
        except:
            print ("IPG record retrieval error: "+ortholog_acc)
            time.sleep(sleepy)
            pass
    
    #from the IPG record, retrieve all the genome accessions from all CDS
    #keeping only accession, location of start and strand, as well as 
    #priority score
    if 'ProteinList' in records["IPGReport"].keys():
        for idprotein in records["IPGReport"]["ProteinList"]:
            if 'CDSList' in idprotein.keys():
                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver'].split(".")[0]
                    if cds_acc == nucleotide:
                        cds_start = cds.attributes['start']
                        cds_stop = cds.attributes['stop']
                        cds_strand = cds.attributes['strand']
                        cds_org = cds.attributes['org'].replace(" ","_")

                        #create and append record
                        cds_rec = {'acc':cds_acc, 'start':cds_start, \
                                   'stop':cds_stop, 'strand':cds_strand,\
                                   'org':cds_org}
                        
                        return (cds_rec)

            else:
                return (None)
    #GenBank record that has no proper IPG record (yes, they exist;
    #see for instance: https://www.ncbi.nlm.nih.gov/protein/RJR51       119.1)
    #in these cases, there should be a CDS within the protein record that 
    #contains the information we want; priority should be lowest
    else:
#        TO BE IMPLEMENTED
#        records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
#                                            rettype='genbank', retmode='xml'))
        return (None)


def env_retrieval(nucid,start,stop,start_adj=5000,stop_adj=5000):
    
    """Download genetic environtment sequences from NCBI RefSeq/GenBank 
       database using the RefSeq/GenBank ID and the coordinates  of 
       gene of interest.
       The start_adj and stop_adj delimit the size.
       Returns a list of protein acc
    """
    
    s_start=int(start)-start_adj
    s_stop=int(stop)+stop_adj
    
    #Fetch and read the annotated GenBank record
    handle = None
    while handle == None:
        try: 
            handle = Entrez.efetch(db="nuccore",id=nucid, seq_start=s_start, seq_stop=s_stop, rettype='gbwithparts', retmode = "text")
        except:
            time.sleep(5)
            print ("Env retrieval error: "+nucid)
            pass
    
    # Find all coding regions in the returned GenBank sequence. 
    cds_aa = []
    records = SeqIO.parse(handle, "gb")
    for record in records:
        for f in record.features:
            if f.type == "CDS":
                try:
                    cds_aa.append(f.qualifiers["protein_id"][0]+"__"+f.qualifiers["product"][0])
                except:
                    pass
    
    return cds_aa


def protein_retrieval(protein_list, output_file):
    
    """ Download protein sequences from NCBI RefSeq/GenBank database and 
        save them into a FASTA file.
        The headers are modified to obtain the following structure:
        >ProteinIdentifier__Species
    """
    
    out_handle = open(output_file, "w")
    handle = Entrez.efetch(db="protein", id=protein_list, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    for record in records:
        out_handle.write(">"+record.id+"__"+record.description.split("[")[1].split("]")[0].replace(" ","_")+"\n")
        out_handle.write(str(record.seq)+"\n")
    out_handle.close()
    handle.close()


def COG_annotation(protein_file, hmm_db, evalue):
    
    """ Annotate a given protein file with single FASTA seq & return the COG id
    """
    
    bashCommand = "hmmscan --tblout cog_temp_arg.txt -E " +str(evalue)+" "+hmm_db+" "+protein_file
    subprocess.call(bashCommand.split())
    cog_file = "cog_temp_arg.txt"
    COG_acc = ""
    with open(cog_file) as fp2:
        line = fp2.readline()
        cnt = 1
        while line:
            line = fp2.readline()
            cnt += 1
            if not line.startswith("#") and line != "":
                COG_acc = line.split(" ")[0]
                break
            
    os.remove("cog_temp_arg.txt")
            
    return COG_acc


#######################################################

#Tell to who I am
Entrez.email = conf["email"]
Entrez.api_key = conf["api_key"]

#Open a file containing deired protein ids and save them into a list
protein = []
with open(conf["input_file"], 'r', encoding='ISO-8859-1') as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter = "\t")
    for row in csvreader:
        hits = row["PBP5_Hits"].split(";")
        protein += hits
        protein = list(set(protein))
protein.remove("")

#Create a dictionary using protein ids as keys and their location, env, etc as values
ortholog_summary_dict = {}
for acc in protein:
    if acc.count(".") >= 2:
        print ("Getting "+acc+" coordinates...")
        protein, nucc = acc.split("_prot_")[-1].split(".")[0], acc.split("_prot_")[0].split("|")[-1].split(".")[0]
        genome_record_retrieval_result = genome_record_retrieval(protein, nucc, 2)
        if genome_record_retrieval_result != None:
            print (genome_record_retrieval_result)
            ortholog_summary_dict[acc] = genome_record_retrieval_result
            ortholog_summary_dict[acc]["env"] = env_retrieval(ortholog_summary_dict[acc]["acc"], ortholog_summary_dict[acc]["start"], ortholog_summary_dict[acc]["stop"])
        
#COG annotation & printing
out_handle_F = open(conf["output_dir"], "w")
out_handle_F.write("Organism_name\tBlastp_hit\tNucleotide_Id\tEnv\n")
for protein_key in ortholog_summary_dict.keys():
    print ("Working with "+protein_key+" env...")
    env_COG = []
    for env_protein in ortholog_summary_dict[protein_key]["env"]:
        protein_retrieval([env_protein.split("__")[0]],"temp_prot.fasta")
        COG_protein = COG_annotation("temp_prot.fasta","/home/miquel/HMMERdbs/COG_database.hmm", 1e-05)
        env_COG.append(env_protein+"__"+COG_protein)
    if "+" in ortholog_summary_dict[protein_key]["strand"]:
        out_handle_F.write(ortholog_summary_dict[protein_key]["org"]+"\t"+protein_key+"\t"+ortholog_summary_dict[protein_key]["acc"]+"\t"+"\t".join(env_COG)+"\n")
    else:
        out_handle_F.write(ortholog_summary_dict[protein_key]["org"]+"\t"+protein_key+"\t"+ortholog_summary_dict[protein_key]["acc"]+"\t"+"\t".join(env_COG[::-1])+"\n")
out_handle_F.close()
