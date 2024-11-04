from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import os
import sys

#open human.fa and make a list of seqrecord objects
#open mouse.fa and create an output file

human_seqs = list(SeqIO.parse("human.fa","fasta"))
mouse_seqs = "mouse.fa"
output_file = "homework-3-results.txt"

#create a blast database with the mouse data if it doesn't already exist

if not os.path.exists("mouse.nsq"):
   makenewdb = f"makeblastdb -in {mouse_seqs} -dbtype prot -input_type fasta -out mouse"
   os.system(makenewdb)

#write output to a text file

with open (output_file, "w") as output:
   
   for record in human_seqs:

      #write the seqrecord to its own temporary fasta file

      temp_file = "temp.fasta"
      
      with open (temp_file, "w") as handle:
         SeqIO.write(record, handle, "fasta")

      #define and run the BLAST program to output to an xml file

      blastp_cline = NcbiblastpCommandline(
         cmd="~/ncbi-blast-2.16.0+/bin/blastp",
         query=temp_file,
         db="mouse",
         evalue=0.001,
         matrix = "BLOSUM80",
         out="out.xml",
         outfmt=5)

      blastp_cline()

      #parse the xml file to find best alignment for the human sequence

      with open ("out.xml") as result_handle:
         blast_records = NCBIXML.parse(result_handle)

         best_alignment = None
         best_ID = None
         best_evalue = 0.001
         best_bitscore = 0

         for br in blast_records:
            
            for alignments in br.alignments:
               
               for hsp in alignments.hsps:
                  
                  if hsp.expect < best_evalue:
                     
                     best_alignment = hsp
                     best_ID = alignments.title
                     best_evalue = hsp.expect
                     best_bitscore = hsp.bits

         output.write("Human Sequence ID: " + str(record.id) + "\n")
         
         #if an alignment with e-value < 0.001 exists, write to output file
         
         if best_alignment:
            output.write("Mouse Sequence ID: " + best_ID + "\n" +
                         "Alignment: " + best_alignment.match + "\n" +
                         "E- value: " + str(best_evalue) + "\n" +
                         "Bitscore: " + str(best_bitscore) + "\n \n")
         else:
            output.write("No homologous sequence with e < 0.001" + "\n")

#get rid of temporary fasta file used before

os.remove(temp_file)

print("Done")
