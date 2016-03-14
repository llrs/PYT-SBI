#This module provides code to work with the WWW version of BLAST provided by the NCBI.
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

def run_BLAST(query):
	blast_record_out = open('blast_record.fa','w')
	record = SeqIO.read(query, format="fasta")
	result_handle=NCBIWWW.qblast('blastp','nr',record.format("fasta"), 
								perc_ident=30, 
								service='psi',  
								hitlist_size=1000) 
	blast_record = NCBIXML.read(result_handle)
	
	blast_record_out.write(">%s\n%s\n" %(record.id, record.seq))
	genere_set=set()
	i=0
	for alignment in blast_record.alignments:
  		for hsp in alignment.hsps:	
  			percentage_identity = 100*float(hsp.identities)/float(len(record.seq))
  			if percentage_identity > 30:
  				if '[' in alignment.title:
					specie=((alignment.title).split('[')[len((alignment.title).split('['))-1]).split(']')[0]
					genere=(specie.split(' '))[0]
					if not set([genere]).issubset(genere_set) and genere.find('Cloning')<0 and genere.find('synthetic')<0:
						blast_record_out.write(">%d-%s_%.2f\n%s\n" %(i, genere, 
																 percentage_identity, 
																 hsp.sbjct))
					genere_set.add(str(genere))
					i+=1 
									

if __name__ == "__main__":
	run_BLAST('1cd8.fasta')