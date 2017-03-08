#!/usr/bin/env python
from __future__ import print_function
import ftplib
import sys

genome_list_fh = open("downloaded_genomes.txt","w")

ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
ftp.connect()
ftp.login()
base_path="/genomes/refseq/fungi"
ftp.cwd(base_path)
#Get species
species_list = ftp.nlst()
#for each species
for species in species_list:
	print("Species: {}".format(species))
	species_path = base_path + "/{}/latest_assembly_versions".format(species)
	#Go into species assembly directory
	try:
		ftp.cwd(species_path)
	except Exception as e:
		print("+ Cannot go into dir {}".format(species_path),file=sys.stderr)
		print(e,file=sys.stderr)
		continue
		
	#Get assemblies
	assemblies = ftp.nlst()
	for asm in assemblies:
		print("Entering assembly {}".format(asm))
		try:
			ftp.cwd(species_path + "/" +  asm)
		except Exception as e:
			print("+ Cannot enter dir {}".format(species_path+"/"+asm),file=sys.stderr)
			print(e,file=sys.stderr)
			continue
		
		asm_file_list = ftp.nlst()
		files_to_dl = [f for f in asm_file_list if (f.endswith("genomic.fna.gz")) and ("cds_from_genomic" not in f) and ("rna_from_genomic" not in f)]
		if not files_to_dl:
			#Try downloading the cds
			files_to_dl = [f for f in asm_file_list if (f.endswith("cds_from_genomic.fna.gz")) ]
		for filename in files_to_dl:
			print("Downloading file {}".format(filename))
			print('{}\t{}\t{}'.format(species,asm,filename),file=genome_list_fh)
			try:
				with open(filename,"wb") as fh:
					ftp.retrbinary("RETR {}".format(filename), fh.write, 8*1024)
			except Exception as e:
				print(">Error downloading {}".format(asm_filename),file=sys.stderr)
				print(e,file=sys.stderr)
			
ftp.close()
