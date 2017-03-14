#!/usr/bin/env python
from __future__ import print_function
import ftplib
import sys
import os.path

genome_list_fh = open("bct_genomes.tsv","a")

ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
ftp.connect()
ftp.login()
base_path="/genomes/refseq/bacteria"
ftp.cwd(base_path)
#Get species
species_list = ftp.nlst()
species_list = species_list[ species_list.index("Mycoplasma_moatsii"):]

#for each species
for species in species_list:
	print("Species: {}".format(species))
	species_path = base_path + "/{}".format(species)

	#Go into species assembly directory
	try:
		ftp.cwd(species_path)
	except Exception as e:
		print("{}\tCannot go into dir {}".format(species,species_path),file=sys.stderr)
		print("{}\t{}"r.,format(species,e),file=sys.stderr)
		continue
		
	assembly_folders = ftp.nlst()

	asm_type = None
	if "representative" in assembly_folders:
		asm_type = "representative"
	elif "reference" in assembly_folders:
		asm_type = "reference"
	elif "latest_assembly_versions" in assembly_folders:
		asm_type =  "latest_assembly_versions"
	elif "all_assembly_versions" in assembly_folders:
		asm_type =  "all_assembly_versions"
		
	
	if not asm_type:
		print("{}\tNo valid assemblies".format(species),file=sys.stderr)
		continue


	print("Entering folder {}".format(asm_type))
	try:
		asm_type_path = "{}/{}".format(species_path,asm_type)
		ftp.cwd(asm_type_path)
	except Exception as e:
		print("{}\tCannot enter dir {}".format(species,asm_type_path),file=sys.stderr)
		print("{}\t{}".format(species,e),file=sys.stderr)
		continue
		

	assemblies = ftp.nlst()
	if "suppressed" in asm_file_list:
		print("{}\tAssembly has been suppressed".format(species,asm_path),file=sys.stderr)
		break

    #Get assemblies
	for asm in assemblies:
		print("Entering assembly {}".format(asm))
		try:
			asm_path = asm_type_path + "/" + asm
			ftp.cwd(asm_path)
		except Exception as e:
			print("{}\tCannot enter dir {}".format(species,asm_path),file=sys.stderr)
			print("{}\t{}".format(species,e),file=sys.stderr)
			continue
		
		asm_file_list = ftp.nlst()
		#Get list of files to download
		files_to_dl = [f for f in asm_file_list if (f.endswith("genomic.fna.gz")) and ("cds_from_genomic" not in f) and ("rna_from_genomic" not in f)]
		if not files_to_dl:
			#Try downloading the cds
			files_to_dl = [f for f in asm_file_list if (f.endswith("cds_from_genomic.fna.gz")) ]

		if not files_to_dl:
			print("{}\tAssembly {} has no valid fasta".format(species,asm_path),file=sys.stderr)
			continue
			

		#Download genome fastas
		for filename in files_to_dl:
			#Verify the file has not been downloaded already
			if not os.path.isfile("fna/{}".format(filename)):
				print("Downloading file {}".format(filename))
				print('{}\t{}\t{}\t{}'.format(species,asm_type,asm,filename),file=genome_list_fh)
				try:
					with open("fna/{}".format(filename),"wb") as fh:
						ftp.retrbinary("RETR {}".format(filename), fh.write, 8*1024)
				except Exception as e:
					print("{}\tError downloading {}".format(filename),file=sys.stderr)
					print("{}\t{}".format(species,e),file=sys.stderr)
			else:
				print("{}\tFile {} already downloaded".format(species,filename),file=sys.stderr)

		#Download only the first assembly that is successfully downloaded
		# unless the assemblies are reference or representative
		if asm_type in ["latest_assembly_versions","all_assembly_versions"]:
			break
			
ftp.close()
genome_list_fh.close()
