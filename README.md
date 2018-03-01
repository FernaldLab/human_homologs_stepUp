# human_homologs_stepUp
Below are basic steps for procedure we used to get human homologs for A.burtoni genes by stepping up through other better-annotated fish species.
Most text is copied from email to ZB who was repeating the process starting with T.guttata (zebra finch)

1) Get fasta files
- z.f. transcriptome (ftp://ftp.ncbi.nlm.nih.gov/genomes/Taeniopygia_guttata/RNA/rna.fa.gz)
- z.f. proteome (ftp://ftp.ncbi.nlm.nih.gov/genomes/Taeniopygia_guttata/protein/protein.fa.gz)
- proteomes for other bird (or whatever) species, we used the proteomes in Ensembl for 5 other fish species

2) Prepare blast databases with Makeblastdb
- protein database with all other bird species, we used cat to combine all the proteome fasta files then ran Makeblastdb once, I'll call this "otherProteinDb"
- protein database for z.f., I'll call this "zfProteinDb"

3) Use blastx to map the z.f. transcriptome to otherProteinDb
For us this took ~1 day per species. Here's what we ran:
blastx  -query H_burtoni_rna.fa \\\
	-db ./db/Drer_Olat_Onil_Trub_ENS_pep \\\
	-out H_burtoni_rna_blastx_FISH_ENS_top1 \\\
	-outfmt '7' \\\
	-max_target_seqs 1

-Note: -outfmt and -max_target_seqs must be set to '7' and '1', respectively

4)  Parse the blastx results
Use 1_parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py, the arguments are described in the script. There will be 5 output files, one with a name ending in "topByTranscript" that I'll call "originalBlastx.topByTranscript". The columns of this file are described in the script, but basically they'll be z.f. gene id, z.f. transcript id, the protein id of the best hit for that transcript from otherProteinDb, full z.f. gene description, and the blastx bit score. 

5) Create a fasta of protein sequences from the hits in originalBlastx.topByTranscript
I don't have one script to do this, but if you make a text file of just the protein ids in column 3 of originalBlastx.topByTranscript you can use 1b_ENSid_toFasta.py to create a protein fasta file. If otherProteinDb contained more than one species you'll need to do this separately for each species then cat the fasta files together into a file I'll call "proteinBestHit.fa".

6) Use blastp to map proteinBestHit.fa to zfProteinDb
If hits from originalBlastx.topByTranscript map back to the proteins made by the original z.f. mRNAs then they're likely homologs. 

7) Parse the blastp results
Use parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py again, one of the outputs will be "reciprocalBlastp.topByTranscript". This file and originalBlastx.topByTranscript will be loaded into R for the next step. You also need a lookup file that maps z.f. gene ids/symbols to transcript and protein ids. I can send code for parsing a gff file into the lookup, but I think I remember you already made one?

8) Load blastx and blastp results into R and try to get human homologs 
Right now everything is still at transcript/protein level but you want to get to gene level to get homologs. Follow what I did in burtoni_annotations_blastx_and_reciprocal_Feb2016.R to combine the blastx and blastp results and make decisions about whether genes had enough transcripts pass reciprocal blasting to confidently call them homologs in the other birds. Then use biomaRt to get human homolog entrez ids for GO analysis.
