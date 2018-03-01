# clear workspace
rm(list=ls());
# avoid character-factor automatic conversion
options(stringsAsFactors=F);
# move to directory that holds gff and transcriptome, need the trailing backslash 
setwd('/Users/abseq/Documents/_annotationsDec2015/');

source('/Volumes/fishstudies-1/_code/gff_cleaning_parsing_makeLookup_functions.R');

####################################################################################################

# process raw gff3 file from ncbi
# need comment.char arg to avoid reading comment lines
gff_rawFile = 'ref_AstBur1.0_scaffolds.gff3';

# set colClasses explicitly
#  if start/end are numeric instead of integer R may convert locations like 100000 to 1e5

colClasses = c(rep('character',3),                               # scaffold, source, feature type
               rep('integer',2),                                 # start, end
               rep('character',4));                              # score, strand, phase, metadata

# read file
gff_raw = read.table(gff_rawFile, colClasses=colClasses, sep='\t', header=F, quote='', comment.char='#');

# double-check that there's nothing unexpected in the first column
table(sapply(strsplit(gff_raw[,1], '_'), function(f) f[1]));
# NW 
# 1294055

# get rid of region feature lines
#  .buildGffLookup will also check for this so is actually unnecessary here
gff = subset(gff_raw, V3!='region');

# parse to lookup table, will be slow
#  if returnReducedGff=F takes ~10 min on 8-core 4GHz i7 iMac
#  if returnReducedGff=T will probably be ~3x longer
#  you can benchmark somewhat using the arg numToTest 
# by default, if write=T will write a file 'gff_lookup.tsv'
lookup0 = .buildGffLookup(gff, featureCol=3, sourceCol=2, metaCol=9, metaDelim=';', verbose=500, numToTest=NULL, write=T, returnReducedGff=T);

lookup = lookup0$lookup;
gffReduced = lookup0$gffReduced;
rm(lookup0); gc();

# make a gene level version of the gff, add a column for number of transcripts
gffReducedGenes = .makeGeneOnlyGff(gffReduced, featureCol=3, metaCol=9, delim=';', pattern='gene=', addTranscripts=T, verbose=500);

