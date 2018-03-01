#! /usr/bin/python
# ========================================================================================
# Jan.8 2016, Austin Hilliard, Stanford University Biology, Fernald lab
# parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py
# ========================================================================================
#
# This script will parse blastx results that include comments and write a new file
#  requires only base Python, see below for example input
#
# This script is an improvement on parseBestHitBlastxResultsWithComments.py
#  overall more flexible, matching is better, might be a little faster
#  parsing results from blastx on A.burtoni transcriptome that kept 3 best hits:
#   ~1 min running time on 4 GHz 8-core i7 iMac
#
# Probably should re-write as function definitions but not going to happen right now
#
# May modify how to decide best hit for a gene when transcripts have different top hits
#  would involve using homolog information across subject species
#  probably would store in dictionary, thus increasing memory overhead 
#  but I don't think it would increase running time 
#   depends on trade-off between how long it would take to read in the homolog info
#    vs time saved by streamlining the multiple top hits decision
#
# ========================================================================================
# Currently only one command line argument:
#
# sys.argv[1] is name of output file from blast run...
#
#	blastx -query H_burtoni_rna.fa \
#		   -db ./db/Drer_Olat_Onil_Trub_ENS_pep \
#		   -out H_burtoni_rna_blastx_FISH_ENS_top1 \
#		   -outfmt '7' \
#		   -max_target_seqs 1
#
#          -outfmt and -max_target_seqs must be set to '7' and '1', respectively
#
#-----------------------------------------------------------------------------------------
# Example portion of blast results:
# 
# # BLASTX 2.2.31+
# # Query: gi|565671690|ref|NM_001287404.1| Haplochromis burtoni corticoliberin-like (LOC102289224), mRNA
# # Database: ../_proteome/dRER_tRUB_oNIL_oLAT_gACU
# # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# # 1 hits found
# gi|565671690|ref|NM_001287404.1|	gi|348512264|ref|XP_003443663.1|	98.80	166	2	0	5	502	1	166	3e-94	282
# # BLASTX 2.2.31+
# # Query: gi|554802573|ref|XM_005912205.1| PREDICTED: Haplochromis burtoni uncharacterized LOC102291469 (LOC102291469), transcript variant X2, mRNA
# # Database: ../_proteome/dRER_tRUB_oNIL_oLAT_gACU
# # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# # 1 hits found
# gi|554802573|ref|XM_005912205.1|	gi|542193782|ref|XP_005471908.1|	93.89	131	8	0	232	624	1	131	8e-84	254
#
#-----------------------------------------------------------------------------------------
# There are 5 output files:
#
# inputFilename.allHits
# inputFilename.topByTranscript
# inputFilename.topByGene
# inputFilename.multiHitGenes
# inputFilename.noHits
#
#-------------------------------
#--inputFilename.allHits
#  Contains all hits from the blastx results
#  Same number of lines as there are non-header lines in the input, e.g.
#   $ wc -l inputFilename.allHits
#    should equal 
#   $ grep -v "^#" inputFilename | wc -l
#  Columns are gene id, query id, hit id, gene description, bit score of hit
#
#-------------------------------
#--inputFilename.topByTranscript
#  Contains only the best scoring hit for each query transcript
#  This is only file needed to move on to reciprocal blasting
#  Same number of lines as there are input queries minus number of no hits, e.g. if
#   $ wc -l inputFilename.noHits
#    returns n
#   $ wc -l inputFilename.topByTranscript
#    should equal
#   $ grep "^# Query" inputFilename | wc -l
#    plus n
#  Columns are same as inputFilename.allHits
#
#-------------------------------
#--inputFilename.topByGene
#  Contains the best hit for each gene
#  If all transcripts have same best hit or there's only one transcript:
#   that hit is reported here
#  If there are different best hits:
#   will check if more than half have the same one
#    if yes, that hit is reported here
#    if no, that gene will be excluded from this file 
#           hits for that gene will be reported in inputFilename.multiHitGenes
#  Columns are gene id, hit id, mean bit score, number of transcripts with hit, 
#               total number of transcripts, ratio of column 4 over column 5
#
#-------------------------------
#--inputFilename.multiHitGenes
#  Contains all hits for genes where at least one transcript had a different best hit
#  This includes genes where >half of the hits were the same
#   i.e. genes that still made it into inputFilename.topByGene
#  Columns are the same as inputFilename.allHits and inputFilename.topByTranscript
#
#-------------------------------
#--inputFilename.noHits
#  Contains transcripts that did not have a blastx hit
#  Columns are gene id, query id, gene description
#
# ========================================================================================
# Currently hardcoded variables that could be moved to arguments
#  outFilename = blastResultsFilename+'.transcriptsAndHitsByGenes'
#  geneNamePattern = 'LOC[0-9]{9}' ### now old, only worked when all ids were LOCs
#  geneNamePattern = '\(([^()]*)\)' ### returns anything enclosed in parentheses
#                                    ### ideally no nested parentheses
#  any other variables set in first couple sections below
#
# ========================================================================================

import sys, re, collections, time

#-----------------------------------------------------------------------------------------
# set variables used for identifying and splitting line types and finding gene names
#-----------------------------------------------------------------------------------------

# check for additional argument 
if len(sys.argv) < 3:
	headerType = 'ncbi'
else:
	headerType = sys.argv[2]
	
# regex to pull gene names from headers
print '---------------------'
if re.match('^-[Ee]', headerType):
	print 'Using header type:     Ensembl'
	geneNamePattern = 'gene:(.*) transcript:'
	queryLongNameField = 4
else:
	# as is will return anything between parentheses and might return multiple
	geneNamePattern = '\(([^()]*)\)'
	queryLongNameField = 3
	print 'Using header type:     NCBI'
print 'Gene name regex:       '+geneNamePattern
print '---------------------'


# rest of this section unlikely to change unless blast changes output format
# delimiters to split header lines with query info and hit lines
queryLineDelim = ' '
hitLineDelim = '\t'

# regex patterns for identifying line type
queryLinePattern = '^# Query'
allHeaderLinePattern = '^#'
noHitLinePattern = '^# 0 hits found'

# fields in query line holding query id and beginning of full gene description
#  assumes line was split on queryLineDelim
queryIdField = 2
#queryLongNameField = 3

# field in hit line containing hit id, assumes line was split on hitLineDelim
hitIdField = 1

#-----------------------------------------------------------------------------------------
# set input/output file names and open for read/write
#-----------------------------------------------------------------------------------------
# set input and output file names
blastResultsFilename = sys.argv[1]
outnameAll = blastResultsFilename+'.allHits'
outnameTopByGene = blastResultsFilename+'.topByGene'
outnameTopByTran = blastResultsFilename+'.topByTranscript'
outnameMultiGene = blastResultsFilename+'.multiHitGenes'
outnameNoHits = blastResultsFilename+'.noHits'

# show an example of how query header lines will be parsed
with open(blastResultsFilename, 'r') as input:
	for line in input:
		if re.match(queryLinePattern, line):
			sLine = line.split(queryLineDelim)
			queryTmp = re.findall(geneNamePattern, line)
			print 'Query parsing example'
			print 'raw header:            '+line
			print 'queryId:               '+sLine[queryIdField]
			print 'queryLongName:         '+' '.join(sLine[queryLongNameField:]).strip()
			print 'queryParentGene:       '+queryTmp[len(queryTmp)-1]
			break

# open files for read/write
blastResults = open(blastResultsFilename, 'r')
outfile = open(outnameAll, 'w')
outfileTopByGene = open(outnameTopByGene, 'w')
outfileTopByTran = open(outnameTopByTran, 'w')
outfileMultiHitGenes = open(outnameMultiGene, 'w')
outfileNoHitGenes = open(outnameNoHits, 'w')

# print message to console
#  include regex that will identify gene names, input filename, current time
#print '---------------------'
#print 'Gene name regex:       '+geneNamePattern
print '---------------------'
print 'Parsing input file:    '+blastResultsFilename
t = time.ctime().split()[1:4]
print t[0]+'.'+t[1]+' '+t[2]
print '---------------------'

#-----------------------------------------------------------------------------------------
# process blastResults one line at a time to build geneDict
#-----------------------------------------------------------------------------------------

# dict to hold info per gene
#  keys will be gene names, as recognized by geneNamePattern
#  values will be lists containing tuples with query/hit ids, gene description, bit score
geneDict = dict()

# count total number of transcripts/queries and how many did not have a blastx hit
tranCount, noHitCount = 0, 0

for line in blastResults:
#	--------------------------------------------------------------------------------------
#	check if line has query info
#	--------------------------------------------------------------------------------------
	if re.match(queryLinePattern, line):
	
		# update transcript count
		tranCount += 1
		
		# split line on space, get query id and gene description
		sLine = line.split(queryLineDelim)
		queryId = sLine[queryIdField]
		queryLongName = ' '.join(sLine[queryLongNameField:]).strip()
		
		# extract gene name from line
		#  regex may hit multiple parts of string, assumes gene name is the last hit
		queryTmp = re.findall(geneNamePattern, line)
		queryParentGene = queryTmp[len(queryTmp)-1]
			
#	--------------------------------------------------------------------------------------				
#	check if current query had no blastx hits
#	--------------------------------------------------------------------------------------	
	elif re.match(noHitLinePattern, line):
	
		# update no hit transcript count
		noHitCount += 1
		
		# write tab-delimited line to noHits output file
		#  includes gene and transcript ids, full gene description
		writeNoHit = queryParentGene+'\t'+queryId+'\t'+queryLongName+'\n'
		outfileNoHitGenes.write(writeNoHit)

#	--------------------------------------------------------------------------------------				
#	if line isn't another header line split on tabs and get info on the hit
#	--------------------------------------------------------------------------------------
	elif not re.match(allHeaderLinePattern, line):
	
		# get hit id and score, assumes score is in final field of split line
		sLine = line.split(hitLineDelim)
		hitId = sLine[hitIdField]
		bitScore = sLine[len(sLine)-1].strip()
		hitInfo = (queryId, hitId, queryLongName, bitScore)
		
		# check whether gene name is already a key in the dict
		#  if yes append hitInfo tuple to existing list
		#  if no make a new key and add new list containing hitInfo
		if queryParentGene in geneDict.keys():
			geneDict[queryParentGene].append(hitInfo)
		else:
			geneDict[queryParentGene] = [hitInfo]
		
#	--------------------------------------------------------------------------------------				
#	--- place to include handling of other header lines ---
#	--------------------------------------------------------------------------------------
#	else:
#
#

# print information to console:
#  names of output files and current time
print 'Writing output files:  '+outnameAll
t = time.ctime().split()[1:4]
print t[0]+'.'+t[1]+' '+t[2]+'        '+outnameTopByTran
print '                       '+outnameTopByGene
print '                       '+outnameMultiGene
print '                       '+outnameNoHits
print '---------------------'

#-----------------------------------------------------------------------------------------
# iterate through geneDict to write output files
#-----------------------------------------------------------------------------------------

# initialize list to keep track of how many transcripts each gene has
geneTranCount = []

# count genes where: 
#  all transcripts had the same best hit
#  at least one transcript had different best hit than others
#  at least one transcript had different best hit than others and >half were the same
numSingleBestGenes, numMultiHitGenes, numMultiHitGenesThatPassRatio = 0, 0, 0

# 'gene' will be a tuple where gene[0] is the key, i.e. gene name
#  gene[1] will be a list containing a tuple for each hit
for gene in geneDict.items():

	# dict to hold info per transcript 
	#  keys will be queryId
	#  values will be lists containing tuples with hitId, queryLongName, bitScore
	isoDict = dict()
	
#	--------------------------------------------------------------------------------------				
#	iterate through hits to write main output file and build isoDict for this gene
#	--------------------------------------------------------------------------------------
	# each hitInfo is a tuple with ids, gene description, bit score
	#  should be one for each gene isoform/transcript
	for hitInfo in gene[1]:
		
#	    ----------------------------------------------------------------------------------			
#	    write line in main output file - all hits
#	    ----------------------------------------------------------------------------------
		# tab-delimited line with gene name and hitInfo
		toWrite = gene[0]+'\t'+'\t'.join(hitInfo)+'\n'
		outfile.write(toWrite)
	
#	    ----------------------------------------------------------------------------------			
#	    build isoDict to store hits for this gene by transcript
#	    ----------------------------------------------------------------------------------
		# get info for each isoform/transcript (queryId)
		isoform = hitInfo[0]
		isoInfo = hitInfo[1:]
		
		# check whether isoform is already a key
		#  if yes append isoInfo to existing list
		#  if no make a new key and add new list containing isoInfo
		if isoform in isoDict.keys():
			isoDict[isoform].append(isoInfo)
		else:
			isoDict[isoform] = [isoInfo]
	
#	--------------------------------------------------------------------------------------		
#	iterate through isoDict to write topByTranscript output file
#   --------------------------------------------------------------------------------------
	# topList will be used in this section 
	# topIds and topScores will only be used later when writing topByGene output file
	#  only if all transcripts for the gene don't have the same top hit
	#  (in A.burtoni this is only ~5% of genes)
	topIds, topScores, topList = [], [], []
	
	# 'iso' will be a tuple where iso[0] is the key, i.e. isoform/transcript/query id
	#  iso[1] will be a list containing a tuple for each hit
	for iso in isoDict.items():
	
		# build list of top hit tuples for transcripts of this gene
		#  isoTopHit is tuple holding hitId, queryLongName, bitScore for top hit
		isoTopHit = iso[1][0]
		topList.append(isoTopHit)
		
		# build separate lists of top hitIds and bitScores (see below)
		#  use hitIdField-1 because tuples are 1 item shorter than in geneDict
		#  assumes bit score is last entry in tuple
		topIds.append(isoTopHit[hitIdField-1])
		topScores.append(isoTopHit[len(isoTopHit)-1])
		
#	    ----------------------------------------------------------------------------------			
#	    write line in topByTranscript output file
#	    ----------------------------------------------------------------------------------
		# tab-delimited line with ids, gene long name, and bit score 
		toWriteTopIso = gene[0]+'\t'+iso[0]+'\t'+'\t'.join(isoTopHit)+'\n'
		outfileTopByTran.write(toWriteTopIso)
	
#	--------------------------------------------------------------------------------------		
#	examine top hits to write topByGene and multiHitGenes output files
#   --------------------------------------------------------------------------------------
	# convert bit scores of top hits to floats before computing mean score
	topScores = map(float, topScores)
	
	# check if only one transcript for this gene or if all had same top hit
	#  if yes write a line in the topByGene output file
	if len(set(topIds)) == 1:
	
		# update count
		numSingleBestGenes += 1
	
		# get mean score and number of isoforms/transcripts 
		meanScore = sum(topScores) / len(topScores)
		meanScore = str(round(meanScore, 2))
		numIso = str(len(topScores))
		
		# numHit is number of isoforms with this top hit, here it's all by definition 
		numHit = numIso
		
#	    ----------------------------------------------------------------------------------			
#	    write line in topByGene output file
#	    ---------------------------------------------------------------------------------- 
		# tab-delimited line with ids, mean bit score, number of isoforms, and 1
		#  the final entry is the ratio numHit/numIso, here it's 1 by definition
		toWriteTopGene = '\t'.join([gene[0], topIds[0], meanScore, numHit, numIso, '1\n'])
		outfileTopByGene.write(toWriteTopGene)

#	--------------------------------------------------------------------------------------		
#	if not all top hits are the same need to decide what to report
#   --------------------------------------------------------------------------------------
	# Currently uses hit counts but may not be necessary with homolog information
	# Instead could check if different hits are:
	#  1: from the same species
	#	   check if they're from same gene 
	#       if yes take that gene as best hit
	#      if not from same gene separate into other file of ambiguous hits
	#  2: not from the same species
	#      check if they're from same gene
	#       if yes take more numerous hit, check if it maps to human
	#        if no check if other id maps to human
	else:
		
		# update count
		numMultiHitGenes += 1
		
#	    ----------------------------------------------------------------------------------		
#	   	count hits for isoforms/transcripts and get scores for most numerous
#	   	----------------------------------------------------------------------------------
		# create hitCounts dictionary
		#  keys are hit ids, values are number of isoforms that had the hit
		hitCounts = collections.Counter(topIds)
		
		# convert the isoform counts to floats to compute ratio
		vals = map(float, hitCounts.values())
		
		# ratio of total top hits that are the most numerous hit
		#  e.g. a gene has 4 isoforms, 3 had topHit1, 1 had topHit2 -> ratio = 0.75
		ratio = max(vals) / sum(vals)
		
#	    ----------------------------------------------------------------------------------		
#	   	check if more than half of isoforms had the same top hit
#	   	----------------------------------------------------------------------------------
		# If yes take that hit as best for the gene and write to topByGene output file
		# Not sure this is best way to filter
		#  probably should check if the different hits are same gene
		#  if yes would take more numerous hit or one that maps to human
		if ratio > 0.5:
			
			# update count
			numMultiHitGenesThatPassRatio += 1
			
			# get key/hitId with value that matches the max count
			topCount = int(max(vals))
			bestHit = [h for h, ct in hitCounts.items() if ct == topCount][0]
			
			# iterate through topIds and topScores (from above section)
			#  need scores for top hit to compute mean
			# may be better to do with dictionary
			#  but would need to be careful about duplicate keys getting collapsed
			bestScores = []
			for sc in range(len(topIds)):
			
				# if current id is the best hit take corresponding score
				if topIds[sc] == bestHit:
					bestScores.append(topScores[sc])
			bestScoresMean = str(sum(bestScores) / len(bestScores))
			
#	    	------------------------------------------------------------------------------		
#	    	write line in topByGene output file
#	    	------------------------------------------------------------------------------
			# convert isoform hit counts to strings for writing
			numIso = str(int(sum(vals)))
			numHit = str(topCount)
			ratio = str(round(ratio, 2))
			
			# tab-delimited line with ids, mean bit score, number and ratio of isoforms 
			toWriteRatio = '\t'.join([gene[0], bestHit, bestScoresMean, numHit, numIso, ratio, '\n'])
			outfileTopByGene.write(toWriteRatio)
			
#	    ----------------------------------------------------------------------------------		
#	   	write line in multiHitGenes output file
#	   	----------------------------------------------------------------------------------
		# multiHitGenes is subset of the main output file with all hits
		# some genes here may not make it into topByGene file
		#  currently depends on more than half of isoforms having the same hit
		for iso in isoDict.items():
			toWriteMulti = gene[0]+'\t'+iso[0]+'\t'+'\t'.join(iso[1][0])+'\n'
			outfileMultiHitGenes.write(toWriteMulti)
	
	# update list to keep track of how many transcripts each gene has
	geneTranCount.append(len(isoDict))

#-----------------------------------------------------------------------------------------
# finish by closing files and getting info about number of transcripts per gene
#-----------------------------------------------------------------------------------------

# close all open files
blastResults.close()
outfile.close()
outfileTopByTran.close()
outfileTopByGene.close()
outfileMultiHitGenes.close()	
outfileNoHitGenes.close()

# convert transcript counts per gene to floats for computing mean
geneTranCount = map(float, geneTranCount)
geneTranCountMean = round(sum(geneTranCount) / len(geneTranCount), 2)

# get median, min, max, number of transcripts per gene
tmp = geneTranCount
tmp.sort()
geneTranCountMedian = int(tmp[len(tmp)/2])
minTran = int(min(geneTranCount))
maxTran = int(max(geneTranCount))

# print information to console
#  numbers of total genes and transcripts in the input, transcripts with no hit
#  stats on numbers of transcripts per gene
#  number of genes where 
#   all transcripts had the same best hit
#   at least one transcript had a different best hit
#   at least one transcript had a different best hit but >half were the same
#  current time
print 'Total input genes:     '+str(len(geneDict))
print '      transcripts:     '+str(tranCount)
print '---------------------'
print 'No hit transcripts:    '+str(noHitCount)
print '---------------------'
print 'Transcripts/gene'
print ' mean, median:         '+str(geneTranCountMean)+', '+str(geneTranCountMedian)
print ' min, max:             '+str(minTran)+', '+str(maxTran)
print '---------------------'
print 'Genes where all'  
print ' transcripts had'
print ' same best hit:        '+str(numSingleBestGenes)
print '---------------------'
print 'Multi-hit genes:       '+str(numMultiHitGenes)
print '---------------------'
print 'Multi-hit genes'
print ' where >half were'
print ' the same:             '+str(numMultiHitGenesThatPassRatio)
print '---------------------'
t = time.ctime().split()[1:4]
print 'Done\n'+t[0]+'.'+t[1]+' '+t[2]
print '---------------------'

##########################################################################################
# Example - blastx results
##########################################################################################
# 
# BIO-985AEBE2A677:_annotationsDec2015 abseq$ /Volumes/fishstudies-1/_scripts/parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out 
# ---------------------
# Using header type:     NCBI
# Gene name regex:       \(([^()]*)\)
# ---------------------
# Query parsing example
# raw header:            # Query: gi|930722691|ref|XM_014337201.1| PREDICTED: Haplochromis burtoni transmembrane protein 56-B-like (LOC102288583), transcript variant X4, mRNA
# 
# queryId:               gi|930722691|ref|XM_014337201.1|
# queryLongName:         PREDICTED: Haplochromis burtoni transmembrane protein 56-B-like (LOC102288583), transcript variant X4, mRNA
# queryParentGene:       LOC102288583
# ---------------------
# Parsing input file:    hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out
# Feb.12 15:34:10
# ---------------------
# Writing output files:  hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.allHits
# Feb.12 15:35:01        hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.topByTranscript
#                        hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.topByGene
#                        hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.multiHitGenes
#                        hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.noHits
# ---------------------
# Total input genes:     25574
#       transcripts:     47807
# ---------------------
# No hit transcripts:    712
# ---------------------
# Transcripts/gene
#  mean, median:         1.84, 1
#  min, max:             1, 33
# ---------------------
# Genes where all
#  transcripts had
#  same best hit:        24418
# ---------------------
# Multi-hit genes:       1156
# ---------------------
# Multi-hit genes
#  where >half were
#  the same:             684
# ---------------------
# Done
# Feb.12 15:35:02
# ---------------------
#
##########################################################################################
# Example - reciprocal blastp results, use -e flag since headers are ensembl format
##########################################################################################
# 
# BIO-985AEBE2A677:_annotationsDec2015 abseq$ /Volumes/fishstudies-1/_scripts/parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py topByTranscript.recpBlastP_query.092915_hBurtoni.out -e
# ---------------------
# Using header type:     Ensembl
# Gene name regex:       gene:(.*) transcript:
# ---------------------
# Query parsing example
# raw header:            # Query: ENSDARP00000136057.1 pep:known chromosome:GRCz10:2:36281866:36282399:-1 gene:ENSDARG00000099434.1 transcript:ENSDART00000163282.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene
# 
# queryId:               ENSDARP00000136057.1
# queryLongName:         chromosome:GRCz10:2:36281866:36282399:-1 gene:ENSDARG00000099434.1 transcript:ENSDART00000163282.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene
# queryParentGene:       ENSDARG00000099434.1
# ---------------------
# Parsing input file:    topByTranscript.recpBlastP_query.092915_hBurtoni.out
# Feb.12 15:42:17
# ---------------------
# Writing output files:  topByTranscript.recpBlastP_query.092915_hBurtoni.out.allHits
# Feb.12 15:42:45        topByTranscript.recpBlastP_query.092915_hBurtoni.out.topByTranscript
#                        topByTranscript.recpBlastP_query.092915_hBurtoni.out.topByGene
#                        topByTranscript.recpBlastP_query.092915_hBurtoni.out.multiHitGenes
#                        topByTranscript.recpBlastP_query.092915_hBurtoni.out.noHits
# ---------------------
# Total input genes:     23256
#       transcripts:     24051
# ---------------------
# No hit transcripts:    23
# ---------------------
# Transcripts/gene
#  mean, median:         1.03, 1
#  min, max:             1, 4
# ---------------------
# Genes where all
#  transcripts had
#  same best hit:        22788
# ---------------------
# Multi-hit genes:       468
# ---------------------
# Multi-hit genes
#  where >half were
#  the same:             5
# ---------------------
# Done
# Feb.12 15:42:45
# ---------------------
