###
rm(list=ls());
options(stringsAsFactors=F);
setwd('/Users/abseq/Documents/_annotationsDec2015/');
library(biomaRt);

#############################################################
### load files
#############################################################

# this should be file generated from parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py
parsedBlastxResults = 'hBurtoni09252015_blastx_dRERtRUBoNILoLATgACU.out.topByTranscript';

# this should be file generated from parseBestHitBlastxResultsWithComments_betterMatchingJan2016.py
reciprocalBlastpResults = 'topByTranscript.recpBlastP_query.092915_hBurtoni.out.topByTranscript';

# this should be file generated from parseGffFileToLookup_Feb2016.R
AburtoniLookup = 'ref_AstBur1.0_scaffolds_translated.gff3_lookup_Feb2016_fromR';

# read in parsed burtoni blastx results
x = read.table(parsedBlastxResults, header=F, sep='\t', colClasses='character', quote='');

# read in reciprocal blastp results 
rx = read.table(reciprocalBlastpResults, header=F, sep='\t', colClasses='character', quote='');

# read in lookup table to translate between burtoni protein, rna, and gene ids
ablookup = read.table(AburtoniLookup, header=T, sep='\t', fill=T, colClasses='character');

#############################################################
### compare results from reciprocal blastp to original blastx
#############################################################

blastHitCol = 3;
blastQueryCol = 2;
abDescriptionCol = 4;

# check burtoni protein hits from reciprocal blastp
#  first check that '^NP|^XP' regex will be ok to use
table(sapply(strsplit(ablookup$protein_id,'_'), function(f) f[1]));
# NP    XP 
# 50 44603
abProtRegex = '^NP|^XP';
rxBurtoniProteinHits = grep(abProtRegex, unlist(strsplit(rx[, blastHitCol], '|', fixed=T)), value=T);
# check that all burtoni protein hits from reciprocal blast are in the lookup table
sum(rxBurtoniProteinHits %in% ablookup$protein_id) == nrow(rx);
# [1] TRUE
rx[, blastHitCol] = rxBurtoniProteinHits;

# check burtoni transcript queries from original blastx
#  first check that '^NM|^XM|^XR' regex will be ok to use
table(sapply(strsplit(ablookup$transcript_id,'_'), function(f) f[1]));
# NM    XM    XR 
# 50 44603  3154
# get burtoni transcripts from original blastx
abTranscriptRegex = '^NM|^XM|^XR';
burtoniTranscripts = grep(abTranscriptRegex, unlist(strsplit(x[, blastQueryCol], '|', fixed=T)), value=T);
# check that all burtoni transcript ids from blastx are in the lookup table
sum(burtoniTranscripts %in% ablookup$transcript_id) == nrow(x);
# [1] TRUE
x[, blastQueryCol] = burtoniTranscripts;

# make new version of ablookup that's reordered and filtered to match the blastx results
xLookup0 = ablookup[match(x[, blastQueryCol], ablookup$transcript_id), ];

# add column with blastx hits
#xLookup = as.data.frame(cbind(xLookup0, oBlastxHit=x[, blastHitCol]));
xLookup = as.data.frame(cbind(xLookup0, oBlastxHit=x[, blastHitCol], abDescription=x[, abDescriptionCol]));


# ### separate genes with ncbi-assigned symbols
# ###  will get human homologs for these genes regardless of reciprocal blast results anyway
# xLookup_sym = subset(xLookup0, !grepl('^LOC', gene));
# xLookup_loc = subset(xLookup0, grepl('^LOC', gene));

# add blastp results
lookup = .addRecipBlastpToBlastxLookup(xLookup, rx);
lookupGeneSplit = split(lookup, lookup$gene);

# make a list with mart for each fish species
host = 'www.ensembl.org';
bmart = 'ENSEMBL_MART_ENSEMBL';
speciesMarts = list(DAR=useMart(biomart=bmart, dataset='drerio_gene_ensembl', host=host),
                    ORL=useMart(biomart=bmart, dataset='olatipes_gene_ensembl', host=host),
                    ONI=useMart(biomart=bmart, dataset='oniloticus_gene_ensembl', host=host),
                    TRU=useMart(biomart=bmart, dataset='trubripes_gene_ensembl', host=host),
                    GAC=useMart(biomart=bmart, dataset='gaculeatus_gene_ensembl', host=host));
attributes = c('ensembl_peptide_id','ensembl_gene_id','description','external_gene_name','hsapiens_homolog_ensembl_gene');
filters = 'ensembl_peptide_id';

# try converting fish protein ids to fish genes, then human genes
#  hsENSdf is a list of dataframes, one per fish species
#   rows are fish proteins, columns are attributes defined above
hsENSdf = .convertMultiSpeciesENStoHuman(unique(lookup$oBlastxHit), speciesMarts, attributes, filters);

# try mapping human ensembl ids to entrez ids
hsENSvec = unique(unlist(sapply(hsENSdf, function(f) f$hsapiens_homolog_ensembl_gene)));
hsENSvec = grep('^ENSG', hsENSvec, value=T);
# hsEntrez is a 2 column dataframe, ensembl ids and entrez ids
# there will be duplicates in each column but downstream functions account for this
hsEntrez = .convertHsEnsemblToEntrez(hsENSvec);

# collapse lookup (where each row is a burtoni transcript) to dataframe where each row is a burtoni gene
an0 = .processLookupToGeneAnnotations(lookupGeneSplit, hsENSdf);

# rows in an0$pass with NAs are genes that I couldn't make a decision on
#  these rows correspond to the elements of an0$leftNAs
# recipPass should be $pass from .processLookupToGeneAnnotations with all NA rows removed
recipPass0 = an0$pass[!is.na(an0$pass$gene), ];
# table(recipPass$passType)
# all_same_blastx  all_same_mapped_gene   same_symbol    same_symbol_strip1  single twothirds_same_blastx 
# 7943             671                    13             5                   11941  14

recipFail = an0$fail;
# table(recipFail$failType)
# multiple_fail_blastp   single_fail_blastp   single_ncRNA  single_no_blastp_hit 
# 1939                   2809                 1423           7 

# add gene ensembl ids, descriptions, symbols for original blastx hits
#  also add homologous human ensembl and entrez ids
recipPass = .addGenesToRecipPass(recipPass0, hsENSdf, hsEntrez);

sum(is.na(recipPass$hsaHomologGene));
# [1] 3855
sum(grepl('^N|A$', recipPass$hsaHomologGene));
# [1] 0
sum(recipPass$hsaHomologGene=='', na.rm=T);
# [1] 0
sum(is.na(recipPass$hsaHomologGene) & !grepl('^LOC', recipPass$gene), na.rm=T);
# [1] 364

sum(is.na(recipPass$hsaHomologEntrez));
# [1] 3893
sum(grepl('^N|A$', recipPass$hsaHomologEntrez));
# [1] 0
sum(recipPass$hsaHomologEntrez=='', na.rm=T);
# [1] 0
sum(is.na(recipPass$hsaHomologEntrez) & !grepl('^LOC', recipPass$gene), na.rm=T);
# [1] 373

# identify genes that passed reciprocal blastp, have a ncbi gene symbol in the burtoni annotations,
#  but no mapped human homolog
# use gene symbol to get human homolog ensembl and entrez ids 
recipPass2 = .addMissingHsHomologsToRecipPassSymbols(recipPass);

sum(is.na(recipPass2$hsaHomologGene));
# [1] 3348
sum(grepl('^N|A$', recipPass2$hsaHomologGene));
# [1] 0
sum(recipPass2$hsaHomologGene=='', na.rm=T);
# [1] 0
sum(is.na(recipPass2$hsaHomologGene) & !grepl('^LOC', recipPass2$gene), na.rm=T);
# [1] 1

sum(is.na(recipPass2$hsaHomologEntrez));
# [1] 3370
sum(grepl('^N|A$', recipPass2$hsaHomologEntrez));
# [1] 0
sum(recipPass2$hsaHomologEntrez=='', na.rm=T);
# [1] 0
sum(is.na(recipPass2$hsaHomologEntrez) & !grepl('^LOC', recipPass2$gene), na.rm=T);
# [1] 15

subset(recipPass2, is.na(hsaHomologEntrez) & !grepl('^LOC', gene))
#          Dbxref    gene  rna_num  transcript_id     protein_id           oBlastxHit     rBlastpHit recipPass       recipType     oBlastxHitGene
# 504   102305428   apoc2 rna30166 XM_005938231.2 XP_005938293.1 ENSTRUP00000007424.1 XP_005938293.1      TRUE          single ENSTRUG00000003166
# 661   102313712   asip2  rna7311 XM_005918395.2 XP_005918457.1 ENSONIP00000004516.1 XP_005918457.1      TRUE          single ENSONIG00000003590
# 1910  102302800  crocc2  rna8910 XM_014332844.1 XP_014188330.1 ENSONIP00000024290.1 XP_014188330.1      TRUE all_same_blastx ENSONIG00000019301
# 3148  102309351  fignl2  rna3340 XM_005914824.1 XP_005914886.1 ENSONIP00000015020.1 XP_005914886.1      TRUE all_same_blastx ENSONIG00000011929
# 3680  102294189    grk1 rna23782 XM_005932677.2 XP_005932739.1 ENSONIP00000003977.1 XP_005932739.1      TRUE          single ENSONIG00000003174
# 4136  102300458  inafm2  rna5703 XM_005916953.1 XP_005917015.1 ENSDARP00000117627.1 XP_005917015.1      TRUE          single ENSDARG00000070423
# 4489  102296014  kif28p  rna3572 XM_005915237.1 XP_005915299.1 ENSONIP00000021611.1 XP_005915299.1      TRUE          single ENSONIG00000017138
# 20614 102308536   muc19  rna4331 XM_014329828.1 XP_014185314.1 ENSONIP00000012378.1 XP_014185314.1      TRUE          single ENSONIG00000009851
# 21202 102306246   ofcc1 rna34901 XM_014340114.1 XP_014195600.1 ENSONIP00000018300.1 XP_014195600.1      TRUE all_same_blastx ENSONIG00000014542
# 21557 102293125   pgbd5  rna2568 XM_005914199.2 XP_005914261.1 ENSONIP00000000692.1 XP_005914261.1      TRUE          single ENSONIG00000000545
# 22693 102312596 rnpepl1 rna41080 XM_005947846.2 XP_005947908.1 ENSONIP00000021791.1 XP_005947908.1      TRUE          single ENSONIG00000017277
# 23172 102311615  shank3 rna25420 XM_014337401.1 XP_014192887.1 ENSONIP00000018731.1 XP_014192887.1      TRUE all_same_blastx ENSONIG00000014881
# 23652 102301214  spata1 rna13105 XM_014333904.1 XP_014189390.1 ENSONIP00000009318.1 XP_014189390.1      TRUE all_same_blastx ENSONIG00000007390
# 24739 102291833   tstd3 rna14156 XM_014334228.1 XP_014189714.1 ENSONIP00000000142.1 XP_014189714.1      TRUE all_same_blastx ENSONIG00000000115
# 24941 106632346   umad1 rna40136 XM_014329122.1 XP_014184608.1 ENSONIP00000024739.1 XP_014184608.1      TRUE all_same_blastx ENSONIG00000019652
#                                                                              oBlastxHitDescription    oBlastxHitSym                  hsaHomologGene hsaHomologEntrez
# 504   Takifugu rubripes apolipoprotein C-II (apoc2), mRNA. [Source:RefSeq mRNA;Acc:NM_001032718]             <NA>                 ENSG00000234906             <NA>
#   661                                                                                         <NA>             <NA>                            <NA>             <NA>
#   1910  ciliary rootlet coiled-coil, rootletin family member 2 [Source:HGNC Symbol;Acc:HGNC:51677]           CROCC2                 ENSG00000226321             <NA>
#   3148                                                                                        <NA>             <NA>                 ENSG00000261308             <NA>
#   3680                   G protein-coupled receptor kinase 1 a [Source:ZFIN;Acc:ZDB-GENE-050823-1]   grk1a (1 of 2) ENSG00000185974,ENSG00000281988             <NA>
#   4136                                            zgc:153157 [Source:ZFIN;Acc:ZDB-GENE-060825-224]       zgc:153157                 ENSG00000259330             <NA>
#   4489                                        si:rp71-84d9.1 [Source:ZFIN;Acc:ZDB-GENE-131121-579]   si:rp71-84d9.1                 ENSG00000223519             <NA>
#   20614                                                                                       <NA>             <NA>                 ENSG00000205592             <NA>
#   21202                          orofacial cleft 1 candidate 1 [Source:HGNC Symbol;Acc:HGNC:21017]            OFCC1                 ENSG00000181355             <NA>
#   21557              piggyBac transposable element derived 5 [Source:ZFIN;Acc:ZDB-GENE-040724-229]            pgbd5                 ENSG00000177614             <NA>
#   22693       arginyl aminopeptidase (aminopeptidase B)-like 1 [Source:HGNC Symbol;Acc:HGNC:10079]          RNPEPL1                 ENSG00000142327             <NA>
#   23172           SH3 and multiple ankyrin repeat domains 3a [Source:ZFIN;Acc:ZDB-GENE-060503-369]          shank3a                 ENSG00000251322             <NA>
#   23652                                       si:dkey-205k8.5 [Source:ZFIN;Acc:ZDB-GENE-141216-11]  si:dkey-205k8.5                 ENSG00000122432             <NA>
#   24739                                     si:dkey-165n16.1 [Source:ZFIN;Acc:ZDB-GENE-100922-224] si:dkey-165n16.1                 ENSG00000228439             <NA>
#   24941                                     si:dkey-119m17.2 [Source:ZFIN;Acc:ZDB-GENE-070705-289] si:dkey-119m17.2                 ENSG00000219545             <NA>






recipFail2 = .addGenesToFailSymbols(recipFail);
# edit colnames so dataframe can be combined with rbind
names(recipPass2)[names(recipPass2)=='passType'] = 'recipType';
names(recipFail2)[names(recipFail2)=='failType'] = 'recipType';
# 
anno = as.data.frame(rbind(recipPass2, recipFail2));

abSyms = !grepl('^LOC', anno$gene);
goodHsaENS = grepl('ENS', anno$hsaHomologGene);
goodHsaEntrez = grepl('[0-9]', anno$hsaHomologEntrez);
toFix = (!goodHsaENS | !goodHsaEntrez) & abSyms;

anno[toFix, ];
#          Dbxref    gene  rna_num  transcript_id     protein_id           oBlastxHit     rBlastpHit recipPass       recipType     oBlastxHitGene
# 504   102305428   apoc2 rna30166 XM_005938231.2 XP_005938293.1 ENSTRUP00000007424.1 XP_005938293.1      TRUE          single ENSTRUG00000003166
# 661   102313712   asip2  rna7311 XM_005918395.2 XP_005918457.1 ENSONIP00000004516.1 XP_005918457.1      TRUE          single ENSONIG00000003590
# 1910  102302800  crocc2  rna8910 XM_014332844.1 XP_014188330.1 ENSONIP00000024290.1 XP_014188330.1      TRUE all_same_blastx ENSONIG00000019301
# 3148  102309351  fignl2  rna3340 XM_005914824.1 XP_005914886.1 ENSONIP00000015020.1 XP_005914886.1      TRUE all_same_blastx ENSONIG00000011929
# 3680  102294189    grk1 rna23782 XM_005932677.2 XP_005932739.1 ENSONIP00000003977.1 XP_005932739.1      TRUE          single ENSONIG00000003174
# 4136  102300458  inafm2  rna5703 XM_005916953.1 XP_005917015.1 ENSDARP00000117627.1 XP_005917015.1      TRUE          single ENSDARG00000070423
# 4489  102296014  kif28p  rna3572 XM_005915237.1 XP_005915299.1 ENSONIP00000021611.1 XP_005915299.1      TRUE          single ENSONIG00000017138
# 20614 102308536   muc19  rna4331 XM_014329828.1 XP_014185314.1 ENSONIP00000012378.1 XP_014185314.1      TRUE          single ENSONIG00000009851
# 21202 102306246   ofcc1 rna34901 XM_014340114.1 XP_014195600.1 ENSONIP00000018300.1 XP_014195600.1      TRUE all_same_blastx ENSONIG00000014542
# 21557 102293125   pgbd5  rna2568 XM_005914199.2 XP_005914261.1 ENSONIP00000000692.1 XP_005914261.1      TRUE          single ENSONIG00000000545
# 22693 102312596 rnpepl1 rna41080 XM_005947846.2 XP_005947908.1 ENSONIP00000021791.1 XP_005947908.1      TRUE          single ENSONIG00000017277
# 23172 102311615  shank3 rna25420 XM_014337401.1 XP_014192887.1 ENSONIP00000018731.1 XP_014192887.1      TRUE all_same_blastx ENSONIG00000014881
# 23652 102301214  spata1 rna13105 XM_014333904.1 XP_014189390.1 ENSONIP00000009318.1 XP_014189390.1      TRUE all_same_blastx ENSONIG00000007390
# 24739 102291833   tstd3 rna14156 XM_014334228.1 XP_014189714.1 ENSONIP00000000142.1 XP_014189714.1      TRUE all_same_blastx ENSONIG00000000115
# 24941 106632346   umad1 rna40136 XM_014329122.1 XP_014184608.1 ENSONIP00000024739.1 XP_014184608.1      TRUE all_same_blastx ENSONIG00000019652
#                                                                              oBlastxHitDescription    oBlastxHitSym                  hsaHomologGene hsaHomologEntrez
# 504   Takifugu rubripes apolipoprotein C-II (apoc2), mRNA. [Source:RefSeq mRNA;Acc:NM_001032718]             <NA>                 ENSG00000234906             <NA>
#   661                                                                                         <NA>             <NA>                            <NA>             <NA>
#   1910  ciliary rootlet coiled-coil, rootletin family member 2 [Source:HGNC Symbol;Acc:HGNC:51677]           CROCC2                 ENSG00000226321             <NA>
#   3148                                                                                        <NA>             <NA>                 ENSG00000261308             <NA>
#   3680                   G protein-coupled receptor kinase 1 a [Source:ZFIN;Acc:ZDB-GENE-050823-1]   grk1a (1 of 2) ENSG00000185974,ENSG00000281988             <NA>
#   4136                                            zgc:153157 [Source:ZFIN;Acc:ZDB-GENE-060825-224]       zgc:153157                 ENSG00000259330             <NA>
#   4489                                        si:rp71-84d9.1 [Source:ZFIN;Acc:ZDB-GENE-131121-579]   si:rp71-84d9.1                 ENSG00000223519             <NA>
#   20614                                                                                       <NA>             <NA>                 ENSG00000205592             <NA>
#   21202                          orofacial cleft 1 candidate 1 [Source:HGNC Symbol;Acc:HGNC:21017]            OFCC1                 ENSG00000181355             <NA>
#   21557              piggyBac transposable element derived 5 [Source:ZFIN;Acc:ZDB-GENE-040724-229]            pgbd5                 ENSG00000177614             <NA>
#   22693       arginyl aminopeptidase (aminopeptidase B)-like 1 [Source:HGNC Symbol;Acc:HGNC:10079]          RNPEPL1                 ENSG00000142327             <NA>
#   23172           SH3 and multiple ankyrin repeat domains 3a [Source:ZFIN;Acc:ZDB-GENE-060503-369]          shank3a                 ENSG00000251322             <NA>
#   23652                                       si:dkey-205k8.5 [Source:ZFIN;Acc:ZDB-GENE-141216-11]  si:dkey-205k8.5                 ENSG00000122432             <NA>
#   24739                                     si:dkey-165n16.1 [Source:ZFIN;Acc:ZDB-GENE-100922-224] si:dkey-165n16.1                 ENSG00000228439             <NA>
#   24941                                     si:dkey-119m17.2 [Source:ZFIN;Acc:ZDB-GENE-070705-289] si:dkey-119m17.2                 ENSG00000219545             <NA>

anno$hsaHomologEntrez[toFix & anno$gene=='apoc2'] = '344';
anno$hsaHomologGene[toFix & anno$gene=='asip2'] = 'ENSG00000101440';
anno$hsaHomologEntrez[toFix & anno$gene=='asip2'] = '434';
anno$hsaHomologEntrez[toFix & anno$gene=='crocc2'] = '728763';
anno$hsaHomologEntrez[toFix & anno$gene=='fignl2'] = '401720';
anno$hsaHomologEntrez[toFix & anno$gene=='grk1'] = '6011';
anno$hsaHomologEntrez[toFix & anno$gene=='inafm2'] = '100505573';
anno$hsaHomologEntrez[toFix & anno$gene=='kif28p'] = '100130097';
anno$hsaHomologEntrez[toFix & anno$gene=='muc19'] = '283463';
anno$hsaHomologEntrez[toFix & anno$gene=='ofcc1'] = '266553';
anno$hsaHomologEntrez[toFix & anno$gene=='pgbd5'] = '79605';
anno$hsaHomologEntrez[toFix & anno$gene=='rnpepl1'] = '57140';
anno$hsaHomologEntrez[toFix & anno$gene=='shank3'] = '85358';
anno$hsaHomologEntrez[toFix & anno$gene=='spata1'] = '100505741';
anno$hsaHomologEntrez[toFix & anno$gene=='tstd3'] = '100130890';
anno$hsaHomologEntrez[toFix & anno$gene=='umad1'] = '729852';

anno[which(toFix), ];

table(substr(anno$gene[is.na(anno$hsaHomologEntrez)], 1, 3))
# LOC 
# 3355

### try more mapping for zfid formatted oBlastxHit gene symbols
speciesVec = unlist(list(DAR='drerio_gene_ensembl',
                         ORL='olatipes_gene_ensembl',
                         ONI='oniloticus_gene_ensembl',
                         TRU='trubripes_gene_ensembl',
                         GAC='gaculeatus_gene_ensembl'));
zfid = .getHomologsFromDARFormatIDsMultipleSpecies(anno, speciesVec);

anno = .addDARFormatHomologs(anno, zfid);
table(substr(anno$gene[is.na(anno$hsaHomologEntrez)], 1, 3))
# LOC 
# 3331





write.table(recipPass, file='passedReciprocalBlastp_wHomologs.tsv', sep='\t', quote=F, col.names=T, row.names=F);
write.table(recipPass2, file='passedReciprocalBlastp_wHomologs_v2.tsv', sep='\t', quote=F, col.names=T, row.names=F);
write.table(recipFail, file='failedReciprocalBlastp.tsv', sep='\t', quote=F, col.names=T, row.names=F);
write.table(recipFail2, file='failedReciprocalBlastpSymbols_wHomologs.tsv', sep='\t', quote=F, col.names=T, row.names=F);
write.table(anno, file='passedReciprocalBlastp_wHomologs_v2_plusHomologsForFailedWithSymbols.tsv', sep='\t', quote=F, col.names=T, row.names=F);


##############################################################
#### functions
#############################################################

.addRecipBlastpToBlastxLookup = function (xLookup, rx, hsENSdf, verbose=1000) {
  blastQueryCol = 2;
  blastHitCol = 3;
  NAvec = rep(NA, nrow(xLookup));
  lookup = as.data.frame(cbind(xLookup, rBlastpHit=NAvec, recipPass=NAvec));
  # skip rows where the blastx hit didn't hit a burtoni protein in the reciprocal blastp
  toSkip = which(!(xLookup$oBlastxHit %in% rx[, blastQueryCol]));
  goodRows = (1:nrow(xLookup))[-toSkip];
  for (tr in goodRows) {
    #cat(tr,'...');
    if (tr %% verbose == 0) { cat(tr,'...') }
    xQuery = xLookup$transcript_id[tr];
    xQueryProtein = xLookup$protein_id[tr];
    if (is.na(xQueryProtein)) { next }
    xHit = xLookup$oBlastxHit[tr];
    pHit = rx[rx[,blastQueryCol]==xHit, blastHitCol];
    lookup$rBlastpHit[tr] = pHit;
    lookup$recipPass[tr] = xQueryProtein == pHit;
  }
  return(lookup);
}

.processLookupToGeneAnnotations = function (lookupGeneSplit, hsENSdf) {
  outdf = as.data.frame(matrix(nrow=length(lookupGeneSplit), ncol=ncol(lookupGeneSplit[[1]])));
  names(outdf) = names(lookupGeneSplit[[1]]);
  
  fails = c();
  passType = c();
  failType = c();
  
  for (g in 1:length(lookupGeneSplit)) {
    if (g %% 1000 == 0) { cat(g,'...') }
    #cat(g,'...')
    this_g = lookupGeneSplit[[g]];
    
    # check if only a single transcript
    #  this should be something like 16180 genes out of 25566
    #   > length(lookupGeneSplit)
    #   [1] 25566
    #   > sum(sapply(lookupGeneSplit, nrow) == 1)
    #   [1] 16180
    if (nrow(this_g) == 1) {
      # check if transcript was non-coding to begin with
      if (is.na(this_g$protein_id)) {
        fails = c(fails, g);
        failType = c(failType, 'single_ncRNA');
        passType = c(passType, NA);
        next;
      }
      # check if transcript's blastx hit didn't hit anything in reciprocal blastp
      if (is.na(this_g$rBlastpHit)) {
        fails = c(fails, g);
        failType = c(failType, 'single_no_blastp_hit');
        passType = c(passType, NA);
        next;
      }
      # check if transcript's blastx hit failed reciprocal blastp
      if (!this_g$recipPass) {
        fails = c(fails, g);
        failType = c(failType, 'single_fail_blastp');
        passType = c(passType, NA);
        next;
        # otherwise it passed
      } else {
        outdf[g, ] = this_g;
        passType = c(passType, 'single');
        next;
      }
    }
    
    # check for all failed reciprocal blastp
    anyRecipPass = any(this_g$recipPass, na.rm=T);
    if (!anyRecipPass) {
      fails = c(fails, g);
      failType = c(failType, rep('multiple_fail_blastp', nrow(this_g)));
      passType = c(passType, NA);
      next;
    }
    
    # remove any non-coding transcripts before continuing
    this_g = this_g[!is.na(this_g$protein_id), ];
    
    # check if all blastx or blastp hits are the same and at least one reciprocal blastp passed
    blastxHitsAllSame = length(unique(this_g$oBlastxHit)) == 1;
    blastpHitsAllSame = length(unique(this_g$rBlastpHit)) == 1;
    if ((blastxHitsAllSame | blastpHitsAllSame) & anyRecipPass) {
      outdf[g, ] = this_g[which(this_g$recipPass), ][1,];       # keep all instead of just first pass?
      passType = c(passType, 'all_same_blastx');
      next;
    } 
    
    # get genes for the different blastx hits
    blastxHitsGeneMap = .getMultiProteinGenes(this_g, hsENSdf);
    geneMaps = .process_getMultiProteinGenesOutput(blastxHitsGeneMap);
    
    # check genes the blastx hits map to
    #  are all of either ensembl ids, external symbol, or human homolog the same?
    if (geneMaps$checks$ens | geneMaps$checks$ext | geneMaps$checks$hsa) {
      outdf[g, ] = this_g[which(this_g$recipPass), ][1,];
      passType = c(passType, 'all_same_mapped_gene');
      next;
    }
    
    this_g_pass = this_g[this_g$recipPass, ];
    #  are any of the external gene symbols the same as the burtoni gene symbol
    #   and if so did the corresponding blastx hit pass reciprocal blastp?
    matchSym = which(geneMaps$gene$ext == names(lookupGeneSplit)[g]);
    if (length(matchSym) > 0) {
      blastx_proteins = geneMaps$protein[matchSym];
      oBlastxHitsStripped = substr(this_g_pass$oBlastxHit, 1, (nchar(this_g_pass$oBlastxHit)-2));
      outdf[g, ] = this_g_pass[match(blastx_proteins, oBlastxHitsStripped), ][1, ];  # use %in% ?
      passType = c(passType, 'same_symbol');
      next;
    }
    # what if last character of external symbol is stripped?
    ext_strip = substr(geneMaps$gene$ext, 1, (nchar(geneMaps$gene$ext)-1))
    matchSym = which(ext_strip == names(lookupGeneSplit)[g]);
    if (length(matchSym) > 0) {
      blastx_proteins = geneMaps$protein[matchSym];
      oBlastxHitsStripped = substr(this_g_pass$oBlastxHit, 1, (nchar(this_g_pass$oBlastxHit)-2));
      outdf[g, ] = this_g_pass[match(blastx_proteins, oBlastxHitsStripped), ][1, ];  # use %in% ?
      passType = c(passType, 'same_symbol_strip1');
      next;
    }
    
    # check if any single blastx hit was > .6 of all hits
    #  .6 forces this to only apply to genes with >2 transcripts
    blastxHitsCounts = table(this_g$oBlastxHit);
    hitEnough = which(blastxHitsCounts > (nrow(this_g)*.6));
    if (length(hitEnough) > 0) {
      outdf[g, ] = this_g_pass[match(names(hitEnough), this_g_pass$oBlastxHit), ][1, ];
      passType = c(passType, 'twothirds_same_blastx');
      next;
    }
    
    passType = c(passType, NA);

  }
  outdf = as.data.frame(cbind(outdf, passType=passType));
  allNAs = which(is.na(outdf[,1]));
  leftNAs = allNAs[!(allNAs %in% fails)];
  outdf = outdf[-fails, ];
  faildf = do.call(rbind, lookupGeneSplit[fails]);
  faildf = as.data.frame(cbind(faildf, failType=failType));
  
  return(list(pass=outdf, fail=faildf, leftNAs=lookupGeneSplit[leftNAs]));
}

# hsENSdf should be list of data frames
#  i.e. run .convertMultiSpeciesENStoHuman with spSplit=NULL
#       or run lapply(hsENS, function(f) do.call(rbind, f))
.getMultiProteinGenes = function (this_g, hsENSdf) {
  blastxHits = unique(this_g$oBlastxHit);
  blastxHits = substr(blastxHits, 1, nchar(blastxHits)-2);
  blastxHitsSpecies = unique(substr(blastxHits, 4, 6));
  hsg = list();
  for (sp in blastxHitsSpecies) {
    this_ens = hsENSdf[[match(sp, names(hsENSdf))]];
    this_hits = grep(sp, blastxHits, value=T);
    hsg[[sp]] = this_ens[this_ens$ensembl_peptide_id %in% this_hits, ];
  }
  return(hsg);
}

# assumes anyRecipPass==T
.process_getMultiProteinGenesOutput = function (hsg) {
  ensgene = unlist(lapply(hsg, function(f) f$ensembl_gene_id));
  extgene = strsplit(unlist(lapply(hsg, function(f) f$external_gene_name)), ' (', fixed=T);
  extgene = tolower(sapply(extgene, function(f) f[1]));
  hsagene = unlist(lapply(hsg, function(f) f$hsapiens_homolog_ensembl_gene));
  ensCheck = length(unique(ensgene)) == 1;
  extCheck = length(unique(extgene)) == 1;
  hsaCheck = length(unique(hsagene)) == 1;
  outgenes = list(ens=ensgene, ext=extgene, hsa=hsagene);
  outchecks = list(ens=ensCheck, ext=extCheck, hsa=hsaCheck);
  proteins = unlist(lapply(hsg, function(f) f$ensembl_peptide_id));
  return(list(gene=outgenes, checks=outchecks, protein=proteins));
}

.convertMultiSpeciesENStoHuman = function (multiSpeciesENSvec, speciesMarts, attributes, filters, spSplit=NULL) {
  speciesTriplets = unique(substr(multiSpeciesENSvec, 4, 6));
  if (!all(speciesTriplets %in% names(speciesMarts))) { stop('All species in ENS ids not represented in speciesMarts') }
  if (any(grepl('.1', multiSpeciesENSvec, fixed=T))) {
    cat(paste0('Stripping final 2 characters from ', deparse(substitute(multiSpeciesENSvec)), '...\n'));
    multiSpeciesENSvec = substr(multiSpeciesENSvec, 1, nchar(multiSpeciesENSvec)-2);
  }
  hsENS = list();
  cat('Working on... ');
  for (species in speciesTriplets) {
    cat(species, '...');
    spGenesInds = which(grepl(species, multiSpeciesENSvec));
    spMart = speciesMarts[[match(species, names(speciesMarts))]];
    hsENS[[species]] = getBM(attributes=attributes, filters=filters, values=multiSpeciesENSvec[spGenesInds], mart=spMart);
    if (is.numeric(spSplit)) {
      hsENS[[species]] = split(hsENS[[species]], hsENS[[species]][,spSplit]);
    }
  }
  return(hsENS);
}

.convertHsEnsemblToEntrez = function (hs_ensemblVec, host='www.ensembl.org', bmart='ENSEMBL_MART_ENSEMBL') {
  mart = useMart(biomart=bmart, dataset='hsapiens_gene_ensembl', host=host);
  attributes = c('ensembl_gene_id', 'entrezgene');
  filters = 'ensembl_gene_id';
  res = getBM(attributes=attributes, filters=filters, values=hs_ensemblVec, mart);
  return(res)
}

.biomartHs = function (values, host='www.ensembl.org', bmart='ENSEMBL_MART_ENSEMBL', attributes=c('ensembl_gene_id','entrezgene'), filters='ensembl_gene_id') {
  mart = useMart(biomart=bmart, dataset='hsapiens_gene_ensembl', host=host);
  res = getBM(attributes=attributes, filters=filters, values=values, mart);
  return(res)
}

.addGenesToRecipPass = function (recipPass, hsENSdf, hsEntrez) {
  blastxHits = recipPass$oBlastxHit;
  geneMaps = as.data.frame(do.call(rbind, hsENSdf));
  fillNAs = rep(NA, nrow(recipPass));
  out = as.data.frame(cbind(recipPass, 
                            oBlastxHitGene=fillNAs, 
                            oBlastxHitDescription=fillNAs,
                            oBlastxHitSym=fillNAs,
                            hsaHomologGene=fillNAs,
                            hsaHomologEntrez=fillNAs))
  for (row in 1:nrow(recipPass)) {
    if (row %% 1000 == 0) {cat(row,'...')}
    tmp = subset(geneMaps, ensembl_peptide_id==substr(blastxHits[row], 1, (nchar(blastxHits[row])-2)));
    if(nrow(tmp) == 0) {
      print(recipPass[row, ]); stop('blastxHit missing from hsENSdf');
    }
    hsa_entrez = hsEntrez[which(hsEntrez$ensembl_gene_id %in% tmp$hsapiens_homolog_ensembl_gene), ];
    if (nrow(tmp) > 1) {
      hsaENS = paste(tmp$hsapiens_homolog_ensembl_gene, collapse=',');
    } else {
      hsaENS = tmp$hsapiens_homolog_ensembl_gene;
    }
    if (nrow(hsa_entrez) > 0) {
      hsa_entrez = paste(as.vector(na.omit(hsa_entrez$entrezgene)), collapse=',');
    } else {
      hsa_entrez = NA;
    }
    out$oBlastxHitGene[row] = unique(tmp$ensembl_gene_id);
    out$oBlastxHitDescription[row] = unique(tmp$description);
    out$oBlastxHitSym[row] = unique(tmp$external_gene_name);
    out$hsaHomologGene[row] = hsaENS;
    out$hsaHomologEntrez[row] = hsa_entrez;
  }
  for (i in 1:ncol(out)) {
    missingvals = out[, i] == '';
    if (any(missingvals, na.rm=T)) {
      out[which(missingvals), i] = NA;
      }    
  }
  return(out);
}

.addMissingHsHomologsToRecipPassSymbols = function (recipPass, strip=NULL) {
  
  # check that every row has valid blastx hit
  goodBlastxHits = grepl('^ENS', recipPass$oBlastxHitGene);
  if (sum(goodBlastxHits) != nrow(recipPass)) { stop() }
  
  # get rows with homolog info
  goodHsaENS = grepl('^ENS', recipPass$hsaHomologGene);
  goodHsaEntrez = grepl('^[0-9]|[0-9]$', recipPass$hsaHomologEntrez);
  
  # get gene symbols from burtoni and blastx hits
  abSyms = grepl('^[a-z]', recipPass$gene);
  oxSyms = grepl('^[a-z]', tolower(recipPass$oBlastxHitSym));
  symRows = abSyms | oxSyms;
  
  # get rows that have a gene symbol but are missing either type of homolog id
  #  likely some will have ensembl but none will have entrez ids
  symRowsToFix = symRows & (!goodHsaENS | !goodHsaEntrez);
  
  syms = unique(c(recipPass$gene[symRowsToFix], recipPass$oBlastxHitSym[symRowsToFix]));
  syms = unique(grep(':|^LOC', syms, invert=T, value=T));
  syms = unique(sapply(strsplit(syms, ' (', fixed=T), function(f) f[1]));
  if (is.numeric(strip)) {
    symsStripped = grep('[0-9][a-z]$', syms, value=T);
    syms = unique(c(syms, substr(symsStripped, 1, nchar(symsStripped)-strip)));
  }
  maps = .biomartHs(values=syms, 
                    attributes=c('ensembl_gene_id','entrezgene','hgnc_symbol'), 
                    filters='hgnc_symbol');
  
  for (r in which(symRowsToFix)) {
    this_sym = unique(c(recipPass$gene[r], recipPass$oBlastxHitSym[r]));
    this_sym = sapply(strsplit(this_sym, ' (', fixed=T), function(f) f[1]);
    this_sym = this_sym[!is.na(this_sym)];
    this_map = subset(maps, tolower(hgnc_symbol) %in% this_sym);
    map_g = NA; map_e = NA;
    if (nrow(this_map) > 0) {
      map_g0 = as.vector(na.omit(unique(this_map$ensembl_gene_id)));
      map_e0 = as.vector(na.omit(unique(this_map$entrezgene)));
      if (length(map_g0) >= 1) { map_g = paste(map_g0, collapse=',') }
      if (length(map_e0) >= 1) { map_e = paste(map_e0, collapse=',') }
    } 
    recipPass$hsaHomologGene[r] = map_g;
    recipPass$hsaHomologEntrez[r] = map_e;
  }
  return(recipPass);
}

.addGenesToFailSymbols = function (recipFail0) {
  tmp = subset(recipFail0, !grepl('^LOC', gene) & !duplicated(gene));
  syms = unique(grep(':|^LOC', unique(tmp$gene), invert=T, value=T));
  failSymsMap = .biomartHs(values=syms, 
                           attributes=c('ensembl_gene_id','entrezgene','hgnc_symbol','description'), 
                           filters='hgnc_symbol');
  fillNAs = rep(NA, nrow(tmp));
  tmp = as.data.frame(cbind(tmp, 
                            oBlastxHitGene=fillNAs, 
                            oBlastxHitDescription=fillNAs,
                            oBlastxHitSym=fillNAs,
                            hsaHomologGene=fillNAs,
                            hsaHomologEntrez=fillNAs));
  for (s in 1:nrow(tmp)) {
    map_g = NA; map_e = NA;
    this_sym = tmp$gene[s];
    this_map = subset(failSymsMap, tolower(hgnc_symbol)==this_sym);
    map_g = unique(this_map$ensembl_gene_id);
    map_e = unique(this_map$entrezgene);
    map_d = unique(this_map$description);
    if (length(map_g) > 1) { map_g = paste(map_g, collapse=',') }
    if (length(map_e) > 1) { map_e = paste(map_e, collapse=',') }
    if (length(map_d) > 1) { print(this_map); print(tmp[s,]); stop('duplicate description') }
    tmp$hsaHomologGene[s] = map_g;
    tmp$hsaHomologEntrez[s] = map_e;
    tmp$oBlastxHitDescription[s] = map_d;
  }
  return(tmp);
}

.convertIDsWithBiomaRt = function (biomart='ENSEMBL_MART_ENSEMBL', host='www.ensembl.org', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', ids) {
  mart = useMart(biomart=biomart, dataset=dataset, host=host);
  return(getBM(attributes=c(fromID, toID), filters=fromID, values=ids, mart=mart));
}

.getHomologsFromDARFormatIDs = function (anno, dataset='drerio_gene_ensembl') {
  cat('Getting ids... ');
  darIDs = subset(anno, grepl(':',oBlastxHitSym) & is.na(hsaHomologEntrez))$oBlastxHitSym;
  darIDs = unique(sapply(strsplit(darIDs, ' (', fixed=T), function(f) f[1]));
  cat('hsa homologs... ');
  zf0 = .convertIDsWithBiomaRt(dataset=dataset, 
                               fromID='external_gene_name', 
                               toID=c('external_gene_name','ensembl_gene_id','hsapiens_homolog_ensembl_gene'), 
                               ids=darIDs);
  zf = subset(zf0, hsapiens_homolog_ensembl_gene!='');
  cat('converting to entrez... ');
  zfentrez = .convertHsEnsemblToEntrez(zf$hsapiens_homolog_ensembl_gene);
  out = list(fish=zf, hsa=zfentrez);
  return(out);
}

.getHomologsFromDARFormatIDsMultipleSpecies = function (anno, speciesVec) {
  out = list();
  for (s in 1:length(speciesVec)) {
    cat('[', names(speciesVec)[s], '] ...');
    out[[names(speciesVec)[s]]] = .getHomologsFromDARFormatIDs(anno=anno, dataset=speciesVec[s]);
  }
  return(out);
}

.addDARFormatHomologs = function (anno, zfid) {
  rows = which(is.na(anno$hsaHomologEntrez) & grepl(':',anno$oBlastxHitSym));
  for (r in rows) {
    this_id = strsplit(anno$oBlastxHitSym[r], ' (', fixed=T)[[1]][1];
    this_species = substr(anno$oBlastxHitGene[r], 4, 6);
    map = zfid[[match(this_species, names(zfid))]];
    check = which(map$fish$external_gene_name == this_id);
    if (length(check) > 0) {
      if (!is.na(anno$hsaHomologGene[r])) { print(anno[r, ]); stop() }
      this_ens = unique(map$fish$hsapiens_homolog_ensembl_gene[check]);
      this_entrez = as.vector(na.omit(map$hsa$entrezgene[match(this_ens, map$hsa$ensembl_gene_id)]));
      anno$hsaHomologGene[r] = paste(this_ens, collapse=',');
      if (length(this_entrez) > 0) {
        anno$hsaHomologEntrez[r] = paste(this_entrez, collapse=',');
      } else {
        anno$hsaHomologEntrez[r] = NA;
      }
    } else {
      next;
    }
  }
  return(anno);
}