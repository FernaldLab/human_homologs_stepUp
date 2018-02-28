####################################################################################################
##### functions
####################################################################################################

####################################################################################################
#### extract desired attributes from gff metadata column
## Args
##  gff: data frame containing at least one valid row of a gff file
##  metaCol: number designating the column to be parsed, will always be 9 in well-formed gffs
##  delim: character delimiter to parse on, will always be ';' in well-formed gffs
##  pattern: name of attribute to be extracted, default is 'Dbxref=GeneID:'
## Output
##  character vector containing the parsed attribute values
## Comment
##  should add options to turn off checks since they add significant running time to parent functions

.parseGffMetaCol = function(gff, metaCol=9, delim=';', pattern='Dbxref=GeneID:') {
  metaSplit = strsplit(gff[, metaCol], delim);
  check = sapply(metaSplit, function(f) any(grepl(pattern, f)))
  if (!all(check) | sum(check)!=nrow(gff)) {
    stop('all lines in gff do not contain arg pattern') 
  }
  res = grep(pattern, unlist(metaSplit), value=TRUE);
  return(gsub(pattern, '', res));
}

####################################################################################################
#### make dataframe containing all geneID-geneSym-rna_num-transcript_id-protein_id mappings for a gene
####  meant to be used in the function .buildGffLookup()
####  requires .parseGffMetaCol()
## Args
##  geneDFgff: data frame containing all relevant lines for a particular gene
##             first line must be a gene feature and there must be some type of rna feature line
##             exon, but not CDS, features are required
##  metaCol: number designating the column to be parsed, will always be 9 in well-formed gffs
##  featureCol: number designating the column that contains feature type, will always be 3 in well-formed gffs
##  delim: character delimiter to parse metadata column on, will always be ';' in well-formed gffs
## Output
##  data frame with 5 columns: Dbxref, gene, rna_num, transcript_id, protein_id
##                  1 row per rna
## Comment
##  this is slow but should be guaranteed correct mappings
##  lookupColNames is basically pointless since these names are hardcoded into multiple parts of the function

.build_mRNAdfForGene = function (geneDFgff, metaCol=9, featureCol=3, delim=';') {
  
  # set some parameters for metadata parsing and output format
  #  probably won't change but wanted them here all in one place
  geneIDPattern = 'Dbxref=GeneID:';
  geneSymbolPattern = 'gene=';
  parentPattern = 'Parent=';
  rnaNamePattern = 'Name=';
  rnaNumPattern = 'ID=';
  transcriptPattern = 'transcript_id=';
  proteinPattern = 'protein_id=';
  lookupColNames = c('Dbxref','gene','rna_num','transcript_id','protein_id');
  
  # check that first row is a 'gene' feature
  if (geneDFgff[1, featureCol] != 'gene') {
    print(unique(geneDFgff[ ,featureCol])); 
    stop('First row should be gene feature');
  }
  
  # get gene ID and symbol, then remove gene feature row
  geneID = .parseGffMetaCol(gff=geneDFgff[1,], metaCol=metaCol, delim=delim, pattern=geneIDPattern);
  geneSymbol = .parseGffMetaCol(gff=geneDFgff[1,], metaCol=metaCol, delim=delim, pattern=geneSymbolPattern);
  geneDFgff = geneDFgff[-1, ];
  
  # split on Parent= field 
  #  first data frame in resulting list will hold the rna feature lines
  #   could be mRNA, ncRNA, transcript
  #  rest of data frames will hold exons and CDSs for rnas, named by rna num
  #  this means exon and CDS features for the same transcript in the same data frames
  #   which should guarantee correct mappings
  parents = .parseGffMetaCol(gff=geneDFgff, metaCol=metaCol, delim=delim, pattern=parentPattern);
  geneParentSplit = split(geneDFgff, parents);
  
  # map rna ids to nums
  #  df is a data frame containing rna feature lines
  #  rnas is a vector where values are rna names (e.g. XM_) and names are rna nums
  df = geneParentSplit[[grep('^gene', names(geneParentSplit))]];
  rnas = .parseGffMetaCol(gff=df, metaCol=metaCol, delim=delim, pattern=rnaNamePattern);
  names(rnas) = .parseGffMetaCol(gff=df, metaCol=metaCol, delim=delim, pattern=rnaNumPattern);
  
  # make output table and add gene info
  out = as.data.frame(matrix(nrow=length(rnas), ncol=length(lookupColNames)));
  names(out) = lookupColNames;
  #  if lookupColNames ever changes this will also need to
  out$Dbxref = rep(geneID, nrow(out));
  out$gene = rep(geneSymbol, nrow(out));
  
  # loop through rnas for this gene
  for (num in 1:length(rnas)) {
    # get current id from df
    id = rnas[num];
    df = geneParentSplit[[match(names(id), names(geneParentSplit))]];
    
    # get rna and protein id from current df
    exrows = df[, featureCol] == 'exon';
    transcript_id = unique(.parseGffMetaCol(gff=df[exrows, ], metaCol=metaCol, delim=delim, pattern=transcriptPattern));
    
    # check if transcript codes for a protein before trying to get protein id
    #  ncRNAs and 'transcript' features will have NA in the protein column of output
    cdsrows = df[, featureCol] == 'CDS';
    if (any(cdsrows)) {
      protein_id = unique(.parseGffMetaCol(gff=df[cdsrows, ], metaCol=metaCol, delim=delim, pattern=proteinPattern));
    } else {
      protein_id = NA;
    }
    
    # make sure all ids are as expected
    if ((transcript_id != id) | (length(transcript_id) > 1) | (length(protein_id) > 1)) {
      print(geneParentSplit); stop('Problem with transcript_id or protein_id');
    }
    
    # add output for this rna
    #  if lookupColNames ever changes this will also need to
    out$rna_num[num] = names(rnas)[num];
    out$transcript_id[num] = transcript_id;
    out$protein_id[num] = protein_id;
  }
  # final output table
  return(out);
}

####################################################################################################
#### make dataframe containing all geneID-geneSym-rna_num-transcript_id-protein_id mappings from the gff
## Args
##  gff: dataframe containing a well-formed gff file (ok to have extra columns)
##       column 3 (feature type) is assumed to contain:
##         C_gene_segment, CDS, exon, gene, mRNA, ncRNA, transcript, V_gene_segment, and possibly region
##  featureCol: number designating the column that contains feature type, will always be 3 in well-formed gffs
##  sourceCol: number designating the column with annotation source, will always be 2 in well-formed gffs
##  metaCol: number designating the column to be parsed, will always be 9 in well-formed gffs
##  metaDelim: character delimiter to parse on, will always be ';' in well-formed gffs
##  verbose: number indicating how often to print a line # update to the console
##  numToTest: number indicating how many genes to process, useful for benchmarking and testing
##  write: logical indicating whether the final lookup table should be written to 'working_dir/gff_lookup.tsv'
##  returnReducedGff: logical indicating whether to return the gff 
##
## Output
##  a 2 element list, containing $lookup and $gffReduced
##   if returnReducedGff=T, $gffReduced is the input gff reduced to lines with 'GeneID' in the metadata column and
##    the following features removed: 
##     -tRNAs (based on 'tRNAscan-SE' in source column)
##     -V_gene_segments and C_gene_segments (based on 'V_gene_segment|C_gene_segment' in feature type column)
##     -pseudo genes (based on lack of 'mRNA|ncRNA|transcript' in feature type column)
##    not just the exact lines with/without these strings are removed, but every associated line as well (e.g. exons)
##    The reduced gff is sorted by gene name, may need to coordinate sort for use in things like IGV
##   if returnReducedGff=F, $gffReduced is NULL
##
## Comment
##  Again, lookupColNames is basically pointless since these names are hardcoded throughout the function
##
##  returnReducedGff=T is extremely slow, it's running do.call(rbind, gffSplitGeneSym) since V/C_gene_segments 
##   and pseudo genes have to be removed after making the gffSplitGeneSym list
##   But, do.call is ~4x faster than a for loop on my computer
##   Alternatively could return the version of the gff that's input to split, but it will still contain 
##    V/C_gene_segments and pseudo genes
##   Or could just return gffSplitGeneSym
##
##  if returnReducedGff=F, function takes ~10min on 8-core 4GHz i7 iMac, ~30min if returnReducedGff=T
##
##  Could add option for user to provide gffSplitGeneSym, to skip splitting step
##  Could add options for writing files/filenames

.buildGffLookup = function (gff, featureCol=3, sourceCol=2, metaCol=9, metaDelim=';', verbose=500, numToTest=NULL, write=F, returnReducedGff=T) {
  
  tRNAsourcePattern = 'tRNAscan-SE';
  geneFilterPattern = 'GeneID';
  geneSymPattern = 'gene=';
  v_c_segPattern = 'V_gene_segment|C_gene_segment';
  noPseudoPattern = 'mRNA|ncRNA|transcript';
  lookupColNames = c('Dbxref','gene','rna_num','transcript_id','protein_id');
  
  # check that input gff is a data frame
  if (!is.data.frame(gff)) { stop('Input gff must be a data frame') }
  
  # get rid of region feature lines if needed
  if (any(gff[, featureCol] == 'region')) {
    cat('Removing region feature lines... ');
    gff = gff[gff[, featureCol] != 'region', ];
  }
  
  # only keep lines with "GeneID" (geneFilterPattern) in the metadata column
  geneIDrows = grep(geneFilterPattern, gff[, metaCol]);
  gff = gff[geneIDrows, ];
  
  # get rid of tRNAs
  tRNArows = gff[, sourceCol] == tRNAsourcePattern;
  if (any(tRNArows)) { gff = gff[!tRNArows, ] }
  
  # table(gff[, featureCol])
  # C_gene_segment  CDS      exon     gene    mRNA    ncRNA   transcript   V_gene_segment 
  # 17              562346   606181   26104   44394   2492    621          72 
  
  # make list by splitting on gene symbols (geneSymPattern, from gene=)
  #  would get same results using geneIDs (from Dbxref=GeneID:) but with one extra step needed
  #   (to strip rna and protein ids from Dbxref field)
  #   gffGeneIDs = .parseGffMetaCol(gff=gff, metaCol=metaCol, delim=metaDelim, pattern='Dbxref=GeneID:');
  #   gffGeneIDs = grep('^Genbank', unlist(strsplit(gffGeneIDs,',')), value=T, invert=T)
  if (is.numeric(verbose)) { cat('Getting genes... ') }
  gffGeneSyms = .parseGffMetaCol(gff=gff, metaCol=metaCol, delim=metaDelim, pattern=geneSymPattern);
  if (is.numeric(verbose)) { cat('Formatting gff... ') } 
  gffSplitGeneSym = split(gff, gffGeneSyms);
  
  # subset the gene list if needed
  if (is.numeric(numToTest)) { 
    cat(paste0('Processing first ', numToTest, ' genes... '));
    gffSplitGeneSym = gffSplitGeneSym[1:numToTest];
  }
  
  # # only types of 'Parent=' are gene, rna, id
  # #  all the Parent=id rows are either V_gene_segment or C_gene_segment features
  # geneFeatureRows = which(gff[, featureCol] == 'gene');
  # xp = .parseGffMetaCol(gff[!geneFeatureRows, ], metaCol=metaCol, delim=metaDelim, pattern='Parent=')
  # table(substr(xp,1,2))
  # ge      id      rn 
  # 48719   222 1167182 
  
  # remove V_gene_segment and C_gene_segment features
  indsVC = sapply(gffSplitGeneSym, function(f) any(grepl(v_c_segPattern, f[, featureCol])));
  if (any(indsVC)) { 
    gffSplitGeneSym = gffSplitGeneSym[-which(indsVC)];
  }
  
  # remove pseudogenes by only keeping genes with an associated rna
  indsNoPseudo = sapply(gffSplitGeneSym, function(f) any(grepl(noPseudoPattern, f[, featureCol])));
  if (any(indsNoPseudo)) {
    gffSplitGeneSym = gffSplitGeneSym[which(indsNoPseudo)];
  }
  
  # intialize lookup table output
  lookup = as.data.frame(matrix(nrow=1, ncol=length(lookupColNames)));
  names(lookup) = lookupColNames;
  
  # iterate through genes to build lookup table
  cat('Building lookup table... ');
  if (is.numeric(verbose)) {
    for (g in 1:length(gffSplitGeneSym)) {
      if (g %% verbose == 0) { cat(g, '...') }
      lookup = rbind(lookup, .build_mRNAdfForGene(gffSplitGeneSym[[g]], metaCol, featureCol, metaDelim));
    }
  } else {
    for (g in 1:length(gffSplitGeneSym)) {
      lookup = rbind(lookup, .build_mRNAdfForGene(gffSplitGeneSym[[g]], metaCol, featureCol, metaDelim));
    }
  }
  
  # remove first row since only NAs from initialization
  lookup = lookup[-1, ];
  
  # write lookup to file if needed
  if (write) {
    cat('Writing output file gff_lookup.tsv to working directory... ');
    write.table(lookup, file='gff_lookup.tsv', col.names=T, quote=F, row.names=F, sep='\t');
  }
  
  # rebuild the gff if needed, will take probably at least 20min
  if (returnReducedGff) {
    cat('Building reduced gff... ');
    gffReduced = do.call(rbind, gffSplitGeneSym);
  } else {
    gffReduced = NULL;
  }
  return(list(lookup=lookup, gffReduced=gffReduced));
}

####################################################################################################
#### make a version of the gff with gene feature lines only
## Args
##  gff: dataframe containing a well-formed gff file (ok to have extra columns)
##       column 3 (feature type) is assumed to contain at least some gene lines
##       typically assume this will be $gffReduced from .buildGffLookup
##  featureCol: number designating the column that contains feature type, will always be 3 in well-formed gffs
##  metaCol: number designating the column to be parsed, will always be 9 in well-formed gffs
##  delim: character delimiter to parse on, will always be ';' in well-formed gffs
##  pattern: name of attribute to be extracted, default is 'gene='
##  addTranscripts: logical indicating whether to add additional column to the output containing the # of transcripts/gene
##                  args metaCol, delim, pattern, and verbose only needed if addTranscripts=T
##  verbose: number indicating how often to print a line # update to the console
## Output
##  dataframe containing only the lines from the input that have 'gene' in the feature column
##   has extra column named 'numTranscripts' if addTranscripts=T
## Comment
##  addTranscripts=T is much slower since otherwise it's just a simple subset on the feature column
##   this should probably just be a separate function

.makeGeneOnlyGff = function (gff, featureCol=3, metaCol=9, delim=';', pattern='gene=', addTranscripts=F, verbose=500) {
  if (addTranscripts) {
    # get gene names to split the gff into a list
    #  there are faster ways to get # of transcripts/gene but this is most trustworthy for me right now
    cat('Getting gene names... ');
    genes = .parseGffMetaCol(gff, metaCol=metaCol, delim=delim, pattern=pattern);
    cat(paste0('Splitting ', deparse(substitute(gff)), ' into gene list... '));
    gffSplit = split(gff, genes);
    
    # add column to the output
    gffGenes = as.data.frame(matrix(nrow=length(gffSplit), ncol=(ncol(gff)+1)));
    names(gffGenes) = c(names(gff), 'numTranscripts');
    
    # iterate through genes to count transcripts
    for (g in 1:length(gffSplit)) {
      if (g %% verbose == 0) { cat(g, '...') }
      this_g = gffSplit[[g]];
      ex = subset(this_g, this_g[,featureCol]=='exon');
      tr = .parseGffMetaCol(ex, metaCol=metaCol, delim=delim, pattern='transcript_id=');
      ntr = length(unique(tr));
      geneRow = subset(this_g, this_g[,featureCol]=='gene');
      gffGenes[g, ] = c(geneRow, ntr);
    }
  } else {
    # if addTranscripts=F just do a quick subset
    cat(paste0('Subsetting ', deparse(substitute(gff)), ' to rows where column ', featureCol, '=="gene"'));
    gffGenes = subset(gff, gff[,featureCol]=='gene');
  }
  return(gffGenes);
}