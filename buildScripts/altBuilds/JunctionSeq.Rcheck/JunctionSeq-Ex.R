pkgname <- "JunctionSeq"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('JunctionSeq')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("JunctionSeqCountSet-class")
### * JunctionSeqCountSet-class

flush(stderr()); flush(stdout())

### Name: JunctionSeqCountSet-class
### Title: Class '"JunctionSeqCountSet"'
### Aliases: JunctionSeqCountSet-class JunctionSeqCountSet
### Keywords: classes

### ** Examples

showClass("JunctionSeqCountSet")



cleanEx()
nameEx("buildAllPlots")
### * buildAllPlots

flush(stderr()); flush(stdout())

### Name: buildAllPlots
### Title: Create and save a full battery of JunctionSeq expression plots.
### Aliases: buildAllPlots

### ** Examples


data(exampleDataSet,package="JctSeqExData2");
buildAllPlots(jscs);

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Generate all plots and the html index
##D #   Save them as pngs to the current directory:
##D buildAllPlots(jscs);
##D 
## End(Not run)



cleanEx()
nameEx("buildAllPlotsForGene")
### * buildAllPlotsForGene

flush(stderr()); flush(stdout())

### Name: buildAllPlotsForGene
### Title: Create and save one or more JunctionSeq expression plots.
### Aliases: buildAllPlotsForGene

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
buildAllPlotsForGene(geneID = "ENSRNOG00000009281", jscs);

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Generate several related plots for the selected gene:
##D buildAllPlotsForGene(geneID = "ENSRNOG00000009281", jscs);
##D 
## End(Not run)



cleanEx()
nameEx("defaultColorList")
### * defaultColorList

flush(stderr()); flush(stdout())

### Name: defaultColorList
### Title: JunctionSeq Color Parameters
### Aliases: JUNCTIONSEQ.DEFAULT.COLOR.LIST defaultColorList
###   junctionSeqColors
### Keywords: datasets

### ** Examples

data(exampleDataSet,package="JctSeqExData2");

#Set a few alternative colors:
buildAllPlotsForGene(geneID = "ENSRNOG00000009281", jscs, 
                     outfile.prefix = "./oddColors.",
                     colorList = list(SIG.FEATURE.COLOR = "red",
                                      SIG.FEATURE.FILL.COLOR = "green",
                                      NOSIG.FEATURE.FILL.COLOR = "blue"
                                      ));



cleanEx()
nameEx("estimateEffectSizes")
### * estimateEffectSizes

flush(stderr()); flush(stdout())

### Name: estimateEffectSizes
### Title: Estimate Effect Sizes, parameter estimates, etc.
### Aliases: estimateEffectSizes

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
jscs <- estimateEffectSizes(jscs);

## Not run: 
##D #Full example (from scratch):
##D 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Advanced Analysis:
##D 
##D #Make a "design" dataframe:
##D design <- data.frame(condition = factor(decoder$group.ID));
##D #Read the QoRTs counts.
##D jscs = readJunctionSeqCounts(countfiles = countFiles,
##D            samplenames = decoder$sample.ID,
##D            design = design,
##D            flat.gff.file = gff.file
##D );
##D #Generate the size factors and load them into the JunctionSeqCountSet:
##D jscs <- estimateJunctionSeqSizeFactors(jscs);
##D #Estimate feature-specific dispersions:
##D jscs <- estimateJunctionSeqDispersions(jscs);
##D #Fit dispersion function and estimate MAP dispersion:
##D jscs <- fitJunctionSeqDispersionFunction(jscs);
##D #Test for differential usage:
##D jscs <- testForDiffUsage(jscs);
##D #Estimate effect sizes and expression estimates:
##D jscs <- estimateEffectSizes( jscs);
##D 
## End(Not run)



cleanEx()
nameEx("estimateJunctionSeqDispersions")
### * estimateJunctionSeqDispersions

flush(stderr()); flush(stdout())

### Name: estimateJunctionSeqDispersions
### Title: JunctionSeq Dispersion Estimation
### Aliases: estimateJunctionSeqDispersions

### ** Examples

## Not run: 
##D #Full example (from scratch):
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Advanced Analysis:
##D 
##D #Make a "design" dataframe:
##D design <- data.frame(condition = factor(decoder$group.ID));
##D #Read the QoRTs counts.
##D jscs = readJunctionSeqCounts(countfiles = countFiles,
##D            samplenames = decoder$sample.ID,
##D            design = design,
##D            flat.gff.file = gff.file
##D );
##D #Generate the size factors and load them into the JunctionSeqCountSet:
##D jscs <- estimateJunctionSeqSizeFactors(jscs);
##D #Estimate feature-specific dispersions:
##D jscs <- estimateJunctionSeqDispersions(jscs);
##D #Fit dispersion function and estimate MAP dispersion:
##D jscs <- fitJunctionSeqDispersionFunction(jscs);
##D #Test for differential usage:
##D jscs <- testForDiffUsage(jscs);
##D #Estimate effect sizes and expression estimates:
##D jscs <- estimateEffectSizes( jscs);
##D 
## End(Not run)



cleanEx()
nameEx("estimateJunctionSeqSizeFactors")
### * estimateJunctionSeqSizeFactors

flush(stderr()); flush(stdout())

### Name: estimateJunctionSeqSizeFactors
### Title: Estimate Size Factors
### Aliases: estimateJunctionSeqSizeFactors writeSizeFactors

### ** Examples


data(exampleDataSet,package="JctSeqExData2");
jscs <- estimateJunctionSeqSizeFactors(jscs);

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Advanced Analysis:
##D 
##D #Make a "design" dataframe:
##D design <- data.frame(condition = factor(decoder$group.ID));
##D #Read the QoRTs counts.
##D jscs = readJunctionSeqCounts(countfiles = countFiles,
##D            samplenames = decoder$sample.ID,
##D            design = design,
##D            flat.gff.file = gff.file
##D );
##D #Generate the size factors and load them into the JunctionSeqCountSet:
##D jscs <- estimateJunctionSeqSizeFactors(jscs);
##D #Estimate feature-specific dispersions:
##D jscs <- estimateJunctionSeqDispersions(jscs);
##D #Fit dispersion function and estimate MAP dispersion:
##D jscs <- fitJunctionSeqDispersionFunction(jscs);
##D #Test for differential usage:
##D jscs <- testForDiffUsage(jscs);
##D #Estimate effect sizes and expression estimates:
##D jscs <- estimateEffectSizes( jscs);
##D 
## End(Not run)



cleanEx()
nameEx("fitJunctionSeqDispersionFunction")
### * fitJunctionSeqDispersionFunction

flush(stderr()); flush(stdout())

### Name: fitJunctionSeqDispersionFunction
### Title: Fit Shared Dispersion Function
### Aliases: fitJunctionSeqDispersionFunction

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
jscs <- fitJunctionSeqDispersionFunction(jscs);

## Not run: 
##D #Full example (from scratch):
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Advanced Analysis:
##D 
##D #Make a "design" dataframe:
##D design <- data.frame(condition = factor(decoder$group.ID));
##D #Read the QoRTs counts.
##D jscs = readJunctionSeqCounts(countfiles = countFiles,
##D            samplenames = decoder$sample.ID,
##D            design = design,
##D            flat.gff.file = gff.file
##D );
##D #Generate the size factors and load them into the JunctionSeqCountSet:
##D jscs <- estimateJunctionSeqSizeFactors(jscs);
##D #Estimate feature-specific dispersions:
##D jscs <- estimateJunctionSeqDispersions(jscs);
##D #Fit dispersion function and estimate MAP dispersion:
##D jscs <- fitJunctionSeqDispersionFunction(jscs);
##D #Test for differential usage:
##D jscs <- testForDiffUsage(jscs);
##D #Estimate effect sizes and expression estimates:
##D jscs <- estimateEffectSizes( jscs);
##D 
## End(Not run)



cleanEx()
nameEx("plotDispEsts")
### * plotDispEsts

flush(stderr()); flush(stdout())

### Name: plotDispEsts
### Title: Plot Fitted and Test-wise Dispersion
### Aliases: plotDispEsts

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
plotDispEsts(jscs);

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ######################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Plot dispersions:
##D plotDispEsts(jscs);
##D 
## End(Not run)



cleanEx()
nameEx("plotJunctionSeqResultsForGene")
### * plotJunctionSeqResultsForGene

flush(stderr()); flush(stdout())

### Name: plotJunctionSeqResultsForGene
### Title: Generate a JunctionSeq expression plot.
### Aliases: plotJunctionSeqResultsForGene

### ** Examples

data(exampleDataSet,package="JctSeqExData2");

plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs);


## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ######################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Make an expression plot for a given gene:
##D plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs);
##D 
##D #Plot normalized read counts for a given gene:
##D plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
##D                     plot.type = "normCounts");
##D 
##D #Plot relative expression for a given gene:
##D plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
##D                     plot.type = "rExpr");
##D 
##D #Plot raw read counts for a given gene:
##D plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
##D                     plot.type = "rawCounts");
##D 
##D #Same thing, but with isoforms shown:
##D plotJunctionSeqResultsForGene(geneID = "ENSRNOG00000009281", jscs,
##D                     plot.type = "rawCounts",
##D                     displayTranscripts = TRUE);
##D 
## End(Not run)



cleanEx()
nameEx("plotMA")
### * plotMA

flush(stderr()); flush(stdout())

### Name: plotMA
### Title: Generate a MA-Plot
### Aliases: plotMA

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
plotMA(jscs);


## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ######################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Plot M-A:
##D plotMA(jscs);
##D 
## End(Not run)



cleanEx()
nameEx("readAnnotationData")
### * readAnnotationData

flush(stderr()); flush(stdout())

### Name: readAnnotationData
### Title: Read junctionSeq annotation files produced by QoRTs.
### Aliases: readAnnotationData

### ** Examples

gff.file <- system.file("extdata/cts/withNovel.forJunctionSeq.gff.gz",
                        package="JctSeqExData2");

#Parse the GFF file:
annoData <- readAnnotationData(gff.file);
head(annoData);



cleanEx()
nameEx("readJunctionSeqCounts")
### * readJunctionSeqCounts

flush(stderr()); flush(stdout())

### Name: readJunctionSeqCounts
### Title: Read junctionSeq count files
### Aliases: readJunctionSeqCounts

### ** Examples

########################################
#Set up example data:
decoder.file <- system.file(
                  "extdata/annoFiles/decoder.bySample.txt",
                  package="JctSeqExData2");
decoder <- read.table(decoder.file,
                  header=TRUE,
                  stringsAsFactors=FALSE);
gff.file <- system.file(
            "extdata/tiny/withNovel.forJunctionSeq.gff.gz",
            package="JctSeqExData2");
countFiles <- system.file(paste0("extdata/tiny/",
     decoder$sample.ID,
     "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
     package="JctSeqExData2");
########################################
#Advanced Analysis:

#Make a "design" dataframe:
design <- data.frame(condition = factor(decoder$group.ID));
#Read the QoRTs counts.
jscs = readJunctionSeqCounts(countfiles = countFiles,
           samplenames = decoder$sample.ID,
           design = design,
           flat.gff.file = gff.file
);




cleanEx()
nameEx("runJunctionSeqAnalyses")
### * runJunctionSeqAnalyses

flush(stderr()); flush(stdout())

### Name: runJunctionSeqAnalyses
### Title: Run a JunctionSeq analysis.
### Aliases: runJunctionSeqAnalyses

### ** Examples

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D 
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D 
## End(Not run)



cleanEx()
nameEx("testForDiffUsage")
### * testForDiffUsage

flush(stderr()); flush(stdout())

### Name: testForDiffUsage
### Title: Test Junctions for Differential Junction Usage
### Aliases: testForDiffUsage

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
jscs <- testForDiffUsage(jscs);

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ########################################
##D #Advanced Analysis:
##D 
##D #Make a "design" dataframe:
##D design <- data.frame(condition = factor(decoder$group.ID));
##D #Read the QoRTs counts.
##D jscs = readJunctionSeqCounts(countfiles = countFiles,
##D            samplenames = decoder$sample.ID,
##D            design = design,
##D            flat.gff.file = gff.file
##D );
##D #Generate the size factors and load them into the JunctionSeqCountSet:
##D jscs <- estimateJunctionSeqSizeFactors(jscs);
##D #Estimate feature-specific dispersions:
##D jscs <- estimateJunctionSeqDispersions(jscs);
##D #Fit dispersion function and estimate MAP dispersion:
##D jscs <- fitJunctionSeqDispersionFunction(jscs);
##D #Test for differential usage:
##D jscs <- testForDiffUsage(jscs);
##D #Estimate effect sizes and expression estimates:
##D jscs <- estimateEffectSizes( jscs);
##D 
## End(Not run)



cleanEx()
nameEx("writeBedTrack")
### * writeBedTrack

flush(stderr()); flush(stdout())

### Name: writeBedTrack
### Title: Write splice junction browser tracks
### Aliases: writeExprBedTrack writeSigBedTrack

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
writeExprBedTrack("test.exonCoverage.bed.gz", jscs, 
                  plot.exons = TRUE, plot.junctions = FALSE)

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ######################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Exon coverage:
##D writeExprBedTrack("test.exonCoverage.bed.gz", jscs, 
##D                   plot.exons = TRUE, plot.junctions = FALSE)
##D #Junction coverage:
##D writeExprBedTrack("test.jctCoverage.bed.gz", jscs, 
##D                   plot.exons = FALSE, plot.junctions = TRUE)
##D #Both Exon and Junction coverage:
##D writeExprBedTrack("test.featureCoverage.bed.gz", jscs)
##D 
##D #p-values of significant features:
##D writeSigBedTrack("test.pvals.bed.gz", jscs)
##D 
## End(Not run)



cleanEx()
nameEx("writeCompleteResults")
### * writeCompleteResults

flush(stderr()); flush(stdout())

### Name: writeCompleteResults
### Title: Produce output data files, given annotation files and DEXSeq
###   exonCountSet object and DEXSeq results data.
### Aliases: writeCompleteResults

### ** Examples

data(exampleDataSet,package="JctSeqExData2");
#Write results tables and browser track files:
writeCompleteResults(jscs, outfile.prefix = "./results.");

## Not run: 
##D ########################################
##D #Set up example data:
##D decoder.file <- system.file(
##D                   "extdata/annoFiles/decoder.bySample.txt",
##D                   package="JctSeqExData2");
##D decoder <- read.table(decoder.file,
##D                   header=TRUE,
##D                   stringsAsFactors=FALSE);
##D gff.file <- system.file(
##D             "extdata/cts/withNovel.forJunctionSeq.gff.gz",
##D             package="JctSeqExData2");
##D countFiles <- system.file(paste0("extdata/cts/",
##D      decoder$sample.ID,
##D      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
##D      package="JctSeqExData2");
##D ######################
##D #Run example analysis:
##D jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
##D            sample.names = decoder$sample.ID,
##D            condition=factor(decoder$group.ID),
##D            flat.gff.file = gff.file,
##D            analysis.type = "junctionsAndExons"
##D );
##D ########################################
##D 
##D #Write results tables and browser track files:
##D writeCompleteResults(jscs, outfile.prefix = "./results.");
##D 
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
