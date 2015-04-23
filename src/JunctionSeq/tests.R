

############################################
### BASE EXAMPLE CODE:

library(JunctionSeq);

sample.files <- system.file(c("extdata/withNovel/C_D1.counts.txt.gz",
                               "extdata/withNovel/C_D2.counts.txt.gz",
                               "extdata/withNovel/C_D3.counts.txt.gz",
                               "extdata/withNovel/C_N1.counts.txt.gz",
                               "extdata/withNovel/C_N4.counts.txt.gz",
                               "extdata/withNovel/C_N5.counts.txt.gz"
                               ),
                               package="JunctionSeq", mustWork=TRUE);

anno.file <- system.file("extdata/anno.rat.exampleGenes.flattened_with_novel.gtf.gz",package="JunctionSeq",mustWork=TRUE);

sample.names <- c("C_D1","C_D2","C_D3","C_N1","C_N2","C_N3" )
condition <- factor(c("D","D","D","N","N","N"));

anno.data <- read.anno.data(anno.file, TRUE);
typeof(anno.data$start);

####Part 1:

rda <- run.dexseq.analyses( sample.files = sample.files,
                            sample.names = sample.names, 
                            outfile.prefix = NULL,
                            condition = condition, 
                            saveState = FALSE,
                            annofile = anno.file,
                            verbose = TRUE,
                            use.splice.sites = TRUE, use.novel.splice.sites = TRUE, use.exons = FALSE);

ecs <- rda[[1]];
res <- rda[[2]];

merged.data <- generate.complete.results(ecs, res, outfile.prefix=NULL, verbose = TRUE);

high.sig.thresh <- 0.05

very.sig <- which(merged.data$padjust < high.sig.thresh);
length(very.sig)
very.sig.genes <- unique(as.character(merged.data$geneID[very.sig]))

FDR=0.05

very.sig.genes

bap <- build.all.plots(gene.list=very.sig.genes,merged.results.data=merged.data, ecs=ecs,  anno.file=anno.file,
   FDR.threshold = FDR, condition=condition, color = NULL, use.vst = TRUE, 
   outfile.prefix = "test.",
   variance.plot=TRUE, ma.plot=TRUE, verbose = TRUE,
   plot.exon.results = TRUE, plot.splice.results = TRUE, plot.novel.splice.results = TRUE);


############################################
### MARKED-UP EXAMPLE CODE:

\examples{
  library("JunctionSeq");

}