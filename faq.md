    > v1.3.4 (Revised Wed May 18 17:03:38 EDT 2016)

#Frequently Asked Questions

JunctionSeq is an R package designed to detect and assess 
differential usage of exons and splice junctions in High-Throughput, Next-Generation RNA-Seq datasets. 
The methodology is based on the methods used by the [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) 
bioconductor package, originally proposed by [Anders, Reyes, and Huber](http://www.ncbi.nlm.nih.gov/pubmed/22722343).

Help, documentation, and the most recent release of JunctionSeq is available on the 
[JunctionSeq github pages](http://hartleys.github.io/JunctionSeq/index.html).

If you do not find an answer to your question here, you can email the developer at JunctionSeq-Contact (at) list.nih.gov.

Note: the current version of JunctionSeq is ONLY compatible with Bioconductor 3.2 or higher.
For older versions of JunctionSeq compatible with Bioconductor 3.0 and 3.1, see JunctionSeq release v0.5.1.

##HELP AND DOCUMENTATION:
For more information see the [JunctionSeq vignette](http://hartleys.github.io/JunctionSeq/doc/JunctionSeq.pdf) or the 
[online reference documentation](http://hartleys.github.io/JunctionSeq/Rhtml/index.html).

There is also a [comprehensive walkthrough](http://hartleys.github.io/JunctionSeq/doc/example-walkthrough.pdf) of 
the entire analysis pipeline, along with a full 
[example dataset](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/QoRTsPipelineWalkthrough.zip) with 
[example bam files](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/bamfiles.zip).

For help with individual R functions in the R utility, use the R 
command:

    > help(functionname);

For a full listing of all help topics for the R utility, use the R 
command: 

    > help(package="JunctionSeq");

##HOW DO I CITE JUNCTIONSEQ?

You can cite the JunctionSeq methods paper, now published in Nucleic Acids Research:

Hartley SW, Mullikin JC. [Detection and visualization of differential splicing in RNA-Seq data with JunctionSeq](http://nar.oxfordjournals.org/content/early/2016/06/07/nar.gkw501.full). 
[Nucleic Acids Research](http://nar.oxfordjournals.org/). 2016 Jun 1. pii: gkw501. doi: [10.1093/nar/gkw501](http://dx.doi.org/10.1093/nar/gkw501). PubMed PMID: [27257077](http://www.ncbi.nlm.nih.gov/pubmed/27257077).

##BASIC INSTALLATION (from source):
JunctionSeq can be installed automatically from source using the R commands:

    source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    JS.install();

##BIOCONDUCTOR INSTALLATION:
If you are using the devel version of bioconductor (v3.3), then you can install the [Bioconductor version of JunctionSeq](http://bioconductor.org/packages/JunctionSeq/):

    source("https://bioconductor.org/biocLite.R")
    biocLite("JunctionSeq")

##MANUAL INSTALLATION:

If you encounter problems with installation, you can install all prerequisite packages manually using the R commands:

    #Install CRAN packages:
    install.packages("statmod")
    install.packages("plotrix")
    install.packages("stringr")
    install.packages("Hmisc")
    install.packages("locfit")
    #Install Bioconductor packages:
    source("http://bioconductor.org/biocLite.R");
    biocLite();
    biocLite("Biobase");
    biocLite("BiocGenerics");
    biocLite("BiocParallel");
    biocLite("GenomicRanges");
    biocLite("IRanges");
    biocLite("S4Vectors");
    biocLite("genefilter");
    biocLite("geneplotter");
    biocLite("SummarizedExperiment");
    biocLite("DESeq2");

You can then install JunctionSeq itself with the command:

    #Install JunctionSeq (from source):
    install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz", 
                       repos = NULL, 
                       type="source")

JunctionSeq NO LONGER REQUIRES COMPILATION. Thus, Rtools and XCode toolsets are no longer necessary to install JunctionSeq.

##Reducing Memory Usage:

JunctionSeq may use large amounts of RAM. You can reduce the memory usage considerably by reducing the number of cores used with the nCores parameter. 
Unfortunetely, BiocParallel duplicates the entire environment whenever it runs in multicore mode, so amount of RAM required is multiplied by the number of cores in 
use.

The exact memory requirements will vary depending on a large number of different factors, such as genome size/complexity, number of replicates, and the number of novel splice junctions.

##LEGAL:
This software package is licensed under the GNU-GPL v3:

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Portions of this software are "United States Government Work" 
under the terms of the United States Copyright Act.  
It was written as part of the authors' official duties for the 
United States Government and thus those portions cannot be 
copyrighted.  Those portions of this software are freely 
available to the public for use without a copyright notice.  
Restrictions cannot be placed on its present or future use.

Although all reasonable efforts have been taken to ensure the 
accuracy and reliability of the software and data, the National 
Human Genome Research Institute (NHGRI) and the U.S. Government 
does not and cannot warrant the performance or results that may 
be obtained by using this software or data.  NHGRI and the U.S. 
Government disclaims all warranties as to performance, 
merchantability or fitness for any particular purpose.
