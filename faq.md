    > v0.99.8 (Revised Thu Jan 21 13:51:03 EST 2016)

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

You can currently cite the JunctionSeq methodology preprint paper as:

Hartley SW, Mullikin JC. [Detection and Visualization of Differential Exon and Splice Junction Usage in RNA-Seq Data with JunctionSeq.](http://arxiv.org/abs/1512.06038) [arXiv](http://arxiv.org) preprint [arXiv:1512.06038](http://arxiv.org/abs/1512.06038). 2015 Dec 18.

##BASIC INSTALLATION (from source):
JunctionSeq can be installed automatically from source using the R commands:

    source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    JS.install();

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
