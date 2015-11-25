    > v0.6.16 (Revised Wed Nov 25 14:52:43 EST 2015)

[JunctionSeq](http://hartleys.github.io/JunctionSeq/) is an R package for detection and visualization of differential usage of 
exons and splice junctions in High-Throughput, Next-Generation RNA-Seq datasets. 
The methodology is heavily based on the [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) 
bioconductor package, originally proposed by [Anders, Reyes, and Huber](http://www.ncbi.nlm.nih.gov/pubmed/22722343). More information
is available on the [JunctionSeq GitHub Pages](http://hartleys.github.io/JunctionSeq/).

The core advantage of JunctionSeq over other similar tools is that it provides a powerful automated tools for
generating readable and interpretable plots and tables to facilitate the interpretation of the results.
An example results report is available [here](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/exampleResults/testForDU.html). 
An example set of browser tracks from this same dataset is available 
[here](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=stephen.hartley&hgS_otherUserSessionName=rn6_pipelineWalkthrough_finalTracks), which uses [this trackhub](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/trackHub/index.html).

Issues, bug reports, or feature requests can be posted to the 
[github issues page](https://github.com/hartleys/JunctionSeq/issues).
The developers can be contacted at JunctionSeq-Contact (at) list.nih.gov.

Note: the current version of JunctionSeq is ONLY compatible with Bioconductor 3.2 or higher.
For older versions of JunctionSeq compatible with Bioconductor 3.0 and 3.1, see JunctionSeq release v0.5.1.

##HELP AND DOCUMENTATION:

* [The JunctionSeq user manual](doc/JunctionSeq.pdf): An introduction to differential usage analysis with JunctionSeq.
* [The online reference documentation](Rhtml/index.html): The complete JunctionSeq R documentation.
* [A comprehensive walkthrough](doc/example-walkthrough.pdf): describes the entire analysis pipeline from aligned reads through analysis. This includes analysis with other tools such as DESeq2, edgeR, and DEXSeq, as well as the creation of plots and browser tracks.
* [Example dataset](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/QoRTsPipelineWalkthrough.zip): used with the walkthrough so you can follow along. The [example bam files](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/bamfiles.zip) are separate (~1gb download).
* [An example results report](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/exampleResults/testForDU.html): An example of the plots and html files produced by JunctionSeq (created using the example dataset).
* An example set of browser tracks from this same dataset is available 
[here](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=stephen.hartley&hgS_otherUserSessionName=rn6_pipelineWalkthrough_finalTracks). The trackhub used to generate this session can be found [here](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/trackHub/index.html)
* [Frequently Asked Questions](faq.html): If you still have questions, check here.

Note: *The example dataset and results are for testing and demonstration purposes only.* The samples and annotation have been heavily modified and down-sampled 
both to test artificial edge cases and to provide smaller and more portable testing files. 
The results should not be taken as an indication of any biological phenomenon.

For help with individual R functions in the R utility, use the R 
command:

    > help(functionname);

For a full listing of all help topics for the R utility, use the R 
command: 

    > help(package="JunctionSeq");

##INSTALLATION:

JunctionSeq is primarily designed for use on high performance linux workstations or cluster environments.

JunctionSeq can be installed automatically in R using the command:

    > source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    > JS.install();


The splice, gene, and exon read-counts required by JunctionSeq can be created 
using the QoRTs software package, available [here](http://hartleys.github.io/QoRTs/index.html).

###Installation on Windows or OSX:

In order to successfully install from source on windows, you must have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 
installed.

Installation on OSX (mac) is not currently supported. In order to successfully install from source on OSX, 
you must have the Xcode tools installed, along with several other R 
development tools, found [here](https://cran.r-project.org/bin/macosx/tools/). According to some reports these
packages are actually out-of-date and may not function properly. See [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)
for a description of some of these problems. JunctionSeq has not been tested on OSX at this time, although testing is currently pending.

##EXAMPLE DATA:
Another example dataset, used in the vignette, is packaged as an R package, and can be installed with the command:
    
    > install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqExData2_LATEST.tar.gz", 
                      repos = NULL, 
                      type="source")

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
