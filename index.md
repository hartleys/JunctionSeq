    > v1.3.4 (Revised Wed May 18 17:03:38 EDT 2016)

[JunctionSeq](http://hartleys.github.io/JunctionSeq/) is a [Bioconductor](https://www.bioconductor.org/) package for detection and visualization of differential usage of 
exons and splice junctions in High-Throughput, Next-Generation RNA-Seq datasets. 
The methodology is heavily based on the [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) 
bioconductor package, originally proposed by [Anders, Reyes, and Huber](http://www.ncbi.nlm.nih.gov/pubmed/22722343).

One major advantage of JunctionSeq over other similar tools is that it provides a powerful automated tools for
generating readable and interpretable plots and tables to facilitate the interpretation of the results.
An example results report is available [here](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/exampleResults/testForDU.html). 
An example set of browser tracks from this same dataset is available 
[here](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=stephen.hartley&hgS_otherUserSessionName=rn6_pipelineWalkthrough_finalTracks), which uses [this trackhub](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/trackHub/index.html).

Issues, bug reports, or feature requests can be posted to the 
[github issues page](https://github.com/hartleys/JunctionSeq/issues).
The developers can be contacted at JunctionSeq-Contact (at) list.nih.gov.

**JunctionSeq is now part of [Bioconductor](https://www.bioconductor.org/), and will be included in the next release (3.3).** 
You can install JunctionSeq using the devel branch of Bioconductor at the [JunctionSeq bioconductor page](http://bioconductor.org/packages/JunctionSeq/). Alternatively, you 
can use the installation instructions below to install the most recent version of JunctionSeq onto the current Bioconductor release (3.2).

##HELP AND DOCUMENTATION:

* [The JunctionSeq user manual/vignette](doc/JunctionSeq.pdf): An introduction to analysis with JunctionSeq. Read this first.
* [The online reference documentation](Rhtml/index.html): The complete JunctionSeq R documentation.
* [An example results report](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/exampleResults/testForDU.html): An example of the plots and html files produced by JunctionSeq (created using the example dataset, see below).
* An example set of browser tracks from this same dataset is available 
[here](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=stephen.hartley&hgS_otherUserSessionName=rn6_pipelineWalkthrough_finalTracks). The trackhub used to generate this session can be found [here](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/trackHub/index.html)
* For Advanced Users: [A comprehensive walkthrough](doc/example-walkthrough.pdf) describing the entire analysis pipeline from aligned reads through analysis. This includes analysis with other tools such as DESeq2, edgeR, and DEXSeq, as well as the creation of plots and browser tracks.
* The [Example dataset](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/QoRTsPipelineWalkthrough.zip) used with the walkthrough, so you can follow along (~280mb download). 
* The [example bam files](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/bamfiles.zip) are separate (~1gb download).
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

##CITING JUNCTIONSEQ:

You can currently cite the JunctionSeq methods paper preprint as:

Hartley SW, Mullikin JC. [Detection and Visualization of Differential Exon and Splice Junction Usage in RNA-Seq Data with JunctionSeq.](http://arxiv.org/abs/1512.06038) [arXiv](http://arxiv.org) preprint [arXiv:1512.06038](http://arxiv.org/abs/1512.06038). 2015 Dec 18.

##INSTALLATION:

JunctionSeq can be installed automatically in R using the command:

    source("http://hartleys.github.io/JunctionSeq/install/JS.install.R");
    JS.install();

See the [FAQ](faq.html) for advanced installation options.

The splice, gene, and exon read-counts required by JunctionSeq can be created 
using the QoRTs software package, available [here](http://hartleys.github.io/QoRTs/index.html).

##BIOCONDUCTOR INSTALLATION:
If you are using the devel version of bioconductor (v3.3), then you can install the [Bioconductor version of JunctionSeq](http://bioconductor.org/packages/JunctionSeq/):

    source("https://bioconductor.org/biocLite.R")
    biocLite("JunctionSeq")

##EXAMPLE DATA:
Another example dataset, used in the vignette, is packaged as an R package, and can be installed with the command:
    
    install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqData_LATEST.tar.gz", 
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
