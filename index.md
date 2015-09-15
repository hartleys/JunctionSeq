    > v0.4.07 (Revised Tue Sep 15 13:39:27 EDT 2015)

JunctionSeq is an R package designed to detect and assess 
differential usage of exons and splice junctions in High-Throughput, Next-Generation RNA-Seq datasets. 
The methodology is based on the methods used by the [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) 
bioconductor package, originally proposed by [Anders, Reyes, and Huber](http://www.ncbi.nlm.nih.gov/pubmed/22722343).

Help, documentation, and the most recent release of JunctionSeq is available on the 
[JunctionSeq github pages](http://hartleys.github.io/JunctionSeq/index.html).

**JunctionSeq is currently in an early alpha stage, and is not yet intended for general use.**

Issues, bug reports, or feature requests can be posted to the 
[github issues page](https://github.com/hartleys/JunctionSeq/issues).

##HELP AND DOCUMENTATION:
For more information see the [JunctionSeq vignette](http://hartleys.github.io/JunctionSeq/helpDocs/doc/JunctionSeq.pdf) or the 
[online reference documentation](http://dl.dropboxusercontent.com/u/103621176/JunctionSeq/helpDocs/Rhtml/index.html).

There is also a [comprehensive walkthrough](http://hartleys.github.io/JunctionSeq/helpDocs/doc/example-walkthrough.pdf) of 
the entire analysis pipeline, along with a full 
[example dataset](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/QoRTsPipelineWalkthrough.zip) with 
[example bam files](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/bamfiles.zip).

For help with individual R functions in the R utility, use the R 
command:

    > help(functionname);

For a full listing of all help topics for the R utility, use the R 
command: 

    > help(package="JunctionSeq");

##INSTALLATION:

JunctionSeq can be installed automatically in R using the command:

    > source("http://hartleys.github.io/JunctionSeq/install/JunctionSeq.install.R");

In order to successfully install R on windows, you must have Rtools installed.

The splice, gene, and exon read-counts required by JunctionSeq can be created 
using the QoRTs software package, available [here](http://hartleys.github.io/QoRTs/index.html).

##EXAMPLE DATA:
The example dataset can be found on the github repository and can be 
installed with the command:
    
    > install.packages(JctSeqExData2_0.4.07.tar.gz", repos = NULL, type="source")

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
