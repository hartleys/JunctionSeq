# JunctionSeq: Detection of Differential Splice Junction Usage in RNA-Seq Data
v0.3.5
Revised Thu Apr 23 16:47:45 EDT 2015

JunctionSeq is an R package designed to detect and assess 
differential splice junction usage in High-Throughput, 
Next-Generation RNA-Seq datasets. 
The methodology is similar to the methods used by the DEXSeq 
bioconductor package.

The most recent release of JunctionSeq is available on the 
[JunctionSeq github page](http://github.com/hartleys/JunctionSeq).
Additional help and documentation can be found 
[here](http://dl.dropboxusercontent.com/u/103621176/JunctionSeq/helpDocs/index.html).
Issues, bug reports, or feature requests can be posted to the 
[github issues page](https://github.com/hartleys/JunctionSeq/issues).

##HELP AND DOCUMENTATION:
For more information, see the [JunctionSeq vignette](doc/JunctionSeq.pdf) or the 
[online reference documentation](Rhtml/index.html) (or as a [pdf](doc/JunctionSeq-reference.pdf)).

Help is also available from within R, accessed via the command:

    > help(functionname);

For a full listing of all help topics, use the R 
command: 

    > help(package="JunctionSeq");

##INSTALLATION:
JunctionSeq is dependent on a number of other R packages. These 
dependencies can be installed using the R commands:

    > install.packages("statmod")
    > install.packages("plotrix")
    > install.packages("stringr")
    > source("http://bioconductor.org/biocLite.R")
    > biocLite()
    > biocLite("Biobase")

JunctionSeq can be installed in R using the command:

    > install.packages("JunctionSeq_0.3.5.tar.gz", repos = NULL, type="source")

The splice-junction counts required by JunctionSeq can be created 
using the QoRTs software package, available 
[here](http://github.com/hartleys/QoRTs).

##EXAMPLE DATA:
The example dataset can be found on the github main page, and can be
installed with the command:
    
    > install.packages(JctSeqExData_0.3.5.tar.gz", repos = NULL, type="source")

##LEGAL:
This software package is licensed under the GNU-GPL v3. A full
copy of the GPL v3 follows:

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
