MapMap: An exact read mapper for very large databases with short update time
===================================================================================

Incorparating Sequence Mappability into the read mapping process using Optimum Search Schemes. As foundation of this tool Dream-Yara was used. 

Software requirements
~~~~~~~~~~~~~~~~~~~~~

**A modern C++11 compiler with OpenMP 3.0 extensions is required to build Yara. If unsure, use GNU G++ 4.9 or newer.**

* Git.
* CMake 3.2 or newer.
* G++ 4.9 or newer.

Download
~~~~~~~~

MapMap sources downloaded by executing:

::

  $ git clone --recurse-submodules https://github.com/svnbgnk/dream_yara.git

Configuration
~~~~~~~~~~~~~

Create a build project by executing CMake as follows:

::
 $ git checkout oss
 $ cd dream_yara/include/seqan/
 $ git remote add upstream https://github.com/svnbgnk/seqan.git
 $ git checkout upstream/mappa 
 $ cd ../../..
 $ mkdir mapmap-build
 $ cd mapmap-build
 $ cmake ../dream_yara

Build
~~~~~

Invoke make as follows:

::

  $ make all


Acquiring hg38.fa
::

 $ DATA/reference
 $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 $ gunzip hg38.fa.gz
 $ mkdir bin
 $ cp hg38.fa bin/0.fa



Usage
-----

Building the Index (requires 200GB of secondary memory)
::
 $ dream_yara_indexer --threads 8 --output-prefix DATA/hg38_N_index/ DATA/reference/bin/*.fa -td /srv/public/svnbngk/tmp/

Computing sequence mappability and bit vectors with TH = 10
::
 $ dream_yara_mappability DATA/hg38_N_index/ -b 1 -K 100 -E 3 -T 10 -s 0 -t 20 -o 35 -v -i -O DATA/hg38_N_index/mappability10E3

All mapping with up to 3 errors
::
 $ dream_yara_mapper DATA/hg38_N_index/ DATA/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 3 -o result.sam -vv 

Stratified all-mapping with strata 2 and 3 errors
::
 $ dream_yara_mapper DATA/hg38_N_index/ DATA/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 2 -o result.sam -vv

Mapping with sequence mappability up to 3 errors
::
 $ dream_yara_mapper DATA/hg38_N_index/ DATA/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 3 -m DATA/hg38_N_index/mappability10E3/ -o result.sam -vv



Output format
^^^^^^^^^^^^^

Output files follow the `SAM/BAM format specification <http://samtools.github.io/hts-specs/SAMv1.pdf>`_.
In addition, Yara generates the following optional tags:

+-----+----------------------------------------------------+
| Tag | Meaning                                            |
+=====+====================================================+
| NM  | Edit distance                                      |
+-----+----------------------------------------------------+
| X0  | Number of co-optimal mapping locations             |
+-----+----------------------------------------------------+
| X1  | Number of sub-optimal mapping locations            |
+-----+----------------------------------------------------+
| XA  | Alternative locations: (chr,begin,end,strand,NM;)* |
+-----+----------------------------------------------------+


Contact
-------

For questions or comments, feel free to contact: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>


References
----------
Dadi, T. H., Siragusa, E., Piro, V. C., Andrusch, A., Seiler, E., Renard, B. Y., & Reinert, K. (2018).
DREAM-Yara: An exact read mapper for very large databases with short update time.
BioRxiv, 256354. https://doi.org/10.1101/256354
