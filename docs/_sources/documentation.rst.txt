
************
Introduction
************

CoverView is a simple quality control tool for assessing the coverage and quality of mapped reads (.bam file)
in a user-specified targeted panel or whole exome (.bed file). CoverView reports QC metrics in multiple output
files with increasing levels of detail from a chromosome level summary to per base profiles. It also flags regions
that do not pass pre-defined quality requirements.


************
Installation
************

Stable releases of CoverView can be downloaded from Github `here <https://github.com/RahmanTeamDevelopment/CoverView/releases>`_
in either ``.zip`` or ``.tar.gz`` format. To unpack these run one of the following commands::

	unzip CoverView_X.X.X.zip

or::

	tar -xvzf CoverView_X.X.X.tar.gz

and then you can install CoverView with the following commands::

    cd CoverView
    ./install.sh

CoverView uses ``virtualenv`` and ``pip`` to manage all its extra dependencies, which means that it will not clutter up your system by installing
things globally. Everthing it installs will go into a sub-directory in the ``CoverView`` directory (specifically, ``CoverView_X.X.X/env``). If
you delete CoverView then everything it has installed will also be deleted. Once the installation script has finished successfully,
CoverView is ready for use. 


Dependencies
============

To install and run CoverView v1.3.0 you will need `Git <https://git-scm.com>`_, `Python <https://www.python.org>`_ (only
version 2.7.X is currently supported), and `Virtualenv <https://virtualenv.pypa.io/en/stable/>`_. 


*****************
Running CoverView
*****************

Once downloaded and correctly installed, CoverView can be run with the following simple command::

    env/bin/coverview -c config.json -i input.bam –b panel.bed -o output

or you can optionally supply a transcript database::

    env/bin/coverview -c config.json -i input.bam –b panel.bed -o output -t transcript_database.gz

By default, CoverView takes four command line arguments: the name of the configuration file (-c), the
name of the input bam file (-i), the name of a bed file (-b) and the output file name prefix (-o). 

* The input bam file (-i) should follow the `BAM format <http://samtools.github.io/hts-specs/SAMv1.pdf>`_ and should contain the mapped reads. The .bai index file should also be present in the same directory.
* The bed file (-b) should follow the `BED format <http://genome.ucsc.edu/FAQ/FAQformat>`_ with each record corresponding to a region of interest (e.g. exon).
* The configuration file (-c) should follow the `JSON <http://www.json.org>`_ format  and should contain the user-specified settings (see possible configuration options in Section 4).
* The transcript database (-t) is optional, but if supplied must be compressed with `bgzip <http://www.htslib.org/doc/tabix.html>`_ and indexed with `tabix <http://www.htslib.org/doc/tabix.html>`_


Running without bed file
========================

The -b command line flag is optional. If a bed file is not specified, only a simplified chromosome level read count summary will be outputted, e.g.::

    env/bin/coverview -c config.json -i input.bam -o output


******************
Configuration File
******************

The configuration file uses the JSON format. An example configuration is shown below.

.. highlight:: json

::

	{
	    "count_duplicate_reads": true,
	    "outputs": 
		{ 
			"regions": true, 
			"profiles": true 
		},
	    "low_bq": 10,
	    "low_mq": 20,
	    "pass": 
		{ 
			"MINQCOV_MIN": 50, 
			"MAXFLMQ_MAX": 0.05, 
			"MAXFLBQ_MAX": 0.15 
		},
	    "transcript":  
		{
			"regions": false, 
			"profiles": false, 
			"poor": true 
		}
	}

The following options may be specified in the configuration file

.. sidebar:: Future Changes to CoverView Configuration
    
    CoverView currently uses JSON (Javascript Object Notation) for configuration of complex input options. Whilst this is a very 
    flexible format, it is not always the easiest for users to write. Therefore we will be switching
    to a new configuration style based on the .ini file format in a future release.


.. csv-table::
    :header: "Option", "Type", "Default Value", "Effect"

    count_duplicate_reads,  Boolean, true, if true then duplicate reads are included in the analysis. CoverView counts reads as duplicates if they have the duplicate bit set in the BAM record. 
    low_bq, Integer, 10, The base quality cut-off used in the FLBQ metrics. Only bases with this value or higher will be counted as high-quality.
    low_mq”,Integer, 20, The mapping quality cut-off used in the FLMQ metrics. Only reads with this mapping quality or higher will be counted as high-quality.
	outputs {"regions"}, Boolean, true,  If this is true then the _regions.txt output file will be written.
	outputs {"profiles"}, Boolean, true,  If this is true then the _profiles.txt and _poor.txt output files will be written.
	direction, Boolean, false, If this is true then summary metrics and profiles are output for forward and reverse-stranded reads separately
	transcript {"regions"}, Boolean, true, If this is true then transcript coordinates are reported in the _regions.txt file (N.B. this options requires that a transcript database be provided)
	transcript {"profiles"}, Boolean, true, If this is true then transcript coordinates are reported in the _profiles.txt file (N.B. this options requires that a transcript database be provided)

 
* “pass” (value = JSON object): requirements a region has to satisfy to be flagged as “passed”. Can contain the following numerical fields:
    *   “MINCOV_MIN":  minimum value allowed for MINCOV metrics
    *   “MEDCOV_MIN":  minimum value allowed for MEDCOV metrics
    *   “MINQCOV_MIN": minimum value allowed for MINQCOV metrics
    *   “MEDQCOV_MIN": minimum value allowed for MEDQCOV metrics
    *   “MAXFLMQ_MAX”: maximum value allowed for MAXFLMQ metrics
    *   “MAXFLBQ_MAX”: maximum value allowed for MAXFLBQ metrics

For instance, if “pass”: {“MINQCOV_MIN”: 15, “MAXFLMQ_MAX”: 0.2} is set, regions with MINQCOV>=15 and MAXFLMQ<=0.2 are going to be flagged as passed.

* “only_fail_profiles” (value = Boolean): if true, the _profiles output file will only include failed regions

Note that a template configuration file (config.json) is provided in the CoverView package.


******
Output
******


Chromosome level summary
========================

The <outputprefix>_summary.txt output file provides chromosome level summary (read counts) and contains the following 4 columns:

* Chromosome name
* Total read count (RC): total number of reads mapped to the chromosome
* Read count in targeted regions (RCIN): total number of reads mapping to the chromosome that overlap targeted regions from the bed file 
* Read count outside of targeted regions (RCOUT): total number of reads mapping to the chromosome that do not overlap targeted regions from the bed file 

In addition to the list of chromosomes, the outputted table also reports the mapped, unmapped and total read counts for the whole dataset.


Per base profiles
=================

The <outputprefix>_profiles.txt output file provides per base profiles for the targeted regions of interest. Each position is described in a 
separate line with the following 8 compulsory columns:

* Chromosome 
* Position
* Coverage (COV): number of reads covering the position 
* Quality coverage (QCOV): number of reads covering the position with a read mapping quality larger than the threshold set by the “low_mq” configuration flag and a base quality larger than the threshold set by the “low_bq” configuration flag
* Median base quality (MEDBQ): median base quality of all read bases mapping to the position 
* Fraction of low base quality (FLBQ): fraction of read bases mapping to the position with a base quality smaller or equal than the cutoff set by the “low_bq” configuration flag
* Median mapping quality (MEDMQ): median mapping quality of all reads covering the position 
* Fraction of low mapping quality (FLMQ): fraction of reads covering the position with a mapping quality smaller or equal than the cutoff set by the “low_mq” configuration flag

If set in the configuration file (see “transcript” key in Section 4), an additional column named “Transcript_coordinate” is included providing
the transcript coordinate of the position with regards to the overlapping transcript. In case the position overlaps with multiple transcripts,
the coordinates in all transcripts are reported separated by commas. Transcripts data are read from the user-specified transcript database (see “transcript_db” key in Section 4).

Finally, if directionality information is requested in the configuration file (see “direction” key in Section 4), 12 additional columns are added to the _profiles.txt file: 

* Columns COV+, QCOV+, MEDBQ+, FLBQ+, MEDMQ+ and FLMQ+ provide the same metrics as COV, QCOV, MEDBQ, FLBQ, MEDMQ and FLMQ defined above, however, considering only forward-stranded reads. 
* Columns COV-, QCOV-, MEDBQ-, FLBQ-, MEDMQ- and FLMQ- provide the same information, considering only reverse-stranded reads.


Summary metrics for targeted regions
====================================

The <outputprefix>_regions.txt output file provides a number of different metrics summarizing the per base profiles of each region.
These summary metrics give information on the overall quality of each region. In addition, regions are marked as ‘PASS’ or ‘FAIL’ based
on the requirements set in the configuration file (see “pass” key in Section 4). Each line in the file corresponds to a region described
by the following 12 columns:

* Region name 
* Chromosome 
* Start position of region 
* End position of region 
* ‘PASS’ or ‘FAIL’: Does the region pass the user-specified requirements?
* Read count (RC): Total number of reads overlapping with the region 
* Median coverage (MEDCOV): Median of coverage (COV) values across all positions in the region 
* Minimum coverage (MINCOV): Minimum of coverage (COV) values across all positions in the region 
* Median quality coverage (MEDQCOV): Median of quality coverage (QCOV) values across all positions in the region 
* Minimum quality coverage (MINQCOV): Minimum of quality coverage (QCOV) values across all positions in the region
* Maximum fraction of low mapping quality (MAXFLMQ): Maximum of FLMQ values across all positions in the region
* Maximum fraction of low base quality (MAXFLBQ): Maximum of FLBQ values across all positions in the region

Note that the MEDCOV, MINCOV, MEDQCOV, MINQCOV, MAXFLMQ and MAXFLBQ values are derived from the per-base COV, QCOV, FLMQ and FLBQ 
profiles defined in Section 5.2. The region name in the first column is taken from the 4th column of the BED file. If there are
multiple regions in the BED file with the same name in their 4th column (e.g. the regions correspond to different exons of the 
same gene), CoverView adds an index to the region names joined by an underscore. For example, multiple regions of the BRCA2 gene
would be referred to as BRCA2_1, BRCA2_2, BRCA2_3, etc.

If set in the configuration file (see “transcript” key in Section 4), two additional columns named “Start_transcript” and
“End_transcript” are included providing the transcript coordinates of the start and end positions of the region with regards to
overlapping transcripts.

Finally, if directionality information is requested in the configuration file (see “direction” key in Section 4), 12 additional columns are added to the _region.txt file: 

* Columns MEDCOV+, MINCOV+, MEDQCOV+, MINQCOV+, MAXFLMQ+ and MAXFLBQ+ provide the same metrics as MEDCOV, MINCOV, MEDQCOV, MINQCOV, MAXFLMQ and MAXFLBQ defined above, however, considering only forward-stranded reads. 
* Columns MEDCOV-, MINCOV-, MEDQCOV-, MINQCOV-, MAXFLMQ- and MAXFLBQ- provide the same information, considering only reverse-stranded reads.


Poor quality ranges
===================

If the _profiles.txt file is outputted (see “output” key in Section 4), an additional file named <outputprefix>_poor.txt is also created.
The _poor.txt file provides a comprehensive list of all continuous ranges within the regions of interest with QCOV<15 for all bases (referred
to as 'poor quality' ranges). Note that multiple such ranges may exist in a single region. Each line in the file corresponds to a 'poor quality'
range with the following 6 columns:

* Region name 
* Chromosome 
* Start position of region 
* End position of region 
* Start coordinate of region in transcript
* End coordinate of region in transcript

In case the start or end position overlaps with multiple transcripts, the coordinates in all transcripts are reported separated by commas.
