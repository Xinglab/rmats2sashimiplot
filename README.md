# rmats2sashimiplot


Requirements
------------

Install Python 2.6.x or Python 2.7.x

setup.py will automatically install matplotlib which is required to run this program.

rmats2sashimiplot is intended to be used in a Unix-based environment. It has
been tested on Mac OS and Linux.


Installation
------------

### Install ###
The package, rmats2sashimiplot is installed by typing:

    python setup.py install

Usage
-----
The following is a detailed description of the options used with rmats2sashimiplot.

### Required Parameters ###
        --s1 s1_rep1.sam[,s1_rep2.sam]	Mapping results for the sample_1 in sam format.
                                        Replicates must be in a comma separated list
                                        (Only if using sam).
        --s2 s2.rep1.sam[,s2.rep2.sam]	Mapping results for the sample_2 in sam format.
                                        Replicates must be in a comma separated list
                                        (Only if using sam).
        --b1 s1_rep1.bam[,s1_rep2.bam]	Mapping results for the sample_1 in bam format.
                                        Replicates must be in a comma separated list
                                        (Only if using bam).
        --b2 s2.rep1.bam[,s2.rep2.bam]	Mapping results for the sample_2 in bam format.
                                        Replicates must be in a comma separated list
                                        (Only if using bam).
        -t eventType	                Type of event from rMATS result used in the analysis.
                                        eventType is 'SE', 'A5SS', 'A3SS', 'MXE' or 'RI'.
                                        'SE' is for skipped exon events, 'A5SS' is for
                                        alternative 5' splice site events, 'A3SS' is for
                                        alternative 3' splice site events, 'MXE' is for
                                        mutually exclusive exons events and 'RI' is for
                                        retained intron events (Only if using rMATS format
                                        result as event file).
        -e eventsFile	                The rMATS output event file (Only if using rMATS
                                        format result as event file).
        -c coordinate:annotaionFile	    The coordinate of genome region and an annotation
                                        of genes and transcripts in GFF3 format. Coordinate
                                        and annotation file must be colon separated
                                        (Only if using coordinate and annotaion file).
        --l1 SampleLabel1	            The label for first sample.
        --l2 SampleLabel2	            The label for second sample.
        -o outDir	                    The output directory.
    
        Optional:
        --exon_s <int>	                The size of scale down exons. The default is 1.
        --intron_s <int>	            The size of scale down introns. For example, if
                                        -intron_s is 5, it means the size of intron is 5:1
                                        (if the real size of intron is 5, the size in the
                                        plot will be scaled down to 1). The default is 1.
        -h                              Print this help message and exit (also --help).

Running with sam files:

```shell
    $rmats2sashimiplot --s1 s1_rep1.sam[,s1_rep2.sam]* --s2 s2.rep1.sam[,s2.rep2.sam]* -t eventType -e eventsFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir
```

Running with bam files:

```shell
    $rmats2sashimiplot --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir
```

### Examples ###
Example of using sam files, drawing sashimiplot by rMATS format event files

```shell
    $rmats2sashimiplot --s1 ./testData/S1.R1.test.sam,./testData/S1.R2.test.sam,./testData/S1.R3.test.sam --s2 ./testData/S2.R1.test.sam,./testData/S2.R2.test.sam,./testData/S2.R3.test.sam -t SE -e ./testData/MATS_output/test_PC3E_GS689.SE.MATS.events.txt --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_events_output
```

Example of using bam files, drawing sashimiplot by user provided coordinates and
gff3 format annotation file

```shell
    $rmats2sashimiplot --b1 ./testData/S1.R1.test.bam,./testData/S1.R2.test.bam,./testData/S1.R3.test.bam --b2 ./testData/S2.R1.test.bam,./testData/S2.R2.test.bam,./testData/S2.R3.test.bam -c chr2:+:10090000:10110000:./testData/ensGene.gff3 --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_coordinate_output
```

### Test Data ###
Please download and untar the test data from: 

http://www.mimg.ucla.edu/faculty/xing/rmats2sashimiplot/testData.tar

### Output ###
All output sashimiplot pdf files are in -o outDir folder.

Contacts and bug reports
------------------------
Yi Xing
yxing@ucla.edu

Zhijie Xie
shiehshiehzhijie@gmail.com

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Zhijie Xie, Yu-Ting Tseng, Yi Xing

Zhijie Xie, Yu-Ting Tseng, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.