# rmats2sashimiplot


Requirements
------------

Install Python 2.6.x or Python 2.7.x
Install Samtools.

setup.py will automatically install matplotlib which is required to run this program.

rmats2sashimiplot is intended to be used in a Unix-based environment. It has
been tested on Mac OS and Linux.

***BAM file must be sorted before visualization/indexing.***


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
        -t eventType	                Type of event from rMATS result used in the 											analysis.
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
        --group-info                    If user want to divide samples into groups,
        								they can specify this parameter with a "*.gf" file.
        								Format specification can be found in following
                                        section.
        --min-counts                    If the junction count is smaller(<) than this float   
        								number, then this junction would be omitted. The 
        								default value is 3. If you want to display all the 
        								numbers, then set it as 0.
        --color 						User can customerize the colors of the plot using a
        								sequence of color. The number of the colors are
                                        supposed to be corresponding to that of bam_files. 
                                        eg: --color #FFCC99,#99CC99,#99CC99
        --font-size                     Change the default font size which equals to 8.
        --no-text-background            Transparent text background.
        --hide-number 					Hide the numbers of junction.
        
        -h 								Print this help message and exit(also --help).

Running with sam files:

    $rmats2sashimiplot --s1 s1_rep1.sam[,s1_rep2.sam]* --s2 s2.rep1.sam[,s2.rep2.sam]* -t eventType -e eventsFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir

Running with bam files:

    $rmats2sashimiplot --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir

Using grouping function:

    $rmats2sashimiplot --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir --group-info gf.gf


### Grouping

By using this function, user can divide their samples into different groups. Our program will calculate average inclusion level, average read depth and average number of junction-spanning reads of each group and display them in sashimi plot.
It's extremely helpful when you need to do comparisons between different groups of samples.


### Examples ###
Example of using sam files, drawing sashimiplot by rMATS format event files

    $rmats2sashimiplot --s1 ./testData/S1.R1.test.sam,./testData/S1.R2.test.sam,./testData/S1.R3.test.sam --s2 ./testData/S2.R1.test.sam,./testData/S2.R2.test.sam,./testData/S2.R3.test.sam -t SE -e ./testData/MATS_output/test_PC3E_GS689.SE.MATS.events.txt --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_events_output

![images](https://github.com/Xinglab/rmats2sashimiplot/blob/master/img/plotwithevent.png)

Example of using bam files, drawing sashimiplot by user provided coordinates and gff3 format annotation file

    $rmats2sashimiplot --b1 ./testData/S1.R1.test.bam,./testData/S1.R2.test.bam,./testData/S1.R3.test.bam --b2 ./testData/S2.R1.test.bam,./testData/S2.R2.test.bam,./testData/S2.R3.test.bam -c chr16:-:24944500:24955500:./testData/ensGene.gff3 --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_coordinate_output

![images](https://github.com/Xinglab/rmats2sashimiplot/blob/master/img/plotwithcoor.png)

Example of using grouping function:

    $rmats2sashimiplot --b1 ./testData/S1.R1.test.bam,./testData/S1.R2.test.bam,./testData/S1.R3.test.bam --b2 ./testData/S2.R1.test.bam,./testData/S2.R2.test.bam,./testData/S2.R3.test.bam -t SE -e ./testData/MATS_output/test_PC3E_GS689.SE.MATS.events.txt --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_events_output --group-info grouping.gf

![images](https://github.com/Xinglab/rmats2sashimiplot/blob/master/img/plotwitheventgf.png)

content of grouping.gf:

```
group1name: 1-2 # That means we group ./testData/S1.R1.test.bam and ./testData/S1.R2.test.bam together.
group2name: 3-6 # That means we group ./testData/S1.R3.test.bam, ./testData/S2.R1.test.bam, ./testData/S2.R2.test.bam and ./testData/S2.R3.test.bam together.

```

**Group-info**

This section describes the format of '*.gf' file.

Each line stand for a group, which consists of *group name* and *index* of bam files.

***Important notes:*** Index starts from 1. And the order of bam files corresponds to the order we specified in --b1/b2/s1/s2, i.e. concatenate --b1 and --b2 (or --s1 and --s2 if you're using them.). User can confirm this order by checking variable `bam_files` in `sashimi_plot_settings.txt`(under Sashimi_index_* folder.)

Index should be seperated by `','`. And use `'-'` to express a sequence.

**Eg:**

```
group1: 1,2
group2: 3
group3: 4-5
group4: 4-5,6

or

group1:1,2
group2:3
group3:4-5
group4:4-5,6
```

*(White space allowed.)*


### Test Data ###

Please download and untar the test data from: 

http://www.mimg.ucla.edu/faculty/xing/rmats2sashimiplot/testData.tar


### Output ###

All output sashimiplot pdf files are in Sashimi_plot folder

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