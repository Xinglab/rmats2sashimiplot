#!/usr/bin/python
from sets import Set

### import necessary libraries
import re,os,sys,logging,time,datetime;

startTime = time.time();

s1=''; 		## sample_1
s2=''; 		## sample_2
l1=''; 		## label_1
l2=''; 		## label_2
exon_s=1; 		## exon_scale
intron_s=1; 		## intron_scale
samDir=''; 	## path of sam files
outDir=''; 	## path of output files
event_type=''; 		## event type
bamFile=0; ## by default, no bam file
eventsFile=0; ## by default, no events file

helpStr = "Required parameters:                                                                         \n" +\
          "-s1 s1_rep1.sam[,s1_rep2.sam]   Mapping results for the sample_1 in sam format.              \n" +\
          "                                Replicates must be in a comma separated list                 \n" +\
          "                                (Only if using sam).                                         \n" +\
          "-s2 s2.rep1.sam[,s2.rep2.sam]   Mapping results for the sample_2 in sam format.              \n" +\
          "                                Replicates must be in a comma separated list                 \n" +\
          "                                (Only if using sam).                                         \n" +\
          "-b1 s1_rep1.bam[,s1_rep2.bam]   Mapping results for the sample_1 in bam format.              \n" +\
          "                                Replicates must be in a comma separated list                 \n" +\
          "                                (Only if using bam).                                         \n" +\
          "-b2 s2.rep1.bam[,s2.rep2.bam]   Mapping results for the sample_2 in bam format.              \n" +\
          "                                Replicates must be in a comma separated list                 \n" +\
          "                                (Only if using bam).                                         \n" +\
          "-t eventType                    Type of event from rMATS result used in the analysis.        \n" +\
          "                                eventType is \'SE\', \'A5SS\', \'A3SS\', \'MXE\' or \'RI\'.  \n" +\
          "                                \'SE\' is for skipped exon events, \'A5SS\' is for           \n" +\
          "                                alternative 5\' splice site events, \'A3SS\' is for          \n" +\
          "                                alternative 3\' splice site events, \'MXE\' is for           \n" +\
          "                                mutually exclusive exons events and \'RI\' is for            \n" +\
          "                                retained intron events (Only if using rMATS format           \n" +\
          "                                result as event file).                                       \n" +\
          "-e eventsFile                   The rMATS output event file (Only if using rMATS             \n" +\
          "                                format result as event file).                                \n" +\
          "-c coordinate:annotaionFile     The coordinate of genome region and an annotation            \n" +\
          "                                of genes and transcripts in GFF3 format. Coordinate          \n" +\
          "                                and annotation file must be colon separated                  \n" +\
          "                                (Only if using coordinate and annotaion file).               \n" +\
          "-l1 SampleLabel1                The label for first sample.                                  \n" +\
          "-l2 SampleLabel2                The label for second sample.                                 \n" +\
          "-o outDir                       The output directory.                                        \n" +\
          "                                                                                             \n" +\
          "Optional:                                                                                    \n" +\
          "-exon_s <int>                   The size of scale down exons. The default is 1.              \n" +\
          "intron_s <int>                  The size of scale down introns. For example, if              \n" +\
          "                                -intron_s is 5, it means the size of intron is 5:1           \n" +\
          "                                (if the real size of intron is 5, the size in the            \n" +\
          "                                plot will be scaled down to 1). The default is 1.            \n" +\
          "-h                              Print this help message and exit (also --help).              \n" +\
          "                                                                                             \n" +\
          "Usage (with sam files):\n" +\
          "rmats2sashimiplot -s1 s1_rep1.sam[,s1_rep2.sam]* -s2 s2.rep1.sam[,s2.rep2.sam]* -t eventType -e eventsFile -l1 SampleLabel1 -l2 SampleLable2 -exon_s exonScale -intron_s intronScale -o outDir  \n\n" +\
          "Example (with sam files):\n" +\
          "rmats2sashimiplot -s1 ./testData/S1.R1.test.sam,./testData/S1.R2.test.sam,./testData/S1.R3.test.sam -s2 ./testData/S2.R1.test.sam,./testData/S2.R2.test.sam,./testData/S2.R3.test.sam -t SE -e ./testData/MATS_output/test_PC3E_GS689.SE.MATS.events.txt -l1 PC3E -l2 GS689 -exon_s 1 -intron_s 5 -o test_events_output  \n\n" +\
          "Usage (with bam files):\n" +\
          "rmats2sashimiplot -b1 s1_rep1.bam[,s1_rep2.bam]* -b2 s2.rep1.bam[,s2.rep2.bam]* -c coordinate:annotaionFile -l1 SampleLabel1 -l2 SampleLable2 -exon_s exonScale -intron_s intronScale -o outDir  \n\n" +\
          "Example (with bam files):\n" +\
          "rmats2sashimiplot -b1 ./testData/S1.R1.test.bam,./testData/S1.R2.test.bam,./testData/S1.R3.test.bam -b2 ./testData/S2.R1.test.bam,./testData/S2.R2.test.bam,./testData/S2.R3.test.bam -c chr2:+:10090000:10110000:./testData/ensGene.gff3 -l1 PC3E -l2 GS689 -exon_s 1 -intron_s 5 -o test_coordinate_output"

def isHelpString(s) :
  norm = s.strip().lower()
  return norm == "help" or norm == "-help" or norm == "--help" or norm == "-h"

### checking out the argument names
validArgList = ['-s1','-b1','-s2','-b2','-t','-e','-c','-l1','-l2','-exon_s','-intron_s','-o','help','-help','--help','-h'];
for argIndex in range(1,len(sys.argv)): ## going through the all parameters
  if(sys.argv[argIndex][0]=='-' and sys.argv[argIndex] not in validArgList): ## incorrect argument
    print ('Not valid argument: %s' % sys.argv[argIndex]);
    print ('Please provide valid arguments.');
    print (helpStr + "\n\n");
    sys.exit();
  elif len(sys.argv) == 0 or (len(sys.argv) == 1 and isHelpString(sys.argv[0])) :
    print (helpStr + "\n\n");
    sys.exit();

for paramIndex in range(1,len(sys.argv)): ## going through the all parameters
  if(sys.argv[paramIndex] == '-s1' or sys.argv[paramIndex] == '-b1'):  ## sample_1
    if (sys.argv[paramIndex] == '-b1'): ## bam file here
      bamFile=1;
    paramIndex += 1;  ## increase index
    s1 = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-s2' or sys.argv[paramIndex] == '-b2'):  ## sample_2
    if (sys.argv[paramIndex] == '-b2'): ## bam file here
      bamFile=1;
    paramIndex += 1;  ## increase index
    s2 = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-t'):  ## event_type
    paramIndex += 1;  ## increase index
    event_type = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-c'  or sys.argv[paramIndex] == '-e'):  ## coordinate or events_file
    if (sys.argv[paramIndex] == '-e'): ## events file here
      eventsFile=1;
    paramIndex += 1;  ## increase index
    events = sys.argv[paramIndex];
  elif(sys.argv[paramIndex] == '-l1'):  ## label_1
    paramIndex += 1;  ## increase index
    l1 = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-l2'):  ## label_2
    paramIndex += 1;  ## increase index
    l2 = sys.argv[paramIndex];
  elif(sys.argv[paramIndex] == '-exon_s'):  ## exon_scale
    paramIndex += 1;  ## increase index
    exon_s = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-intron_s'):  ## intron_scale
    paramIndex += 1;  ## increase index
    intron_s = sys.argv[paramIndex];
  elif (sys.argv[paramIndex] == '-o'):  ## coverage
    paramIndex += 1;  ## increase index
    outDir = sys.argv[paramIndex];

#  else: ### not valid param.. exit
#    print("Not a valid param detected: %s" % sys.argv[paramIndex]);
#    sys.exit();

### checking out the required arguments
if len(sys.argv) == 0 or (len(sys.argv) == 1 and isHelpString(sys.argv[0])) :
    print (helpStr + "\n\n");
    sys.exit();
elif (s1=='' or  s2=='' or events=='' or l1=='' or  l2=='' or outDir==''): ### at least one required param is missing
    print ('Not enough arguments!\n');
    print (helpStr + "\n\n");
    sys.exit();

outPath = os.path.abspath(outDir);
sashimiPath = outPath + '/Sashimi_index';
os.system('mkdir -p '+ sashimiPath);

### setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir+'/log.sashimiPlot_test.'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')

sample_1=s1.split(',');
sample_2=s2.split(',');

if bamFile==1 and ( ((sample_1[0].split('.'))[-1].strip()).upper() !='BAM' or ((sample_2[0].split('.'))[-1].strip()).upper() !='BAM'):
  print "Incorrect file type. Need to provide bam file for -b1 and -b2";
  sys.exit();

if len(sample_1)!=len(sample_2): ## different number of replicates per sample.. wrong!!
    print("Requires the same number of replicates per sample...");
    sys.exit();

if eventsFile==1 and ( ((events.split('.'))[-1].strip()).upper() !='TXT'):
  print "Incorrect events file type. Need to provide rMATS output format txt file for -e";
  sys.exit();

events_name_level = {}


### process input params ####
#
### sam files or bam files
#
###
### 0. convert sam to bam and build index..
###
logging.debug("convert sam to bam and build index..");
logging.debug("################### folder names and associated input files #############");

for fki in range(0,len(sample_1)): ## for each replicate of sample_1
  if bamFile==0:
    os.system('samtools view -Sbh '+sample_1[fki]+' > '+sample_1[fki].replace(".sam","")+'.bam');
    os.system('samtools index '+sample_1[fki].replace(".sam","")+'.bam');
    #logging.debug("sam file is provided"+"\t"+sample_1[fki]);
  else: ## bam file is provided
    os.system('samtools index '+sample_1[fki]);
    #logging.debug("bam file is provided"+"\t"+sample_1[fki]);
  repTempFolder = "SAMPLE_1\REP_"+str(fki+1);
  associatedFile = sample_1[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);

for fki in range(0,len(sample_2)): ## for each replicate of sample_2
  if bamFile==0:
    os.system('samtools view -Sbh '+sample_2[fki]+' > '+sample_2[fki].replace(".sam","")+'.bam');
    os.system('samtools index '+sample_2[fki].replace(".sam","")+'.bam');
    #logging.debug("sam file is provided"+"\t"+sample_2[fki]);
  else: ## bam file is provided
    os.system('samtools index '+sample_2[fki]);
    #logging.debug("bam file is provided"+"\t"+sample_2[fki]);
  repTempFolder = "SAMPLE_2\REP_"+str(fki+1);
  associatedFile = sample_2[fki];
  logging.debug(repTempFolder+"\t"+associatedFile);

logging.debug("#########################################################################\n");

########## functions here... ############

###
### 1. prepare sashimi plot setting file..
###

def prepareSettingFile(gene_no_str): ## get AS events from GTF and SAM files
  logging.debug("prepare sashimi plot setting file..");

  geneSymbol = (gene_no_str.split('_'))[0];
  settingFile = open(outPath + '/Sashimi_index_'+gene_no_str+'/sashimi_plot_settings.txt', 'w');
  settingFile.write("[data]\n");
  settingFile.write("bam_prefix = "+samDir+"\n");
  settingFile.write("miso_prefix = "+samDir+"\n");
  ## setting string for replicates ##
  bam_files_arr1 = [];
  bam_files_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    if bamFile==0:
      bam_files_arr1.append('\"'+sample_1[rr].replace(".sam","")+'.bam\"');
    else: ## bam file is provided
      bam_files_arr1.append('\"'+sample_1[rr]+'\"');
  for rr in range(0,len(sample_2)): ## sample_2
    if bamFile==0:
      bam_files_arr2.append('\"'+sample_2[rr].replace(".sam","")+'.bam\"');
    else: ## bam file is provided
      bam_files_arr2.append('\"'+sample_2[rr]+'\"');
  setting_bam_str = ','.join(bam_files_arr1)+','+','.join(bam_files_arr2);
  settingFile.write("bam_files = ["+setting_bam_str+"]\n");
  settingFile.write("miso_files = ["+setting_bam_str+"]\n");
  settingFile.write("[plotting]\n");
  settingFile.write("fig_width = 8\n");
  if len(sample_1)<5:
    settingFile.write("fig_height = 7\n");
  else:
    settingFile.write("fig_height = 14\n");
  settingFile.write("exon_scale = "+str(exon_s)+"\n");
  settingFile.write("intron_scale = "+str(intron_s)+"\n");
  settingFile.write("logged = False\n");
  settingFile.write("font_size = 8\n");
  settingFile.write("bar_posteriors = False\n");
  settingFile.write("nyticks = 4\n");
  settingFile.write("nxticks = 6\n");
  settingFile.write("show_ylabel = False\n");
  settingFile.write("show_xlabel = True\n");
  settingFile.write("plot_title = \"gene symbol\"\n");
  settingFile.write("plot_label = plot_label\n");
  settingFile.write("show_posteriors = False\n");
  settingFile.write("number_junctions = True\n");
  settingFile.write("resolution = .5\n");
  ## setting string for replicates ##
  colors_arr1 = [];
  colors_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    colors_arr1.append('\"#CC0011\"');
  for rr in range(0,len(sample_2)): ## sample_2
    colors_arr2.append('\"#FF8800\"');
  setting_color_str = ','.join(colors_arr1)+','+','.join(colors_arr2);
  settingFile.write("colors = ["+setting_color_str+"]\n");
  ## setting string for inclusion levels ##
  inclever_str = events_name_level.get(geneSymbol);
  items = inclever_str.split('_');
  inc_level1 = items[0];
  inc_level2 = items[1];
  inc_items1 = inc_level1.split(',');
  inc_items2 = inc_level2.split(',');
  sample_labels_arr1 = [];
  sample_labels_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    file_str1 = sample_1[rr].split('/');
    inc_1 = "{0:.2f}".format(float(inc_items1[rr]));
    sample_labels_arr1.append('\"'+geneSymbol+' '+l1+'-'+str(rr+1)+' IncLevel: '+inc_1+'\"');
  for rr in range(0,len(sample_2)): ## sample_2
    file_str2 = sample_2[rr].split('/');
    inc_2 = "{0:.2f}".format(float(inc_items2[rr]));
    sample_labels_arr2.append('\"'+geneSymbol+' '+l2+'-'+str(rr+1)+' IncLevel: '+inc_2+'\"');
  setting_labels_str = ','.join(sample_labels_arr1)+','+','.join(sample_labels_arr2);
  settingFile.write("sample_labels = ["+setting_labels_str+"]\n");
  settingFile.write("reverse_minus = True");
  settingFile.close();
  logging.debug("done prepare sashimi plot setting file"+sashimiPath+"/sashimi_plot_settings.txt");

  return;
############ end of prepareSettingFile #####

def prepareCoorSettingFile(): ## get AS events from GTF and SAM files
  logging.debug("prepare sashimi plot setting file..");

  settingFile = open(outPath + '/Sashimi_index/sashimi_plot_settings.txt', 'w');
  settingFile.write("[data]\n");
  settingFile.write("bam_prefix = "+samDir+"\n");
  settingFile.write("miso_prefix = "+samDir+"\n");
  ## setting string for replicates ##
  bam_files_arr1 = [];
  bam_files_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    if bamFile==0:
      bam_files_arr1.append('\"'+sample_1[rr].replace(".sam","")+'.bam\"');
    else: ## bam file is provided
      bam_files_arr1.append('\"'+sample_1[rr]+'\"');
  for rr in range(0,len(sample_2)): ## sample_2
    if bamFile==0:
      bam_files_arr2.append('\"'+sample_2[rr].replace(".sam","")+'.bam\"');
    else: ## bam file is provided
      bam_files_arr2.append('\"'+sample_2[rr]+'\"');
  setting_bam_str = ','.join(bam_files_arr1)+','+','.join(bam_files_arr2);
  settingFile.write("bam_files = ["+setting_bam_str+"]\n");
  settingFile.write("miso_files = ["+setting_bam_str+"]\n");
  settingFile.write("[plotting]\n");
  settingFile.write("fig_width = 8\n");
  if len(sample_1)<5:
    settingFile.write("fig_height = 7\n");
  else:
    settingFile.write("fig_height = 14\n");
  settingFile.write("exon_scale = "+str(exon_s)+"\n");
  settingFile.write("intron_scale = "+str(intron_s)+"\n");
  settingFile.write("logged = False\n");
  settingFile.write("font_size = 8\n");
  settingFile.write("bar_posteriors = False\n");
  settingFile.write("nyticks = 4\n");
  settingFile.write("nxticks = 11\n");
  settingFile.write("show_ylabel = False\n");
  settingFile.write("show_xlabel = True\n");
  settingFile.write("plot_title = \"gene symbol\"\n");
  settingFile.write("plot_label = plot_label\n");
  settingFile.write("show_posteriors = False\n");
  settingFile.write("number_junctions = True\n");
  settingFile.write("resolution = .5\n");
  ## setting string for replicates ##
  colors_arr1 = [];
  colors_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    colors_arr1.append('\"#CC0011\"');
  for rr in range(0,len(sample_2)): ## sample_2
    colors_arr2.append('\"#FF8800\"');
  setting_color_str = ','.join(colors_arr1)+','+','.join(colors_arr2);
  settingFile.write("colors = ["+setting_color_str+"]\n");
  sample_labels_arr1 = [];
  sample_labels_arr2 = [];
  for rr in range(0,len(sample_1)): ## sample_1
    sample_labels_arr1.append('\"'+l1+'-'+str(rr+1)+'\"');
  for rr in range(0,len(sample_2)): ## sample_2
    sample_labels_arr1.append('\"'+l2+'-'+str(rr+1)+'\"');
  setting_labels_str = ','.join(sample_labels_arr1)+','+','.join(sample_labels_arr2);
  settingFile.write("sample_labels = ["+setting_labels_str+"]\n");
  settingFile.write("reverse_minus = True");
  settingFile.close();
  logging.debug("done prepare sashimi plot setting file"+sashimiPath+"/sashimi_plot_settings.txt");

  return;
############ end of prepareCoorSettingFile #####

###
### 2. get AS events from rMATS result or user input coordinates and convert to gff3 format file..
###

def drawPlotWithEventsFile(): ## events file is provided
  logging.debug("drawPlotWithEventsFile()");

  fo = open(events,'r');
  w2 = open(outPath+'/Sashimi_index/SE.event.list.txt','w');

  events_no = 0;
  for line in fo:
    geneSymbol = "";
    gene_no_str = "";
    id_str = "";
    if line.startswith('ID'):
      continue;
    events_no += 1;
    items = line.split("\t")
    geneSymbol = items[2]
    geneSymbol = geneSymbol.replace('\"','');
    logging.debug("***** geneSymbol: "+geneSymbol);
    gene_no_str = geneSymbol+'_'+str(events_no);
    sashimiPath = outPath + '/Sashimi_index_'+geneSymbol+'_'+str(events_no);
    os.system('mkdir -p '+ sashimiPath);
    w1 = open(outPath+'/Sashimi_index_'+geneSymbol+'_'+str(events_no)+'/tmp.gff3','w');
    chr = items[3]
    strand = items[4]
    e1st_s = "";
    e1st_e = "";
    e2st_s = "";
    e2st_e = "";
    se_s = "";
    se_e = "";
    up_s = "";
    up_e = "";
    dn_s = "";
    dn_e = "";
    inc_level1 = "";
    inc_level2 = "";
    if event_type!='MXE':
        se_s = str(int(items[5])+1)
        se_e = items[6]
        up_s = str(int(items[7])+1)
        up_e = items[8]
        dn_s = str(int(items[9])+1)
        dn_e = items[10]
        inc_level1 = items[20] ## IncLevel1
        inc_level2 = items[21] ## IncLevel2
    if event_type=='MXE':
        e1st_s = str(int(items[5])+1)
        e1st_e = items[6]
        e2st_s = str(int(items[7])+1)
        e2st_e = items[8]
        up_s = str(int(items[9])+1)
        up_e = items[10]
        dn_s = str(int(items[11])+1)
        dn_e = items[12]
        inc_level1 = items[22] ## IncLevel1
        inc_level2 = items[23] ## IncLevel2
    events_name_level[geneSymbol] = inc_level1 +"_"+ inc_level2
    if strand =='+':
        id_str = chr+":"+up_s+":"+up_e+":"+strand+"@"+chr+":"+se_s+":"+se_e+":"+strand+"@"+chr+":"+dn_s+":"+dn_e+":"+strand
        name_str = geneSymbol+"_"+chr+":"+up_s+":"+up_e+":"+strand+"@"+chr+":"+se_s+":"+se_e+":"+strand+"@"+chr+":"+dn_s+":"+dn_e+":"+strand
        if event_type!='MXE':
            w2.write( "%s\n" % (name_str))
            w1.write( "%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (chr,up_s,dn_e,strand,id_str,name_str))
            w1.write( "%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (chr,se_s,se_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
        if event_type=='MXE':
            w2.write( "%s\n" % (name_str))
            w1.write( "%s\tMXE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (chr,up_s,dn_e,strand,id_str,name_str))
            w1.write( "%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.1st;Parent=%s.A\n" % (chr,e1st_s,e1st_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.2st;Parent=%s.B\n" % (chr,e2st_s,e2st_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
    if strand =='-':
        id_str = chr+":"+dn_s+":"+dn_e+":"+strand+"@"+chr+":"+se_s+":"+se_e+":"+strand+"@"+chr+":"+up_s+":"+up_e+":"+strand
        name_str = geneSymbol+"_"+chr+":"+dn_s+":"+dn_e+":"+strand+"@"+chr+":"+se_s+":"+se_e+":"+strand+"@"+chr+":"+up_s+":"+up_e+":"+strand
        if event_type!='MXE':
            w2.write( "%s\n" % (name_str))
            w1.write( "%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (chr,up_s,dn_e,strand,id_str,name_str))
            w1.write( "%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (chr,se_s,se_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (chr,up_s,up_e,strand,id_str,id_str))
        if event_type=='MXE':
            w2.write( "%s\n" % (name_str))
            w1.write( "%s\tMXE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (chr,up_s,dn_e,strand,id_str,name_str))
            w1.write( "%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (chr,up_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.1st;Parent=%s.A\n" % (chr,e1st_s,e1st_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (chr,up_s,up_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (chr,dn_s,dn_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.2st;Parent=%s.B\n" % (chr,e2st_s,e2st_e,strand,id_str,id_str))
            w1.write( "%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (chr,up_s,up_e,strand,id_str,id_str))
    w1.close()

    logging.debug("start preparing sashimi plot setting files");

    try:
      prepareSettingFile(gene_no_str);
      pass;
    except:
      logging.debug("There is an exception in preparing setting file");
      logging.debug("Exception: %s" % sys.exc_info()[0]);
      logging.debug("Detail: %s" % sys.exc_info()[1]);
      sys.exit(-1);

    logging.debug("done making setting file..");
    logging.debug("use gff3 format file to run index_gff..");

    os.system('index_gff --index '+outPath+'/Sashimi_index_'+geneSymbol+'_'+str(events_no)+'/tmp.gff3 '+outPath+'/Sashimi_index_'+geneSymbol+'_'+str(events_no)+'/');

    logging.debug("output plot from events set..");
    os.system('sashimi_plot --plot-event \"'+id_str+'\" '+outPath+'/Sashimi_index_'+geneSymbol+'_'+str(events_no)+'/ '+outPath+'/Sashimi_index_'+geneSymbol+'_'+str(events_no)+'/sashimi_plot_settings.txt --output-dir '+outPath+'/Sashimi_plot');
    new_str = id_str.replace(":","_");
    os.system('mv '+outPath+'/Sashimi_plot/'+id_str+'.pdf '+outPath+'/Sashimi_plot/'+geneSymbol+'_'+new_str+'.pdf');
    logging.debug("*** output plot "+geneSymbol);

  # Close opend file
  fo.close()
  w2.close()

  return;

##### end of drawPlotWithEventsFile ####

def drawPlotWithCoordinate(): ## coordinate is provided
  logging.debug("drawPlotWithCoordinate()");

  tmp_str = events.split(':');
  in_chr = tmp_str[0];
  in_strand = tmp_str[1];
  in_coor_s = tmp_str[2];
  in_coor_e = int(tmp_str[3])+1;
  id_str = in_chr+":"+in_coor_s+":"+str(in_coor_e)+":"+in_strand; # chr2:10101175:10104171:+
  gff3_file = tmp_str[4];
  fo = open(gff3_file,'r');
  #logging.debug("events " + events);
  #logging.debug("gff3_file " + gff3_file);
  w2 = open(outPath + '/Sashimi_index/SE.event.list.txt','w');
  w1 = open(outPath+'/Sashimi_index/tmp.gff3','w');
  sashimiPath = outPath + '/Sashimi_index';

  w1.write( "%s\tensGene\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (in_chr,in_coor_s,in_coor_e,in_strand,id_str,id_str))
  w1.write( "%s\tensGene\tmRNA\t%s\t%s\t.\t%s\t.\tName=ENST00000000000;Parent=%s;ID=ENST00000000000\n" % (in_chr,in_coor_s,in_coor_e,in_strand,id_str))

  events_no = 0;
  for line in fo:
    if line.startswith('#'):
      continue;
    events_no += 1;
    items = line.split("\t")
    chr = items[0];
    if in_chr!=chr:
      continue;
    type = items[2];
    if (type=='mRNA' or type=='exon'):
      coor_s = items[3];
      coor_e = items[4];
      strand = items[6];
      annot_str = items[8];
      annot_items = annot_str.split(";");
      if strand =='+' and in_strand==strand and int(in_coor_s)<=int(coor_s) and int(coor_e)<=int(in_coor_e):
        ENST_Name_str = annot_items[0];
        ENST_Parent_str = annot_items[1];
        ENST_ID_str = annot_items[2].replace("\n","");
        if type=='mRNA':
          logging.debug(line.replace("\n",""));
          w1.write( "%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\t%s;Parent=%s;%s\n" % (chr,type,coor_s,coor_e,strand,ENST_Name_str,id_str,ENST_ID_str))
        if type=='exon':
          logging.debug(line.replace("\n",""));
          w1.write( "%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\t%s;%s;%s\n" % (chr,type,coor_s,coor_e,strand,ENST_Name_str,ENST_Parent_str,ENST_ID_str))
      if strand =='-' and in_strand==strand and int(in_coor_s)<=int(coor_e) and int(coor_s)<=int(in_coor_e):
        ENST_Name_str = annot_items[0];
        ENST_Parent_str = annot_items[1];
        ENST_ID_str = annot_items[2].replace("\n","");
        if type=='mRNA':
          logging.debug(line.replace("\n",""));
          w1.write( "%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\t%s;Parent=%s;%s\n" % (chr,type,coor_s,coor_e,strand,ENST_Name_str,id_str,ENST_ID_str))
        if type=='exon':
          logging.debug("exon: " + line.replace("\n",""));
          w1.write( "%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\t%s;%s;%s\n" % (chr,type,coor_s,coor_e,strand,ENST_Name_str,ENST_Parent_str,ENST_ID_str))
  w1.close()

  logging.debug("start preparing sashimi plot setting files");

  try:
    prepareCoorSettingFile();
    pass;
  except:
    logging.debug("There is an exception in preparing coordinate setting file");
    logging.debug("Exception: %s" % sys.exc_info()[0]);
    logging.debug("Detail: %s" % sys.exc_info()[1]);
    sys.exit(-2);

  logging.debug("done making setting file..");
  logging.debug("use gff3 format file to run index_gff..");
  os.system('index_gff --index '+outPath+'/Sashimi_index/tmp.gff3 '+outPath+'/Sashimi_index/');

  logging.debug("output plot from events set..");
  os.system('sashimi_plot --plot-event \"'+id_str+'\" '+outPath+'/Sashimi_index/ '+outPath+'/Sashimi_index/sashimi_plot_settings.txt --output-dir '+outPath+'/Sashimi_plot');
  new_str = id_str.replace(":","_");
  os.system('mv '+outPath+'/Sashimi_plot/'+id_str+'.pdf '+outPath+'/Sashimi_plot/'+new_str+'.pdf');
  logging.debug("*** output plot "+id_str+".pdf");

  # Close opend file
  fo.close()
  w2.close()

  return;

##### end of drawPlotWithCoordinate ####

######## end of functions ##############



################## actual process ##############



if eventsFile==0: #user input coordinate
    logging.debug("get user input coordinates and convert to gff3 format file..");
    logging.debug("coordinate is provided"+"\t"+events);
    try:
      drawPlotWithCoordinate();
      pass;
    except:
      logging.debug("There is an exception in drawPlotWithCoordinate()");du
      logging.debug("Exception: %s" % sys.exc_info()[0]);
      logging.debug("Detail: %s" % sys.exc_info()[1]);
      sys.exit(-1);
    logging.debug("done drawPlotWithCoordinate()");
else: ## events file is provided
    logging.debug("get AS events from rMATS result and convert to gff3 format file..");
    logging.debug("events file is provided"+"\t"+events);
    try:
      drawPlotWithEventsFile();
      pass;
    except:
      logging.debug("There is an exception in drawPlotWithEventsFile()");
      logging.debug("Exception: %s" % sys.exc_info()[0]);
      logging.debug("Detail: %s" % sys.exc_info()[1]);
      sys.exit(-1);
    logging.debug("done drawPlotWithEventsFile()");


#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
