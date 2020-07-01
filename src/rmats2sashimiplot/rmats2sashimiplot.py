import os
import sys
import argparse
import subprocess


def convert_sam2bam(options):
    """
    convert sam files to bam files and store the filename in options.b1 & options.b2
    """
    if options.s1:  # input with sam file
        sample_1 = options.s1.split(",")
        for sam in sample_1:
            os.system("samtools view -Sbh " + sam + " > " + sam.replace(".sam", ".bam"))
        options.b1 = options.s1.replace(".sam", ".bam")
        sample_2 = options.s2.split(",")
        for sam in sample_2:
            os.system("samtools view -Sbh " + sam + " > " + sam.replace(".sam", ".bam"))
        options.b2 = options.s2.replace(".sam", ".bam")
    print('\033[0;33;m')  # change the print color as yellow
    for bam in options.b1.split(","):
        if os.path.isfile(bam + '.bai'):  # if the bam file has been indexed.
            print("\'{0}\' seems to be indexed already. Please Check out this index file \'{0}.bai\'".format(bam))
        else:
            print("Indexing \'{0}\'.".format(bam))
            os.system("samtools index " + bam)
    for bam in options.b2.split(","):
        if os.path.isfile(bam + '.bai'):
            print("\'{0}\' seems to be indexed already. Please Check out this index file \'{0}.bai\'".format(bam))
        else:
            print("Indexing \'{0}\'.".format(bam))
            os.system("samtools index " + bam)
    print('\033[0m')  # set the color as default value
    return


def file_check(string, extension):
    """
    check the existence of the files and whether they are with the right extensions

    :param string: original string jointed file names by comma
    :param extension: a string like ".bam" in lowercase
    :return: True or False
    """
    name_arr = string.split(',')
    for name in name_arr:
        if not os.path.isfile(name) or os.path.splitext(name)[1].lower() != extension:  # in lowercase
            return False
    return True


def checkout(parser, options):
    """
    check out the required arguments

    :return: None
    """
    # bam files and sam files are alternative, the same for the case of events_file and coordinate
    # events_file should be provided together with event_type
    if (options.s1 is None and options.b1 is None) or (options.s2 is None and options.b2 is None):
        parser.error("Not enough arguments! Please provide sam or bam files.")
    if (options.events_file is None or options.event_type is None) and options.coordinate is None:
        parser.error("Not enough arguments! Please provide "
                     "1) coordinates with gff3 files. or "
                     "2) events files together with events type.")

    if options.s1 is not None and options.s2 is not None:  # with sam file
        if not (file_check(options.s1, ".sam") and file_check(options.s2, ".sam")):
            parser.error("Incorrect file type. Need to provide with the right sam files for --s1 and --s2")
    elif options.b1 is not None and options.b2 is not None:  # with bam file
        if not (file_check(options.b1, ".bam") and file_check(options.b2, ".bam")):
            parser.error("Incorrect file type. Need to provide with the right bam files for --b1 and --b2")
    else:  # only with s1 and lack s2 or other similar cases
        parser.error("Incorrect file type. Need to provide with the right sam files or bam files")

    if options.events_file and not file_check(options.events_file, ".txt"):
        parser.error("Incorrect file type. Need to provide rMATS output format txt file for -e")


def conf_setting_file(options, gene_no_str=None, gene_symbol=None, events_name_level=None, id_str=None):
    """
    configure the setting files
    the empty of gene_no_str means plotting with events file, otherwise with coordinates
    """
    if gene_no_str is not None:
        setting_file = open(os.path.join(options.out_dir, "Sashimi_index_" + gene_no_str,
                                         "sashimi_plot_settings.txt"), 'w')
    else:
        setting_file = open(os.path.join(options.sashimi_path, "sashimi_plot_settings.txt"), 'w')
    setting_file.write("[data]\n")
    sam_dir = ""  # since bam path is accessible, we don't need to provide a prefix particularly
    setting_file.write("bam_prefix = " + sam_dir + "\n")
    setting_file.write("miso_prefix = " + sam_dir + "\n")
    # setting string for replicates
    bam_files_arr1 = []
    bam_files_arr2 = []
    sample_1 = options.b1.split(',')  # sam files has already been converted into bam files and stored in options.b1&b2
    sample_2 = options.b2.split(',')
    for s in sample_1:  # sample1
        bam_files_arr1.append('\"' + s + '\"')
    for s in sample_2:  # sample2
        bam_files_arr2.append('\"' + s + '\"')
    setting_bam_str = ','.join(bam_files_arr1) + ',' + ','.join(bam_files_arr2)
    len_sample1 = len(sample_1)
    len_sample2 = len(sample_2)

    setting_file.write("bam_files = [{0}]\n".format(setting_bam_str))
    setting_file.write("miso_files = [{0}]\n".format(setting_bam_str))
    setting_file.write("[plotting]\n")
    # use a dict to store the configuration
    setting = {}
    if len_sample1 < 5:
        setting['fig_height'] = 7
    else:
        setting["fig_height"] = 14
    setting["fig_width"] = 8
    setting["exon_scale"] = str(options.exon_s)
    setting["intron_scale"] = str(options.intron_s)
    setting["logged"] = False
    setting["font_size"] = options.font_size
    setting["bar_posteriors"] = False
    setting["nyticks"] = 4
    if gene_no_str is None:
        setting["nxticks"] = 11
    else:
        setting["nxticks"] = 6
    setting["show_ylabel"] = True
    setting["show_xlabel"] = True
    setting["plot_title"] = "\"gene symbol\""
    setting["plot_label"] = "plot_label"
    setting["show_posteriors"] = False
    setting["number_junctions"] = not options.hide_number
    setting["resolution"] = ".5"
    setting["reverse_minus"] = True
    setting["min_counts"] = max(options.min_counts, 0)
    setting["text_background"] = options.text_background
    if options.group_info is None:
        setting["group_info"] = False
    else:
        setting["group_info"] = True
    for item in setting:
        setting_file.write("{0} = {1}\n".format(item, setting[item]))

    # setting color
    setting_color_str = ""
    if options.color is None:
        colors_arr1 = ['\"#CC0011\"'] * len_sample1
        colors_arr2 = ['\"#FF8800\"'] * len_sample2
        setting_color_str = ','.join(colors_arr1) + ',' + ','.join(colors_arr2)
    else:
        colors_arr = ["\"{0}\"".format(c) for c in options.color.split(',')]
        setting_color_str = ','.join(colors_arr)
    setting_file.write("colors = [{0}]\n".format(setting_color_str))
    # setting label
    sample_labels_arr1 = []
    sample_labels_arr2 = []
    if gene_no_str is None:  # the case with coordinate
        for rr in range(0, len_sample1):
            sample_labels_arr1.append('\"{0}-{1}\"'.format(options.l1, str(rr + 1)))
        for rr in range(0, len_sample2):
            sample_labels_arr2.append('\"{0}-{1}\"'.format(options.l2, str(rr + 1)))

    else:  # the case with events file
        inc_level_str = events_name_level.get(gene_symbol)
        items = inc_level_str.split('_')
        inc_level1 = items[0]
        inc_level2 = items[1]
        inc_items1 = inc_level1.split(',')
        inc_items2 = inc_level2.split(',')
        warning_flag = False
        for rr in range(0, len_sample1):
            try:
                inc_1 = "{0:.2f}".format(float(inc_items1[rr]))
            except Exception:
                inc_1 = "{0:.2f}".format(float('nan'))
                warning_flag = True
            sample_labels_arr1.append('\"' + gene_symbol + ' ' + options.l1 + '-' + str(rr + 1) + ' IncLevel: '
                                      + inc_1 + '\"')
        for rr in range(0, len_sample2):
            try:
                inc_2 = "{0:.2f}".format(float(inc_items2[rr]))
            except Exception:
                inc_2 = "{0:.2f}".format(float('nan'))
                warning_flag = True
            sample_labels_arr2.append('\"' + gene_symbol + ' ' + options.l2 + '-' + str(rr + 1) + ' IncLevel: '
                                      + inc_2 + '\"')
        if warning_flag:
            print >> sys.stderr, "Warning: The inclusion levels of Event \'{0}\' contains 'NA' value," \
                                 " which could lead to unexpected output.".format(id_str)
    setting_label_str = ','.join(sample_labels_arr1) + ',' + ','.join(sample_labels_arr2)
    setting_file.write("sample_labels = [{0}]\n".format(setting_label_str))

    setting_file.close()


def parse_gff3_record(record):
    """TODO: Docstring for parse_gff3_record.
    :returns: TODO

    """
    res = {}
    eles = record.split('\t')
    res['seqid'] = eles[0]
    res['source'] = eles[1]
    res['type'] = eles[2]
    res['start'] = eles[3]
    res['end'] = eles[4]
    res['score'] = eles[5]
    res['strand'] = eles[6]
    res['phase'] = eles[7]
    res['attributes'] = eles[8]

    return res


def rm_invalid_record(tmp_gff3_fn, id_str):
    """TODO: Docstring for rm_invalid_record.
    :returns: TODO

    """
    eles = id_str.split(':')
    start = int(eles[1])
    end = int(eles[2])
    filtered = []

    with open(tmp_gff3_fn, 'r') as tmp_gff3_fp:
        for record in tmp_gff3_fp:
            res = parse_gff3_record(record)

            if res['type'] == 'gene':
                pass
            elif res['type'] == 'mRNA':
                if int(res['start']) < start:
                    res['start'] = start
                if int(res['end']) > end:
                    res['end'] = end
            elif res['type'] == 'exon':
                if int(res['start']) < start:
                    continue
                if int(res['end']) > end:
                    continue

            filtered.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (res['seqid'], res['source'], res['type'],
                                                                    res['start'], res['end'], res['score'],
                                                                    res['strand'], res['phase'], res['attributes'],))

    with open(tmp_gff3_fn, 'w') as tmp_gff3_fp:
        tmp_gff3_fp.writelines(filtered)

    return


def plot_c(options, id_str):
    """
    the plot part of the coordinate method
    """
    path_index_gff = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                  'MISO/misopy/index_gff.py')
    path_sashimi_plot = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                     'MISO/misopy/sashimi_plot/sashimi_plot.py')

    # call python index_gff.py
    tmp_str = os.path.join(options.sashimi_path, "tmp.gff3")
    os.system("python {0} --index {1} {2}".format(path_index_gff, tmp_str, options.sashimi_path))

    # call python sashimi_plot.py
    setting_str = os.path.join(options.sashimi_path, "sashimi_plot_settings.txt")
    output_path = os.path.join(options.out_dir, "Sashimi_plot")
    if options.group_info is not None:
        os.system("python {0} --plot-event \"{1}\" {2} {3} "
                  "--output-dir {4} --group-info {5}".format(path_sashimi_plot, id_str, options.sashimi_path,
                                                             setting_str, output_path, options.group_info))
    else:
        os.system("python {0} --plot-event \"{1}\" {2} {3} "
                  "--output-dir {4}".format(path_sashimi_plot, id_str, options.sashimi_path, setting_str, output_path))

    # move pdf file
    new_str = id_str.replace(':', '_')
    old_file = os.path.join(options.out_dir, "Sashimi_plot", id_str + '.pdf')
    new_file = os.path.join(options.out_dir, "Sashimi_plot", new_str + '.pdf')
    os.system("mv {0} {1}".format(old_file, new_file))
    return


def plot_e(options, id_str, gene_symbol, events_no):
    """
    the plot part of the events file method
    """
    path_index_gff = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                  'MISO/misopy/index_gff.py')
    path_sashimi_plot = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                     'MISO/misopy/sashimi_plot/sashimi_plot.py')

    # call python index_gff.py
    out_index = os.path.join(options.out_dir, "Sashimi_index_" + gene_symbol + '_' + str(events_no))
    tmp_str = os.path.join(out_index, "tmp.gff3")
    os.system("python {0} --index {1} {2}".format(path_index_gff, tmp_str, out_index))

    # call python sashimi_plot.py
    setting_str = os.path.join(out_index, "sashimi_plot_settings.txt")
    output_path = os.path.join(options.out_dir, "Sashimi_plot")
    print("python {0} --plot-event \"{1}\" {2} {3} "
          "--output-dir {4}".format(path_sashimi_plot, id_str, out_index, setting_str, output_path))
    if options.group_info is not None:
        os.system("python {0} --plot-event \"{1}\" {2} {3} "
                  "--output-dir {4} --group-info {5}".format(path_sashimi_plot, id_str, out_index, setting_str,
                                                             output_path, options.group_info))
    else:
        os.system("python {0} --plot-event \"{1}\" {2} {3} "
                  "--output-dir {4}".format(path_sashimi_plot, id_str, out_index, setting_str, output_path))

    # move pdf file
    new_str = id_str.replace(':', '_')
    old_file = os.path.join(options.out_dir, "Sashimi_plot", id_str + '.pdf')
    new_file = os.path.join(options.out_dir, "Sashimi_plot",
                            str(events_no) + '_' + gene_symbol + '_' + new_str + '.pdf')
    os.system("mv {0} {1}".format(old_file, new_file))
    return


def plot_with_coordinate(options):
    """
    if the user provides with coordinate, then plot in this way
    """
    try:
        tmp_str = options.coordinate.split(':')
        in_chr = tmp_str[0]
        # if not in_chr.startswith("chr"):  # add 'chr' prefix to the sequence name which is from the input arguement
        #     in_chr = "chr" + in_chr
        in_strand = tmp_str[1]
        in_coor_s = tmp_str[2]
        in_coor_e = int(tmp_str[3]) + 1
        id_str = in_chr + ":" + in_coor_s + ":" + str(in_coor_e) + ":" + in_strand  # chr2:10101175:10104171:+
        gff3_file = tmp_str[4]
        fo = open(gff3_file, 'r')
        w2 = open(os.path.join(options.sashimi_path, "SE.event.list.txt"), 'w')
        w1 = open(os.path.join(options.sashimi_path, "tmp.gff3"), 'w')

        w1.write("%s\tensGene\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" %
                 (in_chr, in_coor_s, in_coor_e, in_strand, id_str, id_str))
        # w1.write("%s\tensGene\tmRNA\t%s\t%s\t.\t%s\t.\tName=ENST00000000000;Parent=%s;ID=ENST00000000000\n" %
        #          (in_chr, in_coor_s, in_coor_e, in_strand, id_str))

        events_no = 0
        for line in fo:
            if line.startswith('#'):
                continue
            events_no += 1
            items = line.split("\t")
            # if items[0].startswith("chr"):
            #     item_chr = items[0]
            # else:  # add 'chr' prefix to the seqence name which is from the gff3 file
            #     item_chr = "chr" + items[0]
            item_chr = items[0]
            if in_chr != item_chr:
                continue
            item_type = items[2]
            if item_type == "mRNA" or item_type == "exon":
                coor_s = items[3]
                coor_e = items[4]
                strand = items[6]
                annot_str = items[8].strip()
                # judge whether the coordinates fit in the item
                if (in_strand == strand) and \
                        ((item_type == 'exon' and int(in_coor_s) <= int(coor_s) and int(coor_e) <= int(in_coor_e)) or \
                                 (item_type == 'mRNA' and int(coor_s) < int(in_coor_e) and int(coor_e) > int(
                                     in_coor_s))):
                    if item_type == 'mRNA':
                        if int(coor_s) < in_coor_s:
                            coor_s = in_coor_s
                        if int(coor_e) > in_coor_e:
                            coor_e = in_coor_e
                    if item_type == "mRNA":
                        annot_str = annot_str.replace('Parent', 'Note')
                        w1.write("%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\tParent=%s;%s\n" %
                                 (item_chr, item_type, coor_s, coor_e, strand, id_str, annot_str))
                    if item_type == "exon":
                        w1.write("%s\tensGene\t%s\t%s\t%s\t.\t%s\t.\t%s\n" %
                                 (item_chr, item_type, coor_s, coor_e, strand, annot_str))
        w1.close()

        try:
            conf_setting_file(options)
        except Exception as e:
            print(e)
            print("There is an exception in preparing coordinate setting file")
            sys.exit(-2)
        plot_c(options, id_str)
        fo.close()

    except Exception as e:
        print(e)
        print("There is an exception in plot_with_coordinate")
        sys.exit(-1)
    return


class EventCoor(object):
    """
    to store the coordinates regarding to the event_type
    """

    def __init__(self, event_type, items):
        if event_type == "MXE":
            self.e1st_s = str(int(items[5]) + 1)
            self.e1st_e = items[6]
            self.e2st_s = str(int(items[7]) + 1)
            self.e2st_e = items[8]
            self.up_s = str(int(items[9]) + 1)
            self.up_e = items[10]
            self.dn_s = str(int(items[11]) + 1)
            self.dn_e = items[12]
            self.inc_level1 = items[22]  # IncLevel1
            self.inc_level2 = items[23]  # IncLevel2
        elif event_type == "SE" or event_type == "RI":
            self.se_s = str(int(items[5]) + 1)
            self.se_e = items[6]
            self.up_s = str(int(items[7]) + 1)
            self.up_e = items[8]
            self.dn_s = str(int(items[9]) + 1)
            self.dn_e = items[10]
            self.inc_level1 = items[20]  # IncLevel1
            self.inc_level2 = items[21]  # IncLevel2
        else:  # A3SS or A5SS
            self.lo_s = str(int(items[5]) + 1)  # long
            self.lo_e = items[6]
            self.sh_s = str(int(items[7]) + 1)  # short
            self.sh_e = items[8]
            self.fl_s = str(int(items[9]) + 1)  # flanking
            self.fl_e = items[10]
            self.inc_level1 = items[20]  # IncLevel1
            self.inc_level2 = items[21]  # IncLevel2

        self.name_str = ''
        self.id_str = ''

    def generate_in_positive_order(self, seq_chr, gene_symbol, strand, event_type):
        if event_type == 'MXE':
            self.id_str = (seq_chr + ":" + self.up_s + ":" + self.up_e + ":" + strand + "@" +
                           seq_chr + ":" + self.e1st_s + ":" + self.e1st_e + ":" + strand + "@" +
                           seq_chr + ":" + self.e2st_s + ":" + self.e2st_e + ":" + strand + "@" +
                           seq_chr + ":" + self.dn_s + ":" + self.dn_e + ":" + strand)
        elif event_type == 'A5SS':
            self.id_str = (seq_chr + ":" + self.sh_s + ":" + self.sh_e + ":" + strand + "@" +
                           seq_chr + ":" + self.lo_s + ":" + self.lo_e + ":" + strand + "@" +
                           seq_chr + ":" + self.fl_s + ":" + self.fl_e + ":" + strand)
        elif event_type == 'A3SS':
            self.id_str = (seq_chr + ":" + self.fl_s + ":" + self.fl_e + ":" + strand + "@" +
                           seq_chr + ":" + self.lo_s + ":" + self.lo_e + ":" + strand + "@" +
                           seq_chr + ":" + self.sh_s + ":" + self.sh_e + ":" + strand)
        elif event_type == 'SE' or event_type == 'RI':
            self.id_str = (seq_chr + ":" + self.up_s + ":" + self.up_e + ":" + strand + "@" +
                           seq_chr + ":" + self.se_s + ":" + self.se_e + ":" + strand + "@" +
                           seq_chr + ":" + self.dn_s + ":" + self.dn_e + ":" + strand)

        self.name_str = gene_symbol + "_" + self.id_str

    def generate_in_reversed_order(self, seq_chr, gene_symbol, strand, event_type):
        if event_type == 'MXE':
            self.id_str = (seq_chr + ":" + self.dn_s + ":" + self.dn_e + ":" + strand + "@" +
                           seq_chr + ":" + self.e2st_s + ":" + self.e2st_e + ":" + strand + "@" +
                           seq_chr + ":" + self.e1st_s + ":" + self.e1st_e + ":" + strand + "@" +
                           seq_chr + ":" + self.up_s + ":" + self.up_e + ":" + strand)
        elif event_type == 'A3SS':  # the same as the positive order in A5SS
            self.id_str = (seq_chr + ":" + self.sh_s + ":" + self.sh_e + ":" + strand + "@" +
                           seq_chr + ":" + self.lo_s + ":" + self.lo_e + ":" + strand + "@" +
                           seq_chr + ":" + self.fl_s + ":" + self.fl_e + ":" + strand)
        elif event_type == 'A5SS':  # the same as the positive order in A3SS
            self.id_str = (seq_chr + ":" + self.fl_s + ":" + self.fl_e + ":" + strand + "@" +
                           seq_chr + ":" + self.lo_s + ":" + self.lo_e + ":" + strand + "@" +
                           seq_chr + ":" + self.sh_s + ":" + self.sh_e + ":" + strand)
        elif event_type == 'SE' or event_type == 'RI':
            self.id_str = (seq_chr + ":" + self.dn_s + ":" + self.dn_e + ":" + strand + "@" +
                           seq_chr + ":" + self.se_s + ":" + self.se_e + ":" + strand + "@" +
                           seq_chr + ":" + self.up_s + ":" + self.up_e + ":" + strand)
        self.name_str = gene_symbol + "_" + self.id_str


def create_chr_aware_events_file(options):
    """
    The *.MATS.*.txt events file from rmats includes the prefix 'chr'
    for chromosomes. If the BAM files do not have the 'chr' prefix
    then remove the prefix from the events file. Consistent presence or
    absence of 'chr' is needed to search for reads in the BAM files.
    """
    orig_events_file_path = options.events_file
    new_events_file_path = os.path.join(options.sashimi_path, 'events_file.txt')
    first_bam_path = options.b1.split(',')[0]
    tmp_sam_file_path = os.path.join(options.sashimi_path, 'tmp_chr_check.sam')
    with open(tmp_sam_file_path, 'wt') as tmp_sam_file_handle:
        subprocess.check_call(['samtools', 'view', first_bam_path],
                              stdout=tmp_sam_file_handle)

    with open(tmp_sam_file_path, 'rt') as tmp_sam_file_handle:
        first_line_of_sam = tmp_sam_file_handle.readline().rstrip('\n')

    os.remove(tmp_sam_file_path)

    sam_columns = first_line_of_sam.split('\t')
    sam_chr_column = sam_columns[2]
    sam_has_chr_prefix = sam_chr_column.startswith('chr')

    with open(orig_events_file_path, 'rt') as orig_handle:
        with open(new_events_file_path, 'wt') as new_handle:
            for i, line in enumerate(orig_handle):
                is_header_line = i == 0
                columns = line.rstrip('\n').split('\t')

                if is_header_line:
                    chr_index = columns.index('chr')
                elif not sam_has_chr_prefix:
                    orig_chr_column = columns[chr_index]
                    if orig_chr_column.startswith('chr'):
                        columns[chr_index] = orig_chr_column[3:]

                new_line = '\t'.join(columns)
                new_handle.write('{}\n'.format(new_line))

    return new_events_file_path


def plot_with_eventsfile(options):
    """
    if the user provides with event files, then plot in this way
    """
    try:
        options.events_file = create_chr_aware_events_file(options)
        fo = open(options.events_file, 'r')
        w2 = open(os.path.join(options.sashimi_path, options.event_type + ".event.list.txt"), 'w')
        events_name_level = {}
        events_no = 0
        for line in fo:
            if line.startswith('ID'):
                continue
            events_no += 1
            items = line.split("\t")
            gene_symbol = items[2]
            gene_symbol = gene_symbol.replace("\"", '')
            gene_no_str = gene_symbol + '_' + str(events_no)
            sashimi_path = os.path.join(options.out_dir, "Sashimi_index_" + gene_symbol + '_' + str(events_no))
            if not os.path.isdir(sashimi_path):
                os.makedirs(sashimi_path)
            w1 = open(os.path.join(sashimi_path, "tmp.gff3"), 'w')
            seq_chr = items[3]
            strand = items[4]
            # construct the events coordinates depending on the event type
            coor = EventCoor(options.event_type, items)
            events_name_level[gene_symbol] = coor.inc_level1 + "_" + coor.inc_level2

            if strand == '+':
                coor.generate_in_positive_order(seq_chr, gene_symbol, strand, options.event_type)
                w2.write("%s\n" % coor.name_str)
                if options.event_type == "SE" or options.event_type == "RI":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.se_s, coor.se_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "A3SS":  # flanking -- long/short TODO: change the ID name
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "A5SS":  # short/long -- flanking TODO: change the ID name
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "MXE":
                    w1.write("%s\tMXE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.1st;Parent=%s.A\n" % (
                        seq_chr, coor.e1st_s, coor.e1st_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.2st;Parent=%s.B\n" % (
                        seq_chr, coor.e2st_s, coor.e2st_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))


            elif strand == '-':
                coor.generate_in_reversed_order(seq_chr, gene_symbol, strand, options.event_type)
                w2.write("%s\n" % coor.name_str)
                if options.event_type == "SE" or options.event_type == "RI":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.se_s, coor.se_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "A5SS":  # flanking -- long/short TODO: change the ID name
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "A3SS":  # short/long -- flanking TODO: change the ID name
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.se;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))

                elif options.event_type == "MXE":
                    w1.write("%s\tMXE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.1st;Parent=%s.A\n" % (
                        seq_chr, coor.e1st_s, coor.e1st_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.2st;Parent=%s.B\n" % (
                        seq_chr, coor.e2st_s, coor.e2st_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
            w1.close()
            try:
                conf_setting_file(options, gene_no_str, gene_symbol, events_name_level, coor.id_str)
            except Exception as e:
                print(e)
                print("There is an exception in preparing coordinate setting file")
                sys.exit(-2)
            plot_e(options, coor.id_str, gene_symbol, events_no)
        fo.close()
        w2.close()

    except Exception as e:
        print(e)
        print("There is an exception in plot_with_eventsfile")
        sys.exit(-1)


def main():
    parser = argparse.ArgumentParser(prog="rmats2sashimiplot",
                                     usage="\n"
                                           "Usage(with sam files):\n"
                                           "%(prog)s --s1 s1_rep1.sam[,s1_rep2.sam]* --s2 s2.rep1.sam[,s2.rep2.sam]*"
                                           " -t eventType -e eventsFile --l1 SampleLabel1 --l2 SampleLable2 --exon_s "
                                           "exonScale --intron_s intronScale -o outDir\n\n"
                                           "Example (with sam files):\n"
                                           "%(prog)s --s1 ./testData/S1.R1.test.sam,./testData/S1.R2.test.sam,"
                                           "./testData/S1.R3.test.sam --s2 ./testData/S2.R1.test.sam,./testData/"
                                           "S2.R2.test.sam,./testData/S2.R3.test.sam -t SE -e ./testData/MATS_output/"
                                           "test_PC3E_GS689.SE.MATS.events.txt --l1 PC3E --l2 GS689 --exon_s 1 "
                                           "--intron_s 5 -o test_events_output  \n\n"
                                           "Usage (with bam files):\n"
                                           "%(prog)s --b1 s1_rep1.bam[,s1_rep2.bam]* --b2 s2.rep1.bam[,s2.rep2.bam]* "
                                           "-c coordinate:annotaionFile --l1 SampleLabel1 --l2 SampleLable2 --exon_s "
                                           "exonScale --intron_s intronScale -o outDir  \n\n"
                                           "Example (with bam files):\n"
                                           "%(prog)s --b1 ./testData/S1.R1.test.bam,./testData/S1.R2.test.bam,"
                                           "./testData/S1.R3.test.bam --b2 ./testData/S2.R1.test.bam,./testData/"
                                           "S2.R2.test.bam,./testData/S2.R3.test.bam -c chr2:+:10090000:10110000:"
                                           "./testData/ensGene.gff3 --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 "
                                           "-o test_coordinate_output\n",
                                     )
    parser.add_argument("-t", dest="event_type", choices=['SE', 'A5SS', 'A3SS', 'MXE', 'RI'],
                        help="Type of event from rMATS result used in the analysis."
                             "eventType is \'SE\', \'A5SS\', \'A3SS\', \'MXE\' or \'RI\'."
                             "\'SE\' is for skipped exon events,"
                             "\'A5SS\' is for alternative 5\' splice site events,"
                             "\'A3SS\' is for alternative 3\' splice site events,"
                             "\'MXE\' is for mutually exclusive exons events "
                             "and \'RI\' is for retained intron events "
                             "(Only if using rMATS format result as event file).")
    parser.add_argument("-e", dest="events_file", help="The rMATS output event file "
                                                       "(Only if using rMATS format result as event file).")
    parser.add_argument("-c", dest="coordinate", help="The coordinate of genome region and an annotation "
                                                      "of genes and transcripts in GFF3 format. Coordinate"
                                                      "and annotation file must be colon separated"
                                                      "(Only if using coordinate and annotation file).")
    parser.add_argument("--l1", dest="l1", help="The label for first sample.", required=True)
    parser.add_argument("--l2", dest="l2", help="The label for second sample.", required=True)
    parser.add_argument("-o", dest="out_dir", help="The output directory.", required=True)

    group_sam = parser.add_argument_group("Sam Files", "Mapping results for the sample_1 & sample_2 in sam format"
                                                       "Replicates must be in a comma separated list"
                                                       "(Only if using sam).")
    group_sam.add_argument("--s1", action="store", dest="s1",
                           help="sample_1 in sam format (s1_rep1.sam[,s1_rep2.sam])")
    group_sam.add_argument("--s2", action="store", dest="s2",
                           help="sample_2 in sam format (s2_rep1.sam[,s2_rep2.sam])")

    group_bam = parser.add_argument_group("Bam Files", "Mapping results for the sample_1 & sample_2 in bam format"
                                                       "Replicates must be in a comma separated list"
                                                       "(Only if using bam).")
    group_bam.add_argument("--b1", action="store", dest="b1",
                           help="sample_1 in bam format(s1_rep1.bam[,s1_rep2.bam])")
    group_bam.add_argument("--b2", action="store", dest="b2",
                           help="sample_2 in bam format(s2_rep1.bam[,s2_rep2.bam])")

    group_optional = parser.add_argument_group("Optional Parameters",
                                               "These parameters have their default values.")
    group_optional.add_argument("--exon_s", dest="exon_s", type=int, default=1,
                                help="The size of scale down exons. The default is 1.")
    group_optional.add_argument("--intron_s", dest="intron_s", type=int, default=1,
                                help="The size of scale down introns. For example,"
                                     "if -intron_s is 5, it means the size of intron is 5:1"
                                     "if the real size of intron is 5, the size in the "
                                     "plot will be scaled down to 1). The default is 1.")
    group_optional.add_argument("--group-info", dest="group_info", default=None,
                                help="If the user wants to divide the bam files manually, "
                                     "you can provide a \'*.gf\' file.")
    group_optional.add_argument("--min-counts", dest="min_counts", default=0,
                                help="If the junction count is smaller than this number, this single junction's count "
                                     "would be omitted in the plot.")
    group_optional.add_argument("--color", dest="color", default=None,
                                help="Set the color in format(\"#CC0011\"[,\"#CC0011\"]). "
                                     "The number of the colors equal to the total number of bam files in "
                                     "different samples.")
    group_optional.add_argument("--font-size", dest="font_size", default=8,
                                help="Set the font size.")
    group_optional.add_argument("--hide-number", dest="hide_number", action="store_true")
    group_optional.set_defaults(hide_number=False)
    group_optional.add_argument("--no-text-background", dest="text_background", action="store_false")
    group_optional.set_defaults(text_background=True)

    options = parser.parse_args()

    out_path = os.path.abspath(os.path.expanduser(options.out_dir))

    checkout(parser, options)  # 0.check out the arguments

    sashimi_path = os.path.join(out_path, "Sashimi_index")
    if not os.path.isdir(sashimi_path):
        os.makedirs(sashimi_path)
    options.out_dir = out_path
    options.sashimi_path = sashimi_path

    convert_sam2bam(options)  # 1.convert sam to bam format

    if options.events_file is None:  # 2.setting and plot
        plot_with_coordinate(options)
    else:
        plot_with_eventsfile(options)


if __name__ == '__main__':
    main()

