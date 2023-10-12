from __future__ import print_function

import os
import sys
import argparse

import pysam


def convert_sams2bams(sams):
    separator = ","
    split_sams = sams.split(separator)
    bams = list()
    for sam in split_sams:
        bam = sam.replace(".sam", ".bam")
        bams.append(bam)
        os.system("samtools view -Sbh {} > {}".format(sam, bam))

    return separator.join(bams)


def index_bams(bams):
    split_bams = bams.split(",")
    for bam in split_bams:
        index_files = [bam + '.bai', bam + '.csi']
        for index_file in index_files:
            if os.path.isfile(index_file):  # if the bam file has been indexed.
                print("'{}' is indexed already: '{}'".format(bam, index_file))
                return

        print("Indexing '{}'.".format(bam))
        os.system("samtools index " + bam)


def prepare_bams(options):
    """
    convert sam files to bam files and store the filename in options.b1 & options.b2.
    Ensure bams are indexed.
    """
    if options.s1:
        options.b1 = convert_sams2bams(options.s1)
    if options.s2:
        options.b2 = convert_sams2bams(options.s2)

    if options.b1:
        index_bams(options.b1)
    if options.b2:
        index_bams(options.b2)


def file_check(string, expected_ext):
    """
    check the existence of the files and whether they are with the right extensions

    :param string: original string jointed file names by comma
    :param extension: a string like ".bam" in lowercase
    :return: error message or None
    """
    name_arr = string.split(',')
    for name in name_arr:
        if not os.path.isfile(name):
            return '{} is not a file'.format(name)

        extension = os.path.splitext(name)[1]
        if extension.lower() != expected_ext.lower():
            return '{} has extension {} but expected {} (ignoring case)'.format(
                name, extension, expected_ext)

    return None


def maybe_read_paths_from_file(file_path):
    extension = os.path.splitext(file_path)[1]
    if extension.lower() in ['.sam', '.bam']:
        return file_path

    if not os.path.isfile(file_path):
        return file_path

    with open(file_path) as handle:
        first_line = handle.readline()

    return first_line.strip()


def checkout(parser, options):
    """
    check out the required arguments

    :return: None
    """
    # bam files and sam files are alternative, the same for the case of events_file and coordinate
    # events_file should be provided together with event_type
    if ((options.s1 is None and options.b1 is None
         and options.s2 is None and options.b2 is None)):
        parser.error("Not enough arguments! Please provide at least one of"
                     " --s1, --b1, --s2, --b2")
    if (((options.s1 is not None or options.s2 is not None)
         and (options.b1 is not None or options.b2 is not None))):
        parser.error("Specify either sam files or bam files not both")
    if (options.events_file is None or options.event_type is None) and options.coordinate is None:
        parser.error("Not enough arguments! Please provide "
                     "1) coordinates with gff3 file. or "
                     "2) events file together with event type.")

    used_sample_1 = False
    used_sample_2 = False
    if options.s1 is not None:
        used_sample_1 = True
        options.s1 = maybe_read_paths_from_file(options.s1)
        file_check_error = file_check(options.s1, ".sam")
        if file_check_error:
            parser.error("Error checking sam files given as --s1: {}".format(
                file_check_error))

    if options.s2 is not None:
        used_sample_2 = True
        options.s2 = maybe_read_paths_from_file(options.s2)
        file_check_error = file_check(options.s2, ".sam")
        if file_check_error:
            parser.error("Error checking sam files given as --s2: {}".format(
                file_check_error))

    if options.b1 is not None:
        used_sample_1 = True
        options.b1 = maybe_read_paths_from_file(options.b1)
        file_check_error = file_check(options.b1, ".bam")
        if file_check_error:
            parser.error("Error checking bam files given as --b1: {}".format(
                file_check_error))

    if options.b2 is not None:
        used_sample_2 = True
        options.b2 = maybe_read_paths_from_file(options.b2)
        file_check_error = file_check(options.b2, ".bam")
        if file_check_error:
            parser.error("Error checking bam files given as --b2: {}".format(
                file_check_error))

    if options.l1 is None:
        if used_sample_1:
            parser.error('Must provide --l1 if using --s1 or --b1')
        else:
            options.l1 = 'DefaultLabel1'

    if options.l2 is None:
        if used_sample_2:
            parser.error('Must provide --l2 if using --s2 or --b2')
        else:
            options.l2 = 'DefaultLabel2'

    if options.events_file:
        file_check_error = file_check(options.events_file, ".txt")
        if file_check_error:
            parser.error("Error checking rMATS output given as -e: {}".format(
                file_check_error))


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
    sample_1 = list()
    sample_2 = list()
    # sam files have already been converted into bam files and stored in options.b1&b2
    if options.b1:
        sample_1 = options.b1.split(',')
    if options.b2:
        sample_2 = options.b2.split(',')

    for s in sample_1:  # sample1
        bam_files_arr1.append('\"' + s + '\"')
    for s in sample_2:  # sample2
        bam_files_arr2.append('\"' + s + '\"')

    setting_bam_str = ','.join(bam_files_arr1 + bam_files_arr2)

    len_sample1 = len(sample_1)
    len_sample2 = len(sample_2)

    setting_file.write("bam_files = [{0}]\n".format(setting_bam_str))
    setting_file.write("miso_files = [{0}]\n".format(setting_bam_str))
    setting_file.write("[plotting]\n")
    # use a dict to store the configuration
    setting = {}
    if options.fig_height is None:
        if len_sample1 < 5:
            setting['fig_height'] = 7
        else:
            setting["fig_height"] = 14
    else:
        setting["fig_height"] = options.fig_height
    setting["fig_width"] = options.fig_width
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
        setting_color_str = ','.join(colors_arr1 + colors_arr2)
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
            print("Warning: The inclusion levels of Event '{}' contains"
                  " 'NA' value, which could lead to unexpected output."
                  .format(id_str), file=sys.stderr)
    setting_label_str = ','.join(sample_labels_arr1 + sample_labels_arr2)
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


def get_python_executable():
    # Try to get the absolute path of the executable for the running
    # Python interpreter.
    python_executable = sys.executable
    if not python_executable:
        # Fallback
        print('Absolute path for current Python interpreter not found.'
              ' Using "python" without a full path to run scripts',
              file=sys.stderr)
        python_executable = 'python'

    return python_executable


def plot_c(options, id_str):
    """
    the plot part of the coordinate method
    """
    path_index_gff = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                  'MISO/misopy/index_gff.py')
    path_sashimi_plot = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                     'MISO/misopy/sashimi_plot/sashimi_plot.py')

    python_executable = get_python_executable()
    # call python index_gff.py
    tmp_str = os.path.join(options.sashimi_path, "tmp.gff3")
    os.system("{} {} --index {} {}".format(python_executable, path_index_gff,
                                           tmp_str, options.sashimi_path))

    # call python sashimi_plot.py
    setting_str = os.path.join(options.sashimi_path, "sashimi_plot_settings.txt")
    output_path = os.path.join(options.out_dir, "Sashimi_plot")
    if options.group_info is not None:
        os.system("{} {} --plot-event \"{}\" {} {} "
                  "--output-dir {} --group-info {}".format(
                      python_executable, path_sashimi_plot, id_str,
                      options.sashimi_path, setting_str, output_path,
                      options.group_info))
    else:
        os.system("{} {} --plot-event \"{}\" {} {} "
                  "--output-dir {}".format(
                      python_executable, path_sashimi_plot, id_str,
                      options.sashimi_path, setting_str, output_path))

    return


def plot_e(options, id_str, gene_symbol, events_no):
    """
    the plot part of the events file method
    """
    path_index_gff = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                  'MISO/misopy/index_gff.py')
    path_sashimi_plot = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                     'MISO/misopy/sashimi_plot/sashimi_plot.py')

    python_executable = get_python_executable()
    # call python index_gff.py
    out_index = os.path.join(options.out_dir, "Sashimi_index_" + gene_symbol + '_' + str(events_no))
    tmp_str = os.path.join(out_index, "tmp.gff3")
    os.system("{} {} --index {} {}".format(python_executable, path_index_gff,
                                           tmp_str, out_index))

    # call python sashimi_plot.py
    setting_str = os.path.join(out_index, "sashimi_plot_settings.txt")
    output_path = os.path.join(options.out_dir, "Sashimi_plot")
    print("{} {} --plot-event \"{}\" {} {} "
          "--output-dir {}".format(python_executable, path_sashimi_plot, id_str,
                                   out_index, setting_str, output_path))
    if options.group_info is not None:
        os.system("{} {} --plot-event \"{}\" {} {} "
                  "--output-dir {} --group-info {}".format(
                      python_executable, path_sashimi_plot, id_str, out_index,
                      setting_str, output_path, options.group_info))
    else:
        os.system("{} {} --plot-event \"{}\" {} {} "
                  "--output-dir {}".format(
                      python_executable, path_sashimi_plot, id_str, out_index,
                      setting_str, output_path))

    # move pdf file
    old_file = os.path.join(options.out_dir, "Sashimi_plot", id_str + '.pdf')
    new_file = os.path.join(options.out_dir, "Sashimi_plot",
                            str(events_no) + '_' + gene_symbol + '_' + id_str + '.pdf')
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
        id_str = in_chr + "_" + in_coor_s + "_" + str(in_coor_e) + "_" + in_strand  # chr2_10101175_10104171_+
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
            is_mrna_or_transcript = item_type in ["mRNA", "transcript"]
            if is_mrna_or_transcript or item_type == "exon":
                coor_s = items[3]
                coor_e = items[4]
                strand = items[6]
                annot_str = items[8].strip()
                # judge whether the coordinates fit in the item
                if (in_strand == strand
                    and ((item_type == 'exon'
                          and int(in_coor_s) <= int(coor_s)
                          and int(coor_e) <= int(in_coor_e))
                         or (is_mrna_or_transcript
                             and int(coor_s) < int(in_coor_e)
                             and int(coor_e) > int(in_coor_s)))):
                    if is_mrna_or_transcript:
                        if int(coor_s) < int(in_coor_s):
                            coor_s = in_coor_s
                        if int(coor_e) > int(in_coor_e):
                            coor_e = in_coor_e

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
            raise

        plot_c(options, id_str)
        fo.close()

    except Exception as e:
        print(e)
        print("There is an exception in plot_with_coordinate")
        raise

    return


class EventCoor(object):
    """
    to store the coordinates regarding to the event_type
    """

    def __init__(self, event_type, items):
        if event_type == "MXE":
            self.e1st_s = str(int(items[5]) + 1)
            self.e1st_e = items[6]
            self.e2nd_s = str(int(items[7]) + 1)
            self.e2nd_e = items[8]
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
            self.id_str = (seq_chr + "_" + self.up_s + "_" + self.up_e + "_" + strand + "@" +
                           seq_chr + "_" + self.e1st_s + "_" + self.e1st_e + "_" + strand + "@" +
                           seq_chr + "_" + self.e2nd_s + "_" + self.e2nd_e + "_" + strand + "@" +
                           seq_chr + "_" + self.dn_s + "_" + self.dn_e + "_" + strand)
        elif event_type == 'A5SS':
            self.id_str = (seq_chr + "_" + self.sh_s + "_" + self.sh_e + "_" + strand + "@" +
                           seq_chr + "_" + self.lo_s + "_" + self.lo_e + "_" + strand + "@" +
                           seq_chr + "_" + self.fl_s + "_" + self.fl_e + "_" + strand)
        elif event_type == 'A3SS':
            self.id_str = (seq_chr + "_" + self.fl_s + "_" + self.fl_e + "_" + strand + "@" +
                           seq_chr + "_" + self.lo_s + "_" + self.lo_e + "_" + strand + "@" +
                           seq_chr + "_" + self.sh_s + "_" + self.sh_e + "_" + strand)
        elif event_type == 'SE' or event_type == 'RI':
            self.id_str = (seq_chr + "_" + self.up_s + "_" + self.up_e + "_" + strand + "@" +
                           seq_chr + "_" + self.se_s + "_" + self.se_e + "_" + strand + "@" +
                           seq_chr + "_" + self.dn_s + "_" + self.dn_e + "_" + strand)

        self.name_str = gene_symbol + "_" + self.id_str

    def generate_in_reversed_order(self, seq_chr, gene_symbol, strand, event_type):
        if event_type == 'MXE':
            self.id_str = (seq_chr + "_" + self.dn_s + "_" + self.dn_e + "_" + strand + "@" +
                           seq_chr + "_" + self.e2nd_s + "_" + self.e2nd_e + "_" + strand + "@" +
                           seq_chr + "_" + self.e1st_s + "_" + self.e1st_e + "_" + strand + "@" +
                           seq_chr + "_" + self.up_s + "_" + self.up_e + "_" + strand)
        elif event_type == 'A3SS':  # the same as the positive order in A5SS
            self.id_str = (seq_chr + "_" + self.sh_s + "_" + self.sh_e + "_" + strand + "@" +
                           seq_chr + "_" + self.lo_s + "_" + self.lo_e + "_" + strand + "@" +
                           seq_chr + "_" + self.fl_s + "_" + self.fl_e + "_" + strand)
        elif event_type == 'A5SS':  # the same as the positive order in A3SS
            self.id_str = (seq_chr + "_" + self.fl_s + "_" + self.fl_e + "_" + strand + "@" +
                           seq_chr + "_" + self.lo_s + "_" + self.lo_e + "_" + strand + "@" +
                           seq_chr + "_" + self.sh_s + "_" + self.sh_e + "_" + strand)
        elif event_type == 'SE' or event_type == 'RI':
            self.id_str = (seq_chr + "_" + self.dn_s + "_" + self.dn_e + "_" + strand + "@" +
                           seq_chr + "_" + self.se_s + "_" + self.se_e + "_" + strand + "@" +
                           seq_chr + "_" + self.up_s + "_" + self.up_e + "_" + strand)
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
    if options.b1:
        first_bam_path = options.b1.split(',')[0]
    else:
        first_bam_path = options.b2.split(',')[0]

    with pysam.AlignmentFile(first_bam_path, 'rb') as bam_file:
        for alignment in bam_file.fetch(until_eof=True):
            ref_name = alignment.reference_name
            sam_has_chr_prefix = ref_name.startswith('chr')
            break  # only check first alignment

    remove_chr_prefix = (not options.keep_event_chr_prefix
                         and (options.remove_event_chr_prefix
                              or not sam_has_chr_prefix))
    with open(orig_events_file_path, 'rt') as orig_handle:
        with open(new_events_file_path, 'wt') as new_handle:
            for i, line in enumerate(orig_handle):
                is_header_line = i == 0
                columns = line.rstrip('\n').split('\t')

                if is_header_line:
                    chr_index = columns.index('chr')
                elif remove_chr_prefix:
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
                if options.event_type == "SE":
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
                elif options.event_type == "RI":
                    # [se_s, se_e] is the whole retained intron isoform (not just the intron)
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.in;Parent=%s.A\n" % (
                        seq_chr, coor.se_s, coor.se_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                elif options.event_type == "A3SS":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                elif options.event_type == "A5SS":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
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
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.2nd;Parent=%s.B\n" % (
                        seq_chr, coor.e2nd_s, coor.e2nd_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))

            elif strand == '-':
                coor.generate_in_reversed_order(seq_chr, gene_symbol, strand, options.event_type)
                w2.write("%s\n" % coor.name_str)
                if options.event_type == "SE":
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
                elif options.event_type == "RI":
                    # [se_s, se_e] is the whole retained intron isoform (not just the intron)
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.in;Parent=%s.A\n" % (
                        seq_chr, coor.se_s, coor.se_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                elif options.event_type == "A5SS":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.fl_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                elif options.event_type == "A3SS":
                    w1.write("%s\tSE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.sh_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.lo_s, coor.lo_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.sh_s, coor.sh_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tSE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.fl_s, coor.fl_e, strand, coor.id_str, coor.id_str))
                elif options.event_type == "MXE":
                    # On the negative strand, the inclusion isoform includes the 2nd exon
                    w1.write("%s\tMXE\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.name_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.A;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s.B;Parent=%s\n" % (
                        seq_chr, coor.up_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.up;Parent=%s.A\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.2nd;Parent=%s.A\n" % (
                        seq_chr, coor.e2nd_s, coor.e2nd_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.A.dn;Parent=%s.A\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.up;Parent=%s.B\n" % (
                        seq_chr, coor.up_s, coor.up_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.1st;Parent=%s.B\n" % (
                        seq_chr, coor.e1st_s, coor.e1st_e, strand, coor.id_str, coor.id_str))
                    w1.write("%s\tMXE\texon\t%s\t%s\t.\t%s\t.\tID=%s.B.dn;Parent=%s.B\n" % (
                        seq_chr, coor.dn_s, coor.dn_e, strand, coor.id_str, coor.id_str))

            w1.close()
            try:
                conf_setting_file(options, gene_no_str, gene_symbol, events_name_level, coor.id_str)
            except Exception as e:
                print(e)
                print("There is an exception in preparing coordinate setting file")
                raise

            plot_e(options, coor.id_str, gene_symbol, events_no)
        fo.close()
        w2.close()

    except Exception as e:
        print(e)
        print("There is an exception in plot_with_eventsfile")
        raise


def main():
    parser = argparse.ArgumentParser(prog="rmats2sashimiplot")

    required_group = parser.add_argument_group('Required')
    required_group.add_argument("-o", dest="out_dir", required=True,
                                help="The output directory.")

    label_group = parser.add_argument_group('Labels')
    label_group.add_argument("--l1", dest="l1",
                             help="The label for the first sample.")
    label_group.add_argument("--l2", dest="l2",
                             help="The label for the second sample.")

    rmats_group_str = 'rMATS event input'
    coord_group_str = 'Coordinate and annotation input'
    rmats_group = parser.add_argument_group(
        rmats_group_str,
        'Use either ({}) or ({})'.format(rmats_group_str, coord_group_str))
    coordinate_group = parser.add_argument_group(
        coord_group_str,
        'Use either ({}) or ({})'.format(coord_group_str, rmats_group_str))

    rmats_group.add_argument(
        "--event-type", dest="event_type", choices=['SE', 'A5SS', 'A3SS', 'MXE', 'RI'],
        help=("Type of event from rMATS result used in the analysis."
              " 'SE': skipped exon,"
              " 'A5SS': alternative 5' splice site,"
              " 'A3SS' alternative 3' splice site,"
              " 'MXE': mutually exclusive exons,"
              " 'RI': retained intron."
              " (Only if using " + rmats_group_str + ")"))
    rmats_group.add_argument(
        "-e", dest="events_file",
        help=("The rMATS output event file (Only if using "
              + rmats_group_str + ")"))

    coordinate_group.add_argument(
        "-c", dest="coordinate",
        help=("The genome region coordinates and a GFF3 (not GTF) annotation"
              " file of genes and transcripts. The format is"
              " -c {chromosome}:{strand}:{start}:{end}:{/path/to/gff3}"
              " (Only if using " + coord_group_str + ")"))

    sam_bam_group_str_template = '{} Files'
    sam_bam_group_desc_template = (
        'Mapping results for sample_1 & sample_2 in {0} format.'
        ' Replicates must be in a comma separated list.'
        ' A path to a file containing the comma separated list can also be given.'
        ' (Only if using {0})')
    sam_bam_sample_arg_desc_template = (
        'sample_{num} {kind} files: s{num}_rep1.{kind}[,s{num}_rep2.{kind}]')
    group_sam = parser.add_argument_group(
        sam_bam_group_str_template.format('SAM'),
        sam_bam_group_desc_template.format('SAM'))
    group_sam.add_argument(
        "--s1", dest="s1",
        help=sam_bam_sample_arg_desc_template.format(num=1, kind='sam'))
    group_sam.add_argument(
        "--s2", dest="s2",
        help=sam_bam_sample_arg_desc_template.format(num=2, kind='sam'))

    group_bam = parser.add_argument_group(
        sam_bam_group_str_template.format('BAM'),
        sam_bam_group_desc_template.format('BAM'))
    group_bam.add_argument(
        "--b1", dest="b1",
        help=sam_bam_sample_arg_desc_template.format(num=1, kind='bam'))
    group_bam.add_argument(
        "--b2", dest="b2",
        help=sam_bam_sample_arg_desc_template.format(num=2, kind='bam'))

    optional_group = parser.add_argument_group('Optional')
    optional_group.add_argument(
        "--exon_s", dest="exon_s", type=int, default=1,
        help="How much to scale down exons. Default: %(default)s")
    optional_group.add_argument(
        "--intron_s", dest="intron_s", type=int, default=1,
        help=("How much to scale down introns. For example, --intron_s 5"
              " results in an intron with real length of 100 being plotted as"
              " 100/5 = 20. Default: %(default)s"))
    optional_group.add_argument(
        "--group-info", dest="group_info",
        help=('The path to a *.gf file which groups the replicates. One'
              ' sashimi plot will be generated for each group instead of'
              ' the default behavior of one plot per replicate'))
    optional_group.add_argument(
        "--min-counts", dest="min_counts", default=0, type=int,
        help=("Individual junctions with read count below --min-counts will"
              " be omitted from the plot. Default: %(default)s"))
    optional_group.add_argument(
        "--color", dest="color",
        help=('Specify a list of colors with one color per plot. Without'
              ' grouping there is one plot per replicate. With grouping there'
              ' is one plot per group: --color \'#CC0011[,#FF8800]\''))
    optional_group.add_argument(
        "--font-size", dest="font_size", default=8,
        help="Set the font size. Default: %(default)s")
    optional_group.add_argument(
        "--fig-height", dest="fig_height",
        help=('Set the output figure height (in inches). Default is 7'
              ' if sample size < 5 and 14 if sample size is 5 or more'))
    optional_group.add_argument(
        "--fig-width", dest="fig_width", default=8,
        help="Set the output figure width (in inches). Default: %(default)s")
    optional_group.add_argument(
        "--hide-number", dest="hide_number", action="store_true",
        help='Do not display the read count on the junctions')
    optional_group.add_argument(
        "--no-text-background", dest="text_background", action="store_false",
        help='Do not put a white box behind the junction read count')
    optional_group.add_argument(
        "--keep-event-chr-prefix", action="store_true",
        help='force the contig name in the provided events file to be used')
    optional_group.add_argument(
        "--remove-event-chr-prefix", action="store_true",
        help=('remove any leading "chr" from contig names in the provided'
              ' events file'))

    options = parser.parse_args()
    out_path = os.path.abspath(os.path.expanduser(options.out_dir))
    checkout(parser, options)  # 0.check out the arguments

    sashimi_path = os.path.join(out_path, "Sashimi_index")
    if not os.path.isdir(sashimi_path):
        os.makedirs(sashimi_path)
    options.out_dir = out_path
    options.sashimi_path = sashimi_path

    prepare_bams(options)  # 1.convert sam to bam format

    if options.events_file is None:  # 2.setting and plot
        plot_with_coordinate(options)
    else:
        plot_with_eventsfile(options)


if __name__ == '__main__':
    main()
