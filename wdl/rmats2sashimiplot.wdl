version 1.0

workflow rmats2sashimiplot {
  input {
    String? l1
    String? l2
    String? event_type
    File? e
    String? c_chrom
    String? c_strand
    Int? c_start
    Int? c_end
    File? c_gff3
    Array[File] s1 = []
    Array[File] s2 = []
    Array[File] b1 = []
    Array[File] b2 = []
    Int? exon_s
    Int? intron_s
    File? group_info
    Int? min_counts
    String? color
    Int? font_size
    Int? fig_height
    Int? fig_width
    Boolean hide_number = false
    Boolean no_text_background = false
    Boolean keep_event_chr_prefix = false
    Boolean remove_event_chr_prefix = false
    Int mem_gb = 4
    Int disk_gb = 10
    Boolean should_use_ssd = false
    String rmats2sashimiplot_docker = "xinglab/rmats2sashimiplot:v3.0.0"
  }
  call rmats2sashimiplot_task {
    input:
    l1 = l1,
    l2 = l2,
    event_type = event_type,
    e = e,
    c_chrom = c_chrom,
    c_strand = c_strand,
    c_start = c_start,
    c_end = c_end,
    c_gff3 = c_gff3,
    s1 = s1,
    s2 = s2,
    b1 = b1,
    b2 = b2,
    exon_s = exon_s,
    intron_s = intron_s,
    group_info = group_info,
    min_counts = min_counts,
    color = color,
    font_size = font_size,
    fig_height = fig_height,
    fig_width = fig_width,
    hide_number = hide_number,
    no_text_background = no_text_background,
    keep_event_chr_prefix = keep_event_chr_prefix,
    remove_event_chr_prefix = remove_event_chr_prefix,
    mem_gb = mem_gb,
    disk_gb = disk_gb,
    should_use_ssd = should_use_ssd,
    rmats2sashimiplot_docker = rmats2sashimiplot_docker
  }

  output {
    Array[File] output_plots = rmats2sashimiplot_task.output_plots
  }

  meta {
    description: "https://github.com/Xinglab/rmats2sashimiplot"
  }
}

task rmats2sashimiplot_task {
  input {
    String? l1
    String? l2
    String? event_type
    File? e
    String? c_chrom
    String? c_strand
    Int? c_start
    Int? c_end
    File? c_gff3
    Array[File] s1 = []
    Array[File] s2 = []
    Array[File] b1 = []
    Array[File] b2 = []
    Int? exon_s
    Int? intron_s
    File? group_info
    Int? min_counts
    String? color
    Int? font_size
    Int? fig_height
    Int? fig_width
    Boolean hide_number = false
    Boolean no_text_background = false
    Boolean keep_event_chr_prefix = false
    Boolean remove_event_chr_prefix = false
    Int mem_gb = 4
    Int disk_gb = 10
    Boolean should_use_ssd = false
    String rmats2sashimiplot_docker = "xinglab/rmats2sashimiplot:v3.0.0"
  }

  String l1_param = if defined(l1) then "--l1" else ""
  String l2_param = if defined(l2) then "--l2" else ""
  String event_type_param = if defined(event_type) then "--event-type" else ""
  String e_param = if defined(e) then "-e" else ""
  Boolean any_c_defined = defined(c_chrom) || defined(c_strand) || defined(c_start) || defined(c_end) || defined(c_gff3)
  String c_param = if any_c_defined then "-c" else ""
  String c_sep = if any_c_defined then ":" else ""

  Int s1_len = length(s1)
  String s1_file_path = "s1.txt"
  String s1_args = if s1_len > 0 then "--s1 " + s1_file_path else ""

  Int s2_len = length(s2)
  String s2_file_path = "s2.txt"
  String s2_args = if s2_len > 0 then "--s2 " + s2_file_path else ""

  Int b1_len = length(b1)
  String b1_file_path = "b1.txt"
  String b1_args = if b1_len > 0 then "--b1 " + b1_file_path else ""

  Int b2_len = length(b2)
  String b2_file_path = "b2.txt"
  String b2_args = if b2_len > 0 then "--b2 " + b2_file_path else ""

  String exon_s_param = if defined(exon_s) then "--exon_s" else ""
  String intron_s_param = if defined(intron_s) then "--intron_s" else ""
  String group_info_param = if defined(group_info) then "--group-info" else ""
  String min_counts_param = if defined(min_counts) then "--min-counts" else ""
  String color_param = if defined(color) then "--color" else ""
  String font_size_param = if defined(font_size) then "--font-size" else ""
  String fig_height_param = if defined(fig_height) then "--fig-height" else ""
  String fig_width_param = if defined(fig_width) then "--fig-width" else ""
  String hide_number_param = if hide_number then "--hide-number" else ""
  String no_text_background_param = if no_text_background then "--no-text-background" else ""
  String keep_event_chr_prefix_param = if keep_event_chr_prefix then "--keep-event-chr-prefix" else ""
  String remove_event_chr_prefix_param = if remove_event_chr_prefix then "--remove-event-chr-prefix" else ""

  String out_dir = "out_dir"

  command <<<
    if [[ ~{s1_len} -ne 0 ]]; then
      S1_LINES=~{write_lines(s1)}
      IS_FIRST_S1=1
      while read LINE || [[ -n "${LINE}" ]]; do
        if [[ "${IS_FIRST_S1}" -eq 1 ]]; then
          echo -n "${LINE}" > ~{s1_file_path}
        else
          echo -n ",${LINE}" >> ~{s1_file_path}
        fi
        IS_FIRST_S1=0
      done < ${S1_LINES}
    fi

    if [[ ~{s2_len} -ne 0 ]]; then
      S2_LINES=~{write_lines(s2)}
      IS_FIRST_S2=1
      while read LINE || [[ -n "${LINE}" ]]; do
        if [[ "${IS_FIRST_S2}" -eq 1 ]]; then
          echo -n "${LINE}" > ~{s2_file_path}
        else
          echo -n ",${LINE}" >> ~{s2_file_path}
        fi
        IS_FIRST_S2=0
      done < ${S2_LINES}
    fi

    if [[ ~{b1_len} -ne 0 ]]; then
      B1_LINES=~{write_lines(b1)}
      IS_FIRST_B1=1
      while read LINE || [[ -n "${LINE}" ]]; do
        if [[ "${IS_FIRST_B1}" -eq 1 ]]; then
          echo -n "${LINE}" > ~{b1_file_path}
        else
          echo -n ",${LINE}" >> ~{b1_file_path}
        fi
        IS_FIRST_B1=0
      done < ${B1_LINES}
    fi

    if [[ ~{b2_len} -ne 0 ]]; then
      B2_LINES=~{write_lines(b2)}
      IS_FIRST_B2=1
      while read LINE || [[ -n "${LINE}" ]]; do
        if [[ "${IS_FIRST_B2}" -eq 1 ]]; then
          echo -n "${LINE}" > ~{b2_file_path}
        else
          echo -n ",${LINE}" >> ~{b2_file_path}
        fi
        IS_FIRST_B2=0
      done < ${B2_LINES}
    fi

    rmats2sashimiplot \
      -o ~{out_dir} \
      ~{l1_param} ~{l1} \
      ~{l2_param} ~{l2} \
      ~{event_type_param} ~{event_type} \
      ~{e_param} ~{e} \
      ~{c_param} ~{c_chrom}~{c_sep}~{c_strand}~{c_sep}~{c_start}~{c_sep}~{c_end}~{c_sep}~{c_gff3} \
      ~{s1_args} \
      ~{s2_args} \
      ~{b1_args} \
      ~{b2_args} \
      ~{exon_s_param} ~{exon_s} \
      ~{intron_s_param} ~{intron_s} \
      ~{group_info_param} ~{group_info} \
      ~{min_counts_param} ~{min_counts} \
      ~{color_param} ~{color} \
      ~{font_size_param} ~{font_size} \
      ~{fig_height_param} ~{fig_height} \
      ~{fig_width_param} ~{fig_width} \
      ~{hide_number_param} \
      ~{no_text_background_param} \
      ~{keep_event_chr_prefix_param} \
      ~{remove_event_chr_prefix_param}
  >>>

  runtime {
    docker: rmats2sashimiplot_docker
    cpu: 1
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_gb + if should_use_ssd then " SSD" else " HDD"
  }

  output {
    Array[File] output_plots = glob(out_dir + "/Sashimi_plot/*")
  }

  parameter_meta {
    l1: "The label for the first sample."
    l2: "The label for the second sample."
    event_type: "{SE,A5SS,A3SS,MXE,RI} Type of event from rMATS result used in the analysis. 'SE': skipped exon, 'A5SS': alternative 5' splice site, 'A3SS' alternative 3' splice site, 'MXE': mutually exclusive exons, 'RI': retained intron. (Only if using rMATS event input)"
    e: "The rMATS output event file (Only if using rMATS event input)"
    c_chrom: "The genome region chromosome (Only if using Coordinate and annotation input)"
    c_strand: "The genome region strand (Only if using Coordinate and annotation input)"
    c_start: "The genome region start coordinate (Only if using Coordinate and annotation input)"
    c_end: "The genome region end coordinate (Only if using Coordinate and annotation input)"
    c_gff3: "The GFF3 (not GTF) annotation file of genes and transcripts (Only if using Coordinate and annotation input)"
    s1: "sample_1 sam files: s1_rep1.sam[,s1_rep2.sam]"
    s2: "sample_2 sam files: s2_rep1.sam[,s2_rep2.sam]"
    b1: "sample_1 bam files: s1_rep1.bam[,s1_rep2.bam]"
    b2: "sample_2 bam files: s2_rep1.bam[,s2_rep2.bam]"
    exon_s: "How much to scale down exons. Default: 1"
    intron_s: "How much to scale down introns. For example, --intron_s 5 results in an intron with real length of 100 being plotted as 100/5 = 20. Default: 1"
    group_info: "The path to a *.gf file which groups the replicates. One sashimi plot will be generated for each group instead of the default behavior of one plot per replicate"
    min_counts:  "Individual junctions with read count below --min-counts will be omitted from the plot. Default: 0"
    color: "Specify a list of colors with one color per plot. Without grouping there is one plot per replicate. With grouping there is one plot per group: --color '#CC0011[,#FF8800]'"
    font_size: "Set the font size. Default: 8"
    fig_height: "Set the output figure height (in inches). Default is 7 if sample size < 5 and 14 if sample size is 5 or more"
    fig_width: "Set the output figure width (in inches). Default: 8"
    hide_number: "Do not display the read count on the junctions"
    no_text_background: "Do not put a white box behind the junction read count"
    keep_event_chr_prefix: "force the contig name in the provided events file to be used"
    remove_event_chr_prefix: "remove any leading 'chr' from contig names in the provided events file"
    mem_gb: "GB of RAM"
    disk_gb: "GB of disk"
    should_use_ssd: "Set to use SSD instead of HDD"
    rmats2sashimiplot_docker: "The docker image"
  }
}
