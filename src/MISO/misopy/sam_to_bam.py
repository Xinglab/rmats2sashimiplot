# -*- mode: python; -*-
##
## Convert SAM to indexed, sorted BAM file with headers
##
import os
import time

def sam_to_bam(sam_filename, output_dir,
               header_ref=None):
    # Convert to BAM
    print "Converting SAM to BAM..."
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    t1 = time.time()
    bam_filename = os.path.join(output_dir,
                                "%s.bam" %(os.path.basename(sam_filename).split(".sam")[0]))
    cmd = "samtools view -Sbh %s " %(sam_filename)
    if header_ref != None:
        cmd += " -t %s" %(header_ref)
    cmd += " > %s" %(bam_filename)
    print "  - Executing: %s" %(cmd)
    os.system(cmd)

    # Sort
    print "Sorting BAM file..."
    sorted_filename = "%s.sorted" %(bam_filename.split(".bam")[0])
    cmd = "samtools sort %s %s" %(bam_filename,
                                  sorted_filename)
    print "  - Executing: %s" %(cmd)
    os.system(cmd)

    # Index
    final_filename = "%s.bam" %(sorted_filename)
    print "Indexing BAM..."
    cmd = "samtools index %s" %(final_filename)
    print "  - Executing: %s" %(cmd)
    os.system(cmd)

    t2 = time.time()
    print "Conversion took %.2f minutes." %((t2 - t1)/60.)

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--convert", dest="convert", nargs=2, default=None,
                      help="Convert given SAM file to indexed, sorted BAM file "
                      "with headers. Takes SAM filename and output directory.")
    parser.add_option("--ref", dest="ref", nargs=1, default=None,
                      help="References file to use to get chromosome lengths.")
    (options, args) = parser.parse_args()

    if options.convert != None:
        ref = None

        if options.ref != None:
            ref = os.path.abspath(os.path.expanduser(options.ref))
            print "Using ref: %s" %(ref)
            
        sam_filename = os.path.abspath(os.path.expanduser(options.convert[0]))
        output_dir = os.path.abspath(os.path.expanduser(options.convert[1]))

        sam_to_bam(sam_filename, output_dir, header_ref=ref)
        
    else:
        print "Need --convert to convert SAM to BAM."

if __name__ == "__main__":
    main()
    

    
