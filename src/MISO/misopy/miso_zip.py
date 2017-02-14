# -*- mode: python; -*-
##
## Compress/uncompress directories containing MISO output.
##
## Reduces the number of raw output files (*.miso files) by
## storing them in a database, and then compressing using
## zip compression to reduce the space footprint.
##
import os
import sys
import time

import zipfile
import sqlite3
import shutil
import fnmatch
import glob

import misopy
import misopy.misc_utils as utils
import misopy.miso_utils as miso_utils
import misopy.miso_db as miso_db


class MISOCompressor:
    """
    Compressor/uncompressor of MISO output-containing directories.

    The compressor:

    (1) copies the original directory containing MISO directories
    (2) creates MISO databases (.miso_db files)
        from the *.miso-containing directories
    (3) then zip's up the resulting directory

    The compressor works on the copy and leaves the original
    directory unmodified.
    """
    def __init__(self):
        self.input_dir = None
        self.output_dir = None
        # Extension for compressed directory that contains
        # any sort of MISO output within it
        self.comp_ext = ".misozip"
        

    def compress(self, output_filename, miso_dirnames):
        """
        Takes a set of MISO input directories and compresses them
        into 'output_filename'. This involves making SQL databases
        for all the MISO directories and then additionally compressing the
        results as a zip file.
        """
        if os.path.isfile(output_filename):
            print "Error: %s already exists. Please delete to overwrite." \
                  %(output_filename)
        output_dir = "%s%s" %(output_filename, miso_db.MISO_DB_EXT)
        if os.path.isdir(output_dir):
            print "Error: Intermediate compressed directory %s " \
                  "exists. Please delete to overwrite." %(output_dir)
            sys.exit(1)
        for miso_dirname in miso_dirnames:
            print "Processing: %s" %(miso_dirname)
            if not os.path.isdir(miso_dirname):
                print "Error: %s not a directory." %(miso_dirname)
                sys.exit(1)
            if os.path.isfile(output_filename):
                print "Output file %s already exists, aborting. " \
                      "Please delete the file if you want " \
                      "compression to run."
                sys.exit(1)
            self.miso_dirs_to_compress = []
            print "Copying source directory tree.."
            shutil.copytree(miso_dirname, output_dir,
                            ignore=self.collect_miso_dirs)
            for dir_to_compress in self.miso_dirs_to_compress:
                rel_path = os.path.relpath(dir_to_compress, miso_dirname)
                comp_path = os.path.join(output_dir, rel_path)
                # Remove the place holder directory
                os.rmdir(comp_path)
                comp_path = "%s%s" %(comp_path, miso_db.MISO_DB_EXT)
                miso_db.miso_dir_to_db(dir_to_compress, comp_path)
        # Zip directory using conventional zip
        print "Zipping compressed directory with standard zip..."
        t1 = time.time()
        zipper(output_dir, output_filename)
        print "Deleting intermediate directory: %s" %(output_dir)
        shutil.rmtree(output_dir)
        t2 = time.time()
        print "  - Standard zipping took %.2f minutes." \
              %((t2 - t1)/60.)
        print "To access the SQLite representation of raw MISO output "
        print "(*.miso) files, simply unzip with the .miso_zip file "
        print "with standard unzip utility:\n"
        print "  unzip %s" %(output_filename)
        

    def uncompress(self, compressed_filename, output_dir):
        """
        Takes a compressed MISO file 'compressed_filename' and
        uncompresses it into 'output_dir'.
        """
        if not os.path.isfile(compressed_filename):
            print "Error: Cannot find %s, aborting." \
                  %(compressed_filename)
        if not os.path.basename(compressed_filename).endswith(self.comp_ext):
            print "Warning: %s does not end in %s. Are you sure it is " \
                  "a file compressed by miso_zip.py?" \
                  %(compressed_filename, self.comp_ext)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        print "Uncompressing %s into %s" %(compressed_filename,
                                           output_dir)
        # First unzip the file using conventional unzip
        unzipped_files = unzipper(compressed_filename, output_dir)
        # Remove the original .zip
        if os.path.isfile(compressed_filename):
            print "Removing the compressed file %s" %(compressed_filename)
            if os.path.isfile(compressed_filename):
                os.remove(compressed_filename)
                
                
    def collect_miso_dirs(self, path, dirnames):
        """
        Collect raw MISO output directories
        """
        if fnmatch.filter(dirnames, "*.miso"):
            self.miso_dirs_to_compress.append(path)
            return dirnames
        return []



def compress_miso(output_filename, input_dirs,
                  comp_ext=".misozip"):
    """
    Compress a directory containing MISO files.

    Traverse directories, one by one, and look for directories
    that contain
    """
    output_filename = utils.pathify(output_filename)
    for input_dir in input_dirs:
        if not os.path.isdir(input_dir):
            print "Error: Cannot find directory %s" %(input_dir)
            sys.exit(1)
    if not os.path.basename(output_filename).endswith(comp_ext):
        print "Error: Compressed output filename must end in %s" \
              %(comp_ext)
        sys.exit(1)
    if os.path.isfile(output_filename):
        print "Error: Output filename exists. Please delete %s to overwrite." \
              %(output_filename)
        sys.exit(1)
    t1 = time.time()
    miso_comp = MISOCompressor()
    miso_comp.compress(output_filename, input_dirs)
    t2 = time.time()
    print "Compression took %.2f minutes." %((t2 - t1)/60.)

    
def uncompress_miso(compressed_filename, output_dir):
    """
    Uncompress MISO directory.
    """
    if not os.path.isfile(compressed_filename):
        print "Error: Zip file %s is not found." \
              %(compressed_filename)
        sys.exit(1)
    t1 = time.time()
    miso_comp = MISOCompressor()
    miso_comp.uncompress(compressed_filename, output_dir)
    t2 = time.time()
    print "Uncompression took %.2f minutes." %((t2 - t1)/60.)


def zipper(dir, zip_file):
    """
    Zip a directory 'dir' recursively, saving result in
    'zip_file'.
    
    by Corey Goldberg.
    """
    # Enable Zip64 to allow creation of large Zip files
    zip = zipfile.ZipFile(zip_file, 'w',
                          compression=zipfile.ZIP_DEFLATED,
                          allowZip64=True)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


def unzipper(zip_file, outdir):
    """
    Unzip a given 'zip_file' into the output directory 'outdir'.
    
    Return the names of files in the archive.
    """
    zf = zipfile.ZipFile(zip_file, "r",
                         allowZip64=True)
    filenames = zf.namelist()
    zf.extractall(outdir)
    return filenames


def greeting():
    print "Compress/uncompress MISO output. Usage:\n"
    print "To compress a directory containing MISO files \'inputdir\', use: "
    print "   miso_zip --compress outputfile.misozip inputdir"
    print "To uncompress back into a directory \'outputdir\', use: "
    print "   miso_zip --uncompress outputfile.misozip outputdir"
    print "\nNote: compressed filename must end in \'.misozip\'"
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compress", dest="compress", nargs=2, default=None,
                      help="Compress a directory containing MISO output. "
                      "Takes as arguments: (1) the output filename of the "
                      "compressed file, (2) a comma-separated list of "
                      "directory names to be compressed. "
                      "Example: --compress output.misozip dirname1,dirname2")
    parser.add_option("--uncompress", dest="uncompress", nargs=2, default=None,
                      help="Uncompress a file generated by compress_miso. "
                      "Takes as arguments: (1) the filename to be "
                      "uncompressed, and (2) the directory to place the "
                      "uncompressed representation into. "
                      "Example: --uncompress output.misozip outputdir")
    (options, args) = parser.parse_args()

    if (options.compress is None) and (options.uncompress is None):
        greeting()
        sys.exit(1)
    elif (options.compress is not None) and (options.uncompress is not None):
        # Can't be given both.
        greeting()
        print "Error: Cannot process --compress and --uncompress at same time."
        sys.exit(1)

    if options.compress is not None:
        output_filename = utils.pathify(options.compress[0])
        input_dirs = [utils.pathify(d) \
                      for d in options.compress[1].split(",")]
        compress_miso(output_filename, input_dirs)

    if options.uncompress is not None:
        compressed_filename = utils.pathify(options.uncompress[0])
        output_dir = utils.pathify(options.uncompress[1])
        uncompress_miso(compressed_filename, output_dir)


if __name__ == "__main__":
    main()
    
