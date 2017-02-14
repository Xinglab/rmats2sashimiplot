# -*- mode: python; -*-
##
## Pack MISO output into an SQL database
##
import os
import csv
import time
import sys
import subprocess
import shutil
from collections import defaultdict

import pysam

import misopy
import misopy.gff_utils as gff_utils
import misopy.as_events as as_events
import misopy.run_miso as run_miso
import misopy.misc_utils as misc_utils
import misopy.run_events_analysis as run_events
import misopy.miso_db as miso_db
from misopy.settings import Settings, load_settings
from misopy.settings import miso_path as miso_settings_path
import misopy.cluster_utils as cluster_utils

miso_path = os.path.dirname(os.path.abspath(__file__))
manual_url = "http://genes.mit.edu/burgelab/miso/docs/"

class MISOPacker:
    """
    Pack MISO directories. Traverses the directory that contains
    MISO output and converts any *.miso-containing directories
    into MISO database files (*.miso_db).
    """
    def __init__(self, dirs_to_pack):
        self.dirs_to_pack = dirs_to_pack


    def pack_dirs(self, miso_dirnames):
        """
        Takes a set of MISO input containing directories and packs them.

        'miso_dirnames' are directory that have MISO output *somewhere*
        in them -- could be arbitrarily nested within the directory.

        This traverses the directory structure and converts raw *.miso text
        containing directory into *.miso_db files that are much more compact
        and are less of a burden for the filesystem.
        """
        t1 = time.time()
        for miso_dirname in miso_dirnames:
            print "Processing: %s" %(miso_dirname)
            if not os.path.isdir(miso_dirname):
                print "Error: %s not a directory." %(miso_dirname)
                sys.exit(1)
            for dir_to_compress, subdirs, curr_fnames in os.walk(miso_dirname):
                if miso_db.is_miso_unpacked_dir(dir_to_compress):
                    # It's a *.miso containing directory, so pack it
                    # into a MISO database
                    chrom_basename = os.path.basename(dir_to_compress)
                    if len(chrom_basename) == 0:
                        print "Error: Failed to pack MISO directory %s" \
                              %(miso_dirname)
                        raise Exception, "Basename for %s is empty!" \
                              %(dir_to_compress)
                    db_fname = \
                      os.path.join(os.path.dirname(dir_to_compress),
                                   "%s%s" %(chrom_basename,
                                             miso_db.MISO_DB_EXT))
                    # If packed file exists, move on
                    if os.path.isfile(db_fname):
                        continue
                    status = miso_db.miso_dir_to_db(dir_to_compress, db_fname)
                    # If packing was successful, delete the input directory
                    # containing the *.miso file
                    shutil.rmtree(dir_to_compress)
        t2 = time.time()
        print "Packing took %.2f minutes" %((t2 - t1)/60.)
        

def greeting(parser=None):
    print "MISO (Mixture of Isoforms model)"
    print "Pack the MISO output into an SQL database."
    print "Use --help argument to view options.\n"
    print "Example usage:\n"
    print "miso_pack --pack mydir"
    if parser is not None:
        parser.print_help()


def pack_miso_output(dirs_to_pack_as_str):
    """
    Pack MISO output into an SQL database.
    """
    dirs_to_pack = dirs_to_pack_as_str.split(",")
    # Make into absolute paths
    dirs_to_pack = map(misc_utils.pathify, dirs_to_pack)
    miso_packer = MISOPacker(dirs_to_pack)
    miso_packer.pack_dirs(dirs_to_pack)


def view_miso_db(db_fname):
    db_fname = misc_utils.pathify(db_fname)
    if not os.path.isfile(db_fname):
        print "Error: %s does not exist." %(db_fname)
        sys.exit(1)
    curr_db = miso_db.MISODatabase(db_fname)
    event_names = curr_db.get_all_event_names()
    num_events = len(event_names)
    print "Database contains %d events" %(num_events)
    for event in event_names:
        print event


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--pack", dest="pack",
                      nargs=1, default=None,
                      help="Pack a MISO output containing dir(s). Takes as input " \
                      "a directory or a comma-separated set of directories " \
                      "that contain MISO output.")
    parser.add_option("--view", dest="view",
                      nargs=1, default=None,
                      help="View a MISO database (.miso_db file).")
    (options, args) = parser.parse_args()

    if (options.pack is None) and (options.view is None):
        greeting()
        sys.exit(1)

    if options.pack is not None:
        pack_miso_output(options.pack)

    if options.view is not None:
        view_miso_db(options.view)

        
if __name__ == "__main__":
    main()
        
