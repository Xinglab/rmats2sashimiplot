##
## Misc. utilities
##
import os
import sys
import time
import shelve

import time
from time import strftime

COMPRESS_PREFIX = "misocomp"


def get_timestamp():
    """
    Return filename-friendly time stamp.
    """
    tstamp = strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    return tstamp


def inv_dict(mydict):
    """
    Reverse key -> val into val -> key.
    """
    new_dict = {}
    for k in mydict:
        new_dict[mydict[k]] = k
    return new_dict


def load_compressed_ids_to_genes(compressed_filename):
    """
    Load mapping from compressed IDs to genes.
    """
    if not os.path.exists(compressed_filename):
        print "Error: %s compressed file does not exist." \
              %(compressed_filename)
        sys.exit(1)
    compressed_ids_to_genes = {}
    # Load mapping from gene IDs to their hashes
    compressed_ids_to_genes = shelve.open(compressed_filename)
    return compressed_ids_to_genes


def is_compressed_name(event_name):
    return str(event_name).startswith(COMPRESS_PREFIX)


def is_compressed_index(index_filename):
    """
    Check if the given index filename uses a compressed (hash)
    ID or not.
    """
    basename = os.path.basename(index_filename)
    if is_compressed_name(basename):
        return True
    return False


def make_dir(dirpath):
    if os.path.isfile(dirpath):
        print "Error: %s is a file!" %(dirpath)
        sys.exit(1)
    # Try to make the directory
    try:
        os.makedirs(dirpath)
    except OSError:
        pass


def pathify(f):
    return os.path.abspath(os.path.expanduser(f))


def which(program):
    """
    Check if program exists on path.
    """
    def is_exe(fpath):
        if not os.path.isfile(fpath):
            return False
        elif not os.access(fpath, os.X_OK):
            # If the file exists but is not executable, warn
            # the user
            print "WARNING: Found %s but it is not executable." %(fpath)
            print "Please ensure %s is executable." %(fpath)
            print "On Unix, use something like: "
            print "  chmod +x %s" %(fpath)
            time.sleep(10)
            return False
        return True
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def main():
    pass


if __name__ == "__main__":
    main()
