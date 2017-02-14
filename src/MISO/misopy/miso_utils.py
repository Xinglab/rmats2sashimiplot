##
## Utilities for parsing MISO output files
##
import os
import glob

import misopy
import misopy.misc_utils as utils


def get_miso_files_from_dir(dirname):
    """
    Return MISO output files from a directory.
    """
    miso_basename_files = []
    if not os.path.isdir(dirname):
        print "Error: %s not a directory." \
              %(dirname)
        return miso_basename_files
    miso_files = glob.glob(os.path.join(dirname, "*.miso"))
    # return basenames
    miso_basename_files = [os.path.basename(f) for f in miso_files]
    return miso_basename_files


def is_miso_rawdir(dirname, miso_ext=".miso"):
    """
    Return true if the given directory is one which contains
    raw MISO output (i.e. if it contains MISO chromosome directories.)
    """
    dirname = utils.pathify(dirname)
    if not os.path.isdir(dirname):
        return False
    filenames = glob.glob(os.path.join(dirname, "*.miso"))
    # If all of the subdirectories of the current directory
    # contain raw MISO (*.miso) files, consider it a raw MISO
    # output directory.
    if len(filenames) == 0:
        return False
    for fname in filenames:
        if not fname.endswith(miso_ext):
            return False
    return True


def get_miso_output_files(event_name, chrom, settings):
    """
    Get MISO output files, in order of 'miso_files'

    Look recursively in subdirectories of MISO prefix.
    """
    miso_filenames = []
    
    # Apply MISO prefix path if given
    if "miso_prefix" in settings:
        miso_prefix = \
            os.path.abspath(os.path.expanduser(settings["miso_prefix"]))
    else:
        miso_prefix = ""

    print "miso_prefix: %s" %(miso_prefix)

    if "miso_files" not in settings:
        print "Error: need \'miso_files\' to be set in settings file in " \
              "order to plot MISO estimates."
        return miso_filenames

    miso_files = settings['miso_files']

    miso_sample_paths = \
        [os.path.abspath(os.path.expanduser(os.path.join(miso_prefix, f))) \
         for f in miso_files]

    event_with_miso_ext = "%s.miso" %(event_name)

    for curr_sample_path in miso_sample_paths:
        event_found = False
        print "Searching for MISO files in: %s" %(curr_sample_path)
        print "  - Looking for chromosome %s directories" %(chrom)

        if event_with_miso_ext in get_miso_files_from_dir(curr_sample_path):
            # Allow the event to be in a top-level directory outside of a
            # chromosome folder
            event_found = True
            event_filename = os.path.join(curr_sample_path,
                                          event_with_miso_ext)
            miso_filenames.append(event_filename)
            print "Found %s MISO file in top-level directory." %(event_name)
            print "  - Location: %s" %(event_filename)
            print "Please try to keep MISO event files in their chromosome "\
                  "directory."
            break

        for root, dirs, files in os.walk(curr_sample_path):
            # First check if the file is in the current directory
            ### TODO FILL ME IN

            # If there's a directory named after the event's chromosome,
            # see if the MISO file is in there
            if chrom in dirs:
                chrom_dirname = os.path.abspath(os.path.join(root, chrom))
                print "Looking for MISO files in: %s" %(chrom_dirname)
                # Fetch MISO files, if any
                curr_miso_files = get_miso_files_from_dir(chrom_dirname)

                # Is the event in there?
                if event_with_miso_ext in curr_miso_files:
                    # Found relevant event
                    event_found = True
                    # Add to list
                    event_filename = os.path.join(root, chrom,
                                                  event_with_miso_ext)
                    print "Found %s MISO file." %(event_name)
                    print "  - Location: %s" %(event_filename)
                    miso_filenames.append(event_filename)
                    break

        if not event_found:
            # If we're here, it means we couldn't find the MISO
            # output files for the current sample
            print "Error: Could not find MISO output files for " \
                  "sample %s (after searching in %s and its subdirectories). " \
                  "Are you sure MISO output files are present in that " \
                  "directory?" %(os.path.basename(curr_sample_path),
                                 curr_sample_path)
            # Include empty path for this sample
            miso_filenames.append('')

    # Number of event filenames retrieved must equal
    # the number of samples to be plotted
    if len(miso_filenames) != len(miso_files):
        print "WARNING: Could not find MISO files for all samples."
        print "  - miso_filenames: ", miso_filenames
        print "  - miso_samples to be plotted: ", miso_files

    return miso_filenames
