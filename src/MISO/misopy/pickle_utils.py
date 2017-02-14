##
## Wrappers for pickling
##
import misopy
import os
import cPickle as pickle

def load_pickled_file(pickled_filename):
    if os.access(pickled_filename, os.F_OK):
        pickled_file = open(pickled_filename, 'r')
        loaded_obj = pickle.load(pickled_file)
        pickled_file.close()
        return loaded_obj
    return None

def write_pickled_file(obj_to_pickle, pickled_filename):
    pickled_file = file(pickled_filename, 'w')
    pickle.dump(obj_to_pickle, pickled_file, -1)
    #del pickled_file
    pickled_file.close()
