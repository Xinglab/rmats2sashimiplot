##
## Utilities for running scripts on cluster
##
import os
import time
import string
import subprocess

import misopy
import misopy.settings as settings
import misopy.misc_utils as misc_utils
from settings import Settings, load_settings

def write_cluster_preface(file_handle):
    module_preface = \
"""
module add numpy;		   
module add scipy;
module add matplotlib;\n
"""
    file_handle.write(module_preface)

def chunk_list(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0
  if len(seq) < num:
      return [seq]
  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg
  return out


def make_bash_script(filename, cmd, crate_dir=None):
    """
    Make an executable bash script out of the given command.
    """
#    os.system('ls %s' %(filename))
    if crate_dir == None:
        crate_dir = \
            os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
    f = open(filename, 'w')
    f.write("#!/bin/bash\n")
    f.write("export PATH=$PATH:%s\n" %(crate_dir))
    f.write("source ~/.bash_profile\n")
    f.write("cd %s\n" %(crate_dir))
    #write_cluster_preface(f)
    f.write(cmd + "\n")
    f.close()
    os.system('chmod +x \"%s\"' %(filename))

    
def valid_cluster_name(name):
    """
    Return valid qsub ID by removing semicolons, converting
    them into underscores.
    """
    name = name.replace(';', '_')
    return name


def run_SGEarray_cluster(arg_list, argfile, cluster_output_dir,
                         queue_type="long",
                         cluster_scripts_dir=None,
                         chunk=2500,
                         settings=None,
                         cmd_name="qsub",
                         job_name="miso_job"):
    """
    Run MISO jobs on cluster using SGE.

    Function contributed by Michael Lovci, UCSD.
    """
    misc_utils.make_dir(cluster_output_dir)
    # Create arguments file to pass on to job
    f = open(argfile, 'w')
    nargs = len(arg_list)
    if nargs % chunk == 0:
        njobs = nargs/chunk
    else:
        njobs = 1+(nargs/chunk)
    
    for args in arg_list:
        f.write(args[0] + "\n")
    f.close()
        
    if cluster_scripts_dir == None:
	cluster_scripts_dir = os.path.join(cluster_output_dir,
                                           'cluster_scripts')
    misc_utils.make_dir(cluster_scripts_dir)
    scripts_output_dir = os.path.join(cluster_output_dir,
                                      'scripts_output')
    misc_utils.make_dir(scripts_output_dir)
    scripts_output_dir = os.path.abspath(scripts_output_dir)
    script_error = os.path.join(scripts_output_dir,
                                string.join([job_name, "err"], "."))
    script_out = os.path.join(scripts_output_dir,
                              string.join([job_name, "out"], "."))
    cluster_script = os.path.join(cluster_scripts_dir, "run_miso.sh")

    if settings != None:
        load_settings(settings)
        cmd_name = Settings.get_cluster_command()

    if queue_type == "long":
        queue_name = Settings.get_long_queue_name()
    elif queue_type == "short":
        queue_name = Settings.get_short_queue_name()
    else:
        raise Exception, "Unknown queue type: %s" %(queue_type)

    if queue_type == None:
        print "  - queue: unspecified"
    else:
        print "  - queue: %s, using queue name %s" %(queue_type,
                                                     queue_name)
    cs = open(cluster_script, 'w')
    cs.write("#!/bin/sh" + "\n")
    cs.write("#$ -N %s\n" %(job_name))
    cs.write("#$ -S /bin/sh\n")
    cs.write("#$ -p -1023\n")
    cs.write("#$ -o %s\n" %(script_out))
    cs.write("#$ -e %s\n" %(script_error))
    cs.write("#$ -t 1-%s\n" %(njobs))
    
    ##execute from current working directory    
    cs.write("#$ -cwd\n")
    
    ## import environment variables
    cs.write("#$ -V\n") 
    if queue_name:
        cs.write("#$ -l %s\n" %(queue_name))
    cs.write("echo \"hostname is:\"\n")
    cs.write("hostname\n")
    cs.write("ARGFILE=%s\n" %argfile)
    cs.write("SEQ=/usr/bin/seq\n")
    cs.write("index=0\n")
    cs.write("lastindex=0\n")
    cs.write("let \"index = $SGE_TASK_ID * %s\"\n" %(chunk))
    chunk2 = chunk-1
    cs.write("let \"lastindex = $index - %s\"\n" %(chunk2))
    if chunk2 > 0:
        cs.write("for i in `$SEQ $lastindex $index`\n")
    else:
        cs.write("for i in $index\n") # if user chooses 1 for chunk size
    cs.write("do\n")
    cs.write("  line=$(cat $ARGFILE | head -n $i | tail -n 1)\n")
    cs.write("  eval $line\n")
    cs.write("done\n")
    cs.close()

    # Make script executable
    os.system('chmod +x \"%s\"' %(cluster_script))
    qsub_cmd = cmd_name + ' \"%s\"' %(cluster_script)

    os.system(qsub_cmd)


def run_on_cluster(cmd, job_name, cluster_output_dir,
                   cluster_scripts_dir=None,
                   queue_type=None,
                   cmd_name="qsub",
                   settings_fname=None):
    print "Submitting job: %s" %(job_name)
    queue_name = None

    # Load command name from settings file
    if settings_fname != None:
        load_settings(settings_fname)
        cmd_name = Settings.get_cluster_command()

    if queue_type == "long":
        queue_name = Settings.get_long_queue_name()
    elif queue_type == "short":
        queue_name = Settings.get_short_queue_name()
    else:
        print "Warning: Unknown queue type: %s" %(queue_type)
        queue_name = queue_type
    
    if queue_type is None:
        print "  - queue type: unspecified"
    else:
        print "  - queue type: %s" %(queue_type)
    if queue_name is None:
        print " - queue name unspecified"
    else:
        print " - queue name: %s" %(queue_name)
        
    misc_utils.make_dir(cluster_output_dir)
    if cluster_scripts_dir == None:
        cluster_scripts_dir = os.path.join(cluster_output_dir,
                                           'cluster_scripts')
        misc_utils.make_dir(cluster_scripts_dir)
    scripts_output_dir = os.path.join(cluster_output_dir,
                                      'scripts_output')
    misc_utils.make_dir(scripts_output_dir)
    scripts_output_dir = os.path.abspath(scripts_output_dir)
    cluster_call = '%s -o \"%s\" -e \"%s\"' %(cmd_name,
                                              scripts_output_dir,
                                              scripts_output_dir)
    # Add queue type if given one
    if queue_name != None:
        cluster_call += ' -q \"%s\"' %(queue_name)
        
    script_name = \
        valid_cluster_name(os.path.join(cluster_scripts_dir,
                                     '%s_time_%s.sh' \
                                     %(job_name,
                                       time.strftime("%m-%d-%y_%H:%M:%S"))))
    make_bash_script(script_name, cmd)
    cluster_cmd = cluster_call + ' \"%s\"' %(script_name)
    job_id = launch_job(cluster_cmd, cmd_name)
    return job_id


def launch_job(cluster_cmd, cmd_name):
    """
    Execute cluster_cmd and return its job ID if
    it can be fetched.
    """
    print "Executing: %s" %(cluster_cmd)
    proc = subprocess.Popen(cluster_cmd, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE)
    # Read the job ID if it's a known cluster
    # submission system
    output = proc.communicate()
    job_id = None
    if cmd_name == "qsub":
        if "." in output[0][:-1] and ">" not in output[0]:
            job_id = int(output[0].split(".")[0])
    elif cmd_name == "bsub":
        if "is submitted to" in output[0]:
            job_id = int(output[0].strip().split()[1][1:-1])                
    return job_id


def wait_on_job(job_id, cluster_cmd,
                delay=60):
    """
    Wait on a job given a job id.
    """
    if cluster_cmd == "qsub":
        # Handle qsub
        while True:
            output = \
                subprocess.Popen("qstat %i" %(job_id),
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()
            if "Unknown Job" in output[1]:
                break
            time.sleep(delay)
        time.sleep(delay)
    elif cluster_cmd == "bsub":
        # Handle bsub
        while True:
            output = subprocess.Popen("bjobs %i" %(job_id),
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
            if len(output[0]) > 0:
                status = output[0].split()[10]
                if status == "DONE":
                    break
            else:
                # No jobs available
                break
            time.sleep(delay)
        time.sleep(delay)
        
        
def wait_on_jobs(job_ids, cluster_cmd,
                 delay=120,
                 supported_cmds=["qsub",
                                 "bsub"]):
    """
    Wait on a set of job IDs.
    """
    if len(job_ids) == 0:
        return
    if cluster_cmd not in supported_cmds:
        return 
    num_jobs = len(job_ids)
    print "Waiting on a set of %d jobs..." %(num_jobs)
    curr_time = time.strftime("%x, %X")
    t_start = time.time()
    print "  - Starting to wait at %s" %(curr_time)
    completed_jobs = {}
    for job_id in job_ids:
        if job_id in completed_jobs:
            continue
        wait_on_job(job_id, cluster_cmd)
        print "  - Job ", job_id, " completed."
        completed_jobs[job_id] = True
    curr_time = time.strftime("%x, %X")
    t_end = time.time()
    print "Jobs completed at %s" %(curr_time)
    duration = ((t_end - t_start) / 60.) / 60.
    print "  - Took %.2f hours." %(duration)

    
