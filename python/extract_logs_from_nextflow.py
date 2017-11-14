#!/usr/bin/env python

import os
import sys
import pandas as pd
from glob import glob
import subprocess

def main():
    if len(sys.argv) != 3:
        print("Usage: extract_logs_from_nextflow.py <nextflow_trace_file> <task_tag>:<stderr|stdout|both>",file=sys.stderr)
        
    # Load Nextflow trace
    nextflow_trace = pd.read_csv(sys.argv[1],sep="\t")
    #Get task you are interested in and whether you want stderr/stdout/both
    # Format <task>:<stderr|stdout|both>
    task_id,output_type = sys.argv[2].split(":")

    #Create dirs
    out_dir="nextflow_logs/{}".format(task_id)
    os.makedirs(out_dir)

    # Subset tasks of interest
    my_tasks = list(nextflow_trace[ (nextflow_trace.process == task_id) ][["hash","tag","status"]].itertuples(index=False,name=None))

    if len(my_tasks) == 0:
        print("No tasks were found",file=sys.stderr)

    # Iterate through tasks
    for t_hash,t_tag,t_status in my_tasks:
        task_dir= get_hash_directory(t_hash)
        if not task_dir:
            print("Error: work/{}* directory was not found".format(t_hash))
            continue
        print("{}: {}".format(t_tag,task_dir))

        out_prefix="{}_{}_{}".format(t_tag,t_status[0],t_hash.replace("/","_"))

        if output_type != "stderr":
            copy_file_into_dir("{}/.command.out".format(task_dir),out_dir,prefix=out_prefix)
        if output_type != "stdout":
            copy_file_into_dir("{}/.command.err".format(task_dir),out_dir,prefix=out_prefix)

# Helping functions
def get_hash_directory( h ):
    my_task_dir = None
    matching_dirs = glob("work/{}*".format(h))
    if len(matching_dirs) == 1:
        my_task_dir = matching_dirs[0]
    return my_task_dir

def copy_file_into_dir(my_file,my_dir,prefix=""):
    print("\t{}".format(my_file))
    subprocess.check_call(["cp",my_file,"{}/{}.{}".format(my_dir,prefix,my_file[-3:])])

if __name__ == '__main__':
    main()
