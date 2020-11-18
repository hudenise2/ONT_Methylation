#! /usr/bin/python3

__author__ = "Hubert DENISE (2020)"

import sys
sys.path.append("/home/had38/.local/lib/python3.6/site-packages")

import numpy
numpy.set_printoptions(threshold=sys.maxsize)

from ont_fast5_api.fast5_interface import get_fast5_file
if len(sys.argv) <2:
    print ("Please provide the path to a fast5 file")
    print ("Usage: extract_mod.py <path to fast5 file> > outfile")
    sys.exit()
else:
    fast5_filepath=sys.argv[1]

    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read_id in f5.get_read_ids():
            read = f5.get_read(read_id)
            latest_basecall = read.get_latest_analysis('Basecall_1D')
            mod_base_table = read.get_analysis_dataset(
            latest_basecall, 'BaseCalled_template/ModBaseProbs')
            table_path = '{}/BaseCalled_template/ModBaseProbs'.format(latest_basecall)
            metadata = read.get_analysis_attributes(table_path)
            print(read_id, "\n", mod_base_table, metadata['modified_base_long_names'])

