"""
@File       :  nuc.py
@Copyright  :  Yuanting
@Author     :  Yuantiing Ma
@Date       :  2023/7/22
"""
from argparse import ArgumentParser
import sys
import os

PROG_NAME = 'nuc_process'
VERSION = '1.3.2'
DESCRIPTION = 'Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools'
RE_CONF_FILE = 'enzymes.conf'
epilog = 'Note %s can be edited to add further restriction enzyme cut-site defintions. ' % RE_CONF_FILE
arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                           epilog=epilog, prefix_chars='-', add_help=True)

arg_parse.add_argument(nargs='+', metavar='FASTQ_FILE', dest='i',
                       help='Input paired-read FASTQ files to process. Accepts wildcards' \
                            ' that match paired files. If more than two files are input,' \
                            ' processing will be run in batch mode using the same parameters.')
arg_parse.add_argument('-n', '--num-cpu', default=1, metavar='CPU_COUNT', dest='n',
                         type=int, help='Number of CPU cores to use in parallel')

def pair_fastq_files(fastq_paths, pair_tags=('r_1','r_2'), err_msg='Could not pair FASTQ read files.'):

  if len(fastq_paths) != len(set(fastq_paths)):
    msg = '%s Repeat file path present.'
    print(msg % (err_msg))

  t1, t2 = pair_tags

  paths_1 = []
  paths_2 = []

  for path in fastq_paths:
    dir_name, base_name = os.path.split(path)

    if (t1 in base_name) and (t2 in base_name):
      msg = '%s Tags "%s" and "%s" are ambiguous in file %s'
      print(msg % (err_msg, t1, t2, base_name))

    elif t1 in base_name:
      paths_1.append((path, dir_name, base_name))

    elif t2 in base_name:
      paths_2.append((path, dir_name, base_name))

    else:
      msg = '%s File name %s does not contain tag "%s" or "%s"'
      print(msg % (err_msg, base_name, t1, t2))
pair_fastq_files(['SRR18579584r_1.fq', 'SRR18579584_r_2.fq'])

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
    args = vars(arg_parse.parse_args(argv))
    fastq_paths = args['i']

    print(fastq_paths)
    fastq_pairs = fastq_paths
    pair_fastq_files(fastq_pairs)

main()