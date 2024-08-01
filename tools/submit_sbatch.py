import os, sys, re, argparse,glob
import numpy as np


sbatch_memory        = 3000  # Memory in MB per default
sbatch_time          = '0-20:00:00'
sbatch_email         = 'dfranco@in2p3.fr'

parser = argparse.ArgumentParser(description='g4ds submitter')
parser.add_argument('-f', '--folder',         default=None,       help='Output folder. If not defined, dir name is equal to tag name')
parser.add_argument('-t', '--tag',            required=True,      help='Output tag name.')
parser.add_argument('-i', '--input',          required=True,      help='Input mac file name.')
parser.add_argument('-j', '--number-of-jobs', default=1,          help='Number of jobs')
parser.add_argument('-n', '--nevents',        default=0,          help='Number of events. If 0, default in mac file is not changed')
parser.add_argument('-s', '--seed',           default=12345678,   help='first seed number. Default: 12345678')
parser.add_argument('-e', '--exe',            default=False,      action="store_true", help='Execute')
parser.add_argument('-r', '--rooter',         default=False,      action="store_true", help='Run rooter at the end of the process')
parser.add_argument('-bd', '--biasdir',        default=False,  help='folder containing input files')
args = vars(parser.parse_args())  # Convert to dictionary.



def get_job_file(job_name, log_name, mac_name, pwd, memory=sbatch_memory, time=sbatch_time, email=sbatch_email):
  context = {
   'job_name': job_name,
   'log_name': log_name,
   'mac_name': mac_name,
   'pwd':      pwd,
   'memory':   memory,
   'time':     time,
   'email':    email
  }
  template = """#!/bin/sh
# SLURM options:
#SBATCH --job-name={job_name}     # Job name
#SBATCH --output={log_name}       # Standard output and error log
#SBATCH --partition htc               # Partition choice
#SBATCH --ntasks 1                    # Run a single task (by default tasks == CPU)
#SBATCH --mem {memory}         # Memory in MB per default
#SBATCH --time {time}          # 7 days by default on htc partition
#SBATCH --mail-user={email}      # Where to send mail
#SBATCH --mail-type=FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)

#module load /sps/nusol/software/miniconda3/bin/conda
source activate
#conda activate /sps/nusol/software/miniconda3/envs/ds
#conda env list
source /sps/nusol/software/set_environment.sh

cd {pwd}
./g4ds {mac_name}
"""
  return template.format(**context)


def change_event_numbers(lines, n):
  outs = lines
  for i,l in enumerate(outs):
    m = re.match('.*beamOn\s+(\d+)', l)
    if m:  outs[i] = f'/run/beamOn {n}\n'
  return outs

def change_seed(lines, s):
  outs = lines
  for i,l in enumerate(outs):
    m = re.match('.*heprandomseed\s+(\d+)', l)
    if m:  outs[i] = f'/run/heprandomseed {s}\n'
  return outs

def change_output(lines, s):
  outs = lines
  for i,l in enumerate(outs):
    m = re.match('.*filename.*', l)
    if m:  outs[i] = f'/run/filename {s}\n'
  return outs

def change_input(lines, s):
  outs = lines
  for i,l in enumerate(outs):
    m = re.match('.*/ds/generator/bias/file.*', l)
    if m:  outs[i] = f'/ds/generator/bias/file {s}\n'
  return outs



pwd = os.getenv('PWD')
dir = args['folder']
tag = args['tag']

njobs = int(args['number_of_jobs'])
if(args['biasdir']):
  input_dir = args['biasdir']
  print(input_dir)
  input_fils = glob.glob(f'{input_dir}/*.fil')
  njobs = int(len(input_fils))



nevents = int(float(args['nevents']))
seed = int(args['seed'])
exe  = bool(args['exe'])
rooter  = bool(args['rooter'])



if dir is None: dir = tag
if not os.path.exists(dir):
  os.makedirs(dir)
  print(f'The {dir} directory is created')

dir = f'{pwd}/{dir}'


with  open(args['input'], 'r') as f: lines = f.readlines()

if nevents >0: lines = change_event_numbers(lines, nevents)

for k in np.arange(njobs):

  if(args['biasdir']):
    fil = input_fils[k]
    j = int(fil[:-4].split('_')[-1])
  else:
    j = k
  
  job_tag_name = f'{dir}/job_{tag}_{j}.sh'
  out_tag_name = f'{dir}/out_{tag}_{j}'
  mac_tag_name = f'{dir}/mac_{tag}_{j}.mac'
  log_tag_name = f'{dir}/log_{tag}_{j}.log'

  lines = change_seed(lines, seed+j)
  lines = change_output(lines, out_tag_name)
  if(args['biasdir']): lines = change_input(lines, fil)

  with open(mac_tag_name, 'w') as f:
    for l in lines: f.write(l)

  job = get_job_file(job_name=job_tag_name, log_name=log_tag_name, mac_name=mac_tag_name, pwd=pwd)
  with open(job_tag_name, 'w') as f:
    for l in job: f.write(l)
    if rooter:
      f.write(f'{pwd}/../tools/g4rooter {out_tag_name}.fil\n')

  if exe:
    print(job_tag_name)
    os.system(f'sbatch -L sps {job_tag_name}')
