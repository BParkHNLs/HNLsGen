'''
Very unsmart helper to handle resubmission from step2 on 
'''


import sys
import os
import subprocess

from python.common import Point
from python.decays import Decays

def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Resubmission helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)
  parser.add_argument('--iter', type=int, dest='iteration', help='iteration of resubmission', default=1)

  return parser.parse_args()

if __name__ == "__main__":


  opt = getOptions()
  print('===> Resubmission for production {}'.format(opt.ver))
  os.chdir(opt.ver)
  
  if opt.iteration==1:
    f = open('./jobs.txt', 'r')
  else:
    f = open('./resubjobs_{j}.txt'.format(j=opt.iteration-1), 'r')  
  
  fo= open('./resubjobs_{i}.txt'.format(i=opt.iteration), 'w')
  for line in f: 
    #### Operations per Job
    name,job = line.split()
    command = 'sacct -j {job} --format="JobID,JobIDRaw,State,ExitCode,CPUTime,MaxRSS,AveRSS,Priority,Reason,Timelimit,Submit,TimelimitRaw" > status.txt'.format(job=job)
    out = subprocess.check_output(command, shell=True) 
    failed_arrs = []

    fstatus = open('./status.txt', 'r')
    for l in fstatus:
      #### Operations per Task in each Job Array
      if 'JobID' in l or '------------' in l or 'batch' in l: continue
      #if 'batch' not in l: continue
      data = filter(lambda x: x is not '', l.split(' '))
      task_id_full = data[0]
      task_state = data[2]
      if task_state == 'FAILED':
        failed_arrs.append(task_id_full.split('_')[1])
    fstatus.close()

    print('  JobName : {}, Job : {}'.format(name,job)) 
    print('  failed elements of array: {}'.format(len(failed_arrs)))

    fslurm = open('slurm_{name}_prod.sh'.format(name=name), 'r')
    a = fslurm.readlines()
    fslurm.close()
    arrl = '#SBATCH --array={}\n'.format(','.join(failed_arrs))
    copl = [
      'cd $WORKDIR\n', 
      'pwd\n', 
      'echo "Going to copy step1 output to work directory"\n', 
      'xrdcp $OUTSEPREFIX/$SERESULTDIR/step1_nj$SLURM_ARRAY_TASK_ID".root" $WORKDIR/BPH-step1.root\n',
      'if [ $? -eq 0 ]; then echo "Successfully copied step1.root file"; else exit $?; fi\n',
      '\n',]
    args  = a[67].split()
    step2_command = args[0] + ' ' + args[1] + ' ' + args[2] + ' ' + args[3] + ' ' + args[4] + ' ' + 'inputFile=BPH-step1.root' + ' ' + args[6] + ' ' + args[7] + '\n'  
    new_lines = a[0:8] + [arrl] + a[9:34] + copl + a[60:67] + [step2_command] + a[68:]
        
    with open('resubmit_{i}_{name}_prod.sh'.format(i=opt.iter,name=name), 'w') as fresub:
      fresub.writelines(new_lines)

    print('  Written resubmitter')
    
    command = 'sbatch resubmit_{i}_{name}_prod.sh'.format(i=opt.iteration,name=name)
    if opt.dosubmit:
      out = subprocess.check_output(command, shell=True)
      jobn = out.split('\n')[0].split(' ')[-1]      
      fo.write('{name} {j}\n'.format(name,j=jobn))
      print('  Resubmission sent')
  
  fo.close()
  f.close()

