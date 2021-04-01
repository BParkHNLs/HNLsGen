'''
Script to retrieve for each mass-V2 point:
- generator filter efficiency => and average over the files
- number of events generated  => and sum them over the files
'''
import sys
import os
import subprocess
import glob

from python.common import Point

#import ROOT


def getOptions():

  # convention: no capital letters

  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('--points', type=str, dest='pointFile', help='name of file contaning information on scan to be run', default='points.py')

  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  ps = __import__(opt.pointFile.split('.py')[0])
  points = ps.points

  #filterEffs = {}
  nGenEvents = {}
  nTotEvents = {}
  timeEvents = {}
  tot_filterEffs = {}
  tot_nGenEvents = {}
  tot_nTotEvents = {}
  tot_timeEvents = {}
  tot_timeEvent = {}

  final_lines  = []
  #table = []

  for p in points:

    print('\n===> Processing new point')
    p.stamp_simpli()

    #filterEffs[p.name]=[]
    nGenEvents[p.name]=[]
    nTotEvents[p.name]=[]
    timeEvents[p.name]=[]

    logs = './{v}/logs/prod_mass{m}_ctau{ctau}_*.log'.format(v=opt.ver,m=p.mass,ctau=p.ctau)
    if len(glob.glob(logs)) == 0: continue
    for log in glob.glob(logs):
      #print(log)
      
      filter_eff = None
      n_acc = None
      time = None
    
      with open(log, 'r') as f:
        for line in f:
          if ('Filter efficiency' in line) and ('TO BE USED IN MC' in line):
            #print('****************** FILTER')
            filter_eff = float(line.split(' = ')[1].split(' +- ')[0])
            #print('********', filter_eff)
            n_acc = float(line.split('= (')[1].split(') / (')[0])
            n_tot = float(line.split('= (')[1].split(') / (')[1].split(') = ')[0])

          #if ('Wallclock running time step1: ' in line)  or ('   step2: ' in line) or ('   step3: ' in line) or  ('   step4: ' in line) or ('   tot: ' in line):
          #  times[].append(float(line.split(': ')[1].split(' s')[0]))

          if ('Wallclock running time: ' in line):
            time = float(line.split(': ')[1].split(' s')[0])

      #if filter_eff:
      #  filterEffs[p.name].append(filter_eff)
      if n_acc is not None and n_tot is not None and time is not None:
        #print(n_acc, n_tot)
        nGenEvents[p.name].append(n_acc)
        nTotEvents[p.name].append(n_tot)
        timeEvents[p.name].append(time)


      tot_filterEffs[p.name] = sum(nGenEvents[p.name])/sum(nTotEvents[p.name]) if sum(nTotEvents[p.name]) != 0 else 0
      tot_nGenEvents[p.name] = sum(nGenEvents[p.name]) 
      tot_nTotEvents[p.name] = sum(nTotEvents[p.name]) 
      tot_timeEvents[p.name] = sum(timeEvents[p.name])/len(timeEvents[p.name])
      tot_timeEvent[p.name]  = sum(timeEvents[p.name])/len(timeEvents[p.name])/sum(nGenEvents[p.name])


    # summary print out
    this_line = '{:12.1f} {:12.1e} {:12.2e} {:12.1f} {:12.1f} {:12.1f}  {:12.1f}'.format(p.mass,p.vv,tot_filterEffs[p.name], tot_nGenEvents[p.name], tot_nTotEvents[p.name], 
                                                                               tot_timeEvents[p.name],tot_timeEvent[p.name])
    print this_line
    final_lines.append(this_line)

  print('\nSummary table')
  print('\n{:12s} {:12s} {:12s} {:12s} {:12s} {:12s} {:12s} ').format('Mass', 'VV', 'Avg Filter Eff', 'NGen', 'NTot', 'Avg Time (s)', 'Avg Time / evt (s)')
  print('\n'.join(final_lines))



