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

  parser.add_argument('--pl', type=str, dest='ver', help='production label, e.g. V00_v00', default='V00_v00')
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
  sizeEvents = {}
  tot_filterEffs = {}
  tot_nGenEvents = {}
  tot_nTotEvents = {}
  tot_timeGenEvents = {}
  tot_timeGenEvent = {}
  tot_timeTotEvent = {}
  tot_sizeGenEvents = {}
  tot_sizeGenEvent = {}

  final_lines  = []
  final_lines_for_table = []
  #table = []

  for p in points:

    print('\n===> Processing new point')
    p.stamp_simpli()

    #filterEffs[p.name]=[]
    nGenEvents[p.name]=[]
    nTotEvents[p.name]=[]
    timeEvents[p.name]=[]
    sizeEvents[p.name]=[]

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

          if ('Wallclock running time' in line): #Wallclock running time: 10500 s
            #print line
            #print line.split(': ')[1].split(' s')
            time = float(line.split(': ')[1].split(' s')[0])

      #if filter_eff:
      #  filterEffs[p.name].append(filter_eff)
      #print(n_acc, n_tot, time)
      if n_acc is not None and n_tot is not None and time is not None:
        #print(n_acc, n_tot)
        nGenEvents[p.name].append(n_acc)
        nTotEvents[p.name].append(n_tot)
        timeEvents[p.name].append(time)


    tot_filterEffs[p.name]    = sum(nGenEvents[p.name])/sum(nTotEvents[p.name]) if sum(nTotEvents[p.name]) != 0 else 0
    tot_nGenEvents[p.name]    = sum(nGenEvents[p.name]) 
    tot_nTotEvents[p.name]    = sum(nTotEvents[p.name]) 
    tot_timeGenEvents[p.name] = sum(timeEvents[p.name])/len(timeEvents[p.name]) if len(timeEvents[p.name]) !=0 else 0
    tot_timeGenEvent[p.name]  = sum(timeEvents[p.name])/len(timeEvents[p.name])/(sum(nGenEvents[p.name])/len(nGenEvents[p.name])) if tot_timeGenEvents[p.name]!=0 else 0
    tot_timeTotEvent[p.name]  = sum(timeEvents[p.name])/len(timeEvents[p.name])/(sum(nTotEvents[p.name])/len(nTotEvents[p.name])) if tot_timeGenEvents[p.name]!=0 else 0


    # get the size from the root files
    rootfiles = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/{v}/mass{m}_ctau{ctau}/step1*.root'.format(v=opt.ver,m=p.mass,ctau=p.ctau)
    if len(glob.glob(rootfiles)) == 0: print 'WARNING: no root files!'
    for rootfile in glob.glob(rootfiles):
      sizeKB = os.path.getsize(rootfile)/(1024.)
      sizeEvents[p.name].append(sizeKB)
    tot_sizeGenEvents[p.name] = sum(sizeEvents[p.name])
    tot_sizeGenEvent[p.name] = sum(sizeEvents[p.name])/tot_nGenEvents[p.name]

    # summary print out
    this_line = '{:12.1f} {:12.1e} {:12.2e} {:12.1f} {:12.1f} {:12.1f} {:12.0f} {:12.3f} {:12.1f}'.format(p.mass,p.ctau,tot_filterEffs[p.name], tot_nGenEvents[p.name], 
                tot_nTotEvents[p.name], tot_timeGenEvents[p.name],tot_timeGenEvent[p.name], tot_timeTotEvent[p.name], tot_sizeGenEvent[p.name])

    this_line_for_table = '({m:.1f}, {ct:8.1f}, {eff:.2e}, {time:.2f}, {size:.1f}),'.format(m=p.mass,ct=p.ctau,eff=tot_filterEffs[p.name],
                                                                                            time=tot_timeGenEvent[p.name],size=tot_sizeGenEvent[p.name]) 
    #(0.5,100000.,1.00e-04,800),

    print this_line
    final_lines.append(this_line)
    final_lines_for_table.append(this_line_for_table)

  print('\nSummary table')
  print('\n{:12s} {:12s} {:12s} {:12s} {:12s} {:12s} {:12s} {:12s}').format('Mass', 'ctau(mm)', 'Avg Filter Eff', 'NGen', 'NTot', 'Avg Time (s)', 'AvgTime/evt (s)', 'AvgSize/evt (MB)')
  print('\n'.join(final_lines))

  print('\nTable for point file')
  print('\n'.join(final_lines_for_table))

