import os
import glob
import ROOT
import math


def getFirstIndex(line, idx):
  return line.rfind(' ', 0, idx-1)

def getLastIndex(idx):
  return idx

def getItem(line, idx):
  idx_ini = getFirstIndex(line, idx)
  idx_fin = getLastIndex(idx)
  qte = line[idx_ini: idx_fin+1]
  return idx_ini, qte

def getQuantity(num):
  index = line.rfind(' ')
  for i in range(0, num):
    idx, qte = getItem(line, index)
    index = idx-1
  return qte

def convertToFloat(qte):
  power = qte[qte.find('E')+1:]
  return float(qte) * pow(1, float(power))

def getPx(line):
  return convertToFloat(getQuantity(6))

def getPy(line):
  return convertToFloat(getQuantity(5))

def getPz(line):
  return convertToFloat(getQuantity(4))

def getEnergy(line):
  return convertToFloat(getQuantity(3))

print 'will analyse LHE file'

#path = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/LHE_files/Bc_LHE_2M'
#path = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/LHE_files/Bc_LHE_600M'
path = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/LHE_files/Bc_LHE_2M_NoCut2022'

max_events = 2000000

# get the files
lhe_files = [f for f in glob.glob('{}/*.lhe'.format(path))]

hist_pt = ROOT.TH1D('hist_pt', 'hist_pt', 80, 0, 100)
hist_y = ROOT.TH1D('hist_y', 'hist_y', 80, -5, 5)

hist_num = ROOT.TH1D('hist_num', 'hist_num', 2, 0, 2)
hist_deno = ROOT.TH1D('hist_deno', 'hist_deno', 2, 0, 2)

count = 0
do_skip = False
for lhe_file in lhe_files:
  the_lhe_file = open(lhe_file)
  lines = the_lhe_file.readlines()
  for line in lines:
    if ' 541 ' not in line: continue

    count += 1
    if max_events != -1 and count > max_events: 
      do_skip = True
      break

    px = getPx(line)
    py = getPy(line)
    pz = getPz(line)
    energy = getEnergy(line)

    p4 = ROOT.TLorentzVector()
    p4.SetPxPyPzE(px, py, pz, energy)

    hist_pt.Fill(p4.Pt())
    hist_y.Fill(p4.Rapidity())

    # for efficiency
    hist_deno.Fill(1.)

    if p4.Pt() > 8. and abs(p4.Rapidity()) < 2.5:
      hist_num.Fill(1.)

  if do_skip: break

canv_pt = ROOT.TCanvas('canv_pt', 'canv_pt', 900, 800)
hist_pt.Draw('hist')
canv_pt.SaveAs('./plots/lhe_Bc_pt.png')

canv_y = ROOT.TCanvas('canv_y', 'canv_y', 900, 800)
hist_y.Draw('hist')
canv_y.SaveAs('./plots/lhe_Bc_y.png')

efficiency = hist_num.Integral() / hist_deno.Integral()
err = efficiency * (math.sqrt(hist_num.Integral()) / hist_num.Integral() + math.sqrt(hist_deno.Integral()) / hist_deno.Integral())
print 'efficiency = {} / {} = ({} +- {})%'.format(hist_num.Integral(), hist_deno.Integral(), efficiency*100, err*100)

print 'Done'






