import os
import sys
import time
import subprocess
import argparse
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path
projectDir = Path(__file__).parent.absolute()
sys.path.append(os.path.join(projectDir, 'batch_jobs'))
from processDict import production

parser = argparse.ArgumentParser(description='Run QQbarAnalysis')
parser.add_argument('--process', type=str, required=True,
                    help='Production process (P2f_z_h, P4f_ww_h, P4_zz_h, Pe1e1h)')
parser.add_argument('--chiral', type=str, required=True,
                    help='Polarization of beam (eLpR or eRpL or eLpL or eRpR)')
parser.add_argument('--mergeMode', type=int, required=True,
                    help='Merge mode (0: tmp -> stage, 1: stage->out)')
parser.add_argument('--rmStage', action='store_true',
                    help='Remove staged directory and create a new one')

args  = parser.parse_args()
if args.chiral not in ["eLpR", "eRpL", "eLpL", "eRpR"]:
  sys.exit("Error: chiral must be eLpR or eRpL or eLpL or eRpR")
if args.process not in ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"]:
  sys.exit("Error: process must be P2f_z_h, P4f_ww_h, P4f_zz_h, or Pe1e1h")

chiralDot = "eL.pR" if args.chiral == "eLpR" else "eR.pL"

model = "l5"
processID  = production[args.process, args.chiral][0]
prodIDList = production[args.process, args.chiral][1]


tmpDir   = os.path.join(projectDir, 'rootfiles/tmp_root')
outDir   = os.path.join(projectDir, 'rootfiles/merged')
outROOT  = os.path.join(outDir, f'rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{args.process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root')
stageDir = os.path.join(outDir, 'staged')
dirCheck = os.path.isdir(tmpDir) or os.path.isdir(outDir) or os.path.isdir(stageDir)
if not dirCheck:
  sys.exit("Error: directory not found")

def stageROOTFileName(prodID):
  return os.path.join(stageDir, f'stage_{args.process}_{args.chiral}_{processID}_{prodID}.root')

def mergeStage(prodID):
    stageROOT = stageROOTFileName(prodID)
    tmpROOT   = os.path.join(tmpDir, f'*{processID}*{args.process}*{chiralDot}*{prodID}*.root')

    condor = f"bsub -q s -J merge_{processID}_{prodID}"
    stageCommand = f"hadd -f -v 0 -j 24 {stageROOT} {tmpROOT}"
    # subprocess.run(f"{condor} {stageCommand}",shell=True)
    subprocess.run(f"{stageCommand}",shell=True)

if __name__ == '__main__':
  
  if args.mergeMode == 0:
    if args.rmStage:
      print("Removing staged directory...")
      subprocess.run(['rm', '-rf', stageDir])
      print("Creating staged directory...")
      subprocess.run(['mkdir', '-p', stageDir])

    for prodID in prodIDList:
      start_time = time.time()
      print(f"Staging {prodID}...")
      mergeStage(prodID)
      print(f"Staging {prodID} took {time.time() - start_time} seconds")
  
  else:
    stageROOTFileList = [stageROOTFileName(prodID) for prodID in prodIDList]
    stageROOTFiles = " ".join(stageROOTFileList)
    haddCommand = f"hadd -f -j 8 {outROOT} {stageROOTFiles}"
    print(haddCommand)
    subprocess.run(haddCommand,shell=True)