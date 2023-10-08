import os
import sys
import subprocess
import argparse
from pathlib import Path
projectDir = Path(__file__).parent.absolute()
sys.path.append(os.path.join(projectDir, 'batch_jobs'))
print(projectDir)
from processDict import production

parser = argparse.ArgumentParser(description='Run QQbarAnalysis')
parser.add_argument('--process', type=str, required=True,
                    help='Production process (P2f_z_h, P4f_ww_h, P4_zz_h, Pqqh)')
parser.add_argument('--chiral',  type=str, required=True,
                    help='Polarization of beam (eLpR or eRpL or eLpL or eRpR)')

args  = parser.parse_args()
if args.chiral not in ["eLpR", "eRpL", "eLpL", "eRpR"]:
  sys.exit("Error: chiral must be eLpR or eRpL or eLpL or eRpR")
if args.process not in ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pqqh"]:
  sys.exit("Error: process must be P2f_z_h, P4f_ww_h, P4f_zz_h, or Pqqh")

chiralDot = "eL.pR" if args.chiral == "eLpR" else "eR.pL"

model = "l5"
processID  = production[args.process, args.chiral][0]
prodIDList = production[args.process, args.chiral][1]

inDir    = os.path.join(projectDir, 'rootfiles/tmp_root')
outDir   = os.path.join(projectDir, 'rootfiles/merged')
stageDir = os.path.join(outDir, 'staged')
dirCheck = os.path.isdir(inDir) or os.path.isdir(outDir) or os.path.isdir(stageDir)
if not dirCheck:
  sys.exit("Error: directory not found")

print("Removing staged directory...")
subprocess.run(['rm', '-rf', stageDir])
print("Creating staged directory...")
subprocess.run(['mkdir', '-p', stageDir])

for prodID in prodIDList:
  print(f"Processing {args.process} {args.chiral} {processID} {prodID}")
  subprocess.run(f"hadd -f -j 8 {stageDir}/stage_{args.process}_{args.chiral}_{processID}_{prodID}.root \
                  {inDir}/*{processID}*{args.process}*{chiralDot}*{prodID}*.root",shell=True)

# # subprocess.run(f"rm -rf {indir}",shell=True)
# # subprocess.run(f"mkdir -p {indir}",shell=True)

haddCommand = f"hadd -f -j 8 {outDir}/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I{processID}.{args.process}.{chiralDot}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.all.root \
                {stageDir}/stage_*.root"
print(haddCommand)
subprocess.run(haddCommand,shell=True)
