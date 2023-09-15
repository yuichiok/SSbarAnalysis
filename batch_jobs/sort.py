import os
import sys
import glob
import subprocess
import argparse
from processDict import production
from pathlib import Path

projectDir = Path(__file__).parent.absolute()

parser = argparse.ArgumentParser(description='Run QQbarAnalysis')
parser.add_argument('--process', type=str, required=True,
                    help='Production process (P2f_z_h, P4f_ww_h, P4_zz_h, Pe1e1h)')
parser.add_argument('--chiral',  type=str, required=True,
                    help='Polarization of beam (eLpR or eRpL)')

args  = parser.parse_args()

model = "l5"
processID  = production[args.process, args.chiral][0]
prodIDList = production[args.process, args.chiral][1]

print(args.process, args.chiral, processID)
print(prodIDList)

for prodID in prodIDList:
  rootOut  = f"/group/ilc/users/yokugawa/QQbar250/{model}/{args.process}/{args.chiral}/{processID}/{prodID}/dEdx_corr/QQbarProcessor_out/"
  
  if not os.path.isdir(rootOut):
    sys.exit(f'Directory {rootOut} does not exist.')

  file_list = []
  try:
    pattern = os.path.join(rootOut, f'*I{processID}*{prodID}*.root')
    file_list = sorted(glob.glob(pattern))
  except subprocess.CalledProcessError as e:
    print(f"Error: {e}")
    file_list = []

  nf = len(file_list)
  print("number of files for process:", nf)
