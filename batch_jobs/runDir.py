import os
import sys
import glob
import subprocess
import argparse
from processDict import production
from pathlib import Path

batchDir = Path(__file__).parent.absolute()
projectDir = os.path.dirname(batchDir)

runROOT = os.path.join(batchDir, 'runROOT.py')

nfirst    = 1    # first file
nlast_set = -1   # -1: all files
nrun      = 50  # number of runs per job

isAll = False
if nlast_set <= 0:
  isAll = True

parser = argparse.ArgumentParser(description='Run QQbarAnalysis')
parser.add_argument('--process', type=str, required=True,
                    help='Production process (P2f_z_h, P4f_ww_h, P4_zz_h, Pe1e1h)')
parser.add_argument('--chiral',  type=str, required=True,
                    help='Polarization of beam (eLpR or eRpL or eLpL or eRpR)')
parser.add_argument('--copy',type=int, required=True,
                    help='Copy executable and setting files. Set to false when running consecutive jobs.')

args  = parser.parse_args()
if args.chiral not in ["eLpR", "eRpL", "eLpL", "eRpR"]:
  sys.exit("Error: chiral must be eLpR or eRpL or eLpL or eRpR")
if args.process not in ["P2f_z_h", "P4f_ww_h", "P4f_zz_h", "Pe1e1h"]:
  sys.exit("Error: process must be P2f_z_h, P4f_ww_h, P4f_zz_h, or Pe1e1h")
if args.copy:
  print("Copy exe and setting files...")
  subprocess.run(['cp', os.path.join(projectDir, 'main.exe'), batchDir])
  subprocess.run(f'cp {os.path.join(projectDir, "etc")}/* {batchDir}/etc',shell=True)


model = "l5"
processID  = production[args.process, args.chiral][0]
prodIDList = production[args.process, args.chiral][1]

print(args.process, args.chiral, processID)
print(prodIDList)

# for processID in processIDs:
for prodID in prodIDList:
  data_input = f"/group/ilc/users/yokugawa/QQbar250/{model}/{args.process}/{args.chiral}/{processID}/{prodID}/dEdx_corr/QQbarProcessor_out/"
  if not os.path.isdir(data_input):
    sys.exit(f'Directory {data_input} does not exist.')

  try:
    pattern = os.path.join(data_input, f'*I{processID}*{prodID}*.root')
    file_list = sorted(glob.glob(pattern))
  except subprocess.CalledProcessError as e:
    print(f"Error: {e}")
    file_list = []

  nf = len(file_list)
  print("number of files for process:", nf)

  nlast = nlast_set
  if isAll or nlast_set >= nf:
    nlast = nf
    
  njobs = (nlast - nfirst) // nrun + 1
  print("processID {}, prodID {}: nlast {}, nfirst {}, nrun {}, njobs {}".format(processID,prodID, nlast, nfirst, nrun, njobs))
  seqlist = range(1, njobs + 1)
  for seq in seqlist:
    nrun0 = nfirst + nrun * (seq - 1)
    nrun1tmp = nrun0 + nrun - 1

    if nrun1tmp <= nlast:
      nrun1 = nrun1tmp
    else:
      nrun1 = nlast

    print(prodID, nrun0, nrun1, nrun, njobs)

    sub_flist = file_list[nrun0 - 1:nrun1]
    arg_flist = ",".join(sub_flist)

    subprocess.run(['bsub', '-q', 's', '-J', f'ana_{prodID}_{seq:03d}', f'python3', runROOT, '--flist', arg_flist])
    # subprocess.run([f'python3', runROOT, '--flist', arg_flist])
