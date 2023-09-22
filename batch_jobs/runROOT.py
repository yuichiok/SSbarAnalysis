import os
import subprocess
import argparse
from pathlib import Path

batchDir = Path(__file__).parent.absolute()
exe = os.path.join(batchDir, 'main.exe')

parser = argparse.ArgumentParser(description='Start and end of iteration')
parser.add_argument('--flist',  type=str, required=True,
                    help='file list with full path separated by comma')

args  = parser.parse_args()
flist = args.flist.split(",")

logDir = "./log/"
for ifile in flist:
  filename = os.path.splitext(os.path.basename(ifile))[0]
  # log = os.path.join(logDir, f"{filename}.log")
  # subprocess.run(f"./main.exe {ifile} > {log} 2>&1",shell=True)
  print(exe, ifile)
  subprocess.run([exe, ifile])