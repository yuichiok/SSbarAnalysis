import os
import glob
import subprocess

# chiral
# chiral = "eLpR"
chiral = "eRpL"

# processID
processIDs = []
if chiral == "eLpR":
  processIDs = [15162,15271,15273,15315,15319,15351,15353,15355,15357,15463,15465,15467,15469,15471,15473,15475,15477,15479,15499,15501,15503,15505,15507,15509,15511,15513,15515,15517,15539,15541,15549,15551,15553,15555,15557,15559,15561,15563,15585,15587,15589,15591,15593,15595,15597,15599,15601,15603,15629,15633]
elif chiral == "eRpL":
  processIDs = [15165,15275,15277,15279,15359,15361,15363,15365,15367,15369,15481,15483,15485,15487,15489,15491,15493,15495,15497,15519,15521,15523,15525,15527,15529,15531,15533,15535,15537,15565,15567,15569,15571,15573,15575,15577,15579,15581,15583,15605,15607,15609,15611,15613,15616,15618,15620,15625,15627,15631]

subprocess.run('cp ../main.exe ./',shell=True)
subprocess.run('cp ../etc/* ./etc',shell=True)

nfirst    = 1    # first file
nlast_set = -1   # -1: all files
nrun      = 100  # number of runs per job

isAll = False
if nlast_set <= 0:
  isAll = True

for processID in processIDs:
  data_input = f"/group/ilc/users/yokugawa/QQbar250/l5/{chiral}/{processID}/dEdx_corr/QQbarProcessor_out/"

  file_list = []
  try:
    pattern = os.path.join(data_input, f'*{processID}*.root')
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
  print("processID {}: nlast {}, nfirst {}, nrun {}, njobs {}".format(processID, nlast, nfirst, nrun, njobs))
  seqlist = range(1, njobs + 1)
  for seq in seqlist:
    nrun0 = nfirst + nrun * (seq - 1)
    nrun1tmp = nrun0 + nrun - 1

    if nrun1tmp <= nlast:
      nrun1 = nrun1tmp
    else:
      nrun1 = nlast

    print(processID, nrun0, nrun1, nrun, njobs)

    sub_flist = file_list[nrun0 - 1:nrun1]
    arg_flist = ",".join(sub_flist)

    if processID == 15271 or processID == 15275:
      log = f"./sublog/{processID}_{seq:03d}.log"
      subprocess.run(['bsub', '-q', 's', '-J', f'ana_{processID}_{seq:03d}', '-o', log, f'python3', 'runROOT.py', '--flist', arg_flist])
    else:
      subprocess.run(['bsub', '-q', 's', '-J', f'ana_{processID}_{seq:03d}', f'python3', 'runROOT.py', '--flist', arg_flist])
