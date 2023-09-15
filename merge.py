import os
import sys
import subprocess

# udsProcess
udsProcess = "uu"

# chiral
# chiral = "eLpR"
chiral = "eRpL"

# processID
chiral_node = ""
processIDs = []
if chiral == "eLpR":
  chiral_node = "eL.pR"
  processIDs = [15162,15271,15273,15315,15319,15351,15353,15355,15357,15463,15465,15467,15469,15471,15473,15475,15477,15479,15499,15501,15503,15505,15507,15509,15511,15513,15515,15517,15539,15541,15549,15551,15553,15555,15557,15559,15561,15563,15585,15587,15589,15591,15593,15595,15597,15599,15601,15603,15629,15633]
elif chiral == "eRpL":
  chiral_node = "eR.pL"
  processIDs = [15165,15275,15277,15279,15359,15361,15363,15365,15367,15369,15481,15483,15485,15487,15489,15491,15493,15495,15497,15519,15521,15523,15525,15527,15529,15531,15533,15535,15537,15565,15567,15569,15571,15573,15575,15577,15579,15581,15583,15605,15607,15609,15611,15613,15616,15618,15620,15625,15627,15631]

indir    = "./rootfiles/tmp_root/"
outdir   = "./rootfiles/merged/"
stagedir = f"{outdir}/staged/"
dirCheck = os.path.isdir(indir) or os.path.isdir(outdir) or os.path.isdir(stagedir)
if not dirCheck:
  sys.exit("Error: directory not found")

subprocess.run(f"rm -rf {stagedir}/*",shell=True)

for processID in processIDs:
  print(f"Merging {processID}")
  subprocess.run(f"hadd -f -j 8 {stagedir}/stage_{processID}.root {indir}/*{processID}*.root",shell=True)

# subprocess.run(f"rm -rf {indir}",shell=True)
# subprocess.run(f"mkdir -p {indir}",shell=True)

haddCommand = f"hadd -f -j 8 {outdir}/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.{chiral_node}.{udsProcess}.KPiLPFO.dedxPi.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root {stagedir}/stage_*.root"
print(haddCommand)
subprocess.run(haddCommand,shell=True)
