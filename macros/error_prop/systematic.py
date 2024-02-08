import math

def systematic_error(err_stat, err_bkg):
  return math.sqrt(abs(err_stat**2 - err_bkg**2))

def main():

  # ud eLpR
  # values      = [8.53913e+01, -2.33598e+01] # S, A
  # error_stats = [1.34244e-01, 3.83239e-01] # S, A
  # error_bkgs  = [1.73475e-01, 4.90390e-01] # S, A

  # ud eRpL
  # values      = [2.94525e+01, -3.61040e+01] # S, A
  # error_stats = [7.13550e-02, 1.97118e-01] # S, A
  # error_bkgs  = [8.72346e-02, 2.44241e-01] # S, A

  # ss eLpR
  # values      = [1.84637e-02, 3.45007e-02] # S, A
  # error_stats = [8.83454e-05, 2.63519e-04] # S, A
  # error_bkgs  = [1.01899e-04, 2.77443e-04] # S, A

  # ss eRpL
  # values      = [1.83485e-02, 1.36361e-02] # S, A
  # error_stats = [1.65735e-04, 5.43070e-04] # S, A
  # error_bkgs  = [2.07100e-04, 5.81419e-04] # S, A

  ## AFB ##

  # ud eLpR
  # values      = [-0.0994781]
  # error_stats = [0.00011812]
  # error_bkgs  = [0.00948006]

  # ud eLpR p40
  values      = [-0.0940253]
  error_stats = [0.000280467]
  error_bkgs  = [0.00947992]

  # ud eRpL
  # values      = [-0.464087]
  # error_stats = [0.000176716]
  # error_bkgs  = [0.014481]

  # ss eLpR
  # values      = [0.702]
  # error_stats = [0.000100154]
  # error_bkgs  = [0.00167207]

  # ss eRpL
  # values      = [0.272]
  # error_stats = [0.000284787]
  # error_bkgs  = [0.00411347]

  # Calculate and print systematic errors for S and A
  for i in range(len(values)):
    sys_error = systematic_error(error_stats[i], error_bkgs[i])
    print(f'Systematic error for {values[i]}: {sys_error}')

if __name__ == "__main__":
  main()