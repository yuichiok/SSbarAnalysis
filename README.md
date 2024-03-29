# SSbarAnalysis (Light Quark Analysis)
[![GitHub commits](https://img.shields.io/github/last-commit/yuichiok/SSbarAnalysis)](https://GitHub.com/yuichiok/SSbarAnalysis/commit)
[![GitHub license](https://img.shields.io/github/license/yuichiok/SSbarAnalysis)](https://github.com/yuichiok/SSbarAnalysis/blob/main/LICENSE)
<!---![CMake](https://img.shields.io/badge/CMake-%23008FBA.svg?style=for-the-badge&logo=cmake&logoColor=white)--->

The primary objective of these codes was to analyze $s\bar{s}$ events, while the analysis is now extended to other light flavour quark events i.e. $u\bar{u}$ and $d\bar{d}$ processes.
This project expects input root file processed by [QQbarAnalysis](https://github.com/QQbarAnalysis/QQbarAnalysis) project.
List of parameters used by this analysis can be found in [TreeStructures.hh](https://github.com/yuichiok/SSbarAnalysis/blob/main/SSbarLibrary/include/TreeStructures.hh).

Contirbution to this project is always welcome, one just needs to follow [the contribution rules](https://github.com/yuichiok/SSbarAnalysis/blob/main/CONTRIBUTING.md).

## Requirements

We use [CMake](https://cmake.org) to build the project.
Minimum of requirements are:
|           | Version              |
|-----------|----------------------|
| CMake     | V3.0                 |
| ROOT      | V6.24/00             |

## Build
In order to build this project, you only need to execute following line.
```
./build.sh
```
If this gives an error, this may be caused by the version incompatibility from the building of ROOT.
Try changing the `CMAKE_CXX_STANDARD` version in the [CMakeLists.txt](https://github.com/yuichiok/SSbarAnalysis/blob/main/CMakeLists.txt#L12-L13).

## Execusion

After the build, execution file `main.exe` should be created.

### Configuration file

Configuration file determines the cut parameters for this analysis.
These are:
 - [**Legacy**] Types of generated particles. (currently use all types so don't care about this)
 - JET_btag_max, JET_ctag_max: b-tagging and c-tagging
 - JET_nvtx_max: Maximum number of jet vertex
 - PFO_TPCHits_min: Minimum number of TPC hits for a PFO
 - PFO_offset_max: Maximum offset distance for a PFO
 - PFO_p_min, PFO_p_max: Minimum and Maximum momentum for a PFO
 - LPFO_p_min, LPFO_p_max: Minimum and Maximum momentum for a PFO

Then you need to change the of `SSbarAnalysisConfig_example.ini` to `SSbarAnalysisConfig_default.ini`.
As a matter of fact, you can change the input file name in the [`main.cc`](https://github.com/yuichiok/SSbarAnalysis/blob/main/main.cc#L42).
In such case, you need to build the project again. (see above)

### Execute the file

One can execute this file by:
```
./main.exe [input_root_file]
```
replacing `[input_root_file]` by the path to the input root file from [QQbarAnalysis](https://github.com/QQbarAnalysis/QQbarAnalysis) project.
