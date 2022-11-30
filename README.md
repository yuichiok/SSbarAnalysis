# SSbarAnalysis
[![GitHub commits](https://img.shields.io/github/last-commit/yuichiok/SSbarAnalysis)](https://GitHub.com/yuichiok/SSbarAnalysis/commit)
[![GitHub license](https://img.shields.io/github/license/yuichiok/SSbarAnalysis)](https://github.com/yuichiok/SSbarAnalysis/blob/main/LICENSE)
<!---![CMake](https://img.shields.io/badge/CMake-%23008FBA.svg?style=for-the-badge&logo=cmake&logoColor=white)--->

Analysis dedicated to $e^+e^- \rightarrow s\bar{s}$ process at Inernational Linear Collider (ILC) at 250 GeV running scenario.
The primary objective is to identify $s$ and $\bar{s}$ quarks out from among the other quark flavor backgrounds, using the leading Kaons as its signature.

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
Then you need to change the of `SSbarAnalysisConfig_example.ini` to `SSbarAnalysisConfig_default.ini`.
As a matter of fact, you can change the input file name in the [`main.cc`](https://github.com/yuichiok/SSbarAnalysis/blob/main/main.cc#L42).
In such case, you need to build the project again. (see above)

### Configuration file

These are the initialization parameters for running the analysis.
Details will be explained in the future commit for this README.md file.

One can execute this file by:
```
./main.exe [input_root_file]
```
replacing `[input_root_file]` by the path to the input root file.
