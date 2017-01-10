# 2DOrthogonal_Polynomial_Decomposition
OPD Matlab implementation with application to a IR sequence of a simulated CFRP sample and two real CFRP and GFRP samples with teflon insertions. Recreates the results presented in Paper "Characterization of defects of pulsed thermography inspections by orthogonal polynomial decomposition"

Because of the large amount of space required by the thermal sequences, and because github does not support big data files, the link to download IR sequences (in OneDrive) is: https://1drv.ms/f/s!AgtrBjgDiWxVhkikDqqTeU-d76eJ

The following files allow reproducing the results:
  1. mainSyntheticCFRP.m and data file syntheticCFRP.mat
  2. mainCFRPsample.m and data file CFRPsmaple.mat
  3. mainGFRPsample.m and data file GFRPsample.mat

The firrst m-file reproduces Figures 9-12 and Table 1; the second m-fle Figures 15-18 and Table 2; the last m-file
Figures 19-22 and Table 3.

Consider that the mat-fles and the m-fles must be in the same directory as well the functions folder. That
directory can be anywhere in your machine, you only need to add the path when Matlab requests it (by clicking
the Add to Path button after clicking the Run button).  You need the signal processing toolbox in order to run the m-files without problems (findpeaks() and resample() functions from the toolbox are needed).




