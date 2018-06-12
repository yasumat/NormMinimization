# General Norm Approximation Solver

written by Yasuyuki Matsushita (yasumat@ist.osaka-u.ac.jp) at Osaka University.

### What is Norm Approximation?

In various scientific computing tasks, there arises a need for minimizing some vector norm, 
or a combination of different vector norms. A generalized norm approximation problem can be written as

<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5Csum_%7Bk%3D1%7D%5EK%20%5Clambda_k%20%5C%7C%5Cmathbf%7BA%7D_k%20%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D_k%20%5C%7C_%7Bp_k%7D%5E%7Bp_k%7D"/>

where <img src="http://latex.codecogs.com/png.latex?k%3D%5C%7B1%2C%20%5Cldots%2C%20K%5C%7D"> is the term index, 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D_k%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bm_k%20%5Ctimes%20n%7D"/>
and 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bb%7D_k%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bm_k%7D"/>
are the design matrix and vector that define the 
<img src="http://latex.codecogs.com/gif.latex?k"/>-th objective. 
The overall objective function is defined as a linear combination of 
<img src="http://latex.codecogs.com/gif.latex?p_k"/>-th power of 
<img src="http://latex.codecogs.com/gif.latex?%5Cell_%7Bp_k%7D"/>-norm weighted by 
<img src="http://latex.codecogs.com/gif.latex?%5Clambda_k"/>. The goal is to determine 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bx%7D%20%5Cin%20%5Cmathbb%7BR%7D%5En"/>. 
This software package provides a solution method for this problem.


For example, in compressive sensing, 
an unconstrained form of Lasso (least absolute shrinkage and selection operator) can be written as

<img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7Bx%7D%7D%20%5C%7C%5Cmathbf%7BA%7D%20%5Cmathbf%7Bx%7D%20-%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda%20%5C%7C%5Cmathbf%7Bx%7D%5C%7C_1"/>

which corresponds to the case where
<img src="http://latex.codecogs.com/gif.latex?k%3D2"/> and
<img src="http://latex.codecogs.com/gif.latex?%5Cleft%28%5Cmathbf%7BA%7D_1%2C%20%5Cmathbf%7BA%7D_2%2C%20%5Cmathbf%7Bb%7D_1%2C%20%5Cmathbf%7Bb%7D_2%2C%20%5Clambda_1%2C%20%5Clambda_2%2C%20p_1%2C%20p_2%20%5Cright%29%3D%5Cleft%28%5Cmathbf%7BA%7D%2C%20%5Cmathbf%7BI%7D%2C%20%5Cmathbf%7Bb%7D%2C%20%5Cmathbf%7B0%7D%2C%201%2C%20%5Clambda%2C%202%2C%201%5Cright%29"/>
 in the general form described above.
 
As another example, Tikhonov regularization or ridge regression 
 
<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5C%7C%5Cmathbf%7BA%7D%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda%20%5C%7C%5Cmathbf%7B%5CGamma%7D%20%5Cmathbf%7Bx%7D%20%5C%7C_2%5E2"/>
 
that appears in image restoration, super-resolution, and image deblurring. This case corresponds to the case where
<img src="http://latex.codecogs.com/gif.latex?k%3D2"/> and 
<img src="http://latex.codecogs.com/png.latex?%5Cleft%28%5Cmathbf%7BA%7D_1%2C%20%5Cmathbf%7BA%7D_2%2C%20%5Cmathbf%7Bb%7D_1%2C%20%5Cmathbf%7Bb%7D_2%2C%20%5Clambda_1%2C%20%5Clambda_2%2C%20p_1%2C%20p_2%20%5Cright%20%29%3D%5Cleft%28%5Cmathbf%7BA%7D%2C%20%5Cmathbf%7B%5CGamma%7D%2C%20%5Cmathbf%7Bb%7D%2C%20%5Cmathbf%7B0%7D%2C%201%2C%20%5Clambda%2C%202%2C%202%20%5Cright%20%29"/>.
 
These objective functions can be further augmented by additional 
<img src=""http://latex.codecogs.com/gif.latex?%5Cell_%7Bp_k%7D""/>-norm terms that represent constraints.
For example, elastic net is defined as
 
<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5C%7C%5Cmathbf%7BA%7D%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Calpha%20%5C%7C%20%5Cmathbf%7Bx%7D%20%5C%7C_1%20&plus;%20%5Cbeta%20%5C%7C%5Cmathbf%7Bx%7D%5C%7C_2%5E2"/>
 
which corresponds to  <img src="http://latex.codecogs.com/gif.latex?k%3D3"/> and 
<img src="http://latex.codecogs.com/png.latex?%5Cleft%28%5Cmathbf%7BA%7D_1%2C%20%5Cmathbf%7BA%7D_2%2C%20%5Cmathbf%7BA%7D_3%2C%20%5Cmathbf%7Bb%7D_1%2C%20%5Cmathbf%7Bb%7D_2%2C%20%5Cmathbf%7Bb%7D_3%2C%20%5Clambda_1%2C%20%5Clambda_2%2C%20%5Clambda_3%2C%20p_1%2C%20p_2%2C%20p_3%20%5Cright%20%29%3D%5Cleft%28%5Cmathbf%7BA%7D%2C%20%5Cmathbf%7BI%7D%2C%20%5Cmathbf%7BI%7D%2C%20%5Cmathbf%7Bb%7D%2C%20%5Cmathbf%7B0%7D%2C%20%5Cmathbf%7B0%7D%2C%201%2C%20%5Calpha%2C%20%5Cbeta%2C%202%2C%201%2C%202%20%5Cright%20%29"/>.
 
### How to use?
 
 
### Conditions of use

This package is distributed under the GNU General Public License. For
information on commercial licensing, please contact the authors at the
contact address below. If you use this code for a publication, please cite the following paper:

    @inproceedings{FGNA2016,
	  	title={Fast General Norm Approximation via Iteratively Reweighted Least Squares},
	  	author={Masaki Samejima and Yasuyuki Matsushita},
	  	booktitle={Proceedings of Asian Conference on Computer Vision Workshops (ACCVW)},
	  	year={2016}
	}
 
### Acknowledgements




### Contact information

Questions / Comments? Bug reports? Please contact Yasuyuki Matsushita at yasumat@ist.osaka-u.ac.jp.







