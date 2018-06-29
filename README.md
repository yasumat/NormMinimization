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
 
As another example, Tikhonov regularization or ridge regression that appears in image restoration, super-resolution, 
and image deblurring can be written as
 
<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5C%7C%5Cmathbf%7BA%7D%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda%20%5C%7C%5Cmathbf%7B%5CGamma%7D%20%5Cmathbf%7Bx%7D%20%5C%7C_2%5E2"/>.
 
This case corresponds to the case where
<img src="http://latex.codecogs.com/gif.latex?k%3D2"/> and 
<img src="http://latex.codecogs.com/png.latex?%5Cleft%28%5Cmathbf%7BA%7D_1%2C%20%5Cmathbf%7BA%7D_2%2C%20%5Cmathbf%7Bb%7D_1%2C%20%5Cmathbf%7Bb%7D_2%2C%20%5Clambda_1%2C%20%5Clambda_2%2C%20p_1%2C%20p_2%20%5Cright%20%29%3D%5Cleft%28%5Cmathbf%7BA%7D%2C%20%5Cmathbf%7B%5CGamma%7D%2C%20%5Cmathbf%7Bb%7D%2C%20%5Cmathbf%7B0%7D%2C%201%2C%20%5Clambda%2C%202%2C%202%20%5Cright%20%29"/>.
 
These objective functions can be further augmented by additional 
<img src="http://latex.codecogs.com/gif.latex?%5Cell_%7Bp_k%7D"/>-norm terms that represent constraints.
For example, elastic net is defined as
 
<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5C%7C%5Cmathbf%7BA%7D%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Calpha%20%5C%7C%20%5Cmathbf%7Bx%7D%20%5C%7C_1%20&plus;%20%5Cbeta%20%5C%7C%5Cmathbf%7Bx%7D%5C%7C_2%5E2"/>
 
which corresponds to  <img src="http://latex.codecogs.com/gif.latex?k%3D3"/> and 
<img src="http://latex.codecogs.com/png.latex?%5Cleft%28%5Cmathbf%7BA%7D_1%2C%20%5Cmathbf%7BA%7D_2%2C%20%5Cmathbf%7BA%7D_3%2C%20%5Cmathbf%7Bb%7D_1%2C%20%5Cmathbf%7Bb%7D_2%2C%20%5Cmathbf%7Bb%7D_3%2C%20%5Clambda_1%2C%20%5Clambda_2%2C%20%5Clambda_3%2C%20p_1%2C%20p_2%2C%20p_3%20%5Cright%20%29%3D%5Cleft%28%5Cmathbf%7BA%7D%2C%20%5Cmathbf%7BI%7D%2C%20%5Cmathbf%7BI%7D%2C%20%5Cmathbf%7Bb%7D%2C%20%5Cmathbf%7B0%7D%2C%20%5Cmathbf%7B0%7D%2C%201%2C%20%5Calpha%2C%20%5Cbeta%2C%202%2C%201%2C%202%20%5Cright%20%29"/>.

This implementation can take arbitrary numbers of norm terms, each of them is linear but can have different norm, constraint matrix, and vector. 
### How to use?

Call a function ``solve`` defined in ``normapprox.py`` by appropriately forming matrices 
<img src="http://latex.codecogs.com/gif.latex?%5C%7B%5Cmathbf%7BA%7D_k%5C%7D"/> and 
vectors <img src="http://latex.codecogs.com/gif.latex?%5C%7B%5Cmathbf%7Bb%7D_k%5C%7D"/> as well as 
lists of weights <img src="http://latex.codecogs.com/gif.latex?%5C%7B%5Clambda_k%5C%7D"/>
and norm specifiers <img src="http://latex.codecogs.com/gif.latex?%5C%7Bp_k%5C%7D"/> in the form of python lists, and 
pass them to ``A_list``, ``b_list``, ``lambda_list``, and ``p_list``, respectively. 

```
def solve(A_list=None, b_list=None, lambda_list=None, p_list=None, max_iter=10000, tol=1.0e-8):
```

See examples for more.

### Example 1: Polynomial fitting with regularization - ex01.py

We take an example of polynomial fitting. Suppose we have data points in the x-y plane and wish to fit 
a polynomial function of degree <img src="http://latex.codecogs.com/gif.latex?d"/> to the data points:

<img src="http://latex.codecogs.com/gif.latex?y%20%3D%20w_1%20x%5Ed%20&plus;%20w_2%20x%5E%7Bd-1%7D%20&plus;%20...%20&plus;%20w_d%20x%20&plus;%20x_%7Bd&plus;1%7D"/>

Here parameters to estimate are coefficients 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bw%7D%20%3D%20%5C%7Bw_1%2C%20...%2C%20w_%7Bd&plus;1%7D%5C%7D"/>.
This problem can be written in a matrix form by a proper matrix 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D"/> and a vector
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bb%7D"/>, as 
<img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D%20%5Cmathbf%7Bw%7D%20%5Capprox%20%5Cmathbf%7Bb%7D"/> 
(See ``ex01.py`` for how <img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7BA%7D"/> 
and <img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bb%7D"/> are formed).
The goal is to estimate a *good* <img src="http://latex.codecogs.com/gif.latex?%5Cmathbf%7Bw%7D"/> that best approximates the equality.
For example, with an L2 metric, this can be formulated as:

Conventional least-squares regression (Case 1): 
<img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7Bw%7D%7D%20%5C%7C%5Cmathbf%7BA%7D%20%5Cmathbf%7Bw%7D%20-%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20"/>

<BR>
 
Over-fitting is a common issue in this type of problems when a proper polynomial degree is unkonwn. 
A conventional approach to avoiding this issue is to *regularize* the esimates by some additional constraints. 
This example penalizes the estimates that have high values using a *regularizer*:

Least-squares regression with a regularizer (Case 2):
<img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7Bw%7D%7D%20%5C%7C%5Cmathbf%7BA%7D%20%5Cmathbf%7Bw%7D%20-%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda_2%20%5C%7C%5Cmathbf%7Bw%7D%5C%7C_2^2"/>,
where 
<img src="http://latex.codecogs.com/gif.latex?%5Clambda_2%20%3D%201.0%5Cmathrm%7Be%7D%5E%7B-4%7D"/> in this particular example.

<BR>

These problems Case 1 and Case 2 can be solved by the following code:


    import normapprox as na
    ...
    # Case 1: Conventional least-squares fitting
    w1, residue1, ite1 = na.solve(A_list=[A], b_list=[b], lambda_list=[1], p_list=[2])
    # Case 2: Least-squares fitting with L2 regularization (special form of Ridge regression)
    w2, residue2, ite2 = na.solve(A_list=[A, np.identity(n)], b_list=[b, np.zeros(n)],
                                  lambda_list=[1, 1e-4], p_list=[2, 2])


Below shows the result of Case 1 and Case 2 fittings. With a regularizer (Case 2), it can be seen that the issue of
over-fitting is alleviated.

<img src="./ex01.png"/> 
 
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

### Dependencies
The code is written in Python 3.6 but should be able to adapt it to Python 2.x if needed.
You might need the following Python packages installed:
* `numpy` (main computation depends on matrix operations)
* `matplotlib` (for running example codes, but not used in the main computation code)


### Contact information

Questions / Comments? Bug reports? Please contact Yasuyuki Matsushita at yasumat@ist.osaka-u.ac.jp.







