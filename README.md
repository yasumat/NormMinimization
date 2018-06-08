# General Norm Approximation Solver

written by Yasuyuki Matsushita (yasumat@ist.osaka-u.ac.jp) at Osaka University.

# What is Norm Approximation?

In various scientific computing tasks, there arises a need for minimizing some vector norm, 
or a combination of different vector norms. A generalized norm approximation problem can be written as
<p align="center">
<img src="http://latex.codecogs.com/gif.latex?%5Cmin_%5Cmathbf%7Bx%7D%20%5Csum_%7Bk%3D1%7D%5EK%20%5Clambda_k%20%5C%7C%5Cmathbf%7BA%7D_k%20%5Cmathbf%7Bx%7D%20-%20%5Cmathbf%7Bb%7D_k%20%5C%7C_%7Bp_k%7D%5E%7Bp_k%7D">
</p>
where <img src="http://latex.codecogs.com/gif.latex?k%3D%5Cleft%5C%7B1%2C%20%5Cldots%2C%20K%5Cright%5C%7D"> is the term index,
 

For example, in compressive sensing, 
an unconstrained form of Lasso (least absolute shrinkage and selection operator) can be written as
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7Bx%7D%7D%20%5C%7C%5Cmathbf%7BA%7D%20%5Cmathbf%7Bx%7D%20-%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda%20%5C%7C%5Cmathbf%7Bx%7D%5C%7C_1">
</p>



