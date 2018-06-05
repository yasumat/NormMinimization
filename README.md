# General Norm Approximation Solver

written by Yasuyuki Matsushita (yasumat@ist.osaka-u.ac.jp) at Osaka University.

# What is Norm Approximation?

In various scientific computing tasks, there arises a need for minimizing some vector norm, or a combination of different vector norms. For example, in compressive sensing, an unconstrained form of Lasso (least absolute shrinkage and selection operator) can be written as

![equation](https://latex.codecogs.com/gif.latex?%5Cmin_%7B%5Cmathbf%7Bx%7D%7D%20%5C%7C%5Cmathbf%7BA%7D%20%5Cmathbf%7Bx%7D%20-%5Cmathbf%7Bb%7D%5C%7C_2%5E2%20&plus;%20%5Clambda%20%5C%7C%5Cmathbf%7Bx%7D%5C%7C_1)

