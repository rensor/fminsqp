# fminsqp
Matlab based optimizer framework using Sequential Quadratic Programming (SQP)

The implementation is based on the fminslp framework. However, instead of relying on Sequential Linear Programming (SLP), the fminsqp framework relies on Sequential Quadratic Programming (SQP). 
In order to ensure stable convergence, a global convergence filter by [1] is applied. The hessian of the objective function can be approximated by two different methods i.e., BFGS or DFP.

[1] Fletcher R, Leyffer S, Toint PL (2000): On the global convergence of a filter-sqp algorithm. Numerical Analysis Report NA/197,
Department of Mathematics, University of Dundee, Scotland, UK

More documentation to come.
