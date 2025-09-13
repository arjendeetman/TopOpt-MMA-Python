# TopOpt-MMA-Python

Example application of the [GCMMA-MMA-Python](https://github.com/arjendeetman/GCMMA-MMA-Python) library in topology optimization. The original toplogy optimization code is written by [Niels Aage and Villads Egede Johansen (Technical University of Denmark)](http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python). The python code is the equivalent of the efficient 88 lines MATLAB code with a minor difference: instead of a volume equality constraint it implements a volume inequality constraint that enforces the volume fraction to be smaller or equal to a preset volume fraction. The original python code can be downloaded [here](http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python). To use the modified Python code with the MMA optimizer the user needs to install the `mmapy` package with `pip install mmapy`. More information about the package can be found [here](https://github.com/arjendeetman/GCMMA-MMA-Python).

## Intended use

This Python code is based on prior work developed by [Niels Aage and Villads Egede Johansen](http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python), who made it accessible for use in engineering education. It is specifically tailored for students and newcomers to the field of structural optimization, offering a clear and accessible way to explore optimization concepts. The code may be freely used in academic courses and educational settings.

## References

Aage,  N.,  Johansen,  V.E.  (2013).  A  165 LINE  TOPOLOGY  OPTIMIZATION  CODE.  Retrieved  November  2,  2019  from 
http://www.topopt.mek.dtu.dk/Apps-and-software/Topology-optimization-codes-written-in-Python 

[Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., Sigmund, O. (2011). Efficient topology optimization in MATLAB 
using 88 lines of code. Structural and Multidisciplinary Optimization 43. 1-16. doi:10.1007/s00158-010-0594-7](https://link.springer.com/article/10.1007/s00158-010-0594-7)

[Liu, K., Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. Structural and Multidisciplinary Optimization 50. 1175–1196. doi:10.1007/s00158-014-1107-x](https://doi.org/10.1007/s00158-014-1107-x)

[Svanberg, K. (1987). The Method of Moving Asymptotes – A new method for structural optimization. International Journal 
for Numerical Methods in Engineering 24, 359-373. doi:10.1002/nme.1620240207](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620240207)

Svanberg, K. (n.d.). MMA and GCMMA – two methods for nonlinear optimization. Retrieved August 3, 2017 from  
https://people.kth.se/~krille/mmagcmma.pdf