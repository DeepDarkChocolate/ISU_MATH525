# ISU_MATH525
 2021 spring ISU MATH525

This program requires a lot of statistical computations like the followings: solving linear equations, Obtaining a Gauss-Hermite quadrature for numerical integration, generating Normal random variates, computing the probability density functions for Normal random variables, and sampling packages. For the first computation, I used LaPACK routines and for the three next ones, I used open license codes available in https://people.sc.fsu.edu/~jburkardt/f_src/f_src.html that are publicly open(pdflib.f90,  rnglib.f90, and hermite\_rule.f90). For the last one(sampling), I used my own code by generating Uniform random variable from \textit{rnglib.f90}. In order to use these library routines, we compile the program by

mpiifort pdflib.f90 rnglib.f90 hermite_rule.f90 
ftns.f90 main.f90 -llapack
