SBD (Sparse Bessel Decomposition) toolbox
Developped by François Alouges, Matthieu Aussal and Martin AVERSENG
Copyright Greengard for the NUFFT algo. 
Contains functions of the MyBEM toolbox developped by Matthieu Aussal and François Alouges. 
Based on the ideas of the SCSD method by Alouges and Aussal 
François Alouges and Matthieu Aussal. 2015. The sparse cardinal sine decomposition and its application for fast numerical convolution. 
Numer. Algorithms 70, 2 (October 2015), 427-448. DOI=http://dx.doi.org/10.1007/s11075-014-9953-6 

Compute operators \tilde{A} that approximate 
f -> A f
where A is a matrix of the form 
A(i,j) = G(|x(i) - y(j)|)

How to : 
Given two clouds of points in \mathbb{R}^2, X and Y of size [Nx,2] and [Ny,2]

- First create an instance 'kernel' of a Kernel object, containing the info about G. You have several possibilities : 
 	* Either G is x -> log(R*x) (for some positive R),  x -> J0(R*x)  or x -> Y0(R*x) then use the specific subclasses of the class Kernel 
	for example, for log, kernel = LogKernel(R). 
	* Either G is an arbitrary function, then use kernel = Kernel(func,der) where func = @(x)(G(x)) and der = @(x)(G'(x)). 
	* In this case, you can alse improve the speed of the code by creating a specific kernel object and creating a method 
	for the dilatation operator, setting as the scalFunc property the function a,b,rho -> \int_{a}^b G'(x) J1(rho*x)dx (J1 
	bessel function of first kind order 1) if some explicit form is known, and setting as the normFunc property the 
	function a,b -> \int_{a}^b x G'(x)^2 dx. You can use the LogKernel subclass as a model. 
- Then choose a parameter 'a'. This parameter is the cutoff that splits the interactions into close and far. The optimal value for 'a' depends
on the data. To help choosing the value, try to have the number of close interactions (i.e. number of pairs (i,j) such that |x(i) - y(j)|<a Rmax, 
where Rmax is the diameter of the union of the two clouds) scaling as 1/a^2. 
- The parameter 'b' should be always set to 1, for standard cases (G defined as one of the specific kernels). If the radial quadrature of the kernel 
has a lot of terms or doesn't converge, try setting b = 1-a, although stability is not guaranteed. 
- Choose a tolerance 'tol' (for example 1e-3)

Then you can create the operator by calling A = Op(X,Y,kernel,a,b,tol);

You can also add Op(...,'noCloseField',true) when you have two distant clouds X and Y (and you know in advance that the set of close interactions is empty) 
This will save a huge time of assembling since no range search will be performed. 

Then you can apply the operator on a vector f of dimensions (Ny,1) by doing 
q = A*f

you can also validate a subset of the entries q in a time tVal by calling 
err = A.validate(f,q,tVal,tol)

This will compute as many entries of the vector 'q' as possible doing the full product during the time tVal and return the maximal error detected. 
It will issue a warning if an error greater than tol was detected. Before validating, make sure the 'f' vector has a unit 'l^1' norm. 

you can check the performance of the operator A by calling 
A.showPerf;

This displays the memory load (and its repartition in close and far contribution). You can use this info to change 'a' for better performance. It also shows the time spent assembling, 
the time of a Matrix-Vector product, and an estimation of the speed-up due to the method. 
It alsos creates a graph representing the radial quadrature. The left panel shows the radial function r -> G(r) (blue) and the Fourier-Bessel approximation (red dashed)
The middle panel represents the error in log log scale. A vertical dashed line shows the position of parameter 'a' on the r axis, and a horizontal dashed line shows the user-defined
tolerance. If everything went well, the curve should be outside the top right section of the graph. (error(r > a) < tol). The right panel compares the weights derived by the algorithm
(blue circles) with those of the truncated Fourier-Bessel series if the code was able to estimate them (red circles) . You should observe a faster decay for the blue spectrum. 

The radial quadrature related functions are located in the folder SBD/Op/RadialQuad
The 2D quadrature related functions are located in the folder SBD/Op/Quad2D

The code that creates the close and far operator is located in the methods quad2D and B of the object Kernel. 




