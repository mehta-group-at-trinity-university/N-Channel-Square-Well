
MCBound.x:	 Bsplines.f matrix_stuff.f modules_qd.o besselnew.o MCBound.o Bsplines.o matrix_stuff.o 
	gfortran  matrix_stuff.o Bsplines.o modules_qd.o  besselnew.o  -L/usr/local/opt/lapack/lib/ -llapack -lblas -L /Users/mehtan/Code/ARPACK/ARPACK/ -larpack_OSX -C MCBound.o -o MCBound.x

MCBound.o:	MCBound.f90
	gfortran -ffixed-line-length-132 -c -C MCBound.f90

matrix_stuff.o:	matrix_stuff.f
	gfortran -ffixed-line-length-132 -c -C matrix_stuff.f

Bsplines.o:	Bsplines.f
	gfortran -ffixed-line-length-132 -c -C Bsplines.f	

besselnew.o :	besselnew.f
	gfortran -ffixed-line-length-132 -c -C besselnew.f	

modules_qd.o :	modules_qd.f90
	gfortran -ffixed-line-length-132 -c -C modules_qd.f90	