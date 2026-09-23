
LoopOverCalcBound.x:	LoopOverCalcBound.f Bsplines.f matrix_stuff.f modules_qd.o besselnew.o LoopOverCalcBound.o Bsplines.o matrix_stuff.o ~/bin/bspline90_22.o
	g95 LoopOverCalcBound.o matrix_stuff.o Bsplines.o modules_qd.o ~/bin/bspline90_22.o besselnew.o -framework vecLib -L/opt/local/lib/ -larpack -C -o LoopOverCalcBound.x
#	g95 LoopOverCalcBound.o matrix_stuff.o Bsplines.o 1Dpot.o modules_qd.o bspline90_22.o besselnew.o  -framework vecLib -larpack_OSX -lg2c -C -o LoopOverCalcBound.x

#1Dpot.o:	1Dpot.f
#	g95 -ffixed-line-length-132 -c -C 1Dpot.f

LoopOverCalcBound.o:	LoopOverCalcBound.f
	g95 -ftrace=full -ffixed-line-length-132 -c -C LoopOverCalcBound.f

matrix_stuff.o:	matrix_stuff.f
	g95 -ffixed-line-length-132 -c -C matrix_stuff.f

Bsplines.o:	Bsplines.f
	g95 -ftrace=full -ffixed-line-length-132 -c -C Bsplines.f	

besselnew.o :	besselnew.f
	g95 -ffixed-line-length-132 -c -C besselnew.f	

modules_qd.o :	modules_qd.f90
	g95 -ffixed-line-length-132 -c -C modules_qd.f90	