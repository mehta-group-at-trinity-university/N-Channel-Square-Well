FC = gfortran
FFLAGS = -fdefault-real-8 -fdefault-double-8
LAPACKDIR = /opt/homebrew/opt/openblas/lib
ARPACKDIR = /opt/homebrew/opt/arpack/lib
LIBS = -L$(LAPACKDIR) -llapack -lblas -L$(ARPACKDIR) -larpack

# squarenew.f90 defines modules (modb, potmod, DataStructures) that besselnew.f
# and square_routines.f depend on via `use`; it must be compiled first.
# minpack.f90 is linked by dos.f90 in a sibling program, but squarenew.f90
# itself never calls into it, so it is intentionally omitted here. ARPACK is
# still needed at link time: matrix_stuff.f's MyDsband/dsband (unused by
# squarenew.f90, called only from dos.f90) reference dsaupd_/dseupd_.

ssnew.x: squarenew.o besselnew.o square_routines.o GOE.o rgnf_lux.o matrix_stuff.o
	$(FC) squarenew.o GOE.o rgnf_lux.o matrix_stuff.o besselnew.o square_routines.o $(FFLAGS) -o ssnew.x $(LIBS)

squarenew.o: squarenew.f90
	$(FC) -fcheck-new -fcheck=bounds -Wargument-mismatch -Winteger-division -Wsurprising -Wintrinsic-shadow $(FFLAGS) -c squarenew.f90

besselnew.o: besselnew.f squarenew.o
	$(FC) $(FFLAGS) -ffixed-line-length-none -fno-range-check -c besselnew.f

square_routines.o: square_routines.f squarenew.o
	$(FC) $(FFLAGS) -c square_routines.f

GOE.o: GOE.f
	$(FC) -ffixed-line-length-132 -c GOE.f

rgnf_lux.o: rgnf_lux.f
	$(FC) -ffixed-line-length-132 -c rgnf_lux.f

matrix_stuff.o: matrix_stuff.f
	$(FC) -ffixed-line-length-132 -c matrix_stuff.f

clean:
	rm -f *.o *.mod ssnew.x
