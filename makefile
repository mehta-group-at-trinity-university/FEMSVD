4BodySVD.x:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90
	ifort -O4 4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90 -mkl -L/opt/ARPACK -larpack_Intel -parallel -no-wrap-margin 
#-debug full -traceback -check bounds
#-qopt-report=3 
#-debug full -traceback -check bounds


#4BodySVD.exe:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90
#	ifort 4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90 ARPACK.lib -Qmkl /wrap-margin- /traceback /check:all /gen-interfaces /warn:interfaces /fpe:0 /debug:full

clean:
	rm *.optrpt *.exe *.out *.o *.obj *.un~ *.pdb *.mod *__*.f90

