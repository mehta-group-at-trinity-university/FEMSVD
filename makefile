CMP     = gfortran
OPTFLAG = -O3
F132FORM = -ffree-form
LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include

4BodySVD.x:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90
	${CMP} ${OPT} ${F132FORM} 4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90 ${LAPACK} ${ARPACK} ${INCLUDE} 
#-debug full -traceback -check bounds
#-qopt-report=3 
#-debug full -traceback -check bounds


#4BodySVD.exe:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90
#	ifort 4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90 ARPACK.lib -Qmkl /wrap-margin- /traceback /check:all /gen-interfaces /warn:interfaces /fpe:0 /debug:full

