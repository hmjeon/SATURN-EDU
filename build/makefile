#
# Makefiles for SATURN-EDU, written by Hyungmin Jun
#
#	-fast
#	Maximizes speed across the entire program.
#	Description:
#	This option maximizes speed across the entire program.
#
#	It sets the following options:
#	- On Linux* systems: -ipo, -O3, -no-prec-div, -static, and -xHost
#	- On OS X* systems: -ipo, -mdynamic-no-pic, -O3, -no-prec-div, and -xHost
#	- On Windows*  systems:  /O3,  /Qipo,  /Qprec-div-,  and /QxHost
#
#	NOTE: Option fast sets some aggressive optimizations that may not be
#   appropriate for all applications. The resulting executable may not run on
#   processor types different from the one on which you compile. You should
#   make sure that you understand the individual optimization options
#   that are enabled by option fast.

objects = Data_Struct.o DISP4.o Solver.o SATURN.o

f90comp = ifort

# Makefile
SATURN-EDU: $(objects)
	$(f90comp) -o SATURN-EDU $(objects)

Data_Struct.mod: Data_Struct.o Data_Struct.f90
	$(f90comp) -c Data_Struct.f90
Data_Geom.o: Data_Struct.f90
	$(f90comp) -c Data_Struct.f90
	
DISP4.mod: DISP4.o DISP4.f90
	$(f90comp) -c DISP4.f90
DISP4.o: DISP4.f90
	$(f90comp) -c DISP4.f90

Solver.mod: Solver.o Solver.f90
	$(f90comp) -c Solver.f90
Data_Mesh.o: Solver.f90
	$(f90comp) -c Solver.f90

SATURN.mod: SATURN.o SATURN.f90
	$(f90comp) -c SATURN.f90
Data_Prob.o: SATURN.f90
	$(f90comp) -c SATURN.f90

# Cleaning everything
clean :
	rm *.o
	rm *.mod
# End of the makefile
