tales: modtypes.o modhelpers.o modsample.o tales.o
	gfortran -o tales modtypes.o modhelpers.o modsample.o tales.o

modtypes.o: modtypes.f90
	gfortran -c $<

modtypes.mod: modtypes.f90 modtypes.o
	@true

modhelpers.o: modhelpers.f90
	gfortran -c $<

modhelpers.mod: modhelpers.f90 modhelpers.o
	@true

modsample.o: modsample.f90 modtypes.mod modhelpers.mod
	gfortran -c $<

modsample.mod: modsample.f90 modsample.o
	@true

tales.o: tales.f90 modsample.mod
	 gfortran -c $<

clean:
	rm *.o *.mod tales
