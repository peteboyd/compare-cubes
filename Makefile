FC = gfortran
FFLAGS = -O3 
all: compare-cubes

compare-cubes: 
	$(FC) $(FFLAGS) compare_cubes.f -o compare-cubes

clean:
	rm -f compare-cubes
