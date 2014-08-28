FC = gfortran $(TFLAGS)  

p1 = run_IMSRG
p2 = run_HO_energies

FFLAGS =  -O3 -fbounds-check -fopenmp
TFLAGS =  -pg -g

LIBS =  -L/user/local/lib/ -llapack -lblas

all: anglib.o basic_IMSRG.o HF_mod.o generators.o commutators.o main_IMSRG.o
	${FC} anglib.o basic_IMSRG.o HF_mod.o generators.o commutators.o main_IMSRG.o -o ${p1} ${LIBS}

basic_IMSRG.o: basic_IMSRG.f90
	${FC} -c basic_IMSRG.f90 ${LIBS}

main_IMSRG.o:  main_IMSRG.f90
	${FC} -c  main_IMSRG.f90 ${LIBS}

anglib.o: anglib.f
	${FC} -c anglib.f ${LIBS}

HF_mod.o: HF_mod.f90
	${FC} -c HF_mod.f90 ${LIBS}

generators.o: generators.f90
	${FC} -c generators.f90 ${LIBS}

commutators.o: commutators.f90
	${FC} -c commutators.f90 ${LIBS}

HO: add_sp_energies.f90
	${FC} add_sp_energies.f90 -o ${p2} 


clean:
	rm -f ${p1} $
	rm -f ${p2} $
	rm -f *.o
	rm -f *.mod
	rm -f *~

