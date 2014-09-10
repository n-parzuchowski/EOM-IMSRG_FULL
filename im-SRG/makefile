FC = gfortran $(FFLAGS)  

p1 = run_IMSRG
p2 = run_HO_energies

FFLAGS =  -O3 -fbounds-check -fopenmp
TFLAGS =   -g -fbounds-check

LIBS =  -L/user/local/lib/ -llapack -lblas

all: anglib.o basic_IMSRG.o HF_mod.o generators.o commutators.o adams_ode.o IMSRG_ODE.o IMSRG_MAGNUS.o main_IMSRG.o
	${FC} anglib.o basic_IMSRG.o HF_mod.o generators.o commutators.o adams_ode.o IMSRG_ODE.o IMSRG_MAGNUS.o main_IMSRG.o -o ${p1} ${LIBS}

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

adams_ode.o: adams_ode.f90
	${FC} -c adams_ode.f90 ${LIBS} 

IMSRG_ODE.o: IMSRG_ODE.f90
	${FC} -c IMSRG_ODE.f90 ${LIBS} 

IMSRG_MAGNUS.o: IMSRG_MAGNUS.f90
	${FC} -c IMSRG_MAGNUS.f90 ${LIBS} 


HO: add_sp_energies.f90
	${FC} add_sp_energies.f90 -o ${p2} 


clean:
	rm -f ${p1} $
	rm -f ${p2} $
	rm -f *.o
	rm -f *.mod
	rm -f *~

