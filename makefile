FC = gfortran $(FFLAGS)  

p1 = run_IMSRG

FFLAGS =  -O3 -fbounds-check -fopenmp
TFLAGS =  -g  

LIBS =  -L/user/local/lib/ -llapack -lblas

all: basic_IMSRG.o main_IMSRG.o
	${FC} basic_IMSRG.o main_IMSRG.o -o ${p1} ${LIBS}

basic_IMSRG.o: basic_IMSRG.f90
	${FC} -c basic_IMSRG.f90 ${LIBS}

main_IMSRG.o:  main_IMSRG.f90
	${FC} -c  main_IMSRG.f90 ${LIBS}



clean:
	rm -f ${p1} $
	rm -f *.o
	rm -f *.mod
	rm -f *~

