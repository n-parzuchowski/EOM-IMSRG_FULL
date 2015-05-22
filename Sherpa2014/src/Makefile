COMP = gfortran
OPTIONS = -O2 -fopenmp  # -r8 

OBJ1 = \
modulesSHERPA.o\
sherpamain.o \
dla_invert.o \
eiglib.o \
GChamlibpn.o \
GClib.o \
HFdiag.o \
HFmenu.o \
libgrd.o \
libra.o \
makeSD.o \
mat_subs.o \
quad.o \
ranlib.o \
LinAlgLib.o \
vectorslib.o \
prjgrdlib.o \
RPAmenu.o \
librpa.o \
transLib.o \
rg.o \
OpCorr.o \
librpaPN.o \
pnTrans.o

.f.o:;  ($(COMP) -c $(OPTIONS) $(<F)) 


sherpa14.x: $(OBJ1) 
	$(COMP) -o sherpa14.x $(OPTIONS) $(OBJ1) 

clean:
	rm -r *.o
	rm -r *.mod
	rm -r *.x
