
OPT=-Wall -g# -O3 -funroll-loops -fexpensive-optimizations  #-D_HAVE_MPI #-D_REACTIVE#-Wno-deprecated

INC:=-I../include #this includes the header files for bertini_real, etc   

#machine-specific configuration setup

  #check environment variable $(MACHINENAME), if defined
ifeq ($(MACHINENAME),MACBRAKE)
    #dan brake's machine.  he uses different places for his includes, etc.
    INC:=$(INC) 
    BISON=bison# bison in path already
    LIBS=# brake has libraries in path for linking, via brew.  no need to specify here
    COMP=gcc  
    STATICFLAG=#  on brake's macintosh, do not static compile (get missing lcrt0.o problem)

else ifeq ($(MACHINE),yourmachinenamehere)
    #use something like this to configure for your machine
    INC:=$(INC) #put your -I inclusions 
    BISON=bison# #path to bison
    LIBS=# #your libs
	COMP=gcc  #compiler
    STATICFLAG=-static#  need static?
	
else
    #default configuration.  for now, appears to be jon's computer.
    INC:=$(INC) -I/opt/mpfr/3.0.0/include -I/opt/gmp/5.0.1/include -I/opt/mpich2/include -I$(PETSC_DIR)/linux-gnu-c-debug/include -I$(PETSC_DIR)/include -I/usr/X11R6/include -I/bin/mpich2-install/include
    BISON=/home/jhauenst/bison/bison23/bin/bison
    LIBS=-L/opt/mpfr/3.0.0/lib -L/opt/gmp/5.0.1/lib /opt/mpfr/3.0.0/lib/libmpfr.a /opt/gmp/5.0.1/lib/libgmp.a
    COMP=gcc #-cc=gcc4 
    STATICFLAG=-static
endif

  #note: you can export an env variable of MACHINE to the string of your choice, and do configuration like above.  


ARGS=$(INC) $(OPT) 
OBJFILES=memory.o isosingular.o Curve_NotSelfConj.o Curve_SelfConj.o witnessSetParse.o programStartup.o checkSelfConjugate.o
POBJFILES=lex.partitionParse.o 
FLEX=flex


all : $(OBJFILES) $(POBJFILES) bertini_real ;

bertini_real : bertini_real.c $(OBJFILES) $(POBJFILES) ;
	$(COMP) $(ARGS) -o bertini_real bertini_real.c $(OBJFILES) $(POBJFILES) $(STATICFLAG) $(LIBS) -lm -lmpfr -lgmp

.c.o : 
	$(COMP) $(ARGS) -c $*.c 

lex.partitionParse.o:	lex.partitionParse.c
	$(COMP) $(ARGS) -c lex.partitionParse.c

lex.partitionParse.c:	partitionParse.l
	$(FLEX) -PpartitionParse partitionParse.l

clean : 
	rm -f $(OBJFILES) $(POBJFILES) bertini_real  
	rm -f par.out paramDerivs.out arr.out eval.out eval2.out const.out finalFile.out jacV.out jacP.out num.out
	rm -f tempBertiniFile lex.partitionParse.c 
	rm -rf *.dSYM

