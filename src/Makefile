
#choose some options.  by default, compiler optimization is turned off
OPT= -Dverbose  -g  -Wall # -Dprintpathlinprod  -O3 -funroll-loops -fexpensive-optimizations



###
# machine-specific configuration setup
###


#initialize a some variables
INSTALLPATH:=UNDEFINED
BERTINIHEADERS:=UNDEFINED
BERTINILIB:=UNDEFINED

  #check environment variable $(MACHINENAME), if defined
ifeq ($(MACHINENAME),MACBRAKE)
    #dan brake's machine.  he uses different places for his includes, etc.
    INC:=-I/usr/local/include  $(INC) 
    BISON=bison# bison in path already
    LIBS=-lm -lmpfr -lgmp# brake has libraries in path for linking, via brew.  no need to specify here
    COMP=gcc 
		COMPpp = clang++ -std=c++11 
    STATICFLAG=#  on macintosh, do not static compile (get missing lcrt0.o problem)
		BERTINILIB:=$(HOME)/lib/libbertini.a
	  INSTALLPATH:=$(HOME)/bin
		BERTINIHEADERS:=$(HOME)/bertini_source/include
endif

ifeq ($(MACHINENAME),wenrui)
    INC:= $(INC) -I/opt/mpfr/3.0.0/include -I/opt/gmp/5.0.1/include
    LIBS= -lm -L/usr/lib -L/opt/mpfr/3.0.0/lib -L/opt/gmp/5.0.1/lib /opt/mpfr/3.0.0/lib/libmpfr.a /opt/gmp/5.0.1/lib/libgmp.a
    COMP=gcc
    COMPpp=g++
    STATICFLAG=-static#no?
    BERTINILIB:=$(HOME)/lib/libbertini.a
    INSTALLPATH:=$(HOME)/bin
    BERTINIHEADERS:=/home/whao/bertini1.4/bertini/include
endif

ifeq ($(MACHINE),yourmachinenamehere)
    #use something like this to configure for your machine
    INC:=$(INC) #put your -I inclusions 
    BISON=bison# #path to bison
    LIBS=# #your libs
		COMP=gcc  #compiler
		COMPpp = g++
    STATICFLAG=-static#  need static?
endif

  #note: you can export an env variable of MACHINE to the string of your choice, and do configuration like above.  

###
# end machine specific configuration
###





#a few last definitions
INC:=$(INC) -I../include  -I$(BERTINIHEADERS) #this includes the header files for bertini_real, etc   

ARGS=$(INC) $(OPT) 
OBJFILES= data_type.o	programStartup.o determinant_derivative.o output.o  detjactodetjacSolver.o linprodtodetjacSolver.o lintolinSolver.o fileops.o isosingular.o curveNotSelfConj.o curveSelfConj.o witnessSetParse.o  postProcessing.o	checkSelfConjugate.o 
POBJFILES=lex.partitionParse.o 
FLEX=flex
#memory.o


###
#BEGIN MAKE RULES
###

#all comes first, so is the default rule when you `make'
all : bertini_real data_to_set sampler;    


#makes the bertini_real program
bertini_real : bertini_real.c $(OBJFILES) $(POBJFILES) ;
	$(COMP) $(ARGS) -o bertini_real bertini_real.c $(OBJFILES) $(POBJFILES) $(BERTINILIB)  $(LIBS)  $(STATICFLAG)


#makes the data2set program.
data_to_set: dataToSet.cpp
	$(COMPpp) $(ARGS) -o data2set dataToSet.cpp $(BERTINILIB) $(LIBS)  $(STATICFLAG)
#   

#makes the sampler program.
sampler: sampler.c $(OBJFILES) $(POBJFILES);
	$(COMP) $(ARGS) -o sampler sampler.c  $(OBJFILES) $(POBJFILES) $(BERTINILIB)  $(LIBS)  $(STATICFLAG)
#   


#simply copy the executables to the installpath
install: bertini_real data_to_set sampler
	cp bertini_real data2set sampler $(INSTALLPATH)


.c.o : 
	$(COMP) $(ARGS) -c $*.c

.cpp.o : 
	$(COMPpp) $(ARGS) -c $*.cpp

lex.partitionParse.o:	lex.partitionParse.c
	$(COMP) $(ARGS) -c lex.partitionParse.c

lex.partitionParse.c:	partitionParse.l
	$(FLEX) -PpartitionParse partitionParse.l

clean : 
	rm -f $(OBJFILES) $(POBJFILES) bertini_real data2set sampler
	rm -f par.out paramDerivs.out arr.out eval.out eval2.out const.out finalFile.out jacV.out jacP.out num.out
	rm -f tempBertiniFile lex.partitionParse.c 
	rm -rf *.dSYM

