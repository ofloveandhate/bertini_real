
#choose some options.  by default, compiler optimization is turned off
OPT= -Dverbose  -g  -Wall #-Wdeprecated # -Dprintpathlinprod  -O3 -funroll-loops -fexpensive-optimizations






###
# machine-specific configuration setup
###


#initialize a some variables
INSTALLPATH:=UNDEFINEDinstallpath
BERTINIHEADERS:=UNDEFINEDbertiniheaders
BERTINILIB:=UNDEFINEDbertinilib
COMPcc:=UNDEFINEDcompcc
COMPpp:=UNDEFINEDcomppp
BOOSTLOC := UNDEFINEDboostloc
BOOSTSUFFIX:=UNDEFINEDboostsuffix
MACHINE_SPECIFICS:= 

  #check environment variable $(MACHINENAME), if defined
ifeq ($(MACHINENAME),MACBRAKE)
    #dan brake's machine.  he uses different places for his includes, etc.
    INC:= -I/usr/local/include  -I/usr/local/Cellar/mpich2/3.0.4/include $(INC) 
    BISON=bison# bison in path already
    LIBS=-lm -lmpfr -lgmp# brake has libraries in path for linking, via brew.  no need to specify here
    COMPcc=mpicc -cc=clang 
	COMPpp = mpic++ -cxx=clang++
    STATICFLAG=#  on macintosh, do not static compile (get missing lcrt0.o problem)
	BERTINILIB:=$(HOME)/lib/libbertini.a
	INSTALLPATH:=$(HOME)/bin
	BERTINIHEADERS:=$(HOME)/bertini_source/include
	
	BOOSTLOC := /usr/local/Cellar/boost/1.54.0/lib/
	BOOSTSUFFIX:=-mt
	MACHINESPECIFICS:= 
endif

ifeq ($(MACHINENAME),sam)
    #sam
    INC:=-I/usr/local/include/   $(INC) 
    BISON=bison 
    LIBS=-lm -L/usr/local/lib -lmpfr  -lgmp #
    COMPcc=gcc 
	COMPpp = g++  
    STATICFLAG= -static
	BERTINILIB:=$(HOME)/lib/libbertini.a
	INSTALLPATH:=$(HOME)/bin
	BERTINIHEADERS:= $(HOME)/code/darpa/bertini_source/include
endif

ifeq ($(MACHINENAME),wenrui)
    INC:= $(INC) -I/opt/mpfr/3.0.0/include -I/opt/gmp/5.0.1/include
    LIBS= -lm -L/usr/lib -L/opt/mpfr/3.0.0/lib -L/opt/gmp/5.0.1/lib /opt/mpfr/3.0.0/lib/libmpfr.a /opt/gmp/5.0.1/lib/libgmp.a
    COMPcc=gcc
    COMPpp=g++
    STATICFLAG=-static
    BERTINILIB:=$(HOME)/lib/libbertini.a
    INSTALLPATH:=$(HOME)/bin
    BERTINIHEADERS:=/home/whao/bertini1.4/bertini/include
endif

ifeq ($(MACHINE),yourmachinenamehere)
    #use something like this to configure for your machine
    INC:=$(INC) #put your -I inclusions 
    BISON=bison# #path to bison
    LIBS=# #your libs
	COMPcc = gcc  #compiler
	COMPpp = g++
    STATICFLAG=-static#  need static?
endif

  #note: you can export an env variable of MACHINE to the string of your choice, and do configuration like above.  

###
# end machine specific configuration
###



FLEX=flex

LIBSboost := $(BOOSTLOC)libboost_system$(BOOSTSUFFIX).a $(BOOSTLOC)libboost_filesystem$(BOOSTSUFFIX).a   $(BOOSTLOC)libboost_regex$(BOOSTSUFFIX).a 

LIBS:=$(LIBS) ./libpartitionparse.a $(LIBSboost)

#a few last definitions
INC:=$(INC) -I../include  -I$(BERTINIHEADERS) #this includes the header files for bertini_real, etc   

ARGS= $(OPT) -D_HAVE_MPI $(MACHINESPECIFICS) $(INC)
OBJFILES= parallelism.o derivative_systems.o data_type.o programConfiguration.o determinant_derivative.o output.o fileops.o isosingular.o checkSelfConjugate.o nullspace_left.o surface.o curve.o
SOLVERS= solver.o solver_multilin_projection.o solver_midpoint_tracker.o solver_nullspace_left.o solver_detjactodetjac.o solver_linprodtodetjac.o solver_multilintolin.o postProcessing.o solver_lintolin.o
FLEXFILES=lex.partitionParse.o 



###
#BEGIN MAKE RULES
###

#all comes first, so is the default rule when you `make'
all : libpartitionparse bertini_real data_to_set sampler

libpartitionparse: lex.partitionParse.o
	ar rcs libpartitionparse.a lex.partitionParse.o

#makes the bertini_real program
bertini_real : bertini_real.cpp $(OBJFILES)  $(SOLVERS) $(FLEXFILES)
	$(COMPpp) $(ARGS) -o bertini_real bertini_real.cpp $(OBJFILES)  $(SOLVERS)  $(FLEXFILES) $(BERTINILIB)  $(LIBS)  $(STATICFLAG)


#makes the data2set program.
data_to_set: dataToSet.cpp
	$(COMPpp) $(ARGS) -o data2set dataToSet.cpp $(BERTINILIB) $(LIBS)  $(STATICFLAG)
#   

#makes the sampler program.
sampler: sampler.cpp $(OBJFILES) $(FLEXFILES) $(SOLVERS)
	$(COMPpp) $(ARGS) -o sampler sampler.cpp  $(SOLVERS) $(OBJFILES)   $(BERTINILIB) $(FLEXFILES)  $(LIBS)  $(STATICFLAG)
#   


#simply copy the executables to the installpath
install: bertini_real data_to_set sampler
	cp bertini_real data2set sampler $(INSTALLPATH)




%.c.o : %.h
	$(COMPcc) $(ARGS) -c $*.c

%.o : %.cpp ../include/%.hpp
	$(COMPpp) $(ARGS) -c $*.cpp


lex.partitionParse.o:	lex.partitionParse.c
	$(COMPcc) $(ARGS) -c lex.partitionParse.c

lex.partitionParse.c:	partitionParse.l
	$(FLEX) -PpartitionParse partitionParse.l

clean : 
	rm -f $(OBJFILES) $(FLEXFILES) $(SOLVERS) sampler.o bertini.o data2set.o bertini_real data2set sampler
	rm -f par.out paramDerivs.out arr.out eval.out eval2.out const.out finalFile.out jacV.out jacP.out num.out
	rm -f tempBertiniFile lex.partitionParse.c 
	rm -rf *.dSYM

