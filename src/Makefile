
#choose some options.  by default, compiler optimization is turned off
OPT:= -Dverbose  -g  -Ddebug_compile -Dthresholding # -Dfunctionentry_output # -O3 #-Wdeprecated # -Dprintpathlinprod  -O3 -funroll-loops -fexpensive-optimizations

COMPFLAGS= -Wsign-compare -Wformat -Wall # -Wunused-parameter -Wdocumentation -Wdocumentation-unknown-command
OPT := $(OPT) $(COMPFLAGS)



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
    LIBS=  -lm -L/usr/local/lib/ -lmpfr -lgmp          #-lm -lmpfr -lgmp# brake has libraries in path for linking, via brew.  no need to specify here
    COMPcc = mpicc -cc=clang 
	COMPpp = mpic++ -cxx=clang++ -std=c++11
    STATICFLAG=#  on macintosh, do not static compile (get missing lcrt0.o problem)
	BERTINILIB:=$(HOME)/lib/libbertini_parallel.a
	INSTALLPATH:=$(HOME)/bin
	BERTINIHEADERS:=$(HOME)/bertini_source/include

	BOOSTLOC := /usr/local/Cellar/boost/1.54.0/lib/
	BOOSTSUFFIX:=-mt
	MACHINESPECIFICS:= 
endif

ifeq ($(MACHINENAME),sam)
    #sam
    INC:=-I/usr/local/mpfr/include/ -I/usr/local/include/ -I/usr/local/boost/include $(INC) 
    BISON=bison 
    LIBS=-L/usr/local/lib/ -lm -L/usr/local/mpfr/lib/ -lmpfr -L/usr/local/lib/  -lgmp #
    COMPcc=mpicc -cc=gcc 
	COMPpp = mpic++ -cxx=g++  -std=c++11
    STATICFLAG=
	BERTINILIB:=$(HOME)/lib/libbertini_parallel.a
	INSTALLPATH:=$(HOME)/bin
	BERTINIHEADERS:= $(HOME)/code/bertini_source/include

	BOOSTLOC := /usr/local/boost/lib/
	BOOSTSUFFIX :=#no spaces

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

LIBSboost := $(BOOSTLOC)libboost_system$(BOOSTSUFFIX).a $(BOOSTLOC)libboost_filesystem$(BOOSTSUFFIX).a   $(BOOSTLOC)libboost_regex$(BOOSTSUFFIX).a  $(BOOSTLOC)libboost_timer$(BOOSTSUFFIX).a $(BOOSTLOC)libboost_chrono$(BOOSTSUFFIX).a

LIBS:=$(LIBS) ./libpartitionparse.a $(LIBSboost)

#a few last definitions
INC:=$(INC) -I../include  -I$(BERTINIHEADERS) #this includes the header files for bertini_real, etc   

ARGS= $(OPT) -D_HAVE_MPI $(MACHINESPECIFICS) $(INC)

OBJFILES= bertini_extensions.o data_type.o parallelism.o derivative_systems.o programConfiguration.o fileops.o isosingular.o checkSelfConjugate.o curve.o surface.o nullspace_left.o 

SOLVERS= solver.o solver_midpoint_tracker.o solver_nullspace_left.o solver_multilintolin.o postProcessing.o solver_sphere.o

FLEXFILES=lex.partitionParse.o 



###
#BEGIN MAKE RULES
###

#all comes first, so is the default rule when you `make'
all : libpartitionparse bertini_real sampler

libpartitionparse: lex.partitionParse.o
	ar rcs libpartitionparse.a lex.partitionParse.o

#makes the bertini_real program
bertini_real : bertini_real.cpp $(OBJFILES)  $(SOLVERS) $(FLEXFILES)
	$(COMPpp) $(ARGS) -o bertini_real bertini_real.cpp  $(OBJFILES)  $(SOLVERS)  $(FLEXFILES)  -Bstatic  $(BERTINILIB) $(LIBS) -Bdynamic
	#$(COMPpp) $(ARGS) -o bertini_real bertini_real.cpp $(OBJFILES)  $(SOLVERS)  $(FLEXFILES) -Wl,-Bstatic $(BERTINILIB)  $(LIBS) -Wl,-Bdynamic



#makes the sampler program.
sampler: sampler.cpp $(OBJFILES) $(FLEXFILES) $(SOLVERS)
	$(COMPpp) $(ARGS) -o sampler sampler.cpp  -Bstatic  $(SOLVERS) $(OBJFILES)   $(BERTINILIB) $(FLEXFILES)  $(LIBS)  -Bdynamic
	#$(COMPpp) $(ARGS) -o sampler sampler.cpp  -Wl,-Bstatic  $(SOLVERS) $(OBJFILES)   $(BERTINILIB) $(FLEXFILES)  $(LIBS)  -Wl,-Bdynamic 


#simply copy the executables to the installpath
install: all
	cp bertini_real sampler $(INSTALLPATH)




%.c.o : %.h
	$(COMPcc) $(ARGS) -c $*.c

%.o : %.cpp ../include/%.hpp
	$(COMPpp) $(ARGS) -c $*.cpp


lex.partitionParse.o:	lex.partitionParse.c
	$(COMPcc) $(ARGS) -c lex.partitionParse.c

lex.partitionParse.c:	partitionParse.l
	$(FLEX) -PpartitionParse partitionParse.l

clean : 
	rm -f $(OBJFILES) $(FLEXFILES) $(SOLVERS) sampler.o bertini.o bertini_real sampler
	rm -f par.out paramDerivs.out arr.out eval.out eval2.out const.out finalFile.out jacV.out jacP.out num.out
	rm -f tempBertiniFile lex.partitionParse.c 
	rm -rf *.dSYM

