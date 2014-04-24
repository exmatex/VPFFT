#----------------------------------------------------
#   Makefile
#
#   Author:  Frankie Li
#   e-mail:  li31@llnl.gov
#
#   Purpose:  This is a basic Makefile used to compile
#             VPFFT++.   
#----------------------------------------------------

#----------------------------------------------------
#
#  Project (executable) Name
#
#----------------------------------------------------
EXEC   = VPFFT++

#----------------------------------------------------
#
#  Specifying source path and source files for this
#  project.
#
#----------------------------------------------------
PROJ_PATH = ./Src

PROJ_CPPS   =  main.cpp LinearAlgebra.cpp Solvers.cpp MaterialGrid.cpp UnitTester.cpp  Debug.cpp Config.cpp Parser.cpp



#main.cpp LinearAlgebra.cpp Solvers.cpp MaterialGrid.cpp UnitTester.cpp  Debug.cpp

CPPS =  $(patsubst %.cpp, $(PROJ_PATH)/%.cpp, $(PROJ_CPPS) )

#----------------------------------------------------
#  Specify external template library location (BOOST)
#----------------------------------------------------
EIGEN_PATH = ExtLib/eigen-eigen-6e7488e20373/

OBJS =  $(patsubst %.cpp, %.o, $(CPPS))
OBJS_DEBUG =  $(patsubst %.cpp, %.do, $(CPPS))

#----------------------------------------------------
#
#  Specifying precomile external library locations
#
#----------------------------------------------------
INCL = -I. -I${EIGEN_PATH} -I/opt/local/include/

#----------------------------------------------------
#
#  Putting together library locations
#
#----------------------------------------------------
LIB_ROOT = -L/opt/local/lib/ 

#----------------------------------------------------
#
#----------------------------------------------------

RELEASE_FLAGS =   -Wfatal-errors -O3 -DNDEBUG  \
								  -fopenmp 
DEBUG_FLAGS    =  -fopenmp  -ggdb -DDEBUG_LEVEL_MAX \
									-Wall

RELEASE_FLAGS +=  ${LIB_ROOT}
DEBUG_FLAGS   +=  ${LIB_ROOT}	

LIBS   =  -lm -lfftw3_threads -lfftw3 -lfftw3_mpi -lpthread
CC     = mpic++ # /opt/local/bin/g++-mp-4.5   # g++ #icpc



#---------------------------
#  This is the dumb way
#---------------------------
debug:EXEC_DEBUG

EXEC_DEBUG:$(OBJS_DEBUG)
	$(CC) $(DEBUG_FLAGS) $(OBJS_DEBUG) -o $(EXEC)  $(LIBS)	
%.do : %.cpp
	$(CC) -c $(INCL) $(DEBUG_FLAGS) $< -o $@

#---------------------------
# RELEASE
#---------------------------
$(EXEC):$(OBJS)
	$(CC) $(RELEASE_FLAGS) $(OBJS) -o $(EXEC)  $(LIBS)	

%.o : %.cpp
	$(CC) -c  $(INCL) $(RELEASE_FLAGS) $< -o $@


clean:
	rm -f $(OBJS) $(OBJS_DEBUG) $(EXEC) *.exe *~ $(PROJ_PATH)/*~   # -- taken out temporarily

