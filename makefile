###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = pr_sdp
CPPSRC	= pr_sdp.cpp\
            BlockMatrix.cpp\
            Matrix.cpp\
            Vector.cpp\
            rxTPM.cpp\
            Tools.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = ./include

LIBS= -llapack -lblas -lgsl

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -O3
FFLAGS   = -O3
LDFLAGS	= -O3


# =============================================================================
#   Targets & Rules
# =============================================================================
all: I1

#------------------------------------------------------------------------------
#  Compile with only P and Q conditions activated
#------------------------------------------------------------------------------

I1:
	@echo
	@echo '  +++ Building $(BINNAME) with I1 condition'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__I1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1 condition successfully!'; \
	   echo; \
	 fi

I:
	@echo
	@echo '  +++ Building $(BINNAME) with I1 and I2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__I"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1 and I2 conditions successfully!'; \
	   echo; \
	 fi

I1Q2:
	@echo
	@echo '  +++ Building $(BINNAME) with I1 and Q2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__I1Q2"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1 and Q2 conditions successfully!'; \
	   echo; \
	 fi

IQ2:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2 and Q2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQ2"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2 and Q2 conditions successfully!'; \
	   echo; \
	 fi

IQ1:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2 and Q1 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQ1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2 and Q1 conditions successfully!'; \
	   echo; \
	 fi

IQ:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2, Q1 and Q2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQ"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2, Q1 and Q2 conditions successfully!'; \
	   echo; \
	 fi

IQG1:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2, Q1, Q2 and G1 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQG1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2, Q1, Q2 and G1 conditions successfully!'; \
	   echo; \
	 fi

IQG2:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2, Q1, Q2 and G2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQG2"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2, Q1, Q2 and G2 conditions successfully!'; \
	   echo; \
	 fi

IQG:
	@echo
	@echo '  +++ Building $(BINNAME) with I1, I2, Q1, Q2, G1 and G2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-D__IQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with I1, I2, Q1, Q2, G1 and G2 conditions successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

#-----------------------------------------------------------------------------
# Make the documentation
#----------------------------------------------------------------------------
doc:
	@doxygen doc-config

# ====================== End of file 'makefile.in' ========================== #
