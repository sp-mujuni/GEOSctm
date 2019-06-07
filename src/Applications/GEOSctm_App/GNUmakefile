#
# Makefile for FVAdvStandalone Tracer Advection Test Application
#
# REVISION HISTORY:
#
# 20060922   Sawyer      Taken from fvsa/src/Application/GNUmakefile
# 20061203   Sawyer      Naming scheme revisions
# 20070608   Sawyer      Upgraded for MAPL implementation
# 20071224   Sawyer      Moved FV_Shared files to separate directory
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS = GEOSctm
BIN  = $(THIS).x

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS = GEOSctm
BIN  = $(THIS).x
TEXS = $(THIS)

PROTEX_FLAGS = -bM


esma_install install: $(BIN)
	$(MKDIR) $(ESMABIN) $(ESMAETC) 
	$(CP) -p *.x    $(ESMABIN)

#	$(CP) -p *.rc   $(ESMAETC)

esma_clean esma_distclean clean distclean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


#                  ---------------------------
#                  Availability Driven Options
#                  ---------------------------

LIB_GEOSCHEM = $(wildcard $(ESMALIB)/libGEOSchem_GridComp.a \
                          $(ESMALIB)/libStratChem_GridComp.a\
                          $(ESMALIB)/libSC_GridComp.a\
                          $(ESMALIB)/libGEOSpchem_GridComp.a\
                          $(ESMALIB)/libGEOSachem_GridComp.a\
                          $(ESMALIB)/libGMIchem_GridComp.a\
                          $(ESMALIB)/libGEOSCHEMchem_GridComp.a\
                          $(ESMALIB)/libGAAS_GridComp.a\
                          $(ESMALIB)/libGOCART_GridComp.a\
                          $(ESMALIB)/libCFC_GridComp.a\
                          $(ESMALIB)/libCARMA_GridComp.a\
                          $(ESMALIB)/libCARMAchem_GridComp.a\
                          $(ESMALIB)/libRn_GridComp.a\
                          $(ESMALIB)/libBC_GridComp.a\
                          $(ESMALIB)/libCO2_GridComp.a\
                          $(ESMALIB)/libCO_GridComp.a\
                          $(ESMALIB)/libDU_GridComp.a\
                          $(ESMALIB)/libO3_GridComp.a\
                          $(ESMALIB)/libOC_GridComp.a\
                          $(ESMALIB)/libSS_GridComp.a\
                          $(ESMALIB)/libSU_GridComp.a\
                          $(ESMALIB)/libCH4_GridComp.a\
                          $(ESMALIB)/libNI_GridComp.a\
                          $(ESMALIB)/libMAMchem_GridComp.a\
                          $(ESMALIB)/libTR_GridComp.a\
                          $(ESMALIB)/libMATRIXchem_GridComp.a\
                          $(ESMALIB)/libDNA_GridComp.a \
                          $(ESMALIB)/libHEMCO_GridComp.a  )


INC_GEOSCHEM = $(wildcard $(ESMAINC)/Chem_Base \
                          $(ESMAINC)/Chem_Shared \
                          $(ESMAINC)/GEOSchem_GridComp \
                          $(ESMAINC)/GMIchem_GridComp ) 


LIB_SHARED = $(wildcard \
             $(LIB_ANASUPP) \
             $(ESMALIB)/libGFDL_fms.a \
             $(ESMALIB)/libGMAO_pilgrim.a \
             $(ESMALIB)/libChem_Shared.a \
             $(ESMALIB)/libChem_Base.a \
             $(ESMALIB)/libGEOS_Shared.a \
             $(ESMALIB)/libMAPL_Base.a \
             $(ESMALIB)/libMAPL_cfio_r4.a \
             $(ESMALIB)/libGMAO_gfio_r4.a \
             $(ESMALIB)/libGMAO_ods.a  \
             $(ESMALIB)/libGMAO_gems.a  \
             $(ESMALIB)/libGMAO_hermes.a  \
             $(ESMALIB)/libGMAO_mpeu.a \
             $(ESMALIB)/libBRC_GridComp.a )

LIB_DYN = $(wildcard $(ESMALIB)/libFVdycoreCubed_GridComp.a \
                     $(ESMALIB)/libfvdycore.a )

LIB_OTHERS = $(wildcard $(ESMALIB)/libGEOSctm_GridComp.a \
                     $(ESMALIB)/libCTMpTracers_GridComp.a\
                     $(ESMALIB)/libCTMdiffusion_GridComp.a \
                     $(ESMALIB)/libCTMconvection_GridComp.a)

LIB_SHARED += $(LIB_DYN)

#                  --------------------
#                  User Defined Targets
#                  --------------------

LIB_COMP = $(LIB_OTHERS) $(LIB_DYN) $(LIB_GEOSCHEM) $(LIB_SHARED)

FREAL = $(FREAL4) 
SRCS = $(THIS).F90
OBJS = $(SRCS:.F90=.o)

DIR_GC = $(wildcard $(ESMAINC)/GEOS[a-z]*GridComp $(INC_GEOSCHEM) )
INC_GC = $(foreach FF,$(DIR_GC),$(M)$(FF))

INC_DIRS = . $(INC_GEOS_GCS) $(INC_GMAO_SHARED) $(INC_ESMF) $(INC_MAPL_BASE)
MOD_DIRS = . $(INC_DIRS) $(ESMAINC)/GEOSctm_GridComp

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))
USER_FFLAGS = $(USER_DEFS) $(BIG_ENDIAN)

vpath % $(MOD_DIRS)

$(BIN) bin : $(OBJS)
	echo "should include stubs $(LIB_ESMF)"
	$(LD) -o $(BIN) $(LDFLAGS) $(OBJS) \
                         $(LIB_COMP) $(LIB_ESMF) $(LIB_SDF) \
                         $(LIB_SCI) $(LIB_MPI) $(LIB_SYS)

#                  --------------------
#                      Dependencies
#                  --------------------

$(THIS).o: $(LIB_COMP)

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
