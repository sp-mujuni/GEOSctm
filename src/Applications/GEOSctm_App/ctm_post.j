#!/bin/csh -x

#######################################################################
#                Batch Parameters for Post-Processing Job
#######################################################################

#@BATCH_TIME@POST_T
#@POST_P
#@BATCH_JOBNAME@POST_N_@COLLECTION.@YYYYMM
#@POST_Q
#@BATCH_GROUP
#@BATCH_OUTPUTNAME@POST_O
#@BATCH_JOINOUTERR

#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

@SETENVS

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv ARCH `uname`

setenv SITE             @SITE
setenv GEOSBIN          @GEOSBIN
setenv GEOSUTIL         @GEOSSRC
setenv BATCHNAME       "@POST_N"

source $GEOSBIN/g5_modules
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/${ARCH}/lib

if( $?SLURM_NTASKS ) then
      setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "
      set NCPUS = $SLURM_NTASKS
else if( $?PBS_NODEFILE ) then
      setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "
      set NCPUS = `cat $PBS_NODEFILE | wc -l`
else
      set NCPUS = NULL
endif

#######################################################################
#                      Perform Post Processing
#######################################################################

setenv POST_SCRIPT @EXPDIR/post/ctmpost.script.$SLURM_JOB_ID

/bin/cp  $GEOSUTIL/post/gcmpost.script $POST_SCRIPT

# A few simple edits change the GCM script to a CTM script:
sed -i -e "s/gcm_run.j/ctm_run.j/g"       $POST_SCRIPT
sed -i -e "s/GEOS5.0/GEOSctm/g"           $POST_SCRIPT
sed -i -e "s/AGCM.rc/GEOSCTM.rc/g"        $POST_SCRIPT
sed -i -e "s/AGCM_IM/GEOSctm_IM/g"        $POST_SCRIPT
sed -i -e "s/AGCM_JM/GEOSctm_JM/g"        $POST_SCRIPT
sed -i -e "s/GEOSgcm.x/GEOSctm.x/g"       $POST_SCRIPT
sed -i -e "s/gcm_archive/ctm_archive/g"   $POST_SCRIPT

$POST_SCRIPT -source @EXPDIR -ncpus $NCPUS -collections @COLLECTION -rec_plt @YYYYMM

exit
