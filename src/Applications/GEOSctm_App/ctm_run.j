#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################

#SBATCH --time=@RUN_T
#SBATCH --job-name=@RUN_N
#@RUN_P
#@RUN_Q
#@BATCH_GROUP
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#@PBS -o ctm_run.o@RSTDATE

#######################################################################
#                         System Settings 
#######################################################################

umask 022

limit stacksize unlimited

@SETENVS

#######################################################################
# Configuration Settings
#######################################################################

setenv doIdealizedPT @doIdealizedPT
setenv doGEOSCHEMCHEM @doGEOSCHEMCHEM

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv ARCH `uname`

setenv SITE             @SITE
setenv GEOSDIR          @GEOSDIR 
setenv GEOSBIN          @GEOSBIN 
setenv GEOSUTIL         @GEOSSRC/GMAO_Shared/GEOS_Util
setenv RUN_CMD         "@RUN_CMD"
setenv CTMVER           @CTMVER

source $GEOSBIN/g5_modules
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/${ARCH}/lib

#######################################################################
#             Experiment Specific Environment Variables
#######################################################################


setenv  EXPID   @EXPID
setenv  EXPDIR  @EXPDIR
setenv  HOMDIR  @HOMDIR

setenv  RSTDATE @RSTDATE
setenv  GCMEMIP @GCMEMIP

#######################################################################
#                 Create Experiment Sub-Directories
#######################################################################

if (! -e $EXPDIR/restarts   ) mkdir -p $EXPDIR/restarts
if (! -e $EXPDIR/holding    ) mkdir -p $EXPDIR/holding
if (! -e $EXPDIR/archive    ) mkdir -p $EXPDIR/archive
if (! -e $EXPDIR/post       ) mkdir -p $EXPDIR/post
if (! -e $EXPDIR/plot       ) mkdir -p $EXPDIR/plot

if( $GCMEMIP == TRUE ) then
    if (! -e $EXPDIR/restarts/$RSTDATE ) mkdir -p $EXPDIR/restarts/$RSTDATE
    setenv  SCRDIR  $EXPDIR/scratch.$RSTDATE
else
    setenv  SCRDIR  $EXPDIR/scratch
endif

if (! -e $SCRDIR ) mkdir -p $SCRDIR

#######################################################################
#                   Set Experiment Run Parameters
#######################################################################

set         NX = `grep    "^ *NX:" $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set         NY = `grep    "^ *NY:" $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set GEOSCTM_IM = `grep GEOSctm_IM: $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set GEOSCTM_JM = `grep GEOSctm_JM: $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set GEOSCTM_LM = `grep GEOSctm_LM: $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set    OGCM_IM = `grep    OGCM_IM: $HOMDIR/GEOSCTM.rc | cut -d':' -f2`
set    OGCM_JM = `grep    OGCM_JM: $HOMDIR/GEOSCTM.rc | cut -d':' -f2`

# Check for Over-Specification of CPU Resources
# ---------------------------------------------
  if ($?PBS_NODEFILE) then
     set  NCPUS = `cat $PBS_NODEFILE | wc -l`
     @    NPES  = $NX * $NY
        if( $NPES > $NCPUS ) then
             echo "CPU Resources are Over-Specified"
             echo "--------------------------------"
             echo "Allotted NCPUs: $NCPUS"
             echo "Specified  NX : $NX"
             echo "Specified  NY : $NY"
             exit
        endif
     endif
  endif

#######################################################################
#                       GCMEMIP Setup
#######################################################################

if( $GCMEMIP == TRUE & ! -e $EXPDIR/restarts/$RSTDATE/cap_restart ) then

cd $EXPDIR/restarts/$RSTDATE

set      FILE = strip
/bin/rm $FILE
cat << EOF > $FILE
#!/bin/ksh
/bin/mv \$1 \$1.tmp
touch   \$1
while read line
do
echo \$line >> \$1
done < \$1.tmp
exit
EOF
chmod +x $FILE

@CPEXEC $HOMDIR/CAP.rc .
./strip         CAP.rc

set year  = `echo $RSTDATE | cut -d_ -f1 | cut -b1-4`
set month = `echo $RSTDATE | cut -d_ -f1 | cut -b5-6`

# Copy MERRA-2 Restarts
# ---------------------
@CPEXEC /discover/nobackup/projects/gmao/g6dev/ltakacs/MERRA2/restarts/AMIP/M${month}/restarts.${year}${month}.tar .
@TAREXEC xf  restarts.${year}${month}.tar
/bin/rm restarts.${year}${month}.tar
/bin/rm MERRA2*bin


# Regrid MERRA-2 Restarts
# -----------------------
set RSTID = `/bin/ls *catch*bin | cut -d. -f1`
$GEOSBIN/regrid.pl -np -ymd ${year}${month}10 -hr 21 -grout C${AGCM_IM} -levsout ${AGCM_LM} -outdir . -d . -expid $RSTID -tagin Ganymed-4_0 -oceanin e -i -nobkg -lbl -nolcv -tagout Icarus -rs 3 -oceanout @OCEANOUT
/bin/rm $RSTID.*.bin

     set IMC = $AGCM_IM
if(     $IMC < 10 ) then
     set IMC = 000$IMC
else if($IMC < 100) then
     set IMC = 00$IMC
else if($IMC < 1000) then
     set IMC = 0$IMC
endif

$GEOSBIN/stripname C${AGCM_IM}@OCEANOUT_${RSTID}.
$GEOSBIN/stripname .${year}${month}10_21z.bin.@BCSTAG.@ATMOStag_@OCEANtag
/bin/mv gocart_internal_rst gocart_internal_rst.merra2
$GEOSBIN/gogo.x -s $RSTID.Chem_Registry.rc.${year}${month}10_21z -t $EXPDIR/RC/Chem_Registry.rc -i gocart_internal_rst.merra2 -o gocart_internal_rst -r C${AGCM_IM} -l ${AGCM_LM}


# Create CAP.rc and cap_restart
# -----------------------------
set   nymd = ${year}${month}10
set   nhms = 210000
echo $nymd $nhms > cap_restart

set curmonth = $month
      @ count = 0
while( $count < 4 )
       set date  = `$GEOSBIN/tick $nymd $nhms 86400`
       set nymd  =  $date[1]
       set nhms  =  $date[2]
       set year  = `echo $nymd | cut -c1-4`
       set month = `echo $nymd | cut -c5-6`
       if( $curmonth != $month ) then
        set curmonth  = $month
             @ count  = $count + 1
       endif
end
set oldstring =  `cat CAP.rc | grep END_DATE:`
set newstring =  "END_DATE: ${year}${month}01 210000"
/bin/mv CAP.rc CAP.tmp
cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
/bin/rm CAP.tmp

endif

#######################################################################
#   Move to Scratch Directory and Copy RC Files from Home Directory
#######################################################################

# Driving datasets
setenv DRIVING_DATASETS @DRIVING_DATASETS

cd $SCRDIR
/bin/rm -rf *
                             /bin/ln -sf $EXPDIR/RC/* .
                             @CPEXEC     $EXPDIR/cap_restart .
if ( ${DRIVING_DATASETS} == F515_516 || ${DRIVING_DATASETS} == F5131 ) then
                             @CPEXEC     $EXPDIR/FP_ExtData.rc.tmpl .
else
                             @CPEXEC     $EXPDIR/${DRIVING_DATASETS}_ExtData.rc.tmpl .
endif
                             @CPEXEC -f  $HOMDIR/*.rc .
                             @CPEXEC -f  $HOMDIR/*.nml .
                             @CPEXEC     $GEOSBIN/bundleParser.py .

                             cat fvcore_layout.rc >> input.nml

if( $GCMEMIP == TRUE ) then
    @CPEXEC -f  $EXPDIR/restarts/$RSTDATE/cap_restart .
    @CPEXEC -f  $EXPDIR/restarts/$RSTDATE/CAP.rc .
endif

set END_DATE  = `grep     END_DATE:  CAP.rc | cut -d':' -f2`
set NUM_SGMT  = `grep     NUM_SGMT:  CAP.rc | cut -d':' -f2`
set FSEGMENT  = `grep FCST_SEGMENT:  CAP.rc | cut -d':' -f2`
set USE_SHMEM = `grep    USE_SHMEM:  CAP.rc | cut -d':' -f2`

#######################################################################
#                  Rename GEOS-Chem RC files
#######################################################################

if( ${doGEOSCHEMCHEM} == YES) then

set GEOSCHEM_RCS = ( brc.dat chemga.dat dust.dat FJX_j2j.dat FJX_spec.dat h2so4.dat jv_spec_mie.dat org.dat so4.dat soot.dat ssa.dat ssc.dat input.geos )
foreach FILE ( ${GEOSCHEM_RCS} )
  /bin/mv ${FILE}.rc ${FILE}
end

endif

#######################################################################
#         Create Strip Utility to Remove Multiple Blank Spaces
#######################################################################

set      FILE = strip
/bin/rm $FILE
cat << EOF > $FILE
#!/bin/ksh
/bin/mv \$1 \$1.tmp
touch   \$1
while read line
do
echo \$line >> \$1
done < \$1.tmp
exit
EOF
chmod +x $FILE

#######################################################################
#              Create HISTORY Collection Directories
#######################################################################

set collections = ''
foreach line ("`cat HISTORY.rc`")
   set firstword  = `echo $line | awk '{print $1}'`
   set firstchar  = `echo $firstword | cut -c1`
   set secondword = `echo $line | awk '{print $2}'`

   if ( $firstword == "::" ) goto done

   if ( $firstchar != "#" ) then
      set collection  = `echo $firstword | sed -e "s/'//g"`
      set collections = `echo $collections $collection`
      if ( $secondword == :: ) goto done
   endif

   if ( $firstword == COLLECTIONS: ) then
      set collections = `echo $secondword | sed -e "s/'//g"`
   endif
end

done:
   foreach collection ( $collections )
      if (! -e $EXPDIR/$collection )         mkdir $EXPDIR/$collection
      if (! -e $EXPDIR/holding/$collection ) mkdir $EXPDIR/holding/$collection
   end

#######################################################################
#                        Link Boundary Datasets
#######################################################################

setenv BCSDIR    @BCSDIR
setenv CHMDIR    @CHMDIR
setenv DATELINE  DC
setenv EMISSIONS @EMISSIONS

set             FILE = linkbcs
/bin/rm -f     $FILE
cat << _EOF_ > $FILE
#!/bin/csh -f

/bin/mkdir -p            ExtData
/bin/ln    -sf $CHMDIR/* ExtData

/bin/ln -sf $BCSDIR/Shared/pchem.species.CMIP-5.1870-2097.z_91x72.nc4 species.data
/bin/ln -sf $BCSDIR/Shared/*bin .
/bin/ln -sf $BCSDIR/Shared/*c2l*.nc4 .

#link to Ben's special tile files:
/bin/ln -sf /discover/nobackup/bmauer/tile_files/*/* .

_EOF_


chmod +x linkbcs
@CPEXEC  linkbcs $EXPDIR

#######################################################################
#          Get C2L History weights/index file for Cubed-Sphere
#######################################################################

set C_NPX = `echo $GEOSCTM_IM | awk '{printf "%5.5i", $1}'`
set C_NPY = `echo $GEOSCTM_JM | awk '{printf "%5.5i", $1}'`
set H_NPX = `echo @HIST_IM | awk '{printf "%5.5i", $1}'`
set H_NPY = `echo @HIST_JM | awk '{printf "%5.5i", $1}'`

set c2l_file = "${C_NPX}x${C_NPY}_c2l_${H_NPX}x${H_NPY}.bin"

#######################################################################
#                    Get Executable and RESTARTS 
#######################################################################

@CPEXEC $EXPDIR/GEOSctm.x .

set rst_files      = `cat GEOSCTM.rc | grep "RESTART_FILE"    | grep -v VEGDYN | grep -v "#" | cut -d ":" -f1 | cut -d "_" -f1-2`
set rst_file_names = `cat GEOSCTM.rc | grep "RESTART_FILE"    | grep -v VEGDYN | grep -v "#" | cut -d ":" -f2`

set chk_files      = `cat GEOSCTM.rc | grep "CHECKPOINT_FILE" | grep -v "#" | cut -d ":" -f1 | cut -d "_" -f1-2`
set chk_file_names = `cat GEOSCTM.rc | grep "CHECKPOINT_FILE" | grep -v "#" | cut -d ":" -f2`

# Remove possible bootstrap parameters (+/-)
# ------------------------------------------
set dummy = `echo $rst_file_names`
set rst_file_names = ''
foreach rst ( $dummy )
  set length  = `echo $rst | awk '{print length($0)}'`
  set    bit  = `echo $rst | cut -c1`
  if(  "$bit" == "+" | \
       "$bit" == "-" ) set rst = `echo $rst | cut -c2-$length`
  set rst_file_names = `echo $rst_file_names $rst`
end

# Copy Restarts to Scratch Directory
# ----------------------------------
if( $GCMEMIP == TRUE ) then
    foreach rst ( $rst_file_names )
      if(-e $EXPDIR/restarts/$RSTDATE/$rst ) @CPEXEC $EXPDIR/restarts/$RSTDATE/$rst . &
    end
else
    foreach rst ( $rst_file_names )
      if(-e $EXPDIR/$rst ) @CPEXEC $EXPDIR/$rst . &
    end
endif
wait

# Copy and Tar Initial Restarts to Restarts Directory
# ---------------------------------------------------
set edate = e`cat cap_restart | cut -c1-8`_`cat cap_restart | cut -c10-11`z
set numrs = `/bin/ls -1 ${EXPDIR}/restarts/*${edate}* | wc -l`
if($numrs == 0) then
   foreach rst ( $rst_file_names )
      if( -e $rst & ! -e ${EXPDIR}/restarts/$EXPID.${rst}.${edate}.${CTMVER} ) then
            @CPEXEC $rst ${EXPDIR}/restarts/$EXPID.${rst}.${edate}.${CTMVER} &
      endif
   end
   wait
   cd $EXPDIR/restarts
      @TAREXEC cf  restarts.${edate}.tar $EXPID.*.${edate}.${CTMVER}
     /bin/rm -rf `/bin/ls -d -1     $EXPID.*.${edate}.${CTMVER}`
   cd $SCRDIR
endif

# If any restart is binary, set NUM_READERS to 1 so that
# +-style bootstrapping of missing files can occur in 
# MAPL. pbinary cannot do this, but pnc4 can.
# ------------------------------------------------------
set found_binary = 0

foreach rst ( $rst_file_names )
   if (-e $rst) then
      set rst_type = `/usr/bin/file -Lb --mime-type $rst`
      if ( $rst_type =~ "application/octet-stream" ) then
         set found_binary = 1
      endif
   endif
end

if ($found_binary == 1) then
   /bin/mv GEOSCTM.rc GEOSCTM.tmp
   cat GEOSCTM.tmp | sed -e "/^NUM_READERS/ s/\([0-9]\+\)/1/g" > GEOSCTM.rc
   /bin/rm GEOSCTM.tmp
endif


##################################################################
######
######         Perform multiple iterations of Model Run
######
##################################################################

@ counter    = 1
while ( $counter <= ${NUM_SGMT} )

/bin/rm -f  EGRESS

if( $GCMEMIP == TRUE ) then
    @CPEXEC -f  $EXPDIR/restarts/$RSTDATE/CAP.rc .
else
    @CPEXEC -f $HOMDIR/CAP.rc .
endif

./strip CAP.rc

# Set Time Variables for Current_(c), Ending_(e), and Segment_(s) dates 
# ---------------------------------------------------------------------
set nymdc = `cat cap_restart | cut -c1-8`
set nhmsc = `cat cap_restart | cut -c10-15`
set nymde = `cat CAP.rc | grep END_DATE:     | cut -d: -f2 | cut -c2-9`
set nhmse = `cat CAP.rc | grep END_DATE:     | cut -d: -f2 | cut -c11-16`
set nymds = `cat CAP.rc | grep JOB_SGMT:     | cut -d: -f2 | cut -c2-9`
set nhmss = `cat CAP.rc | grep JOB_SGMT:     | cut -d: -f2 | cut -c11-16`

# Compute Time Variables at the Finish_(f) of current segment
# -----------------------------------------------------------
set nyear   = `echo $nymds | cut -c1-4`
set nmonth  = `echo $nymds | cut -c5-6`
set nday    = `echo $nymds | cut -c7-8`
set nhour   = `echo $nhmss | cut -c1-2`
set nminute = `echo $nhmss | cut -c3-4`
set nsec    = `echo $nhmss | cut -c5-6`
       @ dt = $nsec + 60 * $nminute + 3600 * $nhour + 86400 * $nday

set nymdf = $nymdc
set nhmsf = $nhmsc
set date  = `$GEOSBIN/tick $nymdf $nhmsf $dt`
set nymdf =  $date[1]
set nhmsf =  $date[2]
set year  = `echo $nymdf | cut -c1-4`
set month = `echo $nymdf | cut -c5-6`
set day   = `echo $nymdf | cut -c7-8`

     @  month = $month + $nmonth
while( $month > 12 )
     @  month = $month - 12
     @  year  = $year  + 1
end
     @  year  = $year  + $nyear
     @ nymdf  = $year * 10000 + $month * 100 + $day

if( $nymdf >  $nymde )    set nymdf = $nymde
if( $nymdf == $nymde )    then
    if( $nhmsf > $nhmse ) set nhmsf = $nhmse
endif

set yearc = `echo $nymdc | cut -c1-4`
set yearf = `echo $nymdf | cut -c1-4`

# Select proper MERRA-2 GOCART Emission RC Files
# (NOTE: MERRA2-DD has same transition date)
# ----------------------------------------------
if( ${EMISSIONS} == MERRA2 | \
    ${EMISSIONS} == MERRA2-DD ) then
    set MERRA2_Transition_Date = 20000401

    if( $nymdc < ${MERRA2_Transition_Date} ) then
         set MERRA2_EMISSIONS_DIRECTORY = $GEOSDIR/$ARCH/etc/$EMISSIONS/19600101-20000331
         if( $nymdf > ${MERRA2_Transition_Date} ) then
          set nymdf = ${MERRA2_Transition_Date}
          set oldstring = `cat CAP.rc | grep END_DATE:`
          set newstring = "END_DATE: $nymdf $nhmsf"
          /bin/mv CAP.rc CAP.tmp
                     cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
         endif
    else
         set MERRA2_EMISSIONS_DIRECTORY = $GEOSDIR/$ARCH/etc/$EMISSIONS/20000401-present
    endif

    if( $GEOSCTM_LM == 72 ) then
        @CPEXEC --remove-destination ${MERRA2_EMISSIONS_DIRECTORY}/*.rc .
    else
        set files =      `/bin/ls -1 ${MERRA2_EMISSIONS_DIRECTORY}/*.rc`
        foreach file ($files)
          /bin/rm -f   `basename $file`
          /bin/rm -f    dummy
          @CPEXEC $file dummy
              cat       dummy | sed -e "s|/L72/|/L${GEOSCTM_LM}/|g" | sed -e "s|z72|z${GEOSCTM_LM}|g" > `basename $file`
        end
    endif

endif

# Select the proper ExtData resource file
# when MERRA2 datasets are selected.
#----------------------------------------
if( ${DRIVING_DATASETS} == MERRA2) then
    set startYear = `cat cap_restart | cut -c1-4`
    set oldstring = `cat CAP.rc | grep EXTDATA_CF:`
    set COMPNAME = `grep COMPNAME CAP.rc | awk '{print $2}'`

    if( $startYear > 1979 && $startYear < 1992 ) then
        set sYear  = 1980
        set sMonth = jan79
        set MERRA2type = MERRA2_100
        set data_Transition_Date = 19920101
    else if( $startYear > 1991 && $startYear < 2000 ) then
        set sYear  = 1992
        set sMonth = jan91
        set MERRA2type = MERRA2_200
        set data_Transition_Date = 20000101
    else if( $startYear > 1999 && $startYear < 2010 ) then
        set sYear  = 2000
        set sMonth = jan00
        set MERRA2type = MERRA2_300
        set data_Transition_Date = 20100101
    else if( $startYear > 2009 ) then
        set sYear  = 2010
        set sMonth = jan10
        set MERRA2type = MERRA2_400
        set data_Transition_Date = 20200101
    endif

    set newstring = "EXTDATA_CF: ${COMPNAME}_ExtData_${sYear}.rc"

    # Edit the CAP.rc file
    #---------------------
    /bin/mv CAP.rc CAP.tmp
    cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
    rm -f CAP.tmp

    # Transition into a new decade
    #-----------------------------
    if( $nymdf > ${data_Transition_Date} ) then
        set nymdf = ${data_Transition_Date}
        set oldstring = `cat CAP.rc | grep END_DATE:`
        set newstring = "END_DATE: $nymdf $nhmsf"
        /bin/mv CAP.rc CAP.tmp
        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
        rm -f CAP.tmp
    endif

    # Edit the ExtData.rc file
    #-------------------------
    set tFILE = tmpfile
    rm -f $tFILE
    cat MERRA2_ExtData.rc.tmpl > $tFILE
    set sFILE = sedfile
    rm -f $sFILE
cat << EOF > $sFILE 
s/@sMonth/$sMonth/g
s/@MERRA2type/$MERRA2type/g
EOF

    sed -f $sFILE $tFILE > MAPL_ExtData_${sYear}.rc
    chmod 755  MAPL_ExtData_${sYear}.rc
    rm -f $tFILE
    rm -f $sFILE

# Concatenate required ExtData files
# ----------------------------------
if( ${doIdealizedPT} == YES) then
  @CPEXEC MAPL_ExtData_${sYear}.rc ${COMPNAME}_ExtData_${sYear}.rc
else if (${doGEOSCHEMCHEM} == YES) then
  set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc GOCARTdata_ExtData.rc WSUB_ExtData.rc BC_GridComp_ExtData.rc CO_GridComp_ExtData.rc CO2_GridComp_ExtData.rc DU_GridComp_ExtData.rc NI_GridComp_ExtData.rc OC_GridComp_ExtData.rc SU_GridComp_ExtData.rc GEOSCHEMchem_ExtData.rc`
  cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc

#  sed -i 's/QFED\/NRT/QFED/' ${COMPNAME}_ExtData_${sYear}.rc
#  sed -i 's/v2.5r1_0.1_deg/v2.5r1\/0.1/' ${COMPNAME}_ExtData_${sYear}.rc
else
  set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc *_ExtData.rc`
  cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc
endif

endif

# For MERRA1
#-----------
if( ${DRIVING_DATASETS} == MERRA1) then
    set startYear = `cat cap_restart | cut -c1-4`
    set oldstring = `cat CAP.rc | grep EXTDATA_CF:`
    set COMPNAME = `grep COMPNAME CAP.rc | awk '{print $2}'`

    if( $startYear >= 1979 && $startYear <= 1992 ) then
        set sYear  = 1979
        set MERRA1type = MERRA100
        set data_Transition_Date = 19930101
    else if( $startYear >= 1993 && $startYear <= 2000 ) then
        set sYear  = 1993
        set MERRA1type = MERRA200
        set data_Transition_Date = 20010101
    else if( $startYear >= 2001 ) then
        set sYear  = 2001
        set MERRA1type = MERRA300
        set data_Transition_Date = 20200101
    endif

    set newstring = "EXTDATA_CF: ${COMPNAME}_ExtData_${sYear}.rc"

    # Edit the CAP.rc file
    #---------------------
    /bin/mv CAP.rc CAP.tmp
    cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
    rm -f CAP.tmp

    # Transition into a new date range
    #----------------------------------
    if( $nymdf > ${data_Transition_Date} ) then
        set nymdf = ${data_Transition_Date}
        set oldstring = `cat CAP.rc | grep END_DATE:`
        set newstring = "END_DATE: $nymdf $nhmsf"
        /bin/mv CAP.rc CAP.tmp
        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
        rm -f CAP.tmp
    endif

    # Edit the ExtData.rc file
    #-------------------------
    set tFILE = tmpfile
    rm -f $tFILE
    cat MERRA1_ExtData.rc.tmpl > $tFILE
    set sFILE = sedfile
    rm -f $sFILE
cat << EOF > $sFILE 
s/@MERRA1type/$MERRA1type/g
EOF

    sed -f $sFILE $tFILE > MAPL_ExtData_${sYear}.rc
    chmod 755  MAPL_ExtData_${sYear}.rc
    rm -f $tFILE
    rm -f $sFILE

# Concatenate required ExtData files
# ----------------------------------
set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc *_ExtData.rc`
cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc

endif

# For FPIT
#-----------
if( ${DRIVING_DATASETS} == FPIT) then
    set startYear = `cat cap_restart | cut -c1-4`
    set oldstring = `cat CAP.rc | grep EXTDATA_CF:`
    set COMPNAME = `grep COMPNAME CAP.rc | awk '{print $2}'`

    if( $startYear >= 1979 && $startYear <= 1992 ) then
        set sYear  = 1979
        set FPITtype = UNKNOWN
        set data_Transition_Date = 19930101
    else if( $startYear >= 1993 && $startYear <= 2000 ) then
        set sYear  = 1993
        set FPITtype = UNKNOWN
        set data_Transition_Date = 20010101
    else if( $startYear >= 2001 ) then
        set sYear  = 2001
        set FPITtype = d5124_fpit
        set data_Transition_Date = 20200101
    endif

    set newstring = "EXTDATA_CF: ${COMPNAME}_ExtData_${sYear}.rc"

    # Edit the CAP.rc file
    #---------------------
    /bin/mv CAP.rc CAP.tmp
    cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
    rm -f CAP.tmp

    # Transition into a new date range
    #----------------------------------
    if( $nymdf > ${data_Transition_Date} ) then
        set nymdf = ${data_Transition_Date}
        set oldstring = `cat CAP.rc | grep END_DATE:`
        set newstring = "END_DATE: $nymdf $nhmsf"
        /bin/mv CAP.rc CAP.tmp
        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
        rm -f CAP.tmp
    endif

    # Edit the ExtData.rc file
    #-------------------------
    set tFILE = tmpfile
    rm -f $tFILE
    cat FPIT_ExtData.rc.tmpl > $tFILE
    set sFILE = sedfile
    rm -f $sFILE
cat << EOF > $sFILE 
s/@FPITtype/$FPITtype/g
EOF

    sed -f $sFILE $tFILE > MAPL_ExtData_${sYear}.rc
    chmod 755  MAPL_ExtData_${sYear}.rc
    rm -f $tFILE
    rm -f $sFILE

# Concatenate required ExtData files
# ----------------------------------
set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc *_ExtData.rc`
cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc

endif

#-------------
# For F515_516
#-------------
if( ${DRIVING_DATASETS} == F515_516) then
    set startYear = `cat cap_restart | cut -c1-4`
    set oldstring = `cat CAP.rc | grep EXTDATA_CF:`
    set COMPNAME = `grep COMPNAME CAP.rc | awk '{print $2}'`

    if( $startYear >= 1979 && $startYear <= 1992 ) then
        set sYear  = 1979
        set FPtype = UNKNOWN
        set FPver  = UNKNOWN
        set FPmod  = UNKNOWN
        set data_Transition_Date = 19930101
    else if( $startYear >= 1993 && $startYear <= 2016 ) then
        set sYear  = 1993
        set FPtype = f515_fpp
        set FPver  = 5_15
        set FPmod  = 5.15
        set data_Transition_Date = 20170101
    else if( $startYear >= 2017 ) then
        set sYear  = 2017
        set FPtype = f516_fp
        set FPver  = 5_16
        set FPmod  = 5.16
        set data_Transition_Date = 20200101
    endif

    set newstring = "EXTDATA_CF: ${COMPNAME}_ExtData_${sYear}.rc"

    # Edit the CAP.rc file
    #---------------------
    /bin/mv CAP.rc CAP.tmp
    cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
    rm -f CAP.tmp

    # Transition into a new date range
    #----------------------------------
    if( $nymdf > ${data_Transition_Date} ) then
        set nymdf = ${data_Transition_Date}
        set oldstring = `cat CAP.rc | grep END_DATE:`
        set newstring = "END_DATE: $nymdf $nhmsf"
        /bin/mv CAP.rc CAP.tmp
        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
        rm -f CAP.tmp
    endif

    # Edit the ExtData.rc file
    #-------------------------
    set tFILE = tmpfile
    rm -f $tFILE
    cat FP_ExtData.rc.tmpl > $tFILE
    set sFILE = sedfile
    rm -f $sFILE
cat << EOF > $sFILE 
s/@FPtype/$FPtype/g
s/@FPver/$FPver/g
s/@FPmod/$FPmod/g
EOF

    sed -f $sFILE $tFILE > MAPL_ExtData_${sYear}.rc
    chmod 755  MAPL_ExtData_${sYear}.rc
    rm -f $tFILE
    rm -f $sFILE

# Rename ExtData files that are not needed
# ----------------------------------------
set            SC_TRUE = `grep -i "^ *ENABLE_STRATCHEM *: *\.TRUE\."     GEOS_ChemGridComp.rc | wc -l`
if (          $SC_TRUE == 0 && -e StratChem_ExtData.rc          ) /bin/mv          StratChem_ExtData.rc          StratChem_ExtData.rc.NOT_USED
set           GMI_TRUE = `grep -i "^ *ENABLE_GMICHEM *: *\.TRUE\."       GEOS_ChemGridComp.rc | wc -l`
if (         $GMI_TRUE == 0 && -e GMI_ExtData.rc                ) /bin/mv                GMI_ExtData.rc                GMI_ExtData.rc.NOT_USED
set           GCC_TRUE = `grep -i "^ *ENABLE_GEOSCHEM *: *\.TRUE\."      GEOS_ChemGridComp.rc | wc -l`
if (         $GCC_TRUE == 0 && -e GEOSCHEMchem_ExtData.rc       ) /bin/mv       GEOSCHEMchem_ExtData.rc       GEOSCHEMchem_ExtData.rc.NOT_USED
set         CARMA_TRUE = `grep -i "^ *ENABLE_CARMA *: *\.TRUE\."         GEOS_ChemGridComp.rc | wc -l`
if (       $CARMA_TRUE == 0 && -e CARMAchem_GridComp_ExtData.rc ) /bin/mv CARMAchem_GridComp_ExtData.rc CARMAchem_GridComp_ExtData.rc.NOT_USED
set           DNA_TRUE = `grep -i "^ *ENABLE_DNA *: *\.TRUE\."           GEOS_ChemGridComp.rc | wc -l`
if (         $DNA_TRUE == 0 && -e DNA_ExtData.rc                ) /bin/mv                DNA_ExtData.rc                DNA_ExtData.rc.NOT_USED
set         ACHEM_TRUE = `grep -i "^ *ENABLE_ACHEM *: *\.TRUE\."         GEOS_ChemGridComp.rc | wc -l`
if (       $ACHEM_TRUE == 0 && -e GEOSachem_ExtData.rc          ) /bin/mv          GEOSachem_ExtData.rc          GEOSachem_ExtData.rc.NOT_USED
set   GOCART_DATA_TRUE = `grep -i "^ *ENABLE_GOCART_DATA *: *\.TRUE\."   GEOS_ChemGridComp.rc | wc -l`
if ( $GOCART_DATA_TRUE == 0 && -e GOCARTdata_ExtData.rc         ) /bin/mv         GOCARTdata_ExtData.rc         GOCARTdata_ExtData.rc.NOT_USED

# Alternate syntax:
#set EXT_FILE=StratChem_ExtData.rc
#if(`grep -i "^ *ENABLE_STRATCHEM *: *\.TRUE\."     GEOS_ChemGridComp.rc | wc -l`==0 && -e $EXT_FILE) /bin/mv $EXT_FILE $EXT_FILE.NOT_USED

if( ${doGEOSCHEMCHEM} == YES) then

# Rename ExtData file that conflicts w/ GEOS-Chem
# -----------------------------------------------
/bin/mv  HEMCOgocart_ExtData.rc  HEMCOgocart_ExtData.rc.NOT_USED
EOF

endif


# Concatenate required ExtData files
# ----------------------------------
set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc *_ExtData.rc`
cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc

endif

#-------------
# For F5131
#-------------
if( ${DRIVING_DATASETS} == F5131) then
    set startYear = `cat cap_restart | cut -c1-4`
    set oldstring = `cat CAP.rc | grep EXTDATA_CF:`
    set COMPNAME = `grep COMPNAME CAP.rc | awk '{print $2}'`

    if( $startYear >= 1979 && $startYear <= 1992 ) then
        set sYear  = 1979
        set FPtype = UNKNOWN
        set FPver  = UNKNOWN
        set FPmod  = UNKNOWN
        set data_Transition_Date = 19930101
    else if( $startYear >= 1993 && $startYear <= 2000 ) then
        set sYear  = 1993
        set FPtype = UNKNOWN
        set FPver  = UNKNOWN
        set FPmod  = UNKNOWN
        set data_Transition_Date = 20010101
    else if( $startYear >= 2001 ) then
        set sYear  = 2001
        set FPtype = e5131_fp
        set FPver  = 5_13_1
        set FPmod  = 5.13.1
        set data_Transition_Date = 20200101
    endif

    set newstring = "EXTDATA_CF: ${COMPNAME}_ExtData_${sYear}.rc"

    # Edit the CAP.rc file
    #---------------------
    /bin/mv CAP.rc CAP.tmp
    cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
    rm -f CAP.tmp

    # Transition into a new date range
    #----------------------------------
    if( $nymdf > ${data_Transition_Date} ) then
        set nymdf = ${data_Transition_Date}
        set oldstring = `cat CAP.rc | grep END_DATE:`
        set newstring = "END_DATE: $nymdf $nhmsf"
        /bin/mv CAP.rc CAP.tmp
        cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
        rm -f CAP.tmp
    endif

    # Edit the ExtData.rc file
    #-------------------------
    set tFILE = tmpfile
    rm -f $tFILE
    cat FP_ExtData.rc.tmpl > $tFILE
    set sFILE = sedfile
    rm -f $sFILE
cat << EOF > $sFILE 
s/@FPtype/$FPtype/g
s/@FPver/$FPver/g
s/@FPmod/$FPmod/g
EOF

    sed -f $sFILE $tFILE > MAPL_ExtData_${sYear}.rc
    chmod 755  MAPL_ExtData_${sYear}.rc
    rm -f $tFILE
    rm -f $sFILE

# Concatenate required ExtData files
# ----------------------------------
set EXTDATA_FILES = `/bin/ls -1 MAPL_ExtData_${sYear}.rc *_ExtData.rc`
cat ${EXTDATA_FILES} > ${COMPNAME}_ExtData_${sYear}.rc

endif

# don't read NRT QFED data
sed -i 's/QFED\/NRT/QFED/'             ${COMPNAME}_ExtData_${sYear}.rc
sed -i 's/v2.5r1_0.1_deg/v2.5r1\/0.1/' ${COMPNAME}_ExtData_${sYear}.rc

# Run bundleParser.py
#---------------------
python bundleParser.py

# Link Boundary Conditions for Appropriate Date
# ---------------------------------------------
setenv YEAR $yearc
./linkbcs

# Run GEOSctm.x
# -------------
if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh
       @  NPES = $NX * $NY
$RUN_CMD $NPES ./GEOSctm.x
if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh

set rc =  $status

echo GEOSctm Run Status: $rc

if ( $rc != 0 ) then
   echo 'CTM error: exit'
   exit
endif
 
#######################################################################
#   Rename Final Checkpoints => Restarts for Next Segment and Archive
#        Note: cap_restart contains the current NYMD and NHMS
#######################################################################

set edate  = e`cat cap_restart | cut -c1-8`_`cat cap_restart | cut -c10-11`z

set numrst = `echo $rst_files | wc -w`
set numchk = `echo $chk_files | wc -w`

@ n = 1
@ z = $numrst + 1
while ( $n <= $numchk )
   if ( -e $chk_file_names[$n] ) then
       @ m = 1
       while ( $m <= $numrst )
       if(    $chk_files[$n] == $rst_files[$m] || \
            \#$chk_files[$n] == $rst_files[$m]    ) then

            set   chk_type = `/usr/bin/file -Lb --mime-type $chk_file_names[$n]`
            if ( $chk_type =~ "application/octet-stream" ) then
                  set ext  = bin
            else
                  set ext  = nc4
            endif

           /bin/mv $chk_file_names[$n] $rst_file_names[$m]
           @CPEXEC $rst_file_names[$m] ${EXPDIR}/restarts/$EXPID.${rst_file_names[$m]}.${edate}.${CTMVER}.$ext &
           @ m = $numrst + 999
       else
           @ m = $m + 1
       endif
       end
       wait
       if( $m == $z ) then
           echo "Warning!!  Could not find CHECKPOINT/RESTART match for:  " $chk_files[$n]
           exit
       endif
   endif
@ n = $n + 1
end


# TAR ARCHIVED RESTARTS
# ---------------------
cd $EXPDIR/restarts
    if( $FSEGMENT == 00000000 ) then
         @TAREXEC cf  restarts.${edate}.tar $EXPID.*.${edate}.${CTMVER}.*
        /bin/rm -rf `/bin/ls -d -1     $EXPID.*.${edate}.${CTMVER}.*`
    endif
cd $SCRDIR

#######################################################################
#               Move HISTORY Files to Holding Directory
#######################################################################

# Check for files waiting in /holding
# -----------------------------------
set     waiting_files = `/bin/ls -1 $EXPDIR/holding/*/*nc4`
set num_waiting_files = $#waiting_files

# Move current files to /holding
# ------------------------------
foreach collection ( $collections )
   /bin/mv `/bin/ls -1 *.${collection}.*` $EXPDIR/holding/$collection
end

#######################################################################
#                Submit Post-Processing (if necessary)
#######################################################################

if( $num_waiting_files == 0 ) then

cd   $EXPDIR/post
/bin/rm -f sedfile
cat >      sedfile << EOF
s/@POST_O/ctm_post.${edate}/g
s/@COLLECTION/ALL/g
s/-rec_plt @YYYYMM//
EOF
sed -f sedfile ctm_post.j > ctm_post.jtmp
chmod 755  ctm_post.jtmp
 qsub      ctm_post.jtmp
/bin/rm -f ctm_post.jtmp
/bin/rm -f sedfile
cd   $SCRDIR

endif

#######################################################################
#                         Update Iteration Counter
#######################################################################

set enddate = `echo  $END_DATE | cut -c1-8`
set capdate = `cat cap_restart | cut -c1-8`

if ( $capdate < $enddate ) then
@ counter = $counter    + 1
else
@ counter = ${NUM_SGMT} + 1
endif

end

#######################################################################
#                              Re-Submit Job
#######################################################################

if( $GCMEMIP == TRUE ) then
     foreach rst ( `/bin/ls -1 *_rst` )
        /bin/rm -f $EXPDIR/restarts/$RSTDATE/$rst
     end
        /bin/rm -f $EXPDIR/restarts/$RSTDATE/cap_restart
     foreach rst ( `/bin/ls -1 *_rst` )
       @CPEXEC $rst $EXPDIR/restarts/$RSTDATE/$rst &
     end
     wait
     @CPEXEC cap_restart $EXPDIR/restarts/$RSTDATE/cap_restart
else
     foreach rst ( `/bin/ls -1 *_rst` )
        /bin/rm -f $EXPDIR/$rst
     end
        /bin/rm -f $EXPDIR/cap_restart
     foreach rst ( `/bin/ls -1 *_rst` )
       @CPEXEC $rst $EXPDIR/$rst &
     end
     wait
     @CPEXEC cap_restart $EXPDIR/cap_restart
endif


if ( $rc == 0 ) then
      cd  $HOMDIR
      if( $GCMEMIP == TRUE ) then
          if( $capdate < $enddate ) qsub $HOMDIR/ctm_run.j$RSTDATE
      else
          if( $capdate < $enddate ) qsub $HOMDIR/ctm_run.j
      endif
endif
