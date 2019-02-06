C     
C      #######    #######   ##    ##     ######    #######    ########
C     ##         ##         ##    ##    ##    ##   ##    ##   ##
C     ##         ##         ##    ##    ##    ##   ##    ##   ##
C     ##         #######    ########    ##    ##   #######    #######
C     ##               ##   ##    ##    ##    ##   ##  ##     ##
C     ##               ##   ##    ##    ##    ##   ##   ##    ##
C      #######   #######    ##    ##     ######    ##    ##   ########
C     
C     Nobuhisa Kobayashi:         
C     Center for Applied Coastal Research
C     University of Delaware, Newark, Delaware 19716
C     
C     Cross-shore wave transformation with Brad Johnson in 1998
C     
C     Cross-shore sediment transport with Yukiko Tega in 2003
C     and Andres Payo in 2005
C     
C     Bottom permeability with Lizbeth Meigs in 2004
C     
C     Roller effects with Haoyu Zhao in 2004
C     
C     Wave runup, overtopping and transmission with Paco de los Santos in 2006
C     and with Jill Pietropaolo in 2011
C     
C     Longshore current and sediment transport with Arpit Agarwal in 2005
C     
C     Longshore bedload transport rate and wind stresses with Andres Payo
C     in 2007
C     
C     Wave and current interaction and impermeable and permeable wet/dry zone
C     (no sediment and with sediment)with Ali Farhadzadeh in 2008
C     
C     Sediment transport on hard bottom (limited sediment availability),
C     and new input options for storm surge and wave time series as well as
C     for permeable bottom profile evolution with Ali Farhadzadeh in 2009
C     
C     Calibration, improvement and verification of CSHORE
C     by Brad Johnson, Mark Gravens(mg) and Jens Figlus in 2009
C     
C     Infiltration landward of dune crest and dip effect above still water
C     shoreline with Jens Figlus in 2010
C     
C     Onshore ridge migration into ponded runnel, oblique waves on
C     permeable wet/dry zone, and tidal effect on currents with
C     Jens Figlus in 2010
C     
C     Multiple cross-shore lines for alongshore gradient of longshore sediment
C     transport and its effect on beach profile evolution with Hooyoung Jung and
C     Kideok Do in 2011
C     
C     Vegetation effect on wave overtopping and overwash with Kideok Do,
C     Christine Grahler and Berna Ayat (IVEG=1 and 2) in 2012 and 
C     extension to pile fence with Rebecca Quan in 2013
C     
C     Improvement of CSHORE programming by Brad Johnson (bdj) in 2012 
C     
C     Erosion of grass roots and soil on dikes (IPROFL=2) with Berna Ayat and 
C     Heather Weitzner in 2013
C
C     Numerical wire mesh (ISEDAV=2) to examine stability of different stone
C     sizes on front slope, crest, and back slope with Berna Ayat and Rolando 
C     Garcia in 2013
C
C     Fixed stone structure on movable sand beach (ISTSAN=1) with Heather
C     Weitzner and Rolando Garcia in 2014
C
C     Erosion of sand beach and underlying clay bottom (ICLAY=1) with Heather
C     Weitzner in 2014
C     
C     Dike or dune overflow (IOFLOW=1) with landward SWL (including no standing 
C     water) below dike or dune crest (IWTRAN=1) with Rolando Garcia in 2014
C  
C     Options of IPROFL=2,ISTSAN=1,ICLAY=1 and IOFLOW=1 are still under 
C     development
C     
C     #########################  GENERAL NOTES  ##########################
C     
C     The purpose of each of 22 subroutines arranged in numerical order
C     is described in each subroutine and where it is called.
C     
C     All COMMON statements appear in the Main Program. Description of
C     each COMMON statement is given only in Main Program.
C     
C     #00######################  MAIN PROGRAM  ###########################
C     
C     Main program marches from the offshore boundary node toward the
C     shoreline using subroutines.
C     
C      PROGRAM CSHORE
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
C     NN=maximum number of cross-shore nodes
C     NB=maximum number of offshore wave and water level data
C     NL=maximum number of cross-shore lines
C     
      PARAMETER (NN=5000, NB=30000, NL=100)
      CHARACTER FINMIN*100, VER*70, BASENAME*90 !bdj
      DIMENSION DUMVEC(NN),QTIDE(NB),SMDEDY(NB)
C     
C     ... COMMONs
C     
C     Name    Contents
C     ----------------------------------------------------------------
C     /OPTION/  Computation options and time
C     /PERIOD/  Representative period and input wave angle
C     /SEAWAV/  Input waves and water levels
C     /PREDIC/  Unknown wave variables predicted by CSHORE
C     /BINPUT/  Input bottom geometry
C     /BPROFL/  Discritized bottom geometry
C     /CONSTA/  Constants
C     /LINEAR/  Linear wave values and wave angle quantities
C     /FRICTN/  Dimensionless parameters related to bottom friction
C     /WBREAK/  Wave breaking quantities and constants
C     /CRSMOM/  Terms in cross-shore momentum equation
C     /LOGMOM/  Terms in longshore momentum equation
C     /ENERGY/  Terms in energy or wave action equation
C     /RUNUP/   Parameters for landward computation limit in wet zone
C     /VELOCY/  Mean and standard deviation of horizontal velocities
C     /SEDINP/  Sediment input parameters
C     /SEDOUT/  Sediment output variables
C     /SEDVOL/  Sediment transport volume per unit width
C     /PROCOM/  Beach profile computation variables
C     /ROLLER/  Roller slope,volume flux and related quantities
C     /POROUS/  Porous flow input and output variables
C     /OVERTF/  Wave overtopping and overflow variables
C     /WIND/    Wind speed, direction and shear stresses	
C     /SWASHP/  Swash parameters for wet and dry zone
C     /SWASHY/  Computed swash hydrodynamic variables
C     /WATRAN/  Input landward still water level for IWTRAN=1
C     /COMPAR/  Computational parameters in Subroutines
C     /RRPOND/  Variables for ridge and runnel with ponded water
C     /TIDALC/  Tidal input variables for currents
C     /SERIES/  Time series of wave overtopping and sediment transport rates
C     /VEGETA/  Parameters related to vegetation for IVEG=1 and 2
C     /DIKERO/  Dike erosion variables and parameters for IPROFL=2
C     /WIMESH/  Wire mesh input and variables for ISEDAV=2
C     /STONES/  Variables and input for ISTSAN=1 (stone on sand)
C     /SOCLAY/  Variables and input for ICLAY=1 (sand on clay)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     + ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     + IVEG,ICLAY
      COMMON /PERIOD/ TP,WKPO,ANGLE,WT(NN)
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     + SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /BINPUT/ XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),
     + YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     + SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     + FSX,FSY,FE,QWX,QWY
      COMMON /FRICTN/ GBX(NN),GBY(NN),GF(NN)
      COMMON /WBREAK/ GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN)
      COMMON /CRSMOM/ SXXSTA(NN),TBXSTA(NN)
      COMMON /LOGMOM/ SXYSTA(NN),TBYSTA(NN)
      COMMON /ENERGY/ EFSTA(NN),DFSTA(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /VELOCY/ UMEAN(NN),USTD(NN),USTA(NN),VMEAN(NN),VSTD(NN),
     + VSTA(NN)
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     + EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /SEDOUT/ PS(NN),VS(NN),QSX(NN),QSY(NN),
     + PB(NN),GSLOPE(NN),QBX(NN),QBY(NN),Q(NN)
      COMMON /SEDVOL/ VBX(NN,NL),VSX(NN,NL),VBY(NN,NL),VSY(NN,NL),
     + VY(NN,NL),DZX(NN,NL)
      COMMON /PROCOM/ DELT,DELZB(NN,NL)
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     + WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     + UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /WIND/   W10(NB),WANGLE(NB),WINDCD(NB),TWXSTA(NB),
     + TWYSTA(NB)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     + BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     + USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     + H1,JWD,JDRY
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /COMPAR/ HWDMIN,NPT,NPE
      COMMON /RRPOND/ZW,QD,QM,JXW,JX2,NOPOND
      COMMON /TIDALC/DETADY(NB),DSWLDT(NB)
      COMMON /SERIES/TSQO(NL),TSQBX(NL),TSQSX(NL)
      COMMON /VEGETA/VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
      COMMON /DIKERO/EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN),
     + GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL),GRSR(NN,NL),
     + GRSRD(NN,NL), DEEB, DEEF
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
      COMMON /SOCLAY/EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL),
     + FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL)
C     
C     For iteration convergence 
C     EPS1 = 0.001 for depth (m), height (m) and velocity (m/s)
C     EPS2 = 0.000001 for roller volume flux (m*m/s)
C     MAXITE = 20 for maximum number of iteration
      DATA EPS1, EPS2, MAXITE/1.D-3, 1.D-6, 20/
C     
C     Store the first line of this CSHORE program on ODOC output file
C     ------------------------------------------------------------------
      VER = 'CSHORE USACE version, 2014 last edit 2015-07-06' !bdj 
C      VER = 'CSHORE USACE version, last edit 2015-03-23   ' !bdj
C      VER = 'CSHORE USACE version 2014, merged on 2015-03-12  ' !bdj
C      VER = 'CSHORE USACE version 2011, last edit 2012-08-15 ' !bdj
C      VER = '2014 CSHORE: Rolando ; 2014 February 12'
C     ------------------------------------------------------------------
C     
C     WRITE (*,*) 'Name of Primary Input-Data-File?'
C     READ  (*,5000) FINMIN
C      FINMIN = 'infile'
C bdj 2015-03-23 
C 5000 FORMAT (A12)
      NUM_ARGS = COMMAND_ARGUMENT_COUNT()
      IF (NUM_ARGS.EQ.0) then
         BASENAME = ''
      ELSE
         CALL GET_COMMAND_ARGUMENT(1,BASENAME)
         BASENAME=TRIM(BASENAME)//'.'
      ENDIF
c end bdj 2015-03-23 
C     
C     Subr. 1 OPENER opens input and output files.
C     bdj 2015-03-23  
C     CALL OPENER (FINMIN) 
      CALL OPENER (BASENAME)
C     end bdj 2015-03-23 
C     Subr. 2 INPUT gets input wave and bathymetry information
C     from  the input file, FINMIN.
      CALL INPUT (VER)
C     
C     Subr. 3 BOTTOM computes initial bathymetry at each node.
      CALL BOTTOM
C     Subr. 4 PARAM  calculates constants.
      CALL PARAM
C     
C     ************* START OF TIME MARCHING COMPUTATION ***********
C     
      TIME = 0.D0
      ITIME = 0
C     
C     Subr. 8 OUTPUT stores input before time marching
      DO 1111 L=1,ILINE
        CALL OUTPUT(ITIME,L,0,1)
 1111 CONTINUE
C     
C     NTIME sets of constant wave and water level at the seaward
C     boundary x=0 for all ILINE cross-shore lines
C     
      DO 999 ITIME=1,NTIME
      DO 998 L=1, ILINE
      IF(ANGLE.EQ.AGLINE(L)) THEN
        IANGLE=0
C     normally incident waves along line L
      ELSE
        IANGLE=1
C       IWCINT=0
      ENDIF
C     obliquely incident waves along line L
C     
C     IEND=0 during the constant wave and water level 
C     IEND=1 at the end of each ITIME
C     If IPOND=1, QO=wave overtopping rate at ridge
C     crest and QM=wave overtopping rate at landward end node JMAX
        QO(L)=0.D0
      IF(IPOND.EQ.1) QM=0.D0
C     
 888  IEND=0
C     
C     ..... PREPARE FOR LANDWARD MARCHING COMPUTATION
C     
C     SWLDEP(J,L) = still water depth at node J for the present landward
C     marching computation along cross-shore line L
C     
        ICHECK=0
        DO 90 J=1,JMAX(L)
          SWLDEP(J,L) = SWLBC(ITIME) - ZB(J,L)
          IF(ICHECK.EQ.0) THEN
            IF(SWLDEP(J,L).LT.0.D0)THEN
              JSWL(L) = J
              ICHECK = 1
            ENDIF
          ENDIF
 90     CONTINUE
        IF(ICHECK.EQ.0) JSWL(L)=JMAX(L)
C     If ITIDE=1 and ILAB=0, computed cross-shore tidal 
C     water flux QTIDE at wet node J
        IF(ITIDE.EQ.1) THEN
          DO 91 J=1,JMAX(L)
            SMDEDY(J)=DETADY(ITIME)
C           SMDEDY(J)=SMDEDY(J)*(0.5D0+0.5D0*DTANH((XB(J)-6.D0)/1.D0))
C     where the above transition function is specifically for LSTF
C     pumping system
 91       CONTINUE
          IF(ILAB.EQ.0) THEN
            DO 92 J=1,JMAX(L)
              IF(J.LT.JSWL(L)) THEN
                QTIDE(J)=(XB(JSWL(L))-XB(J))*DSWLDT(ITIME)
              ELSE
                QTIDE(J)=0.D0
              ENDIF
 92         CONTINUE
          ENDIF
        ENDIF
C     
C     If IWTRAN=1 and JCREST(L) is less than JMAX(L), a bay exists landward of
C     the emerged structure or dune crest. The landward still water level
C     is SWLAND(ITIME) for nodes J=JSL,(JSL+1),...,JMAX if ISWLSL=1
        IF(IWTRAN.EQ.1) THEN
          IOFLOW=0
C     If ISWLSL=0, seaward SWL=landward SWL and the entire structure or dune
C     can be submerged below SWL
          IF(ISWLSL.EQ.0.AND.SWLBC(ITIME).GT.RCREST(L)) THEN
            JSL=JMAX(L)
            JSL1=JMAX(L)
            GOTO 94
          ENDIF
C     If ISWLSL=1, landward SWL must be below crest elevation to avoid seaward
C     overflow and landward overflow occurs if seaward SWL is above crest
          IF(ISWLSL.EQ.1) THEN
            IF(SWLAND(ITIME).GE.RCREST(L)) SWLAND(ITIME)=RCREST(L)-1.D-2
            IF(SWLBC(ITIME).GT.RCREST(L)) IOFLOW=1
          ENDIF
C     If ISWLSL=2, no standing water exists lanward of the crest and wet and
C     dry zone extends to the end of the computation domain and overflow occurs
C     if seaward SWL is above crest
          IF(ISWLSL.EQ.2) THEN
            JSL=JMAX(L)
            JSL1=JMAX(L)
            IF(SWLBC(ITIME).GT.RCREST(L)) IOFLOW=1
            GOTO 94
          ENDIF
C     If ISWLSL=0 or 1, landward SWL may intersect landward slope between nodes
C     JSL and JSL1
          IF(JCREST(L).LT.JMAX(L)) THEN
            ICHECK=0
            DO 95 J=(JCREST(L)+1),JMAX(L)
              DUM=SWLAND(ITIME)-ZB(J,L)
              IF(DUM.GT.0.D0) THEN
                SWLDEP(J,L)=DUM
                IF(ICHECK.EQ.0) THEN
                  JSL=J
                  JSL1=JSL-1
                  ICHECK=1
                ENDIF
              ENDIF
 95         CONTINUE
C     If ICHECK=0, no standing water exists lanward of crest
            IF(ICHECK.EQ.0) THEN
              JSL=JMAX(L)
              JSL1=JMAX(L)
            ENDIF
          ENDIF
C
 94       CONTINUE
        ENDIF
C     
C     If IPOND=1, Subr.20 PONDED finds ponded water zone
        IF(IPOND.EQ.1) THEN
          CALL PONDED(L)
        ENDIF
C     
C.....ITERATION TO FIND QO(L) ............................................
C     At beginning of each ITIME, QO(L)=0.0 as specified above.
C     During each ITIME for profile evolution computation
C     with IPROFL=1, converged QO(L) is used as an initial quess
C     for the next profile change computation with ITEMAX=4
C     If IOVER=0, QO(L)=0.0 always and no iteration.
        IF(IOVER.EQ.0) THEN
          ITEMAX=1
        ELSE
          ITEMAX=20
C     Computed wave overtopping rates QO(L) for ITEMAX=10-30 changed
C     very little for fixed coastal structures with IPROFL=0
          IF(IPROFL.GE.1) ITEMAX=4
C     Computed overwashed dune profile evolutions changed very
C     little for ITEMAX=3-4.
        ENDIF
C     
        ITEQO=0
 777    ITEQO=ITEQO+1
C     
        SIGMA(1) = HRMS(1)/SQR8
        H(1) = WSETUP(1) + SWLDEP(1,L)
        SIGSTA(1) = SIGMA(1)/H(1)
C     
C     Subr.5 LWAVE returns linear wave number WKP,phase velocity CP(J)
C     ratio WN(J) of group velocity to phase velocity and 
C     sin STHETA(J) and cos CTHETA(J) of  wave angle for given
C     QDISP=water flux in dispersion relation of linear waves.
C     QDISP=0.0 is assumed for J=1 for simplicity
        QDISP=0.D0
        CALL LWAVE(1,L,H(1),QDISP)
C     
C     Tentatively assume VMEAN(1) = 0.0
        VMEAN(1) = 0.D0
        VSIGT = 0.D0
C     At node J=1, no porous layer
        QWX=QO(L)
        IF(ITIDE.EQ.1.AND.ILAB.EQ.0) QWX=QWX+QTIDE(1)
        IF(IPERM.EQ.1) THEN
          UPMEAN(1)=0.D0
          QP(1)=0.D0
          UPSTD(1)=0.D0
          DPSTA(1)=0.D0
        ENDIF
C     
        IF(IROLL.EQ.1) RQ(1)=0.D0
        SIGMA2 =  SIGMA(1)**2.D0
        QWY=GRAV*SIGMA2*STHETA(1)/CP(1)
        SXXSTA(1) = SIGMA2*FSX
        EFSTA(1) = SIGMA2*FE
        IF(IANGLE.EQ.1) SXYSTA(1) = SIGMA2*FSY
        IF(IWCINT.EQ.1) THEN
          DUM=GRAV*H(1)
          SXXSTA(1)=SXXSTA(1)+QWX**2.D0/DUM
          IF(IANGLE.EQ.1) SXYSTA(1)=SXYSTA(1)+QWX*QWY/DUM
        ENDIF
C     
C     where roller volume flux RQ(1)=0 is assumed for SXXSTA(1),
C     SXYSTA(1), USIGT=UMEAN/SIGT, and QWY
C     
C     Subr.6a GBXAGF returns approximate analytical values for
C     the Gbx and Gf factors used in calculating cross-shore
C     bottom shear stress and energy dissipation.
C     Effect of QWX on USIGT is neglected unless IWCINT=1
C     If bottom friction coefficient is positive,
        IF(FB2(1,L).GT.0.D0) THEN
          USIGT = -SIGSTA(1)*GRAV*H(1)/CP(1)/CP(1)
          IF(IANGLE.EQ.1) USIGT = USIGT*CTHETA(1)
          DUM = SIGSTA(1)*CP(1)
          IF(IWCINT.EQ.1) THEN
            IF (DUM.GT.1.D-10) USIGT = USIGT+QWX/H(1)/DUM !bdj
          ENDIF
          CALL GBXAGF(CTHETA(1),USIGT,STHETA(1),VSIGT,GBX(1),GF(1))
          TBXSTA(1) = FB2(1,L)*GBX(1)*DUM**2.D0/GRAV
          DFSTA(1) = FB2(1,L)*GF(1)*DUM**3.D0/GRAV
          IF(IVEG.GE.1) THEN
            DUM=VEGH(1,L)
            IF(DUM.GT.H(1)) DUM=H(1)
            VEGCV=1.D0+DUM*VEGFB(1,L)
            TBXSTA(1)=VEGCV*TBXSTA(1)
            DFSTA(1)=VEGCV*DFSTA(1) 
          ENDIF
        ELSE
          TBXSTA(1) = 0.D0
          DFSTA(1) = 0.D0
        ENDIF
C     
C     Subr.7 DBREAK computes the fraction of breaking waves and
C     the associated wave energy dissipation and returns DBSTA(1).
        CALL DBREAK(1, L, HRMS(1), H(1))
C     
C     ------------ LANDWARD MARCHING COMPUTATION -----------------------
C     
C     Computation marching landward from seaward boundary, J = 1
C     Compute unknown variables at node JP1=(J+1) along line L.
        J = 0
 100    J = J + 1
        JP1 = J + 1
        ITE=0
C     
        DUM=DFSTA(J)+DBSTA(J)
        IF(IPERM.EQ.1) DUM=DUM+DPSTA(J)
        DUM=DUM*WT(J)
        DUM=(EFSTA(J)-DX*DUM)/FE
        IF(DUM.LE.0.D0) THEN
          WRITE(40,2901) JP1,L,TIME,ITEQO,ITE,QO(L)
C     Accept the computed results up to node J and end landward
C     marching computation (go to 400)
          JP1=JP1-1
          GOTO 400
        ENDIF
 2901   FORMAT(/'END OF LANDWARD MARCHING: '/
     +   'Square of sigma SIGTIE is negative at node ',I6,'Line=',I3 /
     +   ' TIME =', F13.3,' ITEQO=',I2,' ITE=',I2,' QO=',F13.9)
C     
        SIGITE  = DSQRT(DUM)
        SXXSTA(JP1) = FSX*SIGITE**2.D0
        IF(IROLL.EQ.1) SXXSTA(JP1) = SXXSTA(JP1) + RX(J)*RQ(J)
        IF(IWCINT.EQ.1) SXXSTA(JP1)=SXXSTA(JP1)+QWX*QWX/GRAV/H(J)
        WSETUP(JP1) = WSETUP(J)-(SXXSTA(JP1)-SXXSTA(J)+
     +   (TBXSTA(J)-TWXSTA(ITIME))*DX)/H(J)
        HITE = WSETUP(JP1) + SWLDEP(JP1,L)
C     
        IF(HITE.LT.EPS1) THEN
          WRITE(40,2902) JP1,L,TIME,ITEQO,QO(L)
          JP1=JP1-1
          GOTO 400
        ENDIF
 2902   FORMAT(/'END OF LANDWARD MARCHING: '/
     +   'Water depth is less than EPS1 at node ',I6,'Line=',I3 /
     +   ' TIME =',F13.3,' ITEQO =',I2,' QO =',F13.9)
C     
        QWX=QO(L)
        IF(IPERM.EQ.1) QWX=QO(L)-QP(J)
        IF(ITIDE.EQ.1.AND.ILAB.EQ.0) QWX=QWX+QTIDE(JP1)
        IF(IWCINT.EQ.1) THEN
          IF(IANGLE.EQ.0) THEN
            QDISP=QWX
          ELSE
            QWY = HITE*VMEAN(J) + GRAV*SIGITE**2.D0*STHETA(J)/CP(J)
            IF(IROLL.EQ.1) QWY=QWY+RQ(J)*STHETA(J)
            QDISP = QWX*CTHETA(J) + QWY*STHETA(J)
          ENDIF
        ENDIF
        CALL LWAVE(JP1,L,HITE,QDISP)
C     
        IF(IANGLE.EQ.1) THEN
          DUM1 = SIGITE**2.D0
          SXYSTA(JP1) = FSY*DUM1
          IF(IROLL.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+RY(J)*RQ(J)
          IF(IWCINT.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+QWX*QWY/GRAV/HITE
          DUM2 = SXYSTA(JP1) - SXYSTA(J)
          SIGN=STHETA(JP1)*DUM2
          IF(SIGN.GT.0.D0) DUM2=0.D0
          TBYSTA(JP1) = -DUM2/DX + TWYSTA(ITIME)
          IF(ITIDE.EQ.1) TBYSTA(JP1)=TBYSTA(JP1)-HITE*SMDEDY(JP1)
          DUM = SIGITE/HITE
          IF(DUM.GT.SISMAX) DUM = SISMAX
          DUM3 = CP(JP1)*CP(JP1)/GRAV
          GBY(JP1) = TBYSTA(JP1)/FB2(JP1,L)/DUM3/DUM/DUM
          IF(IVEG.GE.1) THEN
            DUMH=VEGH(JP1,L)
            IF(DUM.GT.HITE) DUMH=HITE
            VEGCV=1.D0+DUMH*VEGFB(JP1,L)
            GBY(JP1)=GBY(JP1)/VEGCV
          ENDIF

C     Subr. 6b VSTGBY computes VSIGT for specified GBY, CTHETA, USIGT
C     and STHETA where effect of QWX on USIGT is neglected unless IWCINT=1
          USIGT = -CTHETA(J)*DUM*HITE/DUM3
          IF(IROLL.EQ.1) THEN
            USIGT = USIGT*(1.D0+ (CP(JP1)/GRAV)*RQ(J)/SIGITE**2.D0)
          ENDIF
          SIGT = DUM*CP(JP1)
          IF(IWCINT.EQ.1) USIGT=USIGT+QWX/HITE/SIGT
          CALL VSTGBY(CTHETA(J),USIGT,STHETA(J),VSIGT,GBY(JP1))
          VITE = VSIGT*SIGT
        ENDIF
C     
        IF(IROLL.EQ.1) THEN
          RQITE = RQ(J) + DX*(DBSTA(J)-RBETA(J)*RQ(J))/RE(J)
          IF(RQITE.LT.0.D0) RQITE=0.D0
        ENDIF
C     
C******Begin iteration for improved Euler finite difference method****
C     
        DO 200 ITE = 1, MAXITE
C     
          HRMITE = SIGITE*SQR8
C     
          CALL DBREAK(JP1,L,HRMITE, HITE)
          SIGSTA(JP1) = SIGITE/HITE
          IF(SIGSTA(JP1).GT.SISMAX) SIGSTA(JP1) = SISMAX
C     
          SIGT = CP(JP1)*SIGSTA(JP1)
          IF(IANGLE.EQ.0) THEN
            VSIGT = 0.D0
          ELSE
            VSIGT = VITE/SIGT
          ENDIF
C     
C     If IPERM=1, Subr.9 POFLOW computes porous flow variables.
C     UPMEAN(J) = mean of horizontal discharge velocity UP
C     UPSTD(J) = standard deviation of UP
C     DPSTA(J) = energy dissipation rate of porous flow
          QWX=QO(L)
          IF(IPERM.EQ.1) THEN
            PKHSIG = WKP*HITE*SIGSTA(JP1)
            DEDX = (WSETUP(JP1) - WSETUP(J))/DX
            CALL POFLOW(JP1,L,PKHSIG,DEDX)
            QWX = QO(L) - QP(JP1)
          ENDIF
          IF(ITIDE.EQ.1.AND.ILAB.EQ.0) QWX=QWX+QTIDE(JP1)
C     
          IF(FB2(JP1,L).GT.0.D0) THEN
            DUM = GRAV*HITE/CP(JP1)/CP(JP1)
            USIGT = -CTHETA(JP1)*SIGSTA(JP1)*DUM
            IF(IROLL.EQ.1) THEN
              USIGT = USIGT*(1.D0+(CP(JP1)/GRAV)*RQITE/SIGITE**2.D0)
            ENDIF
            IF(IWCINT.EQ.1) USIGT=USIGT+QWX/HITE/SIGT
            CALL GBXAGF(CTHETA(JP1),USIGT,STHETA(JP1),VSIGT,
     +      GBX(JP1), GF(JP1))
            TBXSTA(JP1) = FB2(JP1,L)*GBX(JP1)*SIGT**2.D0/GRAV
            DFSTA(JP1) = FB2(JP1,L)*GF(JP1)*SIGT**3.D0/GRAV
            IF(IVEG.GE.1) THEN
              DUM=VEGH(JP1,L)
              IF(DUM.GT.HITE) DUM=HITE
              VEGCV=1.D0+DUM*VEGFB(JP1,L)
              TBXSTA(JP1)=VEGCV*TBXSTA(JP1)
              DFSTA(JP1)=VEGCV*DFSTA(JP1)
            ENDIF
          ELSE
            TBXSTA(JP1) = 0.D0
            DFSTA(JP1) = 0.D0
          ENDIF
C     
          DUMD = DFSTA(JP1) + DFSTA(J) + DBSTA(JP1) + DBSTA(J)
          IF(IPERM.EQ.1) DUMD=DPSTA(JP1)+DPSTA(J)+DUMD
          DUMD = DUMD*(WT(J)+WT(JP1))/2.D0
          DUM = (EFSTA(J) - DXD2*DUMD)/FE
          IF(DUM.LE.0.D0) THEN 
            WRITE(40,2901) JP1, L, TIME, ITEQO, ITE, QO(L)
C     Accept the computed results up to node J
            JP1=JP1-1
            GOTO 400
          ELSE
            SIGMA(JP1) = DSQRT(DUM)
          ENDIF
C     
          SXXSTA(JP1) = FSX*SIGMA(JP1)**2.D0
          IF(IROLL.EQ.1) SXXSTA(JP1)=SXXSTA(JP1)+RX(JP1)*RQITE
          IF(IWCINT.EQ.1) SXXSTA(JP1)=SXXSTA(JP1)+QWX*QWX/GRAV/HITE
          WSETUP(JP1) = WSETUP(J) - (2.D0* (SXXSTA(JP1)-SXXSTA(J)) +
     +       DX*(TBXSTA(JP1)+TBXSTA(J)-2.D0*TWXSTA(ITIME)))/
     +       (HITE+H(J))
          H(JP1) = WSETUP(JP1) + SWLDEP(JP1,L)
          SIGSTA(JP1) = SIGMA(JP1)/H(JP1)
          IF(SIGSTA(JP1).GT.SISMAX) SIGSTA(JP1)=SISMAX
C     
          IF(H(JP1).LE.EPS1) THEN 
            WRITE(40,2902) JP1, L, TIME, ITEQO, QO(L)
            JP1 = JP1-1
            GOTO 400
          ENDIF 
C     
          IF(IWCINT.EQ.1) THEN
            IF(IANGLE.EQ.0) THEN
              QDISP = QWX
            ELSE
              QWY=H(JP1)*VITE+
     +          GRAV*SIGMA(JP1)**2.D0*STHETA(JP1)/CP(JP1)
              IF(IROLL.EQ.1) QWY=QWY+RQITE*STHETA(JP1)
              QDISP = QWX*CTHETA(JP1) + QWY*STHETA(JP1)
            ENDIF
          ENDIF
          CALL LWAVE(JP1,L,H(JP1),QDISP)
C     
          IF(IANGLE.EQ.1) THEN
            DUM1 = SIGMA(JP1)**2.D0
            SXYSTA(JP1) = FSY*DUM1
            IF(IROLL.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+RY(JP1)*RQITE
            IF(IWCINT.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+QWX*QWY/GRAV/
     +         H(JP1)
            DUM2 = SXYSTA(JP1) - SXYSTA(J)
            SIGN=STHETA(JP1)*DUM2
            IF(SIGN.GT.0.D0) DUM2=0.D0
            TBYSTA(JP1) = -DUM2/DX + TWYSTA(ITIME)
            IF(ITIDE.EQ.1) TBYSTA(JP1)=TBYSTA(JP1)-H(JP1)*SMDEDY(JP1)
            DUM3 = CP(JP1)*CP(JP1)/GRAV  
            GBY(JP1)=TBYSTA(JP1)/FB2(JP1,L)/DUM3/
     +         SIGSTA(JP1)/SIGSTA(JP1)
            IF(IVEG.GE.1) THEN
              DUM=VEGH(JP1,L)
              IF(DUM.GT.H(JP1)) DUM=H(JP1)
              VEGCV=1.0D0+DUM*VEGFB(JP1,L)
              GBY(JP1)=GBY(JP1)/VEGCV
            ENDIF
            USIGT = -CTHETA(JP1)*SIGSTA(JP1)*H(JP1)/DUM3
            IF(IROLL.EQ.1) THEN
              USIGT=USIGT*(1.D0+(CP(JP1)/GRAV)*RQITE/SIGMA(JP1)**2.D0)
            ENDIF
            SIGT = SIGSTA(JP1)*CP(JP1)
            IF(IWCINT.EQ.1) USIGT=USIGT+QWX/H(JP1)/SIGT
            CALL VSTGBY(CTHETA(JP1),USIGT,STHETA(JP1),VSIGT,GBY(JP1))
            VMEAN(JP1) = VSIGT*SIGT
          ENDIF
C     
          IF(IROLL.EQ.1) THEN
            DUM1 = RE(JP1) + DXD2*RBETA(JP1)
            DUM2 = (RE(J) - DXD2*RBETA(J))*RQ(J) + 
     +        DXD2*(DBSTA(JP1) + DBSTA(J))
            RQ(JP1) = DUM2/DUM1
          ENDIF
C     
C     Check for convergence
C     
          ESIGMA = DABS(SIGMA(JP1) - SIGITE)
          EH = DABS(H(JP1) - HITE)
          IF(IANGLE.EQ.1) EV = DABS(VMEAN(JP1) - VITE)
          IF(IROLL.EQ.1) ERQ = DABS(RQ(JP1) - RQITE)
          IF(ESIGMA.LT.EPS1.AND.EH.LT.EPS1) THEN
            IF(IANGLE.EQ.0) THEN
              GOTO 199
            ELSE
              IF(EV.LT.EPS1) GOTO 199
              GOTO 198
            ENDIF
 199        IF(IROLL.EQ.0) THEN
              GOTO 210
            ELSE
              IF(ERQ.LT.EPS2) GOTO 210
              GOTO 198
            ENDIF
          ENDIF
C     
C     Average new and previous values to accelerate convergence
 198      SIGITE =  0.5D0*(SIGMA(JP1) + SIGITE)
          HITE = 0.5D0*(H(JP1) + HITE)
          IF(IANGLE.EQ.1) VITE = 0.5D0*(VMEAN(JP1) + VITE)
          IF(IROLL.EQ.1) RQITE = 0.5D0*(RQ(JP1)+RQITE)
C     
 200    CONTINUE
C     
C*****End of iteration: DO 200 ITE = 1 to MAXITE*********************
C     
C     The iteration did not converge
        WRITE(40,2903) MAXITE, EPS1, JP1, L, TIME, QO(L)
C     Adopt the last iteration values
        SIGMA(JP1) = SIGITE
        H(JP1) = HITE
        IF(IANGLE.EQ.1) VMEAN(JP1) = VITE
        IF(IROLL.EQ.1) RQ(JP1) = RQITE
 2903   FORMAT(/'WARNING: Convergence was not reached after MAXITE= ',
     +     I4/ ' iterations with relative error EPS1 = ',E17.5/
     +     'at node JP1 = ',I4, ' Line L=',I3, ' TIME= ',F13.3/
     +     'QO(L)=',F13.9)
C     
 210    HRMS(JP1) = SQR8*SIGMA(JP1)
        WSETUP(JP1) = H(JP1) - SWLDEP(JP1,L)
C     
        IF(IWCINT.EQ.1) THEN
          IF(IANGLE.EQ.0) THEN
            QDISP = QWX
          ELSE
            QWY = H(JP1)*VMEAN(JP1) + GRAV*SIGMA(JP1)**2.D0*
     +       STHETA(JP1)/CP(JP1)
            IF(IROLL.EQ.1) QWY=QWY + RQ(JP1)*STHETA(JP1)
            QDISP = QWX*CTHETA(JP1) + QWY*STHETA(JP1)
          ENDIF
        ENDIF
C     
        CALL LWAVE(JP1, L, H(JP1), QDISP)
        CALL DBREAK(JP1, L, HRMS(JP1), H(JP1))
        SIGSTA(JP1) = SIGMA(JP1)/H(JP1)
        IF(SIGSTA(JP1).GT.SISMAX) SIGSTA(JP1) = SISMAX
        SIGT = SIGSTA(JP1)*CP(JP1)
        IF(IANGLE.EQ.0) THEN
          VSIGT = 0.D0
        ELSE
          VSIGT = VMEAN(JP1)/SIGT
        ENDIF
C     
        QWX=QO(L)
        IF(IPERM.EQ.1) THEN
          PKHSIG = WKP*H(JP1)*SIGSTA(JP1)
          DEDX = (WSETUP(JP1) - WSETUP(J))/DX
          CALL POFLOW(JP1,L,PKHSIG,DEDX)
          QWX = QO(L) - QP(JP1)
        ENDIF
        IF(ITIDE.EQ.1.AND.ILAB.EQ.0) QWX=QWX+QTIDE(JP1)
C     
        SIGMA2 = SIGMA(JP1)**2.D0
        SXXSTA(JP1) = SIGMA2*FSX
        IF(IROLL.EQ.1) SXXSTA(JP1)=SXXSTA(JP1)+RX(JP1)*RQ(JP1)
        IF(IWCINT.EQ.1) SXXSTA(JP1)=SXXSTA(JP1)+QWX*QWX/GRAV/H(JP1)
        EFSTA(JP1) = SIGMA2*FE
        DUM3 = CP(JP1)*CP(JP1)/GRAV
        IF(FB2(JP1,L).GT.0.D0) THEN
          USIGT = -CTHETA(JP1)*SIGSTA(JP1)*H(JP1)/DUM3
          IF(IROLL.EQ.1) THEN 
            USIGT = USIGT*(1.D0+(CP(JP1)/GRAV)*RQ(JP1)/SIGMA2)
          ENDIF
          IF(IWCINT.EQ.1) USIGT=USIGT+QWX/H(JP1)/SIGT
          CALL GBXAGF(CTHETA(JP1),USIGT,STHETA(JP1),VSIGT,GBX(JP1),
     +       GF(JP1))
          TBXSTA(JP1)=FB2(JP1,L)*GBX(JP1)*SIGT**2.D0/GRAV
          DFSTA(JP1)=FB2(JP1,L)*GF(JP1)*SIGT**3.D0/GRAV
          IF(IVEG.GE.1) THEN
            DUM=VEGH(JP1,L)
            IF(DUM.GT.H(JP1)) DUM=H(JP1)
            VEGCV=1.0D0+DUM*VEGFB(JP1,L)
            TBXSTA(JP1)=VEGCV*TBXSTA(JP1)
            DFSTA(JP1)=VEGCV*DFSTA(JP1)
          ENDIF
        ELSE
          TBXSTA(JP1) = 0.D0
          DFSTA(JP1) = 0.D0
        ENDIF
C     
        IF(IANGLE.EQ.1) THEN
          SXYSTA(JP1) = FSY*SIGMA2
          IF(IROLL.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+RY(JP1)*RQ(JP1)
          IF(IWCINT.EQ.1) SXYSTA(JP1)=SXYSTA(JP1)+QWX*QWY/GRAV/H(JP1)
          DUM2 = SXYSTA(JP1) - SXYSTA(J)
          SIGN=STHETA(JP1)*DUM2
          IF(SIGN.GT.0.D0) DUM2=0.D0
          TBYSTA(JP1) = -DUM2/DX + TWYSTA(ITIME)
          IF(ITIDE.EQ.1) TBYSTA(JP1)=TBYSTA(JP1)-H(JP1)*SMDEDY(JP1)
          IF(J.EQ.1) THEN
            TBYSTA(J) = TBYSTA(JP1)
            VMEAN(J) = VMEAN(JP1)
          ENDIF
          GBY(JP1) = TBYSTA(JP1)/FB2(JP1,L)/DUM3/SIGSTA(JP1)/SIGSTA(JP1)
          IF(IVEG.GE.1) GBY(JP1)=GBY(JP1)/VEGCV
        ENDIF
C     
        JDUM = JMAX(L)
        IF(RCREST(L).GT.SWLBC(ITIME)) JDUM=JCREST(L)
C     If IWTRAN=1 and IOFLOW=1, overflow occurs on submerged crest
        IF(IWTRAN.EQ.1.AND.IOFLOW.EQ.1) JDUM=JCREST(L)
        IF(H(JP1).LT.EPS1.OR.JP1.EQ.JDUM) GOTO 400
C      
        GOTO 100
C     
C----------------End of LANDWARD MARCHING COMPUTATION -------------
C     
 400    CONTINUE
C     
        JR = JP1
        XR = XB(JR)
        ZR = ZB(JR,L)
C     
C BDJ 2011->2014 on 2014-10-02 
          CALL SRFSP(L)
C end BDJ 2011->2014 on 2014-10-02 

C     If IOVER=1, Subr.10 QORATE computes for cross-shore line L
C     QO(L) = sum of wave overtopping, overflow and seepage rates
C     If IOVER=0, QO(L)=0.0, no iteration and no wet/dry zone
C     If IWTRAN=1 and JR=JMAX(L), no wet and dry zone in computattion domain
C     Assume no water flux at landward end of computation domain
        IF(IOVER.EQ.1) THEN 
          IF(IWTRAN.EQ.1.AND.JR.EQ.JMAX(L)) THEN
            JWD=JR
            JDRY=JR
            QO(L)=0.D0
            GOTO 405
          ELSE
            ICONV = 1
            QOUSED = QO(L)
            CALL QORATE(ITIME,L,ITEQO,ICONV,0)
            IF(ICONV.EQ.0) GOTO 405
            IF(ICONV.EQ.1) WRITE(40,2904) JR,L,TIME,ITEQO,QOUSED,QO(L)
          ENDIF
        ELSE
          JWD = JR
          JDRY = JR
          QO(L)=0.D0
          GOTO 405
        ENDIF
 2904   FORMAT(/'NO CONVERGENCE OF QO ITERATION'/
     +     'Landward end node JR=', I6,' Line=',I3,' TIME=',F13.3/
     +     'Iteration number ITEQO=',I3,'   assumed QO=',F13.9/
     +     'computed QO=',F13.9)
C     
        IF(ITEQO.LT.ITEMAX) GOTO 777
C     
C....................END OF QO ITERATION........................... 
 405    CONTINUE
C     
C     Calculate the standard deviation and mean of the horizontal
C     velocities U and V
C     
        DO 410 I = 1,JR
          SIGT = CP(I)*SIGSTA(I)
          USTD(I) = SIGT*CTHETA(I)
          UMEAN(I)= -USTD(I)*SIGSTA(I)*GRAV*H(I)/CP(I)/CP(I)
          IF(IROLL.EQ.1) UMEAN(I)=UMEAN(I)*(1.D0+(CP(I)/GRAV)*
     +       RQ(I)/SIGMA(I)**2.D0)
          QWX = QO(L)
          IF(IPERM.EQ.1) QWX=QO(L)-HP(I,L)*UPMEAN(I)
          IF(ITIDE.EQ.1.AND.ILAB.EQ.0) QWX=QWX+QTIDE(I)
          UMEAN(I) = UMEAN(I) + QWX/H(I)
          IF(SIGT.GT.1.D-10) THEN !bdj
            USTA(I)=UMEAN(I)/SIGT
          ELSE
            USTA(I)=0.D0
          ENDIF
          IF(IANGLE.EQ.1) THEN
            VSTD(I) = SIGT*DABS(STHETA(I))
            VSTA(I) = VMEAN(I)/SIGT
          ENDIF
 410    CONTINUE
C     
C     If IOVER=1, connect H(J) and UMEAN(J) with J=1 to JR with wet/dry-
C     zone HWD(J) and UMEAWD(J) with J=JWD to JDRY using Subr.17 TRANWD
C     also connect the corresponding standard deviations.
        IF(IOVER.EQ.1) THEN
          PWET(1:JWD)=1.D0
          IF(JDRY.GT.JR) THEN
            CALL TRANWD(H,JR,HWD,JWD,JDRY)
            CALL TRANWD(SIGMA,JR,SIGWD,JWD,JDRY)
            CALL TRANWD(UMEAN,JR,UMEAWD,JWD,JDRY)
            CALL TRANWD(USTD,JR,USTDWD,JWD,JDRY)
            IF(IPERM.EQ.1) CALL TRANWD(UPMEAN,JR,UPMWD,JWD,JDRY)
            IF(IANGLE.EQ.1) THEN
              CALL TRANWD(VMEAN,JR,VMEAWD,JWD,JDRY)
              CALL TRANWD(VSTD,JR,VSTDWD,JWD,JDRY)
            ENDIF
          ELSE
            JDRY=JR
            IF(JWD.LT.JR) THEN
              JDUM=JWD+1
              PWET(JDUM:JR)=1.D0
            ENDIF
          ENDIF
        ENDIF
C     Smooth computed H(J), SIGMA(J), USTD(J), UMEAN(J), USTA(J), DFSTA(J),
C     DBSTA(J),RQ(J), VMEAN(J), VSTD(J) and VSTA(J) using Subr. 14 SMOOTH
        DUMVEC = H
        CALL SMOOTH(JDRY,DUMVEC,H)
        DUMVEC = SIGMA
        CALL SMOOTH(JDRY,DUMVEC,SIGMA)
        DUMVEC = USTD
        CALL SMOOTH(JDRY,DUMVEC,USTD)
        DUMVEC = UMEAN	
        CALL SMOOTH(JDRY,DUMVEC,UMEAN)
        DUMVEC = USTA
        CALL SMOOTH(JR,DUMVEC,USTA)
        DUMVEC = DFSTA
        CALL SMOOTH(JR,DUMVEC,DFSTA)
        IF(IPERM.EQ.1) THEN
          DUMVEC=UPMEAN
          CALL SMOOTH(JDRY,DUMVEC,UPMEAN)
          IF(IOVER.EQ.1) THEN
            DO 420 J=2,JDRY
              DUM=ZP(J,L)
              IF(DUM.LT.SWLBC(ITIME).AND.ZP(J,L).GE.ZP(J-1,L)) 
     +           DUM=SWLBC(ITIME)
              ETAPOR=ZB(J,L)*PWET(J)+DUM*(1.D0-PWET(J))
              QP(J)=UPMEAN(J)*(ETAPOR-ZP(J,L))*PWET(J)
 420        CONTINUE
          ENDIF
        ENDIF
        IF(IROLL.EQ.0) THEN
          DUMVEC=DBSTA
          CALL SMOOTH(JR,DUMVEC,DBSTA)
        ELSE 
          DUMVEC = RQ
          CALL SMOOTH(JR,DUMVEC,RQ)
        ENDIF
        IF(IANGLE.EQ.1) THEN
          DUMVEC = VMEAN
          CALL SMOOTH(JDRY,DUMVEC,VMEAN)
          DUMVEC = VSTD
          CALL SMOOTH(JDRY,DUMVEC,VSTD)
          DUMVEC = VSTA
          CALL SMOOTH(JR,DUMVEC,VSTA)
        ENDIF
C     
C     Subr. 21 WTRANS computes transmitted waves (IWTRAN=1) landward of
C     an emerged structrue or barrier island if entire structure is not 
C     submerged and standing water exists
        IF(IWTRAN.EQ.1.AND.JR.LT.JMAX(L)) THEN
          ICHECK=0
          IF(JSWL(L).EQ.JMAX(L).AND.IOFLOW.EQ.0) ICHECK=1
          JEND=JSL1
          IF(ICHECK.EQ.1) JEND=JMAX(L)
          IF(JDRY.LT.JEND) THEN
            JDUM=JDRY+1
            DO 425 J=JDUM,JEND
              PWET(J)=0.D0
              H(J)=0.D0
              IF(ISWLSL.LE.1) THEN
                IF(IOFLOW.EQ.0.OR.JSL.LT.JMAX(L)) THEN
                  IF (SWLDEP(J,L).GT.0.D0) THEN
                    PWET(J)=1.D0
                    H(J)=SWLDEP(J,L)
                  ENDIF
                ENDIF
              ENDIF
              SIGMA(J)=0.D0
              WSETUP(J)=0.D0
              SIGSTA(J)=0.D0
              UMEAN(J)=0.D0
              USTD(J)=0.D0
              IF(IPERM.EQ.1) THEN
                QP(J)=QO(L)
                IF(ICHECK.EQ.1) QP(J)=0.D0
                IF(HP(J,L).GT.1.D-3) THEN
                  UPMEAN(J)=QP(J)/HP(J,L)
                ELSE
                  UPMEAN(J)=0.D0
                ENDIF
              ENDIF
              IF(IANGLE.EQ.1) THEN
                VMEAN(J)=0.D0
                VSTD(J)=0.D0
                STHETA(J)=0.D0
              ENDIF
 425        CONTINUE
          ENDIF
C     IF IWTRAN=1 and JDRY is less than JSL1=(JSL-1), no water at
C     nodes between JDRY and JSL but water flux QP(J) in permeable
C     layer is assumed constant
          IF(ICHECK.EQ.0.AND.JSL.LT.JMAX(L)) CALL WTRANS(ITIME,L)
        ENDIF
C     
C     Subr.11 SEDTRA computes the cross-shore and longshore 
C     sediment transport rates if IPROFL = 1
C     If a vertical wall exists, IVWALL=2 indicates exposure to
C     wave action along cross-shore line L
        IF(IPROFL.EQ.1) THEN
            IF(IVWALL(L).GE.1) THEN
              IF(ZB(JMAX(L)-1,L).LT.ZP(JMAX(L),L)) THEN
                IVWALL(L)=2
              ELSE
                IVWALL(L)=1
              ENDIF
            ENDIF
          CALL SEDTRA(L)
        ENDIF
C     
C     If IPROFL=1, Subr.12 CHANGE computes the change of the 
C     bottom elevation, DELZB(j), at node j during the time step
C     DELT determined in this subroutine which also checks whether
C     IEND=1 and the end of given ITIME is reached.
C     VY(J,L)=total longshore sediment transport rate integrated from
C     TIMEBC(ITIME) to TIMEBC(ITIME+1) used in Subr.12 CHANGE if IQYDY=1
C     If ISTSAN=1,ZB(J,L) and ZP(J,L) at next time level are computed in Subr. 12 CHANGE
C     If ICLAY=1,Subr. 22 EROSON is called to compute erosion of clay below sand
C
        IF(IPROFL.EQ.1) THEN
          CALL CHANGE(ITIME,L,IEND,1)
        DO 430 J=1,JMAX(L)
          IF(TIME.EQ.TIMEBC(ITIME)) THEN
            VY(J,L)=0.D0
            DZX(J,L)=ZB(J,L)
          ENDIF
          IF(ISTSAN.EQ.0) ZB(J,L)=ZB(J,L)+DELZB(J,L)
          IF(ICLAY.EQ.1) CALL EROSON(ITIME,L,IEND)
          IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
            HP(J,L)=ZB(J,L)-ZP(J,L)
            IF(HP(J,L).LT.0.D0) THEN
              HP(J,L)=0.D0
              ZB(J,L)=ZP(J,L)
            ENDIF
            IF(ISEDAV.EQ.2) THEN
              IF(ZB(J,L).LT.ZMESH(J,L)) ZB(J,L)=ZMESH(J,L)
            ENDIF
          ENDIF
          IF(TIME.EQ.0.D0) THEN
            VBX(J,L)=0.D0
            VSX(J,L)=0.D0
            VBY(J,L)=0.D0
            VSY(J,L)=0.D0
          ENDIF
          VBX(J,L)=VBX(J,L)+DELT*QBX(J)
          VSX(J,L)=VSX(J,L)+DELT*QSX(J)
          IF(IANGLE.EQ.1) THEN
            VBY(J,L)=VBY(J,L)+DELT*QBY(J)
            VSY(J,L)=VSY(J,L)+DELT*QSY(J)
            VY(J,L)=VY(J,L)+DELT*(QBY(J)+QSY(J))
          ENDIF
 430    CONTINUE
          IF(IVWALL(L).GE.1) THEN
            IF(ZB(JMAX(L),L).LT.ZP(JMAX(L),L)) ZB(JMAX(L),L)=ZP(JMAX(L)
     +        ,L)
          ENDIF
        ENDIF
        IF(IPROFL.EQ.0) THEN
          IEND=1
          DELT = 1.D0
        ENDIF
        IF(IQYDY.EQ.1.AND.IEND.EQ.1) THEN
          IF(L.EQ.ILINE) CALL CHANGE(ITIME,L,IEND,2)
        ENDIF
C    
C     If IVEG=1 and IPROFL=1, compute the vegetation height 
C     VEGH(J,L) above the local bottom ZB(J,L) in vegetated zone 
C     with UPROOT(J,L)=1.0. Check whether the vegetation is buried 
C     or uprooted using the fixed upper and lower elevations of the 
C     vegetation denoted as VEGZD(J,L) and VEGZR(J,L)
C     
      IF(IVEG.EQ.1.AND.IPROFL.EQ.1) THEN
        DO 440 J=1,JMAX(L)
          IF(UPROOT(J,L).EQ.1.D0) THEN
            VEGH(J,L)=VEGZD(J,L)-ZB(J,L)
            IF(VEGH(J,L).LT.0.D0) VEGH(J,L)=0.D0
            IF(VEGZR(J,L).GE.ZB(J,L)) THEN
              UPROOT(J,L)=0.D0
              VEGH(J,L)=0.D0
            ENDIF
          ENDIF
 440    CONTINUE
      ENDIF
C     where for buried vegetation, VEGH(J,L)=0.0 and UPROOT(J,L)=1.0,
C     whereas for uprooted vegetation, VEGH(J,L)=0.0 and UPROOT(J,L)=0.0
C
C     If IPROFL=2, Subr.22 EROSON computes erosion of grassed dike 
C     and resulting bottom elevation ZB(J,L) where this subroutine
C     determines time step DELT and whether the end of given ITIME
C     is reached 
      IF(IPROFL.EQ.2) CALL EROSON(ITIME,L,IEND)
C     
C     If IOVER=1, store time series of wave overtopping rate and 
C     sediment transport rates at landward end node JMAX
C     QO(L)=sum of wave overtopping, overflow and seepage rates
C     
        IF(IOVER.EQ.1) THEN
          IF(TIME.EQ.TIMEBC(ITIME)) THEN
            TSQO(L) = 0.D0
            TSQBX(L) = 0.D0
            TSQSX(L) = 0.D0
          ENDIF
          IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
            TSQO(L)=TSQO(L)+DELT*QM
          ELSE
            TSQO(L)=TSQO(L)+DELT*QO(L)
          ENDIF
          TSQBX(L)=TSQBX(L)+DELT*QBX(JMAX(L))
          TSQSX(L)=TSQSX(L)+DELT*QSX(JMAX(L))
        ENDIF
C     
C     Subr.8 OUTPUT stores computed results when IEND=1
C     Put "c" below if no output when TIME = 0
C       IF(TIME.EQ.0.D0) CALL OUTPUT(ITIME,L,ITEQO,ICONV)
        IF(IPROFL.GE.1) TIME=TIME+DELT
        IF(IEND.EQ.1) THEN
          CALL OUTPUT(ITIME,L,ITEQO,ICONV)
        ENDIF
        IF(IPROFL.EQ.0) TIME=TIMEBC(ITIME+1)
C     
C     Compute the BSLOPE at time = (TIME+DELT)
        IF(IPROFL.GE.1) THEN
          DO 501 J=1,JMAX(L)
            IF(J.EQ.1) THEN
              BSLOPE(1,L)=(ZB(2,L)-ZB(1,L))/DX
            ELSE
              IF(J.EQ.JMAX(L)) THEN
                BSLOPE(JMAX(L),L)=(ZB(JMAX(L),L)-ZB(JMAX(L)-1,L))/DX
              ELSE
                BSLOPE(J,L)=(ZB(J+1,L)-ZB(J-1,L))/DX2
              ENDIF
            ENDIF
 501      CONTINUE
C     Compute new bottom RCREST if IPOND=0
          IF(IPOND.EQ.0) THEN
            RCREST(L) = ZB(1,L)
            DO 502 J=2,JMAX(L)
              DUM = ZB(J,L) - RCREST(L)
              IF(DUM.GE.0.D0) THEN
                RCREST(L) = ZB(J,L)
                JCREST(L) = J
              ENDIF
 502        CONTINUE
          ENDIF
        ENDIF
C     
C     If IEND=0, go to 888 for the next landward marching computation
        IF(IEND.EQ.0) GOTO 888
C     
C     If IEND=1 and L is less than ILINE, reset time for next cross-shore line
        IF(L.LT.ILINE) THEN
          TIME=TIMEBC(ITIME)
        ENDIF
C     
C     IEND=1 and specify the seaward input for the next value of 
C     ITIME if ITIME is less than NTIME and L=ILINE.
        IF(ITIME.LT.NTIME.AND.L.EQ.ILINE) THEN
          ITIME1=ITIME+1
          TP=TPBC(ITIME1)
          HRMS(1)=HRMSBC(ITIME1)
C     NPT=integer used in Subr.14 SMOOTH
C     NPE=integer used in Subr.15 EXTRAPO
          NPT=1+NINT(HRMS(1)/DX)
          NPE=1+NINT(HRMS(1)/DX2)
          IF(IPROFL.EQ.1.AND.IPERM.EQ.1) THEN
C           NPT=NPT+NINT(HRMS(1)/DXD2)
            IF(HRMS(1).LT.0.05) NPT=NPT+NINT(HRMS(1)/DX)
          ENDIF
C         IF(IVWALL(L).EQ.2) NPT=NPT+NINT(HRMS(1)/DX)
          WSETUP(1)=WSETBC(ITIME1)
          ANGLE=WANGBC(ITIME1)
C     No wave and current interaction for IANGLE=1
          WKPO=TWOPI*TWOPI/(GRAV*TP*TP)
C     where WKPO is the deep water wave number
        ENDIF
C     
 998  CONTINUE
C     **************** END OF ILINE COMPUTATIAON ***************************
C     
 999  CONTINUE
C     
C     **************** END OF TIME MARCHING COMPUTATION ********************
C     
      END
C     -00-------------------  END OF MAIN PROGRAM  ----------------------
C     #01####################  SUBROUTINE OPENER  ########################
C     
C     This subroutine opens all input and output files
C     
      SUBROUTINE OPENER(BASENAME)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*90 BASENAME
C     
      OPEN(UNIT=11,FILE=trim(basename)//'infile',
     +      STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(UNIT=20,FILE=trim(basename)//'ODOC',  
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=21,FILE=trim(basename)//'OBPROF',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=22,FILE=trim(basename)//'OSETUP',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=23,FILE=trim(basename)//'OPARAM',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=24,FILE=trim(basename)//'OXMOME',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=25,FILE=trim(basename)//'OYMOME',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=26,FILE=trim(basename)//'OENERG',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=27,FILE=trim(basename)//'OXVELO',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=28,FILE=trim(basename)//'OYVELO',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=29,FILE=trim(basename)//'OROLLE',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=30,FILE=trim(basename)//'OBSUSL',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=31,FILE=trim(basename)//'OPORUS',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=32,FILE=trim(basename)//'OCROSS',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=33,FILE=trim(basename)//'OLONGS',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=34,FILE=trim(basename)//'OSWASH',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=35,FILE=trim(basename)//'OSWASE',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=36,FILE=trim(basename)//'OTIMSE',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=37,FILE=trim(basename)//'OCRVOL',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=38,FILE=trim(basename)//'OLOVOL',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=39,FILE=trim(basename)//'ODIKER',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      OPEN(UNIT=40,FILE=trim(basename)//'OMESSG',
     +      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
C     
      RETURN
      END
C     
C     -01-----------------  END OF SUBROUTINE OPENER  --------------------
C     #02#####################  SUBROUTINE INPUT  #######################
C     
C     This subroutine reads data from primary input data file 
C     
      SUBROUTINE INPUT(VER)
C     
      IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000, NB=30000, NL=100)
      CHARACTER COMMEN*70, VER*70 !bdj
      DIMENSION TWAVE(NB),TPIN(NB),HRMSIN(NB),WANGIN(NB),TSURG(NB),
     +     SWLIN(NB),TWIND(NB),WIND10(NB),WINDAN(NB),TSLAND(NB),
     +     SLANIN(NB),TTIDE(NB),DEDYIN(NB),DSDTIN(NB)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +  ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +  IVEG,ICLAY
      COMMON /PERIOD/ TP, WKPO,ANGLE,WT(NN)
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     +  SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /CONSTA/ GRAV, SQR2, SRQ8,PI,TWOPI,SQRG1,SQRG2
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /BINPUT/ XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),
     +  YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +  SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /WBREAK/ GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN)
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +  EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +  WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +  UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /WIND/   W10(NB),WANGLE(NB),WINDCD(NB),TWXSTA(NB),
     +  TWYSTA(NB)
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /COMPAR/ HWDMIN,NPT,NPE
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
      COMMON /TIDALC/ DETADY(NB),DSWLDT(NB)
      COMMON /VEGETA/VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
      COMMON /DIKERO/EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN),
     +  GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL),GRSR(NN,NL),
     +  GRSRD(NN,NL),DEEB,DEEF
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
      COMMON /SOCLAY/EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL),
     + FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL)
C     
C     Gravitational acceleration GRAV = 9.81 (m/s/s)
      GRAV=9.81D0
C     
C.............. COMMENT LINES
C     Write the first line of this CSHORE program on ODOC output file
      WRITE (20,*) VER !bdj
      WRITE (*,*) VER !bdj
C     NLINES = number of comment lines preceding input data
C     READ (11,1110) NLINES
      READ (11,*) DUM !bdj
      NLINES=NINT(DUM) !bdj
C          
      DO 110 I = 1,NLINES
C       READ  (11,1120) (COMMEN(J),J=1,14)
C       WRITE (20,1120) (COMMEN(J),J=1,14)
C       WRITE (*,*) (COMMEN(J),J=1,14)
        READ (11,'(A70)') COMMEN !bdj
        WRITE (20,*) COMMEN !bdj
        WRITE (*,*) COMMEN !bdj
 110  CONTINUE
C     
C.............. INPUT COMPUTATION OPTIONS
C     ILINE = number of cross-shore lines
C     IQYDY=0 to neglect alongshore gradient of longshore sediment transport
C     IQYDY=1 to include alongshore gradient for profile evolution
C
C     READ(11,1110) ILINE
      READ(11,*) DUM      !bdj
      ILINE = NINT(DUM)   !bdj
      IQYDY=0
C     IF(ILINE.GT.2) READ(11,1110) IQYDY
      IF(ILINE.GT.2) THEN 
        READ(11,*) DUM  !bdj
        IQYDY=NINT(DUM) !bdj
      ENDIF
      IF(ILINE.GT.NL) THEN
        WRITE(*,*) '**Error** ILINE is lager than NL'
        STOP
      ENDIF
C     
C     IPROFL=0 for fixed bottom profile(no input for ISEDAV=0)
C     IPROFL=1 for profile evolution computation(input ISEDAV)
C     IPROFL=2 for dike erosion computation (ISEDAV=0)
C     ISEDAV=0 for unlimited bottom sediment
C     ISEDAV=1 for sediment availability limited by hard bottom
C     ISEDAV=2 for stone availability limited by wire mesh
C     For IPROFL=0 and 2, IQYDY=0 is imposed
C     
C     READ(11,1110) IPROFL
      READ(11,*) DUM      !bdj
      IPROFL = NINT(DUM)  !bdj
      IF(IPROFL.EQ.0.OR.IPROFL.EQ.2) IQYDY=0
      ISEDAV=0
C     IF(IPROFL.EQ.1) READ(11,1110) ISEDAV
      IF(IPROFL.EQ.1) THEN
        READ(11,*) DUM     !bdj
        ISEDAV = NINT(DUM) !bdj
      ENDIF
C     
C.............. INPUT COMPUTATION OPTION
C     IPERM=0 for impermeable bottom
C     IPERM=1 for permeable bottom of stone structure
C     
C     READ(11,1110) IPERM
      READ(11,*) DUM    !bdj
      IPERM = NINT(DUM) !bdj
      IF(ISEDAV.EQ.2.AND.IPERM.EQ.0) THEN
        WRITE(*,*) '**Error** For ISEDAV=2, IPERM=1 is required'
        STOP
      ENDIF
C     
C...........INPUT COMPUTATION OPTIONS
C     IOVER=0 for no wave overtopping and overflow(no additional input)
C     IOVER=1 for wave overtopping and overflow on crest(input IWTRAN)
C     where for dike erosion (IPROFL=2), IOVER=1 to predict erosion of
C     landward dike slope
C     IWTRAN=0 for no standing water landward of crest above seaward SWL but
C     lanward SWL=seaward SWL if crest is submerged and no overflow (IOFLOW=0)
C     IWTRAN=1 for wave transmission into landward water level (input) and
C     if crest is below seaward SWL, wave overtopping and overflow on crest
C     (IOFLOW=1)
C     IPOND=0 for no ponding water landward of JSWL node
C     IPOND=1 for ponding water in runnel landward of ridge for IPERM=0
C     INFILT=0 for no infiltration landward of sand dune
C     INFILT=1 for infiltration if IOVER=1, IPERM=0, IPROFL=1
C     
C     READ(11,1110) IOVER
      READ(11, *) DUM   !bdj
      IOVER = NINT(DUM) !bdj
      IF(IPROFL.EQ.2) IOVER=1
      IWTRAN=0
      IOFLOW=0
      IPOND=0
      INFILT=0
C     IF(IOVER.EQ.1) READ(11,1110) IWTRAN
      IF(IOVER.EQ.1) READ(11,*)DUM;IWTRAN=NINT(DUM) !bdj
C     IF(IOVER.EQ.1.AND.IWTRAN.EQ.0) READ(11,1110) IPOND
      IF(IOVER.EQ.1.AND.IWTRAN.EQ.0) THEN 
        READ(11,*)DUM   !bdj
        IPOND=NINT(DUM) !bdj
      ENDIF
      IF(IOVER.EQ.1.AND.IPERM.EQ.0) THEN
C       IF(IPROFL.EQ.1) READ(11,1110) INFILT
        IF(IPROFL.EQ.1) READ(11,*) DUM;INFILT =NINT(DUM) !bdj
      ENDIF
C     
C........... INPUT COMPUTATION OPTION
C     IWCINT=0 for no wave and current interaction
C     IWCINT=1 for wave and current interaction in frequency 
C     dispersion, momentum and wave action equations 
C     
C     READ(11,1110) IWCINT
      READ(11,*) DUM     !bdj
      IWCINT = NINT(DUM) !bdj
C     
C........... INPUT COMPUTATION OPTION
C     IROLL=0 for no roller
C     IROLL=1 for roller effects in governing equations
C     
C     READ(11,1110) IROLL
      READ(11,*) DUM    !bdj
      IROLL = NINT(DUM) !bdj
      IF(IROLL.EQ.1) RBZERO=0.1D0
C     
C........... INPUT COMPUTATION OPTION
C     IWIND=0 for no wind effect
C     IWIND=1 for wind shear stresses on momentum equations
C     
C     READ(11,1110) IWIND
      READ(11,*) DUM    !bdj
      IWIND = NINT(DUM) !bdj
C     
C........... INPUT COMPUTATION OPTION
C     ITIDE=0 for no tidal effect on currents
C     ITIDE=1 for longshore and cross-shore tidal currents
C
C     READ(11,1110) ITIDE
      READ(11,*) DUM   !bdj 
      ITIDE = NINT(DUM)!bdj
C     
C........... INPUT COMPUTATION OPTION
C     IVEG=0 for no vegetation or vegetation represented by increased
C     bottom friction factor FBINP
C     IVEG=1 for vegataion whose density, width, height and root depth 
C     are specified as input. The height and root depth vary with the 
C     bottom elevation change
C     IVEG=2 for vegatation whose constant density, width and height 
C     are specified as input
C     
C     READ(11,1110) IVEG
      READ(11,*) DUM 
      IVEG = NINT(DUM)
C      
C........... INPUT COMPUTATION OPTION 
C     ISTSAN=0 except for fixed stone structure on sand bottom
C     ISTSAN=1 for stone structure (IPERM=1) on deforming bottom
C     (IPROFL=1) of unlimited sand (ISEDAV=0)
C     CPSTON=empirical parameter for sand transport reduction on porous stone structure
      ISTSAN=0
      IF(IPROFL.EQ.1.AND.IPERM.EQ.1) THEN
          IF(ISEDAV.EQ.0) THEN
              ISTSAN=1
              CPSTON=1.0D0
          ENDIF
      ENDIF
C
C........... INPUT COMPUTATION OPTION
C     ICLAY=0 except for eroding sand layer on erodible clay
C     ICLAY=1 for sand layer (ISEDAV=1 and IPERM=0) above clay
C     bottom (eroded by wave action) with no vegetation (IVEG=0)    
      ICLAY=0
      IF(ISEDAV.EQ.1.AND.IPERM.EQ.0) THEN
          IF(IVEG.EQ.0) THEN
              READ(11,*)DUM
              ICLAY=NINT(DUM)
          ENDIF
      ENDIF
C     
C........... COMPUTATIONAL INPUT DATA
C     DX=nodal spacing for input bottom geometry
C     
C     READ(11,1130) DX
      READ(11,*) DX !bdj
C     
C........... BREAKER RATIO PARAMETER GAMMA=0.5-1.0
C     READ(11,1130) GAMMA
      READ(11,*) GAMMA !bdj
C     
C........... SEDIMENT CHARACTERISTICS IF IPROFL=1
C     WF    = sediment fall velocity (m/s)
C     SG    = sediment specific gravity
C     D50   = median sediment diameter (mm) 
C     converted to (m) below	
C     EFFB  = suspension efficiency due to breaking, eB
C     EFFF  = suspension efficiency due to friction, ef
C     SLP   = suspended load parameter
C     SLPOT = suspended load parameter due to wave overtopping if IOVER=1
C     SPORO = sediment porosity (SPORO=0.4 for sand but input SNP used for IPERM=1)
C     SHIELD= critical Shields parameter used if D50 is less than CSEDIA
C     CSTABN= critical stability number (0.6 to 1.1) used if IPERM=1
C     CSEDIA= critical sediment diameter to separate sand and stone
C     TANPHI= tangent (sediment friction angle)
C     BLP   = bedload parameter
C     BLD	= BLP/GRAV/(SG-1) used for bedload transport rate
C     BEDLM = parameter m for bedload reduction factor BRF for hard
C     bottom used for ISEDAV=1	
C     Following default values are specified to reduce input error
      IF(IPROFL.EQ.1) THEN
        SPORO = 0.4D0
        SHIELD = 0.05D0
C mg     EFFF = 0.01D0
C        EFFB = 0.002D0 to 0.01D0
C        SLP = 0.1D0 to 0.4D0
C        SLPOT = 0.1D0 to 3.6D0
c        TANPHI = 0.63D0 for sand
c        BLP = 0.001D0 to 0.004D0
        IF(ISEDAV.GE.1) BEDLM=1.0D0
C     
C       READ (11,1150) D50,WF,SG
        READ (11,*) D50,WF,SG !bdj
C mg - added read for EFFF
C       IF(IOVER.EQ.0) READ (11,1150) EFFB,EFFF,SLP
        IF(IOVER.EQ.0) READ (11,*) EFFB,EFFF,SLP !bdj
C       IF(IOVER.EQ.1) READ (11,1150) EFFB,EFFF,SLP,SLPOT
        IF(IOVER.EQ.1) READ (11,*) EFFB,EFFF,SLP,SLPOT !bdj
C mg
C       READ (11,1150) TANPHI,BLP
        READ (11,*) TANPHI,BLP !bdj
C mg
        IF(EFFF.LT.EFFB) THEN
          WRITE(*,*) ' ** Error ** The suspension efficiency parameter'
          WRITE(*,*) ' due to bottom friction must be greater than or '
          WRITE(*,*) ' equal to the suspension efficiency parameter due'
          WRITE(*,*) ' to wave breaking.' 
          STOP
        ENDIF
C mg
        D50    = D50*1.D-3
        SGM1   = SG - 1.D0
        SPORO1 = 1.D0 - SPORO
        WFSGM1 = WF*SGM1
        GSGM1  = GRAV*SGM1
        IF(IPERM.EQ.0.OR.ISTSAN.EQ.1) THEN
          GSD50S = GSGM1*D50*SHIELD
          CSEDIA=2.D0*D50
        ENDIF
        BLD = BLP/GSGM1
      ENDIF
C     
C.....RUNUP WIRE HEIGHT RWH (in meters) IF IOVER=1
C     IF(IOVER.EQ.1) READ (11,1130) RWH
      IF(IOVER.EQ.1) READ (11,*) RWH !bdj
C     
C.....STONE OR GRAVEL CHARACTERISTICS IF IPERM=1
C     SNP = Stone/gravel porosity in porous layer (SNP can be different from sand
C     porosity=0.4 for ISTSAN=1)
C     SDP = Nominal stone/gravel diameter (m)
C     CSTABN = Critical stability number (0.6 to 1.1) for stone
C     
      IF(IPERM.EQ.1) THEN
C       READ(11,1150) SNP,SDP,CSTABN
        READ(11,*) SNP,SDP,CSTABN
        IF(IPROFL.EQ.1.AND.ISTSAN.EQ.0) THEN
          GSD50S = DSQRT(GSGM1*D50*CSTABN)
          CSEDIA=0.5D0*D50
          SPORO=SNP
          SPORO1=1.D0-SPORO
        ENDIF
      ENDIF
C     
C.....DIKE EROSION EFFICIENCIES IF IPROFL=2
C     DEEB=eB due to breaking waves
C     DEEF=ef due to bottom friction
      IF(IPROFL.EQ.2.OR.ICLAY.EQ.1) THEN
        READ(11,*) DEEB,DEEF
      ENDIF
C     
C     HWDMIN=minimum water depth (m)
C     used in the wet and dry zone in Subr.16 WETDRY
C     D50=median sediment diameter (m)
      IF(IOVER.EQ.1) THEN
        HWDMIN=1.D-6
        IF(IPROFL.EQ.1.AND.IPERM.EQ.0) THEN
          HWDMIN=D50
C         HWDMIN=0.2D0*HWDMIN
          IF(DX.GE.0.05D0) HWDMIN=1.D-3
        ELSE
C         HWDMIN=1.D-5
          IF(IPROFL.GE.1) HWDMIN=1.D-4
        ENDIF
      ENDIF
C     
C     ......... INPUT WAVE AND WATER LEVEL
C     NWAVE      = number of waves at x=0 starting from time=0
C     NSURG      = number of water levels at x=0 from time=0
C     NTIME      = number of waves and water levels at x=0
C     During time= TIMEBC(i) to time=TIMEBC(i+1) if NWAVE=NSURG
C     TIMEBC(i)  = time in seconds at the beginning of the
C     specified wave and water level
C     TPBC(i)    = spectral peak or wave period in seconds
C     HRMSBC(i)  = root mean square wave height in meters
C     WSETBC(i)  = wave setup in meters
C     SWLBC(i)   = still water level in meters above the
C     datum used for the input bottom profile
C     WANGBC(i)  = incident wave angle in degrees from shorenormal if ILINE=1
C     and from reference direcction(e.g.,North) if ILINE=2 or larger
C     For IPROFL = 0, use input TIMEBC(i+1)=1.0,2.0,... to identify
C     each combination of waves and still water level
C     
C mg     
C mg........INPUT WAVE and WATER LEVEL OPTION
C mg  ILAB=0 for field data set - waves and water levels read separately
C mg  ILAB=1 for laboratory data set - waves and water levels read together 
C mg  
C     READ(11,1110) ILAB
      READ(11,*) DUM !bdj
      ILAB=NINT(DUM) !bdj
C mg
C     READ(11,1110) NWAVE
      READ(11,*) DUM  !bdj
      NWAVE=NINT(DUM) !bdj
C     READ(11,1110) NSURG
      READ(11,*) DUM  !bdj
      NSURG=NINT(DUM) !bdj
C mg
      IF(ILAB.EQ.1) THEN
C mg
        NTIME=NWAVE
        TIMEBC(1) = 0.D0
        DO 120 I = 1,NTIME
C         READ(11,1160) TIMEBC(I+1),TPBC(I),HRMSBC(I),WSETBC(I),
C    +       SWLBC(I),WANGBC(I)
          READ(11,*) TIMEBC(I+1),TPBC(I),HRMSBC(I),WSETBC(I), !bdj
     +       SWLBC(I),WANGBC(I)                             !bdj
C         IF(WANGBC(I).LT.-80.D0.OR.WANGBC(I).GT.80.D0) THEN
C           WRITE (*,2800) WANGBC(I)
C 2800      FORMAT(/'Incident Wave Angle=',D11.4,
C    +           'but Angle must be in the range of -80 
C    +to 80 in degree')
C           STOP
C         ENDIF
 120    CONTINUE
      ELSE
C     
C     For field data,wave conditions(NWAVE+1) and water level(NSURG+1)
C     at x=0 vary continously in time starting from time=0 unlike
C     step changes assumed for laboratory data, defined by ILAB=1.
C     Choose number of step changes to approximate time series
        NTIME=MAX0(NWAVE,NSURG)
        NWAVE1=NWAVE+1
        DO 130 I=1,NWAVE1
C         READ(11,1170) TWAVE(I),TPIN(I),HRMSIN(I),WANGIN(I)
          READ(11,*) TWAVE(I),TPIN(I),HRMSIN(I),WANGIN(I) !bdj
          IF(NWAVE.EQ.NTIME) TIMEBC(I)=TWAVE(I)
 130    CONTINUE
 1170   FORMAT(D11.1,3D11.4)
        NSURG1=NSURG+1
        DO 131 I=1,NSURG1
C         READ(11,1180) TSURG(I),SWLIN(I)
          READ(11,*) TSURG(I),SWLIN(I) !bdj
          IF(NSURG.EQ.NTIME) TIMEBC(I)=TSURG(I)
 131    CONTINUE
 1180   FORMAT(D11.1,D11.4)
        IF(TWAVE(1).NE.0.D0.OR.TSURG(1).NE.0.D0) THEN
          WRITE(*,2801)
          STOP
 2801     FORMAT(/'Data input is stopped because the start
     +      time for offshore wave conditions and water
     +      level is NOT ZERO'/)
        ENDIF
        IF(TWAVE(NWAVE1).NE.TSURG(NSURG1)) THEN
          WRITE(*,2802)
          STOP
 2802     FORMAT(/'Data input stopped because the durations
     +      of offshore wave conditions and water level
     +      are NOT SAME'/)
        ENDIF
C     
C     Subr.19 TSINTP interpolates input time series at
C     specified time TIMEBC(I) with I=1,2,...,(NTIME+1) and
C     generates stepped time series corresponding to input format
C     for the case of NWAVE=NSURG
        CALL TSINTP(NWAVE,TWAVE,TPIN,NTIME,TIMEBC,TPBC)
        CALL TSINTP(NWAVE,TWAVE,HRMSIN,NTIME,TIMEBC,HRMSBC)
        CALL TSINTP(NWAVE,TWAVE,WANGIN,NTIME,TIMEBC,WANGBC)
        CALL TSINTP(NSURG,TSURG,SWLIN,NTIME,TIMEBC,SWLBC)
C     Wave setup at x=0 is assumed to be zero
        WSETBC(1:NTIME)=0.D0
C       DO 132 I=1,NTIME
C         IF(WANGBC(I).LT.-80.D0.OR.WANGBC(I).GT.80.D0) THEN
C           WRITE(*,2800) WANGBC(I)
C           STOP
C         ENDIF
C 132   CONTINUE
C     
      ENDIF
C     End of field data input for (ILAB=0)
C     
C     Prepare for ITIME=1 computation
C     IF IPOND=1, ponded water level ZW=SWL at time=0
      TP = TPBC(1)
      HRMS(1) = HRMSBC(1)
      WSETUP(1) = WSETBC(1)
      ANGLE= WANGBC(1)
      IF(IPOND.EQ.1) ZW=SWLBC(1)
C     
C     ......... BOTTOM GEOMETRY and POROUS LAYER BOTTOM if IPERM=1
C     The bottom geometry is divided into segments of
C     different inclination and roughness starting from
C     seaward boundary for cross-shore line L.
C     YLINE(L)  = alongshore coordinate for line L=1,2,...,ILINE 
C     AGLINE(L) = angle of line L from reference direction(e.g.,North) but
C     YLINE(1)=0.0 and AGLINE(1)=0.0 if ILINE=1
C     NBINP(L)  = number of input points of bottom elevation ZB(X)
C     XBINP(J,L)= horizontal distance to input bottom point (J)
C     in meters where XBINP(1,L) = 0 at the seaward boundary
C     ZBINP(J,L)= dimensional vertical coordinate (+ above datum)
C     of input bottom point (J) in meters along cross-shore line L
C     FBINP(J,L)= bottom friction factor for segment between
C     points (J) and (J+1) along cross-shore line L where if IVEG=1,
C     FBINP needs to be positive in vegetated zone
C     WMINP(J,L)=1.0 for wire mesh segment between points (J) and (J+1)
C     along cross-shore line L for ISEDAV=2 where WMINP(J,L)=0.0 for no
C     wire mesh segment
C     NPINP(L)  = number of input points of impermeable hard or clay bottom
C     ZP(X) along cross-shore line L only if IPERM=1 or ISEDAV=1 but
C     for ISTSAN=1,ZP(X) is sand bottom elevation beneath stone structure
C     XPINP(J,L)= horizontal distance of input point J from x=0
C     ZPINP(J,L)= dimensional vertical coordinate in meters of
C     porous layer bottom or hard or clay bottom at point (J) with ZPINP(J)
C     equal to or less than ZBINP(J,L) where ZPINP(1,L)=ZBINP(1,L) imposed
C     If ICLAY=1, clay resistance and sand fraction in clay are input
C     RCINP(J,L)= clay resistance parameter of order of 10 m*m/s/s
C     FCINP(J,L)= sand volume per unit clay volume in range of 0.0 to (1-SPORO)=SPORO1
C     IF ISEDAV = 1, an almost vertical impermeable wall can be specified
C     using two points (NPINP-1) and NPINP where
C     IVWALL(L) = 0 for no vertical wall along cross-shore line L
C     IVWALL(L) = 1 for vertical wall with sediment in front
C     IVWALL(L) = 2 for vertical wall exposed to wave action
      DO 160 L=1, ILINE
      YLINE(L)=0.D0
      AGLINE(L)=0.D0
C     IF(ILINE.GT.1) READ(11,1150) YLINE(L),AGLINE(L)
      IF(ILINE.GT.1) READ(11,*) YLINE(L),AGLINE(L) !bdj
C     READ(11,1110) NBINP(L)
      READ(11,*) DUM       !bdj
      NBINP(L) = NINT(DUM) !bdj
C     IF(IPERM.EQ.1.OR.ISEDAV.GE.1) READ(11,1110) NPINP(L)
      IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
        READ(11,*) DUM       !bdj
        NPINP(L) = NINT(DUM) !bdj
      ENDIF
      IF(NBINP(L).GT.NN) THEN
        WRITE(*,2900) L, NBINP(L), NN
        STOP
      ENDIF
 2900 FORMAT(/'Number of Input Bottom Nodes NBINP(',I3,') = ',I8,'
     +   ;NN = ',I8/'Increase PARAMETER NN.')
C     
C     Point J = 1 has no corresponding friction factor.
C     READ (11,1150) XBINP(1,L), ZBINP(1,L)
      READ (11,*) XBINP(1,L), ZBINP(1,L) !bdj
      XBINP(1,L) = 0.D0
      DO 140 J = 2,NBINP(L)
C       READ(11,1150) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L)
        IF(ISEDAV.LE.1) THEN
          READ(11,*) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L) !bdj
        ELSE
          READ(11,*) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L), WMINP(J-1,L)
        ENDIF
C     IF IANGLE = 1, the bottom friction factor must be positive
        IF(ANGLE.NE.AGLINE(L).OR.IVEG.GE.1) THEN
          IF(FBINP(J-1,L).LE.0.D0) THEN
            WRITE(*,2901) FBINP(J-1,L), (J-1), L
            STOP
          ENDIF
        ENDIF
C     Avoid perfect horizontal bottom for possible numerical difficulty
        IF(ZBINP(J,L).EQ.ZBINP(J-1,L)) ZBINP(J-1,L)=ZBINP(J-1,L)-1.D-4
 140  CONTINUE
      DUM=XBINP(NBINP(L),L)/DX
      IDUM=NINT(DUM)
      DUM=DUM-DBLE(IDUM)
      IF(DUM.LT.1.D-5) XBINP(NBINP(L),L)=XBINP(NBINP(L),L)+1.D-4
 2901 FORMAT(/'Bottom Friction Factor FBINP(J-1,L)=', D11.4,
     +   'for (J-1) =',I4,'and L=',I3/'For obliquely incident
     +    wave or vegetated zone, FBINP must be positive')
      IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
        XPINP(1,L)=0.D0
        ZPINP(1,L)=ZBINP(1,L)
        DO 150 J=2,NPINP(L)
          IF(ICLAY.EQ.0)THEN
C             READ(11,1150) XPINP(J,L),ZPINP(J,L)
              READ(11,*) XPINP(J,L),ZPINP(J,L) !bdj
          ELSE
              READ(11,*) XPINP(J,L),ZPINP(J,L),RCINP(J,L),FCINP(J,L)
          ENDIF
 150    CONTINUE
        IF(XPINP(NPINP(L),L).LT.XBINP(NBINP(L),L)) XPINP(NPINP(L),L)
     +    =XBINP(NBINP(L),L)
      ENDIF
      IF(L.GT.1) DYLINE(L-1)=YLINE(L)-YLINE(L-1)
C     
C.....VEGETATION CHARACTERISTICS IF IVEG=1 or 2
C     VEGCD    = Vegetation drag coefficient of order of unity
C     VEGN(J,L)= number of vegetation (1/m/m) per unit horizontal area
C     for segment J(J=1,2,...,NBINP(L)-1) along cross-shore line L where
c     VEGN(J,L)=0.0 if no vegetation
C     VEGB(J,L)= width(m) of each vegetation stand where VEGB(J,L)=0.0
C     if no vegetation
C     VEGD(J,L)= height(m) of each vegetation stand above sand where 
C     VEGD(J,L)=0.0 if no vegetation
C     VEGRD(J,L)=root depth (m) below sand for no vegetation uprooting
C     where uprooting occurs when erosion reaches this depth
C     (input only for IVEG=1) where VEGRD(J,L)=0.0 if no vegetation
      IF(IVEG.GE.1) THEN
C     READ(11,1130) VEGCD
      READ(11,*) VEGCD
      JDUM=NBINP(L)-1
      DO 170 J=1,JDUM
        IF(IVEG.EQ.1) THEN
          READ(11,*)VEGN(J,L),VEGB(J,L),VEGD(J,L),VEGRD(J,L)
        ELSE
C         READ(11,1150) VEGN(J,L),VEGB(J,L),VEGD(J,L)
          READ(11,*) VEGN(J,L),VEGB(J,L),VEGD(J,L)
        ENDIF
        VEGINP(J,L)=VEGCD*VEGN(J,L)*VEGB(J,L)/FBINP(J,L)
        IF(VEGINP(J,L).LT.0.D0) THEN
          WRITE(*,2902) VEGINP(J,L),J,L
          STOP
        ENDIF
 170    CONTINUE
      ENDIF
 2902 FORMAT(/'Vegetation Input Characteristic VEGINP(J,L)
     +   =',F11.4,'for Segment J=',I4,'and Line L=',I3/'
     +   Vegetation CD,N and B must be positive or zero')
C
C.....DIKE GRASS AND SOIL CHARACTERISTICS IF IPROFL=2
C     GDINP(J,L)=thickness (m) of grass cover for segment J along cross-
C     shore line L where GDINP(J,L)=0.0 if no grass cover
C     GRINP(J,L)=grass surface resistance parameter (m*m/s/s)  
C     GRDINP(J,L)=resistance parameter (m*m/s/s) below grass cover
C     where grass resistance is assumed to decrease linearly downward
C     in grass cover and be (+) constant below grass cover
      IF(IPROFL.EQ.2) THEN
        DO 180 J=1,NBINP(L)-1
          READ(11,*) GDINP(J,L),GRINP(J,L),GRDINP(J,L)
          IF(GDINP(J,L).LT.0.D0) GDINP(J,L)=0.D0
          IF(GRDINP(J,L).LE.0.01D0) GRDINP(J,L)=0.01D0
          IF(GDINP(J,L).EQ.0.D0) GRINP(J,L)=GRDINP(J,L)
          IF(GDINP(J,L).GT.0.D0) THEN
            IF(GRINP(J,L).LE.GRDINP(J,L)) GRINP(J,L)=GRDINP(J,L)+0.01D0
          ENDIF
 180    CONTINUE
      ENDIF
C
 160  CONTINUE
C     End of line L=1,2,...,ILINE
C     
C.....WIND SPEED AND DIRECTION IF IWIND=1
C     During time = TIMEBC(i) to time=TIMEBC(i+1)
C     W10(i) = wind speed (m/s) at 10m elevation above mean sea level
C     WANGLE(i)= wind direction in degrees at 10 m
C     from a cross-shore line(not adjusted for a curved beach)
C     WINDCD(i)= wind drag coefficient based on Large and Pond (1981)
C     TWXSTA(i)= cross-shore wind shear stress/specific water weight
C     TWYSTA(i)= longshore wind shear stress/specific water weight
C     RATIO = specific water weight/specific air weight
C     
C     Wind data time series is read in the same way as
C     field data of waves and water level
      IF(IWIND.EQ.1) THEN
C       READ(11,1110) NWIND
        READ(11,*) DUM    !bdj
        NWIND = NINT(DUM) !bdj
        NWIND1=NWIND+1
        DO 190 I=1,NWIND1
C         READ(11,1190) TWIND(I),WIND10(I),WINDAN(I)
          READ(11,*) TWIND(I),WIND10(I),WINDAN(I) !bdj
 190    CONTINUE
 1190   FORMAT(D11.1,2D11.4)
        IF(TWIND(1).NE.0.D0) THEN
          WRITE(*,2905)
          STOP
 2905     FORMAT(/'Data input is stopped because the start time of 
     +       wind data is NOT ZERO'/)
        ENDIF
        IF(TWIND(NWIND1).NE.TIMEBC(NTIME+1)) THEN
          WRITE(*,2906)
          STOP
 2906     FORMAT(/'Data input is stopped because the end time of
     +       wind data is NOT SAME as the end time of 
     +       wave and water level data'/)
        ENDIF
        CALL TSINTP(NWIND,TWIND,WIND10,NTIME,TIMEBC,W10)
        CALL TSINTP(NWIND,TWIND,WINDAN,NTIME,TIMEBC,WANGLE)
        RATIO = 837.D0
        CONVRT = 3.14159D0/180.D0
        DO 200 I=1,NTIME
          IF(W10(I).GT.25.D0) WRITE(*,2910)
          IF(W10(I).LT.11.D0) THEN
            WINDCD(I) = 1.2D-3
          ELSE
            WINDCD(I)=0.49D-3 + 0.065D-3*W10(I)
          ENDIF
          DUM = (WINDCD(I)/RATIO/GRAV)*W10(I)*W10(I)
          ANG = CONVRT*WANGLE(I)
          TWXSTA(I)=DUM*DCOS(ANG)
          TWYSTA(I)=DUM*DSIN(ANG)
 200    CONTINUE
      ELSE
        DO 201 I=1,NTIME
          TWXSTA(I) = 0.D0
          TWYSTA(I) = 0.D0
 201    CONTINUE
      ENDIF
 2910 FORMAT(/'Wind speed at 10m =',D11.4/
     +   'but wind speed must be less than 25m/s for
     +   Large and Pond(1981)')
C     
C.....LANDWARD STILL WATER LEVEL IF IWTRAN=1
C     During time=TIMEBC(i) to time=TIMEBC(i+1)
C     SWLAND(i)=still water level in meters above datum
C     landward of emerged structure or dune if IOFLOW=0
C     If ISWLSL=0, seaward and landward still water levels
C     are same and SWLAND(i)=SWLBC(i)
C     If ISWLSL=1, read time series of landward still water level
C     SLANIN(I) at time TSLAND(I) with I=1,2,...,(NSLAN+1)
C     If ISWLSL=2, no water landward of structure or dune and overflow
C     (IOFLOW=1) occurs if crest is submerged
      IF(IWTRAN.EQ.1) THEN
C       READ(11,1110) ISWLSL
        READ(11,*) DUM !bdj
        ISWLSL = NINT(DUM) !bdj
        IF(ISWLSL.EQ.0) THEN
          DO 300 I=1,NTIME
            SWLAND(I)=SWLBC(I)
 300      CONTINUE
        ENDIF
        IF(ISWLSL.EQ.1) THEN
C         READ(11,1110) NSLAN
          READ(11,*) DUM !bdj
          NSLAN = NINT(DUM) !bdj
          NSLAN1=NSLAN+1
          DO 301 I=1,NSLAN1
C           READ(11,1180) TSLAND(I),SLANIN(I)
            READ(11,*) TSLAND(I),SLANIN(I) !bdj
 301      CONTINUE
          IF(TSLAND(1).NE.0.D0) THEN
            WRITE(*,2950)
            STOP
 2950       FORMAT(/'Data input is stopped because the start time of
     +landward SWL is NOT ZERO'/)
          ENDIF
          IF(TSLAND(NSLAN1).NE.TIMEBC(NTIME+1)) THEN
            WRITE(*,2951)
            STOP
 2951       FORMAT(/'Data input is stopped because the end time of
     +        landward SWL is NOT SAME as the end time of 
     +        other input time series'/)
          ENDIF
          CALL TSINTP(NSLAN,TSLAND,SLANIN,NTIME,TIMEBC,SWLAND)
        ENDIF
      ENDIF
C
C.....ALONGSHORE WATER LEVEL GRADIENT IF ITIDE=1
C     During time=TIMEBC(i) to time=TIMEBC(i+1)
C     DETADY(i) = alongshore water level gradient for longshore current
C     DSWLDT(i) = rate of input water level change only for ILAB=0
C     
C     Alongshore gradient data time series is read in the same way as 
C     field surge data
      IF(ITIDE.EQ.1) THEN
C       READ(11,1110) NTIDE
        READ(11,*) DUM    !bdj
        NTIDE = NINT(DUM) !bdj
        NTIDE1=NTIDE+1
        DO 400 I=1,NTIDE1
C         READ(11,1195) TTIDE(I),DEDYIN(I)
          READ(11,*) TTIDE(I),DEDYIN(I) !bdj
 400    CONTINUE
 1195   FORMAT(D11.1,D11.7)
        IF(TTIDE(1).NE.0.D0) THEN
          WRITE(*,2961)
          STOP
 2961     FORMAT(/'Data input is stopped because the start time
     +       of tide data is NOT ZERO'/)
        ENDIF
        IF(TTIDE(NTIDE1).NE.TIMEBC(NTIME+1)) THEN
          WRITE(*,2962)
          STOP
 2962     FORMAT(/'Data input is stopped because the end
     +      time of tide data is NOT SAME as the end time of
     +      wave and water level data'/)
        ENDIF
        CALL TSINTP(NTIDE,TTIDE,DEDYIN,NTIME,TIMEBC,DETADY)
      ENDIF
C     
C     If ITIDE=1 and ILAB=0, cross-shore water flux associated
C     with continuous input water level change is accounted
C     for in cross-shore current in wet zone
      IF(ITIDE.EQ.1.AND.ILAB.EQ.0) THEN
        DO 410 I=1,NSURG
          K=I+1
          DSDTIN(I)=(SWLIN(K)-SWLIN(I))/(TSURG(K)-TSURG(I))
  410   CONTINUE
        DSDTIN(NSURG1)=DSDTIN(NSURG)
        CALL TSINTP(NSURG,TSURG,DSDTIN,NTIME,TIMEBC,DSWLDT)
      ENDIF
C     
      CLOSE (11)
C     
 1110 FORMAT (I8)
 1120 FORMAT (14A5)
 1130 FORMAT (D11.4)
 1150 FORMAT (4D11.4)
 1160 FORMAT (D11.1,5D11.4)
C     
      RETURN
      END
C     
C     -02-----------------  END OF SUBROUTINE INPUT  --------------------
C     #03####################  SUBROUTINE BOTTOM  ########################
C     
C     This subroutine calculates the bottom geometry using input
C     DX between two adjacent nodes along ILINE cross-shore lines
C     Smooth input ZB(J,L) to reduce numerical irregularity
C     
      SUBROUTINE BOTTOM
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000, NB=30000, NL=100)
      DIMENSION SLOPE(NN), PSLOPE(NN), ZBRAW(NN), ZPRAW(NN)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     + ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     + IVEG,ICLAY
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     + SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /BINPUT/ XBINP(NN,NL), ZBINP(NN,NL), FBINP(NN,NL),XS(NL),
     + YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     + SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     + EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),
     + RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     + WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     + UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /COMPAR/ HWDMIN,NPT,NPE
      COMMON /VEGETA/ VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
      COMMON /DIKERO/EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN),
     + GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL),GRSR(NN,NL),
     + GRSRD(NN,NL),DEEB,DEEF
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
      COMMON /SOCLAY/EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL),
     + FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL)
C     
C     The structure geometry is divided into segments of different
C     inclination and roughness for each cross-shore line L.
C     NBINP(L) = number of input bottom points
C     For segments starting from the seaward boundary:
C     SLOPE(K)  = slope of segment K(+ upslope, - downslope)
C     FBINP(K,L)  = bottom friction factor
C     XBINP(K,L)  = dimensional horizontal distance from seaward
C     boundary to the seaward end of segment K
C     ZBINP(K,L)  = dimensional vertical coordinate (+ above datum)
C     at the seaward end of segment K
C     PSLOPE(K) = slope of porous layer bottom or hard bottom
C     XPINP(K,L)  = dimensional horizontal distance of porous layer bottom
C     at the seaward end of segment K
C     ZPINP(K,L)  = dimensional vertical coordinate of porous layer bottom
C     at the seaward end of segment K
C     
      DO 100 L = 1,ILINE
      DO 120 K = 1,NBINP(L)-1
        DUM=XBINP(K+1,L)-XBINP(K,L)
        SLOPE(K) = (ZBINP(K+1,L)-ZBINP(K,L))/DUM
 120  CONTINUE
C     No  vertical wall at landward end unless IVWALL=1 or 2
      IVWALL(L)=0 
      IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
        DO 121 K=1,NPINP(L)-1
          DUM=XPINP(K+1,L)-XPINP(K,L)
          PSLOPE(K)=(ZPINP(K+1,L)-ZPINP(K,L))/DUM
 121    CONTINUE
        IF(PSLOPE(NPINP(L)-1).GT.TANPHI) IVWALL(L)=1
      ENDIF
C     
C     ... INITIAL SHORELINE LOCATION AT DATUM Z=0
C     
C     XS(L)= horizontal distance between X=0 and shoreline
      K = 0
 900  CONTINUE
      K=K+1
      IF(K.EQ.NBINP(L)) THEN
        XS(L) = XBINP(NBINP(L),L)
        GOTO 901
      ENDIF
      CROSS = ZBINP(K,L)*ZBINP(K+1,L)
      IF (CROSS.GT.0.D0) GOTO 900
      XS(L)  = XBINP(K+1,L) - ZBINP(K+1,L)/SLOPE(K)
 901  IF(L.EQ.1) THEN
      DXD2 = DX/2.D0
      DX2 = 2.D0*DX
      DXDX = DX*DX
C     
C     NPT= integer used in Subr.14 SMOOTH
C     NPE= integer used in Subr.15 EXTRAPO
C BDJ 2011->2014 on 2014-10-02 
          NPT=1+NINT(maxval(HRMSBC)/DX)
          NPE=1+NINT(maxval(HRMSBC)/DX2)
C      NPT=1+NINT(HRMS(1)/DX)
C      NPE=1+NINT(HRMS(1)/DX2)
C END BDJ 2011->2014 on 2014-10-02 
C     IF(IPROFL.EQ.1.AND.IPERM.EQ.1) NPT=NPT+2*NINT(SDP/DX)
      ENDIF
C     
C     ... CALCULATE BOTTOM GEOMETRY AT EACH NODE 
C     
C     JMAX(L) = landward edge node corresponding to maximum node number
C     XB(J)= horizontal coordinate of node J where XB(1) = 0
C     ZB(J,L)= vertical coordinate of bottom at node J (+ above datum)
C     BSLOPE(J,L) = bottom slope at node J for cross-shore line L
C     SLOPE(K) = tangent of local slope of segment K
C     
      DUM = XBINP(NBINP(L),L)/DX
      JMAX(L)  = NINT(DUM)+1
      DUM=DX*DBLE(JMAX(L)-1)-XBINP(NBINP(L),L)
      IF(DUM.GT.0.D0) JMAX(L)=JMAX(L)-1
      IF(JMAX(L).GT.NN) THEN
        WRITE (*,2910) L,JMAX(L),NN
        STOP
      ENDIF
 2910 FORMAT (/' End Node of Line',I3,':JMAX(L)=',I8,'; NN =',I8/
     +   ' Bottom length is too long.'/
     +   ' Cut it, or change PARAMETER NN.')
C     
C     INTERPOLATION OF BOTTOM POSITION at XB(J)
C     RCREST(L) = crest (highest) elevation above datum Z=0
C     JCREST(L) = nodal location of crest for cross-shore Line L
C     If IPOND=1, JCREST(L)=nodal location of ridge crest computed in Subr.21 PONDED
      IF(L.EQ.1) JDUM=JMAX(L)
      IF(JMAX(L).LT.JDUM) GOTO 130
      JDUM=JMAX(L)
      DO 141 J = 1,JMAX(L)
        XB(J) = DX*DBLE(J-1)
 141  CONTINUE
 130  CONTINUE
      ZBRAW(1) = ZBINP(1,L)
      FB2(1,L)=0.5D0*FBINP(1,L)
      IF(ISEDAV.EQ.2) WMNODE(1,L)=WMINP(1,L)
      IF(IVEG.GE.1) THEN
        VEGFB(1,L)=VEGINP(1,L)
        VEGH(1,L)=VEGD(1,L)
        IF(IVEG.EQ.1) VEGRH(1,L)=VEGRD(1,L)
      ENDIF
      IF(IPROFL.EQ.2) THEN
        GRSD(1,L)=GDINP(1,L)
        GRSR(1,L)=GRINP(1,L)
        GRSRD(1,L)=GRDINP(1,L)
      ENDIF
      IF(ICLAY.EQ.1) THEN
          RCLAY(1,L)=GRAV/RCINP(1,L)
          FCLAY(1,L)=1.D0-FCINP(1,L)/SPORO1
      ENDIF
      RCREST(L) = ZBRAW(1)
      DO 142 J = 2, JMAX(L)
        DO 143 K = 1, NBINP(L)-1
          IF((XB(J).GT.XBINP(K,L)).AND.(XB(J).LE.XBINP(K+1,L))) THEN
            ZBRAW(J) = ZBINP(K,L) + (XB(J)-XBINP(K,L))*SLOPE(K)
            FB2(J,L) = 0.5D0*FBINP(K,L)
            IF(ISEDAV.EQ.2) WMNODE(J,L)=WMINP(K,L)
            IF(IVEG.GE.1) THEN
              VEGFB(J,L)=VEGINP(K,L)
              VEGH(J,L)=VEGD(K,L)
              IF(IVEG.EQ.1) VEGRH(J,L)=VEGRD(K,L)
            ENDIF
            IF(IPROFL.EQ.2) THEN
              GRSD(J,L)=GDINP(K,L)
              GRSR(J,L)=GRINP(K,L)
              GRSRD(J,L)=GRDINP(K,L)
            ENDIF
            GOTO 144
          ENDIF
 143    CONTINUE
 144    DUM = ZBRAW(J) - RCREST(L)
        IF(IPROFL.EQ.0.AND.DUM.GE.0.D0) THEN
          RCREST(L) = ZBRAW(J)
          JCREST(L) = J
        ENDIF
        IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
          IF(J.EQ.2) ZPRAW(1)=ZPINP(1,L)
          DO 145 K=1,NPINP(L)-1
            IF((XB(J).GT.XPINP(K,L)).AND.(XB(J).LE.XPINP(K+1,L))) THEN
              ZPRAW(J)=ZPINP(K,L)+(XB(J)-XPINP(K,L))*PSLOPE(K)
              IF(ICLAY.EQ.1) THEN
                  RCLAY(J,L)=GRAV/RCINP(K,L)
                  FCLAY(J,L)=1.D0-FCINP(K,L)/SPORO1
              ENDIF
              GOTO 142
            ENDIF
 145      CONTINUE
        ENDIF
 142  CONTINUE
      IF(IVEG.GE.1)THEN
        DO 146 J=1,JMAX(L)
          VEGINP(J,L)=FB2(J,L)*VEGFB(J,L)
          IF(IVEG.EQ.1) THEN
            IF(VEGFB(J,L).EQ.0.D0) THEN
              UPROOT(J,L)=0.D0
            ELSE
              UPROOT(J,L)=1.D0
            ENDIF
          ENDIF
 146    CONTINUE
      ENDIF
C     VEGFB(J,L) used in wet zone (Main Program) and VEGINP(J,L) used 
C     in wet and dry zone (Subr.16 WETDRY). UPROOT(J,L)=0.0 in zone 
C     of no vegetation or uprooted vegetation
C     
C     Smooth ZBRAW(J) and ZPRAW(J) J=1-JMAX(L) using Subr.14 SMOOTH
      JMAXL=JMAX(L)
      CALL SMOOTH(JMAXL,ZBRAW,SLOPE)
      IF(IPERM.EQ.1.OR.ISEDAV.GE.1) CALL SMOOTH(JMAXL,ZPRAW,PSLOPE)
      DO 149 J=1,JMAX(L)
        ZB(J,L)=SLOPE(J)
        IF(IPROFL.GE.1) ZB0(J,L)=ZB(J,L)
        IF(IPERM.EQ.1.OR.ISEDAV.GE.1) ZP(J,L)=PSLOPE(J)
        IF(ICLAY.EQ.1) ZP0(J,L)=ZP(J,L)
        IF(ISEDAV.EQ.2) THEN
          IF(WMNODE(J,L).LE.0.D0) THEN
            ZMESH(J,L)=ZP(J,L)
          ELSE
            ZMESH(J,L)=ZB(J,L)
          ENDIF
        ENDIF
 149  CONTINUE
C     Calculate bottom slope and JCREST(if IPROFL=1 or 2) using 
C     smoothed ZB(J)
      BSLOPE(1,L) = (ZB(2,L) - ZB(1,L))/DX
      JMAXM1 = JMAX(L) - 1
      BSLOPE(JMAX(L),L) = (ZB(JMAX(L),L) - ZB(JMAXM1,L))/DX
      DO 150 J=2,JMAXM1
        BSLOPE(J,L) = (ZB(J+1,L) - ZB(J-1,L))/DX2
 150  CONTINUE
      IF(IPROFL.GE.1.AND.IPOND.EQ.0) THEN
        RCREST(L)=ZB(1,L)
        DO 151 J=2,JMAX(L)
          DUM=ZB(J,L)-RCREST(L)
          IF(DUM.GE.0.D0) THEN
            RCREST(L)=ZB(J,L)
            JCREST(L)=J
          ENDIF
 151    CONTINUE
      ENDIF
C     
C     HP(J,L) = vertical thickness of porous or sediment layer
      IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
        DO 210 J=1,JMAX(L)
          HP(J,L) = ZB(J,L) - ZP(J,L)
          IF(HP(J,L).LT.0.D0) THEN
            HP(J,L)=0.D0
            ZB(J,L)=ZP(J,L)
          ENDIF
          IF(ISTSAN.EQ.1) THEN
              ZBSTON(J,L)=ZB(J,L)
              ZPSTON(J,L)=ZP(J,L)
              HPSTON(J,L)=HP(J,L)
          ENDIF
 210    CONTINUE
      ENDIF
C     
C     If IVEG=1, VEGZD(J,L) and VEGZR(J,L) are the upper and lower 
C     elevations of non-uprooted vegetation at node J and line L
      IF(IVEG.EQ.1) THEN
        DO 220 J=1,JMAX(L)
          VEGZD(J,L)=ZB0(J,L)+VEGH(J,L)
          VEGZR(J,L)=ZB0(J,L)-VEGRH(J,L)
 220    CONTINUE
      ENDIF
C     where VEGZD(J,L) and VEGZR(J,L) are the same as the initial bottom 
C     elevation ZB0(J,L) in the zone of no vegetation with VEGH(J,L)=0.0
C     and VEGRH(J,L)=0.0
C     
 100  CONTINUE
C
      RETURN
      END
C     
C     -03----------------  END OF SUBROUTINE BOTTOM  ---------------------
C     #04#####################  SUBROUTINE PARAM  ########################
C     
C     This subroutine calculates parameters used in other subroutines
C     
      SUBROUTINE PARAM
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NL=100)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /PERIOD/ TP, WKPO, ANGLE, WT(NN)
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +   EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
C     
C     ... CONSTANTS and PARAMETER
C     
C     PI      = 3.14159
C     TWOPI   = 2.D0 * PI
C     GRAV    = acceleration due to gravity specified in Subr.2 INPUT
C     SQR2    = Sqrt(2)
C     SQR8    = Sqrt(8)
C     SQRG1   = Sqrt(2/PI)
C     SQRG2   = 2*Sqrt(2/PI)
C     WKPO    = deep water wave number for the representative period
C     
      PI = 3.14159D0
      TWOPI = 2.D0*PI
      SQR2 = DSQRT(2.D0)
      SQR8 = DSQRT(8.D0)
      SQRG1= DSQRT(2.D0/PI)
      SQRG2= 2.D0*SQRG1
      WKPO = (TWOPI)**2.D0/(GRAV*TP**2.D0)
C     
C.....POROUS FLOW RESISTANCE PARAMETERS IF IPERM=1
C     SNP = stone porosity specified in Subr.2 INPUT
C     SDP = nominal stone diameter specified in Subr.2 INPUT
C     WNU = kinematic viscosity of water (m*m/s)
C     WPM = maximum seepage velocity (m/s)
C     If INFILT=1, WPM is computed using SNP=SPORO and SDP=D50 of sand
C     in Subr.2 INPUT
      IF(IPERM.EQ.1.OR.INFILT.EQ.1) THEN
        WNU = 1.0D-6
        A = 1000.D0
        B = 5.D0
        IF(IPERM.EQ.1) THEN
          DUMP=SNP
          DUMD=SDP
        ENDIF
        IF(INFILT.EQ.1) THEN
          DUMP=1.D0-SPORO1
          DUMD=D50
        ENDIF
        C = 1.D0 - DUMP
        ALPHA = A*WNU*C**2.D0/(DUMP*DUMD)**2.D0
        BETA1 = B*C/DUMP**3.D0/DUMD
        BETA2 = 7.5D0*B*C/SQR2/DUMP**2.D0
C     Need to divide BETA2 by WT(J) in Subr.9 POFLOW
        ALSTA  = ALPHA/GRAV
        BESTA1 = BETA1/GRAV
        BESTA2 = BETA2/GRAV
        ALSTA2 = ALSTA*ALSTA
        BE2    = 2.D0*BESTA1
        BE4    = 2.D0*BE2
        WPM    = (DSQRT(ALSTA2+BE4)-ALSTA)/BE2
      ENDIF
C     
C.....SWASH PARAMETERS IN WET AND DRY ZONE IF IOVER=1
C     AWD = swash velocity parameter
C     AWD=2.0 calibrated for structures (IPROFL=0 or IPERM=1)
C     AWD=1.6 calibrated for wave overwash of sand dunes
C     EWD = duration-based exceedance probability for output
C     where AWD has not been calibrated extensively and
C     EWD=0.01-0.02 approximately corresponds to 2-percent exceedance
C     probability based on individual overtopping events.
      IF(IOVER.EQ.1) THEN
        AWD=2.0D0
        IF(IPROFL.EQ.1.AND.IPERM.EQ.0) AWD=1.6D0
        EWD = 0.015D0
        IF(IPERM.EQ.1) EWD=0.01D0
C     The following parameters are constant in Subr.16 WETDRY
        CWD= 0.75D0*DSQRT(PI)
        AQWD = CWD*AWD
        AGWD = AWD*AWD
        AUWD = 0.5D0*DSQRT(PI)*AWD
        BWD = (2.D0-9.D0*PI/16.D0)*AGWD + 1.D0
      ENDIF
C     
      RETURN
      END
C     
C     -04-----------------  END OF SUBROUTINE PARAM  ---------------------
C     #05#####################  SUBROUTINE LWAVE  ########################
C     
C     This subroutine calculates quantities based on linear wave theory
C     
      SUBROUTINE LWAVE(J, L, WD, QDISP)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000, NL=100)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY 
      COMMON /PERIOD/ TP, WKPO, ANGLE, WT(NN)
      COMMON /BINPUT/ XBINP(NN,NL), ZBINP(NN,NL), FBINP(NN,NL),XS(NL),
     +   YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX, FSY, FE, QWX, QWY
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN)
C     
C     ... LINEAR WAVE PARAMETERS 
C     
C     WD     = mean water depth from Main Program
C     QDISP  = water flux affecting wave period
C     TP     = representative wave period specified as input
C     WKP    = wave number at node J
C     WT(J)  = wave period at node J
C     CP(J)  = phase velocity based on input TP at node J
C     WN(J)  = ratio of group velocity to phase velocity at node J
C     
C     Solve linear wave dispersion relation with no current to find WKP
      IF(IWCINT.EQ.0.OR.QDISP.EQ.0.D0) THEN
        D = WD*WKPO
        IF(J.EQ.1) THEN 
          X = D/DSQRT(DTANH(D))
        ELSE
          X = WKP*WD
        ENDIF
 10     COTH = 1.D0/DTANH(X)
        XNEW = X - (X-D*COTH)/(1.D0+D*(COTH**2.D0-1.D0))
        IF (DABS(XNEW - X).GT.1.D-7) THEN
          X = XNEW
          GOTO 10
        ENDIF
        AF = TWOPI/TP
C     
C     Solve linear wave dispersion relation with current to find WKP
      ELSE
        B = TP*QDISP/TWOPI/WD/WD
        D = WD*WKPO
        IF(J.EQ.1) THEN
          X = D/DSQRT(DTANH(D))
        ELSE	
          X = WKP*WD
        ENDIF
 11     COTH = 1.D0/DTANH(X)
        C = 1.D0 - B*X
        F = X - D*C*C*COTH
        FD = 1.D0 + D*C*(2.D0*B*COTH + C*(COTH*COTH - 1.D0))
        XNEW = X - F/FD
        IF (DABS(XNEW - X).GT.1.D-7) THEN
          X = XNEW
          GOTO 11
        ENDIF
        AF = DSQRT(GRAV*XNEW*DTANH(XNEW)/WD)
      ENDIF
C     
      WKP = XNEW/WD
      X2 = X*2.D0
      WN(J) = 0.5D0*(1.D0 + X2/DSINH(X2))
      WT(J) = TWOPI/AF
      CP(J) = AF/WKP
      FSX = 2.D0*WN(J) - 0.5D0
      FSY = 0.D0
      FE  = WN(J)*CP(J)*WT(J)
C     
C     If IANGLE=0, normally incident waves 
      IF(IANGLE.EQ.0) THEN
        STHETA(J) = 0.D0
        CTHETA(J) = 1.D0
        GOTO 100
      ENDIF
C     
C     Otherwise, compute wave angle THETA in radians at node J using
C     Snell's Law where ANGLE = incident wave angle in degrees at
C     node J=1, AGLINE(L) = angle of cross-shore line L, and WKPSIN = constant
C     Wave angle from shorenormal is limited to range of -180 to 180 degrees
C     before imposing range of -80 and 80 degrees
C     
      IF(J.EQ.1) THEN
        DUM=ANGLE-AGLINE(L)
        IF(DUM.GT.180.D0) DUM=DUM-360.D0
        IF(DUM.LT.-180.D0) DUM=DUM+360.D0
        IF(DUM.GT.80.D0) DUM=80.D0
        IF(DUM.LT.-80.D0) DUM=-80.D0
        THETA = DUM*PI/180.D0
        STHETA(1) = DSIN(THETA)
        CTHETA(1) = DCOS(THETA)
        WKPSIN = WKP*STHETA(1)
      ELSE
        STHETA(J) = WKPSIN/WKP
        THETA = DASIN(STHETA(J))
        CTHETA(J) = DCOS(THETA)
      ENDIF
C     
      FSX = FSX - WN(J)*STHETA(J)*STHETA(J)
      FSY = WN(J)*STHETA(J)*CTHETA(J)
      FE = FE*CTHETA(J)
C     
 100  IF(IWCINT.EQ.1) FE=FE+WT(J)*QWX/WD
C     
C     Compute RX, RY and RE related to roller momentum and energy fluxes
C     as well as RBETA =wave-front slope of roller with RBZERO = 0.1
      IF(IROLL.EQ.1) THEN
        IF(IANGLE.EQ.0) THEN
          RX(J)=CP(J)/GRAV
          RE(J)=RX(J)*CP(J)
        ELSE
          DUM=CP(J)*CTHETA(J)/GRAV
          RX(J)=DUM*CTHETA(J)
          RY(J)=DUM*STHETA(J)
          RE(J)=DUM*CP(J)
        ENDIF
        RBETA(J)=RBZERO
        IF(BSLOPE(J,L).GT.0.D0) RBETA(J)=RBETA(J)+BSLOPE(J,L)*CTHETA(J)
      ENDIF
C     
      RETURN
      END
C     
C     -05-----------------  END OF SUBROUTINE  LWAVE ---------------------

C     #06############### SUBROUTINES GBXAGF and VSTGBY ###################
C     *****************SUBROUTINE GBXAGF**********************************
C     This subroutine computes GBX and GF for specified CTHETA, USIGT,
C     STHETA and VSIGT for Gaussian variable R
C     
      SUBROUTINE GBXAGF(CTHETA,USIGT,STHETA,VSIGT,GBX,GF)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NL=100)
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /CONSTA/ GRAV,SQR2, SQR8, PI, TWOPI,SQRG1, SQRG2
C     
C     For obliquelly incident waves, use approximate equations
      IF(IANGLE.EQ.1) THEN
        RM  = -USIGT*CTHETA - VSIGT*STHETA
        AFM = DABS(VSIGT*CTHETA - USIGT*STHETA)
        DUM = USIGT*USIGT + VSIGT*VSIGT
        GBX = SQRG1*(USIGT - RM*CTHETA)+ USIGT*AFM
        GF  = SQRG2 + (1.D0 + DUM)*AFM + SQRG1*(DUM + 2.D0*RM*RM)
      ENDIF
C     
C     For normally incident waves, use analytical 
C     expresions involving complementary error function ERFCC below
      IF(IANGLE.EQ.0) THEN
        C1 = 1.D0-ERFCC(USIGT/SQR2)
        C2 = SQRG1*DEXP(-USIGT*USIGT/2.D0)
        C3 = 1.D0 + USIGT*USIGT
        GBX = C3*C1 + C2*USIGT
        GF = USIGT*(C3 + 2.D0)*C1 + (C3 + 1.D0)*C2
      ENDIF
C     
      RETURN
      END
C     
C     -----------------------END OF SUBROUTINE GBXAGF---------------------
C     --------------------------SUBROUTINE VSTGBY-------------------------
C     This subroutine computes VSIGT= VMEAN/SIGT for specified GBY,
C     CTHETA, USIGT=UMEAN/SIGT, and STHETA but neglects USIGT*STHETA
C     
      SUBROUTINE VSTGBY(CTHETA,USIGT,STHETA,VSIGT,GBY)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
C     
      DUM=USIGT*STHETA
C     which is assumed to be zero
C     
      VSIGT = 0.D0
      IF(GBY.EQ.0.D0) GOTO 100
      B = SQRG1*(1.D0+STHETA*STHETA) 
      C = GBY 
      IF(GBY.GT.0.D0) THEN
        D=B*B + 4.D0*CTHETA*C
        IF(D.GE.0.D0) VSIGT=0.5D0*(DSQRT(D)-B)/CTHETA
        IF(VSIGT.LT.0.D0) VSIGT=0.D0
C     
      ELSE
        D = B*B-4.0D0*CTHETA*C
        IF(D.GE.0.D0) VSIGT=0.5D0*(B-DSQRT(D))/CTHETA
        IF(VSIGT.GT.0.D0) VSIGT=0.D0
      ENDIF
C
 100    CONTINUE
      RETURN
      END
C     
C     -------------------END OF SUBROUTINE VSTGBY-------------------------
C     ********************************************************************
      FUNCTION ERFCC(X)
      DOUBLE PRECISION X, Z, T, ERFCC
      Z=DABS(X)      
      T=1.D0/(1.D0+0.5D0*Z)
      ERFCC=T*DEXP(-Z*Z-1.26551223D0+T*(1.00002368D0+T*(.37409196D0+
     +   T*(.09678418D0+T*(-.18628806D0+T*(.27886807D0+
     +   T*(-1.13520398D0+T*(1.48851587D0+
     +   T*(-.82215223D0+T*.17087277D0)))))))))
      IF (X.LT.0.D0) ERFCC=2.D0-ERFCC
      RETURN
      END
C     *********************************************************************
C     -06------------  END OF SUBROUTINES GBXAGF and VSTGBY  --------------
C     #07#####################  SUBROUTINE DBREAK  ########################
C     
C     This subroutine calculates QBREAK and DBSTA for wave breaking 
C     
      SUBROUTINE DBREAK(J, L, WHRMS, D)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NL=100)
C     
      COMMON /PERIOD/ TP, WKPO, ANGLE, WT(NN)
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX, FSY, FE, QWX, QWY
      COMMON /WBREAK/ GAMMA, QBREAK(NN), DBSTA(NN), SISMAX, ABREAK(NN)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
c mg
C     COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
C mg  
C     Calculate energy dissipation factor ABREAK(J) for steep slope
C     where D = mean water depth from Main Program
      ABREAK(J) = (TWOPI/WKP/D)*BSLOPE(J,L)*CTHETA(J)/3.D0
      IF(ABREAK(J).LT.1.D0) ABREAK(J) = 1.D0
C mg
C mg  Allow for variable gamma
C     IF(GAMMA.LT.0) THEN
C mg  Compute deep water wave height
C       CO = GRAV*TP/TWOPI
C         THETAO=DASIN(CO/CP(1)*STHETA(1))
C       HRMSO = HRMS(1)*DSQRT((CP(1)*WN(1)*CTHETA(1))/
C    +    (0.5D0*CO*DCOS(THETAO)))
C mg  Alex Apotsos et al. 2008, Coastal Engineering 55 (2008) 224-235.  (Eq 23)
C       GAMMA_TEMP = 0.3 + 0.45*TANH(0.9*HRMSO)
C     ELSE
C       GAMMA_TEMP = GAMMA
C     ENDIF
C mg  
C     ... FRACTION OF BREAKING WAVES AND ASSOCIATED DISSIPATION
C     
C     QBREAK(J) = Fraction of breaking waves at node J
C     DBSTA(J)  = Time averaged normalized energy dissipation due to 
C     wave breaking at node J
C mg     
      HM = 0.88D0/WKP*DTANH(GAMMA*WKP*D/0.88D0)
C mg  HM = 0.88D0/WKP*DTANH(GAMMA_TEMP*WKP*D/0.88D0)
C mg    
C     Compute QBREAK = fraction of breaking waves
      B = (WHRMS/HM)**2.D0
C     IF(B.LT.0.99999D0) THEN
      IF(B.LT.0.99999D0.AND.WHRMS.GT.1.D-10) THEN !bdj
        QBOLD = B/2.D0
 10     QBREAK(J) = QBOLD - (1.D0-QBOLD + B*DLOG(QBOLD))/(B/QBOLD-1.D0)
        IF(QBREAK(J).LE.0.D0) QBREAK(J) = QBOLD/2.D0
        IF(DABS(QBREAK(J)-QBOLD).GT.1.D-6) THEN
          QBOLD = QBREAK(J)
          GOTO 10
        ENDIF
      ELSE
        QBREAK(J) = 1.D0
        IF(WHRMS.LE.1.D-10) QBREAK(J)=0.D0
        HM=WHRMS
      ENDIF
C     
      DBSTA(J) = 0.25D0*ABREAK(J)*QBREAK(J)*HM*HM/WT(J)
C     
C     Reduce SIGSTA if WHRMS is larger than GAMMA*D 
C     (used only for wave transmission over submerged breakwater)
C     GAMD = GAMMA*D
C     IF(WHRMS.LE.GAMD) THEN
      SISMAX = 1.D0
C     ELSE
C     SISMAX = DSQRT(GAMMA*WHRMS/D/8.0D0)
C     ENDIF
C     
      RETURN
      END
C     
C     -07-----------------  END OF SUBROUTINE DBREAK  ---------------------
C     #08#####################  SUBROUTINE OUTPUT  ########################
C     
C     This subroutine stores computed and input quantities
C     
      SUBROUTINE OUTPUT(ITIME,L,ITEQO,ICONV)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000, NB=30000,NL=100)
      DIMENSION DUMVEC(NN),EDEPTH(NN)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +  ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +  IVEG,ICLAY
      COMMON /PERIOD/ TP, WKPO, ANGLE, WT(NN)
      COMMON /SEAWAV/ TIMEBC(NB), TPBC(NB), HRMSBC(NB),
     +  WSETBC(NB), SWLBC(NB), WANGBC(NB), NWAVE, NSURG,
     +  NWIND, NTIME 
      COMMON /PREDIC/ HRMS(NN), SIGMA(NN), H(NN), WSETUP(NN), SIGSTA(NN)
      COMMON /BINPUT/ XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),
     +  YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +  SWLDEP(NN,NL), BSLOPE(NN,NL), JMAX(NL), JSWL(NL)
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
      COMMON /LINEAR/ WKP, CP(NN), WN(NN), WKPSIN, STHETA(NN), 
     +  CTHETA(NN), FSX, FSY, FE, QWX, QWY
      COMMON /FRICTN/ GBX(NN), GBY(NN), GF(NN)
      COMMON /WBREAK/ GAMMA, QBREAK(NN), DBSTA(NN), SISMAX, ABREAK(NN)
      COMMON /CRSMOM/ SXXSTA(NN), TBXSTA(NN)
      COMMON /LOGMOM/ SXYSTA(NN), TBYSTA(NN)
      COMMON /ENERGY/ EFSTA(NN), DFSTA(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /VELOCY/ UMEAN(NN), USTD(NN), USTA(NN), VMEAN(NN),  
     +  VSTD(NN), VSTA(NN)
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +  EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /SEDOUT/ PS(NN), VS(NN), QSX(NN), QSY(NN),
     +  PB(NN), GSLOPE(NN), QBX(NN), QBY(NN), Q(NN)
      COMMON /SEDVOL/ VBX(NN,NL),VSX(NN,NL),VBY(NN,NL),VSY(NN,NL),
     +  VY(NN,NL),DZX(NN,NL)  
      COMMON /ROLLER/ RBZERO, RBETA(NN), RQ(NN), RX(NN), RY(NN), RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL), 
     +  WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +  UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /WIND/   W10(NB), WANGLE(NB), WINDCD(NB), TWXSTA(NB),
     +  TWYSTA(NB)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +  BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +  USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +  H1,JWD,JDRY
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
      COMMON /TIDALC/ DETADY(NB),DSWLDT(NB)
      COMMON /SERIES/TSQO(NL),TSQBX(NL),TSQSX(NL)
      COMMON /VEGETA/VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
      COMMON /DIKERO/EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN),
     +  GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL),GRSR(NN,NL),
     +  GRSRD(NN,NL),DEEB,DEEF
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
      COMMON /SOCLAY/EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL),
     + FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL)
C     
C     ......... OUTPUT ONLY WHEN ITIME=0.............................
C     
      IF(ITIME.EQ.0) THEN
      IF(L.GT.1) GOTO 888
C     
C     ......... COMPUTATIONAL OPTION
C     ILINE=number of cross-shore lines
      WRITE(20,890) ILINE,IQYDY
 890  FORMAT('COMPUTATION OPTION ILINE=',I3/
     +    'Alongshore gradient IQYDY=',I3/
     +    'ILINE cross-shore lines are computed together'/)
C
C     IPROFL=0 for fixed bottom profile
C     IPROFL=1 for profile evolution computation
C     IPROFL=2 for dike erosion computation
        IF(IPROFL.EQ.0) THEN
          WRITE(20,900) IPROFL
        ENDIF
 900    FORMAT('COMPUTATION OPTION IPROFL =',I3/
     +     'Bottom profile is fixed and No sediment transport
     +     is computed'/)
C     
        IF(IPROFL.EQ.1) THEN
          WRITE(20,901) IPROFL,TIMEBC(NTIME+1),NTIME
          IF(ISEDAV.EQ.1.AND.ICLAY.EQ.0) WRITE(20,902) ISEDAV,BEDLM
          IF(ICLAY.EQ.1) WRITE(20,904) ICLAY,BEDLM,DEEB,DEEF
          IF(ISEDAV.EQ.2) WRITE(20,905) ISEDAV,BEDLM
        ENDIF
 901    FORMAT('COMPUTATION OPTION IPROFL =',I3/
     +     'Profile evolution is computed from Time = 0.0'/
     +     'to Time = ',F13.1,'  for NTIME = ', I4/)
 902    FORMAT('ISEDAV=',I3,' for hard bottom with', 
     +    'bedload reduction factor BEDLM=',F4.1/)
 904    FORMAT('ICLAY=',I3,'for sand on erodible clay bottom'/
     +    'with bedload reduction factor BEDLM=',F4.1/
     +    'erosion efficiency DEEB=',F6.4/
     +    'erosion efficiency DEEF=',F6.4/)
 905    FORMAT('ISEDAV=',I3,' for wire mesh with', 
     +    'bedload reduction factor BEDLM=',F4.1/)
C     
        IF(IPROFL.EQ.2) THEN
          WRITE(20,903) IPROFL,TIMEBC(NTIME+1),NTIME,DEEB,DEEF
        ENDIF
 903    FORMAT('COMPUTATIONAL OPTION IPROFL=',I3/
     +     'Dike erosion is computed from Time=0.0'/
     +     'to Time=',F13.1,'  for NTIME=',I4/
     +     'Efficiency DEEB=',F6.4/
     +     'Efficiency DEEF=',F6.4/)
C     
        IF(IROLL.EQ.0) WRITE(20,910)
        IF(IROLL.EQ.1) WRITE(20,911) RBZERO
 910    FORMAT('NO ROLLER is included in computation'/)
 911    FORMAT('ROLLER is included in computation'/
     +     'ROLLER slope Betazero =', F9.3/)
C     
        IF(IWCINT.EQ.0) WRITE(20,920)
        IF(IWCINT.EQ.1) WRITE(20,921)
 920    FORMAT('NO wave and current interaction included'/)
 921    FORMAT('WAVE and current interaction included'/)
C     
        IF(IOVER.EQ.0) WRITE(20,930)
        IF(IOVER.EQ.1.AND.IPOND.EQ.0) THEN
          WRITE(20,931) RWH,JCREST(L),RCREST(L),AWD,EWD
        ENDIF
        IF(IOVER.EQ.1.AND.IPOND.EQ.1) WRITE(20,932) RWH,AWD,EWD,ZW
 930    FORMAT('NO wave overtopping, overflow and seepage'/)
 931    FORMAT('WAVE OVERTOPPING, OVERFLOW AND SEEPAGE'/
     +     'Runup wire height (m)               RWH=',F9.3/
     +     'Initial crest location for L=1      JCREST=',I6/
     +     'Initial crest height (m) for L=1    RCREST=',F9.3/
     +     'Swash velocity parameter            AWD=',F9.3/
     +     'Output exceedance probability       EWD=',F9.3/)
 932    FORMAT('PONDED WATER IN RUNNEL'/
     +   'Runup wire height (m)               RWH=',F9.3/
     +   'Swash velocity parameter            AWD=',F9.3/
     +   'Output exceedance probability       EWD=',F9.3/
     +   'Initial ponded water level (m)       ZW=',F9.3/)
C     
        IF(IPERM.EQ.0) WRITE(20,940)
        IF(IPERM.EQ.1) WRITE(20,941) SNP,SDP,CSTABN,WNU,WPM
        IF(ISTSAN.EQ.1) WRITE(20,942) CPSTON
 940    FORMAT('IMPERMEABLE BOTTOM assumed'/)
 941    FORMAT('PERMEABLE BOTTOM consisting of'/
     +     'Stone porosity                       SNP=',F9.3/
     +     'Nominal stone diameter (m)         DN50=',F9.4/
     +     'Critical stability number        CSTABN=',F9.3/
     +     'Water kinematic viscosity(m*m/s)       =',F9.7/
     +     'Maximum seepage velocity (m/s)      WPM=',F9.5/)
 942    FORMAT('ISTSAN=1 for fixed stone structure on sand bottom'/
     +     'Empirical parameter CPSTON=',F5.2/)
C     
        IF(IWIND.EQ.0) WRITE(20,950)
        IF(IWIND.EQ.1) WRITE(20,951) NWIND
 950    FORMAT('NO wind shear stresses included'/)
 951    FORMAT('WIND shear stresses in momentum equations'/
     +     'Number of wind speed and direction input =',I4/)
C     
        IF(IVEG.EQ.0) WRITE(20,955)
        IF(IVEG.EQ.1) WRITE(20,956) VEGCD
        IF(IVEG.EQ.2) WRITE(20,957) VEGCD
 955    FORMAT('NO vegetation in computation domain'/)
 956    FORMAT('VEGETATION whose density, width, height and root depth 
     +    are specified as input'/'Vegetation drag coefficient VEGCD=',
     +    F5.2/)
 957    FORMAT('VEGETATION whose density,width and height are
     +    specified as input'/'Vegetation drag coefficient VEGCD=',
     +    F5.2/) 
C     
        WRITE(20,960) GAMMA
 960    FORMAT('BREAKER RATIO PARAMETER'/
     +     'Input Gamma =',F5.2/)
C     
        IF(IPROFL.EQ.1) WRITE(20,970) D50*1000.D0,WF,SG,EFFB,SLP,
     +     TANPHI,BLP
        IF(IPROFL.EQ.1.AND.IOVER.EQ.1) THEN
          WRITE(20,971) SLPOT
          IF(INFILT.EQ.1) WRITE(20,972) WPM
        ENDIF
 970    FORMAT('SEDIMENT PARAMETERS if IPROFL=1'/
     +     'Median diameter (mm)                D50=   ',F8.2/
     +     'Fall velocity (m/s)                  Wf=    ',F6.4/
     +     'Specific gravity                     Sg=   ',F5.2/
     +     'Suspension efficiency                eB=   ',F6.3/
     +     'Suspended load parameter               =   ',F5.2/
     +     'Tangent(friction angle)                =   ',F5.2/
     +     'Bedload parameter                     b=    ',F6.4)
 971    FORMAT('Susp.load para. (IOVER=1)              =   ',F5.2/)
 972    FORMAT('INFILT=1 and infiltr. velocity (m/s)   =   ',F7.5/)
C     
C     ......... INPUT WAVE AND WATER LEVEL
c mg   nout = 10000
       nout = 10
        WRITE(20,1001) NTIME, nout
        NTIM9=nout+1
        IF(NTIME.GT.(2*nout)) NTIM9=NTIME-(nout-1)
        DO 130 I = 1,NTIME
          IF(I.LE.nout.OR.I.GE.NTIM9) THEN
            WRITE(20,1002) TIMEBC(I+1),TPBC(I),HRMSBC(I),
     +          WSETBC(I),SWLBC(I), WANGBC(I)
          ENDIF
 130    CONTINUE
 1001   FORMAT(/'INPUT WAVE AND WATER LEVEL:NTIME=',I6,' from TIME=0.0'/
     +     'Output lines are limited to first and last',I6,' lines'/
     +     'End Time(sec) Tp (sec)  Hrms(m) Wave Setup(m)',
     +     'Storm tide(m) ANGLE(deg)'/)
 1002   FORMAT(F11.1,5F11.4)
  888  CONTINUE
C     End of L=1 output................................................
C     
C     ......... OUTPUT BOTTOM GEOMETRY
C     The bottom geometry is divided into segments of
C     different inclination and roughness starting from
C     seaward boundary for ILINE cross-shore lines.
C     NBINP(L)  = number of input points for L line
C     XBINP(J,L)  = horizontal distance from seaward boundary
C     to landward-end of segment (J-1) in meters
C     ZBINP(J,L)  = dimensional vertical coordinate (+ above datum)
C     of the landward end of segment (J-1) in meters
C     FBINP(J,L) = bottom friction factor
        WRITE (20,1100) L, YLINE(L), AGLINE(L), 0.D0-ZBINP(1,L),
     +     NBINP(L), DX, XS(L), JMAX(L)
C     
 1100   FORMAT (/'INPUT BEACH AND STRUCTURE GEOMETRY'/
     +     'Cross-shore line number               L=  ',I3/
     +     'Alongshore coordinate             YLINE=  ',F13.4/
     +     'Line angle(degrees)              AGLINE=  ',F13.4/
     +     'Depth at seaward boundary (m)          =  ',F13.6/
     +     'Number of input points            NBINP=  ',I8/
     +     'Output lines are limited to first and last 5 lines'/
     +     'Node spacing (m)                    DX=  ',F13.6/
     +     'Shoreline location (m) of Zb=0       Xs=  ',F13.6/
     +     'Maximum landward node              JMAX=',I8//
     +     'X  (m)      Zb  (m)  Fric.factor  Wire mesh')
        WRITE (20,1200) XBINP(1,L), ZBINP(1,L)
        NBINP4=6
        IF(NBINP(L).GT.10) NBINP4=NBINP(L)-4
        DO 140 J = 2,NBINP(L)
          IF(J.LE.5.OR.J.GE.NBINP4) THEN
            IF(ISEDAV.LE.1) THEN
              WRITE (20,1200) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L)
            ELSE
              WRITE (20,1202) XBINP(J,L), ZBINP(J,L), FBINP(J-1,L),
     +           WMINP(J-1,L)
            ENDIF
          ENDIF
 140    CONTINUE
        IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
          IF(ICLAY.EQ.0) THEN
              WRITE(20,1150) L,NPINP(L)
          ELSE
              WRITE(20,1151) L,NPINP(L)
          ENDIF
          NPINP4=6
          IF(NPINP(L).GT.10) NPINP4=NPINP(L)-4
          DO 141 J=1,NPINP(L)
            IF(J.LE.5.OR.J.GE.NPINP4) THEN
              IF(ICLAY.EQ.0) THEN
                  WRITE(20,1201) XPINP(J,L), ZPINP(J,L)
              ELSE    
                  WRITE(20,1202) XPINP(J,L),ZPINP(J,L),RCINP(J,L),
     +             FCINP(J,L)
              ENDIF 
            ENDIF
 141      CONTINUE
        ENDIF
C     where the number of the output lines is limited to 10 or less
C     to reduce the length of the output file ODOC.
 1150   FORMAT(/'INPUT IMPERMEABLE HARD BOTTOM GEOMETRY'/
     +     'Number of input points for line L=',I3, ' NPINP= ',I5//
     +     'X  (m)        ZP  (m)  ')
 1151   FORMAT(/'INPUT ERODIBLE CLAY BOTTOM ELEVATION'/
     +     'Number of input points for line L=',I3, 'NPINP= ',I5//
     +     'X(m)      ZP(m)    RC(m*m/s/s),sand frac  ') 
 1200   FORMAT(3F10.3)
 1201   FORMAT(2F10.3)
 1202   FORMAT(4F10.3)
C
C.....OUTPUT VEGETATION CHARACTERISTICS FOR IVEG=1 or 2
      IF(IVEG.GE.1) THEN
        IF(IVEG.EQ.1) THEN
          WRITE(20,1161)
        ELSE
          WRITE(20,1160)
        ENDIF
        DO 135 J=2,NBINP(L)
          IF(J.LE.5.OR.J.GE.NBINP4) THEN
            J1=J-1
            IF(IVEG.EQ.1) THEN
              WRITE(20,1203) XBINP(J,L),VEGN(J1,L),VEGB(J1,L),
     +        VEGD(J1,L),VEGRD(J1,L)
            ELSE
              WRITE(20,1202) XBINP(J,L),VEGN(J1,L),VEGB(J1,L),VEGD(J1,L)
            ENDIF
          ENDIF
 135    CONTINUE
      ENDIF
 1160 FORMAT(/'INPUT VEGETATION CHARACTERISITCS'/
     +   'X (m)   DENSITY(1/m/m) WIDTH(m) HEIGHT(m) ')
 1161 FORMAT(/'INPUT VEGETATION CHARACTERISITCS'/
     +   'X(m)  DEN.(1/m/m) WID.(m) HEI.(m) ROD.(m)')
 1203 FORMAT(5F8.3)
C     
C.....OUTPUT DIKE GRASS AND SOIL CHARACTERISTICS FOR IPROFL=2
      IF(IPROFL.EQ.2) THEN
        WRITE(20,1170)
        DO 136 J=2,NBINP(L)
          IF(J.LE.5.OR.J.GE.NBINP4) THEN
            J1=J-1
            WRITE(20,1210) XBINP(J,L),GDINP(J1,L),GRINP(J1,L),
     +         GRDINP(J1,L)
          ENDIF
 136    CONTINUE
      ENDIF
 1170 FORMAT(/'INPUT GRASS AND SOIL CHARACTERISTICS'/
     +   'X (m)   THICKNESS(m)  RO(m*m/s/s)  RD(m*m/s/s)')
 1210 FORMAT(4F11.3)
C     
C.....INPUT WIND SHEAR STRESSES FOR IWIND=1
      IF(L.GT.1) GOTO 889
        IF(IWIND.EQ.1) THEN
          WRITE(20,1370)
          DO 142 I=1,NTIME
            IF(I.LE.10.OR.I.GE.NTIM9) THEN
              WRITE(20,1371) TIMEBC(I),TIMEBC(I+1),W10(I),WANGLE(I),
     +           WINDCD(I)
            ENDIF
 142      CONTINUE
        ENDIF
 1370   FORMAT(/'INPUT WIND SPEED, DIRECTION AND STRESSES'/
     +     'Start & End Time(s) Speed(m/s) Dir(deg) DragCoef'/)
 1371   FORMAT(2F11.1,2F11.2,E11.4)
C     
C.....INPUT LANDWARD STILL WATER LEVEL FOR IWTRAN=1
        IF(IWTRAN.EQ.1) THEN
          IF(ISWLSL.EQ.0) WRITE(20,1380)
          IF(ISWLSL.EQ.1) THEN
            WRITE(20,1381)
            DO 143 I=1,NTIME
              IF(I.LE.10.OR.I.GE.NTIM9) THEN
                WRITE(20,1382) TIMEBC(I),TIMEBC(I+1),SWLAND(I)
              ENDIF
 143        CONTINUE
          ENDIF
          IF(ISWLSL.EQ.2) WRITE(20,1383)
        ENDIF
 1380   FORMAT(/'INPUT LANDWARD STILL WATER LEVEL for IWTRAN=1 and ',
     +     'ISWLSL=0'/'same as input seaward still water level'/)
 1381   FORMAT(/'INPUT LANDWARD STILL WATER LEVEL for IWTRAN=1 and ',
     +     'ISWLSL=1'/'Start & End Time(s) SWL(m) above datum'/)
 1382   FORMAT(2F11.1,F11.4)
 1383   FORMAT(/'IWTRAN=1 but ISWLSL=2 and NO WATER landward of crest'/
     +     'Overflow occurs (IOFLOW=1) if crest is submerged'/)
C     
C.....INPUT ALONGSHORE WATER LEVEL GRADIENT FOR ITIDE=1
        IF(ITIDE.EQ.1) THEN
          WRITE(20,1390)
          DO 144 I=1,NTIME
            IF(I.LE.10.OR.I.GE.NTIM9) THEN
              WRITE(20,1385) TIMEBC(I),TIMEBC(I+1),DETADY(I)
            ENDIF
 144      CONTINUE
 1385     FORMAT(2F11.1,F11.7)
        ENDIF
 1390   FORMAT(/'INPUT ALONGSHORE WATER LEVEL GRADIENT'/
     +          'Start & End Time(s)      DETA/DY alongshore'/)
C     
C     End of L=1 OUTPUT.....................................................
 889    CONTINUE
      ENDIF
C     --------------------- END OF OUTPUT ONLY WHEN ITIME = 0 --------------
C     
C     ------------------- COMPUTED CROSS-SHORE VARIATIONS ------------------
C     For each cross-shore line L of ILINE lines
C     Stored at Time = 0.0 and at the end of constant wave and
C     water level at the seaward boundary if laboratory data (ILAB=1)
C     For field data (ILAB=0), stored at the beginning, end, and
C     every ten storage time levels (GOTO 200 goes to end of this subr.)
C     
      IF(ITIME.EQ.0) THEN
        WRITE(21,1490)L,JMAX(L),TIMEBC(ITIME)
        DO 199 J=1,JMAX(L)
          IF(IPERM.EQ.0.AND.ISEDAV.EQ.0) THEN
            WRITE(21,1500)XB(J),ZB(J,L)
          ELSE
            IF(ISEDAV.EQ.1.OR.IPERM.EQ.1) THEN
              WRITE(21,1500) XB(J),ZB(J,L),ZP(J,L)
            ENDIF
            IF(ISEDAV.EQ.2) THEN
              WRITE(21,1500) XB(J),ZB(J,L),ZMESH(J,L),ZP(J,L)
            ENDIF
          ENDIF
 199  CONTINUE
        GOTO 200
      ENDIF
C     
      TIMOUT = TIME
C mg - explicit declaration of laboratory/field data sets
      IF(ILAB.EQ.0) THEN
C
C mg - ensure output at end of simulation for field data sets
        IF(ITIME.EQ.NTIME) GOTO 201
C mg
C        DUM=DBLE(ITIME)/10.D0-DBLE(ITIME/10)
C        IF(DUM.NE.0.D0) GOTO 200
      ENDIF
 201  CONTINUE
      IF(IPROFL.EQ.0) THEN
        TIMOUT = TIMEBC(ITIME+1)
        WRITE(20,1440) TIMOUT,L
      ELSE
        WRITE(20,1441) TIMOUT,L
      ENDIF
 1440 FORMAT(/'**********COMPUTED CROSS-SHORE VARIATIONS**********'/
     +   'on input bottom profile at TIME =',F11.1, '  Line L=',I3/)
 1441 FORMAT(/'**********COMPUTED CROSS-SHORE VARIATIONS**********'/
     +   'on bottom profile computed at TIME (s) = ', F11.1, 
     +   '  Line L=',I3/) 
C     
      WRITE(20,1450) JR, XB(JR), ZB(JR,L), H(JR)
 1450 FORMAT('LANDWARD WET COMPUTATION LIMIT'/
     +   'Most landward node of wet zone computation    JR=',I8/
     +   'X-coordinate at JR (m)                        XR=  ',F13.6/
     +   'Bottom elevation at JR (m)                    ZR=  ',F13.6/
     +   'Mean water depth at this node (m)          H(JR)=  ',F13.6/)
C     
C     Wave Reflection Coeffiecient at node J=1 only for IOVER=0
      IF(IOVER.EQ.0) THEN
        IF(JR.GT.JSWL(L).AND.JSWL(L).LT.JMAX(L)) THEN
          DUM = SIGMA(JSWL(L))*SIGMA(JSWL(L))*CP(JSWL(L))*WN(JSWL(L))
          DUM = DUM/WN(1)/CP(1)
          SIGREF=DSQRT(DUM)
          IF(IANGLE.EQ.1) SIGREF=SIGREF/DSQRT(CTHETA(1)/CTHETA(JSWL(L)))
          REFCOF=SIGREF/SIGMA(1)
          WRITE(20,1460) REFCOF, JSWL(L)
        ENDIF
      ENDIF
 1460 FORMAT('WAVE REFLECTION COEFFICIENT'/
     +   'Wave reflection coefficient (at x=0) = ',F9.6/
     +   'Still water shoreline node location JSWL =',I5/)
C     
C     Output computed wave overtopping, overflow and seepage rates
C     in Subr.10 QORATE
      IF(IOVER.EQ.1)THEN
        IF(IWTRAN.EQ.0.OR.JR.LT.JMAX(L))THEN
          CALL QORATE(ITIME,L,ITEQO,ICONV,1)
        ENDIF
      ENDIF
      IF(JR.EQ.JMAX(L).AND.IWTRAN.EQ.1)THEN
        DUM=SIGMA(JMAX(L))/SIGMA(1)
        WRITE(20,1461) DUM,JMAX(L),RCREST(L),JCREST(L)
 1461   FORMAT('WAVE TRANSMISSION OVER SUBMERGED STRUCTURE'/
     +  'Wave transmission coefficient =',F9.6/
     +  '    at landward end node JMAX=',I5/
     +  'Structure crest elevation (m),RCREST=',F9.4/
     +  '    at crest node JCREST=',I5/)
      ENDIF
C     
C     Longshore (Suspended Sediment plus Bedload) Transport Rate
      IF(IPROFL.EQ.1.AND.IANGLE.EQ.1) THEN
        DUM = 0.5D0*(QBY(1)+QSY(1))
        DO 145 J = 2,JDRY-1
          DUM = DUM + (QBY(J)+QSY(J))
 145    CONTINUE
        DUM = DUM + 0.5D0*(QBY(JDRY)+QSY(JDRY))
        QLONG = DUM * DX 
        SIGMAX = SIGMA(1)
        JB=1
        DO 150 J=2,JR
          IF(SIGMA(J).GT.SIGMAX) THEN
            SIGMAX = SIGMA(J)
            JB = J
          ENDIF
 150    CONTINUE
        DUM = SIGMA(JB)**2.D0*CP(JB)*WN(JB)*CTHETA(JB)*STHETA(JB)
        CERCK = (SG-1.D0)*QLONG/DUM
        WRITE(20,1470) QLONG,CERCK,STHETA(JB)
      ENDIF
 1470 FORMAT('LONGSHORE SUSPENDED AND BEDLOAD SAND TRANSPORT RATE'/
     +'Transport Rate (m**3/s) =',E14.5/'CERC Formula K=',F11.3/
     +'sin(breaker angle)=',F11.5/)
C     
C     Damage (normalized eroded area) of stone structure
C     EDMAX = normalized maximum vertical erosion depth
      IF(ISTSAN.EQ.0) THEN
          IF(IPROFL.EQ.1.AND.IPERM.EQ.1) THEN
              EDMAX=0.D0
              DO 300 J=1,JMAX(L)
                  EDEPTH(J)=ZB0(J,L)-ZB(J,L)
                  IF(EDEPTH(J).GT.EDMAX) EDMAX=EDEPTH(J)
                  IF(EDEPTH(J).LT.0.D0) EDEPTH(J)=0.D0
 300      CONTINUE
          EDMAX=EDMAX/SDP
          JMAXL=JMAX(L)
          CALL INTGRL(JMAXL,DX,EDEPTH,AREA)
          DAMAGE=AREA/SDP/SDP
          STABNO=SQR2*HRMS(1)/SDP/(SG-1.D0)
          WRITE(20,1480) DAMAGE,EDMAX,STABNO
          ENDIF
      ENDIF
 1480 FORMAT('STONE STRUCTURE DAMAGE'/
     +'Damage S=',F10.3/ 'Normalized erosion depth E=',F10.3/
     +'Stability number Nmo=',F8.3/)
C     
C.........COMPUTED CROSS-SHORE VARIATIONS
C     
C     Indicate the number of lines used for storage at specified time
C     for each cross-shore line L=1,2,...,ILINE
      JSWASH = JDRY - JWD +1
      JDUM = JR
      IF(IOVER.EQ.1) THEN
        JDUM=JDRY
        IF(IWTRAN.EQ.1.AND.JSL.LT.JMAX(L)) JDUM=JMAX(L)
      ENDIF
      WRITE(22,1490) L,JDUM,TIMOUT
      WRITE(23,1490) L,JR,TIMOUT
      WRITE(24,1490) L,JR,TIMOUT
      WRITE(25,1490) L,JR,TIMOUT
      WRITE(26,1490) L,JR,TIMOUT
      WRITE(27,1490) L,JDUM,TIMOUT
      WRITE(28,1490) L,JDUM,TIMOUT
      WRITE(29,1490) L,JR,TIMOUT
      WRITE(30,1490) L,JDUM,TIMOUT
      WRITE(31,1490) L,JR,TIMOUT
      WRITE(32,1490) L,JMAX(L),TIMOUT
      WRITE(33,1490) L,JMAX(L),TIMOUT
      WRITE(37,1490) L,JMAX(L),TIMOUT
      WRITE(38,1490) L,JMAX(L),TIMOUT
      WRITE(39,1490) L,JMAX(L),TIMOUT
      IF(IOVER.EQ.1) THEN
        WRITE(34,1490)L,JDUM,TIMOUT
        WRITE(35,1490)L,JSWASH,TIMOUT
        TIMID=0.5D0*(TIMEBC(ITIME)+TIMEBC(ITIME+1))
        DUM=TIMEBC(ITIME+1)-TIMEBC(ITIME)
        WRITE(36,1491) L,TIMID,(TSQO(L)/DUM),(TSQBX(L)/DUM),
     +     (TSQSX(L)/DUM)
      ENDIF
 1490 FORMAT(2I8,F11.1)
 1491 FORMAT(I8,4F17.9)
C     
      IF(IPROFL.GE.1.AND.L.EQ.ILINE) THEN
        DO 181 LL=1,ILINE     
          WRITE(21,1490) LL,JMAX(LL),TIMOUT
          DO 180 J=1,JMAX(LL)
            IF(IPERM.EQ.0.AND.ISEDAV.EQ.0) THEN
              IF (IVEG.EQ.1) THEN
                WRITE(21,1500) XB(J),ZB(J,LL),UPROOT(J,LL)
              ELSE
                WRITE(21,1500) XB(J),ZB(J,LL)
              ENDIF
            ELSE
              IF(IVEG.EQ.1) THEN
                WRITE(21,1500) XB(J),ZB(J,LL),ZP(J,LL),UPROOT(J,LL)
              ELSE
                IF(ISEDAV.EQ.1.OR.IPERM.EQ.1) THEN
                  IF(ISTSAN.EQ.0) WRITE(21,1500) XB(J),ZB(J,LL),ZP(J,LL)
                ENDIF
                IF(ISEDAV.EQ.2) THEN
                  WRITE(21,1500) XB(J),ZB(J,LL),ZMESH(J,LL),ZP(J,LL)
                ENDIF
                IF(ISTSAN.EQ.1) THEN
                  WRITE(21,1500) XB(J),ZB(J,LL),ZP(J,LL),VDSAND(J)
                ENDIF
              ENDIF
            ENDIF
 180      CONTINUE
 181    CONTINUE
      ENDIF
C     
C     Smooth computed PB(J), VS(J), and PS(J) before storing and plotting
      IF(IPROFL.EQ.1) THEN
        DUMVEC=PB
        CALL SMOOTH(JDUM,DUMVEC,PB)
        DUMVEC=VS
        CALL SMOOTH(JDUM,DUMVEC,VS)
        DUMVEC=PS
        CALL SMOOTH(JDUM,DUMVEC,PS)
      ENDIF
C     
      DO 160 J = 1,JR
        WRITE(22,1500) XB(J),(H(J)+ZB(J,L)),H(J),SIGMA(J)
        WRITE(23,1500) XB(J),WT(J),QBREAK(J),SIGSTA(J)
        WRITE(24,1500) XB(J),SXXSTA(J),TBXSTA(J)
        IF(IANGLE.EQ.1) WRITE(25,1500) XB(J),SXYSTA(J),TBYSTA(J)
        WRITE(26,1500) XB(J),EFSTA(J)/WT(J),DBSTA(J),DFSTA(J)
        IF(IPERM.EQ.0) THEN
          WRITE(27,1500) XB(J),UMEAN(J),USTD(J)
        ELSE
          WRITE(27,1500) XB(J),UMEAN(J),USTD(J),UPMEAN(J)
        ENDIF
        IF(IANGLE.EQ.1) WRITE(28,1500) XB(J),STHETA(J),VMEAN(J),VSTD(J)
        IF(IROLL.EQ.1) WRITE(29,1500) XB(J),RQ(J)
        IF(IPROFL.EQ.1) WRITE(30,1500) XB(J),PB(J),PS(J),VS(J)
        IF(IPERM.EQ.1) WRITE(31,1500) XB(J),UPSTD(J),DPSTA(J)
 160  CONTINUE
      IF(IOVER.EQ.1) THEN
C     Store mean values over wet duration
        IF(JDRY.GE.JR.AND.IOVER.EQ.1) THEN
        DO 170 J=(JR+1),JDUM
          DUM=H(J)+ZB(J,L)
          IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
            IF(JX2.LT.JMAX(L)) THEN
              IF(JXW.LE.J.AND.J.LE.JX2) THEN
                DUM=H(J)+ZW
                PWET(J)=1.D0
              ENDIF
            ENDIF
          ENDIF
          WRITE(22,1500) XB(J),DUM,H(J),SIGMA(J)
          IF(IPERM.EQ.0) THEN
            WRITE(27,1500) XB(J),UMEAN(J),USTD(J)
          ELSE
            WRITE(27,1500) XB(J),UMEAN(J),USTD(J),UPMEAN(J)
          ENDIF
          IF(IANGLE.EQ.1) WRITE(28,1500) XB(J),STHETA(J),
     +        VMEAN(J),VSTD(J)
          IF(IPROFL.EQ.1) WRITE(30,1500) XB(J),PB(J),PS(J),VS(J)
 170    CONTINUE
        ENDIF
C     Where UPMEAN, PB, PS, VS, and QP include effects of PWET.
        DO 171 J=1,JDUM
          IF(IPERM.EQ.0) WRITE(34,1500) XB(J),PWET(J)
          IF(IPERM.EQ.1) WRITE(34,1500) XB(J),PWET(J),QP(J)
 171    CONTINUE
        DO 161 J=JWD,JDRY
          WRITE(35,1500) XB(J),HEWD(J),UEWD(J),QEWD(J)
 161    CONTINUE
      ENDIF
C     
      IF(IPROFL.EQ.1) THEN
C     Smooth computed QBX(J) and QSX(J) before storing and plotting
        JMAXL=JMAX(L)
        DUMVEC = QBX
        CALL SMOOTH(JMAXL,DUMVEC,QBX)
        DUMVEC = QSX
        CALL SMOOTH(JMAXL,DUMVEC,QSX)
C     Smooth computed QBY(J) and QSY(J) if IANGLE=1
        IF(IANGLE.EQ.1) THEN
          DUMVEC=QBY
          CALL SMOOTH(JMAXL,DUMVEC,QBY)
          DUMVEC=QSY
          CALL SMOOTH(JMAXL,DUMVEC,QSY)
        ENDIF
        DO 162 J=1,JMAX(L)
          WRITE(32,1500) XB(J),QBX(J),QSX(J),(QBX(J)+QSX(J))
          IF(IANGLE.EQ.1) WRITE(33,1500) XB(J),QBY(J),QSY(J),
     +       (QBY(J) + QSY(J))
 162    CONTINUE
C     Store sediment transport volume per unit width
C     during Time=0.0 to Time=TIMOUT
        DO 163 J=1,JMAX(L)
          WRITE(37,1500) XB(J),VBX(J,L),VSX(J,L),(VBX(J,L)+VSX(J,L))
          IF(IANGLE.EQ.1) WRITE(38,1500) XB(J),VBY(J,L),VSY(J,L),
     +        (VBY(J,L)+VSY(J,L)) 
 163    CONTINUE
      ENDIF
C     
C     If IPROFL=2 or ICLAY=1, the following variables related to dike  
C     erosion at node J and line L are computed in Subr.22 EROSON
C     EDIKE(J,L)=downward erosion depth (m) from initial (time=0.0)
C     dike surface at time=TIMOUT for IPROFL=2 and ICLAY=0
C     EPCLAY(J)=downward clay erosion depth (m) for ICLAY=1 and IPROFL=1
C     DSTA(J)=variable (m*m/s) related to energy dissipation and dike or clay
C     erosion forcing at Time=TIMOUT
C     DSUM(J)= cumulative forcing (m*m) obtained by integrating DSTA(J)
C     from Time=0.0 to Time=TIMOUT
      IF(IPROFL.EQ.2.OR.ICLAY.EQ.1) THEN
        DO 164 J=1,JMAX(L)
          IF(IPROFL.EQ.2) THEN
              WRITE(39,1500) XB(J),EDIKE(J,L),DSTA(J),DSUM(J)
          ELSE    
              WRITE(39,1500) XB(J),EPCLAY(J,L),DSTA(J)
          ENDIF    
 164    CONTINUE
      ENDIF
C     
 1500 FORMAT(4F17.9)
C     
 200  CONTINUE
      RETURN
      END
C     
C     -08-----------------  END OF SUBROUTINE OUTPUT  ---------------------
C     #09#####################  SUBROUTINE POFLOW  ########################
C     
C     This subroutine computes mean and standard deviation of 
C     porous flow velocity and wave energy dissipation rate
C     DPSTA for given PKHSIG and DEDX at node J in the wet zone
C     
      SUBROUTINE POFLOW(J, L, PKHSIG, DEDX)
C     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NN=5000, NL=100)
C     
      COMMON /PERIOD/ TP,WKPO,ANGLE,WT(NN)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),
     +   CTHETA(NN),FSX,FSY,FE,QWX,QWY
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
C     
C     For porous layer thickness HP(J,L)=0.0, no velocity and dissipation
      IF(HP(J,L).EQ.0.D0) THEN
        UPMEAN(J) = 0.D0
        UPSTD(J) = 0.D0
        DPSTA(J) = 0.D0
      ENDIF
C     
      IF(HP(J,L).GT.0.D0) THEN
        A = 1.9D0*BESTA1
        B2 = BESTA2/WT(J)
        B = ALSTA + 1.9D0*B2
        UPSTD(J) = 0.5D0*(DSQRT(B*B+4.D0*A*PKHSIG)-B)/A
        A = SQRG1*(B2+BESTA1*UPSTD(J))
        C = CTHETA(J)*CTHETA(J)
        UPMEAN(J) = -DEDX/(ALSTA+A*(1.D0+C))
C     
C     To reduce numerical oscillations of UPMEAN(J), adjust
        RATIO = UPMEAN(J)/UPSTD(J)
        IF(RATIO.GT.0.5D0) UPMEAN(J)=0.5D0*UPSTD(J)
        IF(RATIO.LT.-0.5D0) UPMEAN(J)=-0.5D0*UPSTD(J)
        QP(J)=UPMEAN(J)*HP(J,L)
C     
        A2 = UPMEAN(J)*UPMEAN(J)
        B2 = UPSTD(J)*UPSTD(J)
        DPSTA(J) = HP(J,L)*(ALSTA*(A2+B2)+A*(2.D0*B2+A2*(1.D0+2.D0*C)))
C     
      ENDIF
C     
      RETURN
      END
C     ------------------  END OF SUBROUTINE POFLOW  ---------------------
C     #10#####################  SUBROUTINE QORATE  ########################
C     
C     This subroutine computes overtopping, overflow and seepage rates 
C     
      SUBROUTINE QORATE(ITIME,L,ITEQO,ICONV,ICALL)
C     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NN=5000, NB=30000,NL=100)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /SEAWAV/	TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     +   SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /BINPUT/ XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),
     +   YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX,FSY,FE, QWX, QWY
      COMMON /WBREAK/ GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /VELOCY/ UMEAN(NN),USTD(NN),USTA(NN),VMEAN(NN),VSTD(NN),  
     +   VSTA(NN)
      COMMON /ROLLER/ RBZERO, RBETA(NN), RQ(NN), RX(NN), RY(NN), RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/	RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
      DIMENSION       WSET(NN),ZRW(NN),SIGPT(NN)
C     
C     Find overtopping, overflow and seepage rates during ICALL=0
C     ...................... Start of ICALL=0 ...........................
      IF(ICALL.EQ.0) THEN
C     
C     Predict combined wave overtopping and overflow rate QOTF
C     by calling Subr.16 WETDRY for wet and dry zone
        CALL WETDRY(ITIME,L,ITEQO)
C     
C     Compute new combined rate QONEW and check convergency (one percent)
C     Allowable error is increased for QONEW less than 1.D-2(m*m/s)
        QONEW = QOTF
        IF(IPERM.EQ.1) QONEW=QONEW+SPRATE
        IF(QONEW.LT.1.D-5) THEN
C         IF(QO(L).LT.1.D-5) THEN
            ICONV = 0
            QO(L)=QONEW
            GOTO 99
C         ELSE
C           GOTO 98
C         ENDIF
        ENDIF
        DUM = DABS(QONEW-QO(L))/QONEW
        AER=1.D-4/QONEW
        IF(AER.LT.1.D-2) AER=1.D-2
        IF(DUM.LT.AER) THEN
          ICONV = 0
          QO(L)=QONEW
          GOTO 99
        ENDIF
 98     ICONV = 1
C     To avoid numerical oscillation and accelerate convergance
C     use FRACTN of previous value and (1.0-FRACTN) of new value
        FRACTN = 0.5D0 + 0.1D0*ITEQO
        IF(FRACTN.GT.0.9D0) FRACTN=0.9D0
        IF(ITEQO.EQ.10) FRACTN=0.5D0
        IF(ITEQO.EQ.20) FRACTN=0.D0
        QO(L) = FRACTN*QO(L) + (1.D0-FRACTN)*QONEW
        IF(IWTRAN.EQ.1) THEN
          DUM=GRAV*SIGMA(1)*SIGMA(1)*CTHETA(1)/CP(1)
          IF(QO(L).GT.DUM) QO(L)=DUM
        ENDIF
C        
 99     CONTINUE
        IF(IPOND.EQ.1) THEN
          IF(NOPOND.EQ.1) QM=QO(L)
          IF(JCREST(L).EQ.JXW) QM=QO(L)
          IF(ZW.GE.ZB(JMAX(L),L)) QM=QO(L)
        ENDIF
      ENDIF
C.....................................End of ICALL=0........................
C     
C....................................Start of ICALL=1.......................
C     Output computed values in file 20 'ODOC' if ICALL=1 in Subr.8 OUTPUT
      IF(ICALL.EQ.1) THEN
C     
C     Mean (ERMEAN) above datum Z=0 and standard deviation (SIGRUN) of runup
C     WSET(J)=wave setup above datum Z=0.0 during wet+dry duration
C     SIGPT(J)=standard deviation during wet+dry duration
C     ZRW(J)=runup wire elevation (RWH above bottom) at node J above Z=0.0 
C     
        DO 170 J=1,JCREST(L)
          IF(J.LE.JDRY) THEN
            WSET(J)= H(J)*PWET(J)+ZB(J,L)
          ELSE
            WSET(J)=ZB(J,L)
          ENDIF
          SIGPT(J)=SIGMA(J)*PWET(J)
          ZRW(J)=ZB(J,L)+RWH
 170    CONTINUE
C     
C     K=1,2 and 3 correspond to intersections of ZRW with WSET, (WSET-SIGPT)
C     and (WSET+SIGPT), respectively
        DO 100 K=1,3
          J=JDRY
          IF(JDRY.GT.JCREST(L)) J=JCREST(L)
          DUM1=ZRW(J)-WSET(J)
          IF(K.EQ.2) DUM1=DUM1+SIGPT(J)
          IF(K.EQ.3) DUM1=DUM1-SIGPT(J)
          IF(DUM1.LT.0.D0) THEN
            IF(K.EQ.1) THEN
              ETARUN=WSET(J)
              GOTO 100
            ENDIF
            IF(K.EQ.2) THEN
              Z1RUN=WSET(J)-SIGPT(J)
              X1RUN=XB(J)  !bdj
              GOTO 100
            ENDIF
            IF(K.EQ.3) THEN
              Z2RUN=WSET(J)+SIGPT(J)
              X2RUN=XB(J)  !bdj
              IF(X2RUN.LE.X1RUN)X2RUN=X1RUN+DX  !bdj
              GOTO 100
            ENDIF
          ENDIF
 105      J=J-1
          DUM2=ZRW(J)-WSET(J)
          IF(K.EQ.2) DUM2=DUM2+SIGPT(J)
          IF(K.EQ.3) DUM2=DUM2-SIGPT(J)
          IF(DUM2.GT.0.D0) THEN
            DUM1=DUM2
            GOTO 105
          ELSE
            DUM3=DUM1-DUM2
            DUMJ1=-DUM2/DUM3
            DUMJ=DUM1/DUM3
            DUMETA=DUMJ*WSET(J)+DUMJ1*WSET(J+1)
            IF(K.EQ.1) ETARUN=DUMETA
            IF(K.EQ.2) THEN 
              Z1RUN=DUMETA-DUMJ*SIGPT(J)-DUMJ1*SIGPT(J+1)
              X1RUN=DUMJ*XB(J)+DUMJ1*XB(J+1)
            ENDIF
            IF(K.EQ.3) THEN
              Z2RUN=DUMETA+DUMJ*SIGPT(J)+DUMJ1*SIGPT(J+1)
              X2RUN=DUMJ*XB(J)+DUMJ1*XB(J+1)
C BDJ 2011->2014 on 2014-10-02 
              IF((WSET(J+1)-WSET(J))/SIGPT(J).GT.10.*DX) THEN
                DUMETA=WSET(J)  !bdj
                Z2RUN=DUMETA+SIGPT(J) !bdj
                X2RUN=XB(J)     !bdj
                IF(x2run-x1run.le.01D0*DX) THEN
                  Z2RUN=Z1RUN  + .01D0*DX*BSLOPE(J,L)
                  X2RUN=X1RUN + .01D0*DX
                ENDIF
              ENDIF
C end BDJ 2011->2014 on 2014-10-02 
            ENDIF
          ENDIF
 100    CONTINUE
        SIGRUN=(Z2RUN-Z1RUN)/2.D0
        ERMEAN=(Z1RUN+ETARUN+Z2RUN)/3.D0
        SLPRUN=(Z2RUN-Z1RUN)/(X2RUN-X1RUN)

C bdj 2015-03-11  added catch for negative slopes
        SIGRUN=max(0.D0,SIGRUN)
        ERMEAN=max(z1run,ERMEAN)
        SLPRUN=max(0.D0,SLPRUN)
C end bdj 2015-03-11 added catch for negative slopes
C bdj 2015-07-06  added catch for cases where waves are very small
        IF(JR.LT.NINT(JSWL(L)/2.)) THEN
           SIGRUN=0.D0
           ERMEAN=SWLBC(ITIME)
           SLPRUN=0.D0
        ENDIF
C end bdj 2015-07-06  added catch for cases where waves are very small
C
C     R13=significant runup height above Z=0.0
C     R2P=two percent runup height above Z=0.0
C     RKAPPA=Kappa for runup probability distribution
        IF(IPERM.EQ.1) THEN
          R13=ERMEAN+(2.D0+SLPRUN)*SIGRUN
          RSC=(RCREST(L)-ERMEAN)/(R13-ERMEAN)
          RKAPPA=2.0D0+0.5D0/RSC**3.D0
        ELSE
          DUM=4.D0*SLPRUN
C BDJ 2011->2014 on 2014-10-02 
C          IF(DUM.GT.2.D0) DUM=2.D0
          IF(DUM.GT.1.D0) DUM=1.D0
C end BDJ 2011->2014 on 2014-10-02 
          R13=(ERMEAN-SWLBC(ITIME)+2.D0*SIGRUN)*(1.D0+DUM)+SWLBC(ITIME)
          RKAPPA=2.0D0
        ENDIF
        IF(RCREST(L).GT.ERMEAN) THEN
          R2P=ERMEAN+(R13-ERMEAN)*1.4D0**(2.D0/RKAPPA)
          R1P=ERMEAN+(R13-ERMEAN)*1.52D0**(2.D0/RKAPPA)
        ELSE
          RKAPPA=1000.D0
          R2P=R13
        ENDIF
C     
C.....Output swash hydrodynamics computed in Subr.16 WETDRY..........
        IF(JDRY.GE.JCREST(L)) THEN
          POTF=(DTANH(5.D0*PWET(JCREST(L))))**0.8D0
        ELSE
          POTF=0.D0
        ENDIF
C     Depth H, velocity U and discharge Q corresponding to exceedance 
C     probability EWD specified in Subr.04 PARAM
        IF(JWD.LE.JDRY) THEN
          DO 300 J=JWD, JDRY
            DUM = PWET(J)/EWD
            IF(DUM.LT.1.1D0) DUM=1.1D0
            HEWD(J)=(H(J)/PWET(J))*DLOG(DUM)
            DUM=USWD(J)
            IF(DUM.LT.0.D0) DUM=0.D0
            UEWD(J) = AWD*DSQRT(GRAV*HEWD(J))+DUM
            QEWD(J) = HEWD(J)*UEWD(J)
 300      CONTINUE
        ENDIF
C     Where computed HEWD(J), UEWD(J) and QEWD(J) are stored in Subr.8 OUTPUT
C     
        WRITE(20,920) SWLBC(ITIME),L,RCREST(L),JSWL(L),JWD,H1,JDRY,POTF,
     +     (QO(L)-SPRATE), SPRATE, QO(L), ITEQO
 920    FORMAT('COMBINED WAVE OVERTOPPING AND OVERFLOW'/
     +     'Still water level above Z=0 (m)              SWL=  ',F13.6/
     +     'Cross-shore line number                        L=  ',I3/
     +     'Structure or dune creat elevation (m)     RCREST=  ',F13.6/
     +     'Node number at SWL                          JSWL=  ',I6/
     +     'Wet and dry transition node                  JWD=  ',I6/
     +     'Mean water depth H1(m) at node JWD            H1=  ',F13.6/
     +     'End node of wet and dry zone                JDRY=  ',I6/
     +     'Wave overtopping probability at JCREST      POTF=  ',F13.6/
     +     'Comb. overtopping and overflow rate(m*m/s)  QOTF=  ',F13.9/
     +     'Seepage rate(m*m/s) at JCREST                 QP=  ',F13.9/
     +     'Total rate (QOTF+QP)(m*m/s)                     =  ',F13.9/
     +     'QO iteration number                        ITEQO=  ',I3/)
C     
C.........................Output empirical runup.......................
        WRITE(20,900) L,SLPRUN,ERMEAN,SIGRUN,R13,R2P,R1P
C     
 900    FORMAT('EMPIRICAL WAVE RUNUP'/
     +     'Cross-shore line number                        L=  ',I3/
     +     'Swash zone bottome slope                  SLPRUN=  ',F13.6/
     +     'Mean runup elevation above Z=0 (m)        ERMEAN=  ',F13.6/
     +     'Runup standard deviation (m)              SIGRUN=  ',F13.6/
     +     'Significant runup height above Z=0 (m)       R13=  ',F13.6/
     +     '2 percent runup height above Z=0 (m)         R2P=  ',F13.6/
     +     '1 Percent runup height above z=0 (m)         R1P=  ',F13.6/)
C     
        IF(IWTRAN.EQ.1) THEN 
C         IF(JDRY.EQ.JSL1.AND.JSL.LT.JMAX(L)) THEN
           WRITE(20,940)L,JSL,XB(JSL),WSETUP(JSL),SIGMA(JSL),XB(JMAX(L))
     +        ,WSETUP(JMAX(L)),SIGMA(JMAX(L)),(SIGMA(JMAX(L))/SIGMA(1))
C         ELSE
C           WRITE(20,941) JDRY,JSL,JSL1
C         ENDIF
        ENDIF
 940    FORMAT('WAVE TRANSMISSION DUE TO OVERTOPPING'/
     +     'Cross-shore line number                        L=  ',I3/
     +     'Starting node for wave transmission          JSL=  ',I6/
     +     'X-coordinate (m) at node JSL                  XB=  ',F13.6/
     +     'Wave setup (m) above SWL at node JSL      WSETUP=  ',F13.6/
     +     'Standard deviation (m) at node JSL         SIGMA=  ',F13.6/
     +     'X-coordinate (m) at landward end node JMAX      =  ',F13.6/
     +     'Wave setup (m) above SWL at landward end node JMAX= ',F13.6/
     +     'Standard dev. (m) at landward end node JMAX     =  ',F13.6/
     +     'Wave transmission coefficient at JMAX           =  ',F13.6/)
 941    FORMAT(/'IWTRAN=1 BUT NO WAVE TRANSMISSION'/'JDRY=',I6,
     +     '  and JSL=',I6,'  and JSL1=',I6/' because entire structure',
     +     ' is submerged or no wet zone exists landward of crest'/) 
C
        IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
          WRITE(20,960) L,JCREST(L),JXW,JX2,ZW,QD,QM
        ENDIF
 960    FORMAT('PONDED WATER IN RUNNEL'/
     +     'Cross-shore line number                        L=  ',I3/
     +     'Ridge crest node                          JCREST=  ',I6/
     +     'Ponded water nodes from                      JXW=  ',I6/
     +     '                     to                      JX2=  ',I6/
     +     'Ponded water level (m)		                ZW=  ',F13.6/
     +     'Wave-induced volume flux (m*m/s)              QD=  ',F13.6/
     +     'Wave overtopping rate (m*m/s) at JMAX         QM=  ',F13.6/)
C     
C.................................End of ICALL=1...........................
C     
      ENDIF
      RETURN
      END
C     -10-----------------  END OF SUBROUTINE QORATE  ---------------------
C     #11#####################  SUBROUTINE SEDTRA  ########################
C     
C     This subr. calculates cross-shore and longshore sediment transport 
C     
      SUBROUTINE SEDTRA(L)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NB=30000,NL=100)
      DIMENSION QRAW(NN),GSLRAW(NN),ASLRAW(NN),ASLOPE(NN),RS(NN),RB(NN),
     +   PBWD(NN),PSWD(NN),VSWD(NN),QSXWD(NN),QBXWD(NN),QRAWD(NN),
     +   HDIP(NN),QSYWD(NN),QBYWD(NN)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /PERIOD/ TP,WKPO,ANGLE,WT(NN)
      COMMON /PREDIC/ HRMS(NN), SIGMA(NN), H(NN), WSETUP(NN), SIGSTA(NN)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /WBREAK/ GAMMA, QBREAK(NN), DBSTA(NN),SISMAX, ABREAK(NN)
      COMMON /ENERGY/ EFSTA(NN), DFSTA(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /CONSTA/ GRAV, SQR2, SQR8, PI, TWOPI, SQRG1, SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX,FSY,FE,QWX,QWY
      COMMON /VELOCY/ UMEAN(NN),USTD(NN),USTA(NN), VMEAN(NN), VSTD(NN),
     +   VSTA(NN)
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +   EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /SEDOUT/ PS(NN), VS(NN), QSX(NN), QSY(NN), 
     +   PB(NN),GSLOPE(NN),QBX(NN),QBY(NN),Q(NN)
      COMMON /PROCOM/ DELT,DELZB(NN,NL)
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /COMPAR/ HWDMIN,NPT,NPE
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
      COMMON /VEGETA/ VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
C     
C     GSLMAX = Maximum absolute value of GSLOPE function
      DATA GSLMAX,BQCOEFF  /10.D0,8.D0/
C     
C.....Cross-Shore and Longshore Sediment Transport at Node J
C     RB(J) = Sediment movement initiation parameter
C     RS(J) = Sediment suspension initiation parameter
C     PB(J) = bedload probability
C     PS(J) = suspended load probability
C     VS(J) = suspended sediment volume per unit area (m)	
C     GSLOPE(J) = bed slope correction for QBX(J)
C     ASLOPE(J) = bed slope correction for suspended load parameter SLP
C     QBX(J)= Cross-shore  bedload transport rate per unit width (m*m/s)
C     QBY(J)= Longshore  bedload transport rate per unit width (m*m/s)
C     BRF   = Bedload reduction factor for hard bottom (ISEDAV=1)
C     QSX(J)= Cross-shore suspended sediment transport rate (m*m/s)
C     QSY(J)= Longshore suspended sediment transport rate (m*m/s)
C     Q(J)  = total cross-shore sedimet transport rate including void
C     (m*m/s) used for beach profile change computation 
C     
      IF(TIME.EQ.0.D0) THEN
        BSLOP1 = -TANPHI*(GSLMAX-1.D0)/GSLMAX
        BSLOP2 =  TANPHI*(GSLMAX+1.D0)/(GSLMAX+2.D0)
        IGMILD=0
        IF(IPERM.EQ.1.AND.JSWL(L).LT.JMAX(L)) THEN
          DUM=0.5D0*TANPHI
          IF(BSLOPE(JSWL(L),L).LT.DUM.AND.IWTRAN.EQ.0) IGMILD=1
          IF(HP(JSWL(L),L).EQ.0.D0) IGMILD=0
        ENDIF
C     where input bedload parameter is increased in the surf zone in the
C     following if IGMILD=1 (based on two gravel tests MH and MB only)
      ENDIF
      IF(IVWALL(L).EQ.2) THEN
        JDUM=JMAX(L)
        BSLOPE(JDUM,L)=0.D0
        BSLOPE(JDUM-1,L)=(ZB(JDUM-1,L)-ZB(JDUM-2,L))/DX
      ENDIF
C     
      DO 100 J=1,JMAX(L)
        IF(BSLOPE(J,L).LT.0.D0) THEN
          IF(BSLOPE(J,L).GT.BSLOP1) THEN
            GSLRAW(J) = TANPHI/(TANPHI + BSLOPE(J,L))
          ELSE
            GSLRAW(J) = GSLMAX
          ENDIF
        ELSE
          IF(BSLOPE(J,L).LT.BSLOP2) THEN
            GSLRAW(J) = (TANPHI - 2.D0*BSLOPE(J,L))/(TANPHI-BSLOPE(J,L))
          ELSE
            GSLRAW(J) = -GSLMAX
          ENDIF
          IF(IGMILD.EQ.1) THEN
            IF(GSLRAW(J).LT.0.D0) GSLRAW(J)=0.D0
          ENDIF
        ENDIF
        ASLRAW(J) = SLP
C     Add bottom slope effect to suspended load parameter
      IF(BSLOPE(J,L).GT.0.D0) ASLRAW(J)=SLP+(BSLOPE(J,L)/TANPHI)
     +  **0.5D0
C     
C     Vegetation effect is included in DFSTA(J) with J=1,2,...,JR
C     for energy dissipation rate due to bottom friction and vegetation
C     Assume no vegetation effect on sediment transport if IVEG=2
C     IF(J.LE.JR.AND.IVEG.GE.1) THEN
C       DUM=VEGH(J,L)
C       IF(DUM.GT.H(J)) DUM=H(J)
C       VEGCV=1.D0+DUM*VEGFB(J,L)
C       DFSTA(J)=DFSTA(J)/VEGCV
C     ENDIF
C     Vegetation effect is removed from DFSTA(J) above
C     
 100  CONTINUE
C     
C     Smoothing GSLOPE before Q is computed in Subr.14 SMOOTH
      JMAXL=JMAX(L)
      CALL SMOOTH(JMAXL, GSLRAW, GSLOPE)
      CALL SMOOTH(JMAXL, ASLRAW, ASLOPE)
C     
C     Sediment transport rates are computed for normally incident waves
C     in wet zone (IANGLE=0); wet and dry zone (IOVER=1) for IANGLE=0
C     and 1; and obliquelly incident waves in wet zone (IANGLE=1)
      IF(IANGLE.EQ.1) GOTO 888
C     
C.....Normally Incident Waves in wet zone...............................
      IF(IANGLE.EQ.0) THEN
        DO 110 J = 1,JR
          IF(D50.LT.CSEDIA) THEN
            RB(J) = DSQRT(GSD50S/FB2(J,L))/USTD(J)
          ELSE
            RB(J)=GSD50S/USTD(J)
          ENDIF
          RS(J) = WF/USTD(J)/FB2(J,L)**0.3333D0
          US = USTA(J)
          PB(J)=0.5D0*(ERFCC((RB(J)+US)/SQR2 )+ERFCC((RB(J)-US)/SQR2))
          PS(J)=0.5D0*(ERFCC((RS(J)+US)/SQR2 )+ERFCC((RS(J)-US)/SQR2))
          IF(PS(J).GT.PB(J)) PS(J) = PB(J)
          IF(IROLL.EQ.0) THEN
            VS(J) = PS(J)*(EFFF*DFSTA(J) + EFFB*DBSTA(J))/WFSGM1
          ELSE
            VS(J) = PS(J)*(EFFF*DFSTA(J) + EFFB*RBETA(J)*RQ(J))/WFSGM1
          ENDIF
          VS(J) = VS(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
          BQ=BLD
C     Input bedload parameter in Subr.2 INPUT is adjusted to account
C     for QBREAK=fraction(0.0 to 1.0) of breaking waves.
C     Put "C" in front of the next line for no adjustment
          IF(D50.LT.CSEDIA)BQ=BQ*(0.5D0+QBREAK(J))
          IF(IGMILD.EQ.1) THEN
            BQ=BLD*(1.D0+BQCOEFF*QBREAK(J))
          ENDIF
C     BDJ added 2012-10-23 
            DECAYL = MIN(XB(JSWL(L))/4.D0,2.D0*TP*CP(1)) ! The decay length
            JDECAY = NINT(DECAYL/DX)! index of decay intrusion length
C     end BDJ added 2012-10-23             
          QBX(J) = BQ*PB(J)*GSLOPE(J)*USTD(J)**3
          IF(ISEDAV.GE.1.OR.ISTSAN.EQ.1) THEN
            IF(ISEDAV.GE.1) THEN
               DUM=HP(J,L)
               IF(ISEDAV.EQ.2) THEN
                  DUM=ZB(J,L)-ZMESH(J,L)
                  IF(DUM.LT.0.D0) DUM=0.D0
               ENDIF
               IF(DUM.GE.D50) THEN
                  BRF=1.D0
               ELSE
                  BRF=(DUM/D50)**BEDLM
               ENDIF
            ELSE
               BRF=DEXP(-CPSTON*HP(J,L)/D50)
            ENDIF
            VS(J)=BRF*VS(J)
            QBX(J)=BRF*QBX(J)
          ENDIF
          QSX(J) = ASLOPE(J)*UMEAN(J)*VS(J)
C     Add onshore suspended sediment transport due to wave overtopping
          IF(IOVER.EQ.1) THEN
            DUM = H(J)
            IF(DUM.LT.HWDMIN) DUM = HWDMIN
            AO=SLPOT
            DUMQ=QO(L)
            QSX(J)=QSX(J)+AO*VS(J)*DUMQ/DUM
          ENDIF
          QRAW(J) = (QBX(J) + QSX(J))/SPORO1
 110    CONTINUE
C     
C     BDJ added on 2012-10-24         
        QSX(1:JDECAY) = QSX(JDECAY)
        QBX(1:JDECAY) = QBX(JDECAY)
        QRAW(1:JDECAY) = QRAW(JDECAY)
C     end BDJ added on 2012-10-24         
C
C     If IOVER=0 or JDRY.LE.JR, no wet and dry zone and use scarping formula
C     If IOVER=1, compute sediment transport in wet and dry zone
        IF(IOVER.EQ.0.OR.JDRY.LE.JR) THEN
C     
C     Linear extrapolation for scarped slope exceeding TANPHI
C     only if QRAW(JR) is offshore
          JR1 = JR+1
          JE = JR1
          IF(QRAW(JR).LT.0.D0) THEN
 102        IF(BSLOPE(JE,L).GT.TANPHI) THEN
              JE = JE+1
              IF(JE.GE.JMAX(L)) GOTO 103
              GOTO 102
            ENDIF
          ENDIF
 103      JD = JE-JR
          IF(JD.GE.2) THEN
            DO 104 J=JR1,JE-1
              DUM=DBLE(JE-J)/DBLE(JD)
              QRAW(J)=DUM*QRAW(JR)
 104        CONTINUE
          ENDIF
C     Subr.15 EXTRAPO, extrapolates for nodes from J1 to J2
          CALL EXTRAPO(JR1, JMAXL, QBX)
          CALL EXTRAPO(JR1, JMAXL, QSX)
          CALL EXTRAPO(JE, JMAXL, QRAW)
          CALL EXTRAPO(JR1, JMAXL, PB)
          CALL EXTRAPO(JR1, JMAXL, PS)
          CALL EXTRAPO(JR1, JMAXL, VS)
          GOTO 900
        ENDIF 
      ENDIF
C     End of IANGLE=0 in wet zone ..........................................
C     
C     ....... Wet and Dry Zone for IANGLE=0 and 1 ..........................
C     For node J=JWD to JDRY in wet and dry (WD) zone
C     PBWD(J)=bedload probability computed in Subr.18 PROBWD
C     PSWD(J)=suspended load probability computed in Subr.18 PROBWD
C     VSWD(J)=suspended sediment volume per unit area(m)
C     QSXWD(J)=cross-shore suspended sediment transport rate(m*m/s)
C     QBXWD(J)=cross-shore bedload sediment transport rate(m*m/s)
C     where hydrodynamic variables in WD zone are computed in Subr.16 WETDRY
C     HDIP(J)=mean water depth adjusted for dip in wet
C     and dry zone used for suspended sediment transport rate
C     if IVWALL(L)=0 (no vertical wall at landward end)
C     
 999  CONTINUE
      IF(IOVER.EQ.1.AND.JDRY.GT.JR) THEN
        IF(IVWALL(L).EQ.0) THEN
          J=JWD
          HDIP(J)=H(J)
          ZBPEAK=ZB(J,L)
 140      J=J+1
          IF(J.EQ.JDRY) GOTO 142
          IF(J.GT.JDRY) GOTO 145
          IF(ZB(J-1,L).LT.ZB(J,L).AND.ZB(J,L).GE.ZB(J+1,L)) 
     +       ZBPEAK=ZB(J,L)
          DUM=ZBPEAK-ZB(J,L)
          IF(DUM.GT.H(J)) THEN
            HDIP(J)=DUM
          ELSE
            HDIP(J)=H(J)
          ENDIF
          IF(J.LT.JCREST(L)) GOTO 140
 142      J=JDRY
          HDIP(J)=H(J)
          ZBPEAK=ZB(J,L)
 141      J=J-1
          IF(J.LE.JCREST(L)) GOTO 145
          IF(ZB(J-1,L).LT.ZB(J,L).AND.ZB(J,L).GE.ZB(J+1,L))
     +       ZBPEAK=ZB(J,L)
          DUM=ZBPEAK-ZB(J,L)
          IF(DUM.GT.H(J)) THEN
            HDIP(J)=DUM
          ELSE
            HDIP(J)=H(J)
          ENDIF
          GOTO 141
        ENDIF
 145    CONTINUE
C     For gravel tests MH and MB (IGMILD=1), landward extension
C     of bedload was necessary
        IF(IGMILD.EQ.1) JEXT=JWD+NINT(4.2D0*HRMS(1)/DX)
C     
        DO 150 J=JWD,JDRY
          IF(IPERM.EQ.0.AND.INFILT.EQ.0) THEN
            QWX=QO(L)
            IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
              IF(J.GE.JX2) QWX=QM
              IF(J.GT.JXW.AND.J.LT.JX2) THEN
                QWX=QO(L)-(QO(L)-QM)*(XB(J)-XB(JXW))/(XB(JX2)-XB(JXW))
              ENDIF
            ENDIF
          ELSE
            QWX=QO(L)-QP(J)
          ENDIF
          USWD(J)=QWX/H(J)-AQWD*DSQRT(GRAV*H(J)/PWET(J))
          IF(D50.LT.CSEDIA) THEN
            UCB=DSQRT(GSD50S/FB2(J,L))
          ELSE
            UCB=GSD50S
          ENDIF
          PWAGH=PWET(J)/AGWD/GRAV/H(J)
          CALL PROBWD(PWET(J),PWAGH,USWD(J),UCB,PBWD(J))
          UCS=WF/FB2(J,L)**0.333333D0
          CALL PROBWD(PWET(J),PWAGH,USWD(J),UCS,PSWD(J))
          IF(PSWD(J).GT.PBWD(J)) PSWD(J)=PBWD(J)
C     
C     Suspended load VBF and bedload factor BLDS in wet and dry zone  
C     are adjusted so that VS(J)=VSWD(J) and QBX(J)=QBXWD(J) at J=JWD 
          IF(J.EQ.JWD) THEN
            VBF=1.D0
            BLDS=1.D0
          ENDIF
          VSWD(J)=VBF*PSWD(J)
          VSWD(J)=VSWD(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
          QBXWD(J)=BLDS*PBWD(J)*GSLOPE(J)*USTD(J)**3.D0
          IF(J.EQ.JWD) THEN
            IF(VSWD(J).GT.1.D-20)THEN
              VBF=VS(J)/VSWD(J)
            ELSE
              VBF=0.D0
            ENDIF
            VSWD(J)=VS(J)
            IF(DABS(QBXWD(J)).GT.1.D-20)THEN
              BLDS=QBX(J)/QBXWD(J)
            ELSE
              BLDS=BLD
            ENDIF
            QBXWD(J)=QBX(J)
          ENDIF
          IF(ISEDAV.GE.1.OR.ISTSAN.EQ.1) THEN
            IF(ISEDAV.GE.1) THEN
              DUM=HP(J,L)
              IF(ISEDAV.EQ.2) THEN
                  DUM=ZB(J,L)-ZMESH(J,L)
                  IF(DUM.LT.0.D0) DUM=0.D0
              ENDIF
              IF(DUM.GE.D50) THEN
                  BRF=1.D0
              ELSE
                  BRF=(DUM/D50)**BEDLM
              ENDIF
            ELSE
              BRF=DEXP(-CPSTON*HP(J,L)/D50)
            ENDIF
            IF(IVWALL(L).EQ.0) VSWD(J)=BRF*VSWD(J)
            QBXWD(J)=BRF*QBXWD(J)
          ENDIF
          QSXWD(J)=ASLOPE(J)*UMEAN(J)*VSWD(J)
          IF(IOVER.EQ.1) THEN
            DUM = H(J)
            IF(IVWALL(L).EQ.0) DUM=HDIP(J)
            IF(DUM.LT.HWDMIN) DUM = HWDMIN
            AO=SLPOT
            DUMQ=QO(L)
            IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
              IF(J.GE.JCREST(L).AND.J.LT.JX2) DUMQ=QD
              IF(J.GE.JX2) DUMQ=QM
            ENDIF
            QSXWD(J)=QSXWD(J)+AO*VSWD(J)*DUMQ/DUM
          ENDIF
C     
C     If IGMILD=1, adjust QBXWD as follows 
            IF(IGMILD.EQ.1) THEN
              IF(J.LE.JEXT) QBXWD(J)=QBXWD(JWD)
              IF(J.EQ.JWD) THEN
                JWD1=JWD+1
                IF(JWD1.LT.JR) THEN
                  DO 149 JJ=JWD1,JR
                    QBX(JJ)=QBX(JWD)
 149              CONTINUE
                ENDIF
              ENDIF
            ENDIF
            QRAWD(J)=(QBXWD(J)+QSXWD(J))/SPORO1
            IF(IANGLE.EQ.1) THEN
              US=UMEAN(J)/USTD(J)
              DUM=(1.D0+US*US)*VMEAN(J)/USTD(J)+2.D0*US*STHETA(J)
              QBYWD(J)=QBXWD(J)*DUM/GSLOPE(J)
              QSYWD(J)=VMEAN(J)*VSWD(J)
            ENDIF
 150      CONTINUE
C     
C     If IPOND=1 and NOPOND=0, ponded water exists between 
C     nodes J=JXW and JX2. Ponded water is assumed to cause
C     sedimentation where WF=sediment fall velocity and 
C     QD=wave-induced onshore volume flux at ridge crest
C     node JCREST for deposition
          IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
            JDUM=JDRY
            DLEN=(XB(JX2)-XB(JXW))/SLPOT
C           DUM=QD/WF
C           IF(DLEN.LT.DUM) DLEN=DUM
            IF(JDUM.GT.JXW) THEN
              JXW1=JXW+1
              DO 151 J=JXW1, JDUM
                DUM=DEXP(-(XB(J)-XB(JXW))/DLEN)
                PBWD(J)=PBWD(J)*DUM
                VSWD(J)=VSWD(J)*DUM
                PSWD(J)=PSWD(J)*DUM
                QBXWD(J)=QBXWD(J)*DUM
                QSXWD(J)=QSXWD(J)*DUM
                QRAWD(J)=(QBXWD(J)+QSXWD(J))/SPORO1
                IF(IANGLE.EQ.1) THEN
                  QBYWD(J)=QBYWD(J)*DUM
                  QSYWD(J)=QSYWD(J)*DUM
                ENDIF
 151          CONTINUE
            ENDIF
          ENDIF
C     
C     Connect wet variables (J=1 to JR) with WD variables
C     (J=JWD to JDRY) using Subr.17 TRANWD
          IF(JDRY.GT.JR) THEN
            CALL TRANWD(PB,JR,PBWD,JWD,JDRY)
            CALL TRANWD(PS,JR,PSWD,JWD,JDRY)
            CALL TRANWD(VS,JR,VSWD,JWD,JDRY)
            CALL TRANWD(QSX,JR,QSXWD,JWD,JDRY)
            CALL TRANWD(QBX,JR,QBXWD,JWD,JDRY)
            CALL TRANWD(QRAW,JR,QRAWD,JWD,JDRY)
            IF(IANGLE.EQ.1) THEN
              CALL TRANWD(QSY,JR,QSYWD,JWD,JDRY)
              CALL TRANWD(QBY,JR,QBYWD,JWD,JDRY)
            ENDIF
          ENDIF
C     
C     Compute sediment transport in landward wet zone of wave transmission
C     Suspended load VBF and bedload factor BQ are adjusted so that
C     VS(JSL)=VS(JSL1) and QBX(JSL)=QBX(JSL1)
        IF(IWTRAN.EQ.1.AND.JDRY.EQ.JSL1) THEN
        BQ=1.D0
        VBF=1.D0
        IF(JSL.GE.JMAX(L)) GOTO 165
        DO 160 J=JSL,JMAX(L)
          IF(D50.LT.CSEDIA) THEN
            RB(J) = DSQRT(GSD50S/FB2(J,L))/USTD(J)
          ELSE
            RB(J)=GSD50S/USTD(J)
          ENDIF
          RS(J) = WF/USTD(J)/FB2(J,L)**0.3333D0
          US = UMEAN(J)/USTD(J)
          PB(J)=0.5D0*(ERFCC((RB(J)+US)/SQR2 )+ERFCC((RB(J)-US)/SQR2))
          PS(J)=0.5D0*(ERFCC((RS(J)+US)/SQR2 )+ERFCC((RS(J)-US)/SQR2))
          IF(PS(J).GT.PB(J)) PS(J) = PB(J)
          VS(J) = PS(J)*VBF
          VS(J) = VS(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
          QBX(J) = BQ*PB(J)*GSLOPE(J)*USTD(J)**3
          IF(J.EQ.JSL) THEN
            IF(VS(J).GT.1.D-20)THEN
              VBF=VS(JSL1)/VS(J)
            ELSE
              VBF=0.D0
            ENDIF
            VS(J)=VS(JSL1)
            IF(DABS(QBX(J)).GT.1.D-20)THEN
              BQ=QBX(JSL1)/QBX(J)
            ELSE
              BQ=BLD
            ENDIF
            QBX(J)=QBX(JSL1)
          ENDIF
          IF(ISEDAV.GE.1.OR.ISTSAN.EQ.1)THEN
            IF(ISEDAV.GE.1) THEN
              DUM=HP(J,L)
              IF(ISEDAV.EQ.2) THEN
                  DUM=ZB(J,L)-ZMESH(J,L)
                  IF(DUM.LT.0.D0) DUM=0.D0
              ENDIF
              IF(DUM.GE.D50)THEN
                  BRF=1.D0
              ELSE
                  BRF=(DUM/D50)**BEDLM
              ENDIF
            ELSE    
              BRF=DEXP(-CPSTON*HP(J,L)/D50)
            ENDIF
            VS(J)=BRF*VS(J)
            QBX(J)=BRF*QBX(J)
          ENDIF
          QSX(J) = ASLOPE(J)*UMEAN(J)*VS(J)
          QRAW(J) = (QBX(J) + QSX(J))/SPORO1
          IF(IANGLE.EQ.1) THEN
            QBY(J)=0.D0
            QSY(J)=0.D0
          ENDIF
 160    CONTINUE
 165    CONTINUE
        ELSE
C     Connect QSX(J), QBX(J) and QRAW(J)=0.0 landward of JDRY for no wave transmission
C     (IWTRAN=0) or no wave overtopping to landward wet zone even if IWTRAN=1
          IF(JDRY.LT.JMAX(L)) THEN
            JDRY1=JDRY+1
            CALL EXTRAPO(JDRY1,JMAX(L),QSX)
            CALL EXTRAPO(JDRY1,JMAX(L),QBX)
            CALL EXTRAPO(JDRY1,JMAX(L),QRAW)
            CALL EXTRAPO(JDRY1,JMAX(L),PB)
            CALL EXTRAPO(JDRY1,JMAX(L),PS)
            CALL EXTRAPO(JDRY1,JMAX(L),VS)
            IF(IANGLE.EQ.1) THEN
              CALL EXTRAPO(JDRY1,JMAX(L),QSY)
              CALL EXTRAPO(JDRY1,JMAX(L),QBY)
            ENDIF
          ENDIF
        ENDIF
C     
        GOTO 900
C     
      ENDIF
C     End of Wet and Dry Zone for IANGLE=0 and 1 .......................
C     
C.....Obliquely Incident Waves in wet zone .............................
 888  CONTINUE
      IF(IANGLE.EQ.1) THEN
        DO 190 J=1,JR
          SIGT = USTD(J)/CTHETA(J)
          IF(D50.LT.CSEDIA) THEN
            RB(J)= DSQRT(GSD50S/FB2(J,L))/SIGT
          ELSE
            RB(J)=GSD50S/SIGT
          ENDIF
          RS(J)= WF/SIGT/FB2(J,L)**0.3333D0
          WSTA = USTA(J)*CTHETA(J) + VSTA(J)*STHETA(J)
          VCUS = VSTA(J)*CTHETA(J) -USTA(J)*STHETA(J)
          FS = RS(J)*RS(J) - VCUS*VCUS
          IF(FS.LT.0.D0) THEN
            PS(J) = 1.D0
          ELSE
            FS = DSQRT(FS)
            PS(J)= 0.5D0*(ERFCC((FS+WSTA)/SQR2)+ERFCC((FS-WSTA)/SQR2))
          ENDIF
          FB = RB(J)*RB(J) - VCUS*VCUS
          IF(FB.LT.0.D0) THEN
            PB(J) = 1.D0
          ELSE
            FB = DSQRT(FB)
            PB(J)= 0.5D0*(ERFCC((FB+WSTA)/SQR2)+ERFCC((FB-WSTA)/SQR2))
          ENDIF
          IF(PS(J).GT.PB(J)) PS(J)=PB(J)
          IF(IROLL.EQ.0) THEN
            VS(J) = PS(J)*(EFFF*DFSTA(J)+EFFB*DBSTA(J))/WFSGM1 
          ELSE
            VS(J) = PS(J)*(EFFF*DFSTA(J)+EFFB*RBETA(J)*RQ(J))/WFSGM1 
          ENDIF
          VS(J) = VS(J)*DSQRT(1.D0+BSLOPE(J,L)*BSLOPE(J,L))
          VSTA2 = VSTA(J)*VSTA(J)
          TWOS = 2.D0*STHETA(J)
          BQ=BLD
C     Input bedload parameter in Subr.2 INPUT is adjusted to account
C     for QBREAK=fraction(0.0 to 1.0) of breaking waves.
C     Put "C" in front of the next line for no adjustment 
          IF(D50.LT.CSEDIA)BQ=BQ*(0.5D0+QBREAK(J))
          IF(IGMILD.EQ.1) THEN
            BQ=BLD*(1.D0+BQCOEFF*QBREAK(J))
          ENDIF
C     BDJ added 2012-10-23 
            DECAYL = MIN(XB(JSWL(L))/4.D0,2.D0*TP*CP(1)) ! The decay length
            JDECAY = NINT(DECAYL/DX)! index of decay intrusion length
C     end BDJ added 2012-10-23             
          DUM = BQ*PB(J)*(USTD(J)*USTD(J) + VSTD(J)*VSTD(J))**1.5D0
          QBX(J) = DUM*GSLOPE(J)*(1.D0 + USTA(J)*VSTA2 + TWOS*VCUS)
          QBY(J) = DUM*(VSTA(J)*(1.D0 + USTA(J)*USTA(J)+ VSTA2)+
     +       TWOS*WSTA)
          IF(ISEDAV.GE.1.OR.ISTSAN.EQ.1) THEN
            IF(ISEDAV.GE.1) THEN
              DUM=HP(J,L)
              IF(ISEDAV.EQ.2) THEN
                   DUM=ZB(J,L)-ZMESH(J,L)
                   IF(DUM.LT.0.D0) DUM=0.D0
              ENDIF
              IF(DUM.GE.D50) THEN
                   BRF=1.D0
              ELSE
                   BRF=(DUM/D50)**BEDLM
              ENDIF
            ELSE
              BRF=DEXP(-CPSTON*HP(J,L)/D50)
            ENDIF
            VS(J)=BRF*VS(J)
            QBX(J)=BRF*QBX(J)
            QBY(J)=BRF*QBY(J)
          ENDIF
          QSX(J) = ASLOPE(J)*UMEAN(J)*VS(J)
          IF(IOVER.EQ.1) THEN
            DUM = H(J)
            IF(DUM.LT.HWDMIN) DUM = HWDMIN
            AO=SLPOT
            DUMQ=QO(L)
            QSX(J)=QSX(J)+AO*VS(J)*DUMQ/DUM
          ENDIF
          QSY(J) = VMEAN(J)*VS(J)
          QRAW(J) = (QBX(J) + QSX(J))/SPORO1
 190      CONTINUE
C     
C     Scarping extrapolation is included for oblique waves as well
        IF(IOVER.EQ.0.OR.JDRY.LE.JR) THEN
          JR1 = JR+1
          JE = JR1
          IF(QRAW(JR).LT.0.D0) THEN
 202        IF(BSLOPE(JE,L).GT.TANPHI) THEN
              JE = JE+1
              IF(JE.GE.JMAX(L)) GOTO 203
              GOTO 202
            ENDIF
          ENDIF
 203      JD = JE-JR
          IF(JD.GE.2) THEN
            DO 204 J=JR1, JE-1
              DUM=DBLE(JE-J)/DBLE(JD)
              QRAW(J) =DUM*QRAW(JR)
 204        CONTINUE
          ENDIF
C     
C     Subr. 15 EXTRAPO extrapolates for nodes from J1 to J2
          CALL EXTRAPO(JR1, JMAXL, QBX)
          CALL EXTRAPO(JR1, JMAXL, QSX)
          CALL EXTRAPO(JR1, JMAXL, QBY)
          CALL EXTRAPO(JR1, JMAXL, QSY)
          CALL EXTRAPO(JE, JMAXL, QRAW)
          CALL EXTRAPO(JR1, JMAXL, PB)
          CALL EXTRAPO(JR1, JMAXL, PS)
          CALL EXTRAPO(JR1, JMAXL, VS)
          GOTO 900
        ELSE
          GOTO 999
        ENDIF
C     
      ENDIF
C     End of IANGLE=1 in wet zone ..........................................
C     
C     Adjust computed QSX(1) and QBX(1) at node 1 to be consistent with
C     the boundary condition used in Subr.12 CHANGE
  900 CONTINUE
      QSX(1)=QSX(2)
      QBX(1)=QBX(2)
      QRAW(1)=QRAW(2)
      IF(IANGLE.EQ.1)THEN
        QSY(1)=QSY(2)
        QBY(1)=QBY(2)
      ENDIF
C     Adjust sediment transport rates at node JMAX to be consitent with the
C     boundary condition used in subr.12 CHANGE
      JMAX1=JMAX(L)-1
      QSX(JMAXL)=QSX(JMAX1)
      QBX(JMAXL)=QBX(JMAX1)
      QRAW(JMAXL)=QRAW(JMAX1)
C     Smoothing QRAW (before DELZB is computed) using Sub.14 SMOOTH
      CALL SMOOTH(JMAXL,QRAW,Q)
C     
      RETURN
      END
C     
C     --11-----------------  END OF SUBROUTINE SEDTRA  ---------------------
C     #12##################### SUBROUTINE CHANGE ###########################
C     
C     Compute the bottom elevation change using the volume conservation
C     of bottom sediment.
C     
      SUBROUTINE CHANGE(ITIME,L,IEND,ICALL)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000, NB=30000,NL=100)
      DIMENSION DZBDT(NN),CB(NN),R(NN),DELZBRW(NN),DELZBJ(NN),V(NL),
     +   VDY(NL),AVY(NL),ADZX(NL)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     +   SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /BINPUT/ XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),
     +   YLINE(NL),DYLINE(NL),AGLINE(NL),NBINP(NL)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /SEDINP/ WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +   EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /SEDOUT/ PS(NN),VS(NN),QSX(NN),QSY(NN),
     +   PB(NN),GSLOPE(NN),QBX(NN),QBY(NN),Q(NN)
      COMMON /SEDVOL/ VBX(NN,NL),VSX(NN,NL),VBY(NN,NL),VSY(NN,NL),
     +   VY(NN,NL),DZX(NN,NL)
      COMMON /PROCOM/ DELT,DELZB(NN,NL)
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON /WIMESH/WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL)
      COMMON /STONES/ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL),
     + VDSAND(NN),CPSTON,ISTSAN
C     
C     If ICALL=1, alonshore uniformity is assumed for profile change
C     Compute the first-order rate of the bottom elevation change where
C     sediment transport rate Q(J) computed in Subr.11 SEDTRA.
C     The seaward boundary location at node 1 is chosen such that
C     bottom change is negligible seaward of node 1. Also at node JMAX
      IF(ICALL.EQ.1) THEN
      JMAXM1 = JMAX(L) - 1
      DZBDT(1) = 0.D0
      DZBDT(JMAX(L)) = 0.D0
      DO 100 J = 2,JMAXM1
        DZBDT(J) = (Q(J-1)-Q(J+1))/DX2
 100  CONTINUE
      IF(IVWALL(L).EQ.2) DZBDT(JMAXM1) = (Q(JMAX(L)-2) - Q(JMAXM1))/DX
C     where backward finite difference is used at the node next to the wall
C     if the vertical wall is exposed to wave action.
C     
C     Find the time step DELT using the numerical stability criterion
C     but the value of DELT is limited by the end time TIMEBC(ITIME+1) 
C     for given TIME 
C     Compute CB(J)=bottom profile phase velocity
C     CBMAX = 0.001D0
C     Increase of CBMAX tends to reduce DELT and improve numerical stability
C     CBMAX=0.05D0
      CBMAX=0.004D0
      DZBMAX=0.1D0*DX
      DO 115 J=1,JMAX(L)
        IF(J.EQ.1) THEN
          DELQ = Q(2) - Q(1)
          DELZB1 = ZB(2,L) - ZB(1,L)
        ELSEIF(J.EQ.JMAX(L)) THEN
          J1 = JMAXM1
          DELQ = Q(J) - Q(J1)
          DELZB1 = ZB(J,L) - ZB(J1,L)
        ELSE
          JP1 = J+1
          JM1 = J-1 
          DELQ = Q(JP1) - Q(JM1)
          DELZB1 = ZB(JP1,L) - ZB(JM1,L)
        ENDIF
        IF(DABS(DELZB1).GT.DZBMAX) THEN
          CB(J) = DELQ/DELZB1
        ELSE
          CB(J) =0.D0
        ENDIF
        DUMC = DABS(CB(J))
        IF(DUMC.GT.CBMAX) CBMAX=DUMC
 115  CONTINUE
      DELT = DX/CBMAX
      IDUM = ITIME + 1
      DUM = (TIMEBC(IDUM) - TIMEBC(ITIME))/2.D0
      IF(DELT.GT.DUM) DELT=DUM
C     
      DUM = TIME+DELT
      IF(DUM.GE.TIMEBC(IDUM)) THEN
        DELT = TIMEBC(IDUM) - TIME
        IEND = 1
      ENDIF
C     
C     Compute DELZBRW(J)=first-order bottom elevation change
C     before smoothing
      DO 120 J = 1,JMAX(L)
        DELZBRW(J) = DELT*DZBDT(J)
 120  CONTINUE
C     
C     Add second-order correction to DELZBRW(J)
      DTDT = DELT*DELT
      DO 121 J=1,JMAX(L)
        R(J) = DTDT*CB(J)*CB(J)/DXDX
 121  CONTINUE
      DO 122 J=2,JMAXM1
        JP1 = J+1
        JM1 = J-1
        DUM = ZB(JP1,L)*(R(JP1)+R(J))/4.D0-ZB(J,L)*(R(J)/2.D0+
     +     (R(JP1)+R(JM1))/4.D0)+ZB(JM1,L)*(R(J)+R(JM1))/4.D0
        DELZBRW(J) = DELZBRW(J)+ DUM
 122  CONTINUE
C
C     If ISTSAN=1, erosion is limited by available sand in stone structure
      IF(ISTSAN.EQ.1) THEN
          DO 123 J=1,JMAX(L)
              IF(HPSTON(J,L).GT.0.D0) THEN
                  DUM=ZP(J,L)-ZBSTON(J,L)
                  IF(DUM.GT.0.D0) THEN
                      VDSAND(J)=DUM+SNP*HPSTON(J,L)
                  ELSE
                      VDSAND(J)=SNP*(ZP(J,L)-ZPSTON(J,L))
                  ENDIF
                  DUM=VDSAND(J)+DELZBRW(J)
                  IF(DUM.LT.0.D0) DELZBRW(J)=-VDSAND(J)
               ELSE
                  VDSAND(J)=0.D0
               ENDIF
 123  CONTINUE
      ENDIF
C           
C     Smoothing DELZBRW using Subr.15 SMOOTH
      JMAXL=JMAX(L)
      CALL SMOOTH(JMAXL,DELZBRW,DELZBJ)
      IF(ISEDAV.EQ.2) THEN
          DO 125 J=1, JMAXL
              DUM=DELZBJ(J)+ZB(J,L)
              IF(DUM.LT.ZMESH(J,L)) DELZBJ(J)=ZMESH(J,L)-ZB(J,L)
 125      CONTINUE
      ENDIF
C     
C     Adjust smoothed bottom elevation change DELZB
C     to satisfy the volume conservation between J=1 to JMAX
      DUM = DELT*(Q(1)-Q(JMAX(L)))
      CALL INTGRL(JMAXL,DX,DELZBJ,AREA)
      ADJUST = (DUM-AREA)/(XB(JMAX(L))-XB(1))
      DO 130 J=1,JMAX(L)
        DELZB(J,L) = ADJUST+DELZBJ(J)
 130  CONTINUE
C     If ISTSAN=1, ZB(J,L) and ZP(J,L) at next time level are computed as follows
      IF(ISTSAN.EQ.1) THEN
      DO 140 J=1,JMAX(L)
              IF(HPSTON(J,L).LE.0.D0) THEN
                  ZB(J,L)=ZB(J,L)+DELZB(J,L)
                  ZP(J,L)=ZB(J,L)
              ELSE
                  IF(DELZB(J,L).GE.0.D0) THEN
                      DUM=SNP*HP(J,L)
                      IF(DELZB(J,L).GE.DUM) THEN
                          ZB(J,L)=ZB(J,L)+DELZB(J,L)-DUM
                          ZP(J,L)=ZB(J,L)
                      ELSE
                          ZB(J,L)=ZBSTON(J,L)
                          ZP(J,L)=ZP(J,L)+DELZB(J,L)/SNP
                      ENDIF
                  ELSE
                      IF(HP(J,L).GT.0.D0) THEN
                          ZB(J,L)=ZBSTON(J,L)
                          ZP(J,L)=ZP(J,L)+DELZB(J,L)/SNP
                      ELSE
                          DUM=ZB(J,L)-ZBSTON(J,L)+DELZB(J,L)
                          IF(DUM.GE.0.D0) THEN
                              ZB(J,L)=ZB(J,L)+DELZB(J,L)
                              ZP(J,L)=ZB(J,L)
                          ELSE
                              ZB(J,L)=ZBSTON(J,L)
                              ZP(J,L)=ZBSTON(J,L)+DUM/SNP
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
 140  CONTINUE
      ENDIF
C
      ENDIF
C..........End of ICALL=1..................................................
C     
C     If ICALL=2 from Main Program, the profile change due to alongshore
C     gradient of longshore sediment transport is included if IQYDY=1 and
C     computed when IEND=1 and L=ILINE
      IF(ICALL.EQ.2) THEN
      DO 200 LL=1, ILINE
        JMAXL=JMAX(LL)
      DO 210 J=1,JMAXL
        R(J)=VY(J,LL)
        DELZBRW(J)=DABS(ZB(J,LL)-DZX(J,LL))
 210  CONTINUE
      CALL SMOOTH(JMAXL,R,CB)
      CALL SMOOTH(JMAXL,DELZBRW,DELZBJ)
      CALL INTGRL(JMAXL,DX,CB,AREA)
      AVY(LL)=AREA
      CALL INTGRL(JMAXL,DX,DELZBJ,AREA)
      ADZX(LL)=AREA
      DO 211 J=1,JMAXL
        VY(J,LL)=CB(J)
        DZX(J,LL)=DELZBJ(J)
 211  CONTINUE
 200  CONTINUE
      CALL SMOOTH(ILINE,AVY,V)
      ILINE1=ILINE-1
      DO 220 LL=1, ILINE1
        VDY(LL)=(V(LL+1)-V(LL))/DYLINE(LL)
 220  CONTINUE
      AVY(1)=VDY(1)
      AVY(ILINE)=VDY(ILINE1)
      DO 230 LL=2, ILINE1
C     Use upstream finite difference method
      DUM=V(LL)*DYLINE(LL)
      IF(DUM.GE.0.D0) THEN
        AVY(LL)=VDY(LL-1)
      ELSE
        AVY(LL)=VDY(LL)
      ENDIF
 230  CONTINUE
      DO 240 LL=1, ILINE
        IF(ADZX(LL).LT.1.D-6) THEN
          AVY(LL)=AVY(LL)/1.D-6/SPORO1
        ELSE
          AVY(LL)=AVY(LL)/ADZX(LL)/SPORO1
        ENDIF
 240  CONTINUE
      DO 250 LL=1, ILINE
        JMAXL=JMAX(LL)
C       DUM=XB(JMAXL)*SPORO1
        DO 260 J=1, JMAXL
          DELZBRW(J)=-DZX(J,LL)*AVY(LL)
C         DELZBRW(J)=-AVY(LL)/DUM
 260  CONTINUE
      CALL SMOOTH(JMAXL,DELZBRW,DELZBJ)
      DO 270 J=1, JMAXL
        ZB(J,LL)=DELZBJ(J)+ZB(J,LL)
        IF(IPERM.EQ.1.OR.ISEDAV.GE.1) THEN
          HP(J,LL)=ZB(J,LL)-ZP(J,LL)
          IF(HP(J,LL).LT.0.D0) THEN
            HP(J,LL)=0.D0
            ZB(J,LL)=ZP(J,LL)
          ENDIF
          IF(ISEDAV.EQ.2) THEN
            IF(ZB(J,LL).LT.ZMESH(J,LL)) ZB(J,LL)=ZMESH(J,LL)
          ENDIF
        ENDIF
 270  CONTINUE
 250  CONTINUE
      ENDIF
C     
C..........End of ICALL=2..................................................
C     
      RETURN
      END
C     
C     -12-----------------  END OF SUBROUTINE CHANGE  ---------------------
C     #13##################### SUBROUTINE INTGRL ##########################
      SUBROUTINE INTGRL(NUM,DEL,F,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000)
      DIMENSION F(NN)
C     
C     NUM can be even or odd integer
      IMAX = (NUM-1)/2
      DUM = DBLE(NUM-1)/2.D0
      IF(DBLE(IMAX).LT.DUM) THEN
        NEND = NUM - 1
        NEVEN = 1
      ELSE
        NEND = NUM
      ENDIF
      SE = F(2)
      SO = 0.D0
      DO 500 I=2,IMAX
        SE = SE + F(I*2)
        SO = SO + F(I*2-1)
 500  CONTINUE
      G = DEL/3.D0*(F(1) + 4.D0*SE + 2.D0*SO + F(NEND))
      IF(NEVEN.EQ.1) G=G+(F(NEND)+F(NUM))*DEL/2.D0
      RETURN
      END
C     -13-----------------  END OF SUBROUTINE INTGRL  ---------------------
C     #14##################### SUBROUTINE SMOOTH ##########################
C     Smooth the vector RAW using NPT computed in Subr.3 BOTTOM
C     where NPFS = (2NPT+1) = number of points for smoothing
      SUBROUTINE SMOOTH(NUM,RAW,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000)
      DIMENSION RAW(NN),F(NN)
      COMMON /COMPAR/ HWDMIN,NPT,NPE
C     
      DO 201 J = 1,NUM
        JSTA = J-NPT
        JEND = J+NPT
        IF(JSTA.LT.1) THEN
          JSTA = 1
          JEND = 2*J-1
        ENDIF
        IF(JEND.GT.NUM) THEN
          JSTA = 2*J - NUM
          JEND = NUM
        ENDIF
        NPFS = JEND-JSTA+1
        TOTJ = DBLE(NPFS)
        SUM1 = 0.D0
        DO 202 JJ = JSTA, JEND
          SUM1 = SUM1+ RAW(JJ)
 202    CONTINUE
        F(J) = SUM1/TOTJ
 201  CONTINUE
C     
      RETURN
      END
C     -14-----------------  END OF SUBROUTINE SMOOTH  ---------------------
C     #15####################### SUBROUTINE EXTRAPO #######################
C     Extrapolate vector F from node J1 to node J2
C     where values of F at nodes J=1 to (J1-1) are computed 
      SUBROUTINE EXTRAPO(J1,J2,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000)
      DIMENSION F(NN)
      COMMON /COMPAR/ HWDMIN,NPT,NPE
C     
C     NPE = number of points (computed in Subr.3 BOTTOM) extrapolated
C     to avoid a sudden jump from computed F(J1-1) to zero
      JJ=J1+NPE
      IF(JJ.LE.J2) THEN
      JM = J1-1
      Y= F(JM)
      DELY = Y/DBLE(NPE+1)
      DO 100 J=1,NPE
        F(JM+J) = Y-DELY*DBLE(J)
 100  CONTINUE
      F(JJ:J2) = 0.D0
      ELSE
        IF(J1.LE.J2) F(J1:J2)=0.D0
      ENDIF
C     
      RETURN
      END
C     -15---------------- END OF SUBROUTINE EXTRAPO ------------------------
C     #16####################### SUBROUTINE WETDRY #########################
C     Compute swash hydrodynamics in wet/dry zone (possibly multiple bottom
C     peaks) and combined wave overtopping and overflow rate QOTF.
C     
      SUBROUTINE WETDRY(ITIME,L,ITEQO)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NB=30000,NL=100)
      DIMENSION G(NN), DG(NN), ETA(NN),ETAP(NN)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     +   SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX,FSY,FE,QWX,QWY
      COMMON /WBREAK/ GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /POROUS/ XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN), 
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON /WATRAN/ SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      COMMON /COMPAR/ HWDMIN,NPT,NPE
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
      COMMON /VEGETA/ VEGCD,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL),
     + VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL),VEGRH(NN,NL),
     + VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL)
C     
C     Compute swash variables for node J=JWD to JDRY
C     JWD= wet and dry transition node
C     JR= landward limit of wet computation in Main Program
C     JDRY= landward limit of wet and dry zone
C     PWET(J)= wet probability at node J where PWET(JWD)=1.0
C     HWD(J)= mean water depth where HWD(JWD)=H1
C     USWD(J)= steady swash velocity component
C     UMEAWD(J)= mean velocity in wet and dry zone
C     QP(J)= water flux in permeable layer
C     If INFILT=1, QP(J)=infiltration rate between dune crest and
C     landward node J where QP(J)=0.0 assumed seaward of dune crest
C     UPMWD(J)= mean discharge velocity in permeable layer
      IF(ITEQO.LE.2) THEN
        JWD = JSWL(L)
        IF(IWTRAN.EQ.1.AND.IOFLOW.EQ.1) JWD=JCREST(L)
C     If IPROFL=1 and JWD=JR, no wet and dry zone above SWL 
C     shoreline
C     To avoid possible overwash under small waves, no marching
        IF(IPROFL.EQ.1.AND.JWD.GT.JR) THEN
          IF(IPOND.EQ.0) THEN
            JWD=JR
            JDRY=JR
            H1=H(JR)
            GOTO 110
          ENDIF
        ENDIF
        IF(JWD.GT.JR) JWD=JR
        H1 = H(JWD)
      ENDIF
      HWD(JWD) = H1
      BGH3=BWD*GRAV*H1*H1*H1
C BDJ 2011->2014 on 2014-10-02 
      SSP_50 = 1.D0 ! The value of SSP such that the CORRECT
!     is .5*(1+GAMMA/SQR8)
      A = 1.D0      ! Dictates the steepness of blending curve near SSP_50
      CORRECT = GAMMA/SQR8
      CORRECT = 0.5D0*(1.D0+CORRECT)+ 
     +     0.5D0*(1.D0-CORRECT)*tanh(a*(SSP-SSP_50));
!      CORRECT = 1.D0  ! comment this line out to use correction
      SIGWD(JWD) = CORRECT*H1
!      SIGWD(JWD) = H1
C end C BDJ 2011->2014 on 2014-10-02 
      PMG1=AWD/DSQRT(GRAV)
      PMGH1=PMG1/H1

      PWET(JWD) = 1.D0
      IF(IPERM.EQ.0) THEN
        QWX=QO(L)
        IF(INFILT.EQ.1) QP(JWD)=0.D0
      ELSE
        QWX=QO(L)-QP(JWD)
        UPMWD(JWD)=UPMEAN(JWD)
        ETA(JWD)=H1+ZB(JWD,L)
        ETAP(JWD)=ZB(JWD,L)
      ENDIF
      QS = QWX-AQWD*H1*DSQRT(GRAV*H1)
      IF(QS.GT.0.D0) QS=0.D0
      USWD(JWD) = QS/H1
      UMEAWD(JWD) = AUWD*DSQRT(GRAV*H1) + USWD(JWD)
      DUM=AGWD*GRAV*H1 - (UMEAWD(JWD)-USWD(JWD))**2.D0
      IF(DUM.LT.0.D0) THEN
        JDRY=JWD
        GOTO 110
      ENDIF
      USTDWD(JWD)=DSQRT(DUM)
      A = QWX*QWX/BGH3
      A1=A
C     
C     Empirical formula for wet probability parameter n=WDN
      WDN=1.01D0+0.98D0*DTANH(QO(L)*QO(L)/BGH3)**0.3D0
      W1=WDN-1.D0
      BNWD=BWD*(1.D0+A1)*(2.D0-WDN)/(WDN-1.D0)
C     
C----------------LANDWARD MARCHING COMPUTATION ---------------
C     If IWTRAN=1, the landward wet zone starts from node J=JSL
      JEND=JMAX(L)-1
      IF(IWTRAN.EQ.1.AND.JSL.LT.JMAX(L)) THEN
        JEND=JSL1-1
      ENDIF
C     LSTART=1 indicates beginning of upslope computation
      LSTART=1
      IF(JWD.GT.JEND) THEN
        JDRY=JR
        GOTO 110
      ENDIF
      DO 100 J=JWD,JEND
        JP1 = J+1
C     
C------------------ BOTTOM ELEVATION INCREASING LANDWARD -------------------
C     On the seaward upslope and crest(J<JCREST but J<JMAX if IPOND=1)
C     use an empirical formula for wet probability PWET and
C     compute mean depth and return flow velocity USWD for IUPSLP=1
        IUPSLP=0
        JDUM=JCREST(L)
        IF(IPOND.EQ.1) JDUM=JMAX(L)
        IF(J.LT.JDUM.AND.ZB(JP1,L).GE.ZB(J,L)) IUPSLP=1
C     For J=JWD, IUPSLP=1 for any slope
        IF(J.EQ.JWD) IUPSLP=1
C     If IPOND=1 and NOPOND=0, ponded water zone is treated like 
C     downslope zone with IUPSLP=0
        IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
          IF(JP1.GE.JXW.AND.JP1.LE.JX2) IUPSLP=0
        ENDIF
        IF(IUPSLP.EQ.1) THEN
          IF(LSTART.EQ.1) THEN
            H2=HWD(J)
            BN12=BNWD*(H1/H2)**W1
            D = BN12 - ZB(J,L)/H1
            AH = AGWD/H1
            G(J) = 0.D0
            DUM=QWX-QS
            IF(DABS(DUM).LT.1.D-3) THEN
              R=0.D0
            ELSE
              R=CWD*QS/DUM
            ENDIF
            IF(IVEG.EQ.0)THEN
              DG(J) = AH*DX*FB2(J,L)*GBWD(R)
C BDJ 2011->2014 on 2014-10-02 
C  added to kill friction in wetdry
            DG(J) = 0.D0
C end BDJ 2011->2014 on 2014-10-02 

            ELSE
              DUM=HWD(J)/PWET(J)
              X=VEGH(J,L)/DUM
              DG(J)=AH*DX*(FB2(J,L)*GBWD(R)+VEGINP(J,L)*DUM*GDWD(R,X))
            ENDIF
C     Functions GBWD(R) and GDWD(R,X) are specified below this subroutine
            LSTART=0
          ENDIF
          CX = D + ZB(JP1,L)/H1
          IF(IPERM.EQ.0) THEN
            WPGH=0.D0
          ELSE
            IF(HP(JP1,L).LE.0.D0) THEN
              WPGH=0.D0
            ELSE
              DUM=0.5D0*(HP(J,L)+HP(JP1,L))/SDP
C             IF(DUM.GT.10.D0) DUM=10.D0
              PMGH=PMGH1*DUM**0.3D0
              WPGH=PMGH*WPM*PWET(J)*DX/DSQRT(HWD(J))
            ENDIF
          ENDIF
          DGJP1 = DG(J)+WPGH
          DO 200 ITEH=1,20
            G(JP1) = G(J) + DGJP1
            C = CX + G(JP1)
            IF(C.LE.0.D0) THEN
              JDRY = J
              GOTO 110
            ELSE
              Y = (C/BN12)**(1.D0/W1)
            ENDIF
            HWD(JP1)=H2/Y
            Y=H1/HWD(JP1)
            DUM = (1.D0 + A1)*Y**WDN - A*Y*Y*Y
C           IF(DUM.LE.0.D0) THEN
C             JDRY=J
C             GOTO 110
C           ENDIF
            IF(DUM.LT.1.D0) THEN
              PWET(JP1) = PWET(J)
            ELSE
              PWET(JP1) = 1.D0/DUM
              IF(PWET(JP1).GT.PWET(J)) PWET(JP1)=PWET(J)
            ENDIF
            QWAVE=AQWD*HWD(JP1)*DSQRT(GRAV*HWD(JP1)/PWET(JP1))
C     
C     Compute QP and UPMWD in permeable layer if IPERM=1
C     ETAP(JP1)=mean water level above datum inside permeable layer
C     where ETA(JP1) and ETAP(JP1) are mean water levels above datum
            IF(IPERM.EQ.1) THEN
              IF(HP(JP1,L).LE.0.D0) THEN
                UPMWD(JP1)=0.D0
                QP(JP1)=0.D0
                WPGH=0.D0
              ELSE
                ETA(JP1)=HWD(JP1)+ZB(JP1,L)
                DUM=ZP(JP1,L)
                IF(DUM.LT.SWLBC(ITIME).AND.ZP(JP1,L).GE.ZP(J,L)) 
     +             DUM=SWLBC(ITIME)
                ETAP(JP1)=ZB(JP1,L)*PWET(JP1)+DUM*(1.D0-PWET(JP1))
                IF(ETAP(JP1).LT.ZP(JP1,L)) ETAP(JP1)=ZP(JP1,L)
                C=(ETA(JP1)-ETA(J))/DX
                DUM=DSQRT(ALSTA2+BE4*DABS(C))
                UPMWD(JP1)=(DUM-ALSTA)/BE2
                IF(C.GT.0.D0) UPMWD(JP1)=-UPMWD(JP1)
                QP(JP1)=UPMWD(JP1)*(ETAP(JP1)-ZP(JP1,L))*PWET(JP1)
                DUM=QO(L)-QWAVE
                IF(QP(JP1).LT.DUM.AND.DABS(QP(JP1)).GT.1.D-5) THEN
                  UPMWD(JP1)=UPMWD(JP1)*DUM/QP(JP1)
                  QP(JP1)=DUM
                ENDIF
                DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                WPGH=PMGH*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
              ENDIF
              QWX=QO(L)-QP(JP1)
              A=QWX*QWX/BGH3
            ENDIF
            IF(INFILT.EQ.1) QP(JP1)=0.D0
C     
            QS = QWX-QWAVE
C     QS=return flux must be zero or negative for J<JCREST
            IF(QS.GT.0.D0) QS=0.D0
            USWD(JP1) = QS/HWD(JP1)
            DUM=QWX-QS
            IF(DABS(DUM).LT.1.D-3) THEN
              R=0.D0
            ELSE
              R = CWD*QS/DUM
            ENDIF
            IF(IVEG.EQ.0)THEN
              DG(JP1) = AH*DX*FB2(JP1,L)*GBWD(R)
C BDJ 2011->2014 on 2014-10-02 
C  added to kill friction in wetdry
              DG(JP1) = 0.D0
C end BDJ 2011->2014 on 2014-10-02 
            ELSE
              DUM=HWD(JP1)/PWET(JP1)
              X=VEGH(JP1,L)/DUM
              DG(JP1)=AH*DX*(FB2(JP1,L)*GBWD(R)+VEGINP(JP1,L)*DUM*GDWD(R
     +        ,X))
            ENDIF
            DUM=0.5D0*(DG(J)+DG(JP1))+WPGH
            IF(DABS(DUM-DGJP1).GT.1.D-5) THEN
              DGJP1 = DUM
              GOTO 200
            ELSE
              G(JP1)=G(J)+DUM
C     LSTART=2 indicates that bottom elevation is peaked at node JC
              IF(JP1.LT.JMAX(L)) THEN
                IF(ZB(J+2,L).LT.ZB(JP1,L)) THEN
                  JC=JP1
                  HC=HWD(JC)
                  PC = PWET(JC)
                  QWC=QWX
                  LSTART=2
                ELSE
                  IF(J.EQ.JWD) LSTART=1
                ENDIF
              ENDIF
              GOTO 220
            ENDIF
 200      CONTINUE
        ELSE
C     
C---------------------- BOTTOM ELEVATION DECREASING LANDWARD OR J>JCREST --------------
C     On the landward slope (J>JCREST) or downslope zone for J<JCREST or ponded water
C     zone, PWET=constant on impermeable bottom and compute HWD and USWD for IUPSLP=0
          IF(LSTART.EQ.2) THEN
            PCI=1.D0/PC
            QWC2=QWC*QWC
            BG=BWD*GRAV
            CPC = 0.5D0*PC/BWD/HC
            AB = 0.25D0*PC*QWC2/(BG*HC*HC*HC)
            G(J) = 0.D0
            QS=USWD(J)*HWD(J)
            DUM=QWC-QS
            IF(DABS(DUM).LT.1.D-3) THEN
              R=0.D0
            ELSE
              R=CWD*QS/DUM
            ENDIF
            IF(IVEG.EQ.0)THEN
              DG(J)=AGWD*DX*FB2(J,L)*GBWD(R)
            ELSE
              DUM=HWD(J)/PWET(J)
              X=VEGH(J,L)/DUM
              DG(J)=AGWD*DX*(FB2(J,L)*GBWD(R)+VEGINP(J,L)*DUM*GDWD(R,X))
            ENDIF
            LSTART=0
          ENDIF
          DZB=ZB(JC,L)-ZB(JP1,L)
          IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
            IF(JP1.GE.JXW.AND.JP1.LE.JX2) THEN
              IF(JX2.LT.JMAX(L)) DZB=ZB(JC,L)-ZW
              QWX=QO(L)-(QO(L)-QM)*(XB(JP1)-XB(JXW))/(XB(JX2)-XB(JXW))
              A=QWX*QWX/BGH3
            ENDIF
          ENDIF
          IF(IPERM.EQ.0) THEN
            WPGH=0.D0
            IF(INFILT.EQ.1) THEN
              WPGH=PMG1*WPM*PWET(J)*DX/DSQRT(HWD(J))
            ENDIF
          ELSE
            IF(HP(JP1,L).LE.0.D0) THEN
              WPGH=0.D0
            ELSE
              DUM=0.5D0*(HP(J,L)+HP(JP1,L))/SDP
C             IF(DUM.GT.10.D0) DUM=10.D0
              PMG=PMG1*DUM**0.3D0
              WPGH=PMG*WPM*PWET(J)*DX/DSQRT(HWD(J))
            ENDIF
          ENDIF
          DGJP1 = DG(J)+WPGH
          DO 210 ITEH=1,20
            G(JP1)= G(J)+DGJP1
            C=CPC*(DZB-G(JP1))
            IF(C.LT.0.D0) C=0.D0
            IF(HC.GT.1.D-6) THEN
              Y= HWD(J)/HC
            ELSE
              JDRY=J
              GOTO 110
            ENDIF
 205        DUM=1.D0/Y/Y
            F=Y-1.D0+AB*(DUM-1.D0)-C
            DF=1.D0-2.D0*AB*DUM/Y
            IF(DABS(DF).LT.1.D-6) THEN
              JDRY=J
              GOTO 110
            ENDIF
            YNEW=Y-F/DF
            IF(DABS(YNEW-Y).GT.1.D-6) THEN
              Y=YNEW
              GOTO 205
            ENDIF
            HWD(JP1)=YNEW*HC
            IF(HWD(JP1).LT.1.D-6) THEN
              JDRY=J
              GOTO 110
            ENDIF
            IF(HWD(JP1).GT.HWD(J)) HWD(JP1)=HWD(J)
            IF(IPERM.EQ.0.AND.INFILT.EQ.0) THEN
              PWET(JP1)=PC
            ELSE
              DUM=PCI+(QWC2-QWX*QWX)/BG/HWD(JP1)**3.D0
C             IF(DUM.LE.0.D0) THEN
C               JDRY=J
C               GOTO 110
C             ENDIF
              IF(DUM.LT.PCI) DUM=PCI
              PWET(JP1)=1.D0/DUM
              IF(PWET(JP1).GT.PWET(J)) PWET(JP1)=PWET(J)
            ENDIF
            QWAVE=AQWD*HWD(JP1)*DSQRT(GRAV*HWD(JP1)/PWET(JP1))
C     
C     Compute QP and UPMWD in permeable layer if IPERM=1 as above
            IF(IPERM.EQ.1) THEN
              IF(HP(JP1,L).LE.0.D0) THEN
                UPMWD(JP1)=0.D0
                QP(JP1)=0.D0
                WPGH=0.D0
              ELSE
                ETA(JP1)=HWD(JP1)+ZB(JP1,L)
                DUM=ZP(JP1,L)
                IF(DUM.LT.SWLBC(ITIME).AND.ZP(JP1,L).GE.ZP(J,L))
     +             DUM=SWLBC(ITIME)
                ETAP(JP1)=ZB(JP1,L)*PWET(JP1)+DUM*(1.D0-PWET(JP1))
                IF(ETAP(JP1).LT.ZP(JP1,L)) ETAP(JP1)=ZP(JP1,L)
                C=(ETA(JP1)-ETA(J))/DX
                DUM=DSQRT(ALSTA2+BE4*DABS(C))
                UPMWD(JP1)=(DUM-ALSTA)/BE2
                IF(C.GT.0.D0) UPMWD(JP1)=-UPMWD(JP1)
                QP(JP1)=UPMWD(JP1)*(ETAP(JP1)-ZP(JP1,L))*PWET(JP1)
C               DUM=QO(L)-QWAVE
C               IF(J.GE.JCREST(L).AND.QP(JP1).GT.DUM) THEN 
C                 UPMWD(JP1)=UPMWD(JP1)*DUM/QP(JP1)
C                 QP(JP1)=DUM
C               ENDIF
                DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                WPGH=PMG*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
              ENDIF
              QWX=QO(L)-QP(JP1)
              A=QWX*QWX/BGH3
            ENDIF
C     
C     Compute QP(J)=infiltration rate landward of dune crest if INFILT=1
            IF(INFILT.EQ.1) THEN
              IF(JP1.GT.JCREST(L)) THEN
                QP(JP1)=QP(J)+0.5D0*DX*WPM*(PWET(J)+PWET(JP1))
                DUM=WPM*DX*0.5D0*(PWET(JP1)+PWET(J))
                WPGH=PMG1*DUM/DSQRT(0.5D0*(HWD(JP1)+HWD(J)))
                QWX=QO(L)-QP(JP1)
                A=QWX*QWX/BGH3
              ELSE
                QP(JP1)=0.D0
                WPGH=0.D0
              ENDIF
            ENDIF
            QS = QWX - QWAVE
C           QS= steady flux on landward slope (J>JCREST) must be zero or positive
            IF(IPOND.EQ.0) THEN
              IF(J.GE.JCREST(L).AND.QS.LT.0.D0) QS=0.D0
            ENDIF
            IF(HWD(JP1).LT.1.D-3.AND.QS.GT.1.D-3) QS=1.D-3
            USWD(JP1) = QS/HWD(JP1)
            DUM=QWX-QS
            IF(DABS(DUM).LT.1.D-3) THEN
              R=0.D0
            ELSE
              R=CWD*QS/DUM
            ENDIF
            IF(IVEG.EQ.0)THEN
              DG(JP1)=AGWD*DX*FB2(JP1,L)*GBWD(R)
            ELSE
              DUM=HWD(JP1)/PWET(JP1)
              X=VEGH(JP1,L)/DUM
              DG(JP1)=AGWD*DX*(FB2(JP1,L)*GBWD(R)+VEGINP(JP1,L)*DUM*GDWD
     +        (R,X))
            ENDIF
            DUM=0.5D0*(DG(J)+DG(JP1))+WPGH
            IF(DABS(DUM-DGJP1).GT.1.D-5) THEN
              DGJP1 = DUM
              GOTO 210
            ELSE
              G(JP1)=G(J)+DUM
C     LSTART=1 indicates beginning of upslope computation
            IF(IPOND.EQ.0.OR.NOPOND.EQ.1) THEN
              IF(JP1.LT.JCREST(L)) THEN
                IF(ZB(J+2,L).GE.ZB(JP1,L)) LSTART=1
              ENDIF
              ELSE
                IF(JP1.EQ.JX2) THEN
                  LSTART=1
                  QWX=QM
                  A=QWX*QWX/BGH3
                ENDIF
              ENDIF
              GOTO 220
            ENDIF
 210      CONTINUE
        ENDIF
C----------------- END OF BOTTOM ELEVATION INCREASING OR DECREASING ----------------
C     
C     Compute mean velocity UMEAWD and standard deviations
C     SIGWD and USTDWD in wet and dry zone
 220    UMEAWD(JP1)=AUWD*DSQRT(GRAV*PWET(JP1)*HWD(JP1))
     +     +PWET(JP1)*USWD(JP1)
C BDJ 2011->2014 on 2014-10-02 
C        SIGWD(JP1)=HWD(JP1)*DSQRT(2.D0/PWET(JP1)-2.D0+PWET(JP1))
        SIGWD(JP1)=CORRECT*HWD(JP1)*DSQRT(2.D0/PWET(JP1)-2.D0+PWET(JP1))
C end BDJ 2011->2014 on 2014-10-02 
        DUM = UMEAWD(JP1) - USWD(JP1)
        DUM1 = PWET(JP1)*DUM**2.D0-2.D0*DUM*(UMEAWD(JP1)-
     +     PWET(JP1)*USWD(JP1))
        DUM = AGWD*GRAV*HWD(JP1)+DUM1
        IF(DUM.GT.0.D0) THEN
          USTDWD(JP1) = DSQRT(DUM)
        ELSE
          JDRY=J
          GOTO 110
        ENDIF
        IF(IANGLE.EQ.1) THEN
          STHETA(JP1)=STHETA(JWD)
          VMEAWD(JP1)=AUWD*DSQRT(GRAV*PWET(JP1)*HWD(JP1))*STHETA(JP1)
          DUM=1.D0-0.25D0*PI*PWET(JP1)*(2.D0-PWET(JP1))
          VSTDWD(JP1)=AWD*DSQRT(GRAV*HWD(JP1)*DUM)*DABS(STHETA(JP1))
        ENDIF
C     Mean water depth limited by HWDMIN specified in Subr.2 INPUT
C     Horizontal distance of wet and dry zone is limited because of 
C     assumed alongshore uniformity
C       DUM=(XB(JP1)-XB(JWD))/H1
C       IF(HWD(JP1).LT.HWDMIN.OR.DUM.GT.1000D0) THEN
        IF(HWD(JP1).LT.HWDMIN) THEN
          JDRY = JP1
C       IF(DUM.GT.1000D0.AND.JP1.GT.JCREST(L)) JMAX(L)=JP1
          GOTO 110
        ENDIF
C     
        IF(J.EQ.JEND) JDRY=JP1
 100  CONTINUE  
C-------------------END OF LANDWARD MARCHING --------------------------
C     
 110  CONTINUE
C     
C     QOTF=Combined overtopping and overflow rate QOTF
C     SPRATE=seepage rate through permeable layer predicted by modified
C     formula of Kobayashi and de los Santos(2007) for no overtopping
C     where USWD=0.0(unidirectional flow) at JCREST is assumed
      QOTF=0.D0
      SPRATE=0.D0
      JDAM=JSWL(L)
      IF(IPROFL.EQ.2) JDAM=JCREST(L)
      IF(JDRY.GE.JCREST(L).AND.JDAM.LT.JMAX(L)) THEN
        J=JCREST(L)
        IF(JWD.EQ.JMAX(L)) J=JMAX(L)
        QOTF = AQWD*HWD(J)*DSQRT(GRAV*HWD(J)/PWET(J))
        IF(IPERM.EQ.1) THEN
          IF(JDRY.EQ.JMAX(L).OR.IWTRAN.EQ.1) THEN
            SPRATE=QP(JCREST(L))
            IF(SPRATE.LT.0.D0) SPRATE=0.D0
C         ELSE
C           IF(IWTRAN.EQ.0) QOTF=0.D0
          ENDIF
        ENDIF
      ENDIF
      IF(IPOND.EQ.1.AND.NOPOND.EQ.0) THEN
        QD=QOTF
        IF(JDRY.EQ.JMAX(L)) THEN
          IF(ZW.LT.ZB(JMAX(L),L)) THEN
            QM=AQWD*HWD(JMAX(L))*DSQRT(GRAV*HWD(JMAX(L))/PWET(JMAX(L)))
            IF(QM.GT.QOTF) QM=QOTF
          ELSE
            QM=QOTF
          ENDIF
        ELSE
          QM=0.D0
        ENDIF
        IF(JCREST(L).EQ.JXW) QOTF=QM
      ENDIF
      IF(IPERM.EQ.1.AND.QOTF.EQ.0.D0) THEN
C     Find node JDUM for highest and most landward point of ZP(J,L)
        JSEP=JR
        IF(ZB(JSEP,L).LE.SWLBC(ITIME)) GOTO 301
        IF(IWTRAN.EQ.0) THEN
          JDUM=JSEP
          DUM=ZP(JSEP,L)
          DO 300 J=(JSEP+1),JMAX(L)
            IF(ZP(J,L).GE.DUM) THEN
              DUM=ZP(J,L)
              JDUM=J
            ENDIF
 300      CONTINUE
          DETA=ZB(JSEP,L)-ZP(JDUM,L)
        ELSE
          JDUM=JSL1
          DETA=ZB(JSEP,L)-ZB(JDUM,L)
        ENDIF
        IF(DETA.GT.0.D0) THEN
          DUM=XB(JDUM)-XB(JSEP)
          IF(DUM.LE.0.D0) GOTO 301
          SPRATE=0.2D0*DETA**1.5D0/DSQRT(BESTA1*DUM)
        ENDIF
 301    CONTINUE
      ENDIF
C     
      RETURN
      END
C     
C     -16---------------- END OF SUBROUTINE WETDRY -------------------------
C     **********************************************************************
C     This function related to bottom shear stress in wet and dry zone 
C     is called from Subr.16 WETDRY
      FUNCTION GBWD(R)
      DOUBLE PRECISION R,R2,GBWD,ERFCC
      IF(R.GE.0.D0) THEN
        GBWD = 1.D0 + 1.77245D0*R + R*R
      ELSE
        R2 = R*R
        GBWD = 2.D0*DEXP(-R2)-R2-1.D0+1.77245D0*R*(3.D0-2.D0*ERFCC(R))
      ENDIF
C     Complementary error function ERFCC below Subr.06 GBXAGF
      RETURN
      END
C     ***********************************************************************
C     This function related to vegetation drag force in wet and and dry zone
C     is called from Subr.16 WETDRY
      FUNCTION GDWD(R,X)
      DOUBLE PRECISION GDWD,R,X,EX,SX,FX,R2,ER2,FR,C,ERFCC
C     IF(X.LE.0.D0) THEN
C       GDWD = 0.D0
C     ELSE
C       GDWD=2.D0-(X+2.D0)*DEXP(-X)
        EX=DEXP(-X)
        SX=DSQRT(X)
        FX=1.D0-ERFCC(SX)
        R2 = R*R
        IF(R.GE.0.D0)THEN
          GDWD=2.D0-(X+2.D0)*EX+R*(1.77245D0*X-3.D0*SX*EX+1.77245D0*(1.5
     +    D0-X)*FX)+R2*(1.D0-EX)
        ELSE
          ER2=DEXP(-R2)
          FR=1.D0-ERFCC(R)
          C=(X+2.D0+R2+3.D0*R*SX)*EX
          IF(X.LE.R2)THEN
           GDWD=2.D0*X*ER2-2.D0-R2+1.77245D0*R*((X-1.5D0)*FX+2.D0*X*FR+X
     +     )+C
          ELSE
           GDWD=4.D0*ER2-2.D0-R2+5.317362D0*R*FR+1.77245D0*R*(X+(1.5D0-X
     +     )*FX)-C
          ENDIF
        ENDIF
C     ENDIF
      RETURN
      END
C     ***********************************************************************
C     #17####################### SUBROUTINE TRANWD #########################
C     Connect vector F1(J) with J=1 to JR with vector F2(J)
C     with J=JS to JE where JE is not less than JR
C     
      SUBROUTINE TRANWD(F1,JR,F2,JS,JE)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NN=5000)
      DIMENSION F1(NN),F2(NN)
C     
      IF(JR.GE.JS) THEN
        DO 100 J=JS,JR
          F1(J) = 0.5D0*(F1(J)+F2(J))
 100    CONTINUE
        DO 105 J=(JR+1),JE
          F1(J) = F2(J)
 105    CONTINUE
      ENDIF
      IF(JR.LT.JS) THEN
        DSR=DBLE(JS-JR)
        DO 200 J=(JR+1),JS
          DUM=DBLE(JS-J)/DSR
          F1(J)=F1(JR)*DUM+F2(JS)*(1.D0-DUM)
 200    CONTINUE
        DO 205 J=(JS+1),JE
          F1(J)=F2(J)
 205    CONTINUE
      ENDIF
C     
      RETURN
      END
C     -17---------------------- END OF SUBROUTINE TRANWD ----------------------
C     #18########################## SUBROUTINE PROBWD ##########################
C     Compute bedload probability PBWD(J) and suspended load
C     probability PSWD(J) in wet and dry zone in Subr.11 SEDTRA
C     
      SUBROUTINE PROBWD(PW,A,US,UC,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      IF(DABS(US).LE.UC) THEN
        P=PW*DEXP(-A*(UC-US)**2)
      ELSE
        IF(US.GT.UC) THEN
          P=PW
        ELSE
          P=PW*(1.D0-DEXP(-A*(UC+US)**2)+DEXP(-A*(UC-US)**2))
        ENDIF
      ENDIF
C     
      RETURN
      END
C     -18----------------------- END OF SUBROUTINE PROBWD ---------------------
C     #19########################### SUBROUTINE TSINTP ########################
C     
C     This subroutine interpolates time series W1(I) specified
C     at time levels T1(I) with I=1,2,...,(N1+1), 
C     obtain time series F(K) at time levels T2(K)
C     with K=1,2,...,(N2+1) and compute average value W2(K)
C     with K=1,2,...,N2 between time levels T2(K) and T2(K+1)
C     
      SUBROUTINE TSINTP(N1,T1,W1,N2,T2,W2)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NB=30000)
      DIMENSION T1(NB),W1(NB),T2(NB),W2(NB),F(NB)
C     
      F(1)=W1(1)
      F(N2+1)=W1(N1+1)
      IF(N2.GE.2) THEN
        DO 100 K=2,N2
          DO 200 I=1,N1
            I1=I+1
            IF(T1(I).LE.T2(K).AND.T2(K).LT.T1(I1)) THEN
              DUM=(T2(K)-T1(I))/(T1(I1)-T1(I))
              F(K)=(1.D0-DUM)*W1(I)+DUM*W1(I1)
              GOTO 100
            ENDIF
 200      CONTINUE
 100    CONTINUE
      ENDIF
C     
      DO 300 K=1,N2
        W2(K)=0.5D0*(F(K)+F(K+1))
 300  CONTINUE
C     
      RETURN
      END
C     
C     -19------------------ END OF SUBROUTINE TSINTP --------------------------
C     #20##################### SUBROUTINE PONDED ##############################
C     
C     This subroutine computes ponded water level and zone if IPOND=1
C     
      SUBROUTINE PONDED(L)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NL=100)
C     
      COMMON /OPTION/ TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /PROCOM/ DELT,DELZB(NN,NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /RRPOND/ ZW,QD,QM,JXW,JX2,NOPOND
C   
C     Compute the following quantities for known bottom 
C     profile ZB(J,L) above datum at XB(J) for node J for line L 
C     ZW = ponded water level at time=TIME
C     JCREST(L) = ridge crest node landward of SWL node JSWL
C     JXW = seaward end node of ponded water zone
C     JX2 = landward end node of ponded water zone
C     
C     For TIME=0.0, ZW=SWLBC(1) as specified in Subr.2 INPUT
C     For TIME>0, compute ZW at present time level
      IF(TIME.GT.0.D0) THEN
        IF(JX2.GT.JXW) THEN
          ZW=ZW+DELT*(QO(L)-QM)/(XB(JX2)-XB(JXW))
        ENDIF
        IF(ZW.GT.ZB(JMAX(L),L)) ZW=ZB(JMAX(L),L)
      ENDIF
C     
C     NOPOND=0 for ponded water in runnel
C     NOPOND=1 for submerged ridge and runnel seaward  
C     of node JSWL(L) or dry runnel with no ponded water
      JRUN=JMAX(L)
      JPEAK=JMAX(L)
      DO 100 J=(JSWL(L)+1),(JMAX(L)-1)
        IF(ZB(J-1,L).GE.ZB(J,L).AND.ZB(J,L).LT.ZB(J+1,L)) THEN
          IF(ZB(J,L).LT.ZB(JRUN,L).AND.ZW.GT.ZB(J,L)) JRUN=J
        ENDIF
        IF(ZB(J,L).GE.ZB(JPEAK,L)) JPEAK=J
 100  CONTINUE
      IF(JRUN.EQ.JMAX(L)) THEN
        NOPOND=1
        JCREST(L)=JPEAK
        RCREST(L)=ZB(JCREST(L),L)
        JXW=JSWL(L)
        JX2=JMAX(L)
        ZW=ZB(JSWL(L),L)
        GOTO 200
      ELSE
        NOPOND=0
      ENDIF
C     For NOPOND=1, node JCREST(J) is highest bottom elevation and
C     water level ZW is set to be still water level.
C     
C     JCREST(L) = node of ridge crest located between 
C     nodes JSWL(L) and JRUN if NOPOND=0
      JCREST(L)=JSWL(L)
      DO 110 J=(JSWL(L)+1),(JRUN-1)
        IF(ZB(J-1,L).LE.ZB(J,L).AND.ZB(J,L).GT.ZB(J+1,L)) THEN
          IF(ZB(J,L).GT.ZB(JCREST(L),L)) JCREST(L)=J
        ENDIF
 110  CONTINUE
C     
      IF(JCREST(L).EQ.JSWL(L)) THEN
        NOPOND=1
        JCREST(L)=JPEAK
        RCREST(L)=ZB(JCREST(L),L)
        JXW=JSWL(L)
        JX2=JMAX(L)
        ZW=ZB(JSWL(L),L)
        GOTO 200
      ENDIF
      RCREST(L)=ZB(JCREST(L),L)
C     If ponded water in runnel is full landward of ridge
C     crest, lower ZW to ridge crest elevation
      IF(ZW.GT.ZB(JCREST(L),L)) ZW=ZB(JCREST(L),L)
C     
C     Find nodes JXW and JX2 at water level ZW
      J=JCREST(L)
 120  IF(ZB(J,L).LE.ZW) THEN
        JXW=J
        GOTO 121
      ELSE
        J=J+1
        IF(J.EQ.JRUN) THEN
          JXW=JRUN-1
        GOTO 121
        ENDIF
        GOTO 120
      ENDIF
 121  J=JRUN
 125  IF(ZB(J,L).GT.ZW) THEN
        JX2=J-1
        GOTO 200
      ELSE
        J=J+1
        IF(J.EQ.JMAX(L)) THEN
          JX2=JMAX(L)
          GOTO 200
        ENDIF
        GOTO 125
      ENDIF
C
 200  CONTINUE
      RETURN
      END
C     
C     -20------------------------- END OF SUBROUTINE PONDED -----------------------
C     #21##################### SUBROUTINE WTRANS ##############################
C     
C     This subroutine computes transmitted waves (IWTRAN=1) landward of
C     structure along cross-shore line L if landward standing water exists
C     
      SUBROUTINE WTRANS(ITIME,L)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NB=30000,NL=100)
C     
      COMMON/OPTION/TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON/PREDIC/HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON/BPROFL/DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON/CONSTA/GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON/LINEAR/WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX,FSY,FE,QWX,QWY
      COMMON/VELOCY/UMEAN(NN),USTD(NN),USTA(NN),VMEAN(NN),VSTD(NN),
     +   VSTA(NN)
      COMMON/POROUS/XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON/OVERTF/RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON/SWASHP/AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
      COMMON/SWASHY/PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON/WATRAN/SWLAND(NB),ISWLSL,JSL,JSL1,IOFLOW
      DATA CBREAK/10.0D0/
C     
C     JDRY=landward end node of wet and dry zone computation
C     JSL=most seaward wet node in landward wet zone
C     JSL1=(JSL-1)=most landward node in wet and dry zone
C     If JDRY is less than JSL1, no wave overtopping occurs but 
C     wave transmission through a porous structure is included
C     
C     In landward wet zone, volume flux is constant where computed
C     wave overtopping rate QO(L) along cross-shore line L includes
C     water flux QP in porous layer whose vertical thickness is HP
      IF(IPERM.EQ.1) QPHP=QP(JSL1)/HP(JSL1,L)
C     which is assumed to be constant in landward wet zone
C     
      ICHECK=0
      DUM=CBREAK*H(JSL1)
      IF(DUM.LT.0.01D0) DUM=0.01D0
      DO 100 J=JSL, JMAX(L)
        PWET(J)=1.D0
        IF(ICHECK.EQ.0) THEN
          IF(SWLDEP(J,L).GT.DUM) THEN
            JLIN=J
            ICHECK=1
          ENDIF
        ENDIF
        IF(IPERM.EQ.1) THEN
          UPMEAN(J)=QPHP
          IF(HP(J,L).LE.0.D0) UPMEAN(J)=0.D0
          QP(J)=UPMEAN(J)*HP(J,L)
        ENDIF
        IF(IANGLE.EQ.1) THEN
          STHETA(J)=0.D0
          VMEAN(J)=0.D0
          VSTD(J)=0.D0
        ENDIF
  100   CONTINUE
        IF(ICHECK.EQ.0)JLIN=JMAX(L)
C     
        WSETUP(JSL1)=H(JSL1)+ZB(JSL1,L)-SWLAND(ITIME)
        CSIGMA=SIGMA(JSL1)
        JDUM=JSWL(L)
        IF(CSIGMA.GT.SIGMA(JDUM)) CSIGMA=SIGMA(JDUM)
        IF(CSIGMA.GT.SIGMA(1)) CSIGMA=SIGMA(1)
        DO 110 J=JLIN,JMAX(L)
        H(J)=SWLDEP(J,L)+WSETUP(JSL1)
C       H(J)=SWLDEP(J,L)
        WSETUP(J)=WSETUP(JSL1)
C       WSETUP(J)=0.D0
        SIGMA(J)=CSIGMA
        IF(H(J).LT.1.D-3) H(J)=1.D-3
        CP(J)=DSQRT(GRAV*H(J))
        QWX=QO(L)
        IF(IPERM.EQ.1)QWX=QWX-QP(J)
        IF(SIGMA(J).LE.0.D0) THEN
          IF(QWX.LE.0.D0) THEN
            SIGMA(J)=0.D0
          ELSE
            SIGMA(J)=DSQRT(QWX*H(J)/CP(J))
            IF(SIGMA(J).GT.SIGMA(JSWL(L))) SIGMA(J)=SIGMA(JSWL(L))
            IF(SIGMA(J).GT.SIGMA(1)) SIGMA(J)=SIGMA(1)
          ENDIF
        ENDIF
        SIGSTA(J)=SIGMA(J)/H(J)
        UMEAN(J)=QWX/H(J)-CP(J)*SIGSTA(J)*SIGSTA(J)
C       UMEAN(J)=0.D0
        USTD(J)=CP(J)*SIGSTA(J)
  110   CONTINUE
C     
C     Linear interpolation for transition zone from node JSL1
C     where wet probability is less than unity to wet node JLIN
C     Assume WSETUP(J) above SWL and SIGMA(J) remain constant
C     for nodes J=JSL1,...,JMAX(J) for wave transmission
      DUM=DBLE(JLIN-JSL1)
      IF(DUM.LE.1.D0) GOTO 999
      DO 120 J=JSL,JLIN-1
        DJ=DBLE(J-JSL1)/DUM
        DJ1=1.D0-DJ
        WSETUP(J)=DJ1*WSETUP(JSL1)+DJ*WSETUP(JLIN)
        H(J)=WSETUP(J)+SWLDEP(J,L)
        SIGMA(J)=DJ1*SIGMA(JSL1)+DJ*SIGMA(JLIN)
        UMEAN(J)=DJ1*UMEAN(JSL1)+DJ*UMEAN(JLIN)
        USTD(J)=DJ1*USTD(JSL1)+DJ*USTD(JLIN)
  120 CONTINUE
C     
  999 CONTINUE
      RETURN
      END
C     
C     -21------------------------- END OF SUBROUTINE WTRANS -----------------------
C     #22############################ SUBROUTINE EROSON ###########################
C     
C     This subroutine computes erosion of grassed dike at time level 
C     ITIME and along cross-shore line L if IPROFL=2 and ICLAY=0
C     For ICLAY=1 and IPROFL=1, exposed clay erosion is computed
C     
      SUBROUTINE EROSON(ITIME,L,IEND)
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=5000,NB=30000,NL=100)
      DIMENSION GRS1(NN,NL),GRS2(NN,NL),GRS3(NN,NL),GRS4(NN,NL),
     +   GRS5(NN,NL),FBA3(NN,NL),DFSWD(NN),BSF(NN),DUMVEC(NN)
C     
      COMMON/OPTION/TIME,IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT,
     +   ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY,
     +   IVEG,ICLAY
      COMMON /SEAWAV/ TIMEBC(NB),TPBC(NB),HRMSBC(NB),WSETBC(NB),
     +   SWLBC(NB),WANGBC(NB),NWAVE,NSURG,NWIND,NTIME
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     +   SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /WBREAK/ GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN)
      COMMON /ENERGY/ EFSTA(NN),DFSTA(NN)
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /SEDINP/WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2,
     +   EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA
      COMMON /PROCOM/ DELT,DELZB(NN,NL)
      COMMON /ROLLER/ RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN)
      COMMON /POROUS/XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL),
     +   WNU,SNP,SDP,ALPHA,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN),
     +   UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN),NPINP(NL)
      COMMON /OVERTF/ RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT,JCREST(NL)
      COMMON /SWASHP/ AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,
     +   BE2,BE4
      COMMON /SWASHY/ PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),
     +   USTDWD(NN),VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),
     +   H1,JWD,JDRY
      COMMON /DIKERO/EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN),
     +   GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL),
     +   GRSR(NN,NL),GRSRD(NN,NL),DEEB,DEEF
      COMMON /SOCLAY/EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL),
     + FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL)
C     
C     Dike Erosion efficiencies eB and ef associated with wave breaking 
C     and bottom friction are specified in Subr.2 INPUT
      DATA DELEM,SCP/0.05D0,1.2D0/
C     where DELEM=maximum allowable dike erosion increment (m) and SCP=
C     maximum slope of eroded dike clay soil
C     
C     The following parameters are invariant with time
      IF(TIME.EQ.0.D0)THEN
        DO 100 J=1,JMAX(L)
          FBA3(J,L)=FB2(J,L)*DSQRT(GRAV)*AWD**3.D0
          IF(ICLAY.EQ.1) THEN
              EPCLAY(J,L)=0.D0
              BSF(J)=1.0D0
          ELSE
              GRS3(J,L)=GRAV/GRSRD(J,L)
              GRS4(J,L)=0.5D0*GRSD(J,L)*(GRSR(J,L)-GRSRD(J,L))/
     +         GRSRD(J,L)
              GRS5(J,L)=0.5D0*GRSD(J,L)*(GRSR(J,L)+GRSRD(J,L))/GRAV
              IF(GRSD(J,L).LE.0.D0) THEN
                  GRS1(J,L)=0.D0
                  GRS2(J,L)=0.D0
              ELSE
                  DUM=GRSR(J,L)-GRSRD(J,L)
                  GRS1(J,L)=GRSD(J,L)*GRSR(J,L)/DUM
                  GRS2(J,L)=2.D0*GRAV*DUM/GRSD(J,L)/GRSR(J,L)**2.D0
              ENDIF
              DSUM(J)=0.D0
          ENDIF
 100    CONTINUE
      ENDIF
C     
C     BSF(J)=bottom slope function for dike erosion
C     where computed values are smoothed to obtain BSF(J)
      IF(ICLAY.EQ.0) THEN
          DO 150 J=1,JDRY
              ASB=DABS(BSLOPE(J,L))
              DUM=ASB/SCP
              IF(DUM.GE.0.9D0)THEN
                  DUMVEC(J)=10.D0
              ELSE    
                  DUMVEC(J)=1.D0/(1.D0-DUM)
              ENDIF
 150  CONTINUE
      CALL SMOOTH(JDRY,DUMVEC,BSF)
      ENDIF
C      
C     DSTA(J)=dike erosion forcing at given TIME
C     DSUM(J)=value of DSTA(J) integrated from TIME=0.0
C     DSTA(J,L) is computed for wet nodes (J=1 to JR) and for wet and
C     dry nodes (J=JWD to JDRY) separately
      DO 200 J=1,JR
        IF(IROLL.EQ.0) THEN
          DSTA(J)=DEEB*DBSTA(J)+DEEF*DFSTA(J)
        ELSE
          DSTA(J)=DEEB*RBETA(J)*RQ(J)+DEEF*DFSTA(J)
        ENDIF
        DSTA(J)=DSTA(J)*BSF(J)
 200  CONTINUE
      ED=1.D0
      DO 210 J=JWD,JDRY
        DUM=AQWD*H(J)*DSQRT(GRAV*H(J)/PWET(J))
        IF(DUM.LT.1.D-6) THEN
          RS=0.D0
        ELSE
          RS=CWD*(QO(L)-DUM)/DUM
        ENDIF
        DFSWD(J)=ED*FBA3(J,L)*H(J)*DSQRT(H(J)/PWET(J))*GFDWD(RS)
        DFSWD(J)=DFSWD(J)*BSF(J)
C     Function GFDWD(R) is specified below this subroutine
        IF(J.EQ.JWD)THEN
          ED=DSTA(J)/DFSWD(J)
          DFSWD(J)=DSTA(J)
        ENDIF
 210  CONTINUE
C     
C     Connect DSTA(J) and DFSWD(J) and smooth connected DSTA(J)
C     using Subr.17 TRANWD, Subr.14 SMOOTH and Subr.15 EXTRAPO
      IF(JDRY.GT.JR) THEN
        CALL TRANWD(DSTA,JR,DFSWD,JWD,JDRY)
      ELSE
        JDRY=JR
      ENDIF
      DUMVEC=DSTA
      CALL SMOOTH(JDRY,DUMVEC,DSTA)
      IF(JDRY.LT.JMAX(L)) THEN
        JDRY1=JDRY+1
        CALL EXTRAPO(JDRY1,JMAX(L),DSTA)
      ENDIF
C     
C     Find time step size DELT based on DELEM in DATA for time
C     marching computation of numerical stability if ICLAY=0
      IF(ICLAY.EQ.0) THEN
      DMAX=DSTA(1)*GRS3(1,L)
      DO 300 J=2,JMAX(L)
        DUM=DSTA(J)*GRS3(J,L)
        IF(DUM.GT.DMAX) DMAX=DUM
 300  CONTINUE
      IF(DMAX.LT.1.D-6) DMAX=1.D-6
      DELT=DELEM/DMAX
      IDUM=ITIME+1
      DUM=(TIMEBC(IDUM)-TIMEBC(ITIME))/2.D0
      IF(DELT.GT.DUM) DELT=DUM
      DUM=TIME+DELT
      IF(DUM.GE.TIMEBC(IDUM)) THEN
        DELT=TIMEBC(IDUM)-TIME
        IEND=1
      ENDIF
C     where IEND=1 indicates the end of each ITIME computation 
C     in Main Program
C     
C     EDIKE(J,L)=downward erosion depth (m) from initial (time=0.0)
C     dike surface at time=(TIME+DELT) if ICLAY=0
      DO 400 J=1,JMAX(L)
        DSUM(J)=DSUM(J)+DELT*DSTA(J)
        IF(GRSD(J,L).GT.0.D0) THEN
          IF(DSUM(J).LT.GRS5(J,L)) THEN
            EDIKE(J,L)=GRS1(J,L)*(1.D0-DSQRT(1.D0-GRS2(J,L)*DSUM(J)))
          ELSE
            EDIKE(J,L)=GRS3(J,L)*DSUM(J)-GRS4(J,L)
          ENDIF
        ELSE
          EDIKE(J,L)=GRS3(J,L)*DSUM(J)
        ENDIF
        ZB(J,L)=ZB0(J,L)-EDIKE(J,L)
 400  CONTINUE
      ENDIF
C      
C     ECLAY(J,L)=downward clay erosion depth (m) from initial clay surface
C     below sand layer using DELT computed in Subr. 12 CHANGE if ICLAY=1
      IF(ICLAY.EQ.1) THEN
      DO 500 J=1,JMAX(L)
          IF(HP(J,L).LT.D50) THEN
              DUM=DELT*RCLAY(J,L)*DSTA(J)
              EPCLAY(J,L)=EPCLAY(J,L)+DUM
              ZP(J,L)=ZP0(J,L)-EPCLAY(J,L)
              ZB(J,L)=ZB(J,L)-DUM*FCLAY(J,L)
          ENDIF
 500  CONTINUE
      ENDIF
C     
      RETURN
      END
C     
C     -22------------------------- END OF SUBROUTINE EROSON -----------------------
C     *****************************************************************************
C     This function related to dike erosion forcing in wet and dry 
C     zone is called from Subr.22 EROSON
      FUNCTION GFDWD(R)
      DOUBLE PRECISION GFDWD,R,TR,R2,R3,ERFCC
      TR=3.D0*R
      R2=R*R
      R3=R2*R
      IF(R.GE.0.D0) THEN
        GFDWD=1.32934D0+TR+2.658681D0*R2+R3
      ELSE
        GFDWD=1.32934D0*(1.D0+2.D0*R2)*(2.D0*ERFCC(R)-1.D0)
     +  -TR-R3+(16.D0*R3+9.D0*R)*DEXP(-R2)
      ENDIF
C     Complementary error function ERFCC below Subr.06 GBXAGF
C     
      RETURN
      END
C     *****************************************************************************
C     #23##################### SUBROUTINE SRFSP ##############################
C     
C     This subroutine computes the surf similarity parameter
      SUBROUTINE SRFSP(L)
C     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NN=5000, NL=100)
C     
      COMMON /PERIOD/ TP,WKPO,ANGLE,WT(NN)      
      COMMON /BPROFL/ DXD2,DXDX,DX2,DX,XB(NN),ZB(NN,NL),FB2(NN,NL),
     + SWLDEP(NN,NL),BSLOPE(NN,NL),JMAX(NL),JSWL(NL)
      COMMON /CONSTA/ GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2
      COMMON /LINEAR/ WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN),
     +   FSX, FSY, FE, QWX, QWY
      COMMON /RUNUP/  XR,ZR,SSP,JR
      COMMON /PREDIC/ HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN)

      CO = GRAV*TP/TWOPI
      ARG = ABS(CO/CP(1)*STHETA(1))
      ARG = MIN(ARG,SQR2/2.D0) ! Arbitrary max deep water angle of 45 deg
      THETAO=DASIN(ARG)
      HRMSO = HRMS(1)*DSQRT((CP(1)*WN(1)*CTHETA(1))/
     +     (0.5D0*CO*DCOS(THETAO)))

C First guess at slope uses SWS slope
      TANB = (ZB(JR+1,L)-ZB(JR-1,L))/(XB(JR+1)-XB(JR-1))
      SSP = TANB/DSQRT(SQR2*HRMSO/(TWOPI/WKPO))
C Just to improve slope estimate, estimate Runup with Mase 1989:
C      R2P = SQR2*HRMSO*1.86D0*SSP**0.71D0
      RETURN
      END
C     
C     -23------------------------- END OF SUBROUTINE SRFSP -----------------------
