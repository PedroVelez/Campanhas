C ====================================================================
C                               ANALISIS.F
C ====================================================================
C
C    THIS CODES PERFORMS A UNIVARIATE, 2D SPATIAL ANALYSIS OF STATION
C    DATA BASING ON TWO DIFFERENT INTERPOLATION TECHNIQUES:
C
C   - OPTIMAL STATISTICAL INTERPOLATION (OSI)
C   - SUCCESSIVE CORRECTIONS approaching assimptotically OSI (SC_OSI)
C
C ====================================================================
C           
C                         AUTHOR: Damià Gomis*
C
C       IMEDEA (Institut Mediterrani d'Estudis AvanÇats), a joint 
C     research centre between the Universitat de les Illes Balears 
C    (UIB) and the Spanish Research Council (CSIC). Mallorca (Spain)
C
C    This code has been revised and distributed in the framework of
C    the project REN2000-2599-E, funded by the Spanish Marine Science
C    and Technology subprogram. 
C  
C    * Some of the routines were originally developed by Dr. Mike A. 
C      Pedder (Dep. of Meteorology, University of Reading, UK.)
C
C ====================================================================
C
C   Work structure:
C   - the user selects between carrying out the analysis just at one 
C     level (to be determined by the user) or at the whole set of 
C     levels specified in the domain data file 'grid.dat'
C   - the analysis domain is read from the external file 'grid.dat'
C   - the grid point coordinates are computed 
C   - the level data are read from the corresponding '*.lev' file
C   - the analysis routines are called 
C   - the output values are written in the output file
C
C   External files used:  
C   - INPUT files:
C     * domain info file:                           [../info/grid.dat]
C     * level data files:                        [../lev/NNNNLLLL.lev]
C   - OUTPUT files: 
C     * gridpoint data files:                  [../snp/NNNNLLLLVV.snp] 
C                                           or [../grd/NNNNLLLLVV.grd] 
C
C   Built-in routines:  SUCCOR.F, OSI.F (include lower level routines)
C   External routines:  none
C
C ====================================================================
C
C   ROUTINE STRUCTURE:
C
C                                 ANALISIS
C                ____________________|____________________
C               |                                         |
C               |                                         |
C             SUCCOR _____________ DETREN ______________ OSI..........
C               |                    |                    |          .
C          _____|____            ____|____           _____|____      .
C         |          |          |         |         |          |     .
C   ....SC_OSI     SC_STA     DTRND2   POLINOM    LUBKSB     EIGEN   .
C   .     |          .          |         |       LUDCMP       |     . 
C   .     |          .        PBASE0   ROUT2-7              JACOBI   .
C   .    SCR0        .                                               .
C   .    SCR1        .                                               .
C   .                ........>                                       .
C   .........................> FUNCTION CRMDL <.......................
C
C
C ====================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,NLMAX=100)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      CHARACTER*11 LLEV(2),LNDR(2)
      CHARACTER*8 SYMBOL(NS)
	CHARACTER*8 NAME,STUDA,NA(NLMAX)
      CHARACTER*4 IFOR(2),NA0,NA1
      CHARACTER*2 ICVAR(9)
      CHARACTER*1 NOP,NZON
      DATA IFOR/'.grd','.snp'/
      DATA ICVAR/'tc','sa','st','dh','th','uu','vv','pr','mp'/
      LLEV(1)='     level:'
      LNDR(1)=' ref.level:'
      PI=3.14159
      RD=PI/180.
      XKGLA=60.*1.852
C
    1 WRITE(6,2)
    2 FORMAT(7(/),1X,'####  UNIVARIATE, 2D SPATIAL ANALYSIS OF ',
     .'OCEANOGRAPHIC DATA  ####'///)
C
C ====================================================================
C                           INPUT DATA BLOCK
C ====================================================================
c
      open(unit=99,file='../info/grid.dat',status='old')
      read(99,3)
      read(99,3)
      read(99,*) XLON1,XLAT1
      read(99,3)
      read(99,*) ALF0
      read(99,3) 
      read(99,*) NL,NF,NC
      read(99,3)
      read(99,*) XARM,YARM
      read(99,3)
      read(99,*) P0,PINT
	read(99,3)
	read(99,*) NA0,STUDA
    3 format(1x/1x)
      close(unit=99)
c
      ALF0=ALF0*RD
      NGRI=NF*NC
C     Coordinates of the other corners (degrees):
    5 XLMAX=XARM*REAL(NC-1)
      YLMAX=YARM*REAL(NF-1)
      WRITE(6,7) NC,NF
    7 FORMAT(/1X,'*** The domain, as it has been read from the external'
     .,' file,'/5X,'has been set to ',I3,' X ',I3,' GRID POINTS, with ',
     .'grid cells')
      IF(ABS(ALF0).GT.0.1) THEN
         XLAT4=XLAT1+XLMAX*DSIN(ALF0)/XKGLA
         XKGLO=XKGLA*DCOS(.5*RD*(XLAT1+XLAT4))
         XLON4=XLON1+XLMAX*DCOS(ALF0)/XKGLO
         XLAT3=XLAT4+YLMAX*DSIN(PI/2.+ALF0)/XKGLA
         XKGLO=XKGLA*DCOS(.5*RD*(XLAT3+XLAT4))
         XLON3=XLON4+YLMAX*DCOS(PI/2.+ALF0)/XKGLO
         XLAT2=XLAT1+YLMAX*DSIN(PI/2.+ALF0)/XKGLA
         XKGLO=XKGLA*DCOS(.5*RD*(XLAT1+XLAT2))
         XLON2=XLON1+YLMAX*DCOS(PI/2.+ALF0)/XKGLO
         WRITE(6,9) XARM,YARM,XLAT2,XLON2,XLAT3,XLON3,XLAT1,XLON1,
     .	   XLAT4,XLON4
      ELSE
         XLON4=XLON1+XLMAX
         XLAT4=XLAT1
         XLON3=XLON4
         XLAT3=XLAT4+YLMAX
         XLON2=XLON1
         XLAT2=XLAT1+YLMAX
         WRITE(6,11) XARM,YARM,XLAT2,XLON2,XLAT3,XLON3,XLAT1,XLON1,
     .		XLAT4,XLON4
      ENDIF
    9 FORMAT(5x,'of ',f5.1,'  x ',f5.1,'  km. The LIMITS of the domain',

     .' are:'//2(2(5x,'lat =',f7.3,'  lon =',f8.3)/))
   11 FORMAT(5x,'of',f5.3,'  x ',f5.3,' deg. The LIMITS of the domain',

     .' are:'//2(2(5x,'lat =',f7.3,'  lon =',f8.3)/))
C
      WRITE(6,'(a,$)') 
     .' >>> DO YOU WANT TO USE *THIS* DOMAIN AND GRID ? ([Y]/N) '
      READ(5,13) NOP
   13 FORMAT(A1)
      IF(NOP.NE.'N') GOTO 30
C ----------------------------------------------------------------------
C     DEFINE A NEW ZONE AND A NEW GRID:
   15 WRITE(6,17)
   17 FORMAT(2(/),1X,'***  DEFINITION OF A NEW DOMAIN AND GRID:'/)
   20 WRITE(6,'(a,$)')
     .' >>> Coordinates of the low-left corner  [LAT1,LON1 (deg)]: '
      READ(5,*,ERR=20) XLAT1,XLON1
   23 WRITE(6,'(a,$)')  ' >>> Slope angle  [ALF0 (deg)]: '
      READ(5,*,ERR=23) ALF0
      ALF0=ALF0*RD
   25 WRITE(6,'(a,$)') ' >>> Number of Grid columns and rows  [NC,NR]: '
      READ(5,*,ERR=25) NC,NF
   27 IF(ABS(ALF0).GT.0.1) THEN 
         WRITE(6,'(a,$)') ' >>> Size of grid cells  [XARM,YARM (km)]: '
      ELSE
         WRITE(6,'(a,$)') ' >>> Size of grid cells  [XARM,YARM (deg)]: '
      ENDIF
      READ(5,*,ERR=27) XARM,YARM
      GOTO 5
C ----------------------------------------------------------------------
C     COMPUTES THE GRID POINT COORDINATES
C     The location of grid points bases on a linear transformation
C     between distance and degrees (as given by a rotation matrix). 
C     The common origin is the low-left corner of the domain.
C     The grid point index are: point(1) corresponds to the up-left 
C     corner; the up-right corner is (NC), and so on until the low-
C     right corner (NF*NC=NGRI). 
   30 DO 33 I=1,NF
      DO J=1,NC
      	 K=NC*(NF-I)+J
      	 IF(ABS(ALF0).GT.0.1) THEN
            YL=YARM*REAL(I-1)
            XL=XARM*REAL(J-1)
            D=DSQRT(XL**2+YL**2)
            AL=DATAN2(YL,XL)
            YG(K)=(XLAT1+D*DSIN(AL+ALF0)/XKGLA)*RD
            XKGLO=XKGLA*DCOS(.5*(XLAT1*RD+YG(K)))
            XG(K)=(XLON1+D*DCOS(AL+ALF0)/XKGLO)*RD
         ELSE
            YG(K)=(XLAT1+YARM*REAL(I-1))*RD
            XG(K)=(XLON1+XARM*REAL(J-1))*RD
         ENDIF
      ENDDO
   33 CONTINUE
      IF(NZON.EQ.'Y') GOTO 100
C
C ======================================================================
C                          READING STATION DATA 
C ======================================================================
      WRITE(6,*)
	WRITE(6,31) 
   31 FORMAT('*** The code can operate on a single level (op. 1)
     . or for'/'*** the whole set of levels (op.2)')  
   34 WRITE(6,36)  
   36	FORMAT(/1X,'SELECT OPTION: SINGLE LEVEL(1), ALL LEVELS(2):')
	READ(5,*,ERR=34) REP
	IND=0   
   35	IF(REP.EQ.1) THEN
      WRITE(6,37)
   37 FORMAT(/1X,'*** STATION DATA WILL BE READ FROM AN EXTERNAL FILE'/)
   40 WRITE(6,'(a,$)')' >>> INPUT LEVEL (METERS) OF LEV FILE: [LLLL]: '
      READ(5,42,ERR=40) NA1
   42 FORMAT(A4)
      IL=1
      NA(IL)=NA0//NA1
C
      OPEN(UNIT=2,FILE='../level/'//NA(IL)//'.lev',STATUS='OLD')
      READ(2,44) LLEV(2),LNDR(2)
   44 FORMAT(//2(/30X,A11)//////)
	GOTO 45
      ENDIF

C     For reading all the files
         PMAX=P0+PINT*FLOAT(NL-1)
         DO I=1,NL
            LEV=(P0+PINT*FLOAT(I-1)+0.00001)
            NA1=CHAR(LEV/1000+48)
     .	//CHAR(LEV/100-10*(LEV/1000)+48)//
     .	CHAR(LEV/10-10*(LEV/100)+48)//CHAR(LEV-10*(LEV/10)+48)
          NA(I)=NA0//NA1
         ENDDO
	 
   45	WRITE(6,46)
   46 FORMAT(//1X,'*** VARIABLES:   1: temp   2: sal   3: st    4: dh',
     ./18X,'5: th 6: u      7: v     8: pres  9: mp')  
   47 WRITE(6,'(a,$)') ' >>> CHOOSE variable: '
      READ(5,*,ERR=47) NVAR
	

      IL=1
   48 NND=0
      NSR=0
      NFZ=0
      ND=1
	IF(REP.EQ.1) GOTO 49
      OPEN(UNIT=2,FILE='../level/'//NA(IL)//'.lev',STATUS='OLD')
      READ(2,44) LLEV(2),LNDR(2)
     
   49	DO 75 I=1,NS
         GOTO (50,52,54,56,56,58,60,62,56) NVAR
            GOTO 48
   50       READ(2,80,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   80       FORMAT(1x,A8,1X,F7.3,1X,F8.3,1x,f6.3)
            GOTO 63
   52       READ(2,82,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   82       FORMAT(1x,A8,1X,F7.3,1x,f8.3,8x,f6.3)
            GOTO 63
   54       READ(2,84,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   84       FORMAT(1x,A8,1x,f7.3,1x,f8.3,15x,f6.3)
            GOTO 63
   56       READ(2,86,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   86       FORMAT(1x,A8,1x,f7.3,1x,f8.3,22x,f7.3)
            GOTO 63
   58       READ(2,88,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   88       FORMAT(1x,A8,1X,F7.3,1x,f8.3,30x,f6.1)
            GOTO 63
   60       READ(2,90,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   90       FORMAT(1x,A8,1X,F7.3,1x,f8.3,37X,F6.1)
            GOTO 63
   62       READ(2,92,END=77) SYMBOL(ND),Y(ND),X(ND),V(ND)
   92       FORMAT(1x,A8,1X,F7.3,1x,f8.3,45X,a9)

C        The presence of repeated stations, missing data or located out  
C        of the domian is checked (the latter are used in the analysis, 
C        but they are counted anyway...)
   63    IF(NVAR.LE.3) THEN
C           the value of 'no data' for T, S and St is 99.999 
            IF(V(ND).LT.99.) GOTO 67
         ELSE IF((NVAR.EQ.4).OR.(NVAR.EQ.5).OR.(NVAR.EQ.8)) THEN
C           the value of 'no data' for the geopotencial is 999.999
            IF(V(ND).LT.999.) GOTO 67
         ELSE IF((NVAR.EQ.6).OR.(NVAR.EQ.7)) THEN
C           the value of 'no data' for the current is 9999.9
            IF(V(ND).LT.9999.) GOTO 67
         ELSE IF(NVAR.EQ.8) THEN
C           the value of 'no data' for pressure is -99.99
            IF(V(ND).GT.-99.) GOTO 67
         ENDIF
         NND=NND+1
         GOTO 75           
   67    IF(ABS(ALF0).GT.0.1) THEN
            YL=(Y(ND)-XLAT1)*XKGLA
            XKGLO=XKGLA*DCOS(.5*RD*(XLAT1+Y(ND)))
            XL=(X(ND)-XLON1)*XKGLO
            D=DSQRT(XL**2+YL**2)
            AL=DATAN2(YL,XL)
            YL=D*DSIN(AL-ALF0)
            XL=D*DCOS(AL-ALF0)
            IF((XL.LT.0).OR.(XL.GT.XLMAX)) GOTO 70
            IF((YL.GE.0).AND.(YL.LE.YLMAX)) GOTO 73
         ELSE
            IF((X(ND).LT.XLON1).OR.(X(ND).GT.XLON4)) GOTO 70
            IF((Y(ND).GE.XLAT1).AND.(Y(ND).LE.XLAT2)) GOTO 73
         ENDIF
   70    NFZ=NFZ+1
   73    Y(ND)=Y(ND)*RD
         X(ND)=X(ND)*RD
         ND=ND+1
   75 CONTINUE
   77 NSTA=ND-1
      CLOSE(UNIT=2)
	WRITE(6,93) NA(IL)
      WRITE(6,94) (I-1)
      IF(NOP.NE.'N') WRITE(6,96) NND,NSR,NFZ
      WRITE(6,98) NSTA
   93	FORMAT(//5X,'#######   FILE ',A8,'.LEV          #########')
   94 FORMAT(5X,'#######   READ ',I3,' STATIONS FROM THE FILE  ###')
   96 FORMAT(5X,'### ELIMINATED ',I3,' STATIONS WITH NO DATA   ###'/
     .5X,'### ELIMINATED ',I3,' REPEATED STATIONS       ###'/
     .5X,'### ACCEPTED   ',I3,' STATIONS OUT OF DOMAIN  ###')
   98 FORMAT(5X,'### A TOTAL OF ',I3,' VALUES WILL BE ANALYSED ###'//)

      IF((REP.EQ.2).AND.(NRUT.EQ.1).AND.(IND.EQ.1)) THEN
	CALL OSI(nsta,ngri,rep,ind)
	GOTO 113
	ENDIF

	IF((REP.EQ.2).AND.(NRUT.EQ.2).AND.(IND.EQ.1)) THEN
	CALL SUCCOR(NRUT,NSTA,NGRI,REP,IND)
	GOTO 113
	ENDIF
C
C ======================================================================
C                    CALLING THE INTERPOLATION ROUTINES
C ======================================================================
  100 WRITE(6,102)
  102 FORMAT(//,1X,'*** INTERPOLATION ROUTINES AVAILABLE:'/5X,'1- OPTI',
     .'MAL STATISTICAL INTERPOLATION (OSI)'/5X,
     .'2- SUCCESSIVE CORRECTIONS (->OSI)')
  105 WRITE(6,'(a,$)') ' >>> Choose option: '
      READ(5,*,ERR=105) NRUT
      GOTO (107,110) NRUT
  107    CALL OSI(NSTA,NGRI,REP,IND)
         GOTO 113
  110    CALL SUCCOR(NRUT,NSTA,NGRI,REP,IND)
C      Mean Residuals at Station points:
  113 rms=0.
      DO k=1,nsta
         rms=rms+(v(k)-vs(k))**2
      ENDDO
      rms=sqrt(rms/real(nsta))
      WRITE(6,115) rms
  115 FORMAT(//,1X,'*** INTERPOLATION FINISHED. Rms at Stations: ',f8.4)
C
C ======================================================================
C                       WRITING DATA TO OUTPUT FILE
C ======================================================================
	IF((REP.EQ.2).AND.(IND.EQ.1)) GOTO 126
	WRITE(6,120)
  120 FORMAT(//,1X,'*** WRITING OUTPUT FILES: '/)
  123 WRITE(6,'(a,$)') ' >>> Choose Format option:  1-SURFER  2-SNP : '
      READ(5,*,ERR=123) NSOR
         WRITE(6,125)
  125    FORMAT(1X,'>>> Save complementary information (about the ',
     .   'domain'/5X,'analysis parameters, station values,etc) below')
         WRITE(6,'(a,$)') 
     .'     the data body of the output file (Y/[N]): '
         READ(5,13) NOP
C
  126   OPEN(unit=3,file='../varobs/'//na(IL)//icvar(nvar)//ifor(nsor),
     .                                          status='unknown')
      IF(NSOR.EQ.1) THEN
C        Data plotted in SURFER (.grd) format
         IF(alf0.gt.0.) THEN
            xmin=0.
            ymin=0.
         ELSE
            xmin=xlon1
            ymin=xlat1
         ENDIF
         zmin=1E6
         zmax=-1E6
         DO k=1,ngri
            IF(vg(k).gt.zmax) zmax=vg(k)
            IF(vg(k).lt.zmin) zmin=vg(k)
         ENDDO
         WRITE(3,130) nc,nf,xmin,xlmax+xmin,ymin,ylmax+ymin,zmin,zmax
         DO i=nf,1,-1
            WRITE(3,135) (vg(k),k=(nc*(i-1)+1),(nc*i))
C            WRITE(3,137)
         ENDDO
  130    FORMAT('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f9.4))
  135    FORMAT(10f10.4)
  137    FORMAT(1x)
		IF(NOP.EQ.'Y') THEN
            WRITE(3,145) LLEV,NSTA,LNDR
            WRITE(3,147) XLAT2,XLON2,XLAT3,XLON3,XLAT1,XLON1,XLAT4,
     .                                                       XLON4
            IF (ABS(ALF0).GT.0.1) THEN
               WRITE(3,150) NF,NC,XARM,YARM
            ELSE
               WRITE(3,152) NF,NC,XARM,YARM
            ENDIF
C           Specific information about the method/parameters used:
            IF(NRUT.EQ.1) WRITE(3,155) (PAR(I),I=1,4)
            IF(NRUT.EQ.2) WRITE(3,160) (PAR(I),I=1,5)
C           IF(NRUT.EQ.3) WRITE(3,160) (PAR(I),I=1,5)
C           Values at station points:
            WRITE(3,163)
            DO I=1,NSTA
               WRITE(3,165) SYMBOL(I),Y(I)/RD,X(I)/RD,V(I),VS(I)
            ENDDO
            WRITE(3,167) RMS
		ENDIF
      ELSE IF(NSOR.EQ.2) THEN
c        Data plotted in SNP format 
         DO I=(NF-1)*NC,0,-NC
            WRITE(3,140) (VG(I+J),J=1,NC)
         ENDDO
  140    FORMAT(8F10.4)
C
         IF(NOP.EQ.'Y') THEN
            WRITE(3,145) LLEV,NSTA,LNDR
            WRITE(3,147) XLAT2,XLON2,XLAT3,XLON3,XLAT1,XLON1,XLAT4,
     .                                                       XLON4
            IF (ABS(ALF0).GT.0.1) THEN
               WRITE(3,150) NF,NC,XARM,YARM
            ELSE
               WRITE(3,152) NF,NC,XARM,YARM
            ENDIF
C           Specific information about the method/parameters used:
            IF(NRUT.EQ.1) WRITE(3,155) (PAR(I),I=1,4)
            IF(NRUT.EQ.2) WRITE(3,160) (PAR(I),I=1,5)
C           IF(NRUT.EQ.3) WRITE(3,160) (PAR(I),I=1,5)
C           Values at station points:
            WRITE(3,163)
            DO I=1,NSTA
               WRITE(3,165) SYMBOL(I),Y(I)/RD,X(I)/RD,V(I),VS(I)
            ENDDO
            WRITE(3,167) RMS
  145       FORMAT(///19X,2A11,'  contains ',I3,' valid data'/19X,2A11/)
  147       FORMAT(1X,'Domain:',2(5X,'lat =',F6.3,'  lon =',F7.3)/8X,
     .      2(5X,'lat =',F6.3,'  lon =',F7.3)/)
  150       FORMAT(1X,'Grid:',7X,I2,' X ',I2,' points, with XARM =',
     .      F5.1,' km, YARM =',F5.1,' km')
  152       FORMAT(1X,'Grid:',7X,I2,' X ',I2,' points, with XARM =',
     .      F5.3,' deg, YARM =',F5.3,' deg')
  155       FORMAT(//1X,'Analysis technique: OSI. Parameters used:'/
     .      1X,'POLINOM.DEG=',F2.0,'  SCL=',F4.0,'  GAMMA=',F6.4,
     .      '  XCUT=',F4.0)
  157       FORMAT(//1X,'Analysis technique: SC-Standard. Parameters',
     .      ' used:'/1X,'DETR.DEGREE=',F2.0,'  SCL0=',F4.0,'  G=',
     .      F4.2,'  NITER=',F4.0)
  160       FORMAT(//1X,'Analysis technique: SC->OSI. Parameters ',
     .      'used:'/1X,'DETR.DEGREE=',F2.0,'  SCL=',F4.0,'  GAMMA=',
     .      F6.4,'  XCUT=',F4.0,'  NITER=',F4.0)
  163       FORMAT(/1x,'Data:'//1X,'SYMBOL    LATITUD  LONGITUD   ',
     .      'OBS.DATA   ANA.DATA')
  165       FORMAT(1X,A8,2(3X,F7.3),2(3X,F8.3))
  167       FORMAT(/1X,'Mean Residual at stations: ',f6.3)
         ENDIF
      ELSE
	  GOTO 123
      ENDIF
	WRITE(6,180) NA(IL),ICVAR(NVAR),IFOR(NSOR)
      CLOSE(UNIT=3)
	
C
C ======================================================================
C                           RE-STARTING OPTIONS
C ======================================================================
C     The possibility to restart the program changing some of the 
C     arguments is offered (grid, data file, method,...)

      WRITE(6,175)
  175 FORMAT(//1X,'____________________ END OF THE PROCESS ',
     .'____________________'//)
	
	IF(REP.EQ.2) THEN
	 IL=IL+1 
	 IND=1	 
	 IF(IL.LE.NL) GOTO 48
	ENDIF
	
	IF(REP.EQ.2) GOTO 185
      
	WRITE(6,'(a,$)') ' >>> Change the domain or the grid ?  (Y/[N])'
      READ(5,13) NZON
      IF(NZON.EQ.'Y') GOTO 15
      WRITE(6,*)
      WRITE(6,'(a,$)') ' >>> Change Station data ?  (Y/[N]) '
      READ(5,13) NOP
      IF(NOP.EQ.'Y') GOTO 35
	WRITE(6,*)
	WRITE(6,'(a,$)') ' >>> Repeat the analysis process ?  (Y/[N]) '
      READ(5,13) NOP
      IF(NOP.EQ.'Y') GOTO 100
      WRITE(6,*)
      WRITE(6,'(a,$)') ' >>> Restart the whole program ?  (Y/[N]) '
      READ(5,13) NOP
      IF(NOP.EQ.'Y') GOTO 1
C     TYPE 180,NA(IL),ICVAR(NVAR),IFOR(NSOR)
  180 FORMAT(5(/),8X,'============================================'//
     .18X,'NAME OF THE OUTPUT FILE:'/23X,A8,A2,a4//
     .8X,'============================================',7(/))
  185 END



C ======================================================================
C                     SUBROUTINE SUCCOR (NRUT,NSTA,NGRI)
C
C   Work structure:
C   - detrends the data
C   - calls routines to compute interpolation values (2 options)
C   - returns to the ANÁLISIS.F
C
C   External files used: none
C
C   Built-in routines: SC_STA,SC_OSI
C   External routines: DETREN.F 
C
C ======================================================================
      SUBROUTINE SUCCOR(NRUT,NSTA,NGRI,REP,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),FIN(NS),DIN(NS),SG(ND,NG),TR(NG)
C ----------------------------------------------------------------------
      WRITE(6,3)
    3 FORMAT(10(/),1X,'############# SUCCESSIVE CORRECTIONS ANALYSIS ',
     .' METHOD  #############'//)
C
C     PRELIMINARY STEPS: 
      DO K=1,NSTA
         SS(1,K)=X(K)
         SS(2,K)=Y(K)
         FIN(K)=V(K)
      ENDDO
      DO K=1,NGRI
         SG(1,K)=XG(K)
         SG(2,K)=YG(K)
      ENDDO

	IF((REP.EQ.2).AND.(IND.EQ.1)) THEN
	CALL DETREN(NSTA,NGRI,NP,RMS)
	WRITE(6,11) RMS
	GOTO 12
	ENDIF
C
C     DETRENDING:
      IF(NRUT.EQ.3) THEN
         WRITE(6,5)
    5    FORMAT(1X,'*** The STANDARD VERSION of SC can work without ',
     .   'detrending'/5X,'the data. Go ahead...')
         PAR(1)=99.
      ELSE
         WRITE(6,7)
    7    FORMAT(1X,'*** The version of SC approaching OSI requires a',
     .   ' previous'/5X,'detrending of the data. This will be made e',
     .   'fective by'/5X,'fitting a polinomial to the station data.'/)
    9    WRITE(6,'(a,$)') ' >>> INPUT POLYNOMIAL DEGREE [NP<=5]: '
         READ(5,*,ERR=9) NP
         PAR(1)=NP
   10     CALL DETREN(NSTA,NGRI,NP,RMS)
         WRITE(6,11) RMS
   11    FORMAT(/1X,'*** Detrending finished. Station Residuals:  RMS=',
     .   F6.2/)
      ENDIF
C ----------------------------------------------------------------------
C     INFORMATION ON THE WEIGHT FUNCTIONS:
   12	IF((REP.EQ.2).AND.(IND.EQ.1)) THEN
      CALL SC_OSI(NSTA,NGRI,SCL,GAMMA,FSCL,CF,NITER)
C
C        Adding the trend to the interpolated values:
         DO K=1,NGRI
            VG(K)=VG(K)+TR(K)
         ENDDO
         DO K=1,NSTA
            VS(K)=VS(K)+FIN(K)
         ENDDO
	GOTO 34
	ENDIF

      WRITE(6,13)
   13 FORMAT(/1X,'*** Distance Weight Functions are taken as gausians:'
     ./14X,'W = EXP(-d**2/2*Sn**2)'/5X,'where Sn is the characteristic',
     .'scale used in the'/5X,'n-iteration and d is distance.'/)
   15 WRITE(6,'(a,$)') ' >>> Choose number of iterations  [NITER<100]: '
      READ(5,*,ERR=15) NITER
C
C ======================================================================
C                              SC - STANDARD
C ======================================================================
      IF(NRUT.EQ.3) THEN
         WRITE(6,17)
   17    FORMAT(/1X,'*** For the STANDARD VERSION of SC, Sn is taken ', 
     .   'as'/5X,'Sn=G*S(n-1), with 0<G<1. So you must:')
   20    WRITE(6,'(a,$)') ' >>> Input values of [S1(km),G] : '
         READ(5,*,ERR=20) SCL,G
         PAR(2)=SCL
         PAR(3)=G
         PAR(4)=NITER
         CALL SC_STA(NSTA,NGRI,SCL,G,NITER)
C ======================================================================
C                                 SC-OSI
C ======================================================================
      ELSE IF(NRUT.EQ.2) THEN
         WRITE(6,23)
   23    FORMAT(/1X,'*** For the VERSION of SC that approaches OSI, ',
     .   'Ln are'/5X,'the same for any iteration. So you must only:')
   25    WRITE(6,'(a,$)') ' >>> Input (constant) value of [S(km)] : '
         READ(5,*,ERR=25) SCL
         PAR(2)=SCL
         WRITE(6,28)
   28    FORMAT(/1X,'*** Additionally, you must input:')
   30    WRITE(6,'(a,$)') 
     .   ' >>> The noise to signal variance ratio  [gamma] : '
         READ(5,*,ERR=30) GAMMA
         PAR(3)=GAMMA
         WRITE(6,32)
   32	   FORMAT('>>> The cut-off wavelength [XCUT(km)] of the spatial')
   33    WRITE(6,'(a,$)') 
     .   '     smoothing (set to 0 for no additional smoothing): '
         READ(5,*,ERR=33) XCUT
         FSCL=XCUT/4.
         PAR(4)=XCUT
         PAR(5)=NITER
C        Over-relaxation Parameter set to 1 by default:
         CF=1.
C
        CALL SC_OSI(NSTA,NGRI,SCL,GAMMA,FSCL,CF,NITER)
C
C        Adding the trend to the interpolated values:
         DO K=1,NGRI
            VG(K)=VG(K)+TR(K)
         ENDDO
         DO K=1,NSTA
            VS(K)=VS(K)+FIN(K)
         ENDDO
      ENDIF
C
C ======================================================================
C
   34 RETURN
      END


C ======================================================================
C                 SUBROUTINE SC_STA(NSTA,NGRI,SCL,G,NITER)
C
C   Work structure:
C     - SC scheme with variable weight functions for each iteration. 
C     - Weights are normalized, since the variable mean is NOT zero.
C     - The algorithm is the traditional one, based on the analysis 
C       updating at both, grid points and stations.
C
C   External files used: none
C
C   Built-in routines:  none
C   External routines:  none
C
C   External Functions: CRMDL
C
C ======================================================================
      SUBROUTINE SC_STA(NSTA,NGRI,SCL,G,NITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),FIN(NS),DIN(NS),SG(ND,NG),TR(NG)
      DIMENSION WSS(NS,NS),WGS(NG,NS),DSS(NS),DGS(NG),UPINC(NS)
C
C ======================================================================
C                 Estimate of the initial weight funtions
C ======================================================================
      WRITE(6,5)
    5 FORMAT(/1X,'...Computing initial weight matrices...')
C     No normalized grid-station weight matrix WGS, and station-station 
C     weight matrix WSS:
      DO K=1,NSTA
         WSS(K,K)=1.
         DSS(K)=1.
         DO L=1,K-1
            WSS(K,L)=CRMDL(SS(1,K),SS(1,L),ND,SCL)
            WSS(L,K)=WSS(K,L)
            DSS(K)=DSS(K)+WSS(K,L)
            DSS(L)=DSS(L)+WSS(K,L)
         ENDDO
      ENDDO
C
C ....Weight Grid<--Stations
      DO K=1,NGRI
         DGS(K)=0.
         DO L=1,NSTA
            WGS(K,L)=CRMDL(SG(1,K),SS(1,L),ND,SCL)
            DGS(K)=DGS(K)+WGS(K,L)
         ENDDO
      ENDDO
C
C ======================================================================
C                            Starts iterating
C ======================================================================
C     Iterations process is the traditional one:
C     The analysis after n iteratios is given by:
C                    Vg(n) = Vg(n-1) + Wgs(n)*[Os-Vs(n-1)]
C                    Vs(n) = Vs(n-1) + Wss(n)*[Os-Vs(n-1)]
C     where Os are station observations, and Wgs(n), Wss(n) are the 
C     weight funtion matrices of iteration n.
C
C     Vector UPINC(K) is used to update the increments, and vectors 
C     VG(K) and VS(K) to update analysis:
C
      WRITE(6,*) '...Starting with the iterations...'
      DO K=1,NSTA
         UPINC(K)=V(K)
         VS(K)=0.
      ENDDO
      DO K=1,NGRI
         VG(K)=0.
      ENDDO
C
      FAC=1./G**2
      DO 100 IT=1,NITER
         IF(IT.EQ.1) GOTO 80
C ...... Updating of Wgs, Wss, (Dgs and Dss)
         DO 70 K=1,NSTA
            WSS(K,K)=1.
            DSS(K)=1.
            DO L=1,K-1
               IF(WSS(K,L).GT.0) THEN
                  WSS(K,L)=WSS(K,L)**FAC
                  WSS(L,K)=WSS(K,L)
                  DSS(K)=DSS(K)+WSS(K,L)
                  DSS(L)=DSS(L)+WSS(K,L)
               ENDIF
            ENDDO
   70    CONTINUE
         DO 75 K=1,NGRI
            DGS(K)=0.
            DO L=1,NSTA
               IF(WGS(K,L).GT.0) THEN
                  WGS(K,L)=WGS(K,L)**FAC
                  DGS(K)=DGS(K)+WGS(K,L)
               ENDIF
            ENDDO
   75    CONTINUE
C ...... Analysis update:
C        ATTENTION: NEVER update Vs before Vg, because Vs(n-1) is 
C        used to update Vg(n).
C        Estimate of the analysed fields variances..
   80    VAR=0.
         DO 85 K=1,NGRI
            SUM=0.
            DO L=1,NSTA
               SUM=SUM+WGS(K,L)*UPINC(L)
            ENDDO
            VG(K)=VG(K)+SUM/DGS(K)
            VAR=VAR+(SUM/DGS(K))**2
   85    CONTINUE
         VAR=DSQRT(VAR/REAL(NGRI))
         DO 90 K=1,NSTA
            SUM=0.
            DO L=1,NSTA
               SUM=SUM+WSS(K,L)*UPINC(L)
            ENDDO
            VS(K)=VS(K)+SUM/DSS(K)
   90    CONTINUE
C        Estimate of station residuals:
         RMS=0.
         DO K=1,NSTA
            UPINC(K)=V(K)-VS(K)
            RMS=RMS+UPINC(K)**2
         ENDDO
         RMS=DSQRT(RMS/REAL(NSTA))
         WRITE(6,110) IT,VAR,RMS
  100 CONTINUE
  110 FORMAT(5X,'iteration ',i3,': - variance corr.(grid): Sto=',F6.2
     ./20X,'- station residuals: RMS=',F6.2)
C
      RETURN
      END



C ======================================================================
C           SUBROUTINE SC_OSI(NSTA,NGRI,SCL,GAMMA,FSCL,CF,NITER)
C
C   Work structure:
C     - SC scheme with constant weight functions for each iteration. 
C     - Weights are not normalized (the variable mean is always zero)
C     - The algorithm is such that it approaches OSI after infinite 
C       iterations
C
C   External files used: none
C
C   Built-in routines:  SCR0,SCR1
C   External routines:  none
C
C   External Functions: CRMDL
C
C ======================================================================
      SUBROUTINE SC_OSI(NSTA,NGRI,SCL,GAMMA,FSCL,CF,NITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),FIN(NS),DIN(NS),SG(ND,NG),TR(NG)
      COMMON/SCOSI/B(NS),WK(5*NS),LST(NS*NS),WT(NS*NS)
C
C ======================================================================
      WRITE(6,5)
    5 FORMAT(/1X,'...Calling SCR0...')
      CALL SCR0(NSTA,SCL,CF,GAMMA)
      WRITE(6,*) '...Calling SCR1...'
      DO I=1,NITER
         CALL SCR1(NSTA,GAMMA,RMS)
         WRITE(6,7) I,RMS
      ENDDO
    7 FORMAT(5X,'iteration ',i3,': station residual:  Sto=',F10.6)
C     On exit, the set of parameters B(NSTA) is equivalent to the 
C     product of [C_ss]**_1 * V'
C
      WRITE(6,10)
   10 FORMAT(/1X,'...Computing interpolated values:')
C     Gridpoint values:
      WRITE(6,*) '    - at Grid points...'
      d1=dsqrt(scl**2+fscl**2)
      d2=dsqrt(scl**2+2.*fscl**2)
      c1=(scl/d1)**2
      c2=(scl/d2)**2
      DO k=1,ngri
         vg(k)=0.
         DO j=1,nsta
            weight=2.*c1*CRMDL(SG(1,K),SS(1,J),ND,D1)
     .               -c2*CRMDL(SG(1,K),SS(1,J),ND,D2)
            vg(k)=vg(k)+weight*b(j)
         ENDDO
      ENDDO
C     Station values:
      WRITE(6,*) '    - at Station points...'
      DO k=1,nsta
         vs(k)=0.
         DO j=1,nsta
            weight=2.*c1*CRMDL(SS(1,K),SS(1,J),ND,D1)
     .               -c2*CRMDL(SS(1,K),SS(1,J),ND,D2)
            vs(k)=vs(k)+weight*b(j)
         ENDDO
      ENDDO
C
      RETURN
      END


      SUBROUTINE SCR0(NZ,SCL,CF,EF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/DET/S(ND,NS),FIN(NS),Z(NS),SG(ND,NG),TR(NG)
      COMMON/SCOSI/B(NS),WK(5*NS),LST(NS*NS),WT(NS*NS)
C
      IWK=5*NS
      IWT=NS*NS
      IL=NS*NS
C     Routine to setup weighting factors and initialize parameters
C     ready for successive correction iterations (of SCR1).
C
      IF (IWK.LT.2*NZ) THEN
      WRITE (6,'(1X,A60)') 'SCR0 FAILS . Parameter IWK less than 2*NZ.'
      STOP
      ENDIF
C     IF (IS.LT.ND) THEN
C     WRITE (6,'(1X,A60)') 'SCR0 FAILS . Parameter IS less than ND.'
C     STOP
C     ENDIF
C
      EPS=1.0E-6
      NWT=0
      NWK=0
      DO 10 JK=1,NZ
      WK(JK)=1.0+EF
      JLOC=NWK+1
      NWK=JLOC
      NJ=0
      IF (JK.EQ.1) GOTO 12
      DO 11 IK=1,JK-1
      TEST=CRMDL(S(1,JK),S(1,IK),ND,SCL)
      IF (TEST.LT.EPS) GOTO 11
      NJ=NJ+1
      NWK=NWK+1
      NWT=NWT+1
      IF (NWK.GT.IL.OR.NWT.GT.IWT) GOTO 11
      LST(NWK)=IK
      WT(NWT)=TEST
      WK(JK)=WK(JK)+WT(NWT)
      WK(IK)=WK(IK)+WT(NWT)
   11 CONTINUE
   12 CONTINUE
      IF (JLOC.LE.IL) LST(JLOC)=NJ
   10 CONTINUE
C
      IF (NWK.GT.IL.OR.NWT.GT.IWT) THEN
      WRITE (6,'(1X,A60)') 'Routine SCR0 fails . Insufficient storage al
     .location.'
      WRITE (6,'(1X,A36,I5)') 'Try rerunning with parameter  IL= ',NWK
      WRITE (6,'(1X,A36,I5)') '               and parameter IWT= ',NWT
      STOP
      ENDIF
C
      DO 13 J=1,NZ
      WK(J)=CF/WK(J)
      B(J)=0.0
   13 WK(J+NZ)=Z(J)
      IWT=NWT
      IL=NWK
      RETURN
      END


      SUBROUTINE SCR1(NZ,EF,ZMSR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/DET/S(ND,NS),FIN(NS),Z(NS),SG(ND,NG),TR(NG)
      COMMON/SCOSI/B(NS),WK(5*NS),LST(NS*NS),WT(NS*NS)
C
C     Update algorithm for Bratseth's version of successive correction
C     using weight vectors generated by SCR0
C
      NWT=0
      NWK=0
      DO 10 JS=1,NZ
      B(JS)=B(JS)+WK(NZ+JS)*WK(JS)
      WK(NZ+JS)=Z(JS)-(1.0+EF)*B(JS)
      JLOC=NWK+1
      NJ=LST(JLOC)
      IF (NJ.EQ.0) GOTO 10
      DO 20 I=1,NJ
      NWK=JLOC+I
      IS=LST(NWK)
      NWT=NWT+1
      WK(NZ+JS)=WK(NZ+JS)-WT(NWT)*B(IS)
   20 WK(NZ+IS)=WK(NZ+IS)-WT(NWT)*B(JS)
   10 NWK=JLOC+NJ
C
      ZMSR=0.0
      DO 30 J=1,NZ
      RES=WK(J+NZ)+EF*B(J)
   30 ZMSR=ZMSR+RES**2
      ZMSR=DSQRT(ZMSR/REAL(NZ))
      RETURN
      END


C ======================================================================
C                        SUBRUTINE OSI(NSTA,NGRI)
C
C   Work structure:
C   - detrends the data
C   - computes interpolation values
C   - in case of numerical problems when inverting the matrix: STOP
C   - returns to ANÁLISIS.F
C
C   External files used: none
C
C   Built-in routines:  EIGEN,LUDCMP(C),LUBKSB(C) 
C                       (C) Copr. 1986-92 Numerical Recipes Software
C   External routines:  DETREN.F 
C
C   External Functions: CRMDL
C
C ======================================================================
      SUBROUTINE OSI(nsta,ngri,rep,ind)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),FIN(NS),DIN(NS),SG(ND,NG),TR(NG)
      COMMON/SCOSI/B(NS),WK1(5*NS),LST(NS*NS),WK2(NS,NS)
      COMMON/PROBL/c_ss(ns,ns),c_ss_1(ns,ns)
      dimension indx(ns)
C ----------------------------------------------------------------------
	IF((REP.EQ.2).AND.(IND.EQ.1)) GOTO 4
      WRITE(6,3)
    3 FORMAT(7(/),1X,'########## OPTIMAL STATISTICAL INTERPOLATION ',
     .' METHOD  ##########'//)
C
C     PRELIMINARY STEPS: 
   4  DO K=1,NSTA
         SS(1,K)=X(K)
         SS(2,K)=Y(K)
         FIN(K)=V(K)
      ENDDO
      DO K=1,NGRI
         SG(1,K)=XG(K)
         SG(2,K)=YG(K)
      ENDDO
C
C     DETRENDING:
	IF((REP.EQ.2).AND.(IND.EQ.1)) GOTO 8
      WRITE(6,5)
    5 FORMAT(1X,'*** OSI requires a previous detrending of the data.',
     .' This will be'/5X,'made effective by fitting a polynomial to ',
     .'the station data.'/)
    7 WRITE(6,'(a,$)') ' >>> INPUT POLYNOMIAL DEGREE [NP<=5]: '
      READ(5,*,ERR=7) NP
      PAR(1)=NP
    8 CALL DETREN(NSTA,NGRI,NP,RMS)
      WRITE(6,9) RMS
    9 FORMAT(/1X,'*** Detrending finished. Station Residuals:  RMS=',
     .F6.2/)
C
C ======================================================================
C             INFORMATION AND SETTING UP THE ANALYSIS PARAMETERS
C ======================================================================     
      IF((REP.EQ.2).AND.(IND.EQ.1)) GOTO 24
	WRITE(6,13)
   13 FORMAT(/1X,'*** Correlation Functions are taken as gausians:'
     ./14X,'W = EXP(-d**2/2*S**2)'/5X,'where S is the characteristic',
     .' scale and d is distance.'/)
   15 WRITE(6,'(a,$)') ' >>> Input value of [S(km)] : '
      READ(5,*,ERR=15) SCL
      PAR(2)=SCL
      WRITE(6,18)
   18 FORMAT(/1X,'*** Additionally, you must input:')
   20 WRITE(6,'(a,$)') 
     .' >>> The noise to signal variance ratio  [gamma] : '
      READ(5,*,ERR=20) GAMMA
      PAR(3)=GAMMA
      WRITE(6,*) '>>> The cut-off wavelength [XCUT(km)] of the spatial'
   23 WRITE(6,'(a,$)') 
     .'     smoothing (set to 0 for no additional smoothing): '
      READ(5,*,ERR=23) XCUT
      FSCL=XCUT/4.
      PAR(4)=XCUT
C
C ======================================================================
C                   CONSTRUCTING CORRELATION MATRICES
C ======================================================================     
C.....Constructing matrix C_ss:
   24  do k=1,nsta
         do j=k+1,nsta
            c_ss(k,j)=CRMDL(SS(1,K),SS(1,J),ND,SCL)
            c_ss(j,k)=c_ss(k,j)
         enddo
         c_ss(k,k)=1.+gamma
      enddo
c
c.... Computing the invers of C_ss (saved as C_ss_1):
   25 WRITE(6,*) 
      WRITE(6,26)
   26	format('     ...computing the inverse of the correlation 
     .matrix...')
      do k=1,nsta
         do j=1,nsta
	      wk2(k,j)=c_ss(k,j)
         enddo
      enddo
      do k=1,nsta
	   do j=1,nsta
            c_ss_1(k,j)=0.0
         enddo
         c_ss_1(k,k)=1.0
      enddo
      call ludcmp(wk2,nsta,ns,indx,d)
      do j=1,nsta
         call lubksb(wk2,nsta,ns,indx,c_ss_1(1,j))
      enddo
C
C.....Test sobre la inversa:
   27 WRITE(6,*) '    ...checking inverse accuracy...'
      ncont=0
      do i=1,nsta
         do j=1,nsta
            xx=0.
            do k=1,nsta
               xx=xx+c_ss(i,k)*c_ss_1(k,j)
            enddo
            if(j.eq.i) then
               if(abs(1.-xx).gt.1e-4) then
                  write(6,*) i,i,xx
                  ncont=ncont+1
               endif
            else
               if(abs(xx).gt.1e-4) then
                  write(6,*) i,j,xx
                  ncont=ncont+1
               endif
            endif
         enddo
      enddo
C
      if(ncont.gt.nsta) then
         write(6,30) gamma
   30    format(/1x,'!!! numerical problems when inverting C_ss... Try
     .  to:'/5x,'(1) increase the value of gamma (presently gamma=',
     .  f6.4,')'/5X,'(2) try to solve the problems decomposing C_ss in
     .  eigenvalues'/9X,'(expensive and not easy...)'/5x,'(3) quit'/)
         write(6,'(a,$)') ' >>> Select option: '
         read(5,*) nop
         if(nop.eq.1) then
            write(6,'(a,$)') ' >>> New value for gamma: '
            read(5,*) gamma2
            do k=1,nsta
               c_ss(k,k)=c_ss(k,k)-gamma+gamma2
            enddo
            gamma=gamma2
            goto 25
         else if(nop.eq.2) then 
            call eigen(nsta)
            goto 27
         else
            stop
         endif
      endif
C
C ======================================================================
C                             OPERATING
C ======================================================================
C ....Product  c_ss_1 * V  (saved in B)
      WRITE(6,*)'    ...computing [C_ss]**-1 * [V] ...'
      do j=1,nsta
         b(j)=0.
         do i=1,nsta
            b(j)=b(j)+c_ss_1(j,i)*din(i)
         enddo
      enddo
C
C     Gridpoint values:
      WRITE(6,*) '    ...computing Grid point interpolated values...'
      d1=dsqrt(scl**2+fscl**2)
      d2=dsqrt(scl**2+2.*fscl**2)
      c1=(scl/d1)**2
      c2=(scl/d2)**2
      do k=1,ngri
C        Adding the trend to the interpolated values:
         vg(k)=tr(k)
         do j=1,nsta
            weight=2.*c1*CRMDL(SG(1,K),SS(1,J),ND,D1)
     .               -c2*CRMDL(SG(1,K),SS(1,J),ND,D2)
            vg(k)=vg(k)+weight*b(j)
         enddo
      enddo
C     Station values:
      WRITE(6,*) '    ...computing Station point interpolated values...'
      do k=1,nsta
C        Adding the trend to the interpolated values:
         vs(k)=fin(k)
         do j=1,nsta
            weight=2.*c1*CRMDL(SS(1,K),SS(1,J),ND,D1)
     .               -c2*CRMDL(SS(1,K),SS(1,J),ND,D2)
            vs(k)=vs(k)+weight*b(j)
        enddo
      enddo
C
      return
      end


C ======================================================================
C                        SUBRUTINE EIGEN(NSTA)
C
C   Work structure:
C   - decomposes the correlation matrix in eigenvectors
C   - eliminates the contribution of very small eigenvalues/eigenvectors
C   - re-computes the invers of the correlation matrix
C   - returns to OSI.F
C
C   External files used: none
C
C   Built-in routines:  JACOBI(C)
C                       (C) Copr. 1986-92 Numerical Recipes Software
C   External routines:  none
C
C ======================================================================
      SUBROUTINE EIGEN(nsta)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/SCOSI/EVAL(NS),WK(5*NS),LST(NS*NS),EV(NS,NS)
      COMMON/PROBL/a(ns,ns),a_(ns,ns)
      dimension indx(ns)
C ----------------------------------------------------------------------
C     Attempts spectral decomposition:
      do i=1,nsta
         do j=1,nsta
            a_(i,j)=a(i,j)
         enddo
      enddo
      call jacobi(a_,nsta,nd,eval,ev,nrot)
c     xx=0.
c     do i=1,nsta
c        xx=xx+eval(i)
c        type*,i,eval(i),xx
c     enddo
c     do i=1,nsta
c        do j=1,i
c           xx=0.
c           do k=1,nsta
c              xx=xx+ev(k,i)*ev(k,j)
c           enddo
c           type*,i,j,xx
c        enddo
c     enddo
c
c     do i=1,nsta
c        do j=1,i
c           xx=0.
c           do k=1,nsta
c              xx=xx+ev(i,k)*ev(j,k)*eval(k)
c           enddo
c           type*,i,j,a(i,j)-xx
c        enddo
c     enddo
C
C     Compute inverse C**-1 = EV * EVAL-1 * EV_t  neglecting the contribution 
C     of eigenvalues smaller than input negative power 
      WRITE(6,*) 'input negative power: '
      read(5,*) npow
      nel=0
      do i=1,nsta
         WRITE(6,*) i,eval(i)
         if(eval(i).gt.1./10.**npow) then
            do j=1,nsta
               ev(j,i)=ev(j,i)/dsqrt(eval(i))
            enddo
         else
            WRITE(6,*) '...eliminated'
            nel=nel+1
            do j=1,nsta
               ev(j,i)=0.
            enddo
         endif
      enddo
      WRITE(6,*) 'total eliminated: ',nel
c
      do i=1,nsta
         do j=1,nsta
            a_(i,j)=0.
            do k=1,nsta
               a_(i,j)=a_(i,j)+ev(i,k)*ev(j,k)
            enddo
         enddo
      enddo
C
      RETURN 
      END


      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=500)
      DIMENSION a(np,np),d(np),v(np,np),b(NMAX),z(NMAX)
c*    INTEGER n,np,nrot,NMAX
c*    REAL a(np,np),d(np),v(np,np)
c*    PARAMETER (NMAX=500)
c*    INTEGER i,ip,iq,j
c*    REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
 11     continue
        v(ip,ip)=1.
 12   continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
 13   continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
 14       continue
 15     continue
        WRITE(6,*) i,sm
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     .      .and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+dsqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./dsqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
 16           continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
 17           continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
 18           continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
 19           continue
              nrot=nrot+1
            endif
 21       continue
 22     continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
 23     continue
 24   continue
      pause 'too many iterations in jacobi'
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software 


      SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT REAL*8 (A-H,O-Z)
C     INTEGER n,np,indx(n)
c     REAL a(np,np),b(n)
C     INTEGER i,ii,j,ll
C     REAL sum
      DIMENSION indx(n),a(np,np),b(n)
      ii=0
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
   11       continue
         else if (sum.ne.0.) then
            ii=i
         endif
      b(i)=sum
   12 continue
      do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
   13    continue
      b(i)=sum/a(i,i)
   14 continue
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software 


      SUBROUTINE ludcmp(a,n,np,indx,d)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=5000,TINY=1.0e-20)
C     INTEGER n,np,indx(n),NMAX
C     REAL d,a(np,np),TINY
C     INTEGER i,imax,j,k
C     REAL aamax,dum,sum,vv(NMAX)
      DIMENSION a(np,np),indx(n),vv(NMAX)
      d=1.
      do 12 i=1,n
         aamax=0.
         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
   11    continue
         if (aamax.eq.0.) pause 'singular matrix in ludcmp'
         vv(i)=1./aamax
   12 continue
      do 19 j=1,n
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
   13       continue
            a(i,j)=sum
   14    continue
         aamax=0.
         do 16 i=j,n
            sum=a(i,j)
            do 15 k=1,j-1
               sum=sum-a(i,k)*a(k,j)
   15       continue
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
   16    continue
         if (j.ne.imax)then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
   17       continue
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.)a(j,j)=TINY
         if(j.ne.n)then
            dum=1./a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
   18       continue
         endif
   19 continue
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software



C ======================================================================
C                     FUNCTION CRMDL(S1,S2,ND,SCALE)
C
C     Called by:   SC_STA, SC_OSI, OSI
C
C ======================================================================
      FUNCTION CRMDL(S1,S2,ND,SCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S1(ND),S2(ND)
      SMAX=12.
C
C     To calculate the value of a gaussian correlation function for
C     lag (s1-s2) in ND space.
C
      SPACE=0.0
C***  DO 1 I=1,ND
C***1 SPACE=SPACE+(S1(I)-S2(I))**2
C     Distances between pintsare calculated over the maximum circles on  
C     the globe (ONLY FOR THE BIDIMENSIONAL CASE):
      RT=6370.
      CRMDL=0.
      D=DCOS(S2(2))*DCOS(S1(2))*DCOS(S2(1)-S1(1))
      WAR=D+DSIN(S2(2))*DSIN(S1(2))
      IF(WAR.LT.1) THEN
         ARG=0.5*(RT*DACOS(WAR))**2/SCALE**2
         IF(ARG.LE.SMAX) CRMDL=1.0/EXP(ARG)
      ELSE
         CRMDL=1.
      ENDIF
      RETURN
      END



C ======================================================================
C                     SUBROUTINE DETREN(NSTA,NGRI,NP,RMS)
C
C   Called by:   SUCCOR, OSI
C
C   Work structure:
C   - for a linear detrendind (NP<=1) calls routine DTRN2
C   - otherwise calls routine POLINOM
C   - returns to SUCCOR.F or OSI.F
C
C   External files used: none
C
C   Built-in routines: DTRND2,POLINOM
C   External routines: none
C
C ======================================================================
      SUBROUTINE DETREN(NSTA,NGRI,NP,RMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2)
      COMMON/DET/SS(ND,NS),FIN(NS),DIN(NS),SG(ND,NG),TR(NG)
C ----------------------------------------------------------------------
      IF(NP.LE.1) THEN
         CALL DTRND2(NSTA,NGRI,NP)
      ELSE
C        OPTION 1:  fully bivariate polynomial: (NP+1)*(NP+1) terms
C        OPTION 2:  restrict to the first (NP+1)*(NP+2)/2 terms
         NOPFIT=1
         CALL POLINOM(NSTA,NGRI,NP,NP,NOPFIT)
      ENDIF
C
C     Calculate  residuals in stations...
      rms=0.
      do k=1,nsta
         rms=rms+din(k)**2
      enddo
      rms=dsqrt(rms/real(nsta))
C
      RETURN
      END


C ======================================================================
C                     SUBROUTINE DTRND2(NZ,NGRI,NP)
C
C   External files used: none
C
C   Built-in routines: PBASE0
C   External routines: none
C
C ======================================================================
      SUBROUTINE DTRND2(NZ,NGRI,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2,ID=3)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),ZB(NS),Z(NS),SG(ND,NG),TR(NG)
      DIMENSION D(ID),WK(5*NS)
C
C     Linear detrending routine for two-dimensional fields
C     ZB on input contains NZ observations
C     Z  on output contains detrended observations
C     D on exit contains the (3) parameters of the trend model
C     S on input contains NZ location vectors as columns
C     NP is set equal to 1 by the routine to define trend model order
C
      DO 2 K=1,ID
    2 D(K)=0.0
C
C.....Mean the obs and location variables:
      ZM=0.0
      XM=0.0
      YM=0.0
      CX=0.0
      CY=0.0
      DO 4 I=1,NZ
      ZM=ZM+ZB(I)
      XM=XM+SS(1,I)
    4 YM=YM+SS(2,I)
      ZM=ZM/REAL(NZ)
      XM=XM/REAL(NZ)
      YM=YM/REAL(NZ)
C
C.....Mean field treated as constant if NP=0:
      IF (NP.EQ.0) THEN
      D(1)=ZM
      GOTO 10
      ENDIF
C
C.....Calculate second order moment sums for location variables:
      NP=1
      SX2=0.0
      SY2=0.0
      SXY=0.0
      DO 6 I=1,NZ
      SX2=SX2+(SS(1,I)-XM)**2
      SY2=SY2+(SS(2,I)-YM)**2
    6 SXY=SXY+(SS(1,I)-XM)*(SS(2,I)-YM)
      DEN=SX2*SY2-SXY**2
C
C.....Calculate parameters of trend in mean corrected space:
      DO 8 I=1,NZ
      CX=CX+((SS(1,I)-XM)*SY2 - (SS(2,I)-YM)*SXY)*ZB(I)
    8 CY=CY+((SS(2,I)-YM)*SX2 - (SS(1,I)-XM)*SXY)*ZB(I)
      CX=CX/DEN
      CY=CY/DEN
C
C.....Transform to parameter values in original space and detrend the 
C     observations:
      D(1)=ZM-CX*XM-CY*YM
      D(2)=CX
      D(3)=CY
   10 DO 12 I=1,NZ
      TREND=ZM+(SS(1,I)-XM)*CX+(SS(2,I)-YM)*CY
      Z(I)=ZB(I)-TREND
   12 ZB(I)=TREND
C
C.....Gives the values of the trend at grid points:
      DO K=1,NGRI
         CALL PBASE0(SG(1,K),WK,5*NS,NP,ND,NT)
         TR(K)=0.         
         DO L=1,NT
            TR(K)=TR(K)+D(L)*WK(L)
         ENDDO
      ENDDO
C
      RETURN
      END


      SUBROUTINE PBASE0(S,B,IB,NP,ND,NT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NPOL=20)
      DIMENSION S(ND),B(IB),TERMO(NPOL),TERMN(NPOL)
C
C     Routine to calculate a vector of polynomials complete to order
C     NP in ND coordinate variables sampled by input vector S.
C     The polynomial sequence of NT terms is returned as vector B.
C     Arrays TERM0 and TERMN used as workspace.
C
      NT=1
      B(1)=1.0D0
      IF (NP.EQ.0) RETURN
C
      DO 2 K=1,ND
      TERMN(K)=S(K)
      NT=NT+1
    2 B(NT)=TERMN(K)
      LENN=ND
      IF (NP.EQ.1) RETURN
C
      DO 10 N=2,NP
         LTERM=LENN
         DO 4 I=1,LENN
    4    TERMO(I)=TERMN(I)
         IT=0
         DO K=1,ND
            I0=LENN-LTERM+1
            DO I=I0,LENN
               IT=IT+1
               TERMN(IT)=S(K)*TERMO(I)
               NT=NT+1
               B(NT)=TERMN(IT)
            ENDDO
            IF (K.NE.ND) LTERM=LTERM*(ND-K)/(ND-K+N-1)
         ENDDO
         LENN=LENN*(ND+N-1)/N
   10 CONTINUE
      RETURN
      END


C ======================================================================
C               SUBROUTINE POLINOM(NSTA,NGRI,NPX,NPY,NOPFIT)
c                     original source: Mike Pedder
C                   adapted by: Damia Gomis. May 1999
C                       
C   External files used: none
C
C   Built-in routines: PBASE,BPBASE,PNNORM,ASOLVE,CSOLVE
C                      ARANK (only when cross validation is requested) 
C   External routines: none
C
C ======================================================================
c     Polynomial fitting: Routine 1 out of 7
c     --------------------------------------
      SUBROUTINE POLINOM(NSTA,NGRI,NPX,NPY,NOPFIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=2,ID=50)
      COMMON/INOUT/X(NS),Y(NS),V(NS),VS(NS),XG(NG),YG(NG),VG(NG),PAR(5)
      COMMON/DET/SS(ND,NS),ZB(NS),Z(NS),SG(ND,NG),TR(NG)
      DIMENSION BX(ID),BY(ID),BXY(ID),B(NS,ID),P(NS,ID),AC(ID),CC(ID)
      DIMENSION IR(ID),WK(ID)
C     ----------------------------------------------------------------
c     This is a piece of pseudo-code to show how the method is used
c     in a routine that fits a 2-D polynomial and then interpolates
c     onto a grid point.
c
c     The arrays being used are
c     ss(2,ix)       : ss(1,i)=x, ss(2,i)=y location for observation i
c     bx(ib), by(ib) : 1-d work arrays to store polynomial sequences
c     bxy(ib)        : 1-d work array to store 2-d polynomial sequence
c     b(ix,ib)       : the polynomial design matrix
c     p(ix,ib)       : the orthogonal polynomial design matrix
c     zb(ix)         : contains the observations
c     ac(ib)         : coefficients of orthogonal polynomial functions
c     cc(ib)         : coefficients of original polynomial functions
c
c     The declared array sizes must satisfy:-
c     ix>NSTA, where NSTA is the number of station locations
c     id>=npx*npy, where npx,npy are the orders of 1-d polynomials
c                  given as input parameters when calling this routine
c
c     Calculate the elements of the polynomial design matrix from the
c     station location data, for specified values of npx,npy.
c
c     if no differentiation involved, set:
      ndf=0
      nbxy=0
c
      do k=1,nsta
c        basis sequence for x and y :
         call pbase(ss(1,k),bx,id,npx,ndf,nbx) 
         call pbase(ss(2,k),by,id,npy,ndf,nby) 
c        forms 2-d sequence:
         if(nopfit.eq.2) nbxy=-1
         call bpbase(bx,nbx,by,nby,bxy,id,nbxy) 
c        on output, nbxy is the number of terms in the 2-d sequence:
         do l=1,nbxy
c           filling the design matrix:
            b(k,l)=bxy(l)
         enddo
      enddo
c
c     Orthogonalise the polynomial design matrix:
      call pnnorm(nsta,b,ns,p,nbxy) 
c     orthogonal matrix formed in p
c
c     Note that pnorm need only be called once for problems involving
c     analysis of several sets of observations sampled on the same
c     set of locations. It is the most expensive part of the code!
c
c     Now solve for the coefficients of the orthogonal sequence given
c     one set of observations (in array zb) :
      call asolve(nsta,zb,p,ns,ac,nbxy)
c     result formed in array ac
c
c     Optionally call ARANK to study variance explained:
      nopara=1
      if(nopara.eq.1) then
         kfail=0
         srs=0.0
         call arank(nsta,zb,nbxy,ac,NSIG,SRS,IR,WK,KFAIL)
         var=0.0
         do i=1,nsta
            var=var+zb(i)**2
         enddo
         resvar=var
         do j=1,nbxy
            resvar=resvar-ac(ir(j))**2
c           write(6,*) j,ir(j),var,ac(ir(j))**2
c           write(6,*) j,resvar/real(nsta),resvar/real(nsta-j)
         enddo
c        resvar/real(nsta) is the mean-square residual
c        resvar/real(nsta-j) is the estimated residual variance
c        ir(j) identifies the index of the term in the original polinomial
c              sequence. The order is such that this term always explains 
c              more variance than any of the following terms
      endif
c
c     Calculate corresponding sequence of coefficients for original
c     polynomial basis functions:
      call csolve(nsta,b,ns,p,nbxy,ac,cc)
c     result formed in array cc
c
c     Calculating the fitted field at location xg,yg:
      do k=1,ngri
         call pbase(sg(1,k),bx,id,npx,ndf,nby)
         call pbase(sg(2,k),by,id,npy,ndf,nby)
         if(nopfit.eq.2) nbxy=-1
         call bpbase(bx,nbx,by,nby,bxy,id,nbxy)
         tr(k)=0.0
         do l=1,nbxy
            tr(k)=tr(k)+cc(l)*bxy(l)
         enddo
      enddo
c
      do k=1,nsta
         call pbase(ss(1,k),bx,id,npx,ndf,nby)
         call pbase(ss(2,k),by,id,npy,ndf,nby)
         if(nopfit.eq.2) nbxy=-1
         call bpbase(bx,nbx,by,nby,bxy,id,nbxy)
         trend=0.0
         do l=1,nbxy
            trend=trend+cc(l)*bxy(l)
         enddo
         z(k)=zb(k)-trend
         zb(k)=trend
c	   type*,k,zb(k),z(k)
      ENDDO
C 
      RETURN
      END


c     Polynomial fitting: Routine 2 out of 7
c     --------------------------------------
c     Poynomial fitting in multidimensional space using the orthogonal
c     polynomial approach. Compile with -r8 option to ensure accurate
c     orthogonalization
c
      subroutine pbase(s,b,ib,np,ndf,nt)
      parameter (iwk=20) ! for working arrays index and cf
      implicit real*8 (a-h,o-z)
      integer index(0:iwk)
      dimension b(ib),cf(0:iwk)
c
c     To calculate the sequence of polynomial coefficients of order np in
c     1-d space,
c     with or without differentiation with respect to independent variable s
c     On entry:
c     s is the value of the independent variable
c     np is the order of the polynomial
c     ndf is the order of the differentiation (=0 for none)
c     On exit:
c     nt is the number of terms in the series
c     b(nt) contains the ordered sequence of terms in the series
c
      do k=0,np
         index(k)=k
         cf(k)=1.0
         if (ndf.gt.0) then
            do id=1,ndf
               cf(k)=cf(k)*index(k)
               index(k)=k-id
            enddo
         endif
      enddo
c
      nt=0
      do k=0,np
         nt=nt+1
         b(nt)=0.0
         if (index(k).eq.0) then
            b(nt)=cf(k)
         else if (index(k).gt.0) then
            b(nt)=cf(k)*s**index(k)
         endif
      enddo
c
      return
      end


c     Polynomial fitting: Routine 3 out of 7
c     --------------------------------------
      subroutine bpbase(bx,nx,by,ny,bxy,ib,nxy)
      implicit real*8(a-h,o-z)
      dimension bx(nx),by(ny),bxy(ib)
c
c     routine to generate a multidimensional sequence of polynomial
c     basis functions, as the product of two sequences.
c     On ENTRY
c           bx must contain one sequence already generated,
c           e.g. in x-space (but could be multidimensional)
c           by must contain a second sequence already generated,
c           e.g. in y-space (but could be multidimensional)
c           if nxy<0 on entry, then for nx=ny the routine will
c           construct a sequence containing only terms of order <=nx
c     On EXIT
c           bxy contains the resulting sequence
c           nxy is the total number of terms in this sequence
c
      ig=0
      if(nx.eq.ny.and.nxy.lt.0) ig=1   ! to select subset of terms
      nxy=0 ! counter for terms in bxy
      do ix=1,nx
         do iy=1,ny-ig*(ix-1)
            nxy=nxy+1
            bxy(nxy)=bx(ix)*by(iy)
         enddo
      enddo
      return
      end


c     Polynomial fitting: Routine 4 out of 7
c     --------------------------------------
      SUBROUTINE pnnorm(NZ,B,IB,P,NT)
      implicit real*8 (a-h,o-z)
      dimension B(IB,NT),P(IB,NT)
C     integer nz,ib,nt
C
C     FORM NT ORTHONORMAL VECTORS IN P FROM NT BASIS VECTORS IN B
C
      DOTB=0.0D0
      DO 1 I=1,NZ
      P(I,1)=B(I,1)
    1 DOTB=DOTB+P(I,1)*P(I,1)
      DOTB=1.0D0/SQRT(DOTB)
      DO 2 I=1,NZ
    2 P(I,1)=P(I,1)*DOTB
      DO 9 L=2,NT
      DOTB=0.0D0
      DO 3 I=1,NZ
      P(I,L)=B(I,L)
    3 DOTB=DOTB+B(I,L)*B(I,L)
      DOTB=1.0D0/SQRT(DOTB)
      DO 4 I=1,NZ
    4 P(I,L)=DOTB*P(I,L)
C     ITERATE TO ORTHOGONALITY ON PREVIOUS VECTORS.
      L1=L-1
      DO 5 J=1,L1
      DOTP=0.0D0
      DO 6 I=1,NZ
    6 DOTP=DOTP+P(I,J)*P(I,L)
      DOTB=0.0D0
      DO 7 I=1,NZ
      P(I,L)=P(I,L)-DOTP*P(I,J)
    7 DOTB=DOTB+P(I,L)*P(I,L)
      DOTB=1.0D0/SQRT(DOTB)
      DO 8 I=1,NZ
    8 P(I,L)=DOTB*P(I,L)
    5 CONTINUE
    9 CONTINUE
      RETURN
      END


c     Polynomial fitting: Routine 5 out of 7
c     --------------------------------------
      SUBROUTINE asolve(NZ,Z,P,IP,A,NT)
      implicit real*8 (a-h,o-z)
      dimension P(IP,NT),Z(NZ),A(NT)
c     integer nz,ip,nt
C
C     Estimate the parameters of the orthogonal polynomial
C     approximation to the field observed by data vector Z
C
      DOTA=0.0D0
      DO 1 J=1,NT
      A(J)=0.0D0
      DO 2 I=1,NZ
    2 A(J)=A(J)+P(I,J)*Z(I)
    1 CONTINUE
      RETURN
      END


c     Polynomial fitting: Routine 6 out of 7
c     --------------------------------------
      SUBROUTINE csolve(NZ,B,IB,P,NT,A,C)
      implicit real*8 (a-h,o-z)
      dimension B(IB,NT),P(IB,NT),A(NT),C(NT)
c     integer nz,ib,nt
C
C     RECOVER C PARAMETERS OF BASIS VECTORS B FROM THE A PARAMETERS
C     OF ORTHONORMAL VECTORS P
C
      DOTA=0.0D0
      DO 1 I=1,NZ
    1 DOTA=DOTA+P(I,NT)*B(I,NT)
      C(NT)=A(NT)/DOTA
      IF (NT.EQ.1) RETURN
      N1=NT-1
      DO 2 J=1,N1
      K=NT-J
      DOTB=0.0D0
      DO 3 I=1,NZ
    3 DOTB=DOTB+P(I,K)*B(I,K)
      SUM=0.0D0
      K1=K+1
      DO 4 L=K1,NT
      DOTC=0.0D0
      DO 5 I=1,NZ
    5 DOTC=DOTC+P(I,K)*B(I,L)
    4 SUM=SUM+DOTC*C(L)
    2 C(K)=(A(K)-SUM)/DOTB
      RETURN
      END


c     Polynomial fitting: Routine 7 out of 7
c     --------------------------------------
      SUBROUTINE arank(NZ,Z,NT,A,NSIG,SRS,IR,WK,KFAIL)
      implicit real*8 (a-h,o-z)
      dimension Z(NZ),WK(NT),A(NT)
      dimension ir(nt)
c     integer nz,nt,nsig,ir(nt),kfail
      LOGICAL FLAG
C
C     Optional routine to examine variance explained by the
C     orthogonal components of a polynomial approximation to a field
C     observed by data vector Z.
C     IF KFAIL>0 ON ENTRY, THEN CALCULATE RMS RESIDUAL AND RETURN.
C     OTHERWISE:-
C     RANKS A PARAMETERS OF ORTHONORMAL BASIS ACCORDING TO VARIANCE
C     EXPLAINED. REJECTS REDUNDANT COMPONENTS USING CRITERION
C                RESID VARIANCE=SRS**2 (SRS>0)
C             OR CROSS VALIDATION VARIABLE MINIMIZED (SRS<0)
C
C
      DO 10 J=1,NT
      IR(J)=J
   10 WK(J)=A(J)**2
      VAR=0.0D0
      DO 1 I=1,NZ
    1 VAR=VAR+Z(I)**2
      NSIG=NT
      IF (KFAIL.GT.0) GOTO 20
C
C     ORDER PARAMETERS IN DECREASING ORDER OF A**2
C
      DO 2 IRANK=1,NT
      XMAX=0.0D0
      DO 3 J=1,NT
      IF (WK(J).LT.XMAX) GOTO 3
      XMAX=WK(J)
      K=J
    3 CONTINUE
      IR(IRANK)=K
    2 WK(K)=0.0D0
C
C     TEST OPTIMALITY CRITERION
C
      RSS=VAR
      TEST0=VAR/REAL(NZ)
      FLAG=.FALSE.
      DO 4 K=1,NT
      IF (FLAG) GOTO 4
      NSIG=K
      IF (K.EQ.NZ) THEN
      KFAIL=2
      GOTO 4
      ENDIF
      RSS=RSS-A(IR(K))**2
      TEST1=RSS*REAL(NZ)/REAL(NZ-K)**2
      IF (SRS.GT.0.0D0.AND.RSS/REAL(NZ).LT.SRS**2) THEN
      FLAG=.TRUE.
      GOTO 4
      ELSE IF (SRS.LT.0.0D0.AND.TEST1.GT.TEST0) THEN
      NSIG=K-1
      FLAG=.TRUE.
      GOTO 4
      ENDIF
      TEST0=TEST1
    4 CONTINUE
C
C     SET REDUNDANT COMPONENTS TO ZERO AND GET RESIDUAL VARIANCE
C
   20 SRS=0.0D0
      IF (NSIG.EQ.NZ) RETURN
      IF (NSIG.EQ.NT) KFAIL=1
      DO 9 K=1,NT
      IF (K.GT.NSIG) THEN
      A(IR(K))=0.0D0
      GOTO 9
      ENDIF
      WK(K)=A(IR(K))**2
      VAR=VAR-WK(K)
    9 CONTINUE
      SRS=SQRT(VAR/REAL(NZ-NSIG))
      RETURN
      END
C ===========================================================================


