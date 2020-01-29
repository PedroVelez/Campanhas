C ======================================================================
C                 INVERSION OF THE OMEGA EQUATION
C
C   COMPUTES THE 3D VERTICAL VELOCITY FIELD SOLVING THE OMEGA EQUATION 
C   BY MEANS OF A RELAXATION METHOD. COORDINATES (X,Y,Z) ARE ASSUMED.
C
C                            Damia Gomis
C
C      IMEDEA (Institut Mediterrani d'Estudis Avançats), a joint 
C      centre between the Universitat de les Illes Balears (UIB) and 
C      the Spanish Research Council (CSIC).
C
C      This code has been compiled and distributed in the framework 
C      of the project REN2000-2599-E, funded by the Spanish Marine 
C      Science and Technology subprogram.
C
C ======================================================================
C
C                   Last version: 1 december 2001
C
C ======================================================================
C
C   Work structure:
C   - reads information on the grid/domain (from an external file)
C   - selects a few options and reads some other parameters from an 
C     external file.
C   - reads input files with the grid values of the vertical forcing
C     and the mean level density
C   - starts iterating to solve the equation. Stops when the scheme 
C     converges (i.e., when the effect of subsequent iterations is 
C     less than the prescribed accuracy)
C   - writes the grid values of the vertical velocity field 
C     variable and writes results
C
C ======================================================================
C
C   External files used:  
C
C   - INPUT files:
C     * domain data file [../info/grid.dat]
C     * files with the grid point values of the vertical forcing at 
C       the set of levels conforming the 3D grid [NNNNLLLLqdi.snp or 
C       NNNNLLLqdi.grd]
C     * file with the mean level density [NNNN_st0.dat]
C     * file [ominput.dat] containing the set of parameters needed 
C       to invert the equation. Moreover:
C       -> if NBOUND=2: file [ombc.dat] containing the distributions
C                       of w at the (6) 2D boundaries
C       -> if NBOUND=3: file containing the bathimetry [batim.dat]
C                       at every gridpoint expressed in meters
C
C   - OUPUT files:
C     * files with the grid point values of the vertical velocity at 
C       the set of levels conforming the 3D grid [NNNNLLLLdhw.snp or 
C       NNNNLLLLdhw.grd]
C
C   Built-in routines:  none
C   External routines:  none
C
C ======================================================================
C
C      UNITS OF THE MAIN PARAMETERS USED THROUGH THIS CODE:
C
C      all are expressed in SI units
C
C      UNITS OF THE OUTPUT FIELD:  1E-5 m/day *
C
C      * These units should be suitable for the output formats in most 
C        cases. Nevertheless, if for any reason it happens that output
C        values are too large, they will be scaled to fit within the 
C        output formats. Output units are in any case reported at the 
C        end of each file.
C
C ======================================================================
C ======================================================================
      REAL FORZ(0:100,0:100,0:100),N22(0:100,0:100,0:100),N2(0:100)
      REAL UU(0:100,0:100,0:100),VV(0:100,0:100,0:100)
      REAL W(0:100,0:100,0:100),ZW(0:100,0:100,0:100)
      REAL ROM(0:100),ROMU(0:100),COEF(0:100),C1(0:100),C2(0:100)
      REAL VMIN(0:100),VMAX(0:100)
      REAL N2PARAM
      INTEGER NBAT(0:100,0:100),NBAT_(0:100,0:100)
      CHARACTER*4,NA1
      CHARACTER*4,NA0,ME
      CHARACTER*8,NA(0:100)
      CHARACTER*20,LIN(5)
C
      G=9.81
      PI=3.14159
      RD=PI/180.
      XKGLA=60.*1.852
C
C ======================================================================
C     READING INFO FROM EXTERNAL FILES AND SELECTING OPTIONS...
C
C.....Reading grid information from external file:
      open(unit=99,file='../info/grid.dat',status='old')
      read(99,3)
      read(99,3)
      read(99,*) xlon1,xlat1
      read(99,3)
      read(99,*) alf0
      read(99,3) 
      read(99,*) nl,nf,nc
      read(99,3)
      read(99,*) xarm,yarm
      read(99,3)
      read(99,*) z0,dz
    3 format(1x/1x)
      close(unit=99)
C
C ....Computing the Latitud of the central point of the domain (the 
C     Q-vector formulation assumes a constant value for the Coriolis 
C     parameter) and expressing the grid cell size in meters: 
      ALF0=ALF0*RD
      IF(ALF0.LT.0.01) THEN
         Y00=XLAT1+YARM*REAL(NF-1)/2.
         DY=YARM*XKGLA*1000.
         DX=XARM*XKGLA*1000.*COS(RD*Y00)
      ELSE
         XL=XARM*REAL(NC-1)/2.
         YL=YARM*REAL(NF-1)/2.
         DD=SQRT(XL**2+YL**2)
         AL=ATAN2(YL,XL)
         Y00=XLAT1+DD*SIN(AL+ALF0)/XKGLA
         DY=YARM*1000.
         DX=XARM*1000.
      ENDIF
      W2=1.45444/1E4
      F0=W2*SIN(Y00*RD)
      F02=F0*F0
      DX2=DX*DX
      DY2=DY*DY
      DZ2=DZ*DZ

C     Selecting options:
      WRITE(6,5)
    5 FORMAT(////10X,'##### INVERSION OF THE OMEGA EQUATION #####'//
     .1X,'This code can only operate on a 3D domain, so that vertical'
     ./1X,'forcing grid values must be available on the whole set of'/
     .1X,'levels quoted in the input file GRID.DAT '/)
    8 WRITE(6,'(a,$)') ' >>> ROOT OF THE GRID FILE(S) (NNNN): '
      READ(5,10,ERR=8) NA0
   10 FORMAT(A4)
      ZMAX=Z0+DZ*FLOAT(NL-1)
      DO I=1,NL+1
         LEV=IFIX(Z0+DZ*FLOAT(I-1)+0.00001)
         IF(ZMAX.GT.999) LEV=LEV/10
C         NA1=CHAR(LEV/100+48)//CHAR((LEV-100*(LEV/100))/10+48)//
C     .       CHAR(LEV-10*(LEV/10)+48)

		NA1=CHAR(lev/1000+48)//
     .	CHAR(lev/100-10*(lev/1000)+48)//     
     .	CHAR(lev/10-10*(lev/100)+48)//
     .	CHAR(lev-10*(lev/10)+48)
         NA(I)=NA0//NA1
 
      ENDDO
      NA(0)=NA0//'0000'
   13 WRITE(6,'(a,$)') ' >>> FILE FORMAT:  1: *.GRD   2: *.SNP  : '
      READ(5,*,ERR=13) NOPFOR
      IF(NOPFOR.EQ.1) THEN
         ME='.grd'
      ELSE IF(NOPFOR.EQ.2) THEN
         ME='.snp'
      ELSE
         GOTO 13
      ENDIF
C
C.....Reading parameters from external file:
      OPEN(UNIT=88,FILE='../info/ominput.dat',STATUS='OLD')
      READ(88,*) EPSILON
      READ(88,*) NEX,NEY,NEZ
C      READ(88,*) N2PARAM
	N2PARAM=1
C      READ(88,*) F2PARAM
	F2PARAM=1
      READ(88,*) NBOUND

C
      IF(NEX.LT.2) NEX=2
      IF(NEY.LT.2) NEY=2
      IF((NBOUND.EQ.1).OR.(NBOUND.EQ.3)) THEN
         READ(88,*) WDOW,WUP,WNOR,WSOU,WWES,WEAS
         IF((NBOUND.EQ.3).AND.(WDOW.GT.-98)) THEN
            WRITE(6,15)
   15       FORMAT(/1X,'!!! The bathimetry option selected in the ',
     .      'external file'/1X,'!!! [ominput.dat] implies to assume ',
     .      'zero velocity at the'/1X,'!!! see floor and requires to ',
     .      'have the bathimetry (file'/1X,'!!! [batim.dat]) at every',
     .      ' gridpoint expressed in meters.'/)
   17       WRITE(6,'(a,$)') ' >>> Keep on option 3 or quit ? ([3]/0): '
            READ(5,*,ERR=17) NOP
            IF(NOP.NE.0) THEN
               WDOW=-99
            ELSE
               CLOSE(UNIT=88)        
               STOP
            ENDIF
         ENDIF
      ELSE IF((NBOUND.EQ.0).OR.(NBOUND.EQ.2)) THEN
         WDOW=0.
         WUP=0.
         WNOR=0.
         WSOU=0.
         WWES=0. 
         WEAS=0.
         IF(NBOUND.EQ.2) THEN
            WRITE(6,20)
   20       FORMAT(/1X,'!!! For this option to be applied, there must
     .       be an external'/1X,'!!! file [ombc.snp] containing (6) 
     .       2D distributions of omega'/1X,'!!! values (with units 
     .       quoted at the end.'/)
   23       WRITE(6,'(a,$)') ' >>> Keep on option 2 or quit ? ([2]/0): '
            READ(5,*,ERR=23) NOP
            IF(NOP.EQ.0) THEN
               CLOSE(UNIT=88)        
               STOP
            ENDIF
         ENDIF
      ENDIF
      CLOSE(UNIT=88)
C
C ======================================================================
C     READING DATA FILES: MEAN LEVEL DENSITY AND FORCING GRID VALUES
C
C     Indexes I,J,K run as follows:
C       - I: along parallels, from WEST to EAST 
C       - J: along meridians, from SOUTH to NORTH
C       - K: along the vertical, from BOTTOM to TOP
C
	WRITE(6,*) 'reading mean density file: ../info/'//NA0//'_st0.dat'
      OPEN(UNIT=11,FILE='../info/'//NA0//'_st0.dat',STATUS='OLD')
	DO K=NL,1,-1
         READ(11,*) XX,ROM(K)
      ENDDO
      CLOSE(UNIT=11)
C
      DO K=1,NL
	   WRITE(6,33) NA(NL-K+1)//'qdi'//ME
   33	   FORMAT('reading forcing file: ../varderiv/',A15)
         OPEN(UNIT=12,FILE='../varderiv/'//NA(NL-K+1)//'qdi'//ME,
     .            	                              STATUS='OLD')
         IF(NOPFOR.EQ.1) THEN
            READ(12,32) LIN
   32       FORMAT(A20)
            DO J=1,NF
               READ(12,34) (FORZ(I,J,K),I=1,NC)
c              READ(12,35) 
   34          FORMAT(10F10.4)
   35          FORMAT(1X)
            ENDDO
         ELSE
            DO J=1,NF
               READ(12,36) (FORZ(I,J,K),I=1,NC)
   36          FORMAT(8F10.4)
            ENDDO
c           READ(12,35) 
         ENDIF
         READ(12,38) NESC
   38    FORMAT(//13X,I3)
         CLOSE(UNIT=12)
C        Converting Forcing to SI units:
         DO J=1,NF
            DO I=1,NC
               FORZ(I,J,K)=FORZ(I,J,K)*10.**NESC
            ENDDO
         ENDDO
      ENDDO
C
C ....The external NEX columns, NEY rows, NEZ levels are eliminated.
C     W and the bathymetry are set to zero everywhere by default
      NC=NC-2*NEX
      NF=NF-2*NEY
      NL=NL-2*NEZ
      DO I=0,NC+1
         DO J=0,NF+1
            NBAT_(I,J)=1
            DO K=0,NL+1
               W(I,J,K)=0.
               FORZ(I,J,K)=FORZ(I+NEX,J+NEY,K+NEZ)
            ENDDO
         ENDDO
      ENDDO
      DO K=1,NL
         ROM(K)=ROM(K+NEZ)
      ENDDO
C
C ======================================================================
C     BOUNDARY CONDITIONS: 
C     Although option 3 (case of complex bathimetry) is not allowed, 
C     it has been kept in the code for eventual future applications...
C
C     Boundary values are assigned to (present) index values:
C     I=0,I=NC+1 ; J=0,J=NF+1 ; K=0,K=NL+1 .
C     
      IF((NBOUND.EQ.1).OR.(NBOUND.EQ.3)) THEN
C        W IS SET TO A CONSTANT VALUES AT EACH BOUNDARY
         IF(WUP.GT.-98.) THEN
            DO J=0,NF+1
               DO I=0,NC+1
                  W(I,J,NL+1)=WUP/86400
               ENDDO
            ENDDO
         ENDIF
         IF(WSOU.GT.-98.) THEN
            DO K=0,NL+1
               DO I=0,NC+1
                  W(I,0,K)=WSOU/86400
               ENDDO
           ENDDO
         ENDIF
         IF(WNOR.GT.-98.) THEN
            DO K=0,NL+1
               DO I=0,NC+1
                  W(I,NF+1,K)=WNOR/86400
               ENDDO
            ENDDO
         ENDIF
         IF(WEAS.GT.-98.) THEN
            DO K=0,NL+1
               DO J=0,NF+1
                  W(NC+1,J,K)=WEAS/86400
               ENDDO
            ENDDO
         ENDIF
         IF(WWES.GT.-98.) THEN
            DO K=0,NL+1
               DO J=0,NF+1
                  W(0,J,K)=WWES/86400
               ENDDO
            ENDDO
         ENDIF
C        Finally, the bottom of the domain:
         IF((NBOUND.EQ.1).AND.(WDOW.GT.-98.)) THEN
            DO J=0,NF+1
               DO I=0,NC+1
                  W(I,J,0)=WDOW/86400
               ENDDO
            ENDDO
         ENDIF
C
         IF(NBOUND.EQ.3) THEN
C           Case of Complex Bathimetry: the original number of rows 
C           and columns is first considered. The conversion towards 
C           the restricted domain will be accomplished later on by 
C           re-assigning values.
            OPEN(UNIT=77,FILE='batim.dat',STATUS='OLD')
            NC=NC+2*NEX
            NF=NF+2*NEY
            NL=NL+2*NEZ
            DO J=1,NF
               READ(77,50) (NBAT(I,J),I=1,NC)
            ENDDO
   50       FORMAT(20I4)
            CLOSE(UNIT=77)
C           Defining the limits of the domain:
            DO J=1,NF
               DO I=1,NC
C                 NBAT_ is the first data level above the sea floor:
                  IF(NBAT(I,J).LT.Z0) THEN
                     NBAT_(I,J)=NL+1
                  ELSE IF(NBAT(I,J).LT.ZMAX) THEN
                     NBAT_(I,J)=IFIX((ZMAX-REAL(NBAT(I,J)))/DZ)+2
                  ENDIF
               ENDDO
            ENDDO
            IF(NDOW.GT.-98) THEN
C              The option of defining the bottom value of w depending 
C              on the current-bathimetry interaction is not reliable. 
C              Thus, this option (despite being coded) is never applied.
               DO K=1,NL
                  OPEN(UNIT=11,FILE=NA(NL-K+1)//'uu'//ME,STATUS='OLD')
                  OPEN(UNIT=12,FILE=NA(NL-K+1)//'vv'//ME,STATUS='OLD')
                  IF(NOPFOR.EQ.1) THEN
                     READ(11,32) LIN
                     READ(12,32) LIN
                     DO J=1,NF
                        READ(11,34) (UU(I,J,K),I=1,NC)
C                        READ(11,35) 
                        READ(12,34) (VV(I,J,K),I=1,NC)
C                        READ(12,35) 
                     ENDDO
                  ELSE
                     DO J=1,NF
                        READ(11,36) (UU(I,J,K),I=1,NC)
                        READ(12,36) (VV(I,J,K),I=1,NC)
                     ENDDO
                  ENDIF
                  CLOSE(UNIT=11)
                  CLOSE(UNIT=12)
               ENDDO
               FX=0.5/DX
               FY=0.5/DY
               DO J=1,NF
                  DO I=1,NC
                     IF((NBAT_(I,J).GT.1).AND.(NBAT_(I,J).LE.NL)) THEN
                        WX=UU(I,J,K)*REAL(NBAT(I+1,J)-NBAT(I-1,J))*FX
                        WY=VV(I,J,K)*REAL(NBAT(I,J+1)-NBAT(I,J-1))*FY
                        DO KK=0,K-1
                           W(I,J,KK)=(WX+WY)/100.
C                          Velocity is assumed to be in cm/s (/100)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
C
C           Conversion towards the restricted domain:
            NC=NC-2*NEX
            NF=NF-2*NEY
            NL=NL-2*NEZ
            DO I=0,NC+1
               DO J=0,NF+1
                  NBAT_(I,J)=NBAT_(I+NEX,J+NEX)-NEZ
                  IF(NBAT_(I,J).LT.1) NBAT_(I,J)=1
                  DO K=0,NL+1
                     W(I,J,K)=W(I+NEX,J+NEY,K+NEZ)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
C
      ELSE IF(NBOUND.EQ.2) THEN
C
C        2D DISTRIBUTIONS OF W VALUES ARE READ FOR EACH BOUNDARY
C        -The original number of rows and columns is first considered. 
C         The conversion towards the restricted domain is accomplished 
C         later on by re-assigning values.
C        -At the intersection between faces and in case of conflict, 
C         the values of UP/DOWN dominate over lateral faces, and for 
C         these, NORTH/SOUTH dominate over WEST/EAST.
C
         OPEN(UNIT=77,FILE='ombc.snp',STATUS='OLD')
         DO J=1,NF+2*NEY
            READ(77,36) (W(I,J,0),I=1,NC+2*NEX)
            READ(77,36) (W(I,J,NL+1),I=1,NC+2*NEX)
         ENDDO
         DO K=NL+2*NEZ,1,-1
            READ(77,36) (W(I,0,K),I=1,NC+2*NEX)
            READ(77,36) (W(I,NF+1,K),I=1,NC+2*NEX)
         ENDDO
         DO K=NL+2*NEZ,1,-1
            READ(77,36) (W(0,J,K),J=1,NF+2*NEY)
            READ(77,36) (W(NC+1,J,K),J=1,NF+2*NEY)
         ENDDO
         READ(77,38) NESC
         CLOSE(UNIT=77)
C
C        Converting Omega values to SI units:
c         OMFAC=10.**NESC
         DO J=0,NF+1
            DO I=0,NC+1
               W(I,J,0)=W(I+NEX-1,J+NEY-1,0)/86400
               W(I,J,NL+1)=W(I+NEX-1,J+NEY-1,NL+1)/86400
            ENDDO
         ENDDO
         DO K=0,NL+1
            DO I=0,NC+1
               W(I,0,K)=W(I+NEX-1,0,K+NEZ-1)/86400
               W(I,NF+1,K)=W(I+NEX-1,NF+1,K+NEZ-1)/86400
            ENDDO
         ENDDO
         DO K=0,NL+1
            DO J=0,NF+1
               W(0,J,K)=W(0,J+NEY-1,K+NEZ-1)/86400
               W(NC+1,J,K)=W(NC+1,J+NEY-1,K+NEZ-1)/86400
            ENDDO
         ENDDO
C
      ENDIF
C
C ======================================================================
C     CALCULATES PRELIMINARY ARRAYS APPEARING IN THE DISCRETE FORM OF 
C     OF THE OMEGA EQUATION:
C     > COEF(K) IS THE COEFFICIENT MULTIPLYING THE VERTICAL VELOCITY
C     > C1(K) AND C2(K) ARE OTHER COEFFICIENTS

C     > N2(K) IS THE BRUNT-VAISALA FREQ.
C
      DO K=2,NL-1
        N2(K)=-G*(ROM(K+1)-ROM(K-1))/ROM(K)/(2.*DZ)
      ENDDO
C
C     At the top and bottom levels, N2 cannot be computed. It is taken 
C     equal to the values of the neighbour levels:
      N2(1)=N2(2)
      N2(NL)=N2(NL-1)
C
C     Mean level density is extrapolated to the upper and lower 
C     boundaries using the N2 values:
      ROM(0)=2.*N2(1)*ROM(1)*DZ/G+ROM(2)
      ROM(NL+1)=ROM(NL-1)-2.*N2(NL)*ROM(NL)*DZ/G
C
C*    ROMU(K) IS THE AVERAGE DENSITY FOR LEVEL=K+1/2
C*    DO K=0,NL
C*       ROMU(K)=.5*(ROM(K+1)+ROM(K))
C*    ENDDO
C
      IF (N2PARAM.EQ.0) THEN
         DO K=1,NL
            N2(K)=0.
         ENDDO
      ENDIF
      IF (F2PARAM.EQ.0.) THEN
         F02=0.
      ENDIF
C
      DO K=1,NL
         COEF(K)=2.*(N2(K)*(DX2+DY2)/DX2/DY2+F02/DZ2)
C*       COEF(K)=2.*N2(K)*(DX2+DY2)/DX2/DY2+F02*ROM(K)*(ROMU(K)+
C*   .   ROMU(K-1))/DZ2/ROMU(K)/ROMU(K-1)
         C1(K)=N2(K)/COEF(K)
         C2(K)=F02/COEF(K)/DZ2
      ENDDO
C
      DO K=1,NL
         DO J=1,NF
            DO I=1,NC
               FORZ(I,J,K)=FORZ(I,J,K)/COEF(K)
            ENDDO
         ENDDO
      ENDDO
C
C ======================================================================
C     CALCULATES VERTICAL VELOCITY W(I,J,K)
C
      WRITE(6,*) ' '
      WRITE(6,*) 'Iterating...'
      M=0
   60 M=M+1
      DO I=0,NC+1
         DO J=0,NF+1
            DO K=NBAT_(I,J)-1,NL+1
               ZW(I,J,K)=W(I,J,K)
            ENDDO
         ENDDO
      ENDDO
C
C     Neumann conditions:
      IF(WDOW.LT.-98) THEN
         DO I=1,NC
            DO J=1,NF
               W(I,J,NBAT_(I,J)-1)=W(I,J,NBAT_(I,J)+1)
            ENDDO
         ENDDO
      ENDIF
      IF(WUP.LT.-98) THEN
         DO I=1,NC
            DO J=1,NF
               W(I,J,NL+1)=W(I,J,NL-1) 
            ENDDO
         ENDDO
      ENDIF
      IF(WSOU.LT.-98) THEN
         DO I=1,NC
            DO K=NBAT_(I,J),NL
               W(I,0,K)=W(I,2,K)
            ENDDO
         ENDDO
      ENDIF
      IF(WNOR.LT.-98) THEN
         DO I=1,NC
            DO K=NBAT_(I,J),NL
               W(I,NF+1,K)=W(I,NF-1,K) 
            ENDDO
         ENDDO
      ENDIF
      IF(WWES.LT.-98) THEN
         DO J=1,NF
            DO K=NBAT_(I,J),NL
               W(0,J,K)=W(2,J,K)
            ENDDO
         ENDDO
      ENDIF
      IF(WEAS.LT.-98) THEN
         DO J=1,NF
            DO K=NBAT_(I,J),NL
               W(NC+1,J,K)=W(NC-1,J,K) 
            ENDDO
         ENDDO
      ENDIF
C
      DO I=1,NC
         DO J=1,NF
            DO K=NBAT_(I,J),NL
               WAR=(W(I+1,J,K)+W(I-1,J,K))/DX2+
     .             (W(I,J+1,K)+W(I,J-1,K))/DY2
               WUD=W(I,J,K+1)+W(I,J,K-1)
C*             WUD=ROM(K+1)*W(I,J,K+1)/ROMU(K)+
C*   .             ROM(K-1)*W(I,J,K-1)/ROMU(K-1)
               W(I,J,K)=C1(K)*WAR+C2(K)*WUD-FORZ(I,J,K)
            ENDDO
         ENDDO
      ENDDO
C
C     Checks Convergence...
      DIFWM=0.
      INDEX=0
      NEVE=0
      DO J=1,NF
         DO I=1,NC
            DO K=NBAT_(I,J),NL
               DIFW=ABS(W(I,J,K)-ZW(I,J,K))
               DIFWM=DIFWM+DIFW
               NEVE=NEVE+1
               IF(DIFW.GT.EPSILON) INDEX=1
            ENDDO
         ENDDO
      ENDDO
      DIFWM=DIFWM/REAL(NEVE)
      WRITE (6,*) M,DIFWM
      IF(INDEX.EQ.1) GOTO 60
C
C ======================================================================
C     CONVERGENCE WAS SUCCESFUL... SAVE RESULTS
C
C     The Vertical velocity is expressed in m/day.
C     Boundary conditions are written at the boundaries of the domain:
      NL=NL+2*NEZ
      NC=NC+2*NEX
      NF=NF+2*NEY
      DO K=0,NL+1
         VMIN(K)=1E6
         VMAX(K)=-1E6
         DO I=1,NC
            DO J=1,NF
            IF(NOPFOR.EQ.1) THEN
               ZW(I,J,K)=9999.
            ELSE
               ZW(I,J,K)=0.
            ENDIF
            ENDDO
         ENDDO
      ENDDO
      DO K=NEZ,NL-NEZ+1
         DO I=NEX,NC-NEX+1
            DO J=NEY,NF-NEY+1
               ZW(I,J,K)=W(I-NEX,J-NEY,K-NEZ)*86400.
               IF(ZW(I,J,K).LT.VMIN(K)) VMIN(K)=ZW(I,J,K)
               IF(ZW(I,J,K).GT.VMAX(K)) VMAX(K)=ZW(I,J,K)
            ENDDO
         ENDDO
      ENDDO
C
      WRITE(6,70)
   70 FORMAT(5(/),8X,'============================================'//
     .8X,'NAME OF THE OUTPUT FILES:'/)
      DO K=NL+1,0,-1
         WRITE(6,72) NA(NL-K+1)//'dhw'//ME
   72    FORMAT(23X,A15)
C        Checking max/min values:
         NN=0
         IF(VMAX(K)/10**NN.GT.9999.) THEN
            DO N=1,10
               IF(VMAX(K)/10**(NN+N).LT.9999.) GOTO 75
            ENDDO
   75       NN=NN+N
         ENDIF
         IF(VMIN(K)/10**NN.LT.-999.) THEN
            DO N=1,10
               IF(VMIN(K)/10**(NN+N).GT.-999.) GOTO 78
            ENDDO
   78       NN=NN+N
         ENDIF
C        Re-scaling values:
         IF(NN.GT.0) THEN
            FAC=10.**NN
            VMIN(K)=VMIN(K)/FAC
            VMAX(K)=VMAX(K)/FAC
            DO I=NEX,NC-NEX+1
               DO J=NEY,NF-NEY+1
                  ZW(I,J,K)=ZW(I,J,K)/FAC
               ENDDO
            ENDDO
         ENDIF
C        Writing:
         OPEN(UNIT=34,FILE='../varderiv/'//NA(NL-K+1)//'dhw'//ME,
     .	                                  STATUS='UNKNOWN')
         IF(NOPFOR.EQ.1) THEN
            WRITE(34,32) (LIN(N),N=1,4)
            WRITE(34,80) VMIN(K),VMAX(K)
   80       FORMAT(2(1X,F8.4))
            DO J=1,NF
               WRITE(34,34) (ZW(I,J,K),I=1,NC)
C               WRITE(34,35)
            ENDDO
         ELSE
            DO J=1,NF
               WRITE(34,36) (ZW(I,J,K),I=1,NC)
            ENDDO
         ENDIF
         WRITE(34,82) NN
   82    FORMAT(//1X,'Units: 10**(',I3,') m/day')
         CLOSE(UNIT=34)
      ENDDO
      WRITE(6,84)
   84 FORMAT(/8X,'============================================'///)
C
      END

