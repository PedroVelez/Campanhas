C ======================================================================
C                COMPUTATION OF DERIVED VARIABLES FROM 
C     GRIDDED VALUES OF OBSERVED VARIABLES (THE OUTPUTS OF ANÁLISIS.F)
C                        BY FINITTE DIFERENCES
C
C                             Damia Gomis
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
C   - the derived variable to be computed is selected. Options are:
C     * GEOSTROPHIC VELOCITY and GEOSTROPHIC RELATIVE VORTICITY
C       from the DH field.
C     * RELATIVE VORTICITY and DIVERGENCE from the VELOCITY field.
C     * THE VERTICAL FORCING from the DH and DENSITY fields.
C   - reads input files with the grid values of the corresponding 
C     observed variable, computes the grid values of the derived 
C     variable and writes results
C
C     IT CAN OPERATE EITHER ON A SINGLE LEVEL OR FOR THE WHOLE SET 
C     OF LEVELS CONSTITUTING THE 3D GRID.
C
C ======================================================================
C
C   External files used:  
C   - INPUT files:
C     * domain data file [../info/grid.dat]
C     * file(s) with interpolated grid point values of observed 
C       variables [NNNNLLLLVV.snp or NNNNLLLLVV.grd]
C   - OUPUT files:
C     * file(s) with the computed grid point values of the derived 
C       variable(s) [NNNNLLLLVVV.snp or NNNNLLLLVVV.grd]
C     * file with the mean level density [NNNN_st0.dat] to be used
C       for the computation of the vertical velocity (only for the
C       option of the Vertical Forcing).
C
C   Built-in routines:  POLFIT, FINDIF
C   External routines:  none
C
C ======================================================================
C
C      UNITS OF THE MAIN PARAMETERS USED THROUGH THIS CODE:
C
C      grid cell size:                DX,DY     m
C      Coriolis parameter:            F         s**-1 
C
C      UNITS OF THE OUTPUT FIELDS:*
C
C      geostrophic current:           U,V       1     cm/s 
C      divergence / vorticity:        DIV,VOR   1E-5  s**-1    
C      vorticity advection:           XVA       1E-11 s**-2
C      Q-vector:                      QU,QV     1E-13 s**-3
C      Q-vector divergence:           QDI       1E-17 m**-1 s**-3
C
C      * These units should be suitable for the output formats in most 
C        cases. However, in very energetic regions it might happen that
C        output values are too large, in which case they will be scaled
C        to fit within the output formats. Output units are in any case
C        reported at the end of each file.
C
C ======================================================================
      PARAMETER(NLMAX=100,NCMAX=100,NFMAX=100)
      DIMENSION V11(NFMAX,NCMAX),V12(NFMAX,NCMAX),V31(NFMAX,NCMAX),
     .V32(NFMAX,NCMAX),V33(NFMAX,NCMAX),V34(NFMAX,NCMAX),R0(NLMAX)
	DIMENSION V(9)
      CHARACTER*4,NA0,NA1,ME
      CHARACTER*8,NA(NLMAX)
      CHARACTER*20,LIN(5)
      PI=3.14159
      RD=PI/180.
      XKGLA=60.*1.852
C ----------------------------------------------------------------------
C     Reading info and selecting options...
C
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
      read(99,*) p0,pint
    3 format(1x/1x)
      close(unit=99)
C
C ....Computing the Latitud of the central point of the domain (the 
C     Q-vector formulation (option 5) assumes a constant value for 
C     the Coriolis parameter). DX, DY and F0 will be in MKS units.
      ALF0=ALF0*RD
      IF(ALF0.LT.0.01) THEN
         Y00=XLAT1+YARM*REAL(NF-1)/2.
         DY=YARM*XKGLA*1000.
         DX0=XARM*XKGLA*1000.
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
C
C     Selecting options:
      WRITE(6,2)
    2 FORMAT(////10X,'##### COMPUTATION OF DERIVED VARIABLES #####'/)
      WRITE(6,*) '1. DIVERGENCE FIELD, from actual U,V grids.'
      WRITE(6,*) ' '
      WRITE(6,*) '2. RELATIVE VORTICITY, from actual U,V grids. Option '
      WRITE(6,*) '   to compute also the ADVECTION OF VORTICITY by the'
      WRITE(6,*) '   actual current.'
      WRITE(6,*) ' '
      WRITE(6,*) '3. GEOSTROPHIC VELOCITY, from the DH grid.'
      WRITE(6,*) ' '
      WRITE(6,*) '4. GEOSTROPHIC RELATIVE VORTICITY, from the DH grid.'
      WRITE(6,*) 'Option to compute also the ADVECTION OF GEOSTROPHIC'
      WRITE(6,*) 'VORTICITY BY THE GEOSTROPHIC CURRENT.'
      WRITE(6,*) ' '
      WRITE(6,*)'5. VERTICAL FORCING (right hand side of the OMEGA EQ).'
      WRITE(6,*) 'Expressed in terms of the divergence of the Q_VECTOR'
      WRITE(6,*) 'it can be computed from DH and DENSITY grids of the'
      WRITE(6,*) 'same horizontal level.'
      WRITE(6,*) ' '
    5 WRITE(6,'(a,$)') ' >>> SELECT OPTION: '
      READ(5,*,ERR=5) NOPVAR
      IF((NOPVAR.LT.1).OR.(NOPVAR.GT.5)) GOTO 5
      IF((NOPVAR.EQ.2).OR.(NOPVAR.EQ.4)) THEN
    7    WRITE(6,'(a,$)') 
     .   ' >>> Compute also the ADVECTION OF VORTICITY (SI:1/NO:0) ?: '
         READ(5,*,ERR=7) NOPADV
      ELSE IF(NOPVAR.EQ.5) THEN
    8    WRITE(6,'(a,$)') 
     .   ' >>> Compute also the Q-VECTOR COMPONENTS (SI:1/NO:0) ?: '
         READ(5,*,ERR=8) NOPQVE
      ENDIF
      WRITE(6,*) ' '
      WRITE(6,*)'The code can operate on a single level (op. 1) or for'
      WRITE(6,*)'the whole set of levels conforming the 3D grid (op.2)'
   10 WRITE(6,'(a,$)') ' >>> SELECT OPTION: '
      READ(5,*,ERR=10) NOPLEV
      IF(NOPLEV.EQ.1) THEN
   13    WRITE(6,'(a,$)') ' >>> ROOT OF THE GRID FILE(S) (NNNNLLLL): '
         READ(5,15,ERR=13) NA0,NA1
   15    FORMAT(A4,A4)
   16	   FORMAT(A4)
         NL=1
         NA(1)=NA0//NA1
C	   TYPE*,NA(1)
      ELSE
   17    WRITE(6,'(a,$)') ' >>> ROOT OF THE GRID FILE(S) (NNNN): '
         READ(5,16,ERR=17) NA0
         PMAX=P0+PINT*FLOAT(NL-1)
         DO I=1,NL
            LEV=IFIX(P0+PINT*FLOAT(I-1)+0.00001)
            IF(PMAX.GT.999) LEV=LEV/10
            NA1=CHAR(LEV/1000+48)//CHAR(LEV/100-10*(LEV/1000)+48)//
     .	  CHAR(LEV/10-10*(LEV/100)+48)//CHAR(LEV-10*(LEV/10)+48)
            NA(I)=NA0//NA1
C	      TYPE*,NA(I)
         ENDDO
      ENDIF
   20 WRITE(6,'(a,$)') ' >>> FILE FORMAT:  1: *.GRD   2: *.SNP  : '
      READ(5,*,ERR=20) NOPFOR
      IF(NOPFOR.EQ.1) THEN
         ME='.grd'
      ELSE IF(NOPFOR.EQ.2) THEN
         ME='.snp'
      ELSE
         GOTO 20
      ENDIF
	WRITE(6,*) ' '
C
C     Two more options have been eliminated to simplify the use of 
C     this code, although they remain coded:
C     -> option of computing derivatives by polinomial fitting 
C        (NOPDER=1) instead of by finitte differences (NOPDER=2)
         NOPDER=2
C     -> option of using spatial derivative lags equal to ND times 
C        the grid cell size, in order to smooth the ouput field
         ND=1
         DX0=FLOAT(ND)*DX0
         DX=FLOAT(ND)*DX
         DY=FLOAT(ND)*DY
C
C ----------------------------------------------------------------------
C     Initializing...
      VMIN1=1E6
      VMAX1=-1E6
      VMIN2=1E6
      VMAX2=-1E6
      VMIN3=1E6
      VMAX3=-1E6
      VMIN4=1E6
      VMAX4=-1E6
      DO J=1,NF
         DO K=1,NC
            IF(NOPFOR.EQ.1) THEN
               V11(J,K)=999.
               V12(J,K)=999.
               V31(J,K)=999.
               V32(J,K)=999.
               V33(J,K)=999.
               V34(J,K)=999.
            ELSE
               V11(J,K)=999.
               V12(J,K)=999.
               V31(J,K)=999.
               V32(J,K)=999.
               V33(J,K)=999.
               V34(J,K)=999.
            ENDIF
         ENDDO
      ENDDO
C
      GOTO (30,30,60,60,90) NOPVAR
C
C ----------------------------------------------------------------------
C     OPTIONS 1 AND 2: DIVERGENCE AND RELATIVE VORTICITY:
C
   30 DO I=1,NL
C        Reading Input Files
	   WRITE(6,*) 'Reading  ../varobs/'//NA(I)//'uu'//ME
         OPEN(UNIT=11,FILE='../varobs/'//NA(I)//'uu'//ME,STATUS='OLD')
	   WRITE(6,*) 'Reading  ../varobs/'//NA(I)//'vv'//ME
         OPEN(UNIT=12,FILE='../varobs/'//NA(I)//'vv'//ME,STATUS='OLD')
         IF(NOPFOR.EQ.1) THEN
            READ(11,32) LIN
            READ(12,32) LIN
   32       FORMAT(A20)
            DO J=1,NF
               READ(11,34) (V11(J,K),K=1,NC)
C               READ(11,35) 
               READ(12,34) (V12(J,K),K=1,NC)
C               READ(12,35) 
   34          FORMAT(10F10.4)
   35          FORMAT(1X)
            ENDDO
         ELSE
            DO J=1,NF
               READ(11,36) (V11(J,K),K=1,NC)
               READ(12,36) (V12(J,K),K=1,NC)
   36          FORMAT(8F10.4)
            ENDDO
         ENDIF
         CLOSE(UNIT=11)
         CLOSE(UNIT=12)
C
C        Computing the derived variables
         DO J=ND+1,NF-ND
            IF(ALF0.LT.0.01) DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
            DO K=ND+1,NC-ND
               V(1)=V11(J+ND,K-ND)
               V(2)=V11(J+ND,K)
               V(3)=V11(J+ND,K+ND)
               V(4)=V11(J,K-ND)
               V(5)=V11(J,K)
               V(6)=V11(J,K+ND)
               V(7)=V11(J-ND,K-ND)
               V(8)=V11(J-ND,K)
               V(9)=V11(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DU_DX=V(2)
               DU_DY=V(3)
               V(1)=V12(J+ND,K-ND)
               V(2)=V12(J+ND,K)
               V(3)=V12(J+ND,K+ND)
               V(4)=V12(J,K-ND)
               V(5)=V12(J,K)
               V(6)=V12(J,K+ND)
               V(7)=V12(J-ND,K-ND)
               V(8)=V12(J-ND,K)
               V(9)=V12(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DV_DX=V(2)
               DV_DY=V(3)
C              Units: u,v: cm/s ;
C              --> multiply by 1E3 to obtain output in 1E-5 s**-1
               IF(NOPVAR.EQ.1) THEN
C                 Divergence: 
                  V33(J,K)=1E3*(DU_DX+DV_DY)
               ELSE
C                 Relative Vorticity:
                  V33(J,K)=1E3*(DV_DX-DU_DY)
               ENDIF
               IF(V33(J,K).LT.VMIN3) VMIN3=V33(J,K)
               IF(V33(J,K).GT.VMAX3) VMAX3=V33(J,K)
            ENDDO
         ENDDO
C
C        Writing Output Files
C        Checking max/min values:
         NESC=-5
         NN=0
         IF(VMAX3/10**NN.GT.999.) THEN
            DO N=1,10
               IF(VMAX3/10**(NN+N).LT.999.) GOTO 38
            ENDDO
   38       NN=NN+N
         ENDIF
         IF(VMIN3/10**NN.LT.-999.) THEN
            DO N=1,10
               IF(VMIN3/10**(NN+N).GT.-999.) GOTO 40
            ENDDO
   40       NN=NN+N
         ENDIF
C        Re-scaling values:
         IF(NN.GT.0) THEN
            NESC=NESC+NN
            FAC=10.**NN
            VMIN3=VMIN3/FAC
            VMAX3=VMAX3/FAC
            DO J=ND+1,NF-ND
               DO K=ND+1,NC-ND
                  V33(J,K)=V33(J,K)/FAC
               ENDDO
            ENDDO
         ENDIF
C        Writing:
         IF(NOPVAR.EQ.1) THEN
            OPEN(UNIT=33,FILE='../varderiv/'//NA(I)//'div'//
     .                            	  ME,STATUS='UNKNOWN')
         ELSE
            OPEN(UNIT=33,FILE='../varderiv/'//NA(I)//'rrv'//
     .                        		  ME,STATUS='UNKNOWN')
         ENDIF
         IF(NOPFOR.EQ.1) THEN
            WRITE(33,32) (LIN(L),L=1,4)
            WRITE(33,42) VMIN3,VMAX3
   42       FORMAT(2(1X,F9.4))
            DO J=1,NF
               WRITE(33,34) (V33(J,K),K=1,NC)
C               WRITE(33,35)
            ENDDO
         ELSE
            DO J=1,NF
               WRITE(33,36) (V33(J,K),K=1,NC)
            ENDDO
         ENDIF
         WRITE(33,44) NESC
   44    FORMAT(//1X,'Units: 10**(',I3,') s**-1')
         CLOSE(UNIT=33)
C
C        Optional computation of Vorticity Advection:
         IF(NOPADV.EQ.1) THEN
            DO J=2*ND+1,NF-2*ND
               IF(ALF0.LT.0.01) DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
               DO K=2*ND+1,NC-2*ND
                  V(1)=V33(J+ND,K-ND)
                  V(2)=V33(J+ND,K)
                  V(3)=V33(J+ND,K+ND)
                  V(4)=V33(J,K-ND)
                  V(5)=V33(J,K)
                  V(6)=V33(J,K+ND)
                  V(7)=V33(J-ND,K-ND)
                  V(8)=V33(J-ND,K)
                  V(9)=V33(J-ND,K+ND)
                  IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
                  IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
                  DRV_DX=V(2)
                  DRV_DY=V(3)
C                 Units: u,v: cm/s ; RV: 1E-(5-NN) 1/s ;
C                 --> multiply by 1E(4+NN) to obtain output in 1E-11 s**-2
                  V34(J,K)=10.**(4+NN)*(V11(J,K)*DRV_DX+V12(J,K)*DRV_DY)
                  IF(V34(J,K).LT.VMIN4) VMIN4=V34(J,K)
                  IF(V34(J,K).GT.VMAX4) VMAX4=V34(J,K)
               ENDDO
            ENDDO
C           Checking max/min values:
            NESC=-11
            NN=0
            IF(VMAX4/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX4/10**(NN+N).LT.999.) GOTO 46
               ENDDO
   46          NN=NN+N
            ENDIF
            IF(VMIN4/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN4/10**(NN+N).GT.-999.) GOTO 48
               ENDDO
   48          NN=NN+N
            ENDIF
C           Re-scaling values:
            IF(NN.GT.0) THEN
               NESC=NESC+NN
               FAC=10.**NN
               VMIN4=VMIN4/FAC
               VMAX4=VMAX4/FAC
               DO J=2*ND+1,NF-2*ND
                  DO K=2*ND+1,NC-2*ND
                     V34(J,K)=V34(J,K)/FAC
                  ENDDO
               ENDDO
            ENDIF
C           Writing:
            OPEN(UNIT=34,FILE='../varderiv/'//NA(I)//'rva'//
     .                        		  ME,STATUS='UNKNOWN')
            IF(NOPFOR.EQ.1) THEN
               WRITE(34,32) (LIN(L),L=1,4)
               WRITE(34,42) VMIN4,VMAX4
               DO J=1,NF
                  WRITE(34,34) (V34(J,K),K=1,NC)
C                  WRITE(34,35)
               ENDDO
            ELSE
              DO J=1,NF
                  WRITE(34,36) (V34(J,K),K=1,NC)
               ENDDO
            ENDIF
            WRITE(34,50) NESC
   50       FORMAT(//1X,'Units: 10**(',I3,') s**-2')
            CLOSE(UNIT=34)
         ENDIF
C
      ENDDO
      GOTO 150
C ----------------------------------------------------------------------
C     OPTIONS 3 AND 4: GEOSTROPHIC CURRENT AND RELATIVE VORTICITY:
C
   60 DO I=1,NL
C        Reading Input Files
	   WRITE(6,*) 'Reading  ../varobs/'//NA(I)//'dh'//ME
         OPEN(UNIT=11,FILE='../varobs/'//NA(I)//'dh'//ME,STATUS='OLD')
         IF(NOPFOR.EQ.1) THEN
            READ(11,32) LIN
            DO J=1,NF
               READ(11,34) (V11(J,K),K=1,NC)
C	         READ(11,35)
            ENDDO
         ELSE
            DO J=1,NF
               READ(11,36) (V11(J,K),K=1,NC)
	      ENDDO
         ENDIF
         CLOSE(UNIT=11)
C
C        Computing the derived variables
         DO J=ND+1,NF-ND
            IF(ALF0.LT.0.01) THEN
               F=W2*SIN(RD*(XLAT1+YARM*FLOAT(J-1)))
               DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
            ELSE
               F=W2*SIN(RD*(XLAT1+YARM*FLOAT(J-1)*COS(ALF0)/XKGLA))
            ENDIF
            DO K=ND+1,NC-ND
               V(1)=V11(J+ND,K-ND)
               V(2)=V11(J+ND,K)
               V(3)=V11(J+ND,K+ND)
               V(4)=V11(J,K-ND)
               V(5)=V11(J,K)
               V(6)=V11(J,K+ND)
               V(7)=V11(J-ND,K-ND)
               V(8)=V11(J-ND,K)
               V(9)=V11(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DH_DX=V(2)
               DH_DY=V(3)
               DH_DX2=V(5)
               DH_DY2=V(6)
               IF((NOPVAR.EQ.3).OR.(NOPADV.EQ.1)) THEN
C                 Units:  DH: cm.dyn ;  
C                 --> multiply by 10 to obtain output in cm/s 
                  V31(J,K)=-10.*DH_DY/F
                  V32(J,K)=10.*DH_DX/F
                  IF(V31(J,K).LT.VMIN1) VMIN1=V31(J,K)
                  IF(V31(J,K).GT.VMAX1) VMAX1=V31(J,K)
                  IF(V32(J,K).LT.VMIN2) VMIN2=V32(J,K)
                  IF(V32(J,K).GT.VMAX2) VMAX2=V32(J,K)
               ENDIF
               IF(NOPVAR.EQ.4) THEN
C                 Units: DH: cm.dyn ; 
C                 --> multiply by 10.*100. to obtain output in 1E-5 s**-1 
                  V33(J,K)=1E4*(DH_DX2+DH_DY2)/F
                  IF(V33(J,K).LT.VMIN3) VMIN3=V33(J,K)
                  IF(V33(J,K).GT.VMAX3) VMAX3=V33(J,K)
               ENDIF
            ENDDO
         ENDDO
C
C        Writing Output Files
         IF(NOPVAR.EQ.3) THEN
C           Checking max/min values:
            NESC=0
            NN=0
            IF(VMAX1/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX1/10**(NN+N).LT.999.) GOTO 62
               ENDDO
   62          NN=NN+N
            ENDIF
            IF(VMIN1/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN1/10**(NN+N).GT.-999.) GOTO 64
               ENDDO
   64          NN=NN+N
            ENDIF
            IF(VMAX2/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX2/10**(NN+N).LT.999.) GOTO 66
               ENDDO
   66          NN=NN+N
            ENDIF
            IF(VMIN2/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN2/10**(NN+N).GT.-999.) GOTO 68
               ENDDO
   68          NN=NN+N
            ENDIF
C           Re-scaling values:
            IF(NN.GT.0) THEN
               NESC=NESC+NN
               FAC=10.**NN
               VMIN1=VMIN1/FAC
               VMAX1=VMAX1/FAC
               VMIN2=VMIN2/FAC
               VMAX2=VMAX2/FAC
               DO J=ND+1,NF-ND
                  DO K=ND+1,NC-ND
                     V31(J,K)=V31(J,K)/FAC
                     V32(J,K)=V32(J,K)/FAC
                  ENDDO
               ENDDO
            ENDIF
C           Writing:
            OPEN(UNIT=31,FILE='../varderiv/'//NA(I)//'dhu'//
     .                        		  ME,STATUS='UNKNOWN')
            OPEN(UNIT=32,FILE='../varderiv/'//NA(I)//'dhv'//
     .                        		  ME,STATUS='UNKNOWN')
            IF(NOPFOR.EQ.1) THEN
               WRITE(31,32) (LIN(L),L=1,4)
               WRITE(32,32) (LIN(L),L=1,4)
               WRITE(31,42) VMIN1,VMAX1
               WRITE(32,42) VMIN2,VMAX2
               DO J=1,NF
                  WRITE(31,34) (V31(J,K),K=1,NC)
C	          WRITE(31,35)
                  WRITE(32,34) (V32(J,K),K=1,NC)
C                  WRITE(32,35)
               ENDDO
            ELSE
               DO J=1,NF
                  WRITE(31,36) (V31(J,K),K=1,NC)
                  WRITE(32,36) (V32(J,K),K=1,NC)
               ENDDO
            ENDIF
            WRITE(31,70) NESC
            WRITE(32,70) NESC
   70       FORMAT(//1X,'Units: 10**(',I3,') cm/s') 
            CLOSE(UNIT=31)
            CLOSE(UNIT=32)
         ELSE
C           Checking max/min values:
            NESC=-5
            NN=0
            IF(VMAX3/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX3/10**(NN+N).LT.999.) GOTO 72
               ENDDO
   72          NN=NN+N
            ENDIF
            IF(VMIN3/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN3/10**(NN+N).GT.-999.) GOTO 74
               ENDDO
   74          NN=NN+N
            ENDIF
C           Re-scaling values:
            IF(NN.GT.0) THEN
               NESC=NESC+NN
               FAC=10.**NN
               VMIN3=VMIN3/FAC
               VMAX3=VMAX3/FAC
               DO J=ND+1,NF-ND
                  DO K=ND+1,NC-ND
                     V33(J,K)=V33(J,K)/FAC
                  ENDDO
               ENDDO
            ENDIF
C           Writing:
            OPEN(UNIT=33,FILE='../varderiv/'//NA(I)//'grv'//ME,
     .                                		  STATUS='UNKNOWN')
            IF(NOPFOR.EQ.1) THEN
               WRITE(33,32) (LIN(L),L=1,4)
               WRITE(33,42) VMIN3,VMAX3
               DO J=1,NF
                  WRITE(33,34) (V33(J,K),K=1,NC)
C                  WRITE(33,35)
               ENDDO
            ELSE
               DO J=1,NF
                  WRITE(33,36) (V33(J,K),K=1,NC)
               ENDDO
            ENDIF
            WRITE(33,44) NESC
            CLOSE(UNIT=33)
         ENDIF
C
C        Optional computation of Vorticity Advection:
         IF(NOPADV.EQ.1) THEN
            DO J=2*ND+1,NF-2*ND
               IF(ALF0.LT.0.01) DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
               DO K=2*ND+1,NC-2*ND
                  V(1)=V33(J+ND,K-ND)
                  V(2)=V33(J+ND,K)
                  V(3)=V33(J+ND,K+ND)
                  V(4)=V33(J,K-ND)
                  V(5)=V33(J,K)
                  V(6)=V33(J,K+ND)
                  V(7)=V33(J-ND,K-ND)
                  V(8)=V33(J-ND,K)
                  V(9)=V33(J-ND,K+ND)
                  IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
                  IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
                  DRV_DX=V(2)
                  DRV_DY=V(3)
C                 Units: u,v: cm/s ; RV: 1E-(5-NN) 1/s ;
C                 --> multiply by 1E(4+NN) to obtain output in 1E-11 s**-2
                  V34(J,K)=10.**(4+NN)*(V31(J,K)*DRV_DX+V32(J,K)*DRV_DY)
                  IF(V34(J,K).LT.VMIN4) VMIN4=V34(J,K)
                  IF(V34(J,K).GT.VMAX4) VMAX4=V34(J,K)
               ENDDO
            ENDDO
C           Checking max/min values:
            NESC=-11
            NN=0
            IF(VMAX4/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX4/10**(NN+N).LT.999.) GOTO 76
               ENDDO
   76          NN=NN+N
            ENDIF
            IF(VMIN4/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN4/10**(NN+N).GT.-999.) GOTO 78
               ENDDO
   78          NN=NN+N
            ENDIF
C           Re-scaling values:
            IF(NN.GT.0) THEN
               NESC=NESC+NN
               FAC=10.**NN
               VMIN4=VMIN4/FAC
               VMAX4=VMAX4/FAC
               DO J=2*ND+1,NF-2*ND
                  DO K=2*ND+1,NC-2*ND
                     V34(J,K)=V34(J,K)/FAC
                  ENDDO
               ENDDO
            ENDIF
C           Writing:
            OPEN(UNIT=34,FILE='../varderiv/'//NA(I)//'gva'//ME,
     .		                                STATUS='UNKNOWN')
            IF(NOPFOR.EQ.1) THEN
               WRITE(34,32) (LIN(L),L=1,4)
               WRITE(34,42) VMIN4,VMAX4
               DO J=1,NF
                  WRITE(34,34) (V34(J,K),K=1,NC)
C                  WRITE(34,35)
               ENDDO
            ELSE
               DO J=1,NF
                  WRITE(34,36) (V34(J,K),K=1,NC)
               ENDDO
            ENDIF
            WRITE(34,50) NESC
            CLOSE(UNIT=34)
         ENDIF
C
      ENDDO
      GOTO 150
C ----------------------------------------------------------------------
C     OPTION 5: VERTICAL FORCING
C
   90 OPEN(UNIT=40,FILE='../info/'//NA0//'_st0.dat',STATUS='unknown')
      DO I=1,NL
C        Reading Input Files
	   WRITE(6,*) 'Reading  ../varobs/'//NA(I)//'st'//ME
	   OPEN(UNIT=11,FILE='../varobs/'//NA(I)//'st'//ME,STATUS='OLD')
	   WRITE(6,*) 'Reading  ../varobs/'//NA(I)//'dh'//ME
         OPEN(UNIT=12,FILE='../varobs/'//NA(I)//'dh'//ME,STATUS='OLD')
         IF(NOPFOR.EQ.1) THEN
            READ(11,32) LIN
            READ(12,32) LIN
            DO J=1,NF
               READ(11,34) (V11(J,K),K=1,NC)
C	         READ(11,35)
               READ(12,34) (V12(J,K),K=1,NC)
C               READ(12,35)
            ENDDO
         ELSE
            DO J=1,NF
               READ(11,36) (V11(J,K),K=1,NC)
               READ(12,36) (V12(J,K),K=1,NC)
            ENDDO
         ENDIF
         CLOSE(UNIT=11)
         CLOSE(UNIT=12)
C
C        Preliminary Step: computing the mean level density:

	     LEV=IFIX(P0+PINT*FLOAT(I-1)+0.00001)
C            IF(PMAX.GT.999) LEV=LEV/10
C            NA1=CHAR(LEV/100+48)//CHAR((LEV-100*(LEV/100))/10+48)//
C     .          CHAR(LEV-10*(LEV/10)+48)

	       NA1=CHAR(LEV/1000+48)
     .	//CHAR(LEV/100-10*(LEV/1000)+48)//
     .	CHAR(LEV/10-10*(LEV/100)+48)//CHAR(LEV-10*(LEV/10)+48)

         R0(I)=0.
         DO J=1,NF
            DO K=1,NC
               R0(I)=R0(I)+V11(J,K)
            ENDDO
         ENDDO
         R0(I)=1000.+R0(I)/FLOAT(NF*NC)
         WRITE(40,92) NA1,R0(I)
   92    FORMAT(1X,A7,2X,F10.4)
C
C        Step 1: Computing the Q-Vector:
         DO J=ND+1,NF-ND
            IF(ALF0.LT.0.01) DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
            DO K=ND+1,NC-ND
               V(1)=V11(J+ND,K-ND)
               V(2)=V11(J+ND,K)
               V(3)=V11(J+ND,K+ND)
               V(4)=V11(J,K-ND)
               V(5)=V11(J,K)
               V(6)=V11(J,K+ND)
               V(7)=V11(J-ND,K-ND)
               V(8)=V11(J-ND,K)
               V(9)=V11(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DR_DX=V(2)
               DR_DY=V(3)
               V(1)=V12(J+ND,K-ND)
               V(2)=V12(J+ND,K)
               V(3)=V12(J+ND,K+ND)
               V(4)=V12(J,K-ND)
               V(5)=V12(J,K)
               V(6)=V12(J,K+ND)
               V(7)=V12(J-ND,K-ND)
               V(8)=V12(J-ND,K)
               V(9)=V12(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               D2H_DXY=V(4)
               D2H_DX2=V(5)
               D2H_DY2=V(6)
C              Units:  DH: cm.dyn. ;
C              --> multiply by 10.*1E11. to obtain output in 1E-13 s**-3
               FACT=1E12*9.81/(F0*R0(I))
               V31(J,K)=FACT*(-D2H_DXY*DR_DX+D2H_DX2*DR_DY)
               V32(J,K)=FACT*(-D2H_DY2*DR_DX+D2H_DXY*DR_DY)
               IF(V31(J,K).LT.VMIN1) VMIN1=V31(J,K)
               IF(V31(J,K).GT.VMAX1) VMAX1=V31(J,K)
               IF(V32(J,K).LT.VMIN2) VMIN2=V32(J,K)
               IF(V32(J,K).GT.VMAX2) VMAX2=V32(J,K)
            ENDDO
         ENDDO
C
C        Writing Output Files
         NN=0
         IF(NOPQVE.EQ.1) THEN
C           Checking max/min values:
            NESC=-13
            IF(VMAX1/10**NN.GT.999.) THEN
               DO N=1,10
                  IF(VMAX1/10**(NN+N).LT.999.) GOTO 94
               ENDDO
   94          NN=NN+N
            ENDIF
            IF(VMIN1/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN1/10**(NN+N).GT.-999.) GOTO 96
               ENDDO
   96          NN=NN+N
            ENDIF
            IF(VMAX2/10**NN.GT.9999.) THEN
               DO N=1,10
                  IF(VMAX2/10**(NN+N).LT.9999.) GOTO 98
               ENDDO
   98          NN=NN+N
            ENDIF
            IF(VMIN2/10**NN.LT.-999.) THEN
               DO N=1,10
                  IF(VMIN2/10**(NN+N).GT.-999.) GOTO 100
               ENDDO
  100          NN=NN+N
            ENDIF
C           Re-scaling values:
            IF(NN.GT.0) THEN
               NESC=NESC+NN
               FAC=10.**NN
               VMIN1=VMIN1/FAC
               VMAX1=VMAX1/FAC
               VMIN2=VMIN2/FAC
               VMAX2=VMAX2/FAC
               DO J=ND+1,NF-ND
                  DO K=ND+1,NC-ND
                     V31(J,K)=V31(J,K)/FAC
                     V32(J,K)=V32(J,K)/FAC
                  ENDDO
               ENDDO
            ENDIF
C           Writing:
            OPEN(UNIT=31,FILE='../varderiv/'//NA(I)//'qvu'//ME,
     .		                                STATUS='UNKNOWN')
            OPEN(UNIT=32,FILE='../varderiv/'//NA(I)//'qvv'//ME,
     .        		                        STATUS='UNKNOWN')
            IF(NOPFOR.EQ.1) THEN
               WRITE(31,32) (LIN(L),L=1,4)
               WRITE(32,32) (LIN(L),L=1,4)
               WRITE(31,42) VMIN1,VMAX1
               WRITE(32,42) VMIN2,VMAX2
               DO J=1,NF
                  WRITE(31,34) (V31(J,K),K=1,NC)
C                  WRITE(31,35)
                  WRITE(32,34) (V32(J,K),K=1,NC)
C                  WRITE(32,35)
               ENDDO
            ELSE
               DO J=1,NF
                  WRITE(31,36) (V31(J,K),K=1,NC)
                  WRITE(32,36) (V32(J,K),K=1,NC)
               ENDDO
            ENDIF
            WRITE(31,104) NESC
            WRITE(32,104) NESC
  104       FORMAT(//1X,'Units: 10**(',I3,') s**-3') 
            CLOSE(UNIT=31)
            CLOSE(UNIT=32)
         ENDIF
C
C        Step 2: computing the Q-Vector Divergence:
         DO J=2*ND+1,NF-2*ND
            IF(ALF0.LT.0.01) DX=DX0*COS(RD*(XLAT1+YARM*FLOAT(J-1)))
            DO K=2*ND+1,NC-2*ND
               V(1)=V31(J+ND,K-ND)
               V(2)=V31(J+ND,K)
               V(3)=V31(J+ND,K+ND)
               V(4)=V31(J,K-ND)
               V(5)=V31(J,K)
               V(6)=V31(J,K+ND)
               V(7)=V31(J-ND,K-ND)
               V(8)=V31(J-ND,K)
               V(9)=V31(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DQU_DX=V(2)
               V(1)=V32(J+ND,K-ND)
               V(2)=V32(J+ND,K)
               V(3)=V32(J+ND,K+ND)
               V(4)=V32(J,K-ND)
               V(5)=V32(J,K)
               V(6)=V32(J,K+ND)
               V(7)=V32(J-ND,K-ND)
               V(8)=V32(J-ND,K)
               V(9)=V32(J-ND,K+ND)
               IF(NOPDER.EQ.1) CALL POLFIT(DX,DY,V)
               IF(NOPDER.EQ.2) CALL FINDIF(DX,DY,V)
               DQV_DY=V(3)
C
C              Q-Div Units: QU,QV units: 1E-(13-NN) ; 
C              --> multiply by 1E(4+NN) to obtain output in 1E-17 m**-1 s**-3
               V33(J,K)=10.**(4+NN)*2.*(DQU_DX+DQV_DY)
               IF(V33(J,K).LT.VMIN3) VMIN3=V33(J,K)
               IF(V33(J,K).GT.VMAX3) VMAX3=V33(J,K)
            ENDDO
         ENDDO
C
C        Writing Output Files
C        Checking max/min values:
         NESC=-17
         NN=0
         IF(VMAX3/10**NN.GT.999.) THEN
            DO N=1,10
               IF(VMAX3/10**(NN+N).LT.999.) GOTO 106
            ENDDO
  106       NN=NN+N
         ENDIF
         IF(VMIN3/10**NN.LT.-999.) THEN
            DO N=1,10
               IF(VMIN3/10**(NN+N).GT.-999.) GOTO 108
            ENDDO
  108       NN=NN+N
         ENDIF
C        Re-scaling values:
         IF(NN.GT.0) THEN
            NESC=NESC+NN
            FAC=10.**NN
            VMIN3=VMIN3/FAC
            VMAX3=VMAX3/FAC
            DO J=2*ND+1,NF-2*ND
               DO K=2*ND+1,NC-2*ND
                  V33(J,K)=V33(J,K)/FAC
               ENDDO
            ENDDO
         ENDIF
C        Writing:
         OPEN(UNIT=33,FILE='../varderiv/'//NA(I)//'qdi'//ME,
     .	                                 STATUS='UNKNOWN')
         IF(NOPFOR.EQ.1) THEN
            WRITE(33,32) (LIN(L),L=1,4)
            WRITE(33,42) VMIN3,VMAX3
            DO J=1,NF
               WRITE(33,34) (V33(J,K),K=1,NC)
C               WRITE(33,35)
            ENDDO
         ELSE
            DO J=1,NF
               WRITE(33,36) (V33(J,K),K=1,NC)
            ENDDO
         ENDIF
         WRITE(33,110) NESC
  110    FORMAT(//1X,'Units: 10**(',I3,') m**-1 s**-3') 
         CLOSE(UNIT=33)
C
      ENDDO
      CLOSE(UNIT=40)
      GOTO 150
C ----------------------------------------------------------------------
C     OPTION N: ...
C
C 120 
C ----------------------------------------------------------------------
  150 WRITE(6,152)
  152 FORMAT(5(/),8X,'============================================'//
     .8X,'NAME OF THE OUTPUT FILES:'/)
      DO I=1,NL
         GOTO(155,156,157,158,159) NOPVAR
  155    WRITE(6,162) NA(I)//'div'//ME
         GOTO 160
  156    WRITE(6,162) NA(I)//'rrv'//ME
         IF(NOPADV.EQ.1) WRITE(6,162) NA(I)//'rva'//ME
         GOTO 160
  157    WRITE(6,162) NA(I)//'dhu'//ME
         WRITE(6,162) NA(I)//'dhv'//ME
         GOTO 160
  158    WRITE(6,162) NA(I)//'grv'//ME
         IF(NOPADV.EQ.1) WRITE(6,162) NA(I)//'gva'//ME
         GOTO 160
  159    WRITE(6,162) NA(I)//'qdi'//ME
         IF(NOPQVE.EQ.1) THEN
            WRITE(6,162) NA(I)//'qvu'//ME
            WRITE(6,162) NA(I)//'qvv'//ME
         ENDIF
  160	ENDDO
  162 FORMAT(23X,A15)
  165 FORMAT(/8X,'============================================'///)
      WRITE(6,165)
      END



c     ==================================================================
c                  LOCAL POLINOMIAL FITTING ROUTINE
c
c     Fits a polynomial to a set of 3x3 points and computes first 
c     and second derivatives from the fitted parameters
c
c     ==================================================================
      SUBROUTINE POLFIT(DX,DY,V)
      REAL A(6,6),B(6,1),V(9),WK(36)
      INTEGER M,N,IA,IDGT,IER
C
      N=6
      M=1
      IA=6
      IDGT=5
      DO 5 I=1,6
      DO 5 J=1,6
    5 A(I,J)=0.
      A(1,1)=9.
      A(1,5)=6.*(DX**2)
      A(1,6)=6.*(DY**2)
      A(2,2)=6.*(DX**2)
      A(3,3)=6.*(DY**2)
      A(4,4)=4.*(DX**2)*(DY**2)
      A(5,5)=6.*(DX**4)
      A(5,6)=4.*(DX**2)*(DY**2)
      A(6,6)=6.*(DY**4)
      A(5,1)=A(1,5)
      A(6,1)=A(1,6)
      A(6,5)=A(5,6)
C
      B(1,1)=V(1)+V(2)+V(3)+V(4)+V(5)+V(6)+V(7)+V(8)+V(9)
      B(2,1)=(-V(1)+V(3)-V(4)+V(6)-V(7)+V(9))*DX
      B(3,1)=(V(1)+V(2)+V(3)-V(7)-V(8)-V(9))*DY
      B(4,1)=(-V(1)+V(3)+V(7)-V(9))*DX*DY
      B(5,1)=(V(1)+V(3)+V(4)+V(6)+V(7)+V(9))*DX**2
      B(6,1)=(V(1)+V(2)+V(3)+V(7)+V(8)+V(9))*DY**2
C     TYPE*,V(9),V9
C     -----------------------------------------------------------------
      DO 10 I=1,6
   10 V(I)=B(I,1)
C     Second Derivatives:
      V(5)=2*V(5)
      V(6)=2*V(6)
      DO 15 I=7,9
   15 V(I)=0
      RETURN
      END


c     ==================================================================
C                     FINITE DIFFERENCES SCHEME
c
c     Computes first and second derivatives by applying standard 
c     centred finite differences to a set of 3x3 points
c
c     ==================================================================
      SUBROUTINE FINDIF(DX,DY,V)
      DIMENSION V(9),VV(6)
C
      DO 5 I=1,9
    5 VV(1)=VV(1)+V(I)
      VV(1)=VV(1)/9.
      VV(2)=(V(6)-V(4))/(2*DX)
      VV(3)=(V(2)-V(8))/(2*DY)
      VV(4)=(V(3)+V(7)-V(9)-V(1))/(4*DX*DY)
      VV(5)=(V(4)+V(6)-2*V(5))/DX**2
      VV(6)=(V(2)+V(8)-2*V(5))/DY**2
      DO 10 I=1,6
   10 V(I)=VV(I)
      DO 15 I=7,9
   15 V(I)=0
      RETURN
      END

