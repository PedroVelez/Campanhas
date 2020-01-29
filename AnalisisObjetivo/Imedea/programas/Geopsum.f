C ====================================================================
C                               GEOPSUM.F
C ====================================================================
C
C   THIS CODE COMPUTES THE POTENTIAL FUNCTION (E.G., DYNAMIC HEIGHT)
C   AT EACH LEVEL BY ADDING THE PARTIAL CONTRIBUTIONS OF ALL THE 
C   LAYERS BELOW. IT MUST BE APPLIED TO TRANSFORM THICKNESS VALUES 
C   RESULTING FROM OPTION 2 OF 'C_PRES.F' INTO DYNAMIC HEIGHT VALUES.
C                   
C ====================================================================
C           
C                         AUTHOR: Damià Gomis
C
C       IMEDEA (Institut Mediterrani d'Estudis AvanÇats), a joint 
C    research centre between the Universitat de les Illes Balears 
C    (UIB) and the Spanish Research Council (CSIC). Mallorca (Spain)
C
C    This code has been revised and distributed in the framework of
C    the project REN2000-2599-E, funded by the Spanish Marine Science
C    and Technology subprogram. 
C  
C ====================================================================
C
C   Work structure:
C   - the grid specifications are read from 'grid.dat'
C   - the whole set of TH analysis are read and stored in a 3D matrix
C   - a cummulative sum is applied from the lowest to the highest 
C     level. In this way, DYNH is computed at each level as the 
C     contribution of TH values obtained at all the levels below. 
C   - the whole set of DYNH level files is written on output.
C
C   External files used:  
C   - INPUT files:
C     * domain info file:                           [../info/grid.dat]
C     * gridpoint TH data files:               [../snp/NNNNLLLLth.snp] 
C                                           or [../grd/NNNNLLLLth.grd] 
C   - OUTPUT files:
C     * gridpoint DH data files:               [../snp/NNNNLLLLdh.snp]  
C                                           or [../grd/NNNNLLLLdh.grd] 
C
C   Built-in routines:  none
C   External routines:  none
C
C ====================================================================
      PARAMETER(NLMAX=100,NFMAX=100,NCMAX=100)
      DIMENSION GEOP(NLMAX,NFMAX,NCMAX)
      CHARACTER*4,NA,ME,CAB,NA1
	CHARACTER*8,STUDA
      CHARACTER*15,LIN1
      CHARACTER*63,LIN2
C --------------------------------------------------------------------
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
	read(99,*) NA1,STUDA
    3 format(1x/1x)
      close(unit=99)
c
    4 WRITE(6,5)
    5 FORMAT(5X,'This program generates SNP or GRD files.',/
     .5X,'The program reads thickness files (th)',/
     .5X,'and creates dynamic height (dh) files',//
     .5X,'Input root of thickness files name [NNNN****.th] (a4):')
      READ(5,6,ERR=4) NA
    6 format(a4)
	
    7	WRITE(6,'(a,$)') ' >>> FILE FORMAT:  1: *.GRD   2: *.SNP  : '
      READ(5,*,ERR=7) NOPFOR
      IF(NOPFOR.EQ.1) THEN
         ME='.grd'
      ELSE IF(NOPFOR.EQ.2) THEN
         ME='.snp'
      ELSE
         GOTO 7
      ENDIF

	WRITE(6,8)
    8 FORMAT('Input files:',2(/))
	DO 10 I=1,NL
	
	p00=p0+pint*float(I-1)
	lev=ifix(p00)

	OPEN(UNIT=1,FILE='../varobs/'//NA//
     .	char(lev/1000+48)//
     .	char(lev/100-10*(lev/1000)+48)//     
     .	char(lev/10-10*(lev/100)+48)//
     .	char(lev-10*(lev/10)+48)//'th'//ME,status='OLD')
	
	WRITE(6,*)'>>',NA//char(lev/1000+48)//
     .	char(lev/100-10*(lev/1000)+48)// 
     .	char(lev/10-10*(lev/100)+48)//    
     .	char(lev-10*(lev/10)+48)//'th'//ME
C     Summing the values of all the levels below. A new version of the 
C     input file is created

C	Read data in SURFER (.grd) format

	IF(NOPFOR.EQ.1) THEN

		read(1,26) cab,nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
C		type*,cab,nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
			do J=1,NF
				READ(1,28) (GEOP(I,J,K), K=1,NC)
C				READ(1,29)
			enddo
		CLOSE(UNIT=1)
	goto 10

	ENDIF
   
C	Read data in SNP (.snp) format
	       
	   DO J=1,NF
            READ(1,20) (GEOP(I,J,K),K=1,NC)
         ENDDO

	CLOSE(UNIT=1)

   10 CONTINUE
  
C
      DO I=(NL-1),1,-1
         DO J=1,NF
            DO K=1,NC
               GEOP(I,J,K)=GEOP(I,J,K)+GEOP(I+1,J,K)
            ENDDO
         ENDDO
      ENDDO
C	
	zmax=GEOP(1,1,1)
	zmin=GEOP(1,1,1)
	DO I=(NL-1),1,-1
		DO J=1,NF
			DO K=1,NC
	          IF(GEOP(I,J,K).GT.ZMAX) ZMAX=GEOP(I,J,K)
			  IF(GEOP(I,J,K).LT.ZMIN) ZMIN=GEOP(I,J,K)
			ENDDO
		ENDDO
	ENDDO	  	

C
	WRITE(6,11)
   11 FORMAT(2(/),'Output files:',2(/))
	DO 16 I=1,NL
         p00=p0+pint*float(I-1)
	   lev=ifix(p00)
       
C	-----------------------------------------------------------------
	
	OPEN(unit=3,file='../varobs/'//NA//
     .	char(lev/1000+48)//
     .	char(lev/100-10*(lev/1000)+48)//
     .	char(lev/10-10*(lev/100)+48)//
     .	char(lev-10*(lev/10)+48)//'dh'//ME,status='NEW')
	  
	WRITE(6,*) '>>',NA//char(lev/1000+48)//
     .	char(lev/100-10*(lev/1000)+48)//
     .	char(lev/10-10*(lev/100)+48)//
     .	char(lev-10*(lev/10)+48)//'dh'//ME 
	
C	Formato (*.GRD)


	IF(NOPFOR.EQ.1) THEN
         WRITE(3,27) nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
         DO J=1,NF
            WRITE(3,28) (GEOP(I,J,K),K=1,NC)
C            WRITE(3,29)
         ENDDO

	CLOSE(UNIT=3)
	
	goto 16
	ENDIF

	   DO J=1,NF
            WRITE(3,20) (GEOP(I,J,K),K=1,NC)
         ENDDO

	CLOSE(UNIT=3)	
   16 CONTINUE

   20 FORMAT(8F10.4)
   26	FORMAT(A4/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f8.3))
   27	FORMAT('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f8.3))
   28	FORMAT(10f10.4)
   29	FORMAT(1x)	
      END
