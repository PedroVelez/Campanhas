C ====================================================================
C                             VERTCONS.F
C ====================================================================
C
C   THIS CODE SELECTS A VERTICAL SECTION (ALONG A ROW OR A COLUMN) 
C   FROM A SET OF HORIZONTAL (2D) ANALYSIS.
C
C   THE OUTPUT IS WRITTEN IN THE SAME FORMAT AS THE ORIGINAL (2D)
C   FILES FROM WHICH THE SECTION IS EXTRACTED.
C
C ====================================================================
C           
C            AUTHORS: Rocío Fernández, Simón Ruiz, Damià Gomis
C
C       IMEDEA (Institut Mediterrani d'Estudis AvanÇats), a joint 
C     research centre between the Universitat de les Illes Balears 
C    (UIB) and the Spanish Research Council (CSIC). Mallorca (Spain)
C
C    This code has been revised and distributed in the framework of
C    the project REN2000-2599-E, funded by the Spanish Marine Science
C    and Technology subprogram. 
C  
C ====================================================================
C
C
C   Work structure:
C   - info on the domain is first read
C   - the user selects among the diferent options (location of the
C     vertical section, format of the input/output files,...)
C   - the whole set of horizontal gridpoint data files conforming the  
C     3D domain are read and saved in a 3D matrix
C   - data over the selected section are extracted and written in the
C     output file.
C
C   External files used:  
C   - INPUT files:
C     * domain info file:                           [../info/grid.dat]
C     * set of horizontal gridpoint data files 
C       conforming the 3D domain:              [../snp/NNNNLLLLVV.snp] 
C                                           or [../grd/NNNNLLLLVV.grd] 
C   - OUPUT files: 
C     * vertical section along a row:      [../vertical/NNNNrRRVV.snp]  
C                                       or [../vertical/NNNNrRRVV.grd] 
C                  or along a column:      [../vertical/NNNNcCCVV.snp] 
C                                       or [../vertical/NNNNcCCVV.grd] 
C
C   [For a derived variable, input and output file names would be the 
C    same, but with VVV replaciong VV]
C
C   Built-in routines:  none
C   External routines:  none
C
C ====================================================================
      PARAMETER(NLMAX=150,NCMAX=100,NFMAX=100,NS=750,NG=1500)
	IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VG(NCMAX,NFMAX,NLMAX),VF(NCMAX,NLMAX),VC(NFMAX,NLMAX)	
      CHARACTER*1 KC
	CHARACTER*2 ICVAR(8),LCN
	CHARACTER*3 IVAR(11)
	CHARACTER*4,NA0,NA1,IFOR(2)
	CHARACTER*8,NA(NLMAX)
      DATA IFOR/'.grd','.snp'/
      DATA ICVAR/'tc','sa','st','dh','uu','vv','pr','mp'/
      DATA IVAR/'div','rrv','rva','grv','gva','qvu','qvv','dhu','dhv',
     .                                                    'qdi','dhw'/
      PI=3.14159
      RD=PI/180.
      XKGLA=60.*1.852
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
      read(99,*) p0,pINT
    3 format(1x/1x)
      close(unit=99)
c
      ALF0=ALF0*RD
      NGRI=NF*NC
C     Coordinates of the other corners (degrees):
    5 XLMAX=XARM*REAL(NC-1)
      YLMAX=YARM*REAL(NF-1)
	ZLMAX=PINT*REAL(NL-1)
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
         WRITE(6,9) XARM,YARM,XLAT2,XLON2,XLAT3,XLON3,
     .	   XLAT1,XLON1,XLAT4,XLON4
      ELSE
         XLON4=XLON1+XLMAX
         XLAT4=XLAT1
         XLON3=XLON4
         XLAT3=XLAT4+YLMAX
         XLON2=XLON1
         XLAT2=XLAT1+YLMAX
         WRITE(6,11) XARM,YARM,XLAT2,XLON2,XLAT3,XLON3,
     .	   XLAT1,XLON1,XLAT4,XLON4
      ENDIF
    9 FORMAT(5x,'of ',f5.1,'  x ',f5.1,'  km. The LIMITS of the domain',
     .' are:'//2(2(5x,'lat =',f7.3,'  lon =',f8.3)/))
   11 FORMAT(5x,'of',f5.3,'  x ',f5.3,' deg. The LIMITS of the domain',
     .' are:'//2(2(5x,'lat =',f7.3,'  lon =',f8.3)/))
C
C
   15 WRITE(6,'(a,$)') ' >>> ROOT OF THE GRID FILES (NNNN): '
	READ(5,16,ERR=15) NA0
   16 FORMAT(A4)
	PMAX=P0+PINT*FLOAT(NL-1)
      DO I=1,NL
         LEV=INT(P0+PINT*FLOAT(I-1)+0.00001)
         NA1=CHAR(LEV/1000+48)//CHAR(LEV/100-10*(LEV/1000)+48)//
     .   CHAR(LEV/10-10*(LEV/100)+48)//CHAR(LEV-10*(LEV/10)+48)
         NA(I)=NA0//NA1
      ENDDO
C
   17 WRITE(6,*) ' >>> CHOOSE: 1-OBSERVED VARIABLE  2-DERIVED VARIABLE:'
	READ(5,*,ERR=17) VA
C
	IF(VA.EQ.1) THEN
	   WRITE(6,18)
   18    FORMAT(//1X,'*** VARIABLES:   1: temp   2: sal   3: st   ',
     .   ' 4: dh'/18X,'5: u      6: v     7: pres  8: mp')  
   19    WRITE(6,'(a,$)') ' >>> CHOOSE variable: '
         READ(5,*,ERR=19) NVAR
   25    WRITE(6,'(a,$)') '>>> Input files format: 1: *.GRD  2: *.SNP: '
         READ(5,*) NSOR 
	   IF(NSOR.EQ.1) THEN
C	      Read data in SURFER (.grd) format
		  do i=1,nl
			 WRITE(6,*)'../varobs/'//na(I)//icvar(nvar)//'.grd'
			 open(unit=2,file='../varobs/'//na(I)//icvar(nvar)//
     .		                                '.grd',status='old')
		     read(2,27) nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
			 do j=1,nf
				READ(2,31) (vg(l,j,i), l=1,nc)
C				READ(2,33)
			 enddo
   27	         format('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f8.3))
   31	         format(10f10.4)
   33	         format(1x)
		     close (unit=2)
		  enddo
	   ELSE IF(NSOR.EQ.2) THEN
c           Read data in SNP format 
	      do i= 1,nl
		     WRITE(6,*) '../varobs/'//na(I)//icvar(nvar)//'.snp'
		     open(unit=2,file='../varobs/'//na(I)//icvar(nvar)//
     .		                                '.snp',status='old')
               do j=1,nf
                  READ(2,35) (vg(l,j,i), l=1,nc)
               enddo
   35          format(8F10.4)
               close(unit=2)
	      enddo
 	   ELSE
	      goto 25
	   ENDIF
	ELSE IF(VA.EQ.2) THEN
   22    WRITE(6,24)
   24    FORMAT(//1X,'*** VARIABLES: 1:div  2: rrv  3: rva  4:  grv',
     .   '5: gva  6: gvu  7: gvv  8: dhu  9: dhv  10: qdi  11: dhw')
   26    WRITE(6,'(a,$)') ' >>> CHOOSE variable: '
         READ(5,*,ERR=26) NOPVAR
   32    WRITE(6,'(a,$)') ' >>> Choose Format option: 1-SURFER  2-SNP :'
         READ(5,*) NSOR
	   IF(NSOR.EQ.1) THEN
C	      Read data in SURFER (.grd) format
	 	  do i=1,nl
		 	 WRITE(6,*) '../varderiv/'//na(I)//ivar(nopvar)//'.grd'
			 open(unit=2,file='../varderiv/'//na(I)//ivar(nopvar)//
     .		                                   '.grd',status='old')
	   	     read(2,27) nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
			 do j=1,nf
				READ(2,31) (vg(l,j,i), l=1,nc)
C				READ(2,33)
			 enddo
		     close (unit=2)
		  enddo
	   ELSE IF(NSOR.EQ.2) THEN
c           Read data in SNP format 
	      do i=1,nl
		     WRITE(6,*) '../varderiv/'//na(I)//ivar(nopvar)//'.snp'
		     open(unit=2,file='../varderiv/'//na(I)//ivar(nopvar)//
     .		                                   '.snp',status='old')
               do j=1,nf
                  READ(2,35) (vg(l,j,i), l=1,nc)
               enddo
			 close(unit=2)
	      enddo
	   ELSE
	      goto 32
	   ENDIF
      ENDIF
C
C ======================================================================
C                       EXTRACTING THE VERTICAL SECTION
C ======================================================================
C
   	write(6,'(a,$)') 'DO YOU WANT A VERTICAL SECTION ALONG A ROW (1)  
     .OR ALONG A COLUMN (2) ? : '
	read (5,*) OPT
	IF(OPT.EQ.1) THEN
	   write(6,'(a,$)') 'WHICH ROW number (1=south; NF=north) ? '
	   read(5,*) lfil
	   do i=1,nl	
		  do l=1,nc	
			 vf(l,i)=vg(l,lfil,i)
		  enddo
	   enddo
	   lcn=CHAR(lfil/10-10*(lfil/100)+48)//CHAR(lfil-10*(lfil/10)+48)
	   kc='r'
	ELSE
	   write(6,'(a,$)') 'WHICH COLUMN number (1=west; NC=east) ? '
	   read(5,*) lcol
	   do i=1,nl
		  do j=1,nf
			 vc(j,i)= vg(lcol,j,i)
		  enddo
	   enddo	
	   lcn=CHAR(lcol/10-10*(lcol/100)+48)//CHAR(lcol-10*(lcol/10)+48)
	   kc='c'
	endif
C
C ======================================================================
C                       WRITING DATA TO OUTPUT FILE
C ======================================================================
C
      WRITE(6,43)
   43 FORMAT(//,1X,'*** WRITING OUTPUT FILES: '/)
	if(va.eq.1) then
         open(unit=99,file='../vertical/'//na0//kc//lcn//icvar(nvar)//
     .                                    ifor(nsor),status='unknown')	
	   WRITE(6,*) ''
	   WRITE(6,*) '   OUTPUT FILE:',na0//kc//lcn//icvar(nvar)//
     .                                                 ifor(nsor)
      else
	   open(unit=99,file='../vertical/'//na0//kc//lcn//ivar(nopvar)//
     .                                     ifor(nsor),status='unknown')
	   WRITE(6,*) ''
	   WRITE(6,*) '  OUTPUT FILE:',na0//kc//lcn//ivar(nopvar)//
     .                                                 ifor(nsor)
      endif
C
      IF(NSOR.EQ.1) THEN
C        Data plotted in SURFER (.grd) format
	   IF(OPT.EQ.1) THEN
		  if(alf0.gt.0.) then
			 xmin=0.
			 ymin=0.
		  else
			 xmin=xlon1
			 ymin=p0
		  endif
		  zmin=1E6
	      zmax=-1E6
            do i=1,nl
		     do l=1,nc
	            if(vf(l,i).gt.zmax) then
	               if(va.eq.1) then
		              zmax=vf(l,i)
                     else
		              if(vf(l,i).lt.998) zmax=vf(l,i)
	               endif
	            endif
			    if(vf(l,i).lt.zmin) zmin=vf(l,i)
		     enddo
	      enddo
		  write(99,130) nc,nl,xmin,xmin+xlmax,-(ymin+zlmax),-ymin,
     .                                               	    zmin,zmax
		  do i=nl,1,-1
			 write(99,135) (vf(l,i), l=1,nc)
C			 write(99,137)
		  enddo
  130       format('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f8.2)/2(1x,f9.4))
  135       format(10f10.4)
  137       format(1x)
         ELSE
		  if(alf0.gt.0.) then
			 xmin=0.
			 ymin=0.
		  else
			 xmin=xlat1
			 ymin=p0
		  endif
		  zmin=1E6
	      zmax=-1E6
            do i=1,nl
		     do j=1,nf
	            if(vc(j,i).gt.zmax) then
	               if(va.eq.1) then
		              zmax=vc(j,i)
                     else
		              if(vc(j,i).lt.998) zmax=vc(j,i)
	               endif
	            endif
			    if(vc(j,i).lt.zmin) zmin=vc(j,i)
		     enddo
	      enddo
		  write(99,130) nf,nl,xmin,xmin+ylmax,-(ymin+zlmax),-ymin,
     .                                           		    zmin,zmax
		  do i=nl,1,-1
			 write(99,135) (vc(j,i), j=1,nf)
C			 write(99,137)
		  enddo
         ENDIF
      ELSE IF(NSOR.EQ.2) THEN
C        Data plotted in SNP format 
	   IF(OPT.EQ.1) THEN
		  do i=nl,1,-1
			 write(99,140) (vf(l,i),l=1,nc)
		  enddo
	   ELSE
		  do i=nl,1,-1
			 write(99,140) (vc(j,i), j=1,nf)
		  enddo
  140       FORMAT(8F10.4)
          ENDIF
      ENDIF
C
      CLOSE(UNIT=99)
	end
