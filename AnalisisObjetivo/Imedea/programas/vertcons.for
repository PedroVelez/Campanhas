C ======================================================================
C           
C                              Damia Gomis
C
C      IMEDEA (Institut Mediterrani d'Estudis AvanÇats), a joint 
C      centre between the Universitat de les Illes Balears (UIB) and 
C      the Spanish Research Council (CSIC).
C
C      This code has been compiled and distributed in the framework 
C      of the project REN2000-2599-E, funded by the Spanish Marine 
C      Science and Technology subprogram.
C
C          
C  
C ======================================================================
C
C                       Last version: 1 november 2001
C
C ======================================================================
C
      PARAMETER(NLMAX=100,NCMAX=100,NFMAX=100,NCC=100,NFF=100,NLL=100)
	IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NS=750,NG=1500,ND=20000)
      COMMON/INOUT/X(NS),Y(NS),V(ND),VS(ND),XG(NG),YG(NG),PAR(5)
      DIMENSION VG(NCC,NFF,NLL),VF(NCC,NLL),VC(NFF,NLL)
      CHARACTER*7 NAME,SYMBOL(NS)
      CHARACTER*4 IFOR(2)
      CHARACTER*2 ICVAR(8)
	CHARACTER*3 IVAR(11)
      CHARACTER*1 NOP,KC
	CHARACTER*4,NA0,NA1
	CHARACTER*2 LCN
	CHARACTER*8,NA(NLMAX)
      DATA IFOR/'.grd','.snp'/
      DATA ICVAR/'tc','sa','st','dh','uu','vv','pr','mp'/
      DATA IVAR/'div','rrv','rva','grv','gva','qvu','qvv','dhu','dhv',
     .                                           'qdi','dhw'/
      PI=3.14159
      RD=PI/180.
      XKGLA=60.*1.852
C
C ======================================================================
C                           INPUT DATA BLOCK
C ======================================================================
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
c*      WRITE(6,'(a,$)') 
c*   .' >>> DO YOU WANT TO USE *THIS* DOMAIN AND GRID ? ([Y]/N) '
c*      READ(5,13) NOP
c*   13 FORMAT(A1)
C ======================================================================
C                          READING STATION DATA 
C ======================================================================
   15 WRITE(6,'(a,$)') ' >>> ROOT OF THE GRID FILE(S) (NNNN): '
	READ(5,16,ERR=15) NA0
   16 FORMAT(A4)
	PMAX=P0+PINT*FLOAT(NL-1)
         DO I=1,NL

            LEV=INT(P0+PINT*FLOAT(I-1)+0.00001)

            IF(PMAX.GT.999) LEV=LEV/10
           NA1=CHAR(LEV/1000+48)
     .	//CHAR(LEV/100-10*(LEV/1000)+48)//
     .	CHAR(LEV/10-10*(LEV/100)+48)//CHAR(LEV-10*(LEV/10)+48)
              NA(I)=NA0//NA1
         ENDDO
   17 WRITE(6,*) ' >>> CHOOSE: 1-OBSERVATIONAL OR 2-DERIVATE VARIABLE:'
	READ(5,*,ERR=17) VA
	IF (VA.EQ.1) THEN 
		GOTO 20
	ELSE
		GOTO 22
	ENDIF
   20	WRITE(6,18)
   18 FORMAT(//1X,'*** VARIABLES:   1: temp   2: sal   3: st    4: dh',
     ./18X,'5: u      6: v     7: pres  8: mp')  
   19 WRITE(6,'(a,$)') ' >>> CHOOSE variable: '
      READ(5,*,ERR=19) NVAR
   21 WRITE(6,23)
   23 FORMAT(//1X,'*** STATION DATA WILL BE READ FROM AN EXTENAL FILE'/)
   25 WRITE(6,'(a,$)') ' >>> Choose Format option:  1-SURFER  2-SNP : '
      READ(5,*) NSOR
	IF(NSOR.EQ.1) THEN
		do i=1,nl
			WRITE(6,*)'../varobs/'//na(I)//icvar(nvar)//'.grd'
			open(unit=2, file='../varobs/'//na(I)//icvar(nvar)//'.grd',
     .		                                            status='old')
C	Read data in SURFER (.grd) format
		read(2,27) nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
			do j=1,nf
				READ(2,31) (vg(l,j,i), l=1,nc)
C				READ(2,33)
			enddo
   27	format('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f8.3))
   31	format(10f10.4)
   33	format(1x)
		close (unit=2)
		enddo
	goto 36
	ELSE IF(NSOR.EQ.2) THEN
c        Read data in SNP format 
	do i= 1,nl
		WRITE(6,*) '../varobs/'//na(I)//icvar(nvar)//'.snp'
		open(unit=3, file='../varobs/'//na(I)//icvar(nvar)//'.snp',
     .		                                            status='old')
         DO j=1,nf
             READ(3,35) (vg(l,j,i), l=1,nc)
         ENDDO
   35    FORMAT(8F10.4)
      close(unit=3)
	enddo
	goto 36
	ELSE
	  goto 25
	ENDIF
	goto 36
   22 WRITE(6,24)
   24 FORMAT(//1X,'*** VARIABLES: 1:div  2: rrv  3: rva  4:  grv
     . 5: gva  6: gvu  7: gvv  8: dhu  9: dhv  10: qdi  11: dhw')
   26 WRITE(6,'(a,$)') ' >>> CHOOSE variable: '
      READ(5,*,ERR=26) NOPVAR
      WRITE(6,30)
   30 FORMAT(//1X,'*** STATION DATA WILL BE READ FROM AN EXTENAL FILE'/)
   32 WRITE(6,'(a,$)') ' >>> Choose Format option:  1-SURFER  2-SNP : '
      READ(5,*) NSOR
	
	IF(NSOR.EQ.1) THEN
		do i=1,nl
			WRITE(6,*) '../varderiv/'//na(I)//ivar(nopvar)//'.grd'
			open(unit=2, file='../varderiv/'//na(I)//ivar(nopvar)//
     .		                                     '.grd',status='old')
C	Read data in SURFER (.grd) format
		read(2,27) nc,nf,xmin,xmax,ymin,ymax,zmin,zmax
			do j=1,nf
				READ(2,31) (vg(l,j,i), l=1,nc)
C				READ(2,33)
			enddo
		close (unit=2)
		enddo
	ELSE IF(NSOR.EQ.2) THEN
c        Read data in SNP format 
	do i= 1,nl
		WRITE(6,*) '../varderiv/'//na(I)//ivar(nopvar)//'.snp'
		open(unit=3, file='../varderiv/'//na(I)//ivar(nopvar)//
     .		                                  '.snp', status='old')
         DO j=1,nf
             READ(3,35) (vg(l,j,i), l=1,nc)
         ENDDO
      close(unit=3)
	enddo
	goto 36
	ELSE
	  goto 32
	ENDIF
   36	write(6,'(a,$)') 'DO YOU WANT TO MAKE A PROFILE FOLLOWING ROWS OR 
     .COLUMNS? 1-ROW  2-COLUMN : '
	read (5,*) OPT
	IF(OPT.EQ.1) THEN
		write(6,'(a,$)') 'WHAT ROW?'
		read(5,*) lfil
		do i=1,nl
				do j= lfil,lfil	
					do l=1,nc	
					vf(l,i)= vg(l,j,i)
					enddo
				enddo
		enddo
		lcn=CHAR(lfil/10-10*(lfil/100)+48)//CHAR(lfil-10*(lfil/10)+48)
		kc='r'
	else
		write(6,'(a,$)') 'WHAT COLUMN?'
		read(5,*) lcol
		do i=1,nl
				do j=1,nf
					do l=lcol,lcol
					vc(j,i)= vg(l,j,i)
					enddo
				enddo
		enddo
		
		lcn=CHAR(lcol/10-10*(lcol/100)+48)//CHAR(lcol-10*(lcol/10)+48)
		kc='c'
	endif
	
      
    
C
C ======================================================================
C                       WRITING DATA TO OUTPUT FILE
C ======================================================================
      WRITE(6,43)
   43 FORMAT(//,1X,'*** WRITING OUTPUT FILES: '/)
   45 WRITE(6,'(a,$)') ' >>> Choose Format option:  1-SURFER  2-SNP : '
      READ(5,*,ERR=45) NSOR

	if (va.eq.1) then
		goto 47
	else
		goto 49
	endif

   47	open(unit=99,file='../vertical/'//na0//kc//lcn//icvar(nvar)
     .                    //ifor(nsor),status='unknown')
	
	WRITE(6,*) ' '
	WRITE(6,*) '   OUTPUT FILE:',na0//kc//lcn//icvar(nvar)//ifor(nsor)
     .                    
	                                 
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
				if(vf(l,i).gt.zmax) zmax=vf(l,i)
				if(vf(l,i).lt.zmin) zmin=vf(l,i)
			enddo
	   enddo
		write(99,130) nc,nl,xmin,xlmax+xmin,-(zlmax+ymin),-ymin,
     .                                               	   zmin,zmax
			do i=1,nl
				write(99,135) (vf(l,i), l=1,nc)
C				write(99,137)
			enddo

  130    format('DSAA'/2(1x,i3)/2(1x,f7.2)/2(1x,f7.2)/2(1x,f9.4))
  135    format(10f10.4)
  137    format(1x)
         else
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
				if(vc(j,i).gt.zmax) zmax=vc(j,i)
				if(vc(j,i).lt.zmin) zmin=vc(j,i)
			enddo
	   enddo
		write(99,130) nf,nl,xmin,ylmax+xmin,-(zlmax+ymin),-ymin,
     .                                           		   zmin,zmax
				do i=1,nl
					write(99,135) (vc(j,i), j=1,nf)
C					write(99,137)
				enddo


		endif
      ELSE IF(NSOR.EQ.2) THEN
C        Data plotted in SNP format 
		IF(OPT.EQ.1) THEN
		do i=1,nl
			write(99,140) (vf(l,i),l=1,nc)
		enddo
		else
		do i=1,nl
			write(99,140) (vc(j,i), j=1,nf)
		enddo
  140    FORMAT(8F10.4)
         endif
C
      ELSE
	  goto 45
      ENDIF

    
	goto 50	
   49	open(unit=99,file='../vertical/'//na0//kc//lcn//ivar(nopvar)//
     .                                      ifor(nsor),status='unknown')
	WRITE(6,*) '  OUTPUT FILE:',na0//kc//lcn//ivar(nopvar)//ifor(nsor)

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
				if(vf(l,i).gt.zmax) zmax=vf(l,i)
				if(vf(l,i).lt.zmin) zmin=vf(l,i)
			enddo
	   enddo
		write(99,130) nc,nl,xmin,xlmax+xmin,-(zlmax+ymin),-ymin,
     .                                                   zmin,zmax
			do i=1,nl
				write(99,135) (vf(l,i), l=1,nc)
C				write(99,137)
			enddo
           else
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
				if(vc(j,i).gt.zmax) zmax=vc(j,i)
				if(vc(j,i).lt.zmin) zmin=vc(j,i)
			enddo
	   enddo
		write(99,130) nf,nl,xmin,xlmax+xmin,-(zlmax+ymin),-ymin,
     .                                                   zmin,zmax
				do i=1,nl
					write(99,135) (vc(j,i), j=1,nf)
C					write(99,137)
				enddo
		endif
      ELSE IF(NSOR.EQ.2) THEN
C        Data plotted in SNP format 
		IF(OPT.EQ.1) THEN
		do i=1,nl
			write(99,140) (vf(l,i),l=1,nc)
		enddo
		else
		do i=1,nl
			write(99,140) (vc(j,i), j=1,nf)
		enddo
         endif
C
      ELSE
	  goto 45
      ENDIF
   50 CLOSE(UNIT=99)
C
	end