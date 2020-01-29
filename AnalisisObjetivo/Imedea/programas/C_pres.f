C ====================================================================
C                               C_PRES.F
C ====================================================================
C
C   THIS CODE CONSTRUCTS PRESSURE-LEVEL DATA FILES (*.LEV) FROM THE
C   SMOOTHED STATION-PROFILE DATA FILES (*.SMO). OBSERVATIONS AT THE 
C   SPECIFIED LEVELS ARE SEARCHED FOR IN THE PROFILES AND WRITTEN 
C   ALTOGETHER IN THE OUTPUT FILE.
C   THE CODE ALSO COMPUTES DYNAMIC HEIGHT AT THE SPECIFIED LEVEL WITH
C   RESPECT TO A REFERENCE LEVEL SPECIFIED BY THE USER.
C
C ====================================================================
C           
C                         AUTHOR: Damià Gomis
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
C   Work structure:
C   - the user selects between building up just one level file (to be 
C     determined by the user) or one for each level specified in the 
C     domain data file 'grid.dat'
C   - the CTD profiles to be used are read from the file 'ctd.lis',  
C     and the velccity profiles are read from the file 'veloc.lis'.
C   - for each profile, p is read downwards until reaching the 
C     specified level; at that point p, T, S and DENS are read from 
C     the profile and written in the output file. If a profile does 
C     not contain values exactly at the specified level, these are 
C     computed as a linear interpolation between the closest profile 
C     levels above and below.
C   - the same holds for variables u,v of velocity profiles
C   - Dynamic Height (DYNH) is also computed from profile data and 
C     written in the output file.
C     * Normally, DYNH is computed with respect to a fixed reference 
C       level (option 1 in this code). In this case, if a given 
C       profile is shallower than the selected reference level, DYNH 
C       cannot be computed at that station.
C     * Another option is to use the level immediately below (as 
C       defined in the file 'grid.dat') as a reference level (option 
C       2). In that case (referred as TH after 'thickness'), if a 
C       given profile is shallower than the vertical domain defined in
C       'grid.dat', TH can still be computed at several levels. DYNH 
C       can be retrieved afterwards, once TH has been interpolated 
C       onto a grid, simply by adding the values of TH at all the 
C       levels below [see the help document for further details about 
C       this option.]
C
C   External files used:  
C   - INPUT files:
C     * list of CTD profiles:                        [../info/ctd.lis] 
C     * list of velocity profiles:                 [../info/veloc.lis] 
C     * domain info file:                           [../info/grid.dat] 
C     * revised CTD profiles:                     [../smo/NNNNnnn.smo]
C     * velocity profiles:                        [../vel/NNNNnnn.vel]
C   - OUTPUT files:
C     * level data file(s):                      [../lev/NNNNLLLL.lev]
C
C   Built-in routines:  FUNCTION SVAN
C   External routines:  none
C			      
C ====================================================================
C
	character*3,na2
	character*4,na1,star
	character*1,na3
	character*8,studa
C
C    =================================================================
C                         	input data block
C    =================================================================
C
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
	read(99,3)
	read(99,*) NA1,STUDA
    3 format(1x/1x)
      close(unit=99)
C
	write(6,4)
    4	format(2X,'THIS PROGRAM CONSTRUCTS PRESSURE-LEVEL FILES *.LEV'/
     .2X,'FROM THE SMOOTHED STATION-PROFILE FILES *.SMO AND'/
     .2X,'THE VELOCITY PROFILE FILES *.VEL'/)     
	write(6,5)
    5	format(/'Select option:'//5x,'>>>Individual level (1)'/
     .                          5x,'>>>All levels (2)')                                    
	read(5,*) LE
	if(le.eq.1) then
	   write(6,'(a,$)') 'Input level (meters)? : '
	   read(5,*) p0
	   NL=1
	endif
	write(6,6)
    6	format(/'Select option:'//5X,'>>>Fixed reference level (1)'/
     .           5X,'>>>Level just below (see grid.dat file) (2)')
	read(5,*) RE_RLEV
	if(re_rlev.eq.1) then                                
	   write(6,'(a,$)') 'Depth of reference level ? : '
	   read(5,*) PREF
	endif
C
C    =================================================================
C     Beginning the loop when more than one level is to be constructed
C
     	do 210 k=1,nl
C
      p00=p0+pint*float(k-1)
	lev=ifix(p00)          
	if(re_rlev.eq.2) then
	pref=p00+pINT
	endif
C
C    =================================================================
C                            READING CTD DATA
C    =================================================================
C
C ....opening the output file and designing the file headlines
      indctd=0
	open(unit=99,file='../info/ctd.lis',status='old',err=100)
	read(99,12)
   	read(99,37,end=80) na1,na3,na2
      rewind(99)
   12	format(//)        
C
	open(unit=2,file='../level/'//na1//char(lev/1000+48)//
     . char(lev/100-10*(lev/1000)+48)//char(lev/10-10*(lev/100)+48)//
     . char(lev-10*(lev/10)+48)//'.lev',status='unknown')
	write(2,20) na1,studa
	write(2,25) p00,pref
   20	format(1x,'***     DATA ON ISOBARIC SURFACES READY TO BE ',
     .	'PROCESSED     ***'/1x,'***',5x,'CRUISE CODE: -',a4,
     .  '- .   Region:',a8,5x,'***'/)
   25	format(24x,'level:',f7.2/16x,'pr-ref. level:',
     .  f7.2//22x,'<',12x,'DATA TO BE USED',13x,'>       <REMARKS>'
     .  /'---------',(1x,'-------'),' --------',3(1x,'------'),1x,
     .'-------',2(1x,'------') ,'  ---------'/1x,'  A8 ',2x,(3x,'F7.3')
     .,2x,(3x,'f8.3'),1x,3(3x,'f6.3'),'  dyn. cm',2x,'cm/s',3x,
     .'cm/s',/1x, 'SYMBOL   LATITU   LONGIT  TEMPER SALINI DENSIT 
     .DYN.HEI   U',6x,'V'/'---------',(1x,'-------'),' --------',
     . 3(1x,'------'),1x,'-------',2(1x,'------'), '  ---------')
C
      indctd=1
C    ------------------------------------------------------------------
C ....start reading
C
	read(99,12)
   35	read(99,37,end=80) na1,na3,na2
   37 format(a4,a1,a3)
   38 format(a4,a3)
	write(6,39) p00,na1,na2
   39	format(1x,f6.1,2x,a4,a3,'.smo')
C
	open(unit=1,file='../smo/'//na1//na2//'.smo',status='old') 
      read(1,40) latg0,latm0,lats0,long0,lonm0,lons0
C	write(6,40) latg0,latm0,lats0,long0,lonm0,lons0
   40 format(8x,i3,1x,i2,1x,f5.2,8x,i3,1x,i2,1x,f5.2)  
    	xlats0=0
	xlat0=float(abs(latg0))+float(latm0)/60.+xlats0/3600.
	if(latg0.lt.0.) xlat0=-xlat0
	xlons0=0
	xlon0=float(abs(long0))+float(lonm0)/60.+xlons0/3600.
	if(long0.lt.0) xlon0=-xlon0
	idep=0
C	write(6,89) xlat0,xlon0  
C
C	initial values are read...
   44	read(1,88,end=75) p2,t2,s2,st2,star,sst2
	if(star.eq.' * ') st2=sst2
	if((t2.gt.99.).or.(s2.gt.99.)) goto 44
C	warning: sometimes there are no data just below the surface
	if(p2.gt.p00) goto 75
   45	p1=p2
	t1=t2
	s1=s2
	st1=st2
C	when the chosen level is deeper than the profile, the station  
C	will simply not be included in the data set
   46	read(1,88,end=75) p2,t2,s2,st2,star,sst2
	if(star.eq.' * ') st2=sst2
	if((t2.gt.99).or.(s2.gt.99)) goto 46
	if(p2.lt.p00) goto 45
	van2=svan(s2,t2,p2,sigma)
	t0=t1+(t2-t1)*(p00-p1)/(p2-p1)
	s0=s1+(s2-s1)*(p00-p1)/(p2-p1)
	st0=st1+(st2-st1)*(p00-p1)/(p2-p1)
	van0=svan(s0,t0,p00,sigma)
C	integrating down to the next level...
	sum=(p2-p00)*(van0+van2)
C	warning: both the chosen level and the reference level could 
C	be in between the same two data of the vertical profile
	if(p2.ge.pref) goto 55
   50	p1=p2
	t1=t2
	s1=s2
	van1=van2
   51	read(1,*,end=60) p2,t2,s2
      if((t2.gt.99).or.(s2.gt.99)) goto 51
	van2=svan(s2,t2,p2,sigma)
	sum=sum+(p2-p1)*(van1+van2)
	if(p2.lt.pref) goto 50
   55	tref=t1+(t2-t1)*(pref-p1)/(p2-p1)
	sref=s1+(s2-s1)*(pref-p1)/(p2-p1)
	vanref=svan(sref,tref,pref,sigma)
	sum=sum-(p2-pref)*(van2+vanref)
C
C ....we plot 100*geo0 in order to obtain dynamic centimeters:
	geo0=.5*1E-3*sum
	goto 70
C    ------------------------------------------------------------------
C ....when the reference level is below the bottom, the geopotential
C	is set equal to: 
   60	geo0=999.999
C    ------------------------------------------------------------------
   70	write(2,90) na1,na2,xlat0,xlon0,t0,s0,st0,geo0
   75	close(unit=1)
	goto 35
   80	close(unit=99)
   88 format(1x,f6.1,3f7.3,a3,f7.3)
   89 format(2x,'LATITUD',2x,'LONGITUD',2x,f7.3,2x,f8.3/)
   90	format(1x,a4,'c',a3,1x,f7.3,1x,f8.3,3(1x,f6.3),1x,f7.3,
     .       2(1x,'9999.9'))
C
C    ==================================================================
C                            READING VELOCITY DATA
C    ==================================================================
C
  100	indvel=0
  	open(unit=99,file='../info/veloc.lis',status='old',err=200)
C
C ....opening the output file and designing the file headlines (only in 
C     case this operation was not performed before due to the absence 
C     of CTD profiles).
      if(indctd.eq.0) then
	   read(99,12)
   	   read(99,38,end=180) na1,na2
         rewind(99)
C 
	   open(unit=2,file='../level/'//na1//char(lev/1000+48)//
     .    char(lev/100-10*(lev/1000)+48)//char(lev/10-10*(lev/100)+48)//
     .    char(lev-10*(lev/10)+48)//'.lev',status='unknown')
	   write(2,20) NA1,STUDA
	   write(2,25) p00,pref
      endif
C 
      indvel=1
C    ------------------------------------------------------------------
C
C ....start reading
	read(99,12)
  135	read(99,38,end=180) na1,na2

  136	format(1x,f6.1,2x,a4,a3,'.vel')
	write(6,136) p00,na1,na2
	open(unit=1,file='../vel/'//na1//na2//'.vel',status='old')    
      read(1,40) latg0,latm0,lats0,long0,lonm0,lons0
C	write(6,40) latg0,latm0,lats0,long0,lonm0,lons0
    	xlats0=0
	xlat0=float(abs(latg0))+float(latm0)/60.+xlats0/3600.
	if(latg0.lt.0.) xlat0=-xlat0
	xlons0=0
	xlon0=float(abs(long0))+float(lonm0)/60.+xlons0/3600.
	if(long0.lt.0) xlon0=-xlon0
	idep=0
C	write(6,89) xlat0,xlon0  
C
C	initial values are read...
  144	read(1,188,end=175) p2,u2,v2
	if((u2.gt.9999.).or.(v2.gt.9999.)) goto 144
C	warning: sometimes there are no data just below the surface
	if(p2.gt.p00) goto 175
  145	p1=p2
	u1=u2
	v1=v2
C	when the chosen level is deeper than the profile, the station  
C	will simply not be included in the data set
  146	read(1,188,end=175) p2,u2,v2
  	if((u2.gt.9999).or.(v2.gt.9999)) goto 146
	if(p2.lt.p00) goto 145
	u0=u1+(u2-u1)*(p00-p1)/(p2-p1)
	v0=v1+(v2-v1)*(p00-p1)/(p2-p1)
C
    	write(2,190) na1,na2,xlat0,xlon0,u0,v0
  175	close(unit=1)
	goto 135
  180	close(unit=99)
  188 format(1x,f6.1,2f7.1)
  190	format(1x,a4,'v',a3,1x,f7.3,1x,f8.3,3(1x,'99.999'),1x,'999.999',
     .                                                 2(1x,f6.1))
C    ==================================================================
C
	close(unit=2)
  200	if((indctd.eq.0).and.(indvel.eq.0)) then
         write(6,205)
  205	   format(//,1x,'!!! NEITHER CTD NOR VELOCITY PROFILE NAMES HAVE'
     .   ' BEEN FOUND !!!'/)
         goto 220
	endif        
  210 continue
  220	end


C  ====================================================================
      REAL FUNCTION SVAN(S,T,P0,SIGMA)
C
C  (TESTED IN OCT 92 GIVING EXACT RESULTS ACCORDING TO UNESCO TABLES)
C

C
C  SPECIFIC VOLUME ANOMALY (STERIC ANOMALY) BASED ON 1980 EQUATION
C  OF STATE FOR SEAWATER AND 1978 PRACTICAL SALINITY SCALE.
C  REFERENCES:
C  MILLERO ET AL. (1980) DEEP-SEA RES., 27A, 255-264
C  MILLERO AND POISSON 1981, DEEP-SEA RES. 28A PP 625-629
C  BOTH ABOVE REFERENCES ARE ALSO FOUND IN UNESCO REPORT N0. 38 (1981)
C  UNITS:
C       PRESSURE        P0      DECIBARS
C       TEMPERATURE     T       DEG CELSIUS (IPTS-68)
C       SALINITY        S       (PSS-68)
C       SPEC. VOL. ANO. SVAN    1.0E-8 M**3/KG
C       DENSITY ANO.    SIGMA   KG/M**3
C
C CHECK VALUE: SVAN=981.30210E-8  M**3/KG FOR S=40 (PSS-78),
C              T=40 DEG C, P0=10000 DECIBARS
C CHECK VALUE: SIGMA=59.82037 KG/M**3     FOR S=40 (PSS-78),
C              T=40 DEG C, P0=10000 DECIBARS
C
      REAL P,T,S,SIG,SR,R1,R2,R3,R4
      REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW,K35

C  EQUIV
      EQUIVALENCE (E,D,B1) ,(BW,B,R3),(C,A1,R2)
      EQUIVALENCE (AW,A,R1),(KW,K0,K)

C  DATA
      DATA R3500,R4 /1028.1063, 4.8314E-4/
      DATA DR350 /28.106331/
C  R4 IS REFERRED TO AS C IN MILLERO AND POISSON 1981

C  CONVERT PRESSURE TO BARS AND TAKE SQUARE ROOT SALINITY.
      P=P0/10.
      SR=SQRT(ABS(S))
C
C  PURE WATER DENSITY AT ATMOSPHERIC PRESSURE
C    BIGG P.H. ,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
      R1= ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T
     &    -9.095290E-3)*T + 6.793952E-2)*T-28.263737

C  SEAWATER DENSITY ATM PRESS.
C  COEFFICIENTS INVOLVING SALINITY
C  R2=A  IN NOTATION OF MILLERO AND POISSON 1981
      R2= (((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T
     &   +8.24493E-1

C  R3=B INT NOTATION OF MILLERO AND POISSON 1981
      R3=(-1.6546E-6*T +1.0227E-4)*T -5.72466E-3

C  INTERNATIONAL ONE-ATMOSPHERE EQUATION OF STATE OF SEAWATER
      SIG=(R4*S + R3*SR + R2)*S + R1

C  SPECIFIC VOLUME AT ATMOSPHERIC PRESSURE
      V350P=1.0/R3500
      SVA= -SIG*V350P/(R3500+SIG)
      SIGMA=SIG+DR350

C  SCALE SPECIFIC VOL. ANOMALY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8
      IF (P.EQ.0.0) RETURN

C
C       NEW HIGH PRESSURE EQUATION OF STATE FOR SEAWATER
C
C       MILLERO ET AT, 1980 DSR 27A, PP 255-264
C                CONSTANT NOTATION FOLLOWS ARTICLE
C

C  COMPUTE PRESSURE TERMS
      E =(9.1697E-10*T+2.0816E-8)*T-9.9348E-7
      BW=(5.2787E-8*T -6.12293E-6)*T+3.47718E-5
      B= BW + E*S
C
      D =1.91075E-4
      C =(-1.6078E-6*T-1.0981E-5)*T+2.2838E-3
      AW= ((-5.77905E-7*T+1.16092E-4)*T+1.43713E-3)*T
     &    -0.1194975
      A = (D*SR +C)*S + AW

C    
      B1= (-5.3009E-4*T+1.6483E-2)*T +7.944E-2
      A1= ((-6.1670E-5*T+1.09987E-2)*T -0.603459)*T + 54.6746
      KW= (((-5.155288E-5*T+1.360477E-2)*T-2.327105)*T
     &    + 148.4206)*T -1930.06
      K0= (B1*SR + A1)*S + KW

C  EVALUATE PRESSURE POLYNOMIAL
C
C  K EQUALS THE SECANT BULK MODULUS OF SEAWATER
C  DK=K(S,T,P)-K(35,0,P)
C  K35=K(35,0,P)
C
      DK= (B*P + A)*P + K0
      K35= (5.03217E-5*P + 3.359406)*P + 21582.27
      GAM=P/K35
      PK = 1.0 - GAM
      SVA= SVA*PK + (V350P+SVA)*P*DK/(K35*(K35+DK))

C  SCALE SPECIFIC VOL. ANOMALY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8
      V350P=V350P*PK
C
C  COMPUTE DENSITY ANOMALY WITH RESPECT TO 1000.0 KG/M**3
C  1) DR350: DENSITY ANOMALY AT 35 (PSS-78), 0 DEG. C AND 0 DECIBARS
C  2) DR35P: DENSITY ANOMALY 35 (PSS-78), 0 DEG. C, PRES. VARIATION
C  3) DVAN : DENSITY ANOMALY VARIATIONS INVOLVING SPECIFIC VOL. ANOMALY
C
C  CHECK VALUE: SIGMA = 59.82037 KM/M**3 FOR S=40 (PSS-78),
C  T=40 DEG C, P0=10000 DECIBARS
C
      DR35P=GAM/V350P
      DVAN=SVA/(V350P*(V350P+SVA))
      SIGMA=DR350+DR35P-DVAN

      RETURN
      end
C  =====================================================================
