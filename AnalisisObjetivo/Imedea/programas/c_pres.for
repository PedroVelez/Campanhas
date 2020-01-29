C	=================================================================
C	   THIS PROGRAM CONSTRUCTS  PRESSURE-LEVEL  DATA FILES (*.LEV)
C	      FROM THE SMOOTHED STATION-PROFILE DATA FILES (.SMO).
C
C					Damià Gomis, Version June 2002
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
C	Notes about this code:
C
C	Variables are: TEMP, SAL, DENS  and  DYNH=INTEGRAL[-Van dp],
C                                    'Van' being the Volum Anomaly 
C
C	DYNH is normally computed with respect to a fixed reference  
C     level (option 1 in this code). In this case, if a given profile 
C     is shallower than the selected reference level, DYNH cannot be 
C     computed at that station.
C
C     Another option is to compute DYNH using the level immediately 
C     below (as defined in the file "grid.dat") as a reference level 
C     (option 2). In that case, even if a given profile is shallower 
C     than the vertical domain defined in "grid.dat", DYNH can still 
C     be computed at several levels. On the other hand, this option
C     implies that, once DYHN has been interpolated onto a grid, to
C     retrieve DYNH at a particular level the values of DYNH at all 
C     the levels below must be added. See the help document, for 
C     further details and about the advantadges of this potion.
C
C
C		    STATION DATA FILES ARE READ FROM "CTD.LIS"
C			      
C	=================================================================
	character*4,NA1  
	character*1,NA2
	character*3,NA3
	character*13,name2,name
	character*4,star
	character*8,STUDA 
	real ALF0,lats0,lons0
C	-----------------------------------------------------------------
C .....	input data block
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
	read(99,2) NA1,STUDA 
    2 format(a4,1x,a8)
    3 format(1x/1x)
    
      close(unit=99)
	
	write(6,4)
    4 format(2X,'THIS PROGRAM CONSTRUCTS PRESSURE-LEVEL FILES *.LEV'/
     .2X,'FROM THE SMOOTHED STATION PROFILE FILES *.SMO' /)
	write(6,5)
    5 format(/'Select option:'//
     .         5x,'>>>Individual level (1)'/
     .         5x,'>>>All levels (2)')                                    
	read(5,*) LE
	if(le.eq.1) then
	   write(6,'(a,$)') 'Input level (meters)? : '
	   read(5,*) p0
	   NL=1
	endif
	write(6,6)
    6 format(/'Select option:'// 
     .         5X,'>>>Fixed reference level (1)'/
     .         5X,'>>>Level just below (see grid.dat file) (2)')
	read(5,*) RE_RLEV
	if(re_rlev.eq.1) then                                
	    write(6,'(a,$)') 'Depth of reference level ? : '
	read(5,*) PREF
	endif

    7     do 150 k=1,NL

      p00=p0+pint*float(k-1)
	lev=ifix(p00)          
	if(re_rlev.eq.2) then
	pref=p00+pINT
	endif
C ....start reading
	open(unit=99,file='../info/ctd.lis',status='old')
	read(99,12)
   	read(99,37,end=80) NA1,NA2,NA3
      rewind(99)
   12 format(//)        
C	-----------------------------------------------------------------
C ....designing the output-file headlines
	open(unit=2,file='../level/'//na1//char(lev/1000+48)//
     .char(lev/100-10*(lev/1000)+48)//
     .char(lev/10-10*(lev/100)+48)//
     .char(lev-10*(lev/10)+48)//'.lev',status='unknown')
	write(2,20) NA1,STUDA
	write(2,25) p00,pref
   20 format(1x,'***     DATA ON ISOBARIC SURFACES READY TO BE ',
     .'PROCESSED     ***'/1x,'***',5x,'CRUISE CODE: -',a4,
     .  '- .   Region:',a8,5x,'***'/)
   25	format(24x,'level:',f7.2/16x,'pr-ref. level:',
     .  f7.2//22x,'<',12x,'DATA TO BE USED',13x,'>       <REMARKS>'
     .  /'---------',(1x,'-------'),' --------',3(1x,'------'),1x,
     .'-------',2(1x,'------') ,'  ---------'/1x,'  A8 ',2x,(3x,'F7.3')
     .,2x,(3x,'f8.3'),1x,3(3x,'f6.3'),'  dyn. cm',2x,'cm/s',3x,
     .'cm/s',/1x, 'SYMBOL   LATITU   LONGIT  TEMPER SALINI DENSIT DYN.HE
     .I   U',6x,'V'/'---------',(1x,'-------'),' --------',
     . 3(1x,'------'),1x,'-------',2(1x,'------'), '  ---------')
C	------------------------------------------------------------------
C
C ....start reading
	read(99,12)
   35 read(99,37,end=80) NA1,NA2,NA3
   37 format(a4,a1,a3)
	write(6,38) p00,NA1,NA3
   38 format(1x,f5.1,2x,a4,a3)
C
	open(unit=1,file='../smo/'//NA1//NA3//'.smo',status='old') 
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
C	write(6,119) xlat0,xlon0  
C
C	initial values are read...
   44 read(1,118,end=75) p2,t2,s2,st2,star,sst2
	if(star.eq.' * ') st2=sst2
	if((t2.gt.99.).or.(s2.gt.99.)) goto 44
C	warning: sometimes there are no data just below the surface
	if(p2.gt.p00) goto 75
   45	p1=p2
	t1=t2
	s1=s2
	st1=st2
C	when the chosen level is below the bottom, the station will 
C	simply not be included in the data set
   46	read(1,118,end=75) p2,t2,s2,st2,star,sst2
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
C	-------------------------------------------------------------------
C ....when the reference level is below the bottom, the geopotential
C	is set equal to: 
   60	geo0=999.999
C	------------------------------------------------------------------
   70	write(2,130) na1,na3,xlat0,xlon0,t0,s0,st0,geo0
   75	close(unit=1)
	goto 35
   80	close(unit=99)
	close(unit=2)
C
  118 format(1x,f6.1,3f7.3,a3,f7.3)
  119 format(2x,'LATITUD',2x,'LONGITUD',2x,f7.3,2x,f8.3/)
  120	format(1x,f6.1,3f7.3,a3,f8.5)
  125	format(1x,f6.1,2f7.3)
  130	format(1x,a4,a4,(1x,f7.3),(1x,f8.3),3(1X,F6.3),1X,F7.3,
     .                                  2(1X,'9999.9'),2X,I4)                  
C
  150      continue
C
  	end


C  =====================================================================
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
