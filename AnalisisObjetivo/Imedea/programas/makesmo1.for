c	=================================================================
C               UNSTABLE LAYERS SEARCHER - *SMO FILES GENERATOR
C
C                             Dami� Gomis
C
C      IMEDEA (Institut Mediterrani d'Estudis Avan�ats), a joint 
C      centre between the Universitat de les Illes Balears (UIB) and 
C      the Spanish Research Council (CSIC).
C
C      This code has been compiled and distributed in the framework 
C      of the project REN2000-2599-E, funded by the Spanish Marine 
C      Science and Technology subprogram.
C 
C	Notes about this code:
C
C	This program searchs for unstable layers in the original 
C	vertical profiles at stations and removes them; the way it  does
C	it is simply by making the density profile of the unstable layer
C	to vary linearly from the top to the bottom of the layer.
C
C	The density is computed using the state equation (1980)
C	=================================================================
	parameter(npoin=1000,tolkm=2.)
	real nst(npoin),lzmax
	dimension p(npoin),t(npoin),s(npoin),st(npoin)
	character*8,na
	character*12,name1,name2
	
c	-----------------------------------------------------------------
	rd=3.14159/180.
	tolat=tolkm/(60.*1.852)
c
	write(6,21)
   21 FORMAT(1X,' The program detects jumps on pressure variable.',/
     ./' >>> Introduce the maximum distance (meters) allowed between',/
     ./' consecutive pressure values:')

	read (5,*) lzmax

	open(unit=99,file='../info/ctd.lis',status='old')
	read(99,25)
	write(6,25)
   25 format(//)
c 
   30 read(99,35,end=80) na
   35 format(a8)
	
   	name1=na//'.prs'
	name2=na//'.smo'  
	open(unit=1,file='../prs/'//name1,status='unknown')
	write(*,36) name1
   36 format(2x,'Reading PRS file:',a12)
	open(unit=2,file='../smo/'//name2,status='unknown')
	write(*,37) name2
   37 format(2x,'Creating SMO file:',a12)    
	read(1,115) latg,latm,lats,long,lonm,lons
	write(2,115) latg,latm,lats,long,lonm,lons
c	-------------------------------------------------------------------
C	....... no problem for the moment
		read(1,117)
   43		read(1,*,end=65) p2,t2,s2
C	Case of absent data
		if((t2.gt.99.).or.(s2.gt.99.)) then
		st2=99.99
		goto 48
		endif
          call sigmat(t2,s2,st2)
		if(p2.le.0.) goto 43
		write(2,123) p2,t2,s2,st2
   45		p1=p2
		t1=t2
		s1=s2
		st1=st2
   46		read(1,*,end=65) p2,t2,s2
C		Case of absent data		
		if((t2.gt.99.).or.(s2.gt.99.)) then
		st2=99.99
		write(2,123) p2,t2,s2,st2
		goto 46
		endif
		if(p2-p1.gt.lzmax) write(6,47) p1,p2
   47		format(1x,'!!! Pressure jump:  from ',f6.1,' to ',f6.1)
    		if(p2.gt.float(npoin)) goto 65
		call sigmat(t2,s2,st2)
		if(st2.le.st1) goto 50
   48		write(2,123) p2,t2,s2,st2
		if((t2.gt.99.).or.(s2.gt.99.)) goto 43
          goto 45
c	-------------------------------------------------------------------
c	....... an unstable layer has been detected

   50		ptop=p1
		sttop=st1
		p(1)=p2
		t(1)=t2
		s(1)=s2
		st(1)=st2
			do 53 i=2,npoin
			read(1,*,end=60) p2,t2,s2
			if((t2.gt.99.).or.(s2.gt.99.)) goto 53
			if(p2-p(i-1).gt.lzmax) write(6,47) p1,p2
			if(p2.gt.float(npoin)) goto 60
			call sigmat(t2,s2,st2)
			if(st2.gt.sttop) goto 55
			p(i)=p2
			t(i)=t2
			s(i)=s2
			st(i)=st2
   53		continue
   55		do 57 j=1,(i-1)
			nst(j)=sttop+(st2-sttop)*(p(j)-ptop)/(p2-ptop)
			write(2,125) p(j),t(j),s(j),st(j),nst(j)
   57		continue
		write(2,123) p2,t2,s2,st2
	    goto 45
c	-------------------------------------------------------------------
c	....... end of file detected before the end of the unstable layer
c		the density is set constant in that layer
   60		do 63 j=1,(i-1)
			write(2,125) p(j),t(j),s(j),st(j),sttop
   63		continue
c	-------------------------------------------------------------------
   65		close(unit=1)
		close(unit=2)
	goto 30
   80	close(unit=99)
c
  115	format(8x,i3,1x,i2,1x,f5.2,8x,i3,1x,i2,1x,f5.2)
  117	format(//)
  119	format(2x,'LATITUD',2x,'LONGITUD',2x,'B.DEPTH'/
     .  2x,f7.3,2x,f8.3,2x,f7.1/)
  123	format(1x,f6.1,3f7.3)
  125	format(1x,f6.1,3f7.3,' *',f8.4)
	
	end
c
      subroutine sigmat(t,s,st)
C
C     COMPUTES SIGMA-T ACCORDING TO THE INT. EQ. STATE 1980
C     T= TEMPERATURE IN DEG CC   
C     S= SALINITY
C     ST= SIGMA-T
C
      RW=999.842594+((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)
     ?*T-9.095290E-3)*T+0.06793952)*T
      BT=0.824493+(((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-0.0040899)*T
      CT=-5.72466E-3+(-1.6546E-6*T+1.0227E-4)*T
      SRS=SQRT(S)
      ST=(RW+BT*S+CT*S*SRS+4.8314E-4*S*S)-1000.
      return
      end
