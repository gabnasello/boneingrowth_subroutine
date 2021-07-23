      MODULE GLOBVAR
C
C	 -----------STATE VARIABLES--------------
C
C	----GRANULATION----
C
C	STATEV(1) is the mineralized tissue density
C	of the granulation tissue (-rhotissue-)
C	[g/cm³]
C	STATEV(2) is the Young modulus of the
C	mineralized granulation tissue
C	[Pa]
C	STATEV(3) is the bone deposition rate
C	(-bonerate-) in the UMAT subr
C	[mm³/(mm³*day)]
C	STATEV(4) -cbone- is the normalized cell concentration
C	[adim.]
C
C	----SCAFFOLD----
C
C   STATEV(6) -Dcum- is the cumulative damage
C	[% degradation]
C	STATEV(7) is the Young modulus of the scaffold
C	[Pa]
C   STATEV(8) -Dstate- is the  damage state of the element
C	(0: passive, 1:active, 2:failed)
C	[adim.]
C   STATEV(9) -pit- is the pitting parameter of the element
C	[adim.]
C	STATEV(10...13) contain the connectivity map - i.e the labels of
C	the neighboring elements
C	[adim.]
C
C	 -----------PARAMETERS UPDATED INSIDE UMAT_SCAFFOLD--------------
C	 [system unit mechanical properties - mm, g, s]
C	 [system unit for mechano-driven bone ingrowth - mm, g, day]
C
C	 -Dthresh- is the damage threshold above which the scaffold is replaced
C  by newly mineralized tissue
	 parameter(Dthresh=0.99)
C     -number of elements in the scaffold mesh (used for the definition
C	 of the arrays)
	 parameter(NrElements=100000)
C
C	 -----------PARAMETERS UPDATED INSIDE UMAT_GRANULATION--------------
C	 [system unit mechanical properties - mm, g, s]
C	 [system unit for mechano-driven bone ingrowth - mm, g, day]
C
C	 ---- Strain-based mechanical stimulus ---
C
C	 -nstimul- is the average num. of cycles per timeunit
	 parameter(nstimul=10000)
C	 [cycles/day] (Mohaghegh et al., Comput. Methods Appl. Mech. Engrg.,2014)
C	 -m- an experimental exponent
	 parameter(mpar=4)
C	 [adim.] (Beaupré et al, J of Orth Res, 1990]
C
C	 -refstim- is the reference strain stimulus in the peri-SCAFFOLD region
	 parameter(refstim = 2083.41)
C	 [ustrain]
C
C	 -alpha- is the fraction of reference daily stimulus
	 parameter(alpha = 0.60)
C	 [adim.]
C
C	 --bone deposition rate--
C
C	 -boneratemax- is the maximum bone volume deposition rate
	 parameter(boneratemax=0.04)
C	 [mm³/(mm³*day)]  (adapted from Adachi et al., Phil. trans., 2010)
C
C	 -expemp- bone ingrowth per unitary microstrain stimulus at max cell density
	 parameter(expemp=7.00e-07)
C	 [mm³/(mm³*day)]
C
      END MODULE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	UMAT SUBROUT. IS EXECUTED BEFORE UMATHT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	IF (CMNAME(1:7) .EQ. 'SCAFFOLD') THEN
	 CALL UMAT_SCAFFOLD(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	ELSE IF(CMNAME(1:11) .EQ. 'GRANULATION') THEN
	 CALL UMAT_GRANULATION(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	END IF

     	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	SUBROUTINE UMAT_SCAFFOLD(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C	 -youngscaff- is the scaffold young modulus
C	 [Pa]
C	 -egscaff- is the scaffold shear modulus
C	 [Pa]
C	 -elamscaff- is the scaffold Lamé's first parameter
C	 [Pa]
C	 -Dcum- is the cumulative scaffold degradation
C	 [% degradation]
	   real*8 Dcum, youngscaff, egscaff, elamscaff
C     -Dstate contains the values of the damage state at the start of the increment
      Integer Dstate(NrElements)
C     -pit contains the values of the pitting parameter at the start of the increment      
      real*8 pit(NrElements)
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY
C ----------------------------------------------------------------
C
C
C ----------------------------------------------------------------
C ELASTIC PROPERTIES
C
C	 -younginit- scaffold initial Young modulus
C	 [Pa]
C    -scaffnu- scaffold Possion ratio
C    [adim.]
C	 -degradscaff- mechanical degradation constant
C	 [ %degradation / day ]
C
	scaffnu       =PROPS(12)
	degradscaff   =PROPS(13) 
	younginit     =PROPS(11) 
C
C   -beta is a parameter of the unoiform degradation equation [adim.]
	beta = PROPS(14)
C
C   -deltaU is a characteristic dimension of the corrosion process
C	(Gastaldi et al., 2011) [mm]
	deltaU = PROPS(15)
C      
C   -Dcum- is the cumulative damage constant [% degradation]
	Dcum =  STATEV(6)
C   -Dstate- is the  damage state of the element (0: passive, 1:active, 2:failed)
	Dstate(NOEL) = STATEV(8)
C   -pit- is the pitting parameter (=1 for uniform corrosion)
	pit(NOEL) = STATEV(9)
C
C	NOEL is the label number of the current element
C
C   Dstate = 2: the element is already degraded and replaced by granulation tissue
C
	IF (STATEV(8) .EQ. 2) THEN
	 CALL UMAT_GRANULATION(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
	ELSE
C
C     Dstate = 1: the element is being corroded
C
	IF (STATEV(8) .EQ. 1) THEN
	  Dcum = Dcum + (STATEV(9)*degradcst*du)/CELENT * DTIME
C	CELENT is the characteristic length of the element, it is provided by abaqus
	  IF (Dcum .GE. Dthresh) THEN
C
        STATEV(6) = Dthresh
		STATEV(8) = 2
		CALL UMAT_GRANULATION(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
		
C
C
	  ELSEIF(Dcum .LT. Dthresh) THEN
		STATEV(6) = Dcum
		youngscaff = younginit * (1 - Dcum)
C
C   if one of the nighboring elements is degraded, a fraction of its pitting
C   parameter is transmitted to the current element
        DO i=10,13
          IF (Dstate(STATEV(i)) .EQ. 2) THEN
            STATEV(9) = beta * pit(STATEV(i))
          ENDIF
        ENDDO
	  ENDIF
	ENDIF
C
C     Dstate = 0: the element is not corroded
C
	  IF (STATEV(8) .EQ. 0) THEN
		youngscaff = younginit
C   if one of the nighboring elements is degraded, the current element becomes active
C   and a fraction of its pitting parameter is transmitted to the current element
		DO i=10,13
          IF (Dstate(STATEV(i)) .EQ. 2) THEN
            STATEV(8) = 1
            STATEV(9) = beta * pit(STATEV(i))
          ENDIF
        ENDDO
      ENDIF
C

    !   STATEV(6) = Dcum
    !   youngscaff = younginit * (1 - Dcum)
    !   STATEV(7) = youngscaff
      egscaff = youngscaff/2.0/(1.0 + scaffnu)
	  elamscaff = youngscaff * scaffnu/(1.0 + scaffnu)/(1.0-2.0*scaffnu)
C
C ELASTIC STIFFNESS
C
	DO K1=1,NDI
	  DO K2=1,NDI
	    DDSDDE(K2,K1)=elamscaff
	  ENDDO
	  DDSDDE(K1,K1)=egscaff*2.0+elamscaff
	ENDDO
C
	DO K1=NDI+1,NTENS
	  DDSDDE(K1,K1)=egscaff
	ENDDO
C
C CALCULATE STRESS
C
	DO I1=1,NTENS
	  STRESS(I1)=0.0
	ENDDO
	DO I1=1,NTENS
	  DO I2=1,NTENS
	    STRESS(I1)=STRESS(I1)+DDSDDE(I1,I2)*(STRAN(I2)+DSTRAN(I2))
	  ENDDO
	ENDDO
C
C
	ENDIF
C
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	SUBROUTINE UMAT_GRANULATION(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C	 -youngmod- is the material young modulus
C	 [Pa]
C	 -nu- is the material Poisson's coefficient
C	 [adim.]
C	 -eg- is the material shear modulus
C	 [Pa]
C	 -elam- is the Lamé's first parameter
C	 [Pa]
C	 -tissuestrain- is the effective tissue microstrains
C	 [microstrain]
C	 -mechstimul- is the daily mechanical strain stimulus for bone formation
C	 [microstrain/day]
C	 -localstimul- is the daily mechanical strain stimulus fin the peri-SCAFFOLD region
C	 [microstrain/day]
C	 -bonerate- is the bone volume deposition rate
C	 [mm³/(mm³*day)]
C	 -densityrate- is the bone density rate
C	 [g/cm³/day]
C	 -ENERGY- is the actual strain energy density
C	 [Pa]
C	 -rhotissue- is the mineralized tissue density
C	 map for the granulation tissue
C	 [g/cm³]
      real*8 youngmod, nu, eg, elam, tissuestrain, mechstimul,
     1localstimul, bonerate, densityrate, ENERGY, rhotissue,
     2interceptinterm, ingrowthstrain
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY
C ----------------------------------------------------------------
C
C
C ----------------------------------------------------------------
C ELASTIC PROPERTIES
C
C	 ----Newly formed bone properties---
C
C Material model adapted from Ghouse et al, Appl Mat Today, 2019
C
C	 -rhotrabecular- is the upper limit of
C	 the bone mineral density for trabecular bone
C	 [g/cm³]
C	 -rhocortical- is the lower limit of
C	 the bone mineral density for cortical bone
C	 [g/cm³]
C	 -rhomax- is the maximum bone mineral density
C	 [g/cm³]
C
      rhotrabecular     =PROPS(1)
      rhocortical		=PROPS(2)
      rhomax       	    =PROPS(3)
C
C	 -nu- is the bone Possion ratio
C  [adim.]
      nu       	        =PROPS(4)
C
C	 --trabecular region--
C
C	 -consttrabecular- is the power law constant
C	 in the trabecular range of bone mechanical properties
C  [Pa]
C	 -exptrabecular- is the power law exponent
C	 in the trabecular range of bone mechanical properties
C  [adim.]
C
      consttrabecular	    =PROPS(5)
      exptrabecular      	=PROPS(6)
C
C	 --intermediate region--
C
C	 -slopeinterm- is the linear slope
C	 in the intermediate range of bone mechanical properties
C  [Pa]
C	 -interceptinterm- is the intercept coefficient
C	 in the cortical range of bone mechanical properties
C  [Pa]
C
      slopeinterm        	=PROPS(7)
      interceptinterm    	=PROPS(8)
C
C	 --cortical region--
C
C	 -constcortical- is the power law constant
C	 in the cortical range of bone mechanical properties
C  [Pa]
C	 -exptrabecular- is the power law exponent
C	 in the cortical range of bone mechanical properties
C  [adim.]
C
      constcortical        =PROPS(9)
      expcortical       	 =PROPS(10)
C
C --MODEL--
C
C	STATEV(1) is the mineralized tissue density
C	map for the granulation tissue (-rhotissue-)
	rhotissue = STATEV(1)
C	[g/cm³]
C
C 	UPDATE MECHANICAL PROPERTIES
	IF (rhotissue .LE. rhotrabecular) THEN
	 youngmod = consttrabecular * rhotissue ** exptrabecular
C	 [Pa]
	ELSEIF(rhotissue .GT. rhotrabecular .AND.
     1 rhotissue .LT. rhocortical) THEN
	 youngmod = slopeinterm * rhotissue - interceptinterm
C	 [Pa]
	ELSEIF(rhotissue .GE. rhocortical) THEN
         youngmod = constcortical * rhotissue ** expcortical
C	 [Pa]
	ENDIF
	STATEV(2) = youngmod
C
	eg=youngmod/2.0/(1.0+nu)
	elam=youngmod*nu/(1.0+nu)/(1.0-2.0*nu)
C
C ELASTIC STIFFNESS
C
	DO K1=1,NDI
	 DO K2=1,NDI
	  DDSDDE(K2,K1)=elam
	 ENDDO
	 DDSDDE(K1,K1)=eg*2.0+elam
	ENDDO
C
	DO K1=NDI+1,NTENS
	 DDSDDE(K1,K1)=eg
	ENDDO
C
C CALCULATE STRESS
C
	DO I1=1,NTENS
	 STRESS(I1)=0.0
	ENDDO
	DO I1=1,NTENS
	 DO I2=1,NTENS
	  STRESS(I1)=STRESS(I1)+DDSDDE(I1,I2)*(STRAN(I2)+DSTRAN(I2))
	 ENDDO
	ENDDO
C
C CALCULATE THE STRAIN ENERGY DENSITY AND THE EFFECTIVE STRAIN
	SSE=0.0
	DO I=1,NTENS
	 SSE=SSE+0.5*STRESS(I)*(STRAN(I)+DSTRAN(I))
	ENDDO
C
C	effective tissue microstrains [Pistoia et al. 2002, Bone]
	tissuestrain=SQRT(2*SSE/youngmod)*(1E+06)
	STATEV(5) = tissuestrain
C	 [microstrain]
C
C	tissuestrain reduced by the reference stimulus
	ingrowthstrain = tissuestrain - alpha*refstim
	IF (ingrowthstrain .LT. 0.0) THEN
	 ingrowthstrain = 0.0
	ENDIF
C
C	 -mechstimul- is the effective stress at the tissue level
	mechstimul = (nstimul*ingrowthstrain**mpar)**(1.0/mpar)
C	 [microstrain/day]
C
C	 -bonerate- is the bone deposition rate at max cell density
	bonerate = expemp * mechstimul
C	 [mm³/(mm³*day)]
C
C 	LIMIT BONE DENSITY
       IF (bonerate .GE. boneratemax)THEN
	 bonerate = boneratemax
       ENDIF
	STATEV(3) = bonerate
C
C	 -cbone- is the normalized cell concentration reduced by
C	the threshold val. for bone ingrowth
	cbone = STATEV(4)
C	 [adim.]
C
C	 -densityrate- is the bone density rate
	 densityrate = cbone * bonerate * rhomax
C	 [g/cm³/day]
C
C 	INFERIOR LIMIT DENSITY RATE
       IF (densityrate .LT. 1E-08)THEN
        densityrate = 0
       ENDIF
C
C
C	 -rhotissue- is the mineralized tissue density
	 rhotissue = rhotissue + densityrate * DTIME
C	 [g/cm³]
C
C 	LIMIT BONE DENSITY
       IF (rhotissue .GE. rhomax)THEN
	 rhotissue = rhomax
       ENDIF
	STATEV(1) = rhotissue
C
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
C
	IF (CMNAME(1:7) .EQ. 'SCAFFOLD') THEN
	 CALL UMATHT_SCAFFOLD(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	ELSE IF(CMNAME(1:11) .EQ. 'GRANULATION') THEN
	 CALL UMATHT_GRANULATION(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	END IF
C
     	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATHT_SCAFFOLD(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
C
C	 -Dcum- is the cumulative scaffold degradation
C	 [% degradation]
      real*8 Dcum
C
      Dcum = STATEV(6)
C
	IF (Dcum .GE. Dthresh) THEN
	 CALL UMATHT_GRANULATION(U,DUDT,DUDG,FLUX,DFDT,DFDG,
  	1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
  	2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
  	3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
	ELSEIF(Dcum .LT. Dthresh) THEN
C
C ----------------------------------------------------------------
C UMATHT FOR UNCOUPLED HEAT TRANSFER
C No cell diffusion inside the scaffold if it is not fully degraded.
C ----------------------------------------------------------------
C
C		input specific heat
	DUDT = 0
	DU = 0
	U = 0
C
C		input flux = -[k]*{dtemdx}
C		input isotropic conductivity
	DO I=1, NTGRD
	 FLUX(I) = -0
	 DFDG(I,I) = -0
	END DO
C	
	END IF
C
	RETURN
	END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATHT_GRANULATION(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C	tocall the globvar modul
	USE GLOBVAR
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
C
C	 -cbone- is the normalized cell concentration of the new
C	 increment
C	 [adim.]
      real*8 cbone
c
C ----------------------------------------------------------------
C UMATHT FOR UNCOUPLED HEAT TRANSFER
C ----------------------------------------------------------------
C
	SPECHT = 1
C
C	 -COND- is the cell diffusion constant (from Geris et al., J. Theor. Biol., 2008)
C  [mm²/day]
	COND = PROPS(1)
C
C		input specific heat
	DUDT = SPECHT
	DU = DUDT*DTEMP
	U = U+DU
C
C		input flux = -[k]*{dtemdx}
C		input isotropic conductivity
	DO I=1, NTGRD
	 FLUX(I) = -COND*DTEMDX(I)
	 DFDG(I,I) = -COND
	END DO
C
	cbone = (TEMP+DTEMP)
	STATEV(4) = cbone
	RETURN
	END