
      PROGRAM PBwires_CU111
C     -----------------------------------------------------
C     ****** 3D MONTE CARLO - 3 Pb Nanowires ON  n-LAYER Cu(111) *******
C            Cu-ATOMS MOVABLE
C     ****** Calculation of Vacancies, Holes
C            
C     -----------------------------------------------------
C    
      CHARACTER*30 NAME       
      CHARACTER*17 FNAME,FOLD,COLD,COLDC
      CHARACTER*1  FM1,FM2,FM3,FM4,FN1,FN2,FN3,FN4,FNI
      CHARACTER*8  F1
      CHARACTER*7  F4 
      CHARACTER*2  F6
      CHARACTER*1  F2, F3, F5, F10 
      CHARACTER*5  F7
      CHARACTER(LEN=200) :: CMD
C      
C       SUBSTRAT LATTICE CELL - Cu(111)
C
      PARAMETER(AX = 1.0, AY = 1.7320508, DAY = 1/AY, AZ = 0.816496)

c in units of Cu lattice constant(cub fcc) ::::::ax = 1/sqrt(2);  
c in units of ax::: ay = sqrt(3); az = sqrt(2/3); day = 1/ay
      PARAMETER(AX2 = .5, AY2 = AY/2.,R = 1.987 / 23061.)                                                                    
C
C       POTENTIAL Cu - Cu
C PR B48,22,1993
      PARAMETER(UCU = 1.224, ACU = 0.0855/UCU)
      PARAMETER(RCU = AX, PCU = 10.96, QCU = 4.556)  
C	RCU diametar of Cu atom
C
C       POTENTIAL  Pb - Pb
C
C     PARAMETER(UPB=0.914,APB=0.098/UPB,RPB = 1.36957)  
      PARAMETER(UPB=0.914,APB=0.098/UPB,RPB = 1.05)
C	    RPB = 1.36957 - Initial Value // CHANGE to 1.05, 1.10, etc.

      PARAMETER(PPB = 9.576/RPB, QPB = 7.296/RPB)
C     PRB48,22,1993
C     RPB diameter of Pb atom
C
C       POTENTIAL  Pb - Cu
C
      PARAMETER(UPBCU=1.0577, APBCU=0.09175/UPBCU, RPBCU = (RCU+RPB)/2) !1.1847) !MISTAKE: 1.11847
      PARAMETER(PPBCU = 10.268/RPBCU, QPBCU = 5.926/RPBCU)
C   arithmetic average               
C      
C       SIMULATION BOX
C
      PARAMETER(MX = 76,MY = 44, NML = 2, NML1=NML+1, NZD=(NML+2)*1000)
      PARAMETER(AMX = MX * AX, AMY = MY * AY)
      PARAMETER(AMZPB = RPB/2.+ 25.5)				 !tavan
      PARAMETER(AMZ = NML*AZ+AMZPB, AMZ2 = 2.* AMZ, RAMZ =1./ AMZ)
      PARAMETER(MLE=2*MX*MY, NAT=MLE*NML, NAT1 = 1.5*MLE*NML1)
      PARAMETER(MLE2=2*MLE, NCUAT = NAT + MLE2)
      PARAMETER(NAT2=NAT1*2,NAT3=NAT1*3,NAT4=NAT1*4,NAT5=NAT1*5)
      PARAMETER(CUTOF = 3.0, NAT7 = NAT*7, NATP = NAT1 + MLE2)
      PARAMETER(CUTOFQ= CUTOF*CUTOF)
      PARAMETER(XXIN = 2./AMX,YYIN = 2./AMY,NMLP=4000 * NML1)  
C      
      real dddd1, dddd2, dddd3, V1, V2, Nvacsredno, Nholesredno
      integer Nvac1, Nvac2, Nvac3, Nhole1, Nhole2, Nhole3, TIMENVAC 

      REAL XP(NATP), ZP(NATP), YP(NATP)
      REAL VD(NAT1)
      REAL X(NAT1),Y(NAT1),Z(NAT1)            
      INTEGER*2 LIST(NAT1,NAT1),NRV(NAT1),NNUM(NAT1)
      INTEGER*2 NSTAT(NAT1)
      REAL X1(NAT7),Y1(NAT7),Z1(NAT7)   
      REAL DIX(NAT1), DIY(NAT1), DIZ(NAT1)
      REAL RRANF(NAT5) 
      REAL NMCSTP(1000),ENERPB(1000),ENERCU(1000),ENERAL(1000)  
      REAL DMAX(2),DMAX2(2),SKINQ(2),RPOV(2)
      INTEGER*4 SDF(NZD), PDF(NML1,4000)
      DATA DIX/NAT1*0./,DIY/NAT1*0./,DIZ/NAT1*0./ 
      DATA SDF/NZD*0/,PDF/NMLP*0/ 
      DATA DZZ/1000./ 
C
      INTEGER :: IHR, IMMM, ISEC, I100TH
      INTEGER :: values(8)

C 
      NRUNC = 0
      NRUNT = 0    
      FN1 = CHAR(79)
      FN2 = CHAR(111)
      FN3 = CHAR(77)
      FN4 = CHAR(109)

*******************************************************************
      write(*,*) ''
      write(*,*) ''
      write(*,*) '_________________________________________________'
      write(*,*) '_________________________________________________'
      write(*,*) ''
      write(*,*) ' 3 Pb nanowires on Cu(111)'
      write(*,*) '   each nanowire contains 50 Pb atoms            '
      write(*,*) '_________________________________________________'
      write(*,*) '_________________________________________________'
      write(*,*) ''
      write(*,*) ''
      write(*,*) ''
      write(*,*) 'Left/Right atoms are fixed!!'
      write(*,*) ''
      write(*,*) ''


      write(*,*) 'NAT1', NAT1
      WRITE(*,*)'FILE NUMBER '
      READ(*,121) II9,II1,II2
  121 FORMAT(I1,I1,I1)
      WRITE(*,*)'ONE OR MULTIPLE RUN: One - O OR Multiple - M'
      READ(*,122)FNI
  122 FORMAT(A1) 
      IF ((FNI.EQ.FN1).OR.(FNI.EQ.FN2)) NRUNI = 1
      IF ((FNI.EQ.FN3).OR.(FNI.EQ.FN4)) NRUNI = 2
      IF (NRUNI.GT.1) THEN
       WRITE(*,*)'COVERAGE OR TEMPERATURE STEP, C-1,T-2,CT-3,TC-4'
       READ(*,*)NCT
       IF (NCT.EQ.1) THEN   
        WRITE(*,*)'COVERAGE STEP'
        READ(*,*) DCOV  
        WRITE(*,*)'RUN NUMBER '
        READ(*,*)NRUNC       
       ENDIF
       IF (NCT.EQ.2) THEN 
        WRITE(*,*)'TEMPERATURE STEP'
        READ(*,*) DTEMP 
        WRITE(*,*)'RUN NUMBER '
        READ(*,*)NRUNT
       ENDIF
       IF ((NCT.EQ.3).OR.(NCT.EQ.4)) THEN
        WRITE(*,*)'COVERAGE STEP'
        READ(*,*) DCOV
        WRITE(*,*)'NUMBER OF COVERAGES'
        READ(*,*)NRUNC 
        WRITE(*,*)'TEMPERATURE STEP'
        READ(*,*) DTEMP
        WRITE(*,*)'NUMBER OF TEMPERATURES'
        READ(*,*)NRUNT
       ENDIF 
      ENDIF    
      IF (NRUNI.NE.1) THEN
       IF (NCT.EQ.1) NRUN = NRUNC
       IF (NCT.EQ.2) NRUN = NRUNT 
       IF (NCT.GT.2) NRUN = NRUNT * NRUNC
      ELSE
       NRUN = NRUNI
      ENDIF
      WRITE (*,*) '' 
      WRITE(6,*) 'ISEED = '
      READ(5,*) ISEED1
      ISEED = ISEED1
C ----------------------------
c	write(*,*) ''
c      WRITE(6,*)'number of Pb atoms (for the nanowire)'
c      READ(5,*) NADAT
c      write(*,*) ''
c ----------------------------
c	NE PIPAI
      NADAT = 150 !106 ! 53		! NE PIPAI
c	NE PIPAI
c-----------------------------


c	OPEN(111, file = 'time-Nvac123.dat')

      NCUSUB = NAT    
      WRITE(*,*) ''
      WRITE(6,*)'T = '
      READ(5,*) T1
      RT = R * T1
      WRITE(*,*)''
      WRITE(6,*) 'Number of MC steps = '
      READ(5,*) MCSMA1
      WRITE(6,*) 'Number of steps for equilibration = '
      READ(5,*) NINIT1
      WRITE(*,*)''
      WRITE(*,*)'INITIAL CONFIGURATION:'
      WRITE(*,*)'0 - Random distribution of Pb atoms'
      WRITE(*,*)'1 - (3) NANOWIRE(s) of Pb atoms'
      WRITE(*,*)'2 - DO NOT PRESS !!! '
      WRITE(*,*)''
      READ(*,*)NORDER

      write(*,*) 'NUMBER of MC steps to calculate VACANCIES and HOLES= '
      read(*,*) TIMENVAC

      WRITE(6,*) 'NUMBER of MC steps to be saved ENERGY= '
      READ(5,*) NWCOO1
      NWCOOD = NWCOO1 
      WRITE(6,*) 'NUMBER of MC steps to be saved SNAPSHOTS= '
      READ(5,*) NWSNA1
      NWCOOS = NWSNA1
C                           
      IBG = II9*100 + II1*10 + II2 + 1
      NKT = 1
      NKC = 1
      NRUNK = 0   
      LCT = 0
C      COV1 = 0                     
      COV = COV1
      T = T1    
      TA = T1
      DO IRUN = IBG, IBG + NRUN - 1 
        NRUNK = NRUNK + 1
        DO I = 1, NZD
          SDF(I) = 0
        END DO     
        DO J = 1, NML1   
          DO I = 1, 4000
            PDF(J,I) = 0
          END DO    
        END DO    
        II1 = IRUN - 1 
        IF (II1.GE.100) THEN
         II9 = II1/100
         II1 = II1 - II9*100
        ENDIF               
        II2 = II1 
        IF (II1.LE.9) II1 = 0
        IF (II2.GE.10) THEN
         II1 = II2 / 10
         II2 = II2 - II1*10
        ENDIF       
        FN3 = CHAR(II9 + 48)
        FN1 = CHAR(II1 + 48)
        FN2 = CHAR(II2 + 48)
        FNAME = 'C:\Pb-Cu\RESLG'//FN3//FN1//FN2
        WRITE(*,*)FNAME
C       RESULT = MAKEDIRQQ(FNAME)
C       Construct the command to create the directory
        CMD = 'mkdir ' // TRIM(FNAME)

C       Execute the system command to create the directory
        CALL SYSTEM(CMD)
        NAME = ' '
        ISEED = ISEED1 + 33557
        ISEED1 = ISEED
        IF (NCT.EQ.1) THEN 
         COV   = COV1 + (IRUN - IBG) * DCOV 
         NADAT = MLE * COV
         T = T1
        ENDIF
        IF (NCT.EQ.2) T   = T1  + (IRUN - IBG) * DTEMP 
        IF (NCT.EQ.3) THEN 
         IF (NKT.LE.NRUNT) THEN
          NKT = NKT + 1
          T   = T1  + (IRUN - IBG - LCT*NRUNT) * DTEMP
         ELSE
          LCT = LCT + 1
          T = T1
          COV = COV1 + LCT * DCOV 
          NADAT = MLE * COV
          NKT = 2                     
         ENDIF
        ENDIF
        IF (NCT.EQ.4) THEN
         IF (NKC.LE.NRUNC) THEN
          NKC = NKC + 1
          COV = COV1 + (IRUN - IBG - LCT*NRUNC) * DCOV 
          NADAT = MLE * COV
          T = TA       
         ELSE 
          LCT = LCT + 1
          TA  = T1  + LCT * DTEMP 
          COV = COV1
          NKC = 2
          T = TA   
         ENDIF
        ENDIF
        RT = R * T      
        NWCOOR = NWCOO1
        NWCOOD = NWCOOR 
        NWSNAA = NWSNA1
        NWSOOD = NWSNAA
        NINIT = NINIT1  
        MCSMAX = MCSMA1
C        CALL GETTIM(IHR,IMMM,ISEC,I100TH)

        CALL DATE_AND_TIME(values=values)
        IHR = values(5)    ! Hour
        IMMM = values(6)   ! Minute
        ISEC = values(7)   ! Second
        I100TH = 0         ! Hundredths of seconds (not directly provided by DATE_AND_TIME)

        WRITE(*,*)IHR,':',IMMM,':',ISEC
        MHR = IHR
        MMMM = IMMM
        MSEC = ISEC
        IHR1 = IHR
        MMM1 = IMMM
        MSEC1= ISEC 

      NATALL = NCUSUB + NADAT
 
      UCU2 = UCU*UCU
      UPB2 = UPB*UPB 
      UPC2 = UPBCU*UPBCU
      UACU = UCU*ACU
      UAPB = UPB*APB
      UAPC = UPBCU*APBCU
C
C       CALCULATION OF DMAX FROM THE MORSE POTENTIAL OF CU-CU
C  npar=3 and +9 AT coverage 0.5 
C 
           
      IF (COV.LE.0.4) THEN
       NPAR = 1
       NCU  = 11
      ENDIF
      IF (COV.EQ.0.5) THEN
       NPAR = 3
       NCU  = 9
      ENDIF
      IF (COV.GT.0.5) THEN
       NPAR = 6
       NCU  = 6
      ENDIF 
      DO I = 1,2
        NPAR = NPAR + NCU *(I-1)
        C    = 2.3333 * RT / NPAR 
        SQ   = SQRT(C)
        EA1  = 1 + SQ
        EA2  = 1 - SQ
        R1   = ALOG(EA1)/ 4.
        R2   = ALOG(EA2)/ 4.
        DMAX(I)  = ABS(R1 - R2)
        DMAX2(I) = 2 * DMAX(I)
        SKIN     = 10.*DMAX(I)
        SKINQ(I) = SKIN * SKIN * 0.25
        RM       = CUTOF + SKIN                    
        RMQ      = RM * RM   
      END DO                
C
C       COORDINATES OF HARD SUBSTRATE Cu(111) ATOMS
C
      RCUT  = CUTOF + AX2
      RCUT  = RCUT*RCUT
      N     = 1 
      NL = 0
      NL1 = NL + 1
      DO JZ = -8, 0, 1
        DO JY = -14, 14
          YJ = ABS(FLOAT(JY))
          DDELY = AMOD(YJ,2.)
          DO JX = -12, 12
            X1(N) = JX + DDELY * AX2 
            Y1(N) = JY * AY2 + NL*DAY            
            Z1(N) = JZ * AZ
            RR = (X1(N)-AX2)**2 + (Y1(N)-(AY2+DAY))**2 + Z1(N)**2 
            IF (RR.LE.RCUT) N=N+1
          END DO
        END DO    
        NL = NL + 1
        IF (NL.EQ.3) NL = 0
      END DO                 
      NLPOV = NL 
      NATS = N-1 
C      WRITE(*,*)'NATS= ',NATS                             
C
C --------------------------------------------
C       Cu(111) NML-MOVABLE SURFACE LAYERS
C --------------------------------------------                             
C
      N  = 0 
      DO JZ = 1, NML, 1 
        DO JY = 0,2*MY-1
          YJ = ABS(FLOAT(JY))
          DDELY = AMOD(YJ,2.)
          DO JX = 0,MX-1
            N = N + 1
            X(N) = JX + DDELY * AX2 
            Y(N) = JY * AY2 + NL*DAY            
            Z(N) = JZ * AZ
            NSTAT(N) = 1
          END DO
        END DO    
        NL = NL + 1
        IF (NL.EQ.3) NL = 0
      END DO
      WRITE(*,*)'NCUSUB=',N               
      OPEN(11,FILE='C:\Pb-Cu\CU111.DAT')
       WRITE(11,111) (X(I),Y(I),Z(I),I=1,NCUSUB)
      CLOSE(11)
  111 FORMAT(3(E15.8,1X))
C
C ---------------------------------------------      
C   read unmovable substrate layers for POV-RAY format
C ---------------------------------------------                          
C
      N  = 0 
      IF (NLPOV.NE.0) NLPOV = NLPOV - 1
      DO JZ = 0, -1, -1 
        DO JY = 0,2*MY-1
          YJ = ABS(FLOAT(JY))
          DDELY = AMOD(YJ,2.)
          DO JX = 0,MX-1
            N = N + 1
            XP(N) = JX + DDELY * AX2 
            YP(N) = JY * AY2 + NLPOV*DAY            
            ZP(N) = JZ * AZ
          END DO
        END DO    
        IF (NLPOV.EQ.0) THEN
         NLPOV = NLPOV + 2
        ELSE
         NLPOV = NLPOV - 1
        ENDIF 
      END DO

C      Size of Cu and Pb atoms in pov-ray, according to RCU and RPB  

C        IMPORTANT:

         RPOV(1) = RCU/2.0 ! 0.50   ! Cu - radius za chertane w POV	(1/2 ot RCU)
         RPOV(2) = RPB/2.0 ! 0.68 	! Pb - radius za chertane w POV	(1/2 ot RPB)


C -------------------------------------------------------------
C        AD-LAYER !!!!!!!!!!!!!!!
C -------------------------------------------------------------
C       CHOOSE STARTING POSITIONS OF ATOMS
C -------------------------------------------------------------
       CALL SUN (ISEED,RRANF, NAT5)
C
C  RANDOM Pb-ATOMS CONFIGURATION:
C
      N = 0
      IF (NORDER.EQ.0) THEN 
       DO N = NCUSUB+1,NATALL
         X(N) = RRANF(N) * AMX
         Y(N) = RRANF(NAT + N) * AMY
         Z(N) = NML + RPB*RRANF(NAT2 + N) 
         NSTAT(N) = 2
       END DO
       OPEN(27,FILE='INITIAL.DAT')
        WRITE(27,173) (X(I), Y(I), Z(I), I=NCUSUB+1,NATALL)
       CLOSE(27)
  173  FORMAT(3(E15.8,1X)) 
      ENDIF   
C
C      (3) SINGLE WIRE(s):  
C
      N = 0
      IF (NORDER.EQ.1) THEN
      NNN = 0 
      DO N = NCUSUB+1,NATALL
        NNN = NNN + 1

        if(NNN.LE.50)then
           X(N) = 2.9915 + (NNN-1)*RPB !*2*0.684
        elseif(NNN.LE.100)then
           X(N) = 3.514 + (NNN-50-1)*RPB !*2*0.684
        else
           X(N) = 3.5307 + (NNN-100-1)*RPB !*2*0.684
        endif


        if(NNN.LE.50)then
           Y(N) = 14.9993
        elseif(NNN.LE.100)then
           Y(N) = 38.475
        else
           Y(N) = 60.8887
        endif    


        Z(N) = 2.3586 
        NSTAT(N) = 2
      END DO
      OPEN(27,FILE='INITIAL.DAT')
        WRITE(27,173) (X(I), Y(I), Z(I), I=NCUSUB+1,NATALL)
      CLOSE(27)
C     write(*,*) '1-NCUSUB:', NCUSUB, '2-NATALL:', NATALL
  174  FORMAT(3(E15.8,1X)) 
      ENDIF   
C
C  OLD CONFIGURATION
C
      IF (NORDER.EQ.2) THEN      
       IF (NRUNK.EQ.1) THEN
        WRITE(*,*)'OLD DIRECTORY NUMBER '
        READ(*,121) II8,II3,II4
        FM1 = CHAR(II3 + 48)
        FM2 = CHAR(II4 + 48)
        FM3 = CHAR(II8 + 48)
        FOLD = 'C:\Pb-Cu\RESLG'//FM3//FM1//FM2 
        WRITE(*,*) 'INITIAL CONFIGURATION'
        READ(*,123) II3,II4,II10,II11 
  123   FORMAT(I1,I1,I1,I1)
        FM1  = CHAR(II3 + 48)
        FM2  = CHAR(II4 + 48)
        FM3  = CHAR(II10 + 48)
        FM4  = CHAR(II11 + 48)
        COLD = '\A'//FM1//FM2//FM3//FM4//'.DAT' 
        IF (NRUN.GT.1) THEN 
         WRITE(*,*) 'INITIAL CONFIGURATION OF CHAINED RUNS'
         READ(*,123) II3,II4,II10,II11 
         FM1   = CHAR(II3 + 48)
         FM2   = CHAR(II4 + 48) 
         FM3   = CHAR(II10 + 48)
         FM4   = CHAR(II11 + 48)
         COLDC = '\A'//FM1//FM2//FM3//FM4//'.DAT' 
        ENDIF  
       ELSE                 
        II2 = II2 - 1
        IF(II2.LT.0) THEN
         II1 = II1 - 1
         II2 = 9
        ENDIF 
        FM1 = CHAR(II1 + 48)
        FM2 = CHAR(II2 + 48) 
        FM3 = CHAR(II9 + 48)
        FOLD = 'C:\Pb-Cu\RESLG'//FM3//FM1//FM2
        COLD = COLDC
       ENDIF 
       WRITE(*,*) 'OLD DYRECTORY =',FOLD
       NAME = FOLD//COLD 
       OPEN(18,FILE=NAME)
        DO J = 1, NATALL
          READ(18,190) X(J),Y(J),Z(J)
        END DO
       CLOSE(18)            
      ENDIF   
C
      WRITE(*,*)NAME
      NAME = FNAME//'\INITIAL.DAT'
      WRITE(*,*)NAME
      OPEN(27,FILE=NAME)
       WRITE(27,173) (X(I), Y(I), Z(I), I=1,NATALL)
      CLOSE(27)    
C -----------------------------------------------------
C
C      ADLAYER CALCULATION
C -----------------------------------------------------
C -----------------------------------------------------
C
C UPDATING MOVABLE NML-LAYER_S Cu(111) ATOMS AND Pb - ATOMS
C -----------------------------------------------------
      DO IP = 1,NATALL
        NRV(IP) = 0
      END DO                             
      DO IP = 1,NATALL-1
        DO J = IP+1,NATALL 
          XCC = X(J) - X(IP)
          YCC = Y(J) - Y(IP)
          ZCC = Z(J) - Z(IP)
          XCC = XCC - AMX * AINT(XCC*XXIN)
          YCC = YCC - AMY * AINT(YCC*YYIN)
          RADIUS = XCC*XCC+YCC*YCC + ZCC*ZCC
          IF(RADIUS.LE.CUTOFQ) THEN
           NRV(IP) = NRV(IP) + 1
           NRV(J) = NRV(J) + 1
           LIST(IP,NRV(IP)) = J
           LIST(J,NRV(J)) = IP       
          ENDIF      
        END DO
      END DO            
C ------------------------------------------------        
C -----------------------------------------------------
C
C CALCULATE INITIAL ATOM - SUBSTRATE/ADATOM INTERACTION
C -----------------------------------------------------
      EA = 0.            
      DO I = NCUSUB+1,NATALL
        REP  = 0.
        ATR  = 0.       
        DO J = 1, NRV(I)
          JJ = LIST(I,J)
          XCC = X(JJ) - X(I)
          YCC = Y(JJ) - Y(I)
          ZCC = Z(JJ) - Z(I)
          XCC = XCC - AMX * AINT(XCC*XXIN)
          YCC = YCC - AMY * AINT(YCC*YYIN)
          RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
          IF (RADIUS.LE.CUTOFQ) THEN
           IF (NSTAT(I).EQ.NSTAT(JJ)) THEN
            RADIUS = RPB - SQRT(RADIUS) 
            EX1 = EXP(PPB * RADIUS)
            EX2 = EXP(QPB * RADIUS)
            REP = REP + UAPB*EX1
            ATR = ATR + UPB2*EX2 
           ELSE 
            RADIUS = RPBCU - SQRT(RADIUS) 
            EX11 = EXP(PPBCU * RADIUS)
            EX22 = EXP(QPBCU * RADIUS)
            REP  = REP + UAPC*EX11
            ATR  = ATR + UPC2*EX22  
           ENDIF  
          ENDIF
        END DO           
        DO JJ = 1, NATS        
          XCC = X1(JJ) - AMOD(X(I),AX)
          YCC = Y1(JJ) - AMOD(Y(I),AY)
          ZCC = Z1(JJ) - Z(I)
          RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
          IF (RADIUS.LE.CUTOFQ) THEN
           RADIUS = RPBCU - SQRT(RADIUS)
           EX11 = EXP(PPBCU * RADIUS)
           EX22 = EXP(QPBCU * RADIUS)
           REP  = REP + UAPC*EX11
           ATR  = ATR + UPC2*EX22 
          ENDIF        
        END DO        
        ATR   = SQRT(ATR)
        VD(I) = REP - ATR
        EA    = EA + VD(I)    
      END DO             
      EA = EA/NADAT 
      WRITE(*,*)' PB =', EA
C      
      ES  = 0.       
      DO I = 1, NCUSUB 
        REP  = 0.
        ATR  = 0.       
        DO J = 1, NRV(I)
          JJ = LIST(I,J)
          XCC = X(JJ) - X(I)
          YCC = Y(JJ) - Y(I)
          ZCC = Z(JJ) - Z(I)
          XCC = XCC - AMX * AINT(XCC*XXIN)
          YCC = YCC - AMY * AINT(YCC*YYIN)
          RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
          IF (RADIUS.LE.CUTOFQ) THEN
           IF (NSTAT(I).EQ.NSTAT(JJ)) THEN
            RADIUS = RCU - SQRT(RADIUS) 
            EX1 = EXP(PCU * RADIUS)
            EX2 = EXP(QCU * RADIUS)
            REP = REP + UACU*EX1
            ATR = ATR + UCU2*EX2
           ELSE 
            RADIUS = RPBCU - SQRT(RADIUS) 
            EX11 = EXP(PPBCU * RADIUS)
            EX22 = EXP(QPBCU * RADIUS)
            REP  = REP + UAPC*EX11
            ATR  = ATR + UPC2*EX22
           ENDIF  
          ENDIF
        END DO   
        DO JJ = 1, NATS         
          XCC = X1(JJ) - AMOD(X(I),AX)
          YCC = Y1(JJ) - AMOD(Y(I),AY)
          ZCC = Z1(JJ) - Z(I)
          RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
          IF (RADIUS.LE.CUTOFQ) THEN      
           RADIUS = RCU - SQRT(RADIUS) 
           EX1 = EXP(PCU * RADIUS)
           EX2 = EXP(QCU * RADIUS)
           REP = REP + UACU*EX1
           ATR = ATR + UCU2*EX2
          ENDIF        
        END DO               
        ATR   = SQRT(ATR)
        VD(I) = REP - ATR
        ES    = ES + VD(I)
      END DO    
      ES = ES/NCUSUB
      WRITE(*,*)' CU = ',ES
C
C  INPUT PARAMETERS FILE         
C                             
      NAME = FNAME//'\INPUT.TXT'
      OPEN(11,FILE=NAME)                           
       WRITE(11,131)IHR1,':',MMM1,':',MSEC1,' = BEGIN'
       WRITE(11,*)'ISEED = ',ISEED1 
       WRITE(11,*)'NUMBER OF RUNS ',NRUN,' KIND OF RUNS ',NCT
       WRITE(11,*)'INITIAL ORDER = ',NORDER
       WRITE(11,*)''
      WRITE(11,*)'       SYSTEM of 3 Pb nanowires on  Cu(111)'
      WRITE(11,*)'       ***fixed ends***                    '
       WRITE(11,*)'****************************************'
       WRITE(11,*)'  each nanowire = 50 Pb atoms'
      WRITE(11,*)''
      WRITE(11,*)'Simulation Box Size:'
       WRITE(11,*)'X = ',AMX,'  =  ',MX,' Substrate units'
       WRITE(11,*)'Y = ',AMY,'  =  ',MY,' Substrate units'
       WRITE(11,*)'Z = ',AMZ,'  =  ',NML,'Cu-movable monolayers'
       WRITE(11,*)'MC - Steps  =  ', MCSMAX
       WRITE(11,*)'Steps for equilibration = ',NINIT
      WRITE(11,*)''
       WRITE(11,*)'Temperature = ',T
      WRITE(11,*)''
      WRITE(11,*)'Time interval snapshot = ', NWSNA1
      WRITE(11,*)'' 
      WRITE(11,*)'Time interval calculation vac/holes = ', TIMENVAC
      WRITE(11,*)''
       WRITE(11,*)'Cu - Cu interaction = ', UCU 
      WRITE(11,*)''
      WRITE(11,*)'Pb - Pb interaction = ', UPB
      WRITE(11,*)''
      WRITE(11,*)'Substrate atoms = ',NCUSUB,' Adatom = ',NADAT
      WRITE(11,*)'COVERAGE = ',COV1
      WRITE(11,*)'Total Number of Atoms = ',NATALL
       WRITE(11,*)'Energy per adatom =  ', EA 
       WRITE(11,*)'Energy per substrate atom =  ', ES
      CLOSE(11)
C
 
      NM  = 0
      NM1 = 0        
      INM = 0
      KLM = 0
C ************************************************************      
C ------------------------------------------------------------
C     START MONTE CARLO SIMULATION  <<<<<<<<<<<<<<<<<<<<<<<<
C ------------------------------------------------------------                                                            
C
      write(*,*) ''
      write(*,*) 'START'
      write(*,*) ''

      WRITE(6,*)'EQUILIBRATION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
C
      MCS1 = NM + 1
      MCSTEP = NINIT
      IWRITE = 0     
 3000 IF(NM.LE.NINIT) THEN
       EAVAD  = 0.
       EAVAD2 = 0.
       EAVER  = 0.
       EAVER2 = 0. 
       EAVSU  = 0.
       EAVAL  = 0.
       ENPE   = 0.
       FAVER  = 0.
       PSTEPS = 0.
       IUPDAT = 1
       NACM   = 0
       NACALM = 0
      ENDIF
      IF(NM.GT.NINIT)THEN
       IWTITE = 1
       MCSTEP = MCSMAX
       NM1 = NINIT
       MCS1 = NM - NINIT + 1
      ENDIF



      NAME = FNAME//'\t-v-h.dat'
      open(111, file = NAME)


C --------------------------------------------------
C  MAIN LOOP
C --------------------------------------------------
      DO 2000 MCS = MCS1, MCSTEP! po vremeto

        NM = NM1 + MCS
        CALL SUN(ISEED,RRANF,NAT5)
        IF(IUPDAT.EQ.1) THEN
C-------------------------------------------------------------------
C                 UPDATING OF THE NEIGHBOURS
C-------------------------------------------------------------------
         DO IP = 1, NATALL
           DIX(IP) = 0.
           DIY(IP) = 0.
           DIZ(IP) = 0.
           NRV(IP) = 0
         END DO
         DO IP = 1,NATALL-1                       
           DO J = IP+1,NATALL
             XCC = X(J) - X(IP)
             YCC = Y(J) - Y(IP)
             ZCC = Z(J) - Z(IP)
             XCC = XCC - AMX * AINT(XCC*XXIN)
             YCC = YCC - AMY * AINT(YCC*YYIN)
             RADIUS = XCC*XCC+YCC*YCC + ZCC*ZCC
             IF(RADIUS.LE.RMQ) THEN
              NRV(IP) = NRV(IP) + 1
              NRV(J) = NRV(J) + 1
              LIST(IP,NRV(IP)) = J
              LIST(J,NRV(J)) = IP
             ENDIF
           END DO
         END DO 
         IUPDAT = 0
        ENDIF          
C         
C----------- END UPDATING -----------------
C
C--------------------------------------------------
C--------------------------------------------------
C LOOP THROUGH ALL PARTICLES
C-------------------------------------------------- 
        DDD   = DMAX(2)
        DDD2  = DMAX2(2)
        SKINN = SKINQ(2) 

*****************************************************
C     Izchisliava Nvac1,2,3 i Nhole1,2,3 i osredniava.
C      Zapis v file


      IF (mod(MCS, TIMENVAC).EQ.0) THEN	  ! Ako triabva da izchisliava

C      Nylirane
      dddd1 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_1
      Nvac1 = 0	 ! nylirane predi da skanirame VAC za verijka_1
      Nhole1 = 0	 ! nylirane predi da skanirame HOLES za verijka_1

      dddd2 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_2
      Nvac2 = 0	 ! nylirane predi da skanirame VAC za verijka_2
      Nhole2 = 0	 ! nylirane predi da skanirame HOLES za verijka_2

      dddd3 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_3
      Nvac3 = 0	 ! nylirane predi da skanirame VAC za verijka_3
      Nhole3 = 0	 ! nylirane predi da skanirame HOLES za verijka_3


      V1 = 4.0*RPOV(2)	  ! "dolna" granica za definiciq na vacancy
      V2 = 6.0*RPOV(2)	  !	"gorna" granica za definiciq na vacancy

C	 End nylirane
C      IPP e broiach po olovnite atomi

         DO IPP= NCUSUB+1, NCUSUB+49  ! skanirame po Pb atoms 1 wire

      dddd1=sqrt((X(IPP)-X(IPP+1))**2+(Y(IPP)-Y(IPP+1))**2+
     &(Z(IPP)-Z(IPP+1))**2)


C	write(*,*) dddd1

         if (dddd1.GE.V1.AND.dddd1.LE.V2) then
          Nvac1 = Nvac1 + 1	! za vacancii po purva wire
         endif
         
              if(dddd1.GT.V2)then
               Nhole1 = Nhole1 + 1 ! za dypki po purva wire
      
              endif

         ENDDO	! po Pb atoms 1 wire
******************
         DO IPP= NCUSUB+51, NCUSUB+99  ! skanirame po Pb atoms 2 wire

      dddd2=sqrt((X(IPP)-X(IPP+1))**2+(Y(IPP)-Y(IPP+1))**2+
     &(Z(IPP)-Z(IPP+1))**2)

C	write(*,*) dddd2

         if (dddd2.GE.V1.AND.dddd2.LE.V2) then
          Nvac2 = Nvac2 + 1		! za vacancii po vtora wire
         endif

              if(dddd2.GT.V2)then
               Nhole2 = Nhole2 + 1		! za dypki po vtora wire
              endif

         ENDDO	! po Pb atoms 2 wire
******************
         DO IPP= NCUSUB+101, NCUSUB+149  ! skanirame po Pb atoms 3 wire

      dddd3=sqrt((X(IPP)-X(IPP+1))**2+(Y(IPP)-Y(IPP+1))**2+
     &(Z(IPP)-Z(IPP+1))**2)


C	write(*,*) dddd3

         if (dddd3.GE.V1.AND.dddd3.LE.V2) then
          Nvac3 = Nvac3 + 1		! za vacancii po treta wire
         endif

              if(dddd3.GT.V2)then
               Nhole3 = Nhole3 + 1		! za dypki  po treta wire
              endif

         ENDDO	! po Pb atoms 3 wire
******************

      Nvacsredno = ((Nvac1 + Nvac2 + Nvac3)*1.0)/3.0        ! REAL*4
      Nholesredno = ((Nhole1 + Nhole2 + Nhole3)*1.0)/3.0    ! REAL*4

      write(*,222) MCS, Nvac1, Nvac2, Nvac3, Nhole1, Nhole2, Nhole3,! na ekran
     & Nvacsredno, Nholesredno
      write(*,*) '' 
C	Zapis v file:
      write (111,222) MCS, Nvac1, Nvac2, Nvac3, Nhole1, Nhole2, Nhole3, 
     & Nvacsredno, Nholesredno   

      ENDIF		                           ! Ako triabva da izchisliava
**************************************************
C      Krai na presmiqtaneto na vanacii i dypki
**************************************************





        DO 2050 IP = 1, NATALL  ! po wsichki chastici
          IF (IP.EQ.NCUSUB+1) THEN
           DDD   = DMAX(1)
           DDD2  = DMAX2(1)
           SKINN = SKINQ(1)
          ENDIF 
          XSOLD  = X(IP)
          YSOLD  = Y(IP)
          ZSOLD  = Z(IP)        
C


C       Move particle IP at random ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C

C     Trial move along X axis:
      IF(IP.EQ.NCUSUB+1.OR.IP.EQ.NCUSUB+50.OR.IP.EQ.NCUSUB+51.OR.
     &   IP.EQ.NCUSUB+100.OR.IP.EQ.NCUSUB+101.OR.IP.EQ.NCUSUB+150) THEN
         DX = 0.0
      ELSE
         DX = RRANF(IP)*DDD2 - DDD
         XSNEW = DX + XSOLD
      ENDIF

C     Trial move along Y axis:
      IF(IP.EQ.NCUSUB+1.OR.IP.EQ.NCUSUB+50.OR.IP.EQ.NCUSUB+51.OR.
     &   IP.EQ.NCUSUB+100.OR.IP.EQ.NCUSUB+101.OR.IP.EQ.NCUSUB+150) THEN
         DY = 0.0
      ELSE
         DY = RRANF(NAT+IP)*DDD2 - DDD
         YSNEW = DY + YSOLD
      ENDIF

C     Trial move along Z axis:
      IF(IP .EQ. NCUSUB+1 .OR. IP .EQ. NCUSUB+50 .OR. IP .EQ. NCUSUB+51 
     &   .OR. IP .EQ. NCUSUB+100 .OR. IP .EQ. NCUSUB+101 .OR. IP .EQ. 
     &    NCUSUB+150) THEN
         DZ = 0.0
      ELSE
         DZ = RRANF(NAT4+IP)*DDD2 - DDD
         ZSNEW = DZ + ZSOLD
      ENDIF

C
C       Periodic boundary conditions in X and Y direction:
C       Reflection in Z direction:
C
          XSS = XXIN * XSNEW - 1.
          XSNEW = XSNEW - AMX * AINT(XSS)
          YSS = YYIN * YSNEW - 1.
          YSNEW = YSNEW - AMY * AINT(YSS)
          ZSS = AMZ2 * AINT(ZSNEW*RAMZ)
          ZSNEW = ABS(ZSS - ZSNEW)
          DZ = ZSNEW - ZSOLD
C
C       MOVE A PARTICLE  
C
C       CALCULATE OLD ADSORBATE-ADSORBATE INTERACTIONS 'EADAD'
C
          IF (IP.LE.NCUSUB) THEN
           REP  = 0.
           ATR  = 0.       
           DO J = 1, NRV(IP)
             JJ = LIST(IP,J)
             XCC = X(JJ) - XSOLD
             YCC = Y(JJ) - YSOLD
             ZCC = Z(JJ) - ZSOLD
             XCC = XCC - AMX * AINT(XCC*XXIN)
             YCC = YCC - AMY * AINT(YCC*YYIN)
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
             IF (RADIUS.LE.CUTOFQ) THEN
              IF (NSTAT(IP).EQ.NSTAT(JJ)) THEN
               RADIUS = RCU - SQRT(RADIUS) 
               EX1 = EXP(PCU * RADIUS)
               EX2 = EXP(QCU * RADIUS)
               REP = REP + UACU*EX1
               ATR = ATR + UCU2*EX2
              ELSE 
               RADIUS = RPBCU - SQRT(RADIUS) 
               EX11 = EXP(PPBCU * RADIUS)
               EX22 = EXP(QPBCU * RADIUS)
               REP  = REP + UAPC*EX11
               ATR  = ATR + UPC2*EX22
              ENDIF  
             ENDIF
           END DO                      
           DO JJ = 1, NATS        
             XCC = X1(JJ) - AMOD(XSOLD,AX)
             YCC = Y1(JJ) - AMOD(YSOLD,AY)
             ZCC = Z1(JJ) - ZSOLD
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
             IF (RADIUS.LE.CUTOFQ) THEN
              RADIUS = RCU - SQRT(RADIUS) 
              EX1 = EXP(PCU * RADIUS)
              EX2 = EXP(QCU * RADIUS)
              REP = REP + UACU*EX1
              ATR = ATR + UCU2*EX2 
             ENDIF        
           END DO               
           ATR    = SQRT(ATR)         
           EADAD0  = REP - ATR
          ELSE          
           REP  = 0.
           ATR  = 0.       
           DO J = 1, NRV(IP)
             JJ = LIST(IP,J)
             XCC = X(JJ) - XSOLD
             YCC = Y(JJ) - YSOLD
             ZCC = Z(JJ) - ZSOLD  
             XCC = XCC - AMX * AINT(XCC*XXIN)
             YCC = YCC - AMY * AINT(YCC*YYIN)
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
             IF (RADIUS.LE.CUTOFQ) THEN
              IF (NSTAT(IP).EQ.NSTAT(JJ)) THEN 
               RADIUS = RPB - SQRT(RADIUS) 
               EX1 = EXP(PPB * RADIUS)
               EX2 = EXP(QPB * RADIUS)
               REP = REP + UAPB*EX1
               ATR = ATR + UPB2*EX2
              ELSE 
               RADIUS = RPBCU - SQRT(RADIUS) 
               EX11 = EXP(PPBCU * RADIUS)
               EX22 = EXP(QPBCU * RADIUS)
               REP  = REP + UAPC*EX11
               ATR  = ATR + UPC2*EX22
              ENDIF  
             ENDIF
           END DO                
           DO JJ = 1, NATS        
             XCC = X1(JJ) - AMOD(XSOLD,AX)
             YCC = Y1(JJ) - AMOD(YSOLD,AY)
             ZCC = Z1(JJ) - ZSOLD
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
             IF (RADIUS.LE.CUTOFQ) THEN
              RADIUS = RPBCU - SQRT(RADIUS) 
              EX11 = EXP(PPBCU * RADIUS)
              EX22 = EXP(QPBCU * RADIUS)
              REP  = REP + UAPC*EX11
              ATR  = ATR + UPC2*EX22 
             ENDIF        
           END DO   
           ATR    = SQRT(ATR)         
           EADAD0  = REP - ATR  
          ENDIF
C                                                       
C       CALCULATE NEW ADSORBATE-ADSORBATE INTERACTIONS 'EADAD'
C
          IF (IP.LE.NCUSUB) THEN
           REP  = 0.
           ATR  = 0.       
           DO J = 1, NRV(IP)
             JJ = LIST(IP,J)
             XCC = X(JJ) - XSNEW
             YCC = Y(JJ) - YSNEW
             ZCC = Z(JJ) - ZSNEW
             XCC = XCC - AMX * AINT(XCC*XXIN)
             YCC = YCC - AMY * AINT(YCC*YYIN)
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
             IF (RADIUS.LE.CUTOFQ) THEN
              IF (NSTAT(IP).EQ.NSTAT(JJ)) THEN
               RADIUS = RCU - SQRT(RADIUS) 
               EX1 = EXP(PCU * RADIUS)
               EX2 = EXP(QCU * RADIUS)
               REP = REP + UACU*EX1
               ATR = ATR + UCU2*EX2
              ELSE 
               RADIUS = RPBCU - SQRT(RADIUS) 
               EX11 = EXP(PPBCU * RADIUS)
               EX22 = EXP(QPBCU * RADIUS)
               REP  = REP + UAPC*EX11
               ATR  = ATR + UPC2*EX22
              ENDIF  
             ENDIF
           END DO                      
           DO JJ = 1, NATS        
             XCC = X1(JJ) - AMOD(XSNEW,AX)
             YCC = Y1(JJ) - AMOD(YSNEW,AY)
             ZCC = Z1(JJ) - ZSNEW
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
             IF (RADIUS.LE.CUTOFQ) THEN
              RADIUS = RCU - SQRT(RADIUS) 
              EX1 = EXP(PCU * RADIUS)
              EX2 = EXP(QCU * RADIUS)
              REP = REP + UACU*EX1
              ATR = ATR + UCU2*EX2 
             ENDIF        
           END DO               
           ATR    = SQRT(ATR)         
           EADAD  = REP - ATR
           EDREAL = (EADAD - EADAD0)/RT
          ELSE                   
           REP  = 0.
           ATR  = 0.       
           DO J = 1, NRV(IP)
             JJ = LIST(IP,J)
             XCC = X(JJ) - XSNEW
             YCC = Y(JJ) - YSNEW
             ZCC = Z(JJ) - ZSNEW      
             XCC = XCC - AMX * AINT(XCC*XXIN)
             YCC = YCC - AMY * AINT(YCC*YYIN)
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
             IF (RADIUS.LE.CUTOFQ) THEN
              IF (NSTAT(IP).EQ.NSTAT(JJ)) THEN
               RADIUS = RPB - SQRT(RADIUS) 
               EX1 = EXP(PPB * RADIUS)
               EX2 = EXP(QPB * RADIUS)
               REP = REP + UAPB*EX1
               ATR = ATR + UPB2*EX2 
              ELSE 
               RADIUS = RPBCU - SQRT(RADIUS) 
               EX11 = EXP(PPBCU * RADIUS)
               EX22 = EXP(QPBCU * RADIUS)
               REP  = REP + UAPC*EX11
               ATR  = ATR + UPC2*EX22
              ENDIF  
             ENDIF
           END DO                          
           DO JJ = 1, NATS        
             XCC = X1(JJ) - AMOD(XSNEW,AX)
             YCC = Y1(JJ) - AMOD(YSNEW,AY)
             ZCC = Z1(JJ) - ZSNEW
             RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
             IF (RADIUS.LE.CUTOFQ) THEN
              RADIUS = RPBCU - SQRT(RADIUS) 
              EX11 = EXP(PPBCU * RADIUS)
              EX22 = EXP(QPBCU * RADIUS)
              REP  = REP + UAPC*EX11
              ATR  = ATR + UPC2*EX22 
             ENDIF        
           END DO   
           ATR    = SQRT(ATR)         
           EADAD  = REP - ATR
           EDREAL = (EADAD - EADAD0)/RT
          ENDIF
C          
          NACALM = NACALM + 1
          IF (EDREAL.LE.0.) THEN
           NACM    = NACM + 1
           X(IP)   = XSNEW
           Y(IP)   = YSNEW
           Z(IP)   = ZSNEW 
           DIX(IP) = DIX(IP) + DX
           DIY(IP) = DIY(IP) + DY
           DIZ(IP) = DIZ(IP) + DZ
           RSQUAR  = DIX(IP)*DIX(IP)+DIY(IP)*DIY(IP)+DIZ(IP)*DIZ(IP)
           IF (RSQUAR.GT.SKINN) IUPDAT = 1
          ENDIF
          IF (EDREAL.GT.0.) THEN
           RRR = ALOG(RRANF(NAT2+IP)) + EDREAL
           IF (RRR.LE.0.) THEN
            NACM    = NACM + 1
            X(IP)   = XSNEW
            Y(IP)   = YSNEW
            Z(IP)   = ZSNEW 
            DIX(IP) = DIX(IP) + DX
            DIY(IP) = DIY(IP) + DY
            DIZ(IP) = DIZ(IP) + DZ
            RSQUAR  = DIX(IP)*DIX(IP)+DIY(IP)*DIY(IP)+DIZ(IP)*DIZ(IP)
            IF (RSQUAR.GT.SKINN) IUPDAT = 1
           ENDIF
          ENDIF               
 2050 CONTINUE
C ------------------------------------------------------------------
C  END OF THE LOOP THROUGH ALL PARTICLES
C ------------------------------------------------------------------                                                  
C
C------------------------------------------------------------------
C       CALCULATE ENERGY , SDFunction
C------------------------------------------------------------------
      EG = 0.
      IF(NM1.NE.0) THEN 
       DO L = 1, NML
         J = (L-1)*MLE 
         DO I = J+1,J+MLE
           DO JJ = J+1,J+MLE
             IF (JJ.NE.I) THEN
              XCC = X(JJ) - X(I)
              YCC = Y(JJ) - Y(I)
              ZCC = Z(JJ) - Z(I)
              XCC = XCC - AMX * AINT(XCC*XXIN)
              YCC = YCC - AMY * AINT(YCC*YYIN)
              RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
              IF (RADIUS.LE.16.0) THEN
               NXXX = INT(SQRT(RADIUS) * DZZ)
               PDF(L,NXXX) = PDF(L,NXXX) + 1 
              ENDIF
             ENDIF 
           END DO
         END DO
       END DO       
       DO I = NCUSUB+1, NATALL
         DO JJ = NCUSUB+1, NATALL
           IF (JJ.NE.I) THEN
            XCC = X(JJ) - X(I)
            YCC = Y(JJ) - Y(I)
            ZCC = Z(JJ) - Z(I)
            XCC = XCC - AMX * AINT(XCC*XXIN)
            YCC = YCC - AMY * AINT(YCC*YYIN)
            RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
            IF (RADIUS.LE.16.0) THEN
             NXXX = INT(SQRT(RADIUS) * DZZ)
             PDF(NML1,NXXX) = PDF(NML1,NXXX) + 1 
            ENDIF     
          ENDIF 
         END DO               
       END DO    
       DO IP = 1, NATALL
         NZZZ = INT(Z(IP) * DZZ) + 1
         SDF (NZZZ) = SDF(NZZZ) + 1  
       END DO
      ENDIF                            
C
C  AVERAGE ENERGY CALCULATION
C
      EAL = 0.
      EPB = 0.
      ECU = 0.
      DO I = 1, NATALL
        REP  = 0.
        ATR  = 0.       
        DO J = 1, NRV(I)
          JJ = LIST(I,J)
          XCC = X(JJ) - X(I)
          YCC = Y(JJ) - Y(I)
          ZCC = Z(JJ) - Z(I)
          XCC = XCC - AMX * AINT(XCC*XXIN)
          YCC = YCC - AMY * AINT(YCC*YYIN)
          RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC  
          IF (RADIUS.LE.CUTOFQ) THEN
           IF (NSTAT(I).EQ.NSTAT(JJ)) THEN
            IF (NSTAT(I).EQ.1) THEN
             RADIUS = RCU - SQRT(RADIUS) 
             EX1 = EXP(PCU * RADIUS)
             EX2 = EXP(QCU * RADIUS)
             REP = REP + UACU*EX1
             ATR = ATR + UCU2*EX2
            ELSE                                
             RADIUS = RPB - SQRT(RADIUS)
             EX1 = EXP(PPB * RADIUS)
             EX2 = EXP(QPB * RADIUS)
             REP = REP + UAPB*EX1
             ATR = ATR + UPB2*EX2
            ENDIF 
           ELSE 
            RADIUS = RPBCU - SQRT(RADIUS) 
            EX11 = EXP(PPBCU * RADIUS)
            EX22 = EXP(QPBCU * RADIUS)
            REP  = REP + UAPC*EX11
            ATR  = ATR + UPC2*EX22
           ENDIF  
          ENDIF
        END DO   
        IF (NSTAT(I).EQ.1) THEN
         DO JJ = 1, NATS         
           XCC = X1(JJ) - AMOD(X(I),AX)
           YCC = Y1(JJ) - AMOD(Y(I),AY)
           ZCC = Z1(JJ) - Z(I)
           RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
           IF (RADIUS.LE.CUTOFQ) THEN      
            RADIUS = RCU - SQRT(RADIUS) 
            EX1 = EXP(PCU * RADIUS)
            EX2 = EXP(QCU * RADIUS)
            REP = REP + UACU*EX1
            ATR = ATR + UCU2*EX2
           ENDIF        
         END DO 
        ELSE 
         DO JJ = 1, NATS        
           XCC = X1(JJ) - AMOD(X(I),AX)
           YCC = Y1(JJ) - AMOD(Y(I),AY)
           ZCC = Z1(JJ) - Z(I)
           RADIUS = XCC*XCC + YCC*YCC + ZCC*ZCC
           IF (RADIUS.LE.CUTOFQ) THEN
            RADIUS = RPBCU - SQRT(RADIUS)
            EX11 = EXP(PPBCU * RADIUS)
            EX22 = EXP(QPBCU * RADIUS)
            REP  = REP + UAPC*EX11
            ATR  = ATR + UPC2*EX22
           ENDIF        
         END DO  
        ENDIF                   
C              
        ATR = SQRT(ATR)
        IF (I.LE.NCUSUB) THEN
         ECU = ECU + REP - ATR 
        ELSE
         EPB = EPB + REP - ATR
        ENDIF 
        EAL = EAL + REP - ATR 
      END DO    
      EAL = EAL / NATALL 
      ECU = ECU / NCUSUB 
      EPB = EPB / NADAT
      EG  = EPB / RT
      
C     
      EAVER  = EAVER  + EG
      EAVER2 = EAVER2 + EG * EG 
      EAVSU  = EAVSU  + ECU
      EAVAL  = EAVAL  + EAL
C                    
      EG = EG * RT
C                     
  190  FORMAT(3(E15.8,1X))
  290  FORMAT(5(E15.8,1X))
  390  FORMAT(4(E15.8,1X))
       IF(NM.EQ.NWSNAA) THEN
        NWSNAA = NWSOOD + NM
        IF (NM1.NE.0) THEN 
         INM4 = INM/1000
         INM5 = INM - INM4*1000
         INM3 = INM5/100   
         INM6 = INM5 - INM3*100 
         INM1 = INM6/10 
         INM2 = INM6 - INM1*10
         INM  = INM + 1       
         FM4 = CHAR(INM4 + 48)
         FM3 = CHAR(INM3 + 48)
         FM1 = CHAR(INM1 + 48)
         FM2 = CHAR(INM2 + 48) 

         NAME = FNAME//'\A'//FM4//FM3//FM1//FM2//'.DAT'
         WRITE(*,*)'MCS=',MCS, NAME 
         WRITE(*,*)'T = ',T  
         WRITE(*,*)'COVERAGE = ',COV  
C	    ************************* DATA FILES ***********************                      
C         OPEN(19,FILE=NAME)
C         WRITE(19,290)(X(I),Y(I),Z(I),Y(I),Z(I),I=1,NATALL)
C         CLOSE(19)           
C         
C  SNAPSHOT IN POV-RAY FORMAT         
C 
        F1 = 'sphere{<'
        F2 = ','
        F3 = '>'
        F4 = 'pigment'
        F5 = '{'
        F6 = '}}' 
        F10 = '}' 
         NAME = FNAME//'\A'//FM4//FM3//FM1//FM2//'.POV'  
         DO L=1,NATALL
           K = MLE2 + L
           XP(K) = X(L)
           YP(K) = Y(L)
           ZP(K) = Z(L)
         END DO
C
      OPEN(22,FILE=NAME)
C-----------------------------------------------------
C  POV-RAY HEADER 
C------------------------------------------------------------------

C       WRITE(22,*)'// Persistence Of Vision Raytracer version 2.0 sample file.'
        WRITE(22,*)'#include "colors.inc" '
        WRITE(22,*)'#include "textures.inc" '
        WRITE(22,*)'#include "shapes.inc" '
        WRITE(22,*)'global_settings { assumed_gamma 1.0 }'


        WRITE(22,*)'camera { '
        WRITE(22,*)'location  <-180, 320, -178>'
        WRITE(22,*)'direction <2, 4, 1>'
        WRITE(22,*)'look_at   <35, -1, 35>'
        WRITE(22,*)'}'

        WRITE(22,*)'background {White}'

        WRITE(22,*)'//light_source {<-500, 500, -1500> color White}'


       WRITE(22,*)'light_source { <8, 55, 5> color White}'
       WRITE(22,*)'//light_source { <8, 12, -10> color White}'
C ---------------------------------------------------------------
       WRITE(22,*)'//Cu - UPPER FIX ATOM LAYER ----- '
       F7 = 'Blue'	 

       DO K=1,MLE
         WRITE(22,220)F1,XP(K),F2,ZP(K),F2,YP(K),F3,F2,RPOV(1),
     +    F4,F5,F7,F6
       END DO  

       WRITE(22,*)'//END Cu - FIX ATOMS -------------------- ' 
       WRITE(22,*)'//Cu - LOWER FIX ATOM LAYER----------- '
       F7 = 'Blue'

       DO K= MLE+1,MLE2
         WRITE(22,220)F1,XP(K),F2,ZP(K),F2,YP(K),F3,F2,RPOV(1),
     +    F4,F5,F7,F6
       END DO  

       WRITE(22,*)'//END Cu - UPPER FIX ATOMS -------------------- '
       WRITE(22,*)'//Cu - MOVE ATOMS - UPPER LAYER---------- '
       F7 = 'Green'

       DO K=MLE2+1,NCUAT
         WRITE(22,220)F1,XP(K),F2,ZP(K),F2,YP(K),F3,F2,RPOV(1),
     +    F4,F5,F7,F6
       END DO  

       WRITE(22,*)'//END Cu - MOVE ATOMS - ----------------- ' 
       WRITE(22,*)'//Pb -  ATOMS - ------------------------- '
       F7 = 'Gray'



       DO K= NCUAT+1,NATALL+MLE2



         WRITE(22,221)F1,XP(K),F2,ZP(K),F2,YP(K),F3,F2,RPOV(2),
     +    F4,F5,F7,F10                         
        WRITE(22,*) 'finish {'
        WRITE(22,*) 'ambient .3'
        WRITE(22,*) 'diffuse .8'
        WRITE(22,*) 'phong .65'
        WRITE(22,*) 'phong_size 7 }}'

       END DO  

c	 COUNTING dddd: 
c      Save in file: time, Nvac(1-2-3), Nhol(1-2-3)

c	dddd1 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_1
c	Nvac1 = 0	 ! nylirane predi da skanirame VAC za verijka_1
c	Nhole1 = 0	 ! nylirane predi da skanirame HOLES za verijka_1
c
c	dddd2 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_2
c	Nvac2 = 0	 ! nylirane predi da skanirame VAC za verijka_2
c	Nhole2 = 0	 ! nylirane predi da skanirame HOLES za verijka_2
c
c	dddd3 = 0.0	 ! nylirane predi da skanirame DISTANCES za verijka_3
c	Nvac3 = 0	 ! nylirane predi da skanirame VAC za verijka_3
c	Nhole3 = 0	 ! nylirane predi da skanirame HOLES za verijka_3
c
c
c	V1 = 4.0*RPOV(2)	  ! "dolna" granica za definiciq na vacancy
c	V2 = 8.0*RPOV(2)	  !	"gorna" granica za definiciq na vacancy
*******************
c         DO K= NCUAT+1, NCUAT+49  ! skanirame po Pb atoms 1 wire
c
c	dddd1=sqrt((XP(K)-XP(K+1))**2+(YP(K)-YP(K+1))**2+
c	&(ZP(K)-ZP(K+1))**2)
c
c	   if (dddd1.GE.V1.AND.dddd1.LE.V2) then
c	    Nvac1 = Nvac1 + 1		! za vacancii
c	   elseif(dddd1.GE.V2)then
c	    Nhole1 = Nhole1 + 1		! za dypki
c	   endif
c
c         ENDDO	! po Pb atoms 1 wire
******************
c         DO K= NCUAT+51, NCUAT+99  ! skanirame po Pb atoms 2 wire
c
c	dddd2=sqrt((XP(K)-XP(K+1))**2+(YP(K)-YP(K+1))**2+
c	&(ZP(K)-ZP(K+1))**2)
c
c	   if (dddd2.GE.V1.AND.dddd2.LE.V2) then
c	    Nvac2 = Nvac2 + 1		! za vacancii
c	   elseif(dddd2.GE.V2)then
c	    Nhole2 = Nhole2 + 1		! za dypki
c	   endif
c
c         ENDDO	! po Pb atoms 2 wire
******************
c         DO K= NCUAT+101, NCUAT+149  ! skanirame po Pb atoms 3 wire
c
c	dddd3=sqrt((XP(K)-XP(K+1))**2+(YP(K)-YP(K+1))**2+
c	&(ZP(K)-ZP(K+1))**2)
c
c	   if (dddd3.GE.V1.AND.dddd3.LE.V2) then
c	    Nvac3 = Nvac3 + 1		! za vacancii
c	   elseif(dddd3.GE.V2)then
c	    Nhole3 = Nhole3 + 1		! za dypki
c	   endif
c
c         ENDDO	! po Pb atoms 3 wire
******************
c
c	Nvacsredno = ((Nvac1 + Nvac2 + Nvac3)*1.0)/3.0        ! REAL*4
c	Nholesredno = ((Nhole1 + Nhole2 + Nhole3)*1.0)/3.0    ! REAL*4
c
c	write(*,222) MCS, Nvac1, Nvac2, Nvac3, Nhole1, Nhole2, Nhole3,	! na ekran
c	& Nvacsredno, Nholesredno
c	write(*,*) '' 
C	Zapis v file:
c	write (111,222) MCS, Nvac1, Nvac2, Nvac3, Nhole1, Nhole2, Nhole3, 
c	& Nvacsredno, Nholesredno   



       WRITE(22,*)'//END Pb -  ATOMS - --------------------- '
      CLOSE(22) 
  220 FORMAT(A8,2(F7.4,A1),F7.4,A1,A1,F5.2,1X,A7,1X,A1,A5,A2)
  221 FORMAT(A8,2(F7.4,A1),F7.4,A1,A1,F5.2,1X,A7,1X,A1,A5,A1)
222   FORMAT(I8, 6(I4), F6.2, F6.2)
C END OF SNAPSHOT IN POV-RAY FORMAT   
        ENDIF          
       ENDIF        
       IF(NM.EQ.NWCOOR) THEN
        NWCOOR = NWCOOD + NM
        KLM = KLM + 1
         NAME = FNAME//'\ENERGY.DAT'
         WRITE(*,*)'MCS=',MCS,' Total Energy =',EAL  
C         CALL GETTIM(IHR,IMMM,ISEC,I100TH)
         CALL DATE_AND_TIME(values=values)
         IHR = values(5)    ! Hour
         IMMM = values(6)   ! Minute
         ISEC = values(7)   ! Second
         I100TH = 0         ! No hundredths of seconds available
         WRITE(*,130)IRUN-1,IHR,':',IMMM,':',ISEC,'=',NAME
  130    FORMAT(I5,5X,I2,1A,I2,1A,I2,1A,20A)  
         MHR = IHR
         MMMM = IMMM
         MSEC = ISEC 
         ENERPB(KLM) = EPB 
         ENERCU(KLM) = ECU
         ENERAL(KLM) = EAL
         NMCSTP(KLM) = NM
         OPEN(19,FILE=NAME)
          WRITE(19,390)(NMCSTP(I),ENERPB(I),ENERCU(I),ENERAL(I),I=1,KLM)
         CLOSE(19)             
       ENDIF
C
 2000  CONTINUE                                                     
C---------------------------------------------------
C END OF THE MAIN LOOP
C--------------------------------------------------- 
C
       GOTO (3003) IWRITE
       WRITE(6,*)'COLLECTING AVERAGES<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
       MCSTEP = MCSMAX
       IWRITE = 1
C       NWCOOR = NWCOOR - NINIT
       NM1 = NINIT
       MCS1 = 1
       GOTO 3000



      CLOSE(111)
C
C--------------------------------------------------
C *********** FINAL CALCULATIONS ******************
C--------------------------------------------------
C
 3003  FMCMAX =  FLOAT(MCSMAX)
       EAVERS =  EAVER / FMCMAX * RT
       EAVES2 =  EAVER2 / FMCMAX * RT * RT
       EE = EAVES2 - EAVERS * EAVERS
       SPHEAT = EE / NADAT / RT / T
       SPHEAT = ABS(SPHEAT)
       EAVSU = EAVSU / FMCMAX
       EAVAL = EAVAL / FMCMAX 
C                                           
       WRITE(*,131)IHR1,':',MMM1,':',MSEC1,' = BEGIN'
       WRITE(*,131)IHR,':',IMMM,':',ISEC,' = END  '
  131  FORMAT(5X,I2,1A,I2,1A,I2,A8)        
       NAME = FNAME//'\RESLG.TXT'
       OPEN(11,FILE=NAME)                           
       WRITE(11,131)IHR1,':',MMM1,':',MSEC1,' = BEGIN'
       WRITE(11,131)IHR,':',IMMM,':',ISEC,' = END  '  
       WRITE(11,*)'ISEED = ',ISEED1 
       WRITE(11,*)'NUMBER OF RUNS ',NRUN,' KIND OF RUNS ',NCT
       WRITE(11,*)'INITIAL ORDER = ',NORDER
       WRITE(11,*)'       SYSTEM   Pb  on  Cu(111)'
       WRITE(11,*)'****************************************'
       WRITE(11,*)'Simulation Box Size:'
       WRITE(11,*)'X = ',AMX,'  =  ',MX,' Substrate units'
       WRITE(11,*)'Y = ',AMY,'  =  ',MY,' Substrate units'
       WRITE(11,*)'Z = ',AMZ,'  =  ',NML,' Monolayers'
       WRITE(11,*)'MC - Steps  =  ', MCSMAX
      WRITE(11,*)'****************************************'
       WRITE(11,*)'Steps for equilibration = ',NINIT
       WRITE(11,*)'Temperatur (K) = ',T
       WRITE(11,*)'Trial States Accepted (M) = ',FLOAT(NACM)/NACALM
       WRITE(11,*)'Energy per adatom =  ', EAVERS/NADAT 
       WRITE(11,*)'Energy per substrate atom =  ', EAVSU
       WRITE(11,*)'Energy per atom of the system =  ', EAVAL
       WRITE(11,*)'Substrate atoms = ',NCUSUB,' Adatom = ',NADAT
           WRITE(11,*)'COVERAGE (ML) = ',COV1
           WRITE(11,*)'Specific Heat = ',SPHEAT
           WRITE(11,*)'TB POTENTIAL PARAMETERS:'
           WRITE(11,*)'****************************************'
      WRITE(11,*)'Cu-Cu ksi = ',UCU
      WRITE(11,*)'Cu-Cu p = ',PCU
      WRITE(11,*)'Cu-Cu q = ',QCU
      WRITE(11,*)'Cu-Cu A = ',ACU
      WRITE(11,*)'Pb-Pb ksi = ',UPB
      WRITE(11,*)'Pb-Pb p = ',PPB
      WRITE(11,*)'Pb-Pb q = ',QPB
      WRITE(11,*)'Pb-Pb A = ',APB
      WRITE(11,*)'Pb-Cu ksi = ',UPBCU
      CLOSE(11)
C 
      NAME = FNAME//'\VD.TXT'
      OPEN(13,FILE=NAME)                           
       WRITE(13, '(I5, F10.4)')(I,VD(I),I=1,NATALL)
      CLOSE(13)

C	   SCREEN WRITE NANOWIRE  ENERGY
         WRITE(*,*)'Energy per adatom =  ', EAVERS/NADAT

       NAME = FNAME//'\SDF.DAT'
       OPEN(21,FILE=NAME)
        WRITE(21,219)(I,SDF(I),I=1,NZD)
       CLOSE(21)
 219   FORMAT(I4,2X,I9)
C                              
       NAME = FNAME//'\PDF.DAT'
       OPEN(22,FILE=NAME)
        DO J = 1, NML1
          WRITE(22,219) (I,PDF(J,I),I=1,4000)
        END DO
       CLOSE(22)
      END DO  
c      STOP 
      END
C       ***********************************************************
C     -----------------------------------------------------------------
C
      SUBROUTINE SUN(ISEED,RANF,N)
C
C     Stores N real random numbers in RANF, ISEED serves as
C     startup value.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER ISEED
      REAL RANF(N)
        NUL = 0
C
C                            --------------------
C
      DATA FACTOR /41475557.0D0/, TWO28 /268435456.0D0/, ONE/1.0D0/
C
C                            --------------------
C
      IF (ISEED .GE. NUL) THEN
      R=DBLE(ISEED)/TWO28
      R=DMOD(R*FACTOR,ONE)
      ISEED=-1
      ENDIF
C
      DO 100 I = 1, N
      R=DMOD(R*FACTOR,ONE)
      RANF(I) = SNGL(R)
 100  CONTINUE
C
      RETURN
      END SUBROUTINE SUN
