      MODULE GLOBAL
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NP=4000,NE=4000,NT=1000)
      PARAMETER (NX=80,NY=100)
      PARAMETER (NF=1,NS=3)
      COMMON/ABC/ID(16),MP,ME,MBC,MX,MY,MMX,MMY, IER
      COMMON/BBC/X(NP),Y(NP),JNP(NP),IJK(4,NE),LIST(2,NP)
      COMMON/CBC/KBC(2,NP),LBC(2,NP),PA,PB
      COMMON/DBC/DTIME,TIME,IPRINT(NT),MTIME,ITIME,NPRINT
      COMMON/EBC/VEL(2,NX,NY),GYRA(NX,NY)
      COMMON/FBC/PRE(NX,NY),TTT
      COMMON/GBC/DEN,ERTIA,ETAS,GRA,AGAMA,ALAMDA,AMU,AKAPA
      COMMON/IBC/STIF(2*NP,2*NP),BBB(2*NP)
      COMMON/LBC/CAUCHY(NP,2,2),GSTRESS(NP,3,3)
      COMMON/NBC/IIJJKK(4,NF*NS*NE),XCOOR(NF*NS*NX,NF*NS*NY)&
     &,YCOOR(NF*NS*NX,NF*NS*NY)
      END MODULE
!=======================================================================
!=============================================================
      PROGRAM MICROPOLAR
      USE GLOBAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Program Developer: James Chen                          !
!      This is a program for micropolar fluid                 !
!      It is solved by finite difference method               !
!      The numerical scheme is based on Dr. Fu's              !
!      PhD dissertation in 2008.                              !
!      It is two dimension.                                   !
!      It is unsteady and incompressible                      !
!                                                             !
!      START DATE: April 19, 2010                             !
!                                                             !
!      Update:                                                !
!      Sep. 14, 2010: New BC can be used.                     !
!      Sep. 16, 2010: Pressure Two Way Coupling               !
!                     Pressure two-level correction           !
!      Sep. 17, 2010: Inviscid Flow done                      !
!      Sep. 21, 2010: Cauchy Stress                           !
!      Sep. 21, 2010: Total Velocity Post-process             !
!      Sep. 27, 2010: Vorticity Post-process                  !
!      Oct.  4, 2010: Upwind Scheme                           !
!      Oct.  8, 2010: Projection method for pressure          !
!      Nov.  1, 2010: Program Diet                            !
!      Nov.  5, 2010: General Total Velocity Plot Algorithm   !
!      Dec. 28, 2010: SOR Method for Linear Equations         !
!                                                             !
!      User's Manual:                                         !
!      Only Boundary Conditions are needed for velocity and   !
!      gyration.                                              !
!                                                             !
!      Velocity:                                              !
!      RES1 and STIFFNESS1 for Level 1                        !
!      RES2 and STIFFNESS2 for Level 2                        !
!                                                             !
!      Gyration:                                              !
!      GYRA1 and GSTIFF1 for Level 1                          !
!      GYRA2 and GSTIFF2 for level 2                          !
!                                                             !
!      BC for pressure are controlled by index in CONTROL.    !
!                                                             !
!      VERSION: b5.01                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!====================================================
      OPEN(5, FILE='infile', STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(6, FILE='outfile', STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(7, FILE='tecplot.dat', STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(8, FILE='total.dat', STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN(9, FILE='log.dat', STATUS='UNKNOWN',FORM='FORMATTED')
!====================================================
1     FORMAT(10X,'SOMETHING IS WRONG IN INPUT'/)
2     FORMAT(10X,'SOMETHING IS WRONG IN START'/)
3     FORMAT(10X,'SOMETHING IS WRONG IN INITIAL'/)
4     FORMAT(10X,'SOMETHING IS WRONG IN SOLVE'/)
5     FORMAT(10X,'SOMETHING IS WRONG IN OUTPUT'/)       
!====================================================
100   CONTINUE
      CALL INPUT
!     IER=1
      IF(IER.EQ.0) GO TO 200
      WRITE(6,1)
      GO TO 9000
!==================================================
200   CONTINUE
      CALL INITIAL
!     IER=1
      IF(IER.EQ.0) GO TO 300
      WRITE(6,3)
      GO TO 9000
!====================================================
300   CONTINUE
!====================================================
      VEL=0.0D0
      GYRA=0.0D0
      PRE=0.0D0
      CALL RESTART_READ
!====================================================
      TTT=0.0D0
      DO 1000 ITIME=1,MTIME
!     DO 1000 ITIME=1,5
      TTT=TTT+1.0D0
!====================================================
!====================================================
      IF (ID(8).EQ.0.AND.ID(9).EQ.0) GO TO 9000
      IF (ID(9).EQ.0) GO TO 1111
      CALL RES1
      CALL STIFFNESS1
      CALL SOLVE1
      CALL PREC1
      IF (ID(3).EQ.0) GO TO 1111
      CALL GYRA1
      CALL GSTIFF1
      CALL GSOLVE1
1111  CONTINUE
      IF (ID(8).EQ.0) GO TO 2222
      CALL RES2
      CALL STIFFNESS2
      CALL SOLVE2
      CALL PREC2
      IF (ID(3).EQ.0) GO TO 2222
      CALL GYRA2
      CALL GSTIFF2
      CALL GSOLVE2
2222  CONTINUE
      CALL OUTPUT
!====================================================
!====================================================
!====================================================
!====================================================
1000  CONTINUE
!====================================================
!====================================================
      GO TO 9001
9000  CONTINUE
9192  FORMAT('None of the schemes are used.')
      WRITE(6,9192)
9001  CONTINUE
      CLOSE (5)
      CLOSE (6)
      CLOSE (7)
      CLOSE (8)
      CLOSE (9)
!===================================================
      END
!====================================================
      SUBROUTINE INPUT
      USE GLOBAL
!====================================================
100   CONTINUE
      CALL CONTROL
!     IER=1
      IF(IER.EQ.0) GO TO 200
      GO TO 999
!===================================================
200   CONTINUE
      CALL INFO
!     IER=1
      IF(IER.EQ.0) GO TO 300
      GO TO 999
!====================================================
300   CONTINUE
      CALL DYNA
!     IER=1
      IF(IER.EQ.0) GO TO 400
      GO TO 999
!====================================================
400   CONTINUE
      CALL PRINTING
!     IER=1
      IF(IER.EQ.0) GO TO 500
      GO TO 999
!===================================================
500   CONTINUE
      CALL BCIC
!     IER=1
      IF(IER.EQ.0) GO TO 600
      GO TO 999
!===================================================
600   CONTINUE
      CALL MATERIAL
!     IER=1
      GO TO 999
!===================================================
!====================================================
800   CONTINUE
!====================================================
999   CONTINUE
      RETURN
      END
!====================================================
!====================================================
      SUBROUTINE CONTROL
      USE GLOBAL
!====================================================
1     FORMAT(16I9)
!===================================================
      READ(5,1) (ID(I),I=1,10)
!====================================================
!                                                   !
!      ID(1)=0 : Steady                             !
!            1 : Dynamic                            !
!                                                   !
!      ID(2)=0 : Invisicid                          !
!            1 : Viscostic                          !
!                                                   !
!      ID(3)=0 : Navier Stokes Solver               !
!            1 : Micropolar Solver                  !
!                                                   !
!      ID(4)=0 : Pressure (Homogeneous) at bottom   !
!            1 : Pressure Specified at bottom       !
!                                                   !
!      ID(5)=0 : Pressure (Homogeneous) at top      !
!            1 : Pressure Specified at top          !
!                                                   !
!      ID(6)=0 : Pressure (Homogenous) at LH        !
!            1 : Pressure Specified at LH           !
!                                                   !
!      ID(7)=0 : Pressure (Homogeneous) at RH       !
!            1 : Pressure Specified at RH           !
!                                                   !
!      ID(8)=0 : Full Implicit Scheme               !
!            1 : Semi Implicit Scheme               !
!                                                   !
!      ID(9)=0 : Full Explicit Scheme               !
!            1 : Semi Implicit Scheme               !
!                                                   !
!      ID(10)=0 : Direct Method for Linear System   !                                      
!             1 : Iterative Method for Linear System!
!                                                   !
!===================================================!
      WRITE(6,1) (ID(I),I=1,10)
      RETURN
      END
!====================================================
!====================================================
      SUBROUTINE INFO
      USE GLOBAL
!====================================================
1     FORMAT(16I9)
2     FORMAT(10X/&
     &10X,'Number of Nodes in X-direction:',I5/&
     &10X,'Number of Nodes in Y-direction:',I5/&
     &10X,'Total Number of Nodes: ',I5/&
     &10X,'Total Number of Finite Elements: ',I5/)
3     FORMAT(/'Number of Node',&
     &10X,'Coordinate'/)
4     FORMAT(I9,3F15.8,I7)
5     FORMAT(/'Number of Element',&
     &10X,'Connectivity'/)
6     FORMAT(10I5)
!===================================================
      READ(5,1) MX,MY,MP,ME
      WRITE(6,2) MX,MY,MP,ME 
!===================================================
      L=0
      WRITE(6,3)
      DO 100 J=1,MY
      DO 100 K=1,MX
      L=L+1
      LIST(1,L)=K
      LIST(2,L)=J
      READ(5,4) L,X(L),Y(L)
      WRITE(6,4) L,X(L),Y(L)
100   CONTINUE
      IF (L.EQ.MP) GO TO 101
      IER=1
101   CONTINUE
!===================================================
      L=0
      WRITE(6,5)
      DO 200 I=1,ME
      READ(5,1) L,(IJK(J,I),J=1,4)
      WRITE(6,1) L,(IJK(J,I),J=1,4)
200   CONTINUE
      RETURN
      END
!====================================================
!====================================================
      SUBROUTINE DYNA
      USE GLOBAL
!====================================================
1     FORMAT(I9,F15.8)
2     FORMAT(/&
     &10X,'Total Time Step:',I9/&
     &10X,'One Time Step:',F19.9/)
!===================================================
      READ(5,1) MTIME,TIME
      IF (ID(8).EQ.0.OR.ID(9).EQ.0) GO TO 100
      GO TO 200
100   CONTINUE
      WRITE(6,2) MTIME,TIME/2.0D0
      GO TO 300
200   CONTINUE
      WRITE(6,2) MTIME,TIME
!====================================================
300   CONTINUE
      RETURN
      END
!====================================================
!====================================================
!====================================================
      SUBROUTINE PRINTING
      USE GLOBAL
!====================================================
1     FORMAT(I9,F15.8)
2     FORMAT(10I9)
3     FORMAT(10X,'PRINTING FRAME'/)
4     FORMAT(10X,'IRPINT(',I7,')=',I9/)
5     FORMAT(10X,'It is a STEADY case.&
     & Only one frame will be printed.'/)
6     FORMAT(10X,'It is a DYNAMIC case'&
     &'All of the frames will be printed.'/)
!===================================================
      READ(5,1) NPRINT
      DO 521 I=1,NPRINT/10
      READ(5,2) (IPRINT((I-1)*10+J),J=1,10)
521   CONTINUE
      WRITE(6,3)
      DO 100 J=1,NPRINT
      WRITE(6,4) J,IPRINT(J) 
100   CONTINUE       
!====================================================
!     IF (ID(1).EQ.0) WRITE(6,5)
      IF (ID(1).EQ.1) WRITE(6,6)
!====================================================
!====================================================
      RETURN
      END
!====================================================
!====================================================
      SUBROUTINE BCIC
      USE GLOBAL
!====================================================
1     FORMAT(3I9,4F15.8,2I9)
11    FORMAT(3I5,4F15.8,2I9/)
2     FORMAT(10I5)
5     FORMAT(1X,'There are ',I5 ' boundary conditions.&
     &INDEX=0: GRAD(u)=0, INDEX=1: u=0.'/)
6     FORMAT(1X,'Number',3X,'Node Number',3X,'IC'/)
7     FORMAT(16F15.8)
8     FORMAT(1X'The Pressure at both sides are',2F15.8/)
!===================================================
      WRITE(6,6)
      L=0
      DO 100 J=1,MY
      DO 100 K=1,MX
      L=L+1
      READ(5,1) L,(LBC(I,L),I=1,2),(VEL(I,K,J),I=1,2),GYRA(K,J),&
     &PRE(K,J),(KBC(I,L),I=1,2)
      WRITE(6,11) L,(LBC(I,L),I=1,2),(VEL(I,K,J),I=1,2),GYRA(K,J),&
     &PRE(K,J),(KBC(I,L),I=1,2)
100   CONTINUE       
!====================================================
200   CONTINUE
!====================================================
      RETURN
      END
!=====================================================
!====================================================
!====================================================
      SUBROUTINE MATERIAL
      USE GLOBAL
!====================================================
!====================================================
1     FORMAT(16F15.8)
2     FORMAT(2X,'DENSITY=',F15.8,/,2X,'RADIUS=',F15.8, &
     &2X,'MICROINERTIA=',F15.8,/,2X,'GRAVITY=',F15.8, &
     &2X,'GAMA= ', F15.8,/,2X'LAMDA= ', F15.8, &
     &2X'MU= ', F15.8,2X,'KAPPA= ',F15.8/)
3     FORMAT(1X,'The follwoing are material constants.'/)
!====================================================
      WRITE(6,3)
      READ(5,1) DEN,ERTIA,ETAS,GRA,AGAMA,ALAMDA,AMU,AKAPA 
      WRITE(6,2) DEN,ERTIA,ETAS,GRA,AGAMA,ALAMDA,AMU,AKAPA
!====================================================
!====================================================
!====================================================
      RETURN 
      END
!====================================================
!====================================================
      SUBROUTINE PLOT1
      USE GLOBAL
!====================================================
!===================================================
1     FORMAT('Title="Micropolar Fluid Dynamics"')
2     FORMAT('Variables="X","Y","VEL_X","VEL_Y","PRESSURE","GYRA_Z",&
     &"CSTRESS_XX","CSTRESS_XY","CSTRESS_YX","CSTRESS_YY","GSTRESS_XZ",&
     &"GSTRESS_YZ","GSTRESS_ZX","GSTRESS_ZY"')
3     FORMAT('Zone T="Load Step 0",N=',I5,' E=',I5,' datapacking=point,&
     &zonetype=fequadrilateral')
13    FORMAT('Zone T="Load Step 0",N=',I5,' E=',I5,'')
4     FORMAT(25F15.8)
5     FORMAT(8I10)
6     FORMAT('It is an inviscid flow.')
7     FORMAT('Viscosity is considered.')
8     FORMAT('It is running as a Navier-Stokes Solver')
9     FORMAT('It is running as a Micropolar Solver')
!====================================================
      IF (ID(2).EQ.0) WRITE(6,6)
      IF (ID(2).EQ.1) WRITE(6,7)
      IF (ID(3).EQ.0) WRITE(6,8)
      IF (ID(3).EQ.1) WRITE(6,9)
      WRITE(7,1)
      WRITE(7,2)
      WRITE(7,3) MP,ME
      WRITE(6,13) MP,ME
!====================================================
      DO 100 I=1,MP
      WRITE(7,4) X(I),Y(I),&
     &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
     &PRE(LIST(1,I),LIST(2,I)),&
     &(GYRA(LIST(1,I),LIST(2,I))),&
     &CAUCHY(I,1,1),CAUCHY(I,1,2),CAUCHY(I,2,1),&
     &CAUCHY(I,2,2),GSTRESS(I,1,3),GSTRESS(I,2,3),&
     &GSTRESS(I,3,1),GSTRESS(I,3,2)
      WRITE(6,4) X(I),Y(I),&
     &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
     &PRE(LIST(1,I),LIST(2,I)),&
     &(GYRA(LIST(1,I),LIST(2,I)))
100   CONTINUE
!====================================================
      DO 200 J=1,ME
      WRITE(7,5) (IJK(II,J),II=1,4)
200   CONTINUE
!====================================================
!====================================================
      RETURN
      END
!=====================================================
!====================================================
      SUBROUTINE INITIAL
      USE GLOBAL
!=====================================================
1     FORMAT(I5,F15.8)
!====================================================
!=====================================================
      CALL PLOT1
!     CALL OUTPUT1
!     CALL PLOT3
!=====================================================
      RETURN
      END
!=====================================================
!=====================================================
!=====================================================
      SUBROUTINE RES1
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
!=====================================================
      L=0
      BBB=0.0D0
      IF (ID(3).EQ.0) AKAPA=0.0D0
      DO 199 I=1,MP
      DO 199 J=1,2
!=====================================================
!     L=(I-1)*2+J
      L=L+1
!     IF (KBC(1,I).GT.0.1D0) GO TO 106
      IF (LIST(2,I).EQ.1) GO TO 103
      IF (LIST(2,I).EQ.MY) GO TO 109
      IF (LIST(1,I).EQ.1) GO TO 104
      IF (LIST(1,I).EQ.MX) GO TO 105
!     IF (LIST(1,I).LE.20 .AND. LIST(1,I).GE.10 .AND.&
!    & LIST(2,I).LE.12 .AND. LIST(2,I).GE.9) GO TO 106
      IF (J.EQ.1) GO TO 101
      IF (J.EQ.2) GO TO 102
!=====================================================
101   CONTINUE
      BBB(L)=VEL(1,LIST(1,I),LIST(2,I))
      BBB(L)=BBB(L)+(2.5D-1*TIME*GRA/DEN*0.5D0*DSQRT(2.0D0))
      IF (ID(2).EQ.0) GO TO 111
      IF (ID(3).EQ.0) GO TO 111
      BBB(L)=BBB(L)+(2.5D-1*TIME*AKAPA/DEN*&
     &(GYRA(LIST(1,I),LIST(2,I)+1)-GYRA(LIST(1,I),LIST(2,I)-1))&
     &/(Y(I+MX)-Y(I-MX)))
111  CONTINUE
      GO TO 100
!=====================================================
102   CONTINUE
      BBB(L)=VEL(2,LIST(1,I),LIST(2,I))
      BBB(L)=BBB(L)-(2.5D-1*TIME*GRA/DEN*0.5D0*DSQRT(2.0D0))
      IF (ID(2).EQ.0) GO TO 112
      IF (ID(3).EQ.0) GO TO 112
      BBB(L)=BBB(L)-(2.5D-1*TIME*AKAPA/DEN*&
     &(GYRA(LIST(1,I)+1,LIST(2,I))-GYRA(LIST(1,I)-1,LIST(2,I)))&
     &/(X(I+1)-X(I-1)))
112   CONTINUE
      GO TO 100
!=====================================================
103   CONTINUE
      IF (J.EQ.1) GO TO 123
      IF (J.EQ.2) GO TO 133
123   CONTINUE
      BBB(L)=0.0D0
      BBB(L)=-1.0D0
      GO TO 100
133   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
109   CONTINUE
      IF (J.EQ.1) GO TO 129
      IF (J.EQ.2) GO TO 139
129   CONTINUE
      BBB(L)=0.0D0
      BBB(L)=1.0D0 
      GO TO 100
139   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
104   CONTINUE
      IF (J.EQ.1) GO TO 124
      IF (J.EQ.2) GO TO 134
124   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
134   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
!=====================================================
105   CONTINUE
      IF (J.EQ.1) GO TO 125
      IF (J.EQ.2) GO TO 135
125   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
135   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
!=====================================================
!=====================================================
106   CONTINUE
      IF (J.EQ.1) GO TO 126
      IF (J.EQ.2) GO TO 136
126   CONTINUE
      BBB(L)=5.0D0
!     WRITE(6,1) L,I,J,BBB(L)
      GO TO 100
136   CONTINUE
      IF (KBC(2,I).GT.0) THEN
      BBB(L)=-5.0D0
      ELSE
      BBB(L)=5.0D0
      END IF
!     WRITE(6,1) L,I,J,BBB(L)
      GO TO 100
!=====================================================
100   CONTINUE
!     WRITE(6,1) L,I,J,BBB1(L)
1     FORMAT(3I5,F15.8)
199   CONTINUE
      RETURN
      END
!=====================================================
!=====================================================
!=====================================================
      SUBROUTINE STIFFNESS1
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
!=====================================================
      L=0
      STIF=0.0D0
      IF (ID(3).EQ.0) AKAPA=0.0D0
      DO 100 I=1,MP
      DO 100 J=1,2
!=====================================================
!     L=(I-1)*2+J
      L=L+1
!     IF (KBC(1,I).GT.0.1D0) GO TO 106
      IF (LIST(2,I).EQ.1) GO TO 103
      IF (LIST(2,I).EQ.MY) GO TO 109
      IF (LIST(1,I).EQ.1) GO TO 104
      IF (LIST(1,I).EQ.MX) GO TO 105
!     IF (LIST(1,I).LE.20 .AND. LIST(1,I).GE.10 .AND.&
!    & LIST(2,I).LE.12 .AND. LIST(2,I).GE.9) GO TO 106
      IF (J.EQ.1) GO TO 101
      IF (J.EQ.2) GO TO 101
!=====================================================
!=====================================================
101   CONTINUE
      STIF(L,L)=1.0D0
      IF (VEL(1,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1011
      IF (VEL(1,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1012
1011  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME*VEL(1,LIST(1,I),LIST(2,I))/&
     &(X(I)-X(I-1)))
      STIF(L-2,L)=-2.5D-1*TIME*VEL(1,LIST(1,I),LIST(2,I))/(X(I)-X(I-1))
      GO TO 1013
1012  CONTINUE
      STIF(L,L)=STIF(L,L)-(2.5D-1*TIME*VEL(1,LIST(1,I),LIST(2,I))/&
     &(X(I)-X(I-1)))
      STIF(L+2,L)=2.5D-1*TIME*VEL(1,LIST(1,I),LIST(2,I))/(X(I)-X(I-1))
      GO TO 1013
!=====================================================
1013  CONTINUE
      IF (VEL(2,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1014
      IF (VEL(2,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1015
1014  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME*VEL(2,LIST(1,I),LIST(2,I))/&
     &(Y(I)-Y(I-MX)))
      STIF(L-(2*MX),L)=(-2.5D-1*TIME*VEL(2,LIST(1,I),LIST(2,I))&
     &/(Y(I)-Y(I-MX)))
      GO TO 1016
1015  CONTINUE
      STIF(L,L)=STIF(L,L)-(2.5D-1*TIME*VEL(2,LIST(1,I),LIST(2,I))/&
     &(Y(I)-Y(I-MX)))
      STIF(L+(2*MX),L)=(2.5D-1*TIME*VEL(2,LIST(1,I),LIST(2,I))&
     &/(Y(I)-Y(I-MX)))
      GO TO 1016
!=====================================================
1016  CONTINUE
!=====================================================
      IF (ID(2).EQ.0) GO TO 100
      STIF(L,L)=STIF(L,L)+(0.5D0*TIME*(AMU+AKAPA)/DEN*&
     &((1.0D0/((X(I)-X(I-1))**2.0D0))+(1.0D0/((Y(I)-Y(I-MX))**2.0D0))))
!=====================================================
      STIF(L+2,L)=STIF(L+2,L)-(2.5D-1*TIME*(AMU+AKAPA)/DEN/&
     &((X(I)-X(I-1))**2.0D0))
!=====================================================
      STIF(L+(2*MX),L)=STIF(L+(2*MX),L)-(2.5D-1*TIME*(AMU+AKAPA)/DEN/&
     &((Y(I)-Y(I-MX))**2.0D0))
!=====================================================
      STIF(L-2,L)=STIF(L-2,L)-(2.5D-1*TIME*(AMU+AKAPA)/DEN/&
     &((X(I)-X(I-1))**2.0D0))
!=====================================================
      STIF(L-(2*MX),L)=STIF(L-(2*MX),L)-(2.5D-1*TIME*(AMU+AKAPA)/DEN/&
     &((Y(I)-Y(I-MX))**2.0D0))
!=====================================================
      GO TO 100
!=====================================================
!=====================================================
103   CONTINUE
      IF (J.EQ.1) GO TO 113
      IF (J.EQ.2) GO TO 123
113   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+(2*MX),L)=-1.0D0
      GO TO 100
123   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+(2*MX),L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
109   CONTINUE
      IF (J.EQ.1) GO TO 119
      IF (J.EQ.2) GO TO 129
119   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-(2*MX),L)=-1.0D0
      GO TO 100
129   CONTINUE
!     IF (LIST(1,I).EQ.1) GO TO 1291
!     GO TO 1292
1291  CONTINUE
!     STIF(L,L)=1.0D0
!     STIF(L-2,L)=-1.0D0
!     GO TO 100
1292  CONTINUE
!     STIF(L,L)=1.0D0
!     STIF(L+2,L)=-1.0D0
      STIF(L,L)=1.0D0
!     STIF(L-(2*MX),L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
104   CONTINUE
      IF (J.EQ.1) GO TO 114 
      IF (J.EQ.2) GO TO 124 
114   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+2,L)=-1.0D0
!     STIF(L-(2*MX)+2,L)=-1.0D0
      GO TO 100
124   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+2,L)=-1.0D0
!     STIF(L-(2*MX)+2,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
!=====================================================
105   CONTINUE
      IF (J.EQ.1) GO TO 115 
      IF (J.EQ.2) GO TO 125 
115   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-2,L)=-1.0D0
      GO TO 100
125   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-2,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
!=====================================================
106   CONTINUE
      IF (J.EQ.1) GO TO 116
      IF (J.EQ.2) GO TO 126
116   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
126   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
!=====================================================
!=====================================================
100   CONTINUE
1     FORMAT(I5,50F15.8)
!     DO 199 I=1,L
!     WRITE(6,1) I,STIF1(I,I),STIF1(I+2,I),STIF1(I+(2*MX),I)&
!    &,STIF1(I-2,I),STIF1(I-(2*MX),I)
199   CONTINUE
      RETURN
      END
!=====================================================
!=====================================================
      SUBROUTINE STIFFNESS2
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
!=====================================================
      L=0
      STIF=0.0D0
      IF (ID(3).EQ.0) AKAPA=0.0D0
      DO 100 I=1,MP
      DO 100 J=1,2
!=====================================================
      L=L+1
!     IF (KBC(1,I).GT.0.1D0) GO TO 106
      IF (LIST(2,I).EQ.1) GO TO 103
      IF (LIST(2,I).EQ.MY) GO TO 109
      IF (LIST(1,I).EQ.1) GO TO 104
      IF (LIST(1,I).EQ.MX) GO TO 105
!     IF (LIST(1,I).LE.20 .AND. LIST(1,I).GE.10 .AND.&
!    & LIST(2,I).LE.12 .AND. LIST(2,I).GE.9) GO TO 106
      IF (J.EQ.1) GO TO 101
      IF (J.EQ.2) GO TO 102
!=====================================================
!=====================================================
101   CONTINUE
      STIF(L,L)=1.0D0
      IF (VEL(1,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1011
      IF (VEL(1,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1012
1011  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME/(X(I)-X(I-1))*&
     &(VEL(1,LIST(1,I),LIST(2,I))-VEL(1,LIST(1,I)-1,LIST(2,I))))
      GO TO 1013
1012  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME/(X(I+1)-X(I))*&
     &(VEL(1,LIST(1,I)+1,LIST(2,I))-VEL(1,LIST(1,I),LIST(2,I))))
      GO TO 1013
1013  CONTINUE
      IF (VEL(2,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1014
      IF (VEL(2,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1015
1014  CONTINUE
      STIF(L+1,L)=(2.5D-1*TIME/(Y(I)-Y(I-MX))*&
     &(VEL(1,LIST(1,I),LIST(2,I))-VEL(1,LIST(1,I),LIST(2,I)-1)))
      GO TO 100
1015  CONTINUE
      STIF(L+1,L)=(2.5D-1*TIME/(Y(I+MX)-Y(I))*&
     &(VEL(1,LIST(1,I),LIST(2,I)+1)-VEL(1,LIST(1,I),LIST(2,I))))
      GO TO 100
!=====================================================
!=====================================================
102   CONTINUE
      STIF(L,L)=1.0D0
      IF (VEL(2,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1021
      IF (VEL(2,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1022
1021  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME/(Y(I)-Y(I-MX))*&
     &(VEL(2,LIST(1,I),LIST(2,I))-VEL(2,LIST(1,I),LIST(2,I)-1)))
      GO TO 1023
1022  CONTINUE
      STIF(L,L)=STIF(L,L)+(2.5D-1*TIME/(Y(I+MX)-Y(I))*&
     &(VEL(2,LIST(1,I),LIST(2,I)+1)-VEL(2,LIST(1,I),LIST(2,I))))
      GO TO 1023
!=====================================================
1023  CONTINUE
      IF (VEL(1,LIST(1,I),LIST(2,I)).GE.0.0D0) GO TO 1024
      IF (VEL(1,LIST(1,I),LIST(2,I)).LE.0.0D0) GO TO 1025
1024  CONTINUE
      STIF(L-1,L)=(2.5D-1*TIME/(X(I)-X(I-1))*&
     &(VEL(2,LIST(1,I),LIST(2,I))-VEL(2,LIST(1,I)-1,LIST(2,I))))
      GO TO 100
1025  CONTINUE
      STIF(L-1,L)=(2.5D-1*TIME/(X(I+1)-X(I))*&
     &(VEL(2,LIST(1,I)+1,LIST(2,I))-VEL(2,LIST(1,I),LIST(2,I))))
      GO TO 100
!=====================================================
!=====================================================
103   CONTINUE
      IF (J.EQ.1) GO TO 113
      IF (J.EQ.2) GO TO 123
113   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+(2*MX),L)=-1.0D0
      GO TO 100
123   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+(2*MX),L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
109   CONTINUE
      IF (J.EQ.1) GO TO 119
      IF (J.EQ.2) GO TO 129
119   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-(2*MX),L)=-1.0D0
      GO TO 100
129   CONTINUE
!     IF (LIST(1,I).EQ.1) GO TO 1291
!     GO TO 1292
1291  CONTINUE
!     STIF(L,L)=1.0D0
!     STIF(L-2,L)=-1.0D0
!     GO TO 100
1292  CONTINUE
!     STIF(L,L)=1.0D0
!     STIF(L+2,L)=-1.0D0
      STIF(L,L)=1.0D0
!     STIF(L-(2*MX),L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
104   CONTINUE
      IF (J.EQ.1) GO TO 114
      IF (J.EQ.2) GO TO 124
114   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+2,L)=-1.0D0
!     STIF(L-(2*MX)+2,L)=-1.0D0
      GO TO 100
124   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+2,L)=-1.0D0
!     STIF(L-(2*MX)+2,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
105   CONTINUE
      IF (J.EQ.1) GO TO 115
      IF (J.EQ.2) GO TO 125
115   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-2,L)=-1.0D0
      GO TO 100
125   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-2,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
106   CONTINUE
      IF (J.EQ.1) GO TO 116
      IF (J.EQ.2) GO TO 126
116   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
126   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
!=====================================================
!=====================================================
!     WRITE(6,1) L 
100   CONTINUE
1     FORMAT(50F15.8)
!=====================================================
!     DO 199 I=1,2*MP
!     WRITE(6,1) (STIF1(I,JJJ),JJJ=I,I+3)
199   CONTINUE
!=====================================================
      RETURN
      END
!=====================================================
!=====================================================
      SUBROUTINE GYRA1
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
      L=0
      BBB=0.0D0
!=====================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!=====================================================
      L=L+1
      TUR=-1.0D0
      IF (J.EQ.1 .AND. I.EQ.1) GO TO 1011
      IF (J.EQ.1 .AND. I.EQ.MY) GO TO 1091
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 109
      IF (J.EQ.1) GO TO 102
      IF (J.EQ.MX) GO TO 103
!     IF (J.LE.20 .AND. J.GE.10 .AND.&
!    & I.LE.12 .AND. I.GE.9) GO TO 106
      GO TO 104
!=====================================================
109   CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J,I)-VEL(2,J-1,I))/(x(L)-X(L-1)))-&
     &((VEL(1,J,I)-VEL(1,J,I-1))/(Y(L)-Y(L-MX))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
1091  CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J+1,I)-VEL(2,J,I))/(x(L+1)-X(L)))-&
     &((VEL(1,J,I)-VEL(1,J,I-1))/(Y(L)-Y(L-MX))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
101   CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J,I)-VEL(2,J-1,I))/(x(L)-X(L-1)))-&
     &((VEL(1,J,I+1)-VEL(1,J,I))/(Y(L+MX)-Y(L))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
1011  CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J+1,I)-VEL(2,J,I))/(x(L+1)-X(L)))-&
     &((VEL(1,J,I+1)-VEL(1,J,I))/(Y(L+MX)-Y(L))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
102   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
103   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
104   CONTINUE
      BBB(L)=GYRA(J,I)
      IF (ID(2).EQ.0) GO TO 100
      BBB(L)=BBB(L)+(0.5D0*TIME*AKAPA/DEN/ETAS*&
     &(((VEL(2,J+1,I)-VEL(2,J-1,I))/((X(L+1)-X(L-1))))-&
     &((VEL(1,J,I+1)-VEL(1,J,I-1))/((Y(L+MX)-Y(L-MX))))))
!=====================================================
      GO TO 100
!=====================================================
!=====================================================
106   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
5     FORMAT(15F15.8)
100   CONTINUE
!     DO 999 I=1,MP
!     WRITE(6,5) X(I),Y(I),GBB(I),GYRA(LIST(1,I),LIST(2,I)),&
!    &(VEL(J,LIST(1,I),LIST(2,I)),J=1,2)
999   CONTINUE
!     IF (L.EQ.MP) WRITE(6,*) 'I am in GYRA1 and always right.'
      RETURN
      END
!======================================================
      SUBROUTINE SOLVE1
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(2*NP,2*NP),BB(2*NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=2*MP
      NNPP=2*NP
      MM=1
      MMPP=1
!======================================================
      BB=0.0D0
!======================================================
      DO 100 I=1,2*MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
      AA=0.0D0
!======================================================
      DO 200 I=1,2*MP
      DO 200 J=1,2*MP
!======================================================
      AA(J,I)=STIF(I,J)
!     WRITE(6,*) I,J,STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF (ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF (ID(10).EQ.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      L=0
      DO 300 I=1,MP
      DO 300 J=1,2
      L=L+1
      VEL(J,LIST(1,I),LIST(2,I))=BB(L,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
!======================================================
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=2000)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
!     WRITE(6,*) 'I am here.'
!=====================================================
!     DO I=1,N
!     WRITE(6,*) I,a(I,I),b(I,1)
!     END DO
!=====================================================
      do j=1,n
      ipiv(j)=0
      enddo
      do i=1,n 
      big=0
      do j=1,n 
      if(ipiv(j).ne.1)then
      do k=1,n
      if (ipiv(k).eq.0) then
      if (abs(a(j,k)).ge.big)then
      big=abs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then 
      WRITE(6,*) 'singular matrix in gaussj'
      WRITE(6,*) i,j,a(i,j)
      endif
      enddo
      endif
      enddo
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
      enddo
      do l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
      enddo
      endif
      indxr(i)=irow 
      indxc(i)=icol 
      if (a(icol,icol).eq.0.) WRITE(6,*) 'singular matrix in gaussj'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
      a(icol,l)=a(icol,l)*pivinv
      enddo
      do l=1,m
      b(icol,l)=b(icol,l)*pivinv
      enddo
      do ll=1,n
      if(ll.ne.icol)then 
      dum=a(ll,icol)
      a(ll,icol)=0
      do l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
      enddo
      do l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
      enddo
      endif
      enddo
      enddo
      do l=n,1,-1
      if(indxr(l).ne.indxc(l))then
      do k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
      enddo
      endif
      enddo
      return
      END
!=====================================================
!=====================================================
      SUBROUTINE RES2
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
!=====================================================
      L=0
      IF (ID(3).EQ.0) AKAPA=0.0D0
      BBB=0.0D0
      DO 199 I=1,MP
      DO 199 J=1,2
!=====================================================
      L=L+1
!     IF (KBC(1,I).GT.0.1D0) GO TO 106
      IF (LIST(2,I).EQ.1) GO TO 103
      IF (LIST(2,I).EQ.MY) GO TO 109
      IF (LIST(1,I).EQ.1) GO TO 104
      IF (LIST(1,I).EQ.MX) GO TO 105
!     IF (LIST(1,I).LE.20 .AND. LIST(1,I).GE.10 .AND.&
!    & LIST(2,I).LE.12 .AND. LIST(2,I).GE.9) GO TO 106
      IF (J.EQ.1) GO TO 101
      IF (J.EQ.2) GO TO 102
!=====================================================
101   CONTINUE
      BBB(L)=VEL(1,LIST(1,I),LIST(2,I))
      BBB(L)=BBB(L)+(2.5D-1*TIME*GRA/DEN*0.5D0*DSQRT(2.0D0))
      IF (ID(2).EQ.0) GO TO 111
      IF (ID(3).EQ.0) GO TO 1011
      BBB(L)=BBB(L)+(2.5D-1*TIME*AKAPA/DEN*&
     &(GYRA(LIST(1,I),LIST(2,I)+1)-GYRA(LIST(1,I),LIST(2,I)-1))&
     &/((Y(I+MX)-Y(I-MX))))
1011  CONTINUE
      BBB(L)=BBB(L)+(2.5D-1*TIME*(AMU+AKAPA)/DEN*&
     &(((VEL(1,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+&
     &VEL(1,LIST(1,I)-1,LIST(2,I)))/((X(I)-X(I-1))**2.0D0))&
     &+((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+&
     &VEL(1,LIST(1,I),LIST(2,I)-1))/((Y(I)-Y(I-MX))**2.0D0))))
111  CONTINUE
      GO TO 100
!=====================================================
102   CONTINUE
      BBB(L)=VEL(2,LIST(1,I),LIST(2,I))
      BBB(L)=BBB(L)-(2.5D-1*TIME*GRA/DEN*0.5D0*DSQRT(2.0D0))
      IF (ID(2).EQ.0) GO TO 112
      IF (ID(3).EQ.0) GO TO 1021
      BBB(L)=BBB(L)-(2.5D-1*TIME*AKAPA/DEN*&
     &(GYRA(LIST(1,I)+1,LIST(2,I))-GYRA(LIST(1,I)-1,LIST(2,I)))&
     &/((X(I+1)-X(I-1))))
1021  CONTINUE
      BBB(L)=BBB(L)+(2.5D-1*TIME*(AMU+AKAPA)/DEN*&
     &(((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)))/((X(I)-X(I-1))**2.0D0))&
     &+((VEL(2,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I),LIST(2,I)-1))/((Y(I)-Y(I-MX))**2.0D0))))
112   CONTINUE
      GO TO 100
!=====================================================
103   CONTINUE
      IF (J.EQ.1) GO TO 123
      IF (J.EQ.2) GO TO 133
123   CONTINUE
      BBB(L)=0.0D0
      BBB(L)=-1.0D0
      GO TO 100
133   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
109   CONTINUE
      IF (J.EQ.1) GO TO 129
      IF (J.EQ.2) GO TO 139
129   CONTINUE
      BBB(L)=0.0D0
      BBB(L)=1.0D0
      GO TO 100
139   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
104   CONTINUE
      IF (J.EQ.1) GO TO 124
      IF (J.EQ.2) GO TO 134
124   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
134   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
!=====================================================
105   CONTINUE
      IF (J.EQ.1) GO TO 125
      IF (J.EQ.2) GO TO 135
125   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
135   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
106   CONTINUE
      IF (J.EQ.1) GO TO 126
      IF (J.EQ.2) GO TO 136
126   CONTINUE
      BBB(L)=5.0D0
!     WRITE(6,1) L,I,J,BBB(L)
      GO TO 100
136   CONTINUE
      IF (KBC(2,I).GT.0) THEN
      BBB(L)=-5.0D0
      ELSE
      BBB(L)=5.0D0
      END IF
!     WRITE(6,1) L,I,J,BBB(L)
      GO TO 100
!=====================================================
100   CONTINUE
!     WRITE(6,1) L,BBB1(L)
1     FORMAT(3I5,F15.8)
199   CONTINUE
      RETURN
      END
!=====================================================
!=====================================================
!======================================================
      SUBROUTINE SOLVE2
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(2*NP,2*NP),BB(2*NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=2*MP
      NNPP=2*NP
      MM=1
      MMPP=1
!======================================================
!======================================================
      BB=0.0D0
      BA=0.0D0
      DO 100 I=1,2*MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
!======================================================
      AA=0.0D0
      DO 200 I=1,2*MP
      DO 200 J=1,2*MP
!======================================================
      AA(J,I)=STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF(ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF(ID(10).EQ.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      L=0
      DO 300 I=1,MP
      DO 300 J=1,2
      L=L+1
      VEL(J,LIST(1,I),LIST(2,I))=BB(L,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
!=====================================================
!=====================================================
      SUBROUTINE GYRA2
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
      L=0
      BBB=0.0D0
!=====================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!=====================================================
      L=L+1
      TUR=-1.0D0
      IF (J.EQ.1 .AND. I.EQ.1) GO TO 1011
      IF (J.EQ.1 .AND. I.EQ.MY) GO TO 1091
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 109
      IF (J.EQ.1) GO TO 102
      IF (J.EQ.MX) GO TO 103
!     IF (J.LE.20 .AND. J.GE.10 .AND.&
!    & I.LE.12 .AND. I.GE.9) GO TO 106
      GO TO 104
!=====================================================
101   CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J,I)-VEL(2,J-1,I))/(x(L)-X(L-1)))-&
     &((VEL(1,J,I+1)-VEL(1,J,I))/(Y(L+MX)-Y(L))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
1011  CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J+1,I)-VEL(2,J,I))/(x(L+1)-X(L)))-&
     &((VEL(1,J,I+1)-VEL(1,J,I))/(Y(L+MX)-Y(L))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
109   CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J,I)-VEL(2,J-1,I))/(x(L)-X(L-1)))-&
     &((VEL(1,J,I)-VEL(1,J,I-1))/(Y(L)-Y(L-MX))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
1091  CONTINUE
      PLUS=TUR*&
     &(((VEL(2,J+1,I)-VEL(2,J,I))/(x(L+1)-X(L)))-&
     &((VEL(1,J,I)-VEL(1,J,I-1))/(Y(L)-Y(L-MX))))
      BBB(L)=0.0D0+PLUS
!     BBB(L)=0.0D0
      GO TO 100
!=====================================================
102   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
103   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
104   CONTINUE
      BBB(L)=GYRA(J,I)
      BBB(L)=BBB(L)+(0.5D0*TIME*AKAPA/DEN/ETAS*&
     &(((VEL(2,J+1,I)-VEL(2,J-1,I))/((X(L+1)-X(L-1))))-&
     &((VEL(1,J,I+1)-VEL(1,J,I-1))/((Y(L+MX)-Y(L-MX))))))
!=====================================================
      GO TO 100
!=====================================================
106   CONTINUE
      BBB(L)=0.0D0
      GO TO 100
!=====================================================
!=====================================================
5     FORMAT(15F15.8)
100   CONTINUE
!     DO 999 I=1,MP
!     WRITE(6,5) X(I),Y(I),GBB(I),GYRA(LIST(1,I),LIST(2,I)),&
!    &(VEL(J,LIST(1,I),LIST(2,I)),J=1,2)
999   CONTINUE
      RETURN
      END
!======================================================
!====================================================
!====================================================
      SUBROUTINE PLOT2
      USE GLOBAL
!====================================================
!===================================================
1     FORMAT('Title="Micropolar Fluid Dynamics"')
2     FORMAT('Variables="X","Y","VEL_X","VEL_Y","PRESSURE","GYRA_Z",&
     &"CSTRESS_XX","CSTRESS_XY","CSTRESS_YX","CSTRESS_YY","GSTRESS_XZ",&
     &"GSTRESS_YZ","GSTRESS_ZX","GSTRESS_ZY"')
3     FORMAT('Zone T="Load Step',I9'",N=',I5,' E=',I5,' &
     &datapacking=point,&
     &zonetype=fequadrilateral,connectivitysharezone=1')
13    FORMAT('Zone T="Load Step 0",N=',I5,' E=',I5,'')
4     FORMAT(25F15.8)
5     FORMAT(8I10)
!====================================================
      WRITE(7,3) ITIME,MP,ME
!====================================================
      DO 100 I=1,MP
      WRITE(7,4) X(I),Y(I),&
     &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
     &PRE(LIST(1,I),LIST(2,I)),&
     &(GYRA(LIST(1,I),LIST(2,I))),&
     &CAUCHY(I,1,1),CAUCHY(I,1,2),CAUCHY(I,2,1),&
     &CAUCHY(I,2,2),GSTRESS(I,1,3),GSTRESS(I,2,3),&
     &GSTRESS(I,3,1),GSTRESS(I,3,2)
      WRITE(6,4) X(I),Y(I),&
     &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
     &PRE(LIST(1,I),LIST(2,I)),&
     &(GYRA(LIST(1,I),LIST(2,I)))
100   CONTINUE
!====================================================
!     DO 101 I=1,MP
!     WRITE(7,4) X(I),Y(I),&
!    &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
!    &PRE(LIST(1,I),LIST(2,I)),&
!    &(GYRA(LIST(1,I),LIST(2,I)))
!     WRITE(6,4) X(I),Y(I),&
!    &(VEL(LL,LIST(1,I),LIST(2,I)),LL=1,2),&
!    &PRE(LIST(1,I),LIST(2,I)),&
!    &(GYRA(LIST(1,I),LIST(2,I)))
101   CONTINUE
!====================================================
      RETURN
      END
!=====================================================
!====================================================
!=====================================================
!=====================================================
      SUBROUTINE GSTIFF1
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
      L=0
      STIF=0.0D0
!=====================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!=====================================================
      L=L+1
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 109
      IF (J.EQ.1) GO TO 102
      IF (J.EQ.MX) GO TO 103
!     IF (J.LE.20 .AND. J.GE.10 .AND.&
!    & I.LE.12 .AND. I.GE.9) GO TO 106
      GO TO 104
!=====================================================
101   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+MX,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
109   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-MX,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
102   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+1,L)=-1.0D0
!     STIF(L-MX+1,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
103   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-1,L)=-1.0D0
!     STIF(L-MX+1,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
104   CONTINUE
      STIF(L,L)=1.0D0
!=====================================================
      IF (VEL(1,J,I).GE.0.0D0) GO TO 1041
      IF (VEL(1,J,I).LE.0.0D0) GO TO 1042
1041  CONTINUE
      STIF(L,L)=STIF(L,L)+(0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I))
      STIF(L-1,L)=-0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I)
      GO TO 1045
1042  CONTINUE
      STIF(L,L)=STIF(L,L)-(0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I))
      STIF(L+1,L)=0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I)
      GO TO 1045
!=====================================================
1045  CONTINUE
      IF (VEL(2,J,I).GE.0.0D0) GO TO 1043
      IF (VEL(2,J,I).LE.0.0D0) GO TO 1044
1043  CONTINUE
      STIF(L,L)=STIF(L,L)+(0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I))
      STIF(L-MX,L)=-0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I)
      GO TO 1046
1044  CONTINUE
      STIF(L,L)=STIF(L,L)-(0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I))
      STIF(L+MX,L)=0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I)
      GO TO 1046
1046  CONTINUE
!=====================================================
      IF (ID(2).EQ.0) GO TO 100
      STIF(L,L)=STIF(L,L)+&
     &(TIME*AGAMA/DEN/ETAS*&
     &((1.0D0/(X(L)-X(L-1))**2.0D0)+(1.0D0/((Y(L)-Y(L-MX))**2.0D0))))
      STIF(L,L)=STIF(L,L)+(AKAPA*TIME/DEN/ETAS)
!=====================================================
!=====================================================
      STIF(L+1,L)=STIF(L+1,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((X(L)-X(L-1))**2.0D0))
!=====================================================
!=====================================================
      STIF(L+MX,L)=STIF(L+MX,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((Y(L)-Y(L-MX))**2.0D0))
!=====================================================
!=====================================================
      STIF(L-1,L)=STIF(L-1,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((X(L)-X(L-1))**2.0D0))
!=====================================================
!=====================================================
      STIF(L-MX,L)=STIF(L-MX,L)-(0.5D0*TIME*AGAMA/DEN/ETAS/&
     &((Y(L)-Y(L-MX))**2.0D0))
      GO TO 100
!=====================================================
!=====================================================
106   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
5     FORMAT(F15.8)
100   CONTINUE
      RETURN
      END
!======================================================
!======================================================
      SUBROUTINE GSOLVE1
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(NP,NP),BB(NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=MP
      NNPP=NP
      MM=1
      MMPP=1
!======================================================
      BB=0.0D0
!======================================================
      DO 100 I=1,MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
!======================================================
      AA=0.0D0
      DO 200 I=1,MP
      DO 200 J=1,MP
!======================================================
      AA(J,I)=STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF(ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF(ID(10).EQ.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      DO 300 I=1,MP
      GYRA(LIST(1,I),LIST(2,I))=BB(I,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
!=====================================================
!=====================================================
      SUBROUTINE GSTIFF2
      USE GLOBAL
!=====================================================
!     Checked: September 11, 2010                    =
!=====================================================
!=====================================================
      L=0
      STIF=0.0D0
!=====================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!=====================================================
      L=L+1
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 109
      IF (J.EQ.1) GO TO 102
      IF (J.EQ.MX) GO TO 103
!     IF (J.LE.20 .AND. J.GE.10 .AND.&
!    & I.LE.12 .AND. I.GE.9) GO TO 106
      GO TO 104
!=====================================================
101   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+MX,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
109   CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-MX,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
102   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L+1,L)=-1.0D0
!     STIF(L-MX+1,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
!=====================================================
103   CONTINUE
      STIF(L,L)=1.0D0
      STIF(L-1,L)=-1.0D0
!     STIF(L-MX+1,L)=-1.0D0
      GO TO 100
!=====================================================
!=====================================================
104   CONTINUE
      STIF(L,L)=1.0D0
!=====================================================
      IF (VEL(1,J,I).GE.0.0D0) GO TO 1041
      IF (VEL(1,J,I).LE.0.0D0) GO TO 1042
1041  CONTINUE
      STIF(L,L)=STIF(L,L)+(0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I))
      STIF(L-1,L)=-0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I)
      GO TO 1045
1042  CONTINUE
      STIF(L,L)=STIF(L,L)-(0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I))
      STIF(L+1,L)=0.5D0*TIME/(X(L)-X(L-1))*VEL(1,J,I)
      GO TO 1045
!=====================================================
1045  CONTINUE
      IF (VEL(2,J,I).GE.0.0D0) GO TO 1043
      IF (VEL(2,J,I).LE.0.0D0) GO TO 1044
1043  CONTINUE
      STIF(L,L)=STIF(L,L)+(0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I))
      STIF(L-MX,L)=-0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I)
      GO TO 1046
1044  CONTINUE
      STIF(L,L)=STIF(L,L)-(0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I))
      STIF(L+MX,L)=0.5D0*TIME/(Y(L)-Y(L-MX))*VEL(2,J,I)
      GO TO 1046
1046  CONTINUE
!=====================================================
      IF (ID(2).EQ.0) GO TO 100
      STIF(L,L)=STIF(L,L)+&
     &(TIME*AGAMA/DEN/ETAS*&
     &((1.0D0/(X(L)-X(L-1))**2.0D0)+(1.0D0/((Y(L)-Y(L-MX))**2.0D0))))
      STIF(L,L)=STIF(L,L)+(AKAPA*TIME/DEN/ETAS)
!=====================================================
!=====================================================
      STIF(L+1,L)=STIF(L+1,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((X(L)-X(L-1))**2.0D0))
!=====================================================
      STIF(L+MX,L)=STIF(L+MX,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((Y(L)-Y(L-MX))**2.0D0))
!=====================================================
      STIF(L-1,L)=STIF(L-1,L)-&
     &(0.5D0*TIME*AGAMA/DEN/ETAS/((X(L)-X(L-1))**2.0D0))
!=====================================================
      STIF(L-MX,L)=STIF(L-MX,L)-(0.5D0*TIME*AGAMA/DEN/ETAS/&
     &((Y(L)-Y(L-MX))**2.0D0))
      GO TO 100
!=====================================================
!=====================================================
106   CONTINUE
      STIF(L,L)=1.0D0
      GO TO 100
!=====================================================
!=====================================================
5     FORMAT(F15.8)
100   CONTINUE
      RETURN
      END
!======================================================
!======================================================
!======================================================
!======================================================
      SUBROUTINE GSOLVE2
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(NP,NP),BB(NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=MP
      NNPP=NP
      MM=1
      MMPP=1
!======================================================
!======================================================
      DO 100 I=1,MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
!======================================================
      DO 200 I=1,MP
      DO 200 J=1,MP
!======================================================
      AA(J,I)=STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF(ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF(ID(10).Eq.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      DO 300 I=1,MP
      GYRA(LIST(1,I),LIST(2,I))=BB(I,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
      SUBROUTINE OUTPUT
      USE GLOBAL
!======================================================
1     FORMAT(16I7)
2     FORMAT(10X,'This is',I7,'-th step',I7/)
3     FORMAT(/'Number of Node',&
     &10X,'Coordinate'/)
4     FORMAT(I7,3F15.8,I5)
5     FORMAT(/'Number of Element',&
     &10X,'Connectivity'/)
6     FORMAT(10I7)
!======================================================
      DO 100 I=1,NPRINT
      IF (ITIME.EQ.IPRINT(I)) GO TO 101
      GO TO 199
!======================================================
101   CONTINUE
      WRITE(6,2) ITIME
!     CALL CAUCHYSTRESS
!     CALL MOMENTSTRESS
      CALL PLOT2
!     CALL TOTAL
!     CALL PLOT4
      CALL RESTART_WRITE
      GO TO 100
!======================================================
199   CONTINUE
100   CONTINUE 
!======================================================
      RETURN
      END
!======================================================
!======================================================
      SUBROUTINE PREC1
      USE GLOBAL
!======================================================
!======================================================
      L=0
      STIF=0.0D0
      BBB=0.0D0
!     WRITE(6,*)'I am Here in PREC1.'
!======================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!======================================================
      L=L+1  
      IF(I.EQ.5 .AND. J.EQ.5) GO TO 106   
!     IF (KBC(1,L).GT.0.1D0) GO TO 106 
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 102
      IF (J.EQ.1) GO TO 103
      IF (J.EQ.MX) GO TO 104
      GO TO 105
!======================================================
!======================================================
101   CONTINUE
      BBB(L)=0.0D0
      IF (ID(4).EQ.0) GO TO 1111
      IF (ID(4).EQ.1) GO TO 1112
!======================================================
1111  CONTINUE
      IF (J.EQ.1) GO TO 1011
      IF (J.EQ.MX) GO TO 1012
      GO TO 1013
1011  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
1012  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
1013  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L+1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
!======================================================
!======================================================
1112  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+MX,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
102   CONTINUE
      BBB(L)=0.0D0
      IF (ID(5).EQ.0) GO TO 1121
      IF (ID(5).EQ.1) GO TO 1122
!======================================================
1121  CONTINUE
      IF (J.EQ.1) GO TO 1021
      IF (J.EQ.MX) GO TO 1022
      GO TO 1023
1021  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
1022  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
1023  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
!======================================================
1122  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-MX,L)=-1.0D0
      GO TO 100   
!======================================================
!======================================================
103   CONTINUE
      BBB(L)=0.0D0
      IF (ID(6).EQ.0) GO TO 1131
      IF (ID(6).EQ.1) GO TO 1132
!======================================================
1131  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L+MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
1132  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+1,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
104   CONTINUE
      BBB(L)=0.0D0
      IF (ID(7).EQ.0) GO TO 1141
      IF (ID(7).EQ.1) GO TO 1142
!======================================================
1141  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
1142  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+1,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
105   CONTINUE
      BBB(L)=DEN/2.5D-1/TIME*&
     &(((VEL(1,J+1,I)-VEL(1,J-1,I))/(X(L+1)-X(L-1)))+&
     &((VEL(2,J,I+1)-VEL(2,J,I-1))/(Y(L+MX)-Y(L-MX))))
!======================================================
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=1.0D0/((X(L+1)-X(L))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=1.0D0/((Y(L+MX)-Y(L))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
!======================================================
!======================================================
106   CONTINUE
      BBB(L)=0.0D0
!======================================================
      STIF(L,L)=1.0D0
      GO TO 100
!======================================================
!======================================================
!======================================================
!======================================================
100   CONTINUE
!======================================================
1     FORMAT(2I5,10F15.8)
!     DO 999 K=1,L
!     WRITE(6,1) LIST(1,K),LIST(2,K),SBB(K),SGRE(K,K),SGRE(K+1,K),&
!    &SGRE(K-1,K),SGRE(K+MX,K),SGRE(K-MX,K)&
!    &,VEL(1,LIST(1,K),LIST(2,K))&
!    &,VEL(2,LIST(1,K),LIST(2,K))
999   CONTINUE
!======================================================
      CALL PSOLVE1
!======================================================
      K=0
      DO 200 J=1,MY
      DO 200 I=1,MX
      K=K+1
      IF (J.EQ.1.OR.J.EQ.MY) GO TO 200
      IF (I.EQ.1.OR.I.EQ.MX) GO TO 200
      VEL(1,I,J)=VEL(1,I,J)-(2.5D-1*TIME/DEN*&
     &((PRE(I+1,J)-PRE(I-1,J))/(X(K+1)-X(K-1))))
      VEL(2,I,J)=VEL(2,I,J)-(2.5D-1*TIME/DEN*&
     &((PRE(I,J+1)-PRE(I,J-1))/(Y(K+MX)-Y(K-MX))))
200   CONTINUE
!======================================================
      RETURN
      END
!======================================================
!======================================================
!======================================================
!======================================================
      SUBROUTINE PREC2
      USE GLOBAL
!======================================================
!======================================================
      L=0
      STIF=0.0D0
      BBB=0.0D0
!======================================================
      DO 100 I=1,MY
      DO 100 J=1,MX
!======================================================
      L=L+1      
      IF(I.EQ.5 .AND. J.EQ.5) GO TO 106  
!     IF (KBC(1,L).GT.0.1D0) GO TO 106  
      IF (I.EQ.1) GO TO 101
      IF (I.EQ.MY) GO TO 102
      IF (J.EQ.1) GO TO 103
      IF (J.EQ.MX) GO TO 104
      GO TO 105
!======================================================
!======================================================
101   CONTINUE
      BBB(L)=0.0D0 
      IF (ID(4).EQ.0) GO TO 1111
      IF (ID(4).EQ.1) GO TO 1112
!======================================================
1111  CONTINUE
      IF (J.EQ.1) GO TO 1011
      IF (J.EQ.MX) GO TO 1012
      GO TO 1013
1011  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
1012  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
1013  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L+MX)-Y(L))**2)))
      STIF(L+1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=2.0D0/((Y(L+MX)-Y(L))**2)
      GO TO 100
!======================================================
1112  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+MX,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
102   CONTINUE
      BBB(L)=0.0D0
      IF (ID(5).EQ.0) GO TO 1121
      IF (ID(5).EQ.1) GO TO 1122
!======================================================
!======================================================
1121  CONTINUE
      IF (J.EQ.1) GO TO 1021
      IF (J.EQ.MX) GO TO 1022
      GO TO 1023
1021  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
1022  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
1023  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L-MX,L)=2.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
1122  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-MX,L)=-1.0D0
      GO TO 100   
!======================================================
!======================================================
103   CONTINUE
      BBB(L)=0.0D0
      IF (ID(6).EQ.0) GO TO 1131
      IF (ID(6).EQ.1) GO TO 1132
!======================================================
!======================================================
1131  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L+1)-X(L))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=2.0D0/((X(L+1)-X(L))**2)
      STIF(L+MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
1132  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L+1,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
104   CONTINUE
      BBB(L)=0.0D0
      IF (ID(7).EQ.0) GO TO 1141
      IF (ID(7).EQ.1) GO TO 1142
!======================================================
!======================================================
1141  CONTINUE
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L-1,L)=2.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
1142  CONTINUE
      STIF(L,L)=1.0D0
!     STIF(L-1,L)=-1.0D0
      GO TO 100
!======================================================
!======================================================
105   CONTINUE
      BBB(L)=DEN/2.5D-1/TIME*&
     &(((VEL(1,J+1,I)-VEL(1,J-1,I))/(X(L+1)-X(L-1)))+&
     &((VEL(2,J,I+1)-VEL(2,J,I-1))/(Y(L+MX)-Y(L-MX))))
!======================================================
      STIF(L,L)=-2.0D0*((1.0D0/((X(L)-X(L-1))**2))+&
     &(1.0D0/((Y(L)-Y(L-MX))**2)))
      STIF(L+1,L)=1.0D0/((X(L+1)-X(L))**2)
      STIF(L-1,L)=1.0D0/((X(L)-X(L-1))**2)
      STIF(L+MX,L)=1.0D0/((Y(L+MX)-Y(L))**2)
      STIF(L-MX,L)=1.0D0/((Y(L)-Y(L-MX))**2)
      GO TO 100
!======================================================
!======================================================
!======================================================
!======================================================
106   CONTINUE
      BBB(L)=0.0D0
!======================================================
      STIF(L,L)=1.0D0
      GO TO 100
!======================================================
!======================================================
!======================================================
100   CONTINUE
!======================================================
!======================================================
1     FORMAT(2I5,50F15.8)
!     DO 999 K=1,L
!     WRITE(6,1) LIST(1,K),LIST(2,K),SBB(K),SGRE(K,K),SGRE(K+1,K),&
!    &SGRE(K-1,K),SGRE(K+MX,K),SGRE(K-MX,K)&
!    &,VEL(1,LIST(1,K),LIST(2,K)),VEL(1,LIST(1,K)+1,LIST(2,K)),&
!    &VEL(1,LIST(1,K)-1,LIST(2,K)),VEL(1,LIST(1,K),LIST(2,K)+1),&
!    &VEL(1,LIST(1,K),LIST(2,K)-1)
999   CONTINUE
!======================================================
      CALL PSOLVE2
!======================================================
      K=0
      DO 200 J=1,MY
      DO 200 I=1,MX
      K=K+1
      IF (J.EQ.1.OR.J.EQ.MY) GO TO 200
      IF (I.EQ.1.OR.I.EQ.MX) GO TO 200
      VEL(1,I,J)=VEL(1,I,J)-(2.5D-1*TIME/DEN*&
     &((PRE(I+1,J)-PRE(I-1,J))/(X(K+1)-X(K-1))))
      VEL(2,I,J)=VEL(2,I,J)-(2.5D-1*TIME/DEN*&
     &((PRE(I,J+1)-PRE(I,J-1))/(Y(K+MX)-Y(K-MX))))
200   CONTINUE
!======================================================
      RETURN
      END
!======================================================
!======================================================
!======================================================
!======================================================
      SUBROUTINE PSOLVE2
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(NP,NP),BB(NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=MP
      NNPP=NP
      MM=1
      MMPP=1
!======================================================
!======================================================
      DO 100 I=1,MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
!======================================================
      DO 200 I=1,MP
      DO 200 J=1,MP
!======================================================
      AA(J,I)=STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF(ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF(ID(10).EQ.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      DO 300 I=1,MP
      PRE(LIST(1,I),LIST(2,I))=BB(I,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
!======================================================
!======================================================
      SUBROUTINE PSOLVE1
      USE GLOBAL
!======================================================
!======================================================
      REAL*8 AA(NP,NP),BB(NP,1)
      INTEGER NN,NNPP,MM,MMPP
      NN=MP
      NNPP=NP
      MM=1
      MMPP=1
!======================================================
!======================================================
      DO 100 I=1,MP
      BB(I,1)=BBB(I)
100   CONTINUE
!======================================================
!======================================================
      DO 200 I=1,MP
      DO 200 J=1,MP
!======================================================
      AA(J,I)=STIF(I,J)
!======================================================
200   CONTINUE
!======================================================
      IF(ID(10).EQ.0) CALL gaussj(AA,NN,NNPP,BB,MM,MMPP)
      IF(ID(10).EQ.1) CALL SOR(AA,NN,BB,NNPP)
!======================================================
      DO 300 I=1,MP
      PRE(LIST(1,I),LIST(2,I))=BB(I,1)
300   CONTINUE 
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
      SUBROUTINE CAUCHYSTRESS
      USE GLOBAL
!======================================================
!======================================================
      DO 100 I=1,MP
      IF (LIST(1,I).EQ.1) GO TO 101
      IF (LIST(1,I).EQ.MX) GO TO 102
      IF (LIST(2,I).EQ.1) GO TO 103
      IF (LIST(2,I).EQ.MY) GO TO 104
      GO TO 105
!======================================================
101   CONTINUE
      CAUCHY(I,1,1)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,2,2)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,1,2)=(AKAPA*(((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(VEL(2,LIST(1,I),LIST(2,I))))/&
     &((X(I+1)-X(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(VEL(2,LIST(1,I),LIST(2,I))))/&
     &((X(I+1)-X(I))))))
      CAUCHY(I,2,1)=(AKAPA*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(VEL(2,LIST(1,I),LIST(2,I))))/&
     &((X(I+1)-X(I))))))
      GO TO 100
!======================================================
102   CONTINUE
      CAUCHY(I,1,1)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,2,2)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,1,2)=(AKAPA*(((&
     &(VEL(2,LIST(1,I),LIST(2,I)))-VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I)-X(I-1))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((&
     &(VEL(2,LIST(1,I),LIST(2,I)))-&
     &VEL(2,LIST(1,I)-1,LIST(2,I)-1))/&
     &(2.0D0*(X(I)-X(I-1))))))
      CAUCHY(I,2,1)=(AKAPA*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((&
     &(VEL(2,LIST(1,I),LIST(2,I)))-&
     &VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &((X(I)-X(I-1))))))
      GO TO 100
!======================================================
103   CONTINUE
      CAUCHY(I,1,1)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,2,2)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,1,2)=(AKAPA*(((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(VEL(1,LIST(1,I),LIST(2,I))))/&
     &((Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)-1))/&
     &(2.0D0*(X(I+1)-X(I))))))
      CAUCHY(I,2,1)=(AKAPA*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(VEL(1,LIST(1,I),LIST(2,I))))/&
     &((Y(I+MX)-Y(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(VEL(1,LIST(1,I),LIST(2,I))))/&
     &((Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))))
      GO TO 100
!======================================================
104   CONTINUE
      CAUCHY(I,1,1)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,2,2)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,1,2)=(AKAPA*(((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((&
     &(VEL(1,LIST(1,I),LIST(2,I)))-VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &((Y(I)-Y(I-MX))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)-1))/&
     &(2.0D0*(X(I+1)-X(I))))))
      CAUCHY(I,2,1)=(AKAPA*(((&
     &(VEL(1,LIST(1,I),LIST(2,I)))-VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &((Y(I+MX)-Y(I-MX))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((&
     &(VEL(1,LIST(1,I),LIST(2,I)))-VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &((Y(I)-Y(I-MX))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))))
      GO TO 100
!======================================================
105   CONTINUE
      CAUCHY(I,1,1)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,2,2)=-PRE(LIST(1,I),LIST(2,I))
      CAUCHY(I,1,2)=(AKAPA*(((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)-1))/&
     &(2.0D0*(X(I+1)-X(I))))))
      CAUCHY(I,2,1)=(AKAPA*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))-&
     &GYRA(LIST(1,I),LIST(2,I))))+&
     &(AMU*(((VEL(1,LIST(1,I),LIST(2,I)+1)-&
     &(2.0D0*VEL(1,LIST(1,I),LIST(2,I)))+VEL(1,LIST(1,I),LIST(2,I)-1))/&
     &(2.0D0*(Y(I+MX)-Y(I))))&
     &+((VEL(2,LIST(1,I)+1,LIST(2,I))-&
     &(2.0D0*VEL(2,LIST(1,I),LIST(2,I)))+&
     &VEL(2,LIST(1,I)-1,LIST(2,I)))/&
     &(2.0D0*(X(I+1)-X(I))))))
      GO TO 100
!======================================================
100   CONTINUE
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
      SUBROUTINE OUTPUT1
!     USE COORS
      USE GLOBAL
!======================================================
!===================================================
1     FORMAT('Title="Micropolar Fluid Dynamics-Total Velocity"')
2     FORMAT('Variables="X","Y","VEL_X","VEL_Y"')
3     FORMAT('Zone T="Load Step 0",N=',I5,' E=',I5,' datapacking=point,&
     &zonetype=fequadrilateral')
13    FORMAT('Zone T="Load Step 0",N=',I5,' E=',I5,'')
4     FORMAT(25F15.8)
5     FORMAT(8I10)
6     FORMAT('It is an inviscid flow.')
7     FORMAT('Viscosity is considered.')
!====================================================
!======================================================
      XMIN=X(1)
      YMIN=Y(1)
      XMAX=X(MX)
      YMAX=Y(MP-MX+1)
!     YMAX=Y(MP)
!======================================================
      MMX=NF*NS*(MX-1)
      MMY=NF*NS*(MY-1)
      MME=MMX*MMY
      MMP=(MMX+1)*(MMY+1)
!======================================================
      DDX=(XMAX-XMIN)/(MMX)
      DDY=(YMAX-YMIN)/(MMY)
!======================================================
      WRITE(8,1)
      WRITE(8,2)
      WRITE(8,3) MMP,MME
!======================================================
!     MMX1 = MMX+1
!     MMY1 = MMY+1
!     ALLOCATE(IIJJKK(4,MMX*MMY))
!     ALLOCATE(XCOOR(MMX1,MMY1))
!     ALLOCATE(YCOOR(MMX1,MMY1))
!======================================================
      K=0
      DO 100 IY=1,MMY+1
      YY=(IY-1)*DDY
      DO 100 IX=1,MMX+1
      XX=(IX-1)*DDX
      K=K+1
      XCOOR(IX,IY)=XX
      YCOOR(IX,IY)=YY
      WRITE(8,4) XCOOR(IX,IY),YCOOR(IX,IY),0.0D0,0.0D0
100   CONTINUE
!======================================================
      L=0
      DO 101 IY=1,MMY
      JJ=(IY-1)*(MMX+1)
      DO 101 IX=1,MMX
      L=L+1
      IIJJKK(1,L)=JJ+IX
      IIJJKK(2,L)=IIJJKK(1,L)+1
      IIJJKK(4,L)=IIJJKK(1,L)+MMX+1
      IIJJKK(3,L)=IIJJKK(4,L)+1
      WRITE(8,5) (IIJJKK(II,L),II=1,4)
101   CONTINUE
!======================================================
!======================================================
!======================================================
!======================================================
      RETURN
      END
!======================================================
!======================================================
!======================================================
!====================================================
      SUBROUTINE TOTAL
!     USE COORS
      USE GLOBAL
!====================================================
!===================================================
!===================================================
!===================================================
1     FORMAT('Title="Micropolar Fluid Dynamics-Total Velocity"')
!2     FORMAT('Variables="X","Y","VEL_X","VEL_Y","FVEL_X","FVEL_Y"')
2     FORMAT('Variables="X","Y","VEL_X","VEL_Y"')
3     FORMAT('Zone T="Load Step',I9'",N=',I5,' E=',I5,' &
     &datapacking=point,&
     &zonetype=fequadrilateral,connectivitysharezone=1')
4     FORMAT(25F15.8)
!====================================================
      REAL*8 TVEL(2,NF*NS*NX,NF*NS*NY),FVEL(2,NF*NS*NX,NF*NS*NY)
!====================================================
      WRITE(8,3) ITIME,(MMX+1)*(MMY+1),MMX*MMY
      EERTIA=DSQRT(ERTIA/2.0D0)
      TVEL=0.0D0
      FVEL=0.0D0
!===================================================
      DO 200 IY=1,MMY+1,NS
      DO 200 IX=1,MMX+1,NS
!===================================================
      DO 201 I=1,ME
!===================================================
      DDDX=DABS(X(IJK(2,I))-X(IJK(1,I)))
      DDDY=DABS(Y(IJK(4,I))-Y(IJK(1,I)))
      DISTX=DABS(XCOOR(IX,IY)-X(IJK(1,I)))
      DISTY=DABS(YCOOR(IX,IY)-Y(IJK(1,I)))
      IF (DISTX.GT.DDDX.OR.DISTY.GT.DDDY) GO TO 201
      DISTX=DABS(XCOOR(IX,IY)-X(IJK(2,I)))
      DISTY=DABS(YCOOR(IX,IY)-Y(IJK(2,I)))
      IF (DISTX.GT.DDDX.OR.DISTY.GT.DDDY) GO TO 201
      DISTX=DABS(XCOOR(IX,IY)-X(IJK(3,I)))
      DISTY=DABS(YCOOR(IX,IY)-Y(IJK(3,I)))
      IF (DISTX.GT.DDDX.OR.DISTY.GT.DDDY) GO TO 201
      DISTX=DABS(XCOOR(IX,IY)-X(IJK(4,I)))
      DISTY=DABS(YCOOR(IX,IY)-Y(IJK(4,I)))
      IF (DISTX.GT.DDDX.OR.DISTY.GT.DDDY) GO TO 201
!===================================================
202   CONTINUE
      XXX=-1.0D0+(2.0D0*(XCOOR(IX,IY)-X(IJK(4,I)))/DDDX)
      YYY=-1.0D0+(2.0D0*(YCOOR(IX,IY)-Y(IJK(1,I)))/DDDY)
      S1=0.25D0*(1.0D0-XXX)*(1.0D0-YYY)
      S2=0.25D0*(1.0D0+XXX)*(1.0D0-YYY)
      S3=0.25D0*(1.0D0+XXX)*(1.0D0+YYY)
      S4=0.25D0*(1.0D0-XXX)*(1.0D0+YYY)
      TVEL(1,IX,IY)=(S1*VEL(1,LIST(1,IJK(1,I)),LIST(2,IJK(1,I))))&
     &+(S2*VEL(1,LIST(1,IJK(2,I)),LIST(2,IJK(2,I))))&
     &+(S3*VEL(1,LIST(1,IJK(3,I)),LIST(2,IJK(3,I))))&
     &+(S4*VEL(1,LIST(1,IJK(4,I)),LIST(2,IJK(4,I))))
      TVEL(2,IX,IY)=(S1*VEL(2,LIST(1,IJK(1,I)),LIST(2,IJK(1,I))))&
     &+(S2*VEL(2,LIST(1,IJK(2,I)),LIST(2,IJK(2,I))))&
     &+(S3*VEL(2,LIST(1,IJK(3,I)),LIST(2,IJK(3,I))))&
     &+(S4*VEL(2,LIST(1,IJK(4,I)),LIST(2,IJK(4,I))))
      FVEL(1,IX,IY)=(S1*GYRA(LIST(1,IJK(1,I)),LIST(2,IJK(1,I))))&
     &+(S2*GYRA(LIST(1,IJK(2,I)),LIST(2,IJK(2,I))))&
     &+(S3*GYRA(LIST(1,IJK(3,I)),LIST(2,IJK(3,I))))&
     &+(S4*GYRA(LIST(1,IJK(4,I)),LIST(2,IJK(4,I))))
      FVEL(2,IX,IY)=(S1*GYRA(LIST(1,IJK(1,I)),LIST(2,IJK(1,I))))&
     &+(S2*GYRA(LIST(1,IJK(2,I)),LIST(2,IJK(2,I))))&
     &+(S3*GYRA(LIST(1,IJK(3,I)),LIST(2,IJK(3,I))))&
     &+(S4*GYRA(LIST(1,IJK(4,I)),LIST(2,IJK(4,I))))
      GO TO 201
!===================================================
!===================================================
201   CONTINUE
!===================================================
200   CONTINUE
!===================================================
!===================================================
      DO 100 IY=1,MMY+1
      DO 100 IX=1,MMX+1
!     GO TO 900
!===================================================
      EN=0.0D0
      DO 101 IIY=1,MMY+1,NS
      DO 101 IIX=1,MMX+1,NS
      DIST=DSQRT(((XCOOR(IIX,IIY)-XCOOR(IX,IY))**2.0D0)+&
     &((YCOOR(IIX,IIY)-YCOOR(IX,IY))**2.0D0))
      IF (DIST.GT.EERTIA) GO TO 104
!===================================================
103   CONTINUE
!     EN=EN+1.0D0
      ETA1=XCOOR(IX,IY)-XCOOR(IIX,IIY)
      ETA2=YCOOR(IX,IY)-YCOOR(IIX,IIY)
!     TVEL(1,IX,IY)=TVEL(1,IX,IY)+VEL(1,LIST(1,I),LIST(2,I))&
!    &-(GYRA(LIST(1,I),LIST(2,I))*ETA2)
!     TVEL(2,IX,IY)=TVEL(2,IX,IY)+VEL(2,LIST(1,I),LIST(2,I))&
!    &+(GYRA(LIST(1,I),LIST(2,I))*ETA1)
!     FVEL(1,IX,IY)=FVEL(1,IX,IY)+&
!    &-(GYRA(LIST(1,I),LIST(2,I))*ETA2)
!     FVEL(2,IX,IY)=FVEL(2,IX,IY)+&
!    &+(GYRA(LIST(1,I),LIST(2,I))*ETA1)
      TVEL(1,IX,IY)=TVEL(1,IIX,IIY)&
     &-(FVEL(1,IIX,IIY)*ETA2)
      TVEL(2,IX,IY)=TVEL(2,IIX,IIY)&
     &+(FVEL(1,IIX,IIY)*ETA1)
99    FORMAT(12F15.8)
      GO TO 101
!===================================================
!===================================================
104   CONTINUE
!     WRITE(6,99) X(I),Y(I),XCOOR(IX,IY),YCOOR(IX,IY),DIST,EERTIA
101   CONTINUE
!===================================================
!     IF (EN.EQ.0.0D0) WRITE(6,*) IX,IY
!     TVEL(1,IX,IY)=TVEL(1,IX,IY)/EN
!     TVEL(2,IX,IY)=TVEL(2,IX,IY)/EN
!     FVEL(1,IX,IY)=FVEL(1,IX,IY)/EN
!     FVEL(2,IX,IY)=FVEL(2,IX,IY)/EN
!===================================================
900   CONTINUE
      WRITE(8,4) XCOOR(IX,IY),YCOOR(IX,IY),&
     &(TVEL(LL,IX,IY),LL=1,2)
100   CONTINUE
!===================================================
!     DEALLOCATE (XCOOR)
!     DEALLOCATE (YCOOR)
!     DEALLOCATE (IIJJKK)
!===================================================
      RETURN
      END
!====================================================
!====================================================
!====================================================
      SUBROUTINE SOR(AA,N,BB,NNP)
      USE GLOBAL
      INTEGER N,I,J,K,NNP,ITS,ITERS
      REAL*8 AA(NNP,NNP), BB(NNP,1), XX(N), XNEW(N),DIA,ER,TEMP
      REAL*8 TOL,AME
!1     FORMAT('INPUT NUMBER OF EQUATIONS N=',5I)
!     XX: INITIAL GUESS
!     XNEW: INITIAL GUESS
      XX=0.0D0
      XNEW=0.0D0
!     TOL: TOLENCE
!     ITS: MAX ITERATION
      TOL=1.0E-6
      ITS=5000
!     AME: RELAXATION FACTOR
!     AME=1 -> Gauss-Siedel Algorithm
      AME=1.0D0
!=====================================================
!=====================================================
!==========================================================
!==========================================================
      ITERS=0
4     CONTINUE
      ITERS=ITERS+1
!==========================================================
!==========================================================
      DO 100 I=1,N
      DIA=AA(I,I)
      TEMP=0.0D0
!==========================================================
      DO 200 J=1,I-1
      TEMP=TEMP+(AA(I,J)*XNEW(J))
200   CONTINUE
!==========================================================
      DO 300 K=I+1,N
      TEMP=TEMP+(AA(I,K)*XX(K))
300   CONTINUE
!==========================================================
      XNEW(I)=((1.0D0-AME)*XX(I))+(AME/DIA*(BB(I,1)-TEMP))
100   CONTINUE
!==========================================================
!==========================================================
      ER=0.0D0
      DO 101 I=1,N
      ER=ER+(DABS(XX(I)-XNEW(I)))
101   CONTINUE
!==========================================================
!==========================================================
      DO 111 I=1,N
      XX(I)=XNEW(I)
111   CONTINUE
      IF (ER.LT.TOL.OR.ITERS.GE.ITS) GO TO 5
!==========================================================
!==========================================================
      GO TO 4
5     CONTINUE
!==========================================================
!==========================================================
!==========================================================
9876  FORMAT(1I5'-th step',5X,'The error is',1F15.12&
     &,5X,'The final iteration is',1I8)
      WRITE(9,9876) ITIME, ER, ITERS
!==========================================================
      DO 222 I=1,N
      BB(I,1)=XX(I)
222   CONTINUE
!==========================================================
      RETURN
      END
!==========================================================
!==========================================================
!==========================================================
!==========================================================
      SUBROUTINE  RESTART_WRITE
      USE GLOBAL
      OPEN(10, FILE='restart', STATUS='UNKNOWN',FORM='FORMATTED')
!==========================================================
      DO 100 I=1,MP
      WRITE(10,1) (VEL(K,LIST(1,I),LIST(2,I)),K=1,2),&
     &GYRA(LIST(1,I),LIST(2,I)),PRE(LIST(1,I),LIST(2,I))
100   CONTINUE
!==========================================================
1     FORMAT(5F15.8)
      CLOSE(10)
      RETURN
      END
!==========================================================
!==========================================================
      SUBROUTINE  RESTART_READ
      USE GLOBAL
      OPEN(10, FILE='restart', STATUS='UNKNOWN',FORM='FORMATTED')
!==========================================================
      DO 100 I=1,MP
      READ(10,1) (VEL(K,LIST(1,I),LIST(2,I)),K=1,2),&
     &GYRA(LIST(1,I),LIST(2,I)),PRE(LIST(1,I),LIST(2,I))
100   CONTINUE
!==========================================================
1     FORMAT(5F15.8)
      CLOSE(10)
      RETURN
      END
!==========================================================
!==========================================================

