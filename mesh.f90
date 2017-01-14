      PROGRAM MESH
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NP=10000,NE=10000,NBC=10000,NT=10000)
      COMMON/ABC/ID(16)
      COMMON/BBC/JNP(NP),X(NP),Y(NP),Z(NP),IJK(8,NE)
      COMMON/CBC/LBC(2,NBC),BC(3,NBC)
      COMMON/DBC/ITIME(NT)
      COMMON/EBC/MP,MBC
!=========================================
      OPEN(6, FILE='infile', STATUS='UNKNOWN',FORM='FORMATTED')
!=========================================
1     FORMAT(15I9)
2     FORMAT(I9,3F15.8,I7)
3     FORMAT(2I9,F15.8)
4     FORMAT(3I9,5F15.8)
5     FORMAT(16F15.8)
!=========================================
!=========================================
      ID(1)=1
      ID(2)=1
      ID(3)=1
      ID(4)=0
      ID(5)=0
      ID(6)=0
      ID(7)=0
      ID(8)=1
      ID(9)=1
      ID(10)=1
      WRITE(6,1) (ID(K),K=1,10)
!=========================================
!=========================================
      NX=10
      NY=80
      XL=1.0D0
      YL=2.0D0
      DX=XL/NX
      DY=YL/NY
      MP=(NX+1)*(NY+1)
      ME=NX*NY
      WRITE(6,1) (NX+1),(NY+1),MP,ME
!=========================================
!=========================================
      K=0
      DO 100 IY=1,NY+1
      YY=-1.0+(IY-1)*DY
      DO 100 IX=1,NX+1
      XX=(IX-1)*DX
      K=K+1
      X(K)=XX
      Y(K)=YY
      WRITE(6,2) K,X(K),Y(K)
100   CONTINUE
      L=0
      DO 101 IY=1,NY
      JJ=(IY-1)*(NX+1)
      DO 101 IX=1,NX
      L=L+1
      IJK(1,L)=JJ+IX
      IJK(2,L)=IJK(1,L)+1
      IJK(4,L)=IJK(1,L)+NX+1
      IJK(3,L)=IJK(4,L)+1
      WRITE(6,1) L, (IJK(K,L),K=1,4)
101   CONTINUE
!=========================================
!=========================================
      WRITE(6,2) 10000,1.0D-1
      WRITE(6,1) 100
      DO 210 J=1,10
      DO 200 I=1,10
      ITIME((10*(J-1))+I)=1000*(J-1)+100*I
200   CONTINUE
      WRITE(6,1) (ITIME(10*(J-1)+K),K=1,10)
210   CONTINUE
!=========================================
!=========================================
      PA=1.0D-2
      PB=0.0D0
      L=0
      DO 300 IY=1,NY+1
      YY=(IY-1)*DY
      DO 300 IX=1,NX+1
      XX=(IX-1)*DX
      L=L+1
      WRITE (6,4) L, IX, IY, 0.0D-1, 0.0D-1, 0.0D0,0.0D0
!     WRITE (6,4) L, IX, IY, 0.0D0, 0.0D0, 0.0D0,PA+(PB-PA)/NX*(IX-1)
301   CONTINUE
300   CONTINUE
!=========================================
!=========================================
      DEN=1.0D0
!     ERTIA=2.0D0*&
!    &(DSQRT(((DX/2.0D0)**2.0D0)+((DY/2.0D0)**2.0D0))**2.0D0)
      ETAS=1.0D-3
!     ETAS=ERTIA
      GRA=0D-1
      AGAMA=1.0D-3
      ALAMDA=1.0D-3
      AMU=6.25D-5
      AKAPA=6.25D-5
      AGAMA=((0.5D0*AKAPA)+AMU)*ETAS
!     ETAS=2.0D0*AGAMA*(AMU+AKAPA)/AKAPA/(AKAPA+(2.0D0*AMU))
      WRITE(6,5) DEN,ERTIA,ETAS,GRA,AGAMA,ALAMDA,AMU,AKAPA
!=========================================
!=========================================
      CLOSE (6)
!=========================================
!=========================================
      END
