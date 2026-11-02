!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFWT(YDEPHY,KIDIA,KFDIA,KLON,KLEV,&
 & PUM1,PVM1,PTM1,PQM1,PGEOM1, PAPM1,&
 & PKMFL,PKHFL,PKQFL,&
 & PDHPBL, WDIR_PBL, UPBL, VPBL,WTURN_PBL, MF_PBL)  
!     ------------------------------------------------------------------

!**   *VDFWT* - VDFDPBL with added wind turning and boundary layer 
!                cross isobaric mass flux for diagnostic purposes only

!     Original
!**   *VDFDPBL* - VDFDPBL (Diagnostic PBL height) determines  
!                 PBL height for diagnostic purposes only

!     A.C.M. BELJAARS       E.C.M.W.F.    17/02/1998.

!     PURPOSE
!     -------

!     Determine PBL height

!     INTERFACE
!     ---------

!     *VDFDPBL* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PAPM1*        PRESSURE AT T-1
!     *PKMFL*        KINEMATIC MOMENTUM FLUX                
!     *PKHFL*        KINEMATIC HEAT FLUX                    
!     *PKQFL*        KINEMATIC MOISTURE FLUX           

!     OUTPUT PARAMETERS (REAL):

!     *PDHPBL*        Boundary layer height                  m
!     *WDIR_PBL*      Wind direction at Boundary layer height (PDHPBL) 
!     *UPBL*          U component of wind at PDHPBL
!     *VPBL*          V component of wind at PDHPBL
!     *WTURN_PBL*     Wind turning between lowest layer and PDHPBL
!     *MF_PBL*        Ageostrophic mass flux in boundary layer

!     METHOD
!     ------

!     Troen and Mahrt method using bulk Richardson criterion
!     (see documentation)
!
!     MODIFICATIONS.
!     --------------
!     N. Semane+P.Bechtold     04-10-2012 Add RPLRG/RPLDARE factor for small planet
!     I.Sandu 15/03/2013 new algorithm based on Seidel et al. 2012
!     A. Lewinscahl 05-2025   Wind turning diagmostocs and ageostrophic mass flux
!     ------------------------------------------------------------------



!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RCPD, RETV, RPI,RD
USE PARPHY    , ONLY : REPDU2
USE YOEPHY    , ONLY : TEPHY
USE YOMDYNCORE, ONLY : RPLRG
USE YOMLUN        , ONLY : NULOUT
IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHPBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: WDIR_PBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: UPBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: VPBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: WTURN_PBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: MF_PBL(KLON) 


!*            LOCAL STORAGE
!             ----- -------

LOGICAL :: LLDONE(KLON)
REAL(KIND=JPRB) :: ZRI(KLON),ZDU2(KLON),&
 & ZSVBOT(KLON),ZSVBOTP(KLON)


INTEGER(KIND=JPIM) :: ILEVM1, ITOT, JIT, JK, JL

REAL(KIND=JPRB) :: ZBUST, ZCONS13, ZCONS14, ZCONS15, ZDRORO,&
 & ZEPS, ZPAR, ZPAR1, ZPARZI, ZRICRI, ZRILEV, ZSV, ZREPUST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM)  :: PBLIDX(KLON)
REAL(KIND=JPRB)     :: WDIR_0(KLON) 
REAL(KIND=JPRB)     :: WDIR(KLON,KLEV) 
REAL(KIND=JPRB)     :: AGWIND(KLON,KLEV) 
REAL(KIND=JPRB)     :: ROTANG 
REAL(KIND=JPRB)     :: ZRHO 

#include "surf_inq.h"

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFDPBL',0,ZHOOK_HANDLE)

ZPAR   = 8.5_JPRB
ZPAR1  = 0.6_JPRB
ZCONS14= RCPD*ZPAR
ZCONS15= ZPAR1*RG
ZEPS   = 1.E-10_JPRB
ZBUST  = 100._JPRB
ZRICRI = 0.25_JPRB
ZPARZI = 1000._JPRB/RPLRG

CALL SURF_INQ(YDEPHY%YSURF,PREPUST=ZREPUST)

ZCONS13=1.0_JPRB/3._JPRB

ILEVM1 = KLEV-1
MF_PBL = 0._JPRB

!        1.    PREPARE SURFACE PARAMETERS
!              ------- ------- ----------

DO JL=KIDIA,KFDIA
  ZSVBOT(JL)=RCPD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV))+PGEOM1(JL,KLEV)
  PDHPBL(JL)=ZPARZI
ENDDO

!        2.    DO 1 OR 2 ITERATIONS ON PBL HEIGHT
!              -- - -- - ---------- -- --- ------

!*************************************
JIT=1
!*************************************

!        2.1   UPDATE VELOCITY SCALE AND EXCESS TEMPERATURE
!              ------ -------- ----- --- ------ -----------

  DO JL=KIDIA,KFDIA
    LLDONE(JL)=.FALSE.
    ZRI(JL)=0.0_JPRB
    ZSVBOTP(JL)=ZSVBOT(JL)+ZEPS
  ENDDO

!        2.2    VERTICAL SCAN TO DETERMINE MIXED LAYER DEPTH
!               -------- ---- -- --------- ----- ----- -----

!***
  DO JK=ILEVM1,1,-1
!***

    ITOT=KFDIA-KIDIA+1
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL)) THEN
        ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
        ZDU2(JL)=MAX(REPDU2, PUM1(JL,JK)**2+PVM1(JL,JK)**2) 
        ZDRORO=(ZSV-ZSVBOTP(JL)) &
         & /(ZSV-PGEOM1(JL,JK)) 
!
        ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)
        IF (ZRILEV  >  ZRICRI) THEN
          PDHPBL(JL)=( (ZRILEV-ZRICRI)*PGEOM1(JL,JK+1)&
           & +(ZRICRI-ZRI(JL))*PGEOM1(JL,JK) )/&
           & ((ZRILEV-ZRI(JL))*RG)  

          ! U component of wind at PBL
          UPBL(JL)=( (ZRILEV-ZRICRI)*PUM1(JL,JK+1)&
           & +(ZRICRI-ZRI(JL))*PUM1(JL,JK) )/&
           & (ZRILEV-ZRI(JL))  

          ! V component of wind at PBL
          VPBL(JL)=( (ZRILEV-ZRICRI)*PVM1(JL,JK+1)&
           & +(ZRICRI-ZRI(JL))*PVM1(JL,JK) )/&
           & (ZRILEV-ZRI(JL))  

          ! Wind direction at PBL
          WDIR_PBL(JL) = 270_JPRB - ((180_JPRB / RPI) * atan2(VPBL(JL),UPBL(JL)))
          ! Adjust so that wind direction falls between 0 and 360 and 0=north
          IF (WDIR_PBL(JL) > 360_JPRB) THEN
             WDIR_PBL(JL) = WDIR_PBL(JL) - 360_JPRB
          ENDIF

          LLDONE(JL)=.TRUE.
          ITOT=ITOT-1
          PBLIDX(JL)=JK
        ELSE
          ZRI(JL)=ZRILEV
        ENDIF
      ELSE
        ITOT=ITOT-1
      ENDIF
    ENDDO
    IF (ITOT  <=  0) EXIT
!***
  ENDDO
!***
!        2.3   Wind Turning Boundary Layer
!          -------------------------------

  DO JL=KIDIA,KFDIA !do t_i = 1, ncol

       WDIR_0(JL) = 270_JPRB - ((180_JPRB / RPI) * atan2(PVM1(JL,KLEV),PUM1(JL,KLEV)))
       ! Adjust so that wind direction falls between 0 and 360 according to met definition
       IF (WDIR_0(JL) > 360_JPRB) THEN
          WDIR_0(JL) = WDIR_0(JL) - 360_JPRB
       ENDIF

       WTURN_PBL(JL) = WDIR_PBL(JL) - WDIR_0(JL)    ! in reference to pbl height
       IF(WTURN_PBL(JL) < -180_JPRB) THEN
         WTURN_PBL(JL) = WTURN_PBL(JL) + 360_JPRB
       ELSE IF (WTURN_PBL(JL) > 180_JPRB) THEN
         WTURN_PBL(JL) = WTURN_PBL(JL) - 360_JPRB
       ENDIF
      
  ENDDO


!        2.4  Ageostrophic mass flux in boundary layer 
!        Massflux is calculated from surface to layer below boundarly layer height
!
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PDHPBL       
!       ----------------------------- PBLIDX       -|
!                                                   |
!       ----------------------------- PBLIDX-1      |
!                                                   |
!             ..............                        |
!                                                   | - integrated mass flux
!       ----------------------------- KLEV-1        |
!                                                   |
!       ----------------------------- KLEV          |
!                                                   |
!       _____________________________ surface      -|
!
  DO JL = KIDIA, KFDIA
    ! Calculate rotation angle
    ROTANG = (270_JPRB - WDIR_PBL(JL)) * (RPI / 180_JPRB)

    ! Calculate ageostrophic wind component perpendicular to geostropic wind.
    ! Coordinate reference system rotation for wind vector component normal to geostrophic wind. 
    ! To do: calculate geostropic wind at each height using pressure gradient force
    AGWIND(JL, KLEV) = -(PUM1(JL,KLEV) ) * sin(ROTANG) + (PVM1(JL,KLEV) ) * cos(ROTANG)
    ZRHO = PAPM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    MF_PBL(JL) = MF_PBL(JL) + ZRHO * AGWIND(JL,KLEV) * (PGEOM1(JL,KLEV)/RG)

    DO JK = ILEVM1,PBLIDX(JL)+1,-1  ! what happens at the ground???
         AGWIND(JL, JK) = -(PUM1(JL,JK) ) * sin(ROTANG) + (PVM1(JL,JK) ) * cos(ROTANG)
         ZRHO = PAPM1(JL,JK)/( RD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK)) )
         MF_PBL(JL) = MF_PBL(JL) + ZRHO * AGWIND(JL,JK) *( (PGEOM1(JL,JK)- PGEOM1(JL,JK+1) )/RG)
    ENDDO
  ENDDO



!*************************************
!*************************************

IF (LHOOK) CALL DR_HOOK('VDFWF',1,ZHOOK_HANDLE)
END SUBROUTINE VDFWT
