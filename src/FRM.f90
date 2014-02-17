SUBROUTINE FRM(SEED,NHOPS,NFRAMES,NDIMERS,RATES,DX,DY,DZ,CON, &
               DISTX,DISTY,DISTZ,TIME,NHOP)
  IMPLICIT NONE

! Input arguments
  INTEGER          , DIMENSION(12)            , INTENT(IN ) :: SEED
  DOUBLE PRECISION                            , INTENT(IN ) :: NHOPS
  INTEGER                                     , INTENT(IN ) :: NFRAMES, NDIMERS
  DOUBLE PRECISION, DIMENSION(NFRAMES,NDIMERS), INTENT(IN ) :: RATES, DX, DY, DZ
  INTEGER         , DIMENSION(2      ,NDIMERS), INTENT(IN ) :: CON
! Ouput arguments
  DOUBLE PRECISION, DIMENSION(NFRAMES        ), INTENT(OUT) :: DISTX, DISTY, DISTZ, TIME
  DOUBLE PRECISION, DIMENSION(NFRAMES,NDIMERS), INTENT(OUT) :: NHOP
! Function variables
  INTEGER                                     :: SEED_SIZE, FRAME, HOP, SITE, NSITES
  LOGICAL         , DIMENSION(        NDIMERS):: PATHWAYS
  DOUBLE PRECISION                            :: RANDOM_SMOL
  DOUBLE PRECISION, DIMENSION(        NDIMERS):: ALEAT, TRANSITION_TIMES
  INTEGER         , DIMENSION(1)              :: LOC


! Initialization of the random seed
  CALL RANDOM_SEED(SIZE=SEED_SIZE        )
  CALL RANDOM_SEED(PUT =SEED(1:SEED_SIZE))

! Determination of the number of sites
  NSITES = MAXVAL(CON)

! Initializing the output arguments
  DISTX=0.0E0
  DISTY=0.0E0
  DISTZ=0.0E0
  TIME =0.0E0
  NHOP =0.0D0
! Starting the simulations
  DO FRAME = 1 , NFRAMES
! Chosing randomly the starting molecule. (Random integer from 1 to NSITE)
    CALL RANDOM_NUMBER(RANDOM_SMOL)
    SITE = INT(RANDOM_SMOL*(NSITES)+1)
    HOP=0
    DO WHILE(HOP.LT.NHOPS)
! Chosing random numbers to weight the transition times. (The WHERE statement is use to pass from range [0,1[ to range ]0,1])
      CALL RANDOM_NUMBER(ALEAT)
      WHERE(Aleat.EQ.0.0)
        ALEAT=1.0
      END WHERE
! Finding all the pathways allowing the charge to escape the site.
      PATHWAYS = CON(1,:).EQ.SITE
! Caculating the transition times only for the pathways allowed
      WHERE(PATHWAYS)
        TRANSITION_TIMES = -(1/RATES(FRAME,:))*LOG(ALEAT)
      END WHERE
! Locating the fastest event
      LOC=MINLOC(TRANSITION_TIMES, MASK=PATHWAYS)
! Updating the distance and the time of the simulation
      DISTX(FRAME)=DISTX(FRAME)+DX(FRAME,LOC(1))
      DISTY(FRAME)=DISTY(FRAME)+DY(FRAME,LOC(1))
      DISTZ(FRAME)=DISTZ(FRAME)+DZ(FRAME,LOC(1))
      TIME (FRAME)=TIME (FRAME)+TRANSITION_TIMES(LOC(1))
! Updating the number of hopping event along each pathways
      NHOP (FRAME,LOC(1))=NHOP (FRAME,LOC(1))+1
! Updating the position of the charge
      SITE=CON(2,LOC(1))
      HOP=HOP+1
    END DO
  END DO

END SUBROUTINE FRM
