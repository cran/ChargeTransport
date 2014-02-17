SUBROUTINE BKL(SEED,NHOPS,NFRAMES,NDIMERS,RATES,DX,DY,DZ,CON, &
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
  INTEGER                                     :: SEED_SIZE, FRAME, HOP, SITE, NSITES, I
  LOGICAL         , DIMENSION(        NDIMERS):: PATHWAYS
  DOUBLE PRECISION                            :: RANDOM_SMOL, R
  INTEGER         , DIMENSION(1)              :: LOC
  DOUBLE PRECISION, DIMENSION(        NDIMERS):: CUMSUM


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
! Finding all the pathways allowing the charge to escape the site.
      PATHWAYS = CON(1,:).EQ.SITE
! Calculation of the cumulative sum of the rates for all the pathways allowing to escape the site
      CUMSUM = 0.0E0
      DO I=1,NDIMERS
        IF(PATHWAYS(I)) THEN
          WHERE(PATHWAYS)
            CUMSUM = CUMSUM + RATES(FRAME,I)
          END WHERE
          PATHWAYS(I)=.FALSE.
        END IF
      END DO
! Generating a random number to select the pathway to escape. (0.0 is replace by 1.0)
      CALL RANDOM_NUMBER(R)
      IF(R.EQ.0.0) R=1.0
! Finding the pathway to escape. (First rate for witch: cumsum(i) > r.k_tot )
      LOC=MINLOC(CUMSUM,MASK=CUMSUM.GT.R*MAXVAL(CUMSUM))
! Updating the distance and the time of the simulation
      DISTX(FRAME)=DISTX(FRAME)+DX(FRAME,LOC(1))
      DISTY(FRAME)=DISTY(FRAME)+DY(FRAME,LOC(1))
      DISTZ(FRAME)=DISTZ(FRAME)+DZ(FRAME,LOC(1))
      TIME (FRAME)=TIME (FRAME)-(1/MAXVAL(CUMSUM))*LOG(R)
! Updating the number of hopping event along each pathway
      NHOP (FRAME,LOC(1))=NHOP (FRAME,LOC(1))+1
! Updating the position of the charge
      SITE=CON(2,LOC(1))
      HOP=HOP+1
      END DO
    END DO

END SUBROUTINE BKL
