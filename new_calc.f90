PROGRAM calculator
! by Bucky Badger
! This program calculates meteorological variables.
IMPLICIT NONE

! These are our original variables
REAL :: spd, dir, rh

! These are our computed variables
REAL :: uwind, vwind

! This is our constant
REAL, PARAMETER :: pi=3.14159

! This is our loop variable
INTEGER :: i

! These are our matrices
REAL :: uwindmat(1453), vwindmat(1453), rhmat(1453)

! These are our masks
LOGICAL, DIMENSION(1453) :: ugt0all, ult0all, ugt2all, ult2all
LOGICAL, DIMENSION(363) :: ugt0q1, ult0q1, ugt2q1, ult2q1
LOGICAL, DIMENSION(364) :: ugt0q3, ult0q3, ugt2q3, ult2q3

! These are our averages
REAL :: ugt0allavg, ult0allavg, ugt2allavg, ult2allavg, &
        ugt0q1avg, ult0q1avg, ugt2q1avg, ult2q1avg, &
        ugt0q3avg, ult0q3avg, ugt2q3avg, ult2q3avg

PRINT *, "Program is now starting."

! Read in data
OPEN(UNIT=10, FILE='data/obs.txt', ACTION='read', STATUS='old')
READ(10, *)  ! Skip the first line--column headers

! Write out data
OPEN(UNIT=20, FILE='data/obs_conv.txt', ACTION='write')
WRITE(20, *) "uwind[mph] vwind[mph] rh[%]"

DO i = 1, 1453
  READ(10, '(24X,F8.2,F8.0,F8.0)') rh, dir, spd
!  spd = spd * 1.151        ! this is mph
  CALL ktstomph(spd, spd)       ! this is mph

  dir = dir - 180.0        ! direction wind is going in degrees
  CALL deg2rad(dir, dir)
!  dir = dir * (pi / 180.0) ! direction wind is going in radians

  CALL decompose(spd, dir, uwind, vwind)
!  uwind = spd * SIN(dir)
!  vwind = spd * COS(dir)

  uwindmat(i) = uwind
  vwindmat(i) = vwind
  rhmat(i) = rh

  PRINT *, i, uwind, vwind, rh
  WRITE(20, *), uwind, vwind, rh   ! Write out converted data to file

END DO

CLOSE(UNIT=10)
CLOSE(UNIT=20)

! Construct masks
ugt0all = uwindmat > 0.0
ult0all = uwindmat < 0.0
ugt2all = uwindmat > 2.0
ult2all = uwindmat < 2.0
ugt0q1 = uwindmat(1:363) > 0.0
ult0q1 = uwindmat(1:363) < 0.0
ugt2q1 = uwindmat(1:363) > 2.0
ult2q1 = uwindmat(1:363) < 2.0
ugt0q3 = uwindmat(726:1089) > 0.0
ult0q3 = uwindmat(726:1089) < 0.0
ugt2q3 = uwindmat(726:1089) > 2.0
ult2q3 = uwindmat(726:1089) < 2.0

! Calculate averages of _RH_ using uwind masks
ugt0allavg = SUM(rhmat, MASK=ugt0all) / COUNT(ugt0all)
ult0allavg = SUM(rhmat, MASK=ult0all) / COUNT(ult0all)
ugt2allavg = SUM(rhmat, MASK=ugt2all) / COUNT(ugt2all)
ult2allavg = SUM(rhmat, MASK=ult2all) / COUNT(ult2all)
  ! Note for the sliced matrices, while we indexed the masks before,
  !  we still must index the rhmatrix to match the same area.
ugt0q1avg = SUM(rhmat(1:363), MASK=ugt0q1) / COUNT(ugt0q1)
ult0q1avg = SUM(rhmat(1:363), MASK=ult0q1) / COUNT(ult0q1)
ugt2q1avg = SUM(rhmat(1:363), MASK=ugt2q1) / COUNT(ugt2q1)
ult2q1avg = SUM(rhmat(1:363), MASK=ult2q1) / COUNT(ult2q1)
ugt0q3avg = SUM(rhmat(726:1089), MASK=ugt0q3) / COUNT(ugt0q3)
ult0q3avg = SUM(rhmat(726:1089), MASK=ult0q3) / COUNT(ult0q3)
ugt2q3avg = SUM(rhmat(726:1089), MASK=ugt2q3) / COUNT(ugt2q3)
ult2q3avg = SUM(rhmat(726:1089), MASK=ult2q3) / COUNT(ult2q3)

OPEN(UNIT=30, FILE='data/obs_avg.txt', ACTION='write')
WRITE(30, *) "Here are the RH averages for each of the 12 conditions:"
WRITE(30, *) ugt0allavg, " : for u > 0 over the full year,"
WRITE(30, *) ult0allavg, " : for u < 0 over the full year,"
WRITE(30, *) ugt2allavg, " : for u > 2 over the full year,"
WRITE(30, *) ult2allavg, " : for u < 2 over the full year,"
WRITE(30, *) ugt0q1avg, " : for u > 0 over the first quarter,"
WRITE(30, *) ult0q1avg, " : for u < 0 over the first quarter,"
WRITE(30, *) ugt2q1avg, " : for u > 2 over the first quarter,"
WRITE(30, *) ult2q1avg, " : for u < 2 over the first quarter,"
WRITE(30, *) ugt0q3avg, " : for u > 0 over the third quarter,"
WRITE(30, *) ult0q3avg, " : for u < 0 over the third quarter,"
WRITE(30, *) ugt2q3avg, " : for u > 2 over the third quarter,"
WRITE(30, *) ult2q3avg, " : for u < 2 over the third quarter."
CLOSE(UNIT=30)

PRINT '(T20,A)', "Program has finished."
PRINT *, "Go Bucky!"

END PROGRAM calculator

SUBROUTINE ktstomph(inspeed, outspeed)
! by Bucky Badger
! This subroutine converts a speed in knots to miles per hour.
IMPLICIT NONE

! This is our input and our output variable.
REAL :: inspeed, outspeed

outspeed = inspeed * 1.151

END SUBROUTINE ktstomph

SUBROUTINE deg2rad(degrees, radians)
! by Bucky Badger
! This subroutine converts a degree into radians.
IMPLICIT NONE

! This is our input and our output variable.
REAL :: degrees, radians
REAL, PARAMETER :: pi=3.14159

radians = degrees * (pi / 180.0)

END SUBROUTINE

SUBROUTINE decompose(speed, direction, uwind, vwind)
! by Bucky Badger
! This subroutine splits up a wind speed and direction (going towards)
!  into the zonal and meridional components. This is a simple vector
!  projection onto the x and y axes.
IMPLICIT NONE

! These are the input variables
REAL :: speed, direction

! These are the output variables
REAL :: uwind, vwind

uwind = speed * SIN(direction)
vwind = speed * COS(direction)

END SUBROUTINE decompose
