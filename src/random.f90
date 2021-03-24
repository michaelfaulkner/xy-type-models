!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR.                    C
!  Must be initialized with randinit() before use.                     C
!                                                                      C
!  See F. James, Comp. Phys. Comm. 60, 329 (1990), or                  C
!  G. Marsaglia et al., Stat. Prob. Lett. 9, 35 (1990).                C
!                                                                      C
!----------------------------------------------------------------------C


!----------------------------------------------------------------------C
!                                                                      C
! This is the initialization routine RMARIN for the random number      C
!     generator RANMAR                                                 C
!                                                                      C
! NOTE: The seed variables can have values between:  0 <= IJ <= 31328  C
!                                                    0 <= KL <= 30081  C
!----------------------------------------------------------------------C

      SUBROUTINE randinit(seed)
      IMPLICIT NONE
      INTEGER seed
      INTEGER ij,kl, i,j,k,l, ii,jj, m
      double precision :: s,t
      INTEGER Maxseed
      PARAMETER (Maxseed = 900000000)
      double precision :: u(97), c, cd, cm
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec

      seed = mod(seed,Maxseed)
      ij = seed / 30082
      kl = seed - (30082 * ij)
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      RETURN
      END



!----------------------------------------------------------------------C
!                                                                      C
!  Lagged Fibonacci random number generator RANMAR().                  C
!                                                                      C
!----------------------------------------------------------------------C

      FUNCTION rand()
      IMPLICIT NONE
      double precision :: u(97), c, cd, cm, uni, rand
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec

      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      rand = uni
      RETURN
      END
