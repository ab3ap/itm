C     This program uses NTIA's itm.for to solve the problem presented in
C     section 7.2 of "A Guide to the Use of the ITS Irregular Terrain Model
C     in the Area Prediction Mode," Hufford, Longley, Kissick, NTIA Report
C     82-100, 1982.  Equation numbers in comments below reference those in
C     the report.
C
C     Mike Markowski, mike.ab3ap@gmail.com
C     June 2023

      subroutine coverage(res, Rn, s, rho, r)
      real, dimension(2,2), intent(in) :: res
      real, intent(in) :: Rn, s
      real, intent(out) :: rho, r

C     Linearly interpolate results to find percent coverage.
C
C     res[] is a 2 row list of lists.  Each row is [x, Rx], where x is distance
C     from desired antenna, and Rx is the D/U ratio at x.  The input res looks
C     like
C
C     [[x0, Rx0],
C      [x1, Rx1]]
C
C     where Rn is in interval [Rx0,Rx1] (or [Rx1,Rx0] if Rx1 < Rx0).
C     Rn's fraction f in [Rx0,Rx1] is used
C     to interpolate a value for r in [x0,x1].

C     Inputs:
C       res (float[2][2]): two lines of dist and D/U results.
C       Rn (float): D/U value to interpolate between res[:,1] values.
C       s (float): km, separation between broadcast antennas.
C     Outputs:
C       rho (float): percent interference-free coverage.
C       r (float): km, service range of interference coverage.

      pi = acos(-1.)
      span = res(2,2) - res(1,2)
      f = (Rn - res(1,2))/span ! Rx fraction between Rx's (D/U's).
C     Solve Eq 11 and 17 D/U ratio for x=r, interference-free range.
      span = res(2,1) - res(1,1)
      r = res(1,1) + f*span    ! Interpolate by f.
      rho = 100*2*pi/sqrt(3.)*(r/s)**2 ! Eq 12. % coverage.
      end subroutine

      subroutine deviation(dist_km, qT, qL, Ao, YT, YL, pl_dB)
      real, intent(in) :: dist_km, qT, qL
      real, intent(out) :: Ao, YT, YL

C     After itm.lrProp() has been called, this subroutine is called with
C     2 of the 3 quantiles of attenuation so that deviations from medians are
C     calculated and returned.  Eq 19 is implemented to generate quantiles
C     of individual attenuation.  Service range is based on time and location
C     and is why this subroutine deals with them.  

C     The definition of quantiles of observations are reversed here.  Quantiles
C     of observations are values exceeded for given fractions of time, locations
C     and situation.  See discussion under Eq 3 in 1982 report.  If parameters
C     qT and qL from the calling code represent (1-qT) and (1-qL), as in the
C     case of undesired signals, probabilities are that attenuation does not
C     exceed aVar(qT,qL,qS) for >= 1-qL locations for >= 1-qT of the time.

      common/prop/kwx,aref,mdp,dist,hg(2),wn,dh,ens,gme,zgnd,he(2),dl(2)
     &,the(2)
      common/propv/lvar,sgc,mdvar,klim

      zL = qerfi(qL)
      zT = qerfi(qT)
      zS = 0. ! qerfi(qS=0.5) = 0. qS=0.5 always for interference probs.

C     Make ITM recalculate distance based vals.
      lvar = max(1, lvar)
      call lrProp(dist_km*1e3) ! Longley-Rice propagation calculation.

C     Find time deviation from median for given location.
C     At least qL locations & qT times, atten does not exceed aVar(qT,qL,qS).
      YT1 = aVar(zT, zL, zS) ! Choose T after fixing L and S.
      YT0 = aVar(0., zL, zS) ! Time median.

C     Find location deviation from median.
C     For at least qL of locations, attenuation does not exceed aVar().
      YL1 = aVar(0., zL, zS) ! Choose L after fixing S5
      YL0 = aVar(0., 0., zS) ! Location median.

      YT = YT1 - YT0 ! Deviation from time median.
      YL = YL1 - YL0 ! Deviation from location median.

      Ao = aVar(0., 0., 0.) ! dB, overall median of attenuation.

      mdvarOld = mdvar
      mdvar = 3 ! Broadcast without interference mode.
      lvar = max(lvar, 4) ! New mdvar
      aConf = aVar(0., 0., zT) ! Location median.
      mdvar = 20 + 3 ! Broadcast with interference mode.
      lvar = max(lvar, 4) ! New mdvar
      pl_dB = aConf + 8.685890*log(2*wn*dist_km*1e3)
      end subroutine

      real function du(distD, distU, qT, qL, plD, plU)
      real, intent(in) :: distD, distU, qT, qL

C     Find the D/U ratio for a given distance and quantiles of attenuation.

C     Inputs:
C       distD (float): km, distance receiver is from desired antenna.
C       distU (float): km, distance receiver is from undesired antenna.
C       qT (float): 0 to 1, quantile of time.
C       qL (float): 0 to 1, quantile of location.
C     Output:
C       Rx (float): dB, desired/undesired ratio at distance distD km.

C     Find Rxo = D/U for closer offset stations.
C     Quantiles of attenuation for time and location.
      call deviation(distD,   qT,   qL, AoD, YTD, YLD, plD) ! D desired.
      call deviation(distU, 1-qT, 1-qL, AoU, YTU, YLU, plU) ! U undesired.
      YTR = sign(1., 0.5-qT) * sqrt(YTD**2. + YTU**2.) ! Time.
      YLR = sign(1., 0.5-qL) * sqrt(YLD**2. + YLU**2.) ! Location.
C     Eq 17, D/U desired-to-undesired ratio.
      du = 20*log10(distU/distD) - AoD + AoU + YLR + YTR
      end

      subroutine show(tbl, nRows)
      real, dimension(1000,6), intent(in) :: tbl ! Array, n x 6.
      integer, intent(in) :: nRows

      i = 1
      do j=1,nRows
        if (tbl(j,2) > tbl(i,2)) then
            i = j
        endif
      end do
C     Print table of results.
      rPrev = tbl(1,6) ! First range value.
      write(*,fmt="(a8 a7 a8 a8 a8 a8)") "s km", "cov %", "r km",
     +  "plD dB", "plU dB", "D/U dB"
      do n=1,nRows
        if (n .eq. i) then
          write(*,fmt='(">")',advance="no") ! Optimal s line, max coverage.
        else
          write(*,fmt='(" ")',advance="no")
        end if
        write(*,fmt="(2f7.1,4f8.1)", advance="no")
     +    tbl(n,1),tbl(n,2),tbl(n,3),tbl(n,4),tbl(n,5),tbl(n,6)
        if (rPrev .eq. tbl(n,6)) then
          write(*,fmt='(" ")')
        else
          write(*,fmt='("*")') ! Indicator when D/U changes.
        end if
        rPrev = tbl(n,6)
      end do
      write(*,fmt='(/,"Optimal separation s = ",f5.1," km.")') tbl(i,1)
      write(*,fmt='("Service range      r = ",f5.1," km.")') tbl(i,3)
      write(*,fmt='("Percent coverage   ",f5.1,"%.")') tbl(i,2)
      write(*,*) 
      end subroutine

      subroutine tv(tblN, r)
      real, dimension(1000, 6), intent(out) :: tblN
      integer, intent(out) :: r

      common/prop/kwx,aref,mdp,dist,hg(2),wn,dh,ens,gme,zgnd,he(2),dl(2)
     &,the(2)
      common/propv/lvar,sgc,mdvar,klim
      dimension kst(2)

      real, dimension(1000, 6) :: tblO
      real, dimension(2, 2) :: res
      real, dimension(2) :: cov
      integer s,x

C     Implement problem 7.2 Optimum Television Station Separation, p. 42 in 
C     "A Guide to the Use of the ITS Irregular Terrain Model in the Area
C     Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.

C     Set up values as given in report.
C     Common block vars.
c     prop
      kwx = 0        ! No errors initially.
      hg(1) = 305.   ! m, transmit antenna structure height.
      hg(2) = 9.     ! m, receive antenna structure height.
      dh = 90.       ! m, roughness factor for hilly terrain.
      ens = 301.     ! N-units, continental temperate refractivity.

c     propv
      klim = 5       ! Continental temperate.
      mdvar = 20 + 3 ! Interference + broadcast.
      lvar = 5

      qT = 0.9      ! Quantile of time, Grade A tv service problem.
      qL = 0.7      ! Quantile of location, Grade A tv service problem.
      f_MHz = 193.  ! TV channel 10.
      pol = 0       ! Horizontal antenna polarization.
      eps = 15.     ! F/m, relative ground permittivity.
      sgm = 0.005   ! S/m, relative ground conductivity.
      zSys = 0.     ! m, system elevation above sea level (0 sets Ns to N0).
      eN0 = 301.    ! N-units, sea level refractivity.

C     Required desired/undesired ratios:
      Ro = 28. ! Offset freqs: 28 dB
      R1 = 45. ! Non-offset f: 45 dB.

      lvar = 5
      call qlrps(f_MHz, zSys, eN0, pol, eps, sgm)
      kst(1) = 0 ! Random
      kst(2) = 0 ! Random
      call qlra(kst, klim, mdvar) ! Siting not given in report.

C     Rx=D/U for closest ring of offset freq antennas.
      i = 1
      r = 1
      do s=40,400,5 ! km, broadcast antenna separation.
        do x=1,s ! km, rcvr dist from desired tx antenna.
          i = 1 + mod(i, 2) ! Alternate 1 and 2 index.
          Rx = du(real(x), real(s-x), qT, qL, plD, plU) ! D/U ratio.
          res(i,1) = x
          res(i,2) = Rx ! New results.
          if (Rx < Ro) then ! Stop when Rx=D/U falls below Ro.
            call coverage(res, Ro, real(s), rho, rng) ! Solve Eq 11 and 12.
            ! Build table of results.
            tblO(r,1) = s
            tblO(r,2) = rho
            tblO(r,3) = rng
            tblO(r,4) = plD
            tblO(r,5) = plU
            tblO(r,6) = Ro
            r = r + 1
            exit ! Found Ro for this s, stop looking.
          end if
        end do
      end do

C     Rx=D/U for closest ring of offset freq antennas.
      i = 1
      r = 1
      do s=40,400,5 ! km, broadcast antenna separation.
        do x=1,2*s ! km, rcvr dist from desired tx antenna.
          i = 1 + mod(i, 2) ! Alternate 1 and 2 index.
          Rx = du(real(x), sqrt(3.)*s-x, qT, qL, plD, plU) ! D/U ratio.
          res(i,1) = x
          res(i,2) = Rx ! New results.
          if (Rx < R1) then ! Stop when Rx=D/U falls below R1.
            call coverage(res, R1, real(s), rho, rng) ! Solve Eq 11 and 12.
            ! Build table of results.
            tblN(r,1) = s
            tblN(r,2) = rho
            tblN(r,3) = rng
            tblN(r,4) = plD
            tblN(r,5) = plU
            tblN(r,6) = R1
            r = r + 1
            exit ! Found R1 for this s, stop looking.
          end if
        end do
      end do
      r = r - 1
C     Save row i from whichever table has shortest r, most interfered signal.
C     Each tbl row is [separation, coverage, range, D/U].
      do i=1,r
        if (tblO(i,3) < tblN(i,3)) then
          tblN(i,1) = tblO(i,1)
          tblN(i,2) = tblO(i,2)
          tblN(i,3) = tblO(i,3)
          tblN(i,4) = tblO(i,4)
          tblN(i,5) = tblO(i,5)
          tblN(i,6) = tblO(i,6)
        end if
      end do

      return
      end subroutine

      program tvpack
C     Main program
      real, dimension(1000, 6) :: tbl

      call tv(tbl, n)
      call show(tbl, n)
      end program
