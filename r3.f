c     This program uses itm.for to solve the problem presented in section 7.3
c     of "A Guide to the Use of the ITS Irregular Terrain Model in the Area
c     Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.
c
c     Mike Markowski, mike.ab3ap@gmail.com
c     June 2023

      subroutine predict()

      common/prop/kwx,aref,mdp,dist,hg(2),wn,dh,ens,gme,zgnd,
     x   he(2),dl(2),the(2)
         complex zgnd
      common/propv/lvar,sgc,mdvar,klim
      common/propa/dlsa,dx,ael,ak1,ak2,aed,emd,aes,ems,dls(2),dla,tha
      dimension kst(2)
      integer dist_km

c     Implement problem 7.2 Optimum Television Station Separation, p. 42 in 
c     "A Guide to the Use of the ITS Irregular Terrain Model in the Area
c     Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.

c     Set up values for R3 data set as given in section 7.3 of report.
c     prop
      dh = 126.    ! m, roughness factor for hilly terrain.
      hg(1) = 275. ! m, transmit antenna structure height.
      hg(2) = 6.6  ! m, receive antenna structure height.

c     propv
      lvar = 5
      mdvar = 3    ! Broadcast mode.

      en0 = 300.   ! N-units, sea level refractivity (yields Ns=250).
      eps = 15.    ! F/m, relative ground permittivity.  Average ground.
      f_MHz = 410. ! UHF tv-like.
      ipol = 0     ! Horizontal antenna polarization.
      sgm = 0.005  ! S/m, relative ground conductivity.  Average ground.
      zs = 1700.   ! m, system elevation above sea level.
      zT = 0.      ! qT = 0.5, no time variability.

      kwx = 0
      call qlrps(f_MHz, zs, en0, ipol, eps, sgm)
      if (.true.) then ! Recreate figure 9.
        kst(1) = 0
        kst(2) = 0
      else ! Recreate figure 12.
        kst(1) = 2
        kst(2) = 2
      endif
      call qlra(kst, -1, -1)

c     All curves exactly match report except A(0.1,0.1).
      do i=1,5
        if (i .eq. 1) then
          qL = 0.9
          qC = 0.9
        elseif (i .eq. 2) then
          qL = 0.5
          qC = 0.9
        elseif (i .eq. 3) then
          qL = 0.5
          qC = 0.5
        elseif (i .eq. 4) then
          qL = 0.5
          qC = 0.1
        elseif (i .eq. 5) then
          qL = 0.1
          qC = 0.1
        endif
 100    format(2hA(f4.1, f4.1,1h))
        write(*, 100) qL, qC
        zL = qerfi(qL)   ! Quantile of location.
        zC = qerfi(qC)   ! Quantile of confidence (situation).
        do dist_km=1,100 ! Distance from cliff-top transmitter.
          lvar = max(1, lvar)
          call lrProp(dist_km*1e3)   ! Longley-Rice propagation calculation.
          aRef_dB = aVar(zT, zL, zC) ! Variability of attenuation.
c200      format(i6, f8.1)
c         write(*, 200) dist_km, aRef_dB
 200      format(f8.1,1x ,$)
          write(*, 200) aRef_dB
          if (aRef_dB .gt. 40) then  ! Match max atten in Figures 9 and 12.
            exit
          endif
        enddo
        write(*,*)
      enddo
      endsubroutine

      call predict()
      end
