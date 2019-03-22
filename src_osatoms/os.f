      PROGRAM OS
C this program, OS.F, substitutes the old Opacity Sampling Generation 
C Program (OSGEP), with better line profiles, better programming, etc.
C     PARAMETER (NTEMP=12,MT1=NTEMP+1,nso=6000)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'atomparameter.inc'
      real*4 time_0,time_1,time_2
C     DIMENSION sumop_m(ntemp),sumop_l(ntemp),
      DIMENSION sumop_m(ntemp),opt(ntemp),dadt(ntemp),w(3,ntemp),
     & tmold(ntemp)
     & ,sumop_h2o(0:13),nlin_h2o(0:13),phcn(2)
      character MOLID*4,osfil*60,adum*66

      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      COMMON /COPAC/HALF,sig_nu(ntemp),hckt(ntemp)
     *    ,sqr_pi_m,sig_sqr,sumop_c(ntemp),a17,n16,n17,npr,j
      COMMON/CISO/vkmssp,wgtmola,rel_iso_c(3),pc3(3)
     & ,reliso(5),kiso,nhalf,jderiv,kiso_c
     & ,linesc3_direct,lineshcn_direct,linesc2h2_direct
      COMMON/CMOLSP/MOLID,OSFIL
      COMMON/CSWMX/sgf,sumop_l(ntemp)

      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      common /cph/ph(ndp)
      CHARACTER ATMOS*45

      DATA rmol /8.3143d7/    !Molar gas constant R in erg/(mol*K)
      common /cprespp/prespp(ndp,nspec)           !partial pressures in dyn/cm^2
      common /cabink/abink(ndp,nspec)
      DIMENSION pp_sum(ndp),ropp(ndp)
      common /cwmol/wmol(nspec)

       namelist /vald_inp/ lvald, elm_pick, ion_pick, 
     &     wb, we, lgstel, nsort, ktsort, scalemet
       character elm_pick*3
       common /cvald/ wb, we, lvald, ion_pick, elm_pick
       character elm*2, nameatom*2
       common /cvaldline/ wgtatom, vv0, vexs, vgf,
     &    gam_rad, gam_stark, gam_waal, fac_lande, 
     &    conv_gstel(ntemp),
     &    jatom, jg, ion, lgstel,
     &    elm, nameatom

      namelist /nprofile/lmarcs,vkmssp,nhalf,profile,atmos_voigt
      common/cvoigt/ pgas(mt1),p_el(mt1),p_h(mt1),hwda(ntemp)
     &              ,pemu(mt1),pp(mt1,nspec)
     &              ,profile,atmos_voigt
      character profile*5,atmos_voigt*45
      COMMON /CCHROMOSPHERE/ nchrom,lchrom,itmin

      NAMELIST /INPUTOSMOL/ MOLID, KTEMP, TMOL, NWNOS
     &   ,VKMS, KISO, RELISO,l_per_stellar, lchrom
      character h2olistname*60
      NAMELIST /H2OTEST/ lh2otest, ktemp_h2o, h2o_lim, h2olistname
      data osfil/'os.dat'/
      common /comdata/ BDA,pi,c_vel
      character namesort*2
      dimension  stsort(nso), namesort(nso), ionsort(nso), gfsort(nso),
     &           exssort(nso), v0sort(nso), ppsort(nso)

      BDA = 8.31327E7
      pi =  3.141593
      c_vel = 2.99792E10
C C_VEL: velocity of light in cm/s,
C BOLTZ: Boltzmanns constant in erg*K-1, AMU: atom
C mass unit on 12C scale, BDA=BOLTZ/AMU=1.38042E-16/1.6605E-24

C Open relevant units
C
      OPEN (UNIT=12,FILE='osat.input',STATUS='old')
      READ (12,INPUTOSMOL)
      write(6,612) molid,ktemp,(tmol(kt),kt=1,ktemp)
612   format(' Molid: ',a4,' OS being computed at',i3,' temperatures:',
     &     /12f7.0)
C     WRITE (6,INPUTOSMOL)
      READ (12,H2OTEST)
C     WRITE (6,H2OTEST)
      OPEN (UNIT=18,FILE='osat.output'
     *                                          ,STATUS='unknown')
C
C START TIMETAKING
      call etime(time_0)

C INITIATE OPACITY SAMPLING COMPUTATIONS:
C
      call oswavenumbers
      WRITE(6,812) OSFIL
812   FORMAT (' The OS file will be written on file ',A50)
      OPEN (UNIT=8,FILE=osfil,STATUS='unknown')
      read(12,vald_inp)
C     write(18,vald_inp)
      if(nsort.gt.nso) stop ' increase dimension nso'
      if(lvald.eq.0) lgstel = 0     ! we can only compute per stellar for atoms
C                       in particular conv_gstel(kt) is only defined for atoms and
C                       the opacities will therefore be 0 for molecules if lgstel=1
      l_per_stellar = lgstel
      read(12,nprofile)
C     write(18,nprofile)

      if(profile.eq.'voigt' .or. lvald.eq.1) then

      CALL GEM_INIT        !read gem names and indexes etc
      call molwgt   !compute molecular weights based on molecular name and atomic wgt
      atmos = atmos_voigt
      if(lmarcs.eq.0) then
        call tpgread
      else
        CALL MODEL
      end if


       DO 3151 kd=1,NTAU
       ropp(kd) = 0.
       pp_sum(kd) = 0.
       pe_gem = pg(kd) * abink(kd,1) / (1.d0 - abink(kd,1))
       p_particles = pg(kd) + pe_gem             !pg(input)+pe(computed) in dyn/cm^2
       DO 3153 km=1,nspec
       PRESPP(kd,km) = p_particles * abink(kd,km)  !(pg+pe)*rel.pp = pp in dyn/cm^2
       if(km.gt.1) then
          ropp(kd) = ropp(kd) + prespp(kd,km) * wmol(km)  !sum([dyn/cm2*g/mol]
          pp_sum(kd) = prespp(kd,km) + pp_sum(kd)
       end if
3153   continue
          ph(kd) = prespp(kd,2)
          emu(kd) = ropp(kd)/pg(kd)    !mean molecular weight stellar in g*/mol
          ropp(kd) = ropp(kd) / (rmol * t(kd)) !sum[dyn/cm2*g/mol]/[RT]=sum[g_pp/cm3]
C       if(kd.eq.1 .or. kd.eq.21 .or. kd.eq.41) then
C       write(6,3158) (dlog10(max(1.d-40,prespp(kd,km))),km=1,15)
       modwrite = 0
       if(modwrite.eq.1) write(6,3150)
     & kd,t(kd),pg(kd),pe(kd),ro(kd)
     & ,pp_sum(kd),pe_gem,ropp(kd),ro(kd)/ropp(kd)
     & ,emu(kd)
     & ,log10(prespp(kd,2)/prespp(kd,50))
     & ,prespp(kd,2)+prespp(kd,3)
     & +prespp(kd,50)+prespp(kd,51)
     & +prespp(kd,52)
C       end if
3150   format(i3,f9.2,1p3e9.2,3x,3e9.2,0p3f6.3,1p3e9.2)
3151   continue
3158   format('OS',10f7.2/,3x,10f7.2)


C for voigt profile as well as for computation of a cm2/g* general
C atomic OS we need to know the T,Pg values the computation is to 
C be done for. We therefore read here interpolate T, Pg, EMU to
C the relevant OS temperatures, from the T,Pg form a model atmosphere
C or other T,Pg arrays and the cmemical equilibrium computed here or
C in the model:

      if(profile.eq.'voigt') then
      call voigt_init
      endif


      end if


      if(lvald.eq.1) then
C----------------------------------------------------------------------
C                                                                     V
       call vald_open

C                                         VALD atomic lines           |
C----------------------------------------------------------------------
      end if
C----------------------------------------------------------------------
C  Abs.coef. from a molecule (necessary to call even for atoms only   | 
C  in which all mol_specific will do is to read namelist 'sampling'   V

      call mol_specific
C                                             molecular lines         |
C----------------------------------------------------------------------



      if (ktemp.gt.ntemp) then
	write(18,*) ' ktemp > ntemp; increase dimension ntemp'
	go to 900
      endif
      vkms = vkmssp     !for header in os.dat (just for identification)
      write(18,INPUTOSMOL)

      n16 = 0
      n17 = 0
      a17 = 0.
      sgf = 5.33129e+11      ! conversion from gf to cm/mol

C Calculate data for the doppler motions and the associated line width
	 
C If more than one atomic species is computed, wgtmola will be determined
C for each line, and sig_nu(it) computed for each line, otherwise (only one
C kind of atoms, e.g. Ca) wgtmola is fixed and determined above via the 
C call to vald_open.
	 if(lvald.eq.1 .and. elm_pick.eq.'999') wgtmola = 1.
	 ak2_c = 2.0*BDA/(c_vel**2) ! 2k/c^2
	 ak2_mc = (2.0*BDA/wgtmola)/(c_vel**2) ! 2k/mc^2
	 xic = vkmssp*1.0E5/c_vel ! XI in cm/s, vkmssp in km/s
	 sqr_pi_m = 1.0/sqrt(pi)

	 do it=1,ktemp
	    sumop_c(it) = 0.
	    sumop_m(it) = 0.
	    sumop_l(it) = 0.
	    sig_nu2 = ak2_mc*tmol(it)+xic**2 ! (2kT/m + xi^2)/c^2
	    sig_nu(it) = sqrt(sig_nu2)
	    HCKT(it) = 1.4388/TMOL(IT)
	 enddo

	 half = float(nhalf)    ! max # Gauss halfwidths to follow profile

      call etime(time_1)
      write(6,*) ' time for initiation: ',time_1-time_0
      write(18,*)' time for initiation: ',time_1-time_0
      write(18,14) nwnos,wnos(1),wnos(nwnos)
14    format(' nwnos,wnos(1),wnos(nwnos) = ',i6,4f7.0)
      write(18,*)'Temperature and sig_nu(t):'
      write(18,141) (tmol(it),it=1,ktemp)
141   format(12f6.0)
      write(18,142) (sig_nu(it),it=1,ktemp)
142   format(6f12.8)

C
C computation of the vibrational partition function for the temperature
C for which the Opacity Sampling will be done

C temporarily (1998,05) we read direcvtly from the water line list
C
      rewind(66)
      if (MOLID.eq.'C2H2' .and. kiso_c.gt.1) then
	   p12c = rel_iso_c(1)**2
      end if

       pc3(1) = (rel_iso_c(1))**3
       pc3(2) = (1.-pc3(1)) / 3.0  !there are 3 different 13C12C2 symmetries
       pc3(3) = pc3(2)

       phcn(1) = rel_iso_c(1)
       phcn(2) = rel_iso_c(2)

      if (MOLID .eq. 'H2Os') then
	    open (unit=65,status='old',readonly,
     &      file='/ast/p2/uffegj/2013/h2o/scan_h2o.dat')
c    &      file=h2olistname)
C    &      file='/ste1/uffegj/h2o/spectr1mil.dat')
C    &      file='/ste1/uffegj/h2o/list.sort')
C    &      file='/ste3/uffegj/h2o/h2o.dat')
C     write(6,*) ' I opened unit65=ste3/uffegj/h2o/h2o.dat'
      write(6,*) ' I opened unit65=ste4/uffegj/scan_h2o.dat'
      end if
      if (MOLID .eq. 'H2O ') then
	    open (unit=65,status='old',readonly,
     &      file='/ste1/uffegj/h2o/b2b2.lis_rr9_jmin_dk25')
C    &      file='/ste1/uffegj/h2o/b2b2.lis_rr')
	    nh2ofil = 1
      end if
       do 1998 kh2o=0,13
       nlin_h2o(kh2o) =  0
1998   sumop_h2o(kh2o) = 0.d0

C calculation of OS from Schwenke's 300 mill line list
      if (MOLID .eq. 'swmx') then
	    call sw_h2o
	    go to 119
      end if

      nr = 0
      nrlin = 0
      linesgfv0 = 0
      npick = 0
      do 1459 lsort = 1,nsort
      namesort(lsort) = 'xx'
      ionsort(lsort) = 0
      gfsort(lsort) = 0.d0
      exssort(lsort) = 0.d0
      v0sort(lsort) = 0.d0
      ppsort(lsort) = 0.d0
1459  stsort(lsort) = 0.d0

      do 110 i=1,900000000

      if (MOLID .eq. 'H2Os' .or. MOLID.eq.'H2O ') then
	read(65,*,end=109) v0,s0,exs
	go to 108
109     continue
	if (MOLID .eq. 'H2Os') go to 119      ! this is standard - reading from only one file
C                                               H2O ~ read directly from 3 vibnew.f output files
	if (nh2ofil.le.3) then
	    nh2ofil = nh2ofil + 1
	    close (unit=65)
	    if (nh2ofil.eq.4) go to 119
	    if (nh2ofil.eq.2) open (unit=65,status='old',readonly,
     &      file='/ste2/uffegj/h2o/a1a1.lis_rr9_jmin_dk25')
C    &      file='/ste1/uffegj/h2o/a1a1.lis_rr')
	    if (nh2ofil.eq.3) open (unit=65,status='old',readonly,
     &      file='/ste2/uffegj/h2o/a1b2.lis_rr9_jmin_dk25')
C    &      file='/ste1/uffegj/h2o/a1b2.lis_rr')
	    read(65,*,end=119) v0,s0,exs
	end if
108     continue

C      NAMELIST /H2OTEST/ lh2otest, ktemp_h2o, h2o_lim 
	 IF (lh2otest.eq.1) THEN
	 kth = ktemp_h2o
	 STH2O = s0 
     &     *EXP(-exs*HCKT(kth))*(1.-EXP(-v0*HCKT(kth)))/qvib(kth)
	 kh2o = max( min(log10(1.d-3*sth2o),0.0d0), -13.d0)
	 kh2o = -kh2o   !i.e, st>1.e2~0,  1.e1<st=<1.e2~1, ... st=<1.e-10~13
	 sumop_h2o(kh2o) = sumop_h2o(kh2o) + sth2o   !integrated abs.coef.(kth) in [km/mol]
	 nlin_h2o(kh2o) =  nlin_h2o(kh2o) + 1
	 if(sth2o.lt.h2o_lim) go to 110
	 END IF


	gf = s0*1.87572D-7             !conv from s0 km/mol to gf
      
      else if (MOLID .eq. 'C3  ' .and. linesc3_direct.eq.1) then

      read(66,*,end=119) v0,exs,gf,isoc3
       gf  = pc3(isoc3) * gf

      else if ( (MOLID .eq. 'HCN ' .or. MOLID .eq. 'HCNt')
     & .and. lineshcn_direct.eq.1) then

      IF (MOLID.eq.'HCN ') then 
          read(66,*,end=119) v0,exs,gf,isohcn
          gf  = phcn(isohcn) * gf
      else if (MOLID.eq.'HCNt') then 
          read(66,*,end=119) v0,exs,gf
          isohcn = 1
      END IF

      else if (lvald.eq.1) then
C     read(66,*,end=119) v0,exs,gf,wgtatom
      call vald_line(kadopt)
      if(kadopt.eq.0) go to  110
      if(kadopt.eq.999) go to  119
      gf = vgf
      if(scalemet.ne.1.) gf = gf * scalemet
      v0 = vv0
      exs = vexs
C      if(nrlin.le.10) then 
C           write(6,1162) nrlin,v0,1.d4/v0,exs,gf,wgtatom
C1162   format('Vald lines #/vo/mu/ex/gf/wgt:',
C     &                          i2,3f10.2,1pe10.2,0pf8.2)
C      end if
      if(elm_pick.eq.'999') wgtmola = wgtatom   !if ne 999 wgtmola const in ciso
      npick = npick + 1

      else

      read(66,*,end=119,err=1191) v0,exs,gf
	if (gf .le. 0.D0 .or. v0.le.0.d0) then
	   linesgfv0 = linesgfv0 + 1
	   go to 110
	end if

C       write(6,1098) lw,ws,gf,exs,iel,ivl,ivu,jl,ju,ibr
C1098   format(i3,f11.3,1pe11.3,0pf10.2,i2,2i3,2i4,i3)


      end if

      nrlin = nrlin + 1

      if (MOLID.eq.'C2H2' .and. kiso_c.gt.1) then
	gf = gf * p12c
	liso = 1
      end if 
129   continue      !compute more isotopes not listed

      absl = sgf * gf
C     kt2 = max(1,ktemp/2)
      kt2 = ktsort
      do 139 kt = 1,ktemp
      ST = absl
     &     *EXP(-exs*HCKT(kt))*(1.-EXP(-v0*HCKT(kt)))/qvib(kt)

         if(lvald.eq.1 .and. kt.eq.kt2) then
            do 1450 lsort = 1,nsort
              if(st*conv_gstel(kt2).ge.stsort(lsort)) then
                 do 1451 jsort=nsort,lsort+1,-1
                    stsort(jsort) = stsort(jsort-1)
                    namesort(jsort) = namesort(jsort-1)
                    ionsort(jsort) = ionsort(jsort-1)
                    gfsort(jsort) = gfsort(jsort-1)
                    exssort(jsort) = exssort(jsort-1)
                    v0sort(jsort) = v0sort(jsort-1)
                    ppsort(jsort) = ppsort(jsort-1)
1451             continue
                 stsort(lsort) = st*conv_gstel(kt2)
                 namesort(lsort) = nameatom
                 ionsort(lsort) = ion
                 gfsort(lsort) = gf
                 exssort(lsort) = exs
                 v0sort(lsort) = v0
                 ppsort(lsort) = pp(kt2,jg)/pgas(kt2)
                 go to 1452
              end if
1450        continue
1452        continue
         end if


      sumop_l(kt) = sumop_l(kt) + st
139   continue

C actually qvib is qv*qr so consistency with conv. to gf is ok.
C eg.for C2H2:   gf = 1.87572E-12*s0 * qh296 / exp(-exs*1.4388/296.) /
C    &     (1.-exp(-v0*1.4388/296.) ) * ana
C -- ana = N_A = Avogadro's number = 6.0221367e23 molecules/mol

      if (v0.lt.wnos(1)) go to 110
C     if (v0.gt.wnos(nwnos)) go to 119
C temp we go to 110 because h2o list is not ordered
      if (v0.gt.wnos(nwnos)) go to 110

      call range(v0,nw)   !this subr. identifies nw such that 
C                          wnos(nw) < v0 < wnos(nw+1)

      if (v0.lt.wnos(nw) .or. v0.gt.wnos(nw+1) ) then
	write(18,*) ' nw,nw+1,v0 ???  => stop! '
	write(18,138) nw,v0,wnos(nw),wnos(nw+1)
	go to 900
138     format (' i,nr,nw,v0,wnos{nw,nw+1}=',i5,2i5,3f11.3)
      else
	nr = nr + 1
	if(nr.le.5)write(18,138) i,nr,nw,v0,wnos(nw),wnos(nw+1)
      end if

C if several atoms, compute sig_nu(it) to common for profile of each individual line 
C if only one atom (e.g, elm_pick='Ca+') sig_nu(it) is computed above and only once
	    if(lvald.eq.1 .and. elm_pick.eq.'999') then
	    do 1388 it=1,ktemp
	    sig_nu2 = ak2_c/wgtmola*tmol(it)+xic**2 ! (2kT/m + xi^2)/c^2
	    sig_nu(it) = sqrt(sig_nu2)
1388        continue 
	    end if

	 CALL OPAC_COM(v0,gf,exs,nw)

      if (MOLID.eq.'C2H2' .and. kiso_c.gt.1) then
	 if (liso .eq. 2) go to 110
	   liso = liso + 1
	   gf = gf * (1. - p12c)/p12c
	   v0 = 0.997 * v0    !ws(13)=0.997ws(12), gf(13)=gf(12)
	   go to 129
      end if

	 go to 110
1191     continue
	 backspace(1)
	  read(66,1193) adum
	  write(6,1192) i,adum
1192      format(' read error for i=',i6,2x,a66)
1193      format(a66)

110      continue                      !read next line from input data

119      continue
	 write(6,*) ' we read ',nrlin,' lines from line list'
	 write(6,*) ' we called opac_com ',nr,' times'
	 write(6,*) n16,' lines contributed to opac'
	 write(6,*) ' on the average each line contributed to'
     *     ,float(n17)/float(n16),' opac frequency points'
	 write(6,*) ' the average prof factor was',a17/float(n17)
	 write(6,*) ' there were ',linesgfv0,' lines with gf,0 or v0,0'

         IF (lvald.eq.1) THEN
         if(nsort.gt.npick) then
            nsort = npick
            write(6,*) 
     &      'We cannot sort ',nsort,' lines as requested in input'
            write(6,*) 
     &      'Instead we sort ',npick,' lines only'
         end if
         write(6,1460) nsort, tmol(ktsort)
1460     format(i3,' strongest lines sorted by cm/g* at T=',f8.1)
         write(6,*)' v0[AA]  # elm ion cm/g* gf exs exs_eV pp/pg'
         write(6,1461) (1.d8/v0sort(j),j,namesort(j),ionsort(j),
     *      stsort(j),gfsort(j),exssort(j),exssort(j)*1.2398e-4,
     *      ppsort(j),j=1,nsort)
1461     format(0pf10.2,i3,1x,a2,i1,1p5e10.2)
         END IF

	 IF (lh2otest.eq.1) THEN
	 write(6,*) 
     *  ' h2o line intensity at t = ',ktemp_h2o,' in the ranges'
	 write(6,*) 
     * ' S>100km/mol, 10<S=<100, 1<S<=10, 0.1<S<=1, ... S<=1.e-10km/mol'
	 write(6,1999) 
     *  ( (4-kh2o,3-kh2o,nlin_h2o(kh2o),sumop_h2o(kh2o)), kh2o=0,13 )
C        kh2o = -kh2o   !i.e, st>1.e2~0,  1.e1<st=<1.e2~1, ... st=<1.e-10~13
1999     format(' Lines with 1.e'i4,' -- 1.e',i4,'km/mol:',
     *   ' N: 'i9,' intgr.abs:'1pe11.3,'km/mol')
	 END IF


	   nzf = 0
	 DO 319 JV=1,nwnos
	   OPSUM = 0.
	   DO 317 IT=1,KTEMP
	   OPSUM = OPSUM + OPACOS(IT,JV)
317        CONTINUE
	   if (opsum.eq.0.) go to 319
	   if (opsum.ne.0. .and. nzf.eq.0)
     &     nzf = jv    !first wavenumber with non-zero opacity
	   nzl = jv    !last wavenumber with non-zero opacity
319      continue

         if(nzf.eq.0) then
          write(6,*)' we never got started finding opac.ne.0'
          do 381 jv=1,nwnos
            sumopt1 = sumopt1+opacos(1,jv)
            sumoptk = sumoptk+opacos(ktemp,jv)
381       continue
          write(6,*)' sum of all opacities at T=1 and T=ktemp:'
          write(6,382)sumopt1,sumoptk
382       format(1p2e13.3)
          write(6,*)' ktemp = ',ktemp
         end if

	 write(6,334) wnos(nzf),wnos(nzl),nzf,nzl
334      format(' opacity is non-zero from ',f9.2,
     &      ' to',f9.2,' cm-1 (=point',i5,' to',i6,')')

	nwnos = nzl - nzf + 1
	WRITE(8,INPUTOSMOL)


	 DO 311 jv = nzf,nzl

	 WRITE(8,330) WNOS(JV),(OPACOS(IT,jv),IT=1,KTEMP)


         do 3111 it=1,ktemp
         if(opacos(it,jv).lt.0.) then
            write(6,3110) jv,it,wnos(jv),1.e4/wnos(jv),opacos(it,jv)
3110        format(i7,i3,f10.3,f8.5,1pe12.3)
         end if
3111     continue

	 if (jv.ge.2 .and.jv.le.nwnos-1) then
c         if (jv.ge.nzf+1 .and.jv.le.nzl-1) then
	       dw = (wnos(jv+1)-wnos(jv-1))/2.
	 else if (jv.eq.1) then
	      dw = wnos(2) - wnos(1)
	 else if (jv.eq.nwnos) then
	      dw = wnos(nwnos) - wnos(nwnos-1)
	 end if
	    do 352 it=1,ktemp
	    sumop_m(it) = sumop_m(it) + opacos(it,jv)*dw
352       continue
c         end if

	 IF (jderiv .EQ. 0) go to 311    ! do not list OS derivatives

	 DO 312 IT=1,KTEMP
	 IF (OPACOS(IT,jv).EQ.0) THEN
	   DO 314 ITN=1,KTEMP
	   DADT(ITN)=0.
314        CONTINUE
	   GO TO 315
	 END IF
312      OPT(IT)=LOG( OPACOS(IT,jv) )

C  For each JV, calculate the derivative, DADT(IT)=d(ln(opt)/d(tmol), 
C  at each temperature.

	 CALL TB04AT(NTEMP,KTEMP,TMOL,OPT,DADT,W)

	 DO 313 IT=1,KTEMP
313      DADT(IT) = DADT(IT)*OPACOS(IT,jv)    !dadt = d(opt)/d(tmol)

315      CONTINUE

	 WRITE(8,331) (DADT(IT),IT=1,KTEMP)

311      CONTINUE    !next OS frequency point


	 WRITE(18,*) NWNOS,' lines succesfully written in os.dat '
	 write(18,*) ' From',wnos(nzf),' to',wnos(nzl),' cm^-1'
	 write(18,155) ((sumop_m(it)/1.d5),it=1,ktemp)
155      format('sumop_m(sum OS-opacities*dw)t=1,kt in km/mol:'
     &      /,1p12e10.3)
	 write(18,1551) ((sumop_l(it)/1.d5),it=1,ktemp)
1551     format('sumop_l (i.e.sum of all line intensities on list)'
     &     ,' t=1,kt in km/mol:'/,1p12e10.3)
	 write(18,156) ((sumop_c(it)/1.d5),it=1,ktemp)
156      format('sumop_c(i.e. sum of lines inside OS-intervals) t=1,kt:'
     &      /,1p12e10.3)
	 if(sumop_m(1).ne.0.0)
     &   write(18,*) ' op(ktemp)/op(1):',sumop_m(ktemp)/sumop_m(1)
	 if(sumop_l(1).ne.0.0)
     &   write(18,*) ' op(ktemp)/op(1):',sumop_l(ktemp)/sumop_l(1)
	 write(18,*) ' we read ',nrlin,' lines from line list'
	 write(18,*) ' we called opac_com ',nr,' times'
	 if (n17*n16 .gt. 0) then
	  write(18,*) n16,' lines contributed to opac'
	  write(18,*) ' on the average each line contributed to'
     *      ,float(n17)/float(n16),' opac frequency points'
	  write(18,*) ' the average prof factor was',a17/float(n17)
	 end if
C
C
C330     FORMAT (F10.3,1P12E11.4)
330     FORMAT (F10.3,1P12E11.3)
331     FORMAT (2X,1P12E12.4)

      WRITE(18,*) ' '
      write(18,*) ' *** successfully termination of program *** ' 
      go to 901
900   CONTINUE
      WRITE(18,*) ' ***'
      write(18,*) ' *** the program was excecuted with problems *** ' 
      WRITE(18,*) ' ***'
901   CONTINUE
C
C
      call etime(time_2)

      WRITE(6,906) time_2-time_0
      WRITE(18,906) time_2-time_0
906   FORMAT (' TOTAL TIME SPENT = ',F10.1,' seconds')
C
      STOP
      END

C*************************************************************
C Subroutine to compute opacities for OS:
C*************************************************************

      SUBROUTINE OPAC_COM(WS,gf,elow,nws)
C     PARAMETER (NTEMP=12,MT1=NTEMP+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'atomparameter.inc'

      common /chw/ hwst(ntemp), hwvw(ntemp)
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      COMMON /COPAC/HALF,sig_nu(ntemp),hckt(ntemp)
     *    ,sqr_pi_m,sig_sqr,sumop_c(ntemp),a17,n16,n17,npr,j
      common/cvoigt/ pgas(mt1),p_el(mt1),p_h(mt1),hwda(ntemp)
     &              ,pemu(mt1),pp(mt1,nspec)
     &              ,profile,atmos_voigt
      COMMON/CISO/vkmssp,wgtmola,rel_iso_c(3),pc3(3)
     & ,reliso(5),kiso,nhalf,jderiv,kiso_c
     & ,linesc3_direct,lineshcn_direct,linesc2h2_direct
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      character profile*5,atmos_voigt*45,atmos*45
      common /comdata/ BDA,pi,c_vel
       character elm*2, nameatom*2
       common /cvaldline/ wgtatom, vv0, vexs, vgf,
     &    gam_rad, gam_stark, gam_waal, fac_lande, 
     &    conv_gstel(ntemp),
     &    jatom, jg, ion, lgstel,
     &    elm, nameatom
C
C Initially NW has been determined by subroutine range such that 
C wnos(nw) < v0 < wnos(nw+1); i.e. the lines can affect the opacities
C in opacity sampling points close to nw. Initially nwp=nwq=nwr=nw,
C but while computing the band we keep nwp,nwr,nwq as the one 
C corresponding to the previous line.

         ncall_opac = ncall_opac + 1
C         if(ncall_opac.lt.10) then
C                 write(6,316) jatom,jg,ion,nameatom
C316              format('ID in op: ',3i4,2x,a2)
C                      write(6,309) (qvib(it),it=1,ktemp)
C                      write(6,315) (pp(it,jg),it=1,ktemp)
C                      write(6,314) (conv_gstel(it),it=1,ktemp)
C309                   format ('Qv in op: ',12f8.2)
C315                   format ('Pp in op: ',1p12e12.3)
C314                   format ('conv_gstel: ',1p12e12.3)
C         end if


C GF=1.87572E-12 Sl[cm/mol] QvQr/exp(-[Gv+Fj]hc/kT)(1-e-nuhc/kT)
C GF = 1.87572E-12*S0
C S0 = 5.33129e+11*gf
      sgf = 5.33129d+11
      absl = sgf * gf

1001  continue
      dw = abs(ws-wnos(nws))
      if ( abs(ws-wnos(max(1,nws-1))) .lt. dw ) then
	   nws = nws-1
	   go to 1001
      else if ( abs(ws-wnos(min(nwnos,nws+1))) .lt. dw ) then
	   nws = nws+1
	   go to 1001
      end if
C   nws is now the nearest OS point to ws. Determine the possible contribution 
C   of line absorption absl to this and surrounding OS points.


      if (profile.eq.'voigt') then

C This is a piece from Bernhard's spectrum program, where he has shown that a very
C broad wing is necessary for the strong lines.
C qntau is the partition function for this atom/ion at maximum temperature.
      qntau = 1.
      qntau = qvib(ktemp)

      testvoi = gf * EXP(-elow/10000.0)/qntau !Testvariable, ob Voigt-Profil breit sein ka
      IF (testvoi .GT. 5.0E-2) THEN
         half = float(nhalf) * 300.0
      ELSE IF (testvoi .GT. 5.0E-3) THEN
         half = float(nhalf) * 150.0
      ELSE IF (testvoi .GT. 5.0E-4) THEN
         half = float(nhalf) * 50.0
      ELSE IF (testvoi .GT. 5.0E-5) THEN
         half = float(nhalf) * 30.0
      ELSE IF (testvoi .GT. 5.0E-6) THEN
         half = float(nhalf) * 10.0
      ELSE
         half =  float(nhalf)
      ENDIF

C      write(6,*) ' half in opac_com = ',half
 
C Damping Constants
C
         IF (gam_rad .EQ. 0.0) THEN
            hwlo = 5.9E-13 * ws * ws !Breite des Lorentz-Profils in cm-1, klassisch
         ELSE
            hwlo = (10.0**gam_rad)/(c_vel*pi*4.0) !aus VALD
         ENDIF

         IF (gam_stark .NE. 0.0) THEN
           gam_stark = gam_stark + 15.193259 !-log(k)-1/6log(10000)
         ENDIF
         IF (gam_waal .NE. 0.0) THEN
           gam_waal = gam_waal + 14.659925 !-log(k)-0.3log(10000)
         ENDIF
         IF (elm .EQ. 'H ') THEN ! H --> hier keine p-Effekte
            gam_stark = 0.0
            gam_waal = 0.0
         ENDIF

      do 1201 it = 1,ktemp
            hwst(it) = gam_stark+LOG10(p_el(it))-(5.0*LOG10(t(it))/6.0)
            hwst(it) = (10.0**hwst(it))/(c_vel*pi*4.0)
            hwvw(it) = gam_waal+LOG10(p_h(it))-(0.7*LOG10(t(it)))
            hwvw(it) = (10.0**hwvw(it))/(c_vel*pi*4.0)
            hwda(it) = hwlo + hwvw(it) + hwst(it)
1201  continue 
  

      end if         ! end if profile.eq.'voigt'



      do 1101 kt = 1,ktemp
      nws1 = nws
      nws2 = nws
      dw = abs(ws-wnos(nws))
      if (dw .gt. (half * sig_nu(kt) * ws) ) go to 1101  !if dist to nearest OS
      if (kt.eq.ktemp) n16 = n16 + 1                   !#lines used in OS at ktemp
	 sig_sqr = sig_nu(kt) * ws
	 ST = absl
     &     *EXP(-elow*HCKT(kt))*(1.-EXP(-ws*HCKT(kt)))/qvib(kt)
	 sumop_c(kt) = sumop_c(kt) + st            !integrated abs.coef.(T) in [cm/mol]
1110     continue
	 if(profile.eq.'voigt') then
	      call voigt_entry(kt,prof,ws,dw)
	 else
	 prof = (sqr_pi_m/sig_sqr)*exp(-((dw/sig_sqr)**2))
	 end if
         opadd = st * prof                             !abs.coef. in cm2/mol at nws1 for T=kt
         if (lgstel.eq.1) opadd = opadd * conv_gstel(kt)   !abs.coef. in cm2/g* 
	 OPACOS(kt,nws1) = OPACOS(kt,nws1) + opadd
	 if (kt.eq.ktemp) n17 = n17 + 1 
	 if (kt.eq.ktemp) a17 = a17 + prof
	 nws1 = nws1-1
	 if(nws1.lt.1) go to 1120
	 dw = abs(ws-wnos(nws1))
	 if (dw .le. (half * sig_nu(kt) * ws) ) go to 1110
1120     continue
	 nws2 = nws2 + 1
	 if (nws2.gt.nwnos) go to 1101
	 dw = abs(ws-wnos(nws2))
	 if (dw .gt. (half * sig_nu(kt) * ws) ) go to 1101
	 if(profile.eq.'voigt') then
              call voigt_entry(kt,prof,ws,dw)
         else
         prof = (sqr_pi_m/sig_sqr)*exp(-((dw/sig_sqr)**2))
         end if
         opadd = st * prof                             !abs.coef. in cm2/mol at nws1 for T=kt
         if (lgstel.eq.1) opadd = opadd * conv_gstel(kt)   !abs.coef. in cm2/g* 
         OPACOS(kt,nws2) = OPACOS(kt,nws2) + opadd
         if (kt.eq.ktemp) n17 = n17 + 1 
         if (kt.eq.ktemp) a17 = a17 + prof
         go to 1120
1101  CONTINUE

      return
      end


C*************************************************************
C Subroutine to compute the wavenumbers for OS:
C*************************************************************

      SUBROUTINE OSWAVENUMBERS
C      PARAMETER (NWL=149000,NTEMP=12)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'atomparameter.inc'
      COMMON /COSLIST/ WNB(25),WNSTEP(25),WNEND,INTVOS
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      COMMON/CMOLSP/MOLID,OSFIL
      CHARACTER MOLID*4,OSFIL*60,inwnfil*60
      NAMELIST /INPUTOS/ WNB,WNSTEP,WNEND,INTVOS,OSFIL,listwn,inwnfil
     *   ,losresl,osresl,wnos_first,wnos_last
C
      READ (12,INPUTOS)
      WRITE(18,INPUTOS)
C
C  Calculate (or read, if listwn=1) the wave numbers, WN, for OS opacity 
C  and the total number, NWNOS, of OS frequency points for the OS-tables.
C  If losresl = 1, the os wavenumbers are computed based on a specified 
C  spectral resolution, osresl. Else it is computed in fixed steps inside a
C  number of prespecified intervals.
C

      IF (LISTWN .eq. 1) go to 225   ! read an existing OS - wn list


C  compute an OS - wn list:

      if (losresl.eq.0) go to 228     ! compute list with fixed steps

C  compute OS list with fixed resolution, osresl, through spectrum:
      
      wnos(1) = wnos_first
      step = (wnos_last-wnos_first)/osresl
      do 240 k = 1,nwl-1
      nwnos = k + 1
      wnos(nwnos) = wnos(k) + step
      if (wnos(nwnos) .gt. wnos_last) go to 241
240   continue
      write(6,*) ' wnos_first, wnos_last, nwnos, wnos(nwnos) = '
     &  ,wnos_first, wnos_last, nwnos, wnos(nwnos)
      stop ' error: increase dimension for wnos '
241   continue

      write(6,*) ' We did a resolution based set of os wavenumbers'
      write(6,245) nwnos,osresl,step,wnos(1),wnos(nwnos)
     & ,1.e4/wnos(nwnos),1.e4/wnos(1)
245   format(' There are ',i6,' OS wavenumbers.'
     & ,/' Resolution is',f8.0,' (corresponding to step factor',f8.5,')'
     & ,/' OS interval:',f7.1,'-',f9.1,' cm^-1 (=',f5.2,'-',f5.1,'mu)')
      

      go to 226

228   continue
      L=0
      WNOS(1)=WNB(1)
      WNB(INTVOS+1) = WNEND
      DO 100 I = 1,INTVOS
      WNLAST = WNB(I+1)-WNSTEP(I)
102   L=L+1
      L1=L+1
      WNOS(L1)=WNOS(L)+WNSTEP(I)
      IF ( WNOS(L1).LE.WNLAST ) GO TO 102
100   CONTINUE
      NWNOS = L1

      go to 226

225   continue
C     read an existing OS - wn list:

      open (unit=14,file=inwnfil,status='old')
      do 200 l=1,1000000
      read(14,*,end=201) wnos(l),idum
      nwx = l
200   continue
201   continue
      NWNOS = nwx
      im = 0
      do 210 i=1,nwx
      if (wnos(i).lt.wnb(1)) then
       im = i 
       go to 210
      end if
        iu = i-im
        wnos(iu) = wnos(i)
        nwnos = iu
      if (wnos(iu).gt.wnend) then
       nwnos = iu-1
       go to 211
      end if
210   continue
211   continue


226   continue

      DO 103 I=1,NWNOS
      L=NWNOS-I+1
103   WLOS(L)=1.D8/WNOS(I)

      DO 104 L=2,NWNOS-1
      LP=L+1
      LM=L-1
104   WLOSSTEP(L) = ( WLOS(LP)-WLOS(LM) ) / 2.
      WLOSSTEP(1) = WLOS(2)-WLOS(1)
      WLOSSTEP(NWNOS) = WLOS(NWNOS)-WLOS(NWNOS-1)


      IF (NWNOS.GT.NWL)  STOP ' DIMENSION NWL TOO SMALL'
      WRITE(18,12) OSFIL,NWNOS,WNOS(1),WNOS(NWNOS),
     #             KTEMP,(TMOL(K),k=1,ktemp)
12    FORMAT (' End of call to subroutine OSWAVENUMBERS:',/
     &    ' The OS file will be written on file ',A60,/
     &    ' There will be ',I7,' frequency points from wnos(1)=',
     &    f9.2,' to wnos(nwnos)=',f9.2,' in the sampling file',/
     &    ' and ',I3,' temperatures, which are: ',12f7.0,' K')


      RETURN
      END

C *******************************************************************
      subroutine mol_specific

*  this subroutine computes or reads the molecular specific data
*  such as the partition function, Qv for the temperatures 
*  tmol(k),k=1,ktemp, the isotopes, molecular weights, etc.
C *******************************************************************

C      PARAMETER (NTEMP=12,MT1=NTEMP+1,NQT=200)
      PARAMETER (NQT=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'atomparameter.inc'
      character cu*2,symu*1,cl*2,syml*1,br*1
      DIMENSION ncoi(9),ncni(3),nc2i(3),p_iso(3)
     & ,wgt_iso_h(2),rel_iso_h(2)
     & ,wgt_iso_c(3)
     & ,wgt_iso_n(3),rel_iso_n(3)
     & ,wgt_iso_o(3),rel_iso_o(3)
     & ,wgt_iso_t(5),rel_iso_t(5),ftio(7)
     & ,wgt_iso_s(4),rel_iso_s(4)
     & ,wgt_iso_si(3),rel_iso_si(3)
     & ,wgt_iso_fe(4),rel_iso_fe(4)
     & ,wgt_iso_cr(4),rel_iso_cr(4)
     & ,wgt_iso_mg(4),rel_iso_mg(4)
     & ,tin(nqt),qvin(nqt),qtln(nqt),dqdt(nqt),w(3,nqt)
     & ,wsiso(10),linesgfz(10),nc3(3)
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      COMMON/CMOLSP/MOLID,OSFIL
      COMMON/CISO/vkmssp,wgtmola,rel_iso_c(3),pc3(3)
     & ,reliso(5),kiso,nhalf,jderiv,kiso_c
     & ,linesc3_direct,lineshcn_direct,linesc2h2_direct
      character MOLID*4,osfil*60,aminus*1,adum*66,ftpsite*26
       character elm_pick*3
       common /cvald/ wb, we, lvald, ion_pick, elm_pick
      common/cvoigt/ pgas(mt1),p_el(mt1),p_h(mt1),hwda(ntemp)
     &              ,pemu(mt1),pp(mt1,nspec)
     &              ,profile,atmos_voigt
      character profile*5,atmos_voigt*45
      NAMELIST /SAMPLING/jderiv
     & ,linesc2h2_direct,lhitr_c2h2
     & ,lineshcn_direct,lhitr_hcn
     & ,linesc3_direct
     & ,wgt_iso_h,rel_iso_h,kiso_h
     & ,wgt_iso_c,rel_iso_c,kiso_c,lgvc2,lqc2_red
     & ,wgt_iso_n,rel_iso_n,kiso_n
     & ,wgt_iso_o,rel_iso_o,kiso_o
     & ,wgt_iso_t,rel_iso_t,kiso_t,ftio
     & ,wgt_iso_s,rel_iso_s,kiso_s
     & ,wgt_iso_si,rel_iso_si,kiso_si
     & ,wgt_iso_fe,rel_iso_fe,kiso_fe
     & ,wgt_iso_cr,rel_iso_cr,kiso_cr
     & ,wgt_iso_mg,rel_iso_mg,kiso_mg
C
       ana = 6.022136d23                  !Advogadro's number [molecules/mol]

      READ(12,sampling)
      WRITE(18,sampling)
      if(lvald.eq.1) go to 998

       open (unit=66,status='unknown',
     &   file='/ast/p2/uffegj/lines.tmp')
c      open (unit=66,status='new',
c    &   file='/ast/p2/uffegj/lines.tmp')

      IF (MOLID.NE.'CO  ') go to 901
C----------------------------------------------------------------------
C                                         CO     CO      CO     CO    |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       do 117 kc=1,kiso_c
       do 117 ko=1,kiso_o
          relabu = rel_iso_c(kc)*rel_iso_o(ko)
          wgtmola=(wgt_iso_c(kc)+wgt_iso_o(ko))*relabu + wgtmola
117    continue
       if (kiso_o .eq. 1) then
         do 1171 kc=1,kiso_c
1171     reliso(kc) = rel_iso_c(kc)
       end if
       kiso = kiso_c * kiso_o
       write(18,*) 'mean molecular weight for CO is',wgtmola

       open (unit=65,status='old',readonly,
     &   file='/ast/p2/uffegj/moldat/co/linelist.goorvitch')
       open (unit=67,status='old',readonly,
     &   file='/ast/p2/uffegj/moldat/co/partfunc.goorvitch')

       read(67,110) adum
110    format(a66)
C T[K] 12-16 12-17 12-18 13-16 13-17 13-18 14-16
       kqt = 0
       do 111 k=1,1000
       read(67,*,end=119) tin(k),qvin(k),q1,q2,q3,q4,q5,q6
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
111    continue
119    continue
       write(6,*) ' I read',kqt,' partition function values'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)


       nco=0
       do 118 k=1,9
118    ncoi(k)=0
       write(6,*) ' First 5 lines of line strength data:'
       do 112 i=1,100000000
        read(65,114,end=113) 
     &      v0,exs,gf,ivl,ivu,jl,ju,isoc,aminus,isoo
        if (i.le.5) write(6,115)
     &      v0,exs,gf,ivl,ivu,jl,ju,isoc,aminus,isoo
        idx=(isoc-11)+(isoo-16)*3
        ncoi(idx) = ncoi(idx) + 1
        if (isoc.eq.14) go to 112
        nco = nco + 1
        gfiso = gf * rel_iso_c(isoc-11) * rel_iso_o(isoo-15)
        write(66,116)v0,exs,gfiso
112    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
       go to 999
113    continue
       write(6,*) ' Last line of line strength data:'
       write(6,115)v0,exs,gf,ivl,ivu,jl,ju,isoc,aminus,isoo
       write(6,*)' wrote #lines for 12C16O,13C16O,14C16O,..17O,18O:'
       write(6,109) (ncoi(k),k=1,9),nco
114    format(f11.4,f12.4,e13.5,4i4,i3,a1,i2)
C    3.6759      0.0000  2.10800e-08   0   0   0   1 13-16
115    format(2f11.4,1pe13.5,0p4i3,i3,a1,i2)
116    format(2f11.4,1pe13.5)
109    format(10i8)

C                                                                     ^
C                                         CO     CO      CO     CO    |
C----------------------------------------------------------------------
      GO TO 990
901   CONTINUE

      IF (MOLID.NE.'CN  ' .and. 
     &    MOLID.NE.'CNqc' .and. MOLID.NE.'CNvq') go to 902
C----------------------------------------------------------------------
C                                         CN     CN      CN     CN    |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       do 127 kc=1,kiso_c
       do 127 kn=1,kiso_n
          relabu = rel_iso_c(kc)*rel_iso_n(kn)
          wgtmola=(wgt_iso_c(kc)+wgt_iso_n(kn))*relabu + wgtmola
127    continue
       if (kiso_n .eq. 1) then
         do 1271 kc=1,kiso_c
1271     reliso(kc) = rel_iso_c(kc)
       end if
       kiso = kiso_c * kiso_n
       write(18,*) 'mean molecular weight for CN is',wgtmola
       open (unit=67,status='old',readonly,
     &   file='/ast/p2/uffegj/cn/scan_cn.partfunc')

C
       read(67,110) adum
       kqt = 0
       do 121 k=1,1000
       read(67,*,end=129) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
121    continue
129    continue
       write(6,*) ' I read',kqt,' partition function values for CN'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)


          k_violet = 0
1272   continue

       IF(molid.eq.'CN  ' .or. molid.eq.'CNvq') THEN
       open (unit=65,status='old',readonly,
     &          file='/ast/p2/uffegj/cn/scan_cn.dat')
       ELSE IF(molid.eq.'CNqc') THEN
       open (unit=65,status='old',readonly,
     &          file='/ast/p2/uffegj/cn/qc_cn.dat')
       END IF

       ison = 1
       rel_iso_n(ison) = 1.D0    !971003 only 14N on scan_cn
       ncn=0
       ngfz=0
       do 128 k=1,3
128    ncni(k)=0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0,gf,exs,ivu,ivl,jl,ibr,isoc '
       do 122 i=1,100000000
          IF(molid.eq.'CN  ' .or. molid.eq.'CNvq') THEN
        read(65,*,end=123) v0,gf,exs,ivu,ivl,jl,ibr,isoc
        ncni(isoc) = ncni(isoc) + 1
          ELSE IF(molid.eq.'CNqc') THEN
        read(65,*,end=123) v0,gf,exs,isoc
         if(k_violet.eq.1) then
           if(v0.lt.21100.0) go to 122
         end if
        ncni(isoc) = ncni(isoc) + 1
          ison=1     !..^14N
          if(isoc.eq.3) then
            ison=2   !^12C^15N
            isoc=1
          end if
          END IF
        if (i.le.5) write(6,125) v0,gf,exs,ivu,ivl,jl,ibr,isoc
        ncn = ncn + 1
        gfiso = gf * rel_iso_c(isoc) * rel_iso_n(ison)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 122
        end if
        write(66,116)v0,exs,gfiso
122    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
c      go to 999
123    continue
       write(6,*) ' Last line of line strength data:'
       write(6,125) v0,gf,exs,ivu,ivl,jl,ibr,isoc
       write(6,*)' #lines read: 12CN, 13CN, 12C15N, total; wrote'
       write(6,109) (ncni(k),k=1,3),ncn,ncn-ngfz
125    format(f11.3,1pe13.6,0pf10.2,2i3,i5,i3,i2)
C   333.007 2.045595E-05  24631.21  7 11  -48  5 1

C if we want to add Querci's violet CN system to the Scan_CN red system:
       IF(molid.eq.'CNvq') THEN
          close(65)
          molid = 'CNqc'
          k_violet = 1
          go to 1272
        END IF


C                                                                     ^
C                                         CN     CN      CN     CN    |
C----------------------------------------------------------------------
      GO TO 990
902   CONTINUE


      IF (MOLID.NE.'CH  ') go to 903
C----------------------------------------------------------------------
C                                         CH     CH      CH     CH    |
C                                                                     V
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       do 137 kc=1,kiso_c
       do 137 kh=1,kiso_h
          relabu = rel_iso_c(kc)*rel_iso_h(kh)
          wgtmola=(wgt_iso_c(kc)+wgt_iso_h(kh))*relabu + wgtmola
137    continue
       if (kiso_h .eq. 1) then
         do 1371 kc=1,kiso_c
1371     reliso(kc) = rel_iso_c(kc)
       end if
       kiso = kiso_c * kiso_h
       write(18,*) 'mean molecular weight for CH is',wgtmola

       open (unit=65,status='old',readonly,
     &   file='/ste1/ftp/pub/scan/scan_ch.dat')
       open (unit=67,status='old',readonly,
     &   file='/ast/p2/uffegj/ch/scan_ch.partfunc')

C
       read(67,110) adum
       kqt = 0
       do 131 k=1,1000
       read(67,*,end=139) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
131    continue
139    continue
       write(6,*) ' I read',kqt,' partition function values for CH'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)

       nch=0
       ngfz=0
       isoh = 1
       rel_iso_h(isoh) = 1.d0
       write(6,*) ' Only one H-isotope on Scan-CH line list,'
       write(6,*) ' so we set isoh==1, rel_iso_h(1)==1.0,'
       write(6,*) ' and gfiso = gf * rel_iso_c(isoc)'
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13'
       write(6,*) ' v0,gf,exs,ivu,ivl,jl,ibr,isoc '
       do 132 i=1,100000000
        read(65,*,end=133) iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13
C 1.line: 2  0 11  13 20 4.67100E-11   333.33 24723.21  37.18
        if (i.le.5) write(6,135) iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13
        do 1321 isoc = 1,kiso_c
        nch = nch + 1
        gfiso = gf * rel_iso_c(isoc) * rel_iso_h(isoh)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 1321
        end if
        v0 = v0 + (isoc-1)*dw13
        write(66,116)v0,exs,gfiso
1321   continue
132    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
c      go to 999
133    continue
       write(6,*) ' Last line of line strength data:'
       write(6,135) iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13
       write(6,*)' read 12CH+13CH (total); wrote total, #lines:'
       write(6,109) nch,nch-ngfz
135    format(3i3,i4,i3,1pe12.5,0p2f9.2,f7.2)
C 2  0 11  13 20 4.67100E-11   333.33 24723.21  37.18
C                                                                     ^
C                                         CH     CH      CH     CH    |
C----------------------------------------------------------------------
      GO TO 990
903   CONTINUE


      IF (MOLID.NE.'C2  ' .and. MOLID.NE.'C2H ') go to 904
C----------------------------------------------------------------------
C                                         C2     C2      C2     C2    |
C                                                                     V
C As a test to compute the C2H absorption coefficient we use the C2
C linelist data in the region 2545 - 7360 cm-1 (1.36 - 3.93 microns)
C where the position of the bands are almost identical and the integrated 
C band strengths comparable, when the scaled Querci data are compared with
C the low-resolution plot of Langhoffs data in Goorvitch et al 1983.
C The latter data have been transformed to cm^2?mol and stored in 
C /ste1/uffegj/c2h/gorvlangh.dat. The strongest C2H band is almost identical
C to the corresponding C2 band (in the 2.5-3 micron region), while the two next
C C2 bands are scaled a bit (a factor 20 and 3.8, resp) in order to look more like
C Langhoff's C2H data.
       if (MOLID.eq.'C2H ') then
       lgvc2 = 0
       lqc2_red = 1
       end if

C    Calculate mean molecular weight
          p_iso(1) = rel_iso_c(1)**2    ! 12C12C
          p_iso(3) = rel_iso_c(2)**2    ! 13C13C
          p_iso(2) = 2.D0*rel_iso_c(1)*rel_iso_c(2)  !12C13C
          wgtmola = 2.D0*wgt_iso_c(1)*p_iso(1) +
     &              2.D0*wgt_iso_c(2)*p_iso(3) +
     &    (wgt_iso_c(1)+wgt_iso_c(2))*p_iso(2)
         do 1471 kc=1,kiso_c
1471     reliso(kc) = rel_iso_c(kc)
       kiso = kiso_c
       write(18,*) 'mean molecular weight for C2 is',wgtmola
       write(6,692) wgtmola,p_iso(1),p_iso(2),p_iso(3)
692    format(' C2 mean mol.wgt:',f7.2,' P{12C2,13C2,12C13C}:',
     &       3f8.4)
       if (lgvc2.eq.1)
     & write(6,693) p_iso(1)/2.,p_iso(2)/8.,p_iso(3)/8.
693    format(' C2 gf isotopic mult for 12C2,13C2,12C13C:',
     &       1p3e12.3)

       if (lgvc2.ge.1) then
       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/c2/gvlist.dat')
       else
       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/querci/c2_qc.dat')
       end if

       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/querci/c2_qc.partfunc')

C
       do 140 i=1,5
140    read(67,110) adum
       kqt = 0
       do 141 k=1,1000
       read(67,*,end=149) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
141    continue
149    continue
       write(6,*) ' I read',kqt,' partition function values for C2'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)
       if (lgvc2.eq.1) write(18,*)
     *    ' We used gf C2 data from Goorvitch (OBS: only 2.5mu band)'
       if (lgvc2.eq.2) write(18,*)
     *    ' We used gf C2 data from Goorvitch and multiplied his gf',
     *    ' with 1, 1/2, 1/4 for 12C2, 12C13C, 13C2 as recommended ',
     *    ' in G90 when gf is used with Q of S&T84 or Tatum66'
       if (lgvc2.eq.0 .and. lqc2_red.eq.0) write(18,*)
     *    ' We used gf C2 data from Querci'
       if (lgvc2.eq.0 .and. lqc2_red.eq.1) write(18,*)
     *    ' We used gf C2 data from Querci, reduced with a factor 10',
     *    ' for lam>1.5mu (v<6667cm^-1); f lin incr from 1.15-1.5mu'

       nc2=0
       ngfz=0
       do 148 k=1,3
148    nc2i(k)=0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0,exs,gf,isoc '

       do 142 i=1,100000000
        read(65,*,end=143) v0,gf,exs,isoc
        if (i.le.5) write(6,145) v0,exs,gf,isoc
C       if (lgvc2.eq.1 .and. v0.le.3934. .and. v0.ge.3930. ) 
C    &      write(6,145) v0,exs,gf,isoc
        nc2i(isoc) = nc2i(isoc) + 1
        nc2 = nc2 + 1

C in order to follow the examples given in Table 5 of Goorvitch 1990, which
C however, in my (UGJ) oppinion is inconsistent with the definition of the 
C partition functions in Sauval&Tatum84 or Tatum66 used here (and in G90).
       if (lgvc2.eq.2) then
        if(isoc.eq.1) gf = gf/1.d0
        if(isoc.eq.2) gf = gf/2.d0
        if(isoc.eq.3) gf = gf/4.d0
       end if

        if (lgvc2.eq.0 .and. lqc2_red.eq.1) then
          if (1.d4/v0.ge.1.5d0) then
              fred = 10.d0
          else if (1.d4/v0.ge.1.15d0) then
              fred = 1.d0 + 9.d0*(1.e4/v0-1.15d0)/0.35d0
          else 
              fred = 1.d0
          end if
        gf = gf / fred
        end if 

       IF (MOLID.eq.'C2H ') then
         if(v0.lt.2545.0 .or. v0.gt.7360.0 ) go to 142
         if(v0.gt.4170.0) then
           if(v0.lt.5830) gf = gf * 0.0513
           if(v0.ge.5830) gf = gf * 0.263
         end if
       END IF


        gfiso = gf * p_iso(isoc)

        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 142
        end if
        write(66,116)v0,exs,gfiso
142    continue

       write(6,*)' Stop: more than 1.e8 lines in input file ??'
       go to 999
143    continue
       write(6,*) ' Last line of line strength data:'
       write(6,145) v0,exs,gf,isoc
       write(6,*)' read 12C2, 12C13C, 13C2, total; ',
     &        'wrote total, #lines:'
       write(6,109) (nc2i(k),k=1,3),nc2,nc2-ngfz
145    format(2f11.4,1pe13.5,i4)
C                                                                     ^
C                                         C2     C2      C2     C2    |
C----------------------------------------------------------------------
      GO TO 990
904   CONTINUE


      IF (MOLID.NE.'C3  ') go to 905
C----------------------------------------------------------------------
C                                         C3     C3      C3     C3    |
C                                                                     V

C There are computed 2 isotopes of C3: 12C3, 12C213C, corresponding
C to iso = 1, 2, 3. Iso = 3 are the 12C-12C-13C molecules, iso = 2
C are the 12C-13C-12C molecules, iso = 1 are the 12C3 molecules.
C The strewngth of the 13C3 are include into the gf values of the 
C 12C-13C-12C molecules (iso=2) for symmetry reasons.
C    Calculate mean molecular weight
       write(6,*) ' now ready to do C3 mol-specific'
       pc3(1) = (rel_iso_c(1))**3
       pc3(2) = (1.-pc3(1)) / 3.0  !there are 3 different 13C12C2 symmetries
       pc3(3) = pc3(2)
       wm121212 = 3.*wgt_iso_c(1) * pc3(1)
       wm131212 = (2.*wgt_iso_c(1)+wgt_iso_c(2)) * (1.-pc3(1))
       wgtmola = wm121212 + wm131212
       write(18,*) 'mean molecular weight for C3 is',wgtmola
       if (kiso_c.eq.1) then
         write(18,*) ' we consider only 12C3 lines'
       else
         write(18,*) ' we consider 12C3, 2*12C-12C-13C and 12C-13C-12C'
         write(18,*)' 13C3 and 13C212C lines are put into 12C213C lines'
       end if 
         do 1571 kc=1,kiso_c
1571     reliso(kc) = rel_iso_c(kc)
       kiso = kiso_c
       c1312=1./rel_iso_c(1) -1.
       if (c1312.ne.0) then
       write(18,1572) kiso,(reliso(kc),kc=1,kiso),1./c1312
       else
       write(18,1573) kiso,(reliso(kc),kc=1,kiso)
       end if
1572   format(' #isotopes=',i2,' reliso=',2f7.4,' 12C/13C=',f7.2)
1573   format(' #isotopes=',i2,' reliso=',2f7.4,' Only 12C considered')

       open (unit=65,status='old',readonly,
     &   file='/ast/p2/uffegj/2013/c3/c3linelist.dat')
c    &   file='/ste2/uffegj/c3_linelist.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/c3/c3.partfunc')

C
       read(67,110) adum
       kqt = 0
       do 151 k=1,1000
       read(67,*,end=159) tin(k),q1,q2,qvin(k),q4,q5,q6,q7,q8
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
       B2 = 0.4305
       qvin(k) = (0.69493*tin(k)/B2+0.33333) * qvin(k)
151    continue
159    continue
       write(18,*)' I read',kqt
     &    ,' partition function values for C3: Q_vr, Q_v'
       write(18,*) ' from t=',tin(1),' to t=',tin(kqt)
       do 1591 kt=1,kqt
1591   write(6,'(f8.1,f12.2,f10.2)') tin(kt),qvin(kt)
     &  ,qvin(kt)/(0.69493*tin(kt)/B2+0.33333)

C actually the vibrational partition functions should be multiplied
C with the rotational partition function at each temperature and for
C each bands, as is done in the C3.f program. Here we assume that 
C the B(v) molecular constant can be approximated by B(000) so that
C the sequence:
C      do 1200 kt = 1, ktemp
C      IF(L2.EQ.0) THEN
C        Q_vr(kt) = (0.69493*tmol(kt)/B2+0.33333) * Qvib(kt)
C      ELSE
C        Q_vr(kt) = ( 0.69493*tmol(kt) / B2
C     &    *EXP(-B2*L2*(L2+1)*1.439/tmol(kt))
C     &   +0.33333+L2*(1.-L2*(4.4444E-3+L2/3000.)) ) * Qvib(kt)
C      END IF
C1200  continue
C in the C3.f program can be approximated by Qr = 0.69493*tmol(kt)/B2+0.33333
C B2 = B(V_low) of lower vibrational level in the band transition, is computed as
C      BS=0.
C      DO 132 I=1,3
C132    BS=BS+A(I)*V(I)       --- V(i) -> V(i)+0.5 ?????
C      B(V1,V2,V3)=BE+BS
C      B2 = B(v1l,v2l,v3l)
C where V(i) as well as v1,v2,v3 are the lower vibrational quantum numbers.
C BE = 0.4305;  A{1,2,3} = 0.019, 0.0122, -0.01,   
C so that B(0,0,0) = Be = 0.4305

       if (linesc3_direct.eq.1) then
        write(6,*) ' we read directly from the line list'
        close(unit=65)
        close(unit=66)
        open (unit=66,status='old',readonly,
     &   file='/ste1/uffegj/c3/c3_linelist.dat')
c    &   file='/ste2/uffegj/c3_linelist.dat')
        go to 990
       end if

       do 1521 i=1,3
1521   nc3(i)=0
       ngfz=0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*)
     &    'v0,jl,vl,vm,vn,ll,ju,vi,vj,vk,lu,gf,exs,b2,iso'
       do 152 i=1,100000000
        read(65,*,end=153)
     &     v0,jl,ivl,ivm,ivn,ll,ju,ivi,ivj,ivk,lu,gf,exs,b2,iso
C  746.344 299 0  0 0  0 298 0  1 0  1 6.386E-09 35521.  0.431 3
       nc3(iso) = nc3(iso) + 1
        if (i.le.5) write(6,155)
     &     v0,jl,ivl,ivm,ivn,ll,ju,ivi,ivj,ivk,lu,gf,exs,b2,iso
        gfiso = gf * pc3(iso)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 152
        end if

        write(66,156)v0,exs,gfiso
152    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
c      go to 999
153    continue
       write(6,*) ' Last line of line strength data:'
       write(6,155)
     &     v0,jl,ivl,ivm,ivn,ll,ju,ivi,ivj,ivk,lu,gf,exs,b2,iso
       write(6,*)' I read 12C3, 12C-13C-12C, 12C-12C-13C; wrote total,'
     &   ,' #lines:'
       nnc3 = nc3(1) + nc3(2) + nc3(3)
       write(6,109) nc3(1),nc3(2),nc3(3),nnc3-ngfz
155    format(f10.3,2(1x,i3,i2,i3,i2,i3),1pe10.3,0pf7.0,f7.3,i2)
156    format(2f9.2,1pe11.4)  !larger uncertainty and need for more space

C                                                                     ^
C                                         C3     C3      C3     C3    |
C----------------------------------------------------------------------
      GO TO 990
905   CONTINUE


      IF (MOLID.NE.'C2H2') go to 906
C----------------------------------------------------------------------
C                                         C2H2   C2H2    C2H2   C2H2  |
C                                                                     V
C Detailed C2H2 opacity sampling is best produced directly from the 
C program /ste1/uffegj/c2h2/c2h2.f. This program can also produce a
C line list, and an OS can be produced from that.
C The old C2H2 program produced 13C lines by multiplying the wavenumber
C of the 12C lines by 1.01 and doing nothing else. However, this is
C not correct as concerns the intensity (and of course only a rough 
C approximation concerning the line frequency), so her it is completley
C abandoned for the moment (980701).

C    Calculate mean molecular weight
       wm1212 = 2.*wgt_iso_c(1) * (rel_iso_c(1))**2
       wm1313 = 2.*wgt_iso_c(2) * (rel_iso_c(2))**2
       wm1213 = (wgt_iso_c(1)+wgt_iso_c(2)) * rel_iso_c(1)*rel_iso_c(2)
       wgtmolc = wm1212 + 2.*wm1213 + wm1313
       wm11 = 2.*wgt_iso_h(1) * (rel_iso_h(1))**2
       wm22 = 2.*wgt_iso_h(2) * (rel_iso_h(2))**2
       wm12 = (wgt_iso_h(1)+wgt_iso_h(2)) * rel_iso_h(1)*rel_iso_h(2)
       wgtmolh = wm11 + 2.*wm12 + wm22
       wgtmola = wgtmolc + wgtmolh
       write(6,*) 'mean molecular weight for C2H2 is',wgtmola
       if (kiso_c.eq.1)  then
         write(18,*) ' only 12C2H2 lines are considered'
       else
         write(18,*)' 12C13CH2 lines are considered with v0 = '
         write(18,*)' 1.01*v0(12C2H2) and gf(12C13CH2) = gf(12C2H2)'
       end if
         do 1671 kc=1,kiso_c
1671     reliso(kc) = rel_iso_c(kc)
       kiso = kiso_c

       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/c2h2/linelist.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/c2h2/c2h2.partfunc')

C
       B2 = 1.1765
       temp = 296.
       qh296 = (1.-(EXP(-3294.85*1.4388/temp)))**(-1)
     &  * (1.-(EXP(-1974.32*1.4388/temp)))**(-1)
     &  * (1.-(EXP(-3372.85*1.4388/temp)))**(-1)
     &  * (1.-(EXP(-612.87*1.4388/temp)))**(-2)
     &  * (1.-(EXP(-730.33*1.4388/temp)))**(-2)
       write(6,*) ' Harmonic Q_v(T=296) = ',qh296
       qh296 = 2.*(0.69493*296./B2+0.33333) * qh296
       write(6,*) ' Total Harmonic Q_vr(T=296) = ',qh296

       read(67,110) adum
       kqt = 0
       do 161 k=1,1000
       read(67,*,end=169) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
       qvin(k) = 2.*(0.69493*tin(k)/B2+0.33333) * qvin(k)
161    continue
169    continue
       write(6,*) ' I read',kqt
     &     ,' vibrational partition function values for C2H2'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)
       write(6,*) ' The product of vibrational and rotational '
     &   ,' Q, Q_vr(T), is computed as:'
       do 1691 kt=1,kqt
1691   write(6,'(f8.1,f12.2)') tin(kt),qvin(kt)
C 
C the vibrational partition function is calculated in the program C2H2.f
C as
C      do 200 kt=1,ktemp
C      Qvib(kt) = 0.
C      temp = tmol(kt)
C
C      DO 134 V1=0,M1
C      B1=EXP(-V1*3294.85*1.4388/temp)
C      DO 134 V2=0,M2
C      B2=EXP(-V2*1974.32*1.4388/temp)
C      DO 134 V3=0,M3
C      B3=EXP(-V3*3372.85*1.4388/temp)
C      DO 134 V4=0,M4
C      B4=EXP(-V4*612.87*1.4388/temp)
C      DO 134 V5=0,M5
C      B5=EXP(-V5*730.33*1.4388/temp)
C      DO 134 L4=V4,0,-2
C      DO 134 L5=V5,0,-2
C      FL1=1.
C      L21=IABS(L4-L5)
C      L22=IABS(L4+L5)
C      IF(L21.NE.L22) FL1=2.
C      IF(V4+V5.NE.0) FL1=FL1*2.
C      Qvib(kt)=FL1*B1*B2*B3*B4*B5+Qvib(kt)
C134   CONTINUE
C      WRITE(18,*) ' T,Q_vib(T) = ',tmol(kt),Qvib(kt)
C
C200   Continue
C
C 
C the rotational partition function is calculated in the program C2H2.f
C as
C      do 1200 kt=1,ktemp
C      IF(L2.EQ.0) THEN
C       Q_vr(kt)=2.*(0.69493*tmol(kt)/B2+0.33333) * Qvib(kt)
C      ELSE
C       Q_vr(kt)=( 0.69493*tmol(kt)/B2 * EXP(-B2*L2*(L2+1)
C     &   * 1.4388/tmol(kt))
C     &  + 0.33333 + L2*(1.-L2*(4.444E-3+L2/3000.)) ) * 2. * Qvib(kt)
C      END IF
C1200  continue
C Here we approximate it with:
C       Q_vr(kt)=2.*(0.69493*tmol(kt)/B(00000)+0.33333) * Qvib(kt)
C where the lower level molecular B-value, B(v), is set to B(v=00000).
C In the C2H2 program B(v) is computed for each vibrational level as
C      SUMB=A1*V1+A2*V2+A3*V3+A4*V4+A5*V5
C      B(V1,V2,V3,V4,V5)=B0-SUMB
C      B2=B(VN,VO,VP,VQ,VR)
C      D2=1.E-6*B2*B2*B2
C  with B0,A1,A2,A3,A4,A5 = 1.1765 6.85-3 6.22-3 4.19-3 -1.31-3 -2.15-3

C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
       if (linesc2h2_direct.eq.1 .and. lhitr_c2h2.eq.0) then
        write(6,*) ' we read directly from the line list'
        close(unit=65)
        close(unit=66)
        open (unit=66,status='old',readonly,
     &   file='/ste1/uffegj/c2h2/linelist.dat')
        go to 990
       end if

       IF (lhitr_c2h2.eq.1) THEN
        write(6,*) ' we use Hitran linelist, and conv by Q(296)=',qh296
        close(unit=65)
        open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/c2h2/hitr96_c2h2.dat')
C       read(65,110) adum
       END IF
       nc2h2=0
       ngfz=0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13'
       write(6,*) ' v0,gf,exs,ivu,ivl,jl,ibr,isoc '
       write(6,*) ' ivu,ivl,v0 = '
       do 162 i=1,100000000
        IF (lhitr_c2h2.eq.1) THEN
C        read(65,*,end=163,err=1166) v0, exs, s0, ivu,ivl
         read(65,1167,end=163,err=1166) mol, iso, v0, s0, rdb, gamair,
     &    gammeself, exs, coefair, deltaair, ivu,ivl
C    &    qu, ql, ierr1, ierr2, ierr3, iref1, iref2, iref3
C iso: 1~12C2H2, 2~12C13CH2
1167  format(i2,i1,f12.6,d10.3,d10.3,f5.4,f5.4,f10.4,f4.2,f8.6,
     +i3,i3)
C    +,a9,a9,i1,i1,i1,i2,i2,i2)

C  604.774170      .4461   .292E-25   7   5               P 50
        gf = 1.87572E-12*s0 * qh296 / exp(-exs*1.4388/296.) /
     &     (1.-exp(-v0*1.4388/296.) ) * ana
        if (iso.ne.1) go to 162
        ELSE
        read(65,*,end=163) v0,jl,iv1l,iv2l,iv3l,iv4l,iv5l,l4l,l5l
     &              ,ju,iv1u,iv2u,iv3u,iv4u,iv5u,l4u,l5u,gf,exs,b2
C  734.727  0 0 0 0 0 0 0 0  1 0 0 0 0 1 0 1 7.334E-10 0.000E+00 1.170E+00
       nc2h2 = nc2h2 + 1
        if (i.le.5) write(6,165) v0,jl,iv1l,iv2l,iv3l,iv4l,iv5l,l4l,l5l
     &              ,ju,iv1u,iv2u,iv3u,iv4u,iv5u,l4u,l5u,gf,exs,b2
        END IF
        gfiso = gf 

        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 162
        end if
        write(66,116)v0,exs,gfiso
        go to 162
1166    continue
        if (i.le.10) write(6,*) ivu,ivl,v0
        backspace(1)
        read(65,110) adum
        write(6,1029) i,adum
162    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
c      go to 999
163    continue
       write(6,*) ' Last line of line strength data:'
       write(6,165) v0,jl,iv1l,iv2l,iv3l,iv4l,iv5l,l4l,l5l
     &              ,ju,iv1u,iv2u,iv3u,iv4u,iv5u,l4u,l5u,gf,exs,b2
       write(6,*)' I read 12C2H2; wrote total, #lines:'
       write(6,109) nc2h2,nc2h2-ngfz
165    format(f10.3,2(1x,i4,i2,i3,i2,4i3),1pe10.3,0pf9.2,f7.3)


C                                                                     ^
C                                         C2H2   C2H2    C2H2   C2H2  |
C----------------------------------------------------------------------
      GO TO 990
906   CONTINUE


      IF (MOLID.NE.'HCN ' .and. MOLID.NE.'HCNt') go to 907
C MOLID=HCN is the SCAN list, HCNt is Tennyson's list.
C----------------------------------------------------------------------
C                                         HCN    HCN     HCN    HCN   |
C                                                                     V

C    Calculate mean molecular weight
       wm1 = wgt_iso_h(1) * rel_iso_h(1)
       wm2 = wgt_iso_h(2) * rel_iso_h(2)
       wgtmolh = wm1 + wm2
       wm12 = wgt_iso_c(1) * rel_iso_c(1)
       wm13 = wgt_iso_c(2) * rel_iso_c(2)
       wgtmolc = wm12 + wm13
       wm14 = wgt_iso_n(1) * rel_iso_n(1)
       wm15 = wgt_iso_n(2) * rel_iso_n(2)
       wgtmoln = wm14 + wm15
       wgtmola = wgtmolh + wgtmolc + wgtmoln
       write(6,*) 'mean molecular weight for HCN is',wgtmola
       if (kiso_c.eq.1)
     *    write(18,*) ' only HCN lines of 12C are considered'
       if (kiso_c.eq.2)
     *    write(18,*) ' HCN lines of 12C and 13C are considered'
       if (kiso_n.eq.1)
     *    write(18,*) ' only HCN lines of 14N are considered'
       if (kiso_h.eq.1)
     *    write(18,*) ' only HCN lines of 1H are considered'
         do 1771 kc=1,kiso_c
1771     reliso(kc) = rel_iso_c(kc)
       kiso = kiso_c

      IF (MOLID.EQ.'HCN ') then
       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/hcn/linelist.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/hcn/hcn.partfunc')
      ELSE IF (MOLID.EQ.'HCNt') then
       open (unit=65,status='old',readonly,
     &   file='/ste4/uffegj/hcnt/linelist.dat')
C    &   file='/ste4/uffegj/hcnt/line1.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste4/uffegj/hcnt/hcn.partfunc')
      END IF
C
       B2 = 1.4782      !(from Herzberg)

       read(67,110) adum
       read(67,110) adum
       read(67,110) adum
       kqt = 0
       do 171 k=1,1000
      IF (MOLID.EQ.'HCN ') then
       read(67,*,end=179) tin(k),qhanal,qvin(k),qmol360,qhmm,qmolmm
       qvin(k) = (0.69493*tin(k)/B2+0.33333) * qvin(k)
      ELSE IF (MOLID.EQ.'HCNt') then
       read(67,*,end=179) tin(k),qvin(k)
      END IF
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
171    continue
179    continue
       write(6,*) ' I read',kqt
     &     ,' vibrational partition function values for HCN'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)
       write(6,*) ' The product of vibrational and rotational '
     &   ,' Q, Q_vr(T), is computed as:'
       do 1791 kt=1,kqt
1791   write(6,'(f8.1,f12.2)') tin(kt),qvin(kt)
       qh296 = qvin(1)
C 
C the vibrational partition function is calculated in the program part.f
C while the rotational Q is calculated in HCN.f as:
C     IF(VI.LE.M1+MD1.AND.V2.LE.M2+MD2.AND.V3.LE.M3+MD3) THEN
C     B1=B(VI,VJ,VK)
C     B2=B(VL,VM,VN)
C     Q1=QS(VI,VJ,VK)  
C     Q2=QS(VL,VM,VN)  
C     ELSE 
C      B1=B(0,0,1) 
C      B2=B(0,0,0) 
C      Q1=QS(0,0,1)
C      Q2=QS(0,0,0)
C     ENDIF
C     D1=0.95E-6*B1*B1*B1  
C     D2=0.95E-6*B2*B2*B2  
C  
C     t = tmol(itlim)
C     IF(L2.EQ.0) THEN 
C       CPARTR=1./(0.6951*T/B2+0.33333)
C     ELSE 
C       CPARTR=1./(0.6951*T/B2*EXP(-B2*L2*(L2+1)*1.439/T)+0.33333+
C    &   L2*(1.-L2*(4.4444E-3+L2/3000.)))
C     END IF
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
       if (lineshcn_direct.eq.1 .and. lhitr_hcn.eq.0) then
        write(6,*) ' we read directly from the line list'
        close(unit=65)
        close(unit=66)
           IF (MOLID.EQ.'HCN ') then
            open (unit=66,status='old',readonly,
     &        file='/ste1/uffegj/hcn/linelist.dat')
           ELSE IF (MOLID.EQ.'HCNt') then
            open (unit=66,status='old',readonly,
     &        file='/ste4/uffegj/hcnt/linelist.dat')
C    &        file='/ste4/uffegj/hcnt/line1.dat')
           END IF
        go to 990
       end if

       IF (lhitr_hcn.eq.1) THEN
        write(6,*) ' we use Hitran linelist, and conv by Q(296)=',qh296
        close(unit=65)
        open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/hcn/hitr96_hcn.dat')
C       open (unit=65,status='old',readonly,
C    &   file='/ste1/aringer/molec/hitran96/23_hit96.par')
C       read(65,110) adum
C       read(65,110) adum
       END IF
       nhcn=0
       ngfz=0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' iel,ivu,ivl,jl,ibr,gf,v0,exs,dw13'
       write(6,*) ' v0,gf,exs,ivu,ivl,jl,ibr,isoc '
       write(6,*) ' ivu,ivl,v0 = '
       do 172 i=1,100000000

        IF (lhitr_hcn.eq.1) THEN
C        read(65,*,end=173,err=1176) v0, exs, s0, ivu,ivl
         read(65,1177,end=173,err=1176) mol, iso, v0, s0, rdb, gamair,
     &    gammeself, exs, coefair, deltaair, ivu,ivl
     &    ,br,jhlow
C    &    ,qu, ql, ierr1, ierr2, ierr3, iref1, iref2, iref3
C iso: 1~H12CN, 2~H13CN
1177  format(i2,i1,f12.6,d10.3,d10.3,f5.4,f5.4,f10.4,f4.2,f8.6,
     +i3,i3,13x,a1,i3)
C    +,a9,a9,i1,i1,i1,i2,i2,i2)

C  604.774170      .4461   .292E-25   7   5               P 50
        gf = 1.87572E-12*s0 * qh296 / exp(-exs*1.4388/296.) /
     &     (1.-exp(-v0*1.4388/296.) ) * ana
       if (kiso_c.eq.1) then
          if (iso.ne.1) go to 172
          gfiso = gf
       else if (kiso_c.eq.2) then
          if (iso.ne.1 .and. iso.ne.2) go to 172
          gfiso = gf * reliso(iso)
        end if

        ELSE

        read(65,*,end=173) v0,jl,iv1l,iv2l,iv3l,iv4l,iv5l,l4l,l5l
     &              ,ju,iv1u,iv2u,iv3u,iv4u,iv5u,l4u,l5u,gf,exs,b2
C  734.727  0 0 0 0 0 0 0 0  1 0 0 0 0 1 0 1 7.334E-10 0.000E+00 1.170E+00
        gfiso = gf 
        nhcn = nhcn + 1
        if (i.le.5) write(6,165) v0,gf,exs
        END IF

        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 172
        end if
        write(66,116)v0,exs,gfiso
        go to 172
1176    continue
        if (i.le.10) write(6,*) ivu,ivl,v0
        backspace(1)
        read(65,110) adum
        write(6,1029) i,adum
172    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
c      go to 999
173    continue
       write(6,*) ' Last line of line strength data:'
       write(6,175) v0,gf,exs
       write(6,*)' I read 12HCN; wrote total, #lines:'
       write(6,109) nhcn,nhcn-ngfz
175    format(f10.3,1pe10.3,0pf9.2)

C                                                                     ^
C                                         HCN    HCN     HCN    HCN   |
C----------------------------------------------------------------------
      GO TO 990
907   CONTINUE


      IF(MOLID.NE.'H2O '                       !Scan H2O in a1a1,a1b2,b2b2.dat form
     *    .and.MOLID.NE.'H2Os'                 !Scan H2O in a xxx.dat form
     *    .and.MOLID.NE.'H2Oh'                 !Hitran line list in the form gf_hitlist.dat
     *    .and.MOLID.NE.'H2Ot'                 !Hitran line list in the form gf_hitemp.dat
     *    .and.MOLID.NE.'swmn'                 !small Swenke line list
     *                    .and.molid.ne.'swmx') goto 908                 !large Swenke line list
C----------------------------------------------------------------------
C                                         H2O    H2O     H2O    H2O   |
C                                                                     V

C    Calculate mean molecular weight
          wgtmol_2h = 2.*wgt_iso_h(1)*rel_iso_h(1)**2 +
     &                2.*wgt_iso_h(2)*rel_iso_h(2)**2 +
     &             (wgt_iso_h(1)+wgt_iso_h(2))*rel_iso_h(1)*rel_iso_h(2)
       wgtmola = wgtmol_2h
       do 187 ko=1,kiso_o
          wgtmola=wgt_iso_o(ko)*rel_iso_o(ko) + wgtmola
187    continue
       if (kiso_h .eq. 1) then
         do 1871 ko=1,kiso_o
1871     reliso(ko) = rel_iso_o(ko)
       end if
       kiso = kiso_o
       write(18,*) 'mean molecular weight for H2O is',wgtmola
       write(18,*) 'the relative abundance of H2O, D2O, HDO is:'
       write(18,*) rel_iso_h(1)**2, rel_iso_h(2)**2,
     &                              rel_iso_h(1)*rel_iso_h(2)


C small Swenke line list:
       if(molid.eq.'swmn')
     & open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/h2o/h2osw.dat')
C Scan H2O line list not on a1b2,b2b2,a1a1.dat form:
C      if(molid.eq.'H2O ')
C    & open (unit=65,status='old',readonly,
C    &   file='/ste2/uffegj/scanh2o_vo7.dat')
       open (unit=67,status='old',readonly,
     &   file='/ast/p2/uffegj/2013/h2o/scanh2o_vo7.pfct')

         pi         = 3.1415926536D0
         h_Planck   = 6.6260755D-34
         c_light    = 2.99792458D10
         Bolzm_k    = 1.380658D-23
         Avogadro_N = 6.0221367D23
         hc_div_k   = h_Planck*c_light/Bolzm_k
         A_const    = 8.*pi**3*Avogadro_N*1.D-48/(3.*h_Planck*c_light)
         RC01 = 27.88063
         RC02 = 14.52177
         RC03 =  9.27771

       read(67,110) adum
C T[K] Q-morbid  Q-harmonic  Q-anh_Morbid  Q-anh_30  Q-anh_40
       kqt = 0
       do 181 k=1,1000
       read(67,*,end=189) tin(k),qvin(k),q1,q2,q3,q4
C incl. temp.dependent part of Q_rot; temp.indep. part is incl. in S0 on linelist

C rotational partition function at  T = tin(k)
          temp = tin(k)
          hc_kt = hc_div_k / temp
          Qra = RC01*RC02*RC03
          Qra = sqrt(pi/Qra*(temp/hc_div_k)**3)
          A = RC01
          B = sqrt(RC02*RC03)
          Qrb = B*hc_kt
          B = (1.D0-B/A)*Qrb
          Qrb = Qra*exp(0.25D0*Qrb)*(1.D0+B/12.D0+B*B*7.D0/480.D0)
          Qrb = Qrb + Qrb
      
C      write(6,*)' T, Qv, Qr, Qv*Qr = ',temp,qvin(k),Qrb,qvin(k)*Qrb

C      qvin(k) = qvin(k) * tin(k)**1.5
       qvin(k) = qvin(k) * Qrb
       kqt = kqt + 1
       if(kqt.eq.1) Qrb1 = Qrb
       if (kqt.gt.nqt) stop ' kqt > nqt'
181    continue
189    continue
       write(6,*) ' I read',kqt,' partition function values'
       write(6,1891) tin(1),qvin(1),Qrb1,tin(kqt),qvin(kqt),Qrb
1891   format
     & (' from t/qvr/Qrb=',F7.0,2F11.1,' to t/qvr/Qrb=',F7.0,2F11.1)

C If line list is Scan lines in form a1a1,a1b2,b2b2.dat (H2O) or one file (H2Os), or
C in the form of Swenkes 3.e8 line list (swmx; seperate call to subr.swmx) :
       if (MOLID.eq.'H2O ' .or. MOLID.eq.'H2Os' .or.
     &                               molid.eq.'swmx') go to 990


C If the line list is Hitran in the gf format (gf_hitlist.dat made with hitinf.f):
       if (MOLID.eq.'H2Oh') then
       write(6,*)' we read directly from Hitran gf_hitlist.dat linelist'
        close(unit=66)
        open (unit=66,status='old',
     &   file='/ste1/uffegj/h2o/gf_hitlist.dat')
        go to 990
       end if


C If the line list is Hitemp in the gf format (gf_hitemp.dat made with ../h20/hitinf.f):
       if (MOLID.eq.'H2Ot') then
       write(6,*)' we read directly from Hitran gf_hitemp.dat linelist'
        close(unit=66)
        open (unit=66,status='old',
     &   file='/ste1/uffegj/h2o/gf_hitemp.dat')
        go to 990
       end if


C Read from small Swenke line list or from Scan H2O list,
C i.e.,    molid.eq.'swmn' .or. molid.eq.'H2Os'  :
       nho2=0
       write(6,*) ' First 5 lines of line strength data:'
       do 182 i=1,100000000
        read(65,*,end=183) v0,s0,exs
        if (i.le.5) write(6,185) v0,exs,s0
        gfiso = s0*1.87572D-7    !1997,10 no inf on H2O isotopes
        write(66,116)v0,exs,gfiso
        nh2o = nh2o + 1
182    continue
       write(6,*)' Stop: more than 1.e8 lines in input file ??'
       go to 999
183    continue
       write(6,*) ' Last line of line strength data:'
       write(6,185)v0,exs,s0
185    format(2f11.4,1pe13.5)
       close (unit=65)
C                                                                     ^
C                                         H2O    H2O     H2O    H2O   |
C----------------------------------------------------------------------
      GO TO 990
908   CONTINUE


      IF (MOLID.NE.'TiO ' .and.molid.ne.'TIO '
     &                    .and.molid.ne.'TiOs') go to 909   !s~Schwenke's list
C----------------------------------------------------------------------
C                                         TiO    TiO     TiO    TiO   |
C                                                                     V
C
C    Calculate mean molecular weight
       wgtmola = 0.0
       kiso_o = 1           !only 16O is considered in TiO line list
       do 1097 ko=1,kiso_o
       do 1097 kt=1,kiso_t
          relabu = rel_iso_o(ko)*rel_iso_t(kt)
          wgtmola=(wgt_iso_o(ko)+wgt_iso_t(kt))*relabu + wgtmola
1097   continue
       if (kiso_o .eq. 1) then
         do 1971 kt=1,kiso_t
1971     reliso(kt) = rel_iso_t(kt)
       end if
       kiso = kiso_t
       write(18,*) 'mean molecular weight for TiO is',wgtmola

       open (unit=67,status='old',readonly,
     &   file='/p7/hpg688/hpg688/os/data/tio/tio.partfunc')
C
       read(67,110) adum
       kqt = 0
       do 191 k=1,1000
       read(67,*,end=1099) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
191   continue
1099   continue
       write(6,*) ' I read',kqt,' partition function values for TiO'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)


       write(6,*) ' I adopted the following f_el(nu_0) values for'
     &       ,' the 7 el-systems'
       write(6,*)'--alpha, beta, gamma, gammap, delta, phi, epsilon:'
       write(6,10910) (ftio(itf),itf=1,7)
10910  format (7f7.3)

       if (molid.eq.'TiOs') go to 530

       write(6,*) ' First 5 lines of Scan line list data:'
       write(6,*) ' el.system, v0, gf, exs, id '
C
       DO 300 LW=1,7
       write(6,*) ' TiO electronic system number: ',lw
       IF (LW.EQ.1) THEN
         SGFCONV = ftio(1)/19339.0
         OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/alpha.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.2) THEN
        SGFCONV = ftio(2)/17840.6
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/beta.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.3) THEN
        SGFCONV = ftio(3)/14092.9
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/gamma.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.4) THEN
        SGFCONV = ftio(4)/16147.0
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/gprime.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.5) THEN
        SGFCONV = ftio(5)/11273.3
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/delta.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.6) THEN
        SGFCONV = ftio(6)/9054.0
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/phi.dat',STATUS='OLD',READONLY)
       ELSE IF (LW.EQ.7) THEN
        SGFCONV = ftio(7)/11893.9
        CLOSE(65)
        OPEN(UNIT=65,
     & FILE='/p8/hpg688/os/data/tio/epsilon.dat',STATUS='OLD',READONLY)
       END IF
C
C
       nisotio = kiso_o * kiso_t
       ntio=0
       do 580 i=1,kiso_t
580    linesgfz(i) = 0
       READ (65,15) DUMMY
       READ (65,*)  NLIN
15     FORMAT (A55)


       DO 600 I=1,100000000
C
        READ(65,*,END=699) WS,EXS,GF,iel,ivl,ivu,jl,ju,ibr
     &      ,DW46,DW47,DW49,DW50
        if (i.le.5) 
     &     write(6,1095) lw,ws,gf,exs,iel,ivl,ivu,jl,ju,ibr
        ntio = ntio + 1
        WSISO(1) = WS
        WSISO(2) = WS +DW46
        WSISO(3) = WS +DW47
        WSISO(4) = WS +DW49
        WSISO(5) = WS +DW50
C
C
       DO 610  IWS = 1, kiso_t
       WS = WSISO(IWS)
C
       gfiso = gf * rel_iso_t(iws) * SGFCONV * WS
C
        if (gfiso .eq. 0.D0) then
           linesgfz(iws) = linesgfz(iws) + 1
           go to 610
        end if
        write(66,116)ws,exs,gfiso

610    CONTINUE    ! next isotop
600    CONTINUE    ! next line
699    CONTINUE
       write(6,*) ' Last line of line strength data:'
       write(6,1095) lw,ws,gf,exs,iel,ivl,ivu,jl,ju,ibr
1095   format(i3,f11.3,1pe11.3,0pf10.2,i2,2i3,2i4,i3)

       write(6,*) ' we read ',ntio,' (*5) lines of this el.system'
       write(6,*) ' of them, for each isotopes, gf = 0 for:'
       write(6,19) (linesgfz(iwsw),iwsw=1,kiso_t)
300    CONTINUE    ! consider the next electronic system
C
       go to 531

530    continue


       write(6,*) ' First 5 lines of Schwenke line list:'
       write(6,*) ' el.system, v0, gf, exs, id '
C
       OPEN (unit=65,
     & FILE='/p8/hpg688/os/data/tio/tio.dat',STATUS='OLD',READONLY)

531    continue
19     format(10i8)
C
C                                         TiO    TiO     TiO    TiO   |
C----------------------------------------------------------------------
      GO TO 990
909   CONTINUE


      IF (MOLID.NE.'SiO ' .and.molid.ne.'SIO') go to 910
C----------------------------------------------------------------------
C                                         SiO    SiO     SiO    SiO   |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       kiso_o = 1           !only 16O is considered in SiO line list
       isoo = 1
       do 1107 ko=1,kiso_o
       do 1107 ks=1,kiso_si
          relabu = rel_iso_o(ko)*rel_iso_si(ks)
          wgtmola=(wgt_iso_o(ko)+wgt_iso_si(ks))*relabu + wgtmola
1107   continue
       if (kiso_o .eq. 1) then
         do 11071 ks=1,kiso_si
11071    reliso(ks) = rel_iso_si(ks)
       end if
       kiso = kiso_si
       write(18,*) 'mean molecular weight for SiO is',wgtmola

       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/sio/langhoff_sio.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/sio/langhoff.partfunc')
C
       read(67,110) adum
       kqt = 0
       do 11090 k=1,1000
       read(67,*,end=1109) tin(k),qvin(k),qv29,qv30
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
11090   continue
1109   continue
       write(6,*) ' I read',kqt,' partition function values for SiO'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)


       nsio=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0,gf,exs ivl,ivu, jl, ju, iso: '
       do 1102 i=1,100000000
       read(65,*,end=1103) v0, exs, gf, ivl, ivu, jl, ju, iso
        if (i.le.5) write(6,1105) v0,gf,exs,ivl,ivu,jl,ju,iso
        nsio = nsio + 1
        isos = iso - 27
        gfiso = gf * rel_iso_o(isoo) * rel_iso_si(isos)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 1102
        end if
        write(66,116)v0,exs,gfiso
1102    continue
       write(6,*)' Stop: more than 1.e8 lines in SiO input file ??'
c      go to 999
1103    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1105) v0,gf,exs,ivl,ivu,jl,ju,iso
1105   format(f11.3,1pe12.3,0pf11.2,2i3,2i4,i5)
1101   format(' read total',i6,'; wrote total',i6,' SiO lines.')
       write(6,1101) nsio,nsio-ngfz

C                                                                     ^
C                                         SiO    SiO     SiO    SiO   |
C----------------------------------------------------------------------
      GO TO 990
910   CONTINUE

      IF (MOLID.NE.'OH  ') go to 911
C----------------------------------------------------------------------
C                                         OH     OH      OH     OH    |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       do 1117 ko=1,kiso_o
       do 1117 kh=1,kiso_h
          relabu = rel_iso_o(ko)*rel_iso_h(kh)
          wgtmola=(wgt_iso_o(ko)+wgt_iso_h(kh))*relabu + wgtmola
1117   continue
       if (kiso_h .eq. 1) then
         do 11171 ko=1,kiso_o
11171    reliso(ko) = rel_iso_o(ko)
       end if
       kiso = kiso_o
       write(18,*) 'mean molecular weight for OH is',wgtmola

       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/oh/oh_schw.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/oh/oh.partfunc')

C
       read(67,110) adum
       kqt = 0
       do 1111 k=1,1000
       read(67,*,end=1119) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
1111   continue
1119   continue
       write(6,*) ' I read',kqt,' partition function values for OH'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)

       noh=0
       n2iso = 0    !#lines with iso .ne. 16O1H
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 9 lines of line strength data:'
       write(6,*) ' v0,gf,exs ielu,ivu,ju,ifstru,iparitu ',
     &       ' iell,ivl,jl,ifstrl,iparitl, iso'
       do 1112 i=1,100000000
       read(65,*,end=1113) iso, v0, gf, dip2, exs, aju, ajl
     &     ,ivu, ielu, ifstru,   ivl, iell, ifstrl
     &     ,irootnumf,  irootnumi
     &     ,iparitu,   iparitl
C    irootnumf,  irootnumi   = ju, jl  ?
C    angmomf, angmomi = ju, jl ?
        if (i.le.9) write(6,1115) v0,gf,exs
     &       ,ielu,ivu,aju,ifstru,iparitu
     &       ,iell,ivl,ajl,ifstrl,iparitl, iso
        noh = noh + 1
        if (iso.eq.11) then
           isoo = 1
           isoh = 1
        end if
        if (iso.ne.11) n2iso = n2iso + 1
        gfiso = gf * rel_iso_o(isoo) * rel_iso_h(isoh)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 1112
        end if
        write(66,116)v0,exs,gfiso
1112    continue
       write(6,*)' Stop: more than 1.e8 lines in OH input file ??'
c      go to 999
1113    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1115) v0,gf,exs
     &       ,ielu,ivu,aju,ifstru,iparitu
     &       ,iell,ivl,ajl,ifstrl,iparitl, iso
1115   format(f10.3,1pe10.3,0pf9.2,2(2x,2i3,f5.1,2i3),i4)
1091   format(' read total',i6,'; wrote total',i6,' 16O1H lines.'
     &   ,/' There were',i3,' lines with isotop .ne. 16O1H')
       write(6,1091) noh,noh-ngfz,n2iso

       reliso(1) = 1.d0
       kiso = 1


C                                                                     ^
C                                         OH     OH      OH     OH    |
C----------------------------------------------------------------------
      GO TO 990
911   CONTINUE


      IF (MOLID.NE.'CH4 ' .and. MOLID.NE.'CH4h') go to 912
C----------------------------------------------------------------------
C                                         CH4    CH4     CH4    CH4   |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       kiso_h = 1
       kiso_c = 1
       relabu = 1.
       do 1127 kc=1,kiso_c
       do 1127 kh=1,kiso_h
C         relabu = rel_iso_c(kc)*rel_iso_h(kh)
          wgtmola=(wgt_iso_c(kc)+4*wgt_iso_h(kh))*relabu + wgtmola
1127   continue
C       if (kiso_h .eq. 1) then
C         do 11271 kc=1,kiso_c
C11271    reliso(kc) = rel_iso_c(kc)
C       end if
       kiso = 1
       reliso(1) = 1.d0
       write(18,*) ' mean molecular weight for CH4 is',wgtmola
       write(18,*) ' (only 12C1H4 is available from input) '

      IF (MOLID.EQ.'CH4 ') then
       open (unit=65,status='old',readonly,
     &   file='/ste4/uffegj/ch4/jpclist.dat')
       tstds = 1300.d0  !the temperature at which the line intensity in cm^-2 atm^-1 is given
      ELSE
       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/ch4/12CH4_nils_uffe.linelist')
      END IF
C    &   file='/ste1/uffegj/ch4/hitran.dat')

       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/ch4/ch4.partfunc')

C
C      read(67,110) adum
       kqt = 0
       do 1121 k=1,1000
       read(67,*,end=1129) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
1121   continue
1129   continue
C      IF (MOLID.EQ.'CH4t') then
C        do 11211 k=1,kqt
C11211   qvin(k) = 1.d0
C      END IF
       write(6,*) ' I read',kqt,' partition function values for CH4'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)
C
      IF (MOLID.EQ.'CH4 ') then
      DO 3021 IT=1,KQT
      QTLN(IT)=LOG(qvin(it))
3021  continue
      CALL TB04AT(NQT,KQT,tin,QTLN,DQDT,W)
      II = -1
      TINT=tstds
      qtstds = TG01BT(II,NQT,KQT,tin,QTLN,DQDT,TINT)
      qtstds = EXP(qtstds)
      write(6,3121) molid,tstds,qtstds
3121  format( ' molecule, temp, Q(temp): ',a4,f8.1,1pe11.3)
      end if

       nch4=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0,gf,exs ielu,ivu,ju,ifstru,iparitu '
C  Line list from Hitran data base 1984, 29669 lines, 12C1H
C  wavenumber    E_low    So(296K)
C  [cm-1]       [cm-1]    [cm/molecule]
C 1007.173170  1977.2151   .664E-25   2   1 18  0 F2  2   19  0 F1  5
C 1007.606950  1977.1796   .423E-25   2   1 18  0 F1  2   19  0 F2  5
CJag laeser in foeljande [fra Hitran]:
C          read(10,27,end=99) mo,iso,wave,S,R,agam,sgam,chie,n,d,
C    &      v1,v2,jp,rp,cp,np,symbp,jb,rb,cb,nb,symbb,ierf,
C    &      iers,ierh,ireff,irefs,irefh
C              write(outfil,22) wave,chie,S,v1,v2,jp,rp,cp,np,symp,
C    &                          jb,rb,cb,nb,symb
C22    format(f12.6,x,f10.4,x,e10.3,x,i3,x,i3,x,i2,x,i2,x,a2,x,i2,
C    +x,a1,x,i2,x,i2,x,a2,x,i2,x,a1)
C GF=1.87572E-12 Sl[cm/mol] QvQr/exp(-[Gv+Fj]hc/kT)(1-e-nuhc/kT)
C GF = 1.87572E-12*S0
C The partition function read in above (harmnonic for T<1200K, 16*Irwin81
C for T>1200K) is consistent with the conversion below from Hitran data
C to gf values, such that we get correct opacity in cm2/g* for CH4

      IF (MOLID.EQ.'CH4h') then
       read(65,110) adum
       read(65,110) adum
       read(65,110) adum
      END IF
C1214   format(f12.6,f11.4,e11.3,i4,i4,i3,i3,1x,a2,i3,
C     +        1x,a1,i3,i3,1x,a2,i3,1x,a1)
1214   format(f12.6,f11.4,e11.3,i4,i4,i3)

       do 1122 i=1,100000000
      IF (MOLID.EQ.'CH4 ') then
       read(65,*,end=1123,err=1126) v0, s0, exs
       s00 = s0
       s0 = s0 * 82.056 * tstds
       s0 = s0 * qtstds / exp(-exs*1.4388/tstds) /
     &     (1.-exp(-v0*1.4388/tstds) )
       gf = 1.87572E-12*s0 
      ELSE
       read(65,1214,end=1123,err=1126) v0, exs, s0, ivu,ivl, ju
       gf = 1.87572E-12*s0 * 587. / exp(-exs*1.4388/296.) /
     &     (1.-exp(-v0*1.4388/296.) ) * ana
      END IF
C       read(65,1214,end=1123) v0, exs, s0, ivu,ivl, ju,ru,cu,nu,symu
C     *                     ,jl,rl,cl,nl,syml

1215   format(f14.6,1pe10.3,0pf14.6,1p2e10.3)
        if (i.le.5) write(6,1215) v0, s00, exs, s0, gf
C        if (i.le.5) write(6,1214) 
C     *                      v0, exs, s0, ivu,ivl, ju,ru,cu,nu,symu
C     *                     ,jl,rl,cl,nl,syml
        nch4 = nch4 + 1
C in Hitran we read S0 in units cm/molecule at 296K. We therefore 
C multiply S0 with Q(296K), and then divide with the Boltzman factor and the
C stimulated emission factor to get the temperature independent line strength
C in cm/molecule. We then multiply with N_a (ana) to get in cm/mol, and finally
C GF=1.87572E-12 Sl[cm/mol] QvQr/exp(-[Gv+Fj]hc/kT)(1-e-nuhc/kT)
C       gfiso = gf * rel_iso_c(isoc)
        if (gf .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 1122
        end if
        write(66,116)v0,exs,gf
        go to 1122
1126    continue
        backspace(1)
        read(65,110) adum
        write(6,1029) i,adum
1122    continue
       write(6,*)' Stop: more than 1.e8 lines in CH4 input file ??'
c      go to 999
1123    continue
1029   format(' Problems at line',i3,' :'a66)
       write(6,*) ' Last line of line strength data:'
        write(6,1215) v0, s00, exs, s0, gf
C      write(6,1125) v0,gf,exs
C    &       ,ielu,ivu,aju,ifstru,iparitu
C    &       ,iell,ivl,ajl,ifstrl,iparitl, iso
1125   format(f10.3,1pe10.3,0pf9.2,2(2x,2i3,f5.1,2i3),i4)
1092   format(' read total',i8,'; wrote total',i8,' 12C1H4 lines.')
       write(6,1092) nch4,nch4-ngfz

C                                                                     ^
C                                         CH4    CH4     CH4    CH4   |
C----------------------------------------------------------------------
      GO TO 990
912   CONTINUE



      IF (MOLID.NE.'CS  ') go to 913
C----------------------------------------------------------------------
C                                         CS     CS      CS     CS    |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       sumabu = 0.0
C      kiso_c = 1           !isotopes 12C,13C,32S,33S,34S
       do 1137 kc=1,kiso_c
       do 1137 ks=1,kiso_s
          kss = ks
          if(kc.eq.2) then
            if(ks.ge.2) go to 1137
            kss = 4
          end if
          relabu = rel_iso_c(kc)*rel_iso_s(ks)
          reliso(kss) = relabu
          sumabu = sumabu + relabu
1137   continue
       do 11371 kc=1,kiso_c
       do 11371 ks=1,kiso_s
          kss = ks
          if(kc.eq.2) then
            if(ks.ge.2) go to 11371
            kss = 4
          end if
          reliso(kss) = reliso(kss)/sumabu
          wgtmola=(wgt_iso_c(kc)+wgt_iso_s(ks))*reliso(kss) + wgtmola
11371  continue
       kiso = kss
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
       write(18,*)
     *   'mean molecular weight for',kss,' isotopes of CS is',wgtmola
       write(18,*)
     *  ' reliso(1-4)=',reliso(1),reliso(2),reliso(3),reliso(4)

       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/cs/cs.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/cs/cs.partfunc')
C
       read(67,110) adum
       read(67,110) adum
       read(67,110) adum
       kqt = 0
       do 11390 k=1,1000
       read(67,*,end=1139) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
11390   continue
1139   continue
       write(6,*) ' I read',kqt,' partition function values for CS'
       write(6,*) ' from t=',tin(1),' to t=',tin(kqt)


       ncs=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0,gf,exs ivl,ivu, jl, ju,iso(1-4): '
       do 1132 i=1,100000000
       read(65,*,end=1133) v0, gf, exs, ivl, ivu, jl, ju, iso
        if (i.le.5) write(6,1135) v0,gf,exs,ivl,ivu,jl,ju,iso
        ncs = ncs + 1
        gfiso = gf * reliso(iso)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           go to 1132
        end if
        write(66,116)v0,exs,gfiso
1132    continue
       write(6,*)' Stop: more than 1.e8 lines in CS input file ??'
c      go to 999
1133    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1135) v0,gf,exs,ivl,ivu,jl,ju,iso
1135   format(f11.3,1pe12.3,0pf11.2,2i3,2i4,i5)
1131   format(' read total',i6,'; wrote total',i6,' CS lines.')
       write(6,1131) ncs,ncs-ngfz

C                                                                     ^
C                                         CS     CS      CS     CS    |
C----------------------------------------------------------------------
      GO TO 990
913   CONTINUE




      IF (MOLID.NE.'FeH ' .and. MOLID.NE.'FEH ') go to 914
C----------------------------------------------------------------------
C                                         FeH    FeH     FeH    FeH   |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       sumabu = 0.0

C   isotopes {54,56,57,58}Fe ~ 5.8, 91.7, 2.2, and 0.28%
C Since for the moment we have only absorption coefficient data for 
C 56Fe-1H, we solve the problem by skipping the loop and setting
C wgtmola = ( 55.934939 + 1.00782504 ) * 1.00 = 56.942764
       do 1147 kf=1,kiso_fe
       do 1147 kh=1,kiso_h
          relabu = rel_iso_fe(kf)*rel_iso_h(kh)
          sumabu = sumabu + relabu
1147   continue
          kiso = 0
       do 11471 kf=1,kiso_fe
       do 11471 kh=1,kiso_h
          kiso = kiso + 1
          reliso(kiso) = rel_iso_fe(kf)*rel_iso_h(kh)/sumabu
          wgtmola=(wgt_iso_fe(kf)+wgt_iso_h(kh))*reliso(kiso) + wgtmola
11471  continue
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

       kiso = 1
       reliso(1) = 1.d0
       wgtmola = 56.942764
       write(18,*) ' mean molecular weight for FeH is',wgtmola
       write(18,*) ' (only 56Fe1H is available from input) '


       open (unit=65,status='old',readonly,
     &   file='/ste4/uffegj/feh/dulik.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste4/uffegj/feh/feh.partfunc')
C
       read(67,110) adum
       read(67,110) adum
       read(67,110) adum
       kqt = 0
       do 11490 k=1,1000
       read(67,*,end=1149) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
11490   continue
1149   continue
       write(6,*)'I read',kqt,' partition function values for FeH'
       write(6,*) 'from t=',tin(1),' to t=',tin(kqt)


       nfeh=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0 gf exs ivl ivu jl ju iso: '
       do 1142 i=1,10000000
       read(65,*,end=1143) v0, gf, exs, ivl, ivu, jl, ju, iso
        if (i.le.5) write(6,1145) v0,gf,exs,ivl,ivu,jl,ju,iso
        nfeh = nfeh + 1
        gfiso = gf * reliso(iso)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           if(ngfz.le.5)write(6,1145)v0,gf,exs,ivl,ivu,jl,ju,iso
           go to 1142
        end if
        write(66,116)v0,exs,gfiso
1142    continue
       write(6,*)' Stop: more than 1.e8 lines in FeH input file ??'
c      go to 999
1143    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1145) v0,gf,exs,ivl,ivu,jl,ju,iso
1145   format(f11.3,1pe12.3,0pf11.2,2i3,2i4,i5)
1141   format(' read total',i6,'; wrote total',i6,' FeH lines.')
       write(6,1141) nfeh,nfeh-ngfz

C                                                                     ^
C                                         FeH    FeH     FeH    FeH   |
C----------------------------------------------------------------------
      GO TO 990
914   CONTINUE
 


      IF (MOLID.NE.'CrH ' .and. MOLID.NE.'CrH ') go to 915
C----------------------------------------------------------------------
C                                         CrH    CrH     CrH    CrH   |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       sumabu = 0.0

C   isotopes {50,52,53,54}Cr ~ 4.35, 83.79, 9.50, and 2.36%
C Since for the moment we have only absorption coefficient data for 
C 52Cr-1H, we skip the loop and set
C wgtmola = ( 51.940510 + 1.00782504 ) * 1.00 = 52.948335
       do 1157 kc=1,kiso_cr
       do 1157 kh=1,kiso_h
          relabu = rel_iso_cr(kc)*rel_iso_h(kh)
          sumabu = sumabu + relabu
1157   continue
          kiso = 0
       do 11571 kc=1,kiso_cr
       do 11571 kh=1,kiso_h
          kiso = kiso + 1
          reliso(kiso) = rel_iso_cr(kc)*rel_iso_h(kh)/sumabu
          wgtmola=(wgt_iso_cr(kc)+wgt_iso_h(kh))*reliso(kiso) + wgtmola
11571  continue
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

       kiso = 1
       reliso(1) = 1.d0
       wgtmola = 52.948335
       write(18,*) ' mean molecular weight for CrH is',wgtmola
       write(18,*) ' (only 52Cr1H is available from input) '


       open (unit=65,status='old',readonly,
     &   file='/ste4/uffegj/crh/linelist.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste4/uffegj/crh/crh.partfunc')
C
       read(67,110) adum
       read(67,110) adum
       read(67,110) adum
       kqt = 0
       do 11590 k=1,1000
       read(67,*,end=1159) tin(k),qvin(k)
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
11590   continue
1159   continue
       write(6,*)'I read',kqt,' partition function values for CrH'
       write(6,*) 'from t=',tin(1),' to t=',tin(kqt)


       ncrh=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0 gf exs ivl ivu jl ju iso: '
       do 1152 i=1,10000000
       read(65,*,end=1153) v0, gf, exs, ivl, ivu, jl, ju, iso
        if (i.le.5) write(6,1155) v0,gf,exs,ivl,ivu,jl,ju,iso
        ncrh = ncrh + 1
        gfiso = gf * reliso(iso)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           if(ngfz.le.5)write(6,1155)v0,gf,exs,ivl,ivu,jl,ju,iso
           go to 1152
        end if
        write(66,116)v0,exs,gfiso
1152    continue
       write(6,*)' Stop: more than 1.e8 lines in CrH input file ??'
c      go to 999
1153    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1155) v0,gf,exs,ivl,ivu,jl,ju,iso
1155   format(f11.3,1pe12.3,0pf11.2,2i3,2i4,i5)
1151   format(' read total',i6,'; wrote total',i6,' CrH lines.')
       write(6,1151) ncrh,ncrh-ngfz

C                                                                     ^
C                                         CrH    CrH     CrH    CrH   |
C----------------------------------------------------------------------
      GO TO 990
915   CONTINUE
 


      IF (MOLID.NE.'MgH ' .and. MOLID.NE.'MGH ') go to 916
C----------------------------------------------------------------------
C                                         MgH    MgH     MgH    MgH   |
C                                                                     V
C    Calculate mean molecular weight
       wgtmola = 0.0
       sumabu = 0.0

C Since for the moment we have only absorption coefficient data for 
C 24Mg-1H, we skip the loop and set
C wgtmola = ( 23.985045 + 1.00782504 ) * 1.00 = 24.99287 
       do 1757 kc=1,kiso_mg
       do 1757 kh=1,kiso_h
          relabu = rel_iso_mg(kc)*rel_iso_h(kh)
          sumabu = sumabu + relabu
1757   continue
          kiso = 0
       do 17571 kc=1,kiso_mg
       do 17571 kh=1,kiso_h
          kiso = kiso + 1
          reliso(kiso) = rel_iso_mg(kc)*rel_iso_h(kh)/sumabu
          wgtmola=(wgt_iso_mg(kc)+wgt_iso_h(kh))*reliso(kiso) + wgtmola
17571  continue
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

       kiso = 1
       reliso(1) = 1.d0
       iso = 1
       wgtmola = 24.99287 
       write(18,*) ' mean molecular weight for MgH is',wgtmola
       write(18,*) ' (only 24Mg1H is available from input) '
       write(18,*) ' (so recall that we have FORCED ISO=1)'


       open (unit=65,status='old',readonly,
     &   file='/ste1/uffegj/mgh/linelist.dat')
       open (unit=67,status='old',readonly,
     &   file='/ste1/uffegj/mgh/MgH.partfunc')
C
       read(67,110) adum
       read(67,110) adum
       read(67,110) adum
       kqt = 0
       do 17590 k=1,1000
       read(67,*,end=1759) tin(k),qvin(k), qharm
       kqt = kqt + 1
       if (kqt.gt.nqt) stop ' kqt > nqt'
17590   continue
1759   continue
       write(6,*)'I read',kqt,' partition function values for MgH'
       write(6,*) 'from t=',tin(1),' to t=',tin(kqt)


       nmgh=0
       ngfz=0  !#lines with gf = 0
       write(6,*) ' First 5 lines of line strength data:'
       write(6,*) ' v0 gf exs ivl ivu jl ju iso: '
       do 1752 i=1,10000000
       read(65,*,end=1753) v0, gf, exs, ivl, ivu, jl, ju
        if (i.le.5) write(6,1755) v0,gf,exs,ivl,ivu,jl,ju
        nmgh = nmgh + 1
        gfiso = gf * reliso(iso)
        if (gfiso .eq. 0.D0) then
           ngfz = ngfz + 1
           if(ngfz.le.5)write(6,17551)v0,gf,exs,ivl,ivu,jl,ju
           go to 1752
        end if
        write(66,116)v0,exs,gfiso
1752    continue
       write(6,*)' Stop: more than 1.e8 lines in MgH input file ??'
c      go to 999
1753    continue
       write(6,*) ' Last line of line strength data:'
       write(6,1755) v0,gf,exs,ivl,ivu,jl,ju
1755   format(f11.3,1pe12.3,0pf11.2,2i4,2i4)
17551  format(f11.3,1pe12.3,0pf11.2,2i4,2i4,' (gfiso=0)')
1751   format(' read total',i6,'; wrote total',i6,' MgH lines.')
       write(6,1751) nmgh,nmgh-ngfz

C                                                                     ^
C                                         MgH    MgH     MgH    MgH   |
C----------------------------------------------------------------------
      GO TO 990
916   CONTINUE
 


990   CONTINUE


C
      DO 302 IT=1,KQT
      QTLN(IT)=LOG(qvin(it))
302   continue
      CALL TB04AT(NQT,KQT,tin,QTLN,DQDT,W)
C
      II = -1
      write(6,312) molid
312   format( ' it temp Q(it) for molecule: ',a4)
      DO 303 IT=1,KTEMP
      TINT=TMOL(IT)
      qvib(it) = TG01BT(II,NQT,KQT,tin,QTLN,DQDT,TINT)

      qvib(it) = EXP(qvib(it))
303   CONTINUE
      write(6,308) (tmol(it),it=1,ktemp)
      write(6,309) (qvib(it),it=1,ktemp)
308   format (' T: ',12f9.1)
309   format ('Qv: ',1p12e9.2)

C         DO 411 kt=1,kqt
C411      DADT(kt) = DADT(kt)*Qlog(kt)    !dadt = d(qvr)/d(tqin)

      go to 998
999   CONTINUE
      stop
998   CONTINUE

      return 
      end

C *******************************************************************
C
      SUBROUTINE RANGE (WS,NW)
C
C This subroutine identifies the nearest OS wavenumber to a given frequency
C inside the OS-wavenumber interval; i.e. WNOS(1) < WS < WNOS(NWNOS)
C
C ********************************************************************
C
C      PARAMETER (NWL=149000,NTEMP=12)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'atomparameter.inc'
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
C

      n1 = 1
      n2 = nwnos


      do 100 i=1,nwnos

      if (n2-n1.eq.1) then
          nw = n1
          go to 900
      end if

      khalf = (n2 + n1 ) / 2
      if (ws.lt.wnos(khalf)) then
         n2 = khalf
      else
         n1 = khalf
      end if

100   continue

900   continue
      return 
      end
C
C
      FUNCTION TG01BT(II,ND,N,X,F,D,XX) 
      implicit real*8 (a-h,o-z)
      DIMENSION X(ND),F(ND),D(ND)
      COMMON /TG01BA/I1,IN,KK
C     I1 = 2     ! we fix for linear extrapolation (possible in log)
      I1 = 1
      IN = 1
C i1/in=0~tg01bt=0, =1~set = last value, =2~extrapolate linearly, =3~2.order
C ROUTINE TO CALCULATE VALUE FXX OF SPLINE IN POINT XX WHEN N KNOTS
C (XI,FI) WITH DERIVATIVE DI ARE GIVEN.
C II<0 => SEARCH THE WHOLE RANGE. II >= 0 => FUNCTION HAS PREVIOUSLY 
C BEEN ENTERED WITH A SMALLER VALUE OF XX.
C COMMON VALUES I1 AND IN CONTROLS WHAT TO DO IF XX IS OUTSIDE X INTERVAL.
C
C II NEGATIVE, RESET
      IF(II.LT.0) KK=2  
C   
C CHECK IF OUTSIDE  
      IF(XX.LT.X(1)) GOTO 110   
      IF(XX.GT.X(N)) GOTO 120   
      DO 100 K=KK,N 
      IF(XX.LT.X(K)) GOTO 101   
100   CONTINUE  
      KK=N  
      GOTO 102  
101   KK=K  
C   
C CALCULATE FUNCTION
102   DX=X(KK)-X(KK-1)  
      DF=F(KK)-F(KK-1)  
      P=(XX-X(KK-1))/DX 
      Q=1.-P
      TG01BT=Q*F(KK-1)+P*F(KK)+P*Q*  
     & (Q*(D(KK-1)*DX-DF)-P*(D(KK)*DX-DF))  
      RETURN
C   
C BEFORE X(1)   
110   TG01BT=0.  
      IF(I1.LE.0) RETURN
      TG01BT=F(1)
      IF(I1.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(1))*D(1)
      IF(I1.EQ.2) RETURN
      DX=X(2)-X(1)  
      D2=2.*(3.*(F(2)-F(1))/DX**2-(2.*D(1)+D(2))/DX)
      TG01BT=TG01BT+.5*(XX-X(1))**2*D2
      IF(I1.EQ.3) RETURN
      D36=(D(1)+D(2)-2.*(F(2)-F(1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(1))*(XX-X(1))**2*D36
      RETURN
C   
C AFTER X(N)
120   TG01BT=0.  
      IF(IN.LE.0) RETURN
      TG01BT=F(N)
      IF(IN.EQ.1) RETURN
      TG01BT=TG01BT+(XX-X(N))*D(N)
      IF(IN.EQ.2) RETURN
      DX=X(N)-X(N-1)
      D2=2.*(-3.*(F(N)-F(N-1))/DX**2+(2.*D(N)+D(N-1))/DX)   
      TG01BT=TG01BT+.5*(XX-X(N))**2*D2
      IF(IN.EQ.3) RETURN
      D36=(D(N)+D(N-1)-2.*(F(N)-F(N-1))/DX)/DX**2   
      TG01BT=TG01BT+(XX-X(N))*(XX-X(N))**2*D36
      END
C
C
C
      SUBROUTINE TB04AT(ND,N,X,F,D,W)
      implicit real*8 (a-h,o-z)
C
      DIMENSION X(ND),F(ND),D(ND),W(3,ND)   
C   
C THIS VERSION OF TB04A IS identical to TB04A, except that it can be
C called with arrays which might be dimensioned larger than the part
C of it which is used for the interpolation.
C INPUT ARE X(I), FUNCTION VALUE IN N KNOTS. THEN THE DERIVATIVE IS 
C CALCULATED IN THE KNOTS. W=0 AT SUCCESFULL RETURN, OTHERWISE W=1.
C X(I) SHOULD BE IN STRICTH INCREASING ORDER, X1<X2<...<XN.
C In the call from SUBROUTINE PROFILE, X(1), X(2),...,X(NTAU) are the 
C temperature values at the NTAU optical depth values, and the derivative
C dPg/dT (=D(N))is calculated in the NTAU points, for later use to 
C calculate Pg in the points where the temperature of the
C absorption coefficient is known.  
C         
C FIRST POINT   
      CXB=1./(X(2)-X(1))
      CXC=1./(X(3)-X(2))
      DFB=F(2)-F(1) 
      DFC=F(3)-F(2) 
      W(1,1)=CXB*CXB
      W(3,1)=-CXC*CXC   
      W(2,1)=W(1,1)+W(3,1)  
      D(1)=2.*(DFB*CXB*CXB*CXB-DFC*CXC*CXC*CXC) 
C   
C INTERIOR POINTS   
      N1=N-1
      DO 100 K=2,N1 
      CXA=CXB   
      CXB=1./(X(K+1)-X(K))  
      DFA=DFB   
      DFB=F(K+1)-F(K)   
      W(1,K)=CXA
      W(3,K)=CXB
      W(2,K)=2.*(CXA+CXB)   
      D(K)=3.*(DFB*CXB*CXB+DFA*CXA*CXA) 
100   CONTINUE  
C   
C LAST POINT
      W(1,N)=CXA*CXA
      W(3,N)=-CXB*CXB   
      W(2,N)=W(1,N)+W(3,N)  
      D(N)=2.*(DFA*CXA*CXA*CXA-DFB*CXB*CXB*CXB) 
C   
C ELIMINATE AT FIRST POINT  
      C=-W(3,1)/W(3,2)  
      W(1,1)=W(1,1)+C*W(1,2)
      W(2,1)=W(2,1)+C*W(2,2)
      D(1)=D(1)+C*D(2)  
      W(3,1)=W(2,1) 
      W(2,1)=W(1,1) 
C   
C ELIMINATE AT LAST POINT   
      C=-W(1,N)/W(1,N-1)
      W(2,N)=W(2,N)+C*W(2,N-1)  
      W(3,N)=W(3,N)+C*W(3,N-1)  
      D(N)=D(N)+C*D(N-1)
      W(1,N)=W(2,N) 
      W(2,N)=W(3,N) 
C   
C ELIMINATE SUBDIAGONAL 
      DO 110 K=2,N  
      C=-W(1,K)/W(2,K-1)
      W(2,K)=W(2,K)+C*W(3,K-1)  
      D(K)=D(K)+C*D(K-1)
110   CONTINUE  
C   
C BACKSUBSTITUTE
      D(N)=D(N)/W(2,N)  
      DO 120 KK=2,N 
      K=(N+1)-KK
      D(K)=(D(K)-W(3,K)*D(K+1))/W(2,K)  
120   CONTINUE  
C   
      RETURN
      END
C
      subroutine sw_h2o
C      PARAMETER (NWL=149000,NTEMP=12)
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'
c
c   UGJ990615:  Copy of HSs linexud.f for own format use etc. 
c   this subroutine is meant to work with the OS.f program.
c
c     read files containing ro-vibrational levels
c     and line strengths and produce line list
c     of approximately hitran type format (written to unit 8).
c     for h2o. (proton nuclear spin statistics)
c     unit 1 contains a list of files containing the ro-vibrational
c     levels
c     unit 2 contains a list of files containing the packed line data
c
c     idrv is number of ro-vibrational levels
c     idj is number of j blocks
c     p,q, or r branch.
c
      character*80 fname                                                8d22s96
      character*1 a(80)                                                 8d22s96
      parameter (idrv=180000,idj=60,igrd=8000)
      parameter (nx=38,ncx=1633,ibase=229)
      dimension erov(idrv),iqnrv(5,idrv),iadd(idj,2,2),
     $          pop(idrv),numbb(idj,2,2),
     $          idig(4),ipack(2),idig1(8,10),                           8d22s96
     $          idig2(80)                                               8d22s96
      equivalence (pack,ipack)                                          8d22s96
      equivalence (idig1,idig2)                                         8d23s96
      dimension spec(igrd,2),popj(idj)
      dimension nle(0:55), xinle(0:55)

      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      COMMON /COPAC/HALF,sig_nu(ntemp),hckt(ntemp)
     *    ,sqr_pi_m,sig_sqr,sumop_c(ntemp),a17,n16,n17,npr,j
      COMMON/CSWMX/sgf,sumop_l(ntemp)
C
       data temp,thres/ 1585.d0, 1.d-29/
       tocm=1d0/4.556335d-6                                              5d23s91
       xk=3.16681d-6
       conv=tocm*16.1941d5/6.0221367d23
       todebye=2.5417662d0                                               9d16s94
       todebye2=todebye**2
c
c     temp is temperature in kelvin to evaluate partion function and
c          line strengths
c     thres is threshold (in cm/molecule) for printing out lines.
c           use thres lt zero to just compute partition function.
c     xlow - lower frequcency to consider
c     xhigh - highest frequency to consider
c
c     temp=296d0
c     thres=1d-30
      xlow=0d0
      xhigh=28279d0
      write(6,*)('frequency limits '),xlow,xhigh
      tokm=6.0221367d18
      dx=(xhigh-xlow)/dfloat(igrd-1)
      do 121 i=1,igrd
       spec(i,1)=0d0
       spec(i,2)=0d0
  121 continue
      do 2105 k=0,55
      nle(k) = 0
      xinle(k) = 0.d0
2105  continue
      write(6,*)('temperature = '),temp
      beta=1d0/(xk*temp)
      open(unit=101,file='/ste1/uffegj/sw/H2O/fortE.1'
     *  ,status='old',readonly)      !file with names of energy level files
      open(unit=102,file='/ste1/uffegj/sw/H2O/fortL.2'
     *  ,status='old',readonly)      !file with names of line list files
c
c     data for packing/formatting stuff
c
      write(6,*)('data for packing ')
      write(6,*)('ibase = '),ibase
      write(6,*)('nx = '),nx
      write(6,*)('ncx = '),ncx
      binv=1d0/dfloat(ibase)
      base=dfloat(ibase)
      ibase2=ibase**2
      ibase3=ibase2*ibase
      xnx=dfloat(nx)
      xmat=base/xnx
      write(6,*)('xmat '),xmat
      ixmat=int(xmat)
      xmat=dfloat(ixmat)
      write(6,*)('xmat '),xmat
      xmati=1d0/xmat
      trunc=5d0*base/(xmat*(base**4))
      write(6,*)('truncation error = '),trunc
      rndc=0.5d0*(binv**3)*xmati
      write(6,*)('rndc '),rndc
      ilab4=ncx*ncx
      ilab3=ilab4*60
      ilab2=ilab3*2
      ilab1=ilab2*2
c
c     read in ro-vibrational states and compute partition function
c
      do 1 i=1,idj
       iadd(i,1,1)=0
       iadd(i,2,1)=0
       iadd(i,1,2)=0
       iadd(i,2,2)=0
       numbb(i,1,1)=0                                                   8d22s96
       numbb(i,2,1)=0                                                   8d22s96
       numbb(i,1,2)=0                                                   8d22s96
       numbb(i,2,2)=0                                                   8d22s96
    1 continue
      do 14 i=1,idj
       popj(i)=0d0
   14 continue
      nrv=1
      partf=0d0
      jmax=-1
      jmin=idj                                                          8d22s96
      ehigh=0d0
  200 continue                                                          8d22s96
      read(101,*,end=201)fname                                            8d22s96
      write(6,*)('opening file '),fname
      open(unit=103,file=fname,form='formatted',status='old')             8d22s96
    2 continue
      read(103,3,end=4)jtot,ipar,isym,numb                                8d22s96
    3 format(4i5)
      jmax=max(jmax,jtot)
      jmin=min(jmin,jtot)                                               8d22s96
      if(jtot+1.gt.idj)then
       write(6,*)('too high a jtot. increase idj to '),jtot+1
       stop
      end if
      ipar=mod(ipar,2)+1
      if(isym.eq.0)isym=2
Cugj  write(6,7)jtot,ipar,isym,numb
    7 format('jtot = ',i5,' ipar = ',i1,' isym= ',i1,' numb = ',i4)
      iadd(jtot+1,ipar,isym)=nrv
      numbb(jtot+1,ipar,isym)=numb
      if(nrv+numb-1.gt.idrv)then
       write(6,*)('too many ro-vib levels: increase idrv to '),
     $      nrv+numb-1
       stop
      end if
      xj=dfloat(2*jtot+1)
      if(isym.eq.1)xj=xj*3d0
      sumj=0d0
      do 5 i=1,numb
       read(103,6)erov(nrv),(iqnrv(j,nrv),j=1,5)                          8d22s96
       ehigh=max(ehigh,erov(nrv))
       if(partf.eq.0d0)then
        emin=erov(nrv)
       else
        emin=min(emin,erov(nrv))
       end if
       pop(nrv)=xj*exp(-beta*erov(nrv))
       popj(jtot+1)=popj(jtot+1)+pop(nrv)
       sumj=sumj+pop(nrv)
       partf=partf+pop(nrv)
    6  format(1pe21.14,5i4)
       nrv=nrv+1
    5 continue
      sumj=1d0/sumj
      tsumj=0d0
      plast=pop(nrv-1)*sumj
      do 100 i=1,numb
       tsumj=tsumj+pop(nrv-i)
       if(tsumj*sumj.gt.0.01d0)then
Cugj    write(6,103)numb+1-i,numb,plast
  103   format(1x,'99% of pop at root ',i4,' out of ',i4,' roots.',     8d22s96
     $       ' last level has fraction ',1pe10.2)
        go to 2
       end if
  100 continue
      go to 2
    4 continue
      close(unit=103)                                                     8d22s96
      go to 200                                                         8d22s96
  201 continue                                                          8d22s96
      nrv=nrv-1
      write(6,*)('total no. ro-vibrational levels = '),nrv,
     $     ('   below [au, cm-1]'),ehigh,ehigh*tocm
      write(6,*)('ground state energy [au; cm-1] = '),emin,emin*tocm
      write(6,*)('partition function = '),partf
      partf0=partf*exp(beta*emin)
      write(6,*)('partition function with energy zero at ground'),
     $      (' state energy = '),partf0
      xnorm=1d0/partf
      do 13 i=1,nrv
       pop(i)=pop(i)*xnorm
       erov(i)=erov(i)-emin
   13 continue
      do 15 j=jmin+1,jmax+1                                             8d22s96
       popj(j)=popj(j)*xnorm
       if(mod(j-1,2).eq.0)then                                          8d22s96
        i1a=1                                                           8d22s96
        i1b=2                                                           8d22s96
       else                                                             8d22s96
        i1a=2                                                           8d22s96
        i1b=1                                                           8d22s96
       end if                                                           8d22s96
Cugj   write(6,16)j-1,popj(j),(numbb(j,i1a,i2),numbb(j,i1b,i2),i2=1,2)
   16  format(1x,'fraction in j = ',i3,' is ',f8.6,4i5)
   15 continue
      if(thres.lt.0d0)stop
c
c     read line list
c
      write(6,*)('line intensity threshold at temp = '),thres
      iso=11
      nlinet=0
      nline=0
      nall=0
      xsmall=1d10
      ixx=0
c$$$      write(6,*)('jmax = '),jmax
  300 continue
      read(102,*,end=301)fname
      write(6,*)('reading from '),fname
      open (unit=65,status='old',readonly,file=fname,form='formatted')
      ntotf=0
      omlow=1d10
      omhigh=0d0
    8 continue
       read(65,*,end=9)numb
       ndo=numb/10
       if(ndo*10.ne.numb)ndo=ndo+1
       iread=0
       do 11 ido=1,ndo
        nread=min(10,numb-iread)
        iread=iread+nread
        ird=nread*8
        read(65,402)(a(i),i=1,ird)
  402   format(80a1)
        do 450 i=1,ird
         idig2(i)=ichar(a(i))-27
  450   continue
        do 451 i=1,nread
         ilab=idig1(4,i)+ibase*(idig1(3,i)+ibase*(idig1(2,i)
     $        +ibase*idig1(1,i)))
         jdiff=ilab/ilab1
         ilab=ilab-ilab1*jdiff
         jdiff=jdiff-1
         is=ilab/ilab2
         ilab=ilab-is*ilab2
         is=is+1
         ipi=ilab/ilab3
         ilab=ilab-ilab3*ipi
         ipi=ipi+1
         ji=ilab/ilab4
         ilab=ilab-ji*ilab4
         iqn2=ilab/ncx
         ilab=ilab-iqn2*ncx
         iqn2=iqn2+1
         iqn1=ilab+1
         ipf=3-ipi
         jf=ji+jdiff
         iai=iadd(ji+1,ipi,is)+iqn1-1
         iaf=iadd(jf+1,ipf,is)+iqn2-1
         ei=erov(iai)
         ef=erov(iaf)
         om=ef-ei
         omau=om
         om=om*tocm
         imat=idig1(8,i)/nx
         ix=idig1(8,i)-imat*nx
         xsect=dfloat(imat)*xmati+dfloat(idig1(7,i))
         xsect=xsect*binv+dfloat(idig1(6,i))
         xsect=xsect*binv+dfloat(idig1(5,i))
         xsect=xsect*binv*(1d-1**ix)
         xi=conv*xsect*pop(iai)*omau*(1d0-exp(-omau*beta))
Cugj     gf=1.1295842d12*conv*xsect*omau   == N_A/10^5 * 1.87572d-7 * 16.1941 * xsect * om[cm^-1]
         gf= 3.03756d-6 * xsect * om * (2.d0*ji + 1.d0)
         icent=int((om-xlow)/dx)+1
         if(icent.ge.1.and.icent.le.igrd)spec(icent,1)=spec(icent,1)+xi
         do 2100 k=0,55
         if (xi.gt.1.d0*10.**-k) then
            nle(k) = nle(k) + 1
            xinle(k) = xinle(k) + xi
            go to 2101
         end if
2100     continue
2101     continue
          ef=ef*tocm
          ei=ei*tocm

         if (is.eq.1) gs = 3.d0
         if (is.eq.2) gs = 1.d0

         gf = gs * (2.d0*ji + 1.d0) * 3.03756d-6 * xsect * om
         xsect=xsect*todebye2
         gfnont = gs * (2.d0*ji+1.d0) * xinont * 1.1295842d12 !==1.87572d-12 * 6.0221367d23
Cugj  gf=1.1295842d12*conv*xsect*omau == N_A/10^5 * 1.87572d-7 * 16.1941 * xsect * om[cm-1]

         if(xi.gt.thres.and.om.gt.xlow.and.om.lt.xhigh)then
c
c     iso=11 meaning 1H1H16O
c     om - transition frquence in cm**-1
c     xi - line strength in cm/molecule (for absorbtion)
c     xsect - square of dipole matrix element in Debye**2
c     ei - initial state energy in cm**-1
c     jf - final state total angular momentum quantum number
c     iqnrv(4,iaf) - Kp for final state
c     iqnrv(5,iaf) - Ko for final state
c     ji - initial state total angular momentum quantum number
c     iqnrv(4,iai) - Kp for initialstate
c     iqnrv(5,iai) - Ko for initial state
c     iqnrv(1,iaf) - symmetric stretch quantum number for final state
c     iqnrv(2,iaf) - bending quantum number for final state
c     iqnrv(3,iaf) - asymmetric stretch quantum number for final state
c     iqnrv(1,iai) - symmetric stretch quantum number for iitial state
c     iqnrv(2,iai) - bending quantum number for iitial state
c     iqnrv(3,iai) - asymmetric stretch quantum number for iitial state
c
C      if(xi .gt. 1.e-20) then
C         write(8,116) om, ei, gf
Cugj      write(8,12)iso,om,xi,xsect,ei,jf,iqnrv(4,iaf),iqnrv(5,iaf),
Cugj $         ji,iqnrv(4,iai),iqnrv(5,iai),iqnrv(1,iaf),iqnrv(2,iaf),
Cugj $         iqnrv(3,iaf),iqnrv(1,iai),iqnrv(2,iai),iqnrv(3,iai),
Cugj $         iqn2,iqn1,ipi,is
C      end if
116    format(2f11.4,1pe13.5)
   12     format(i3,f12.5,1p2e10.3,10x,0pf11.5,17x,3i2,3x,3i2,3x,3i2,3x,
     $           3i2,3x,2i4,2i2)



          nline=nline+1
C Lines here are in the interval xlow to xhigh and have intensity above thres at
C temperature temp. These lines are added to the OS.

      v0 = om
      exs = ei
      absl = sgf * gf       !cm/mol
      do 139 kt = 1,ktemp
         ST = absl
     &     *EXP(-exs*HCKT(kt))*(1.-EXP(-v0*HCKT(kt)))/qvib(kt)
         sumop_l(kt) = sumop_l(kt) + st                         !integrated abs.coef in [cm/mol]
139   continue

      if (v0.lt.wnos(1)) go to 451
C     if (v0.gt.wnos(nwnos)) go to 119
C temp we go to 110 because h2o list is not ordered
      if (v0.gt.wnos(nwnos)) go to 451

      call range(v0,nw)   !this subr. identifies nw such that
C                          wnos(nw) < v0 < wnos(nw+1)

      if (v0.lt.wnos(nw) .or. v0.gt.wnos(nw+1) ) then
        write(18,*) ' nw,nw+1,v0 ???  => stop! '
        write(18,138) nw,v0,wnos(nw),wnos(nw+1)
C       go to 900
        stop
138     format (' i,nr,nw,v0,wnos{nw,nw+1}=',i5,2i5,3f11.3)
      else
        nr = nline
        if(nr.le.5)write(18,138) i,nr,nw,v0,wnos(nw),wnos(nw+1)
      end if

         CALL OPAC_COM(v0,gf,exs,nw)


         end if     !endif(xi.gt.thres.and.om.gt.xlow.and.om.lt.xhigh)then

  451   continue
   11  continue
       go to 8
    9 continue
      close(unit=65)
      go to 300
  301 continue
      write(6,*)('total number of lines above threshold = '),nline
      write(6,*) ' #lines with xi in the interval 1.e-k<xi<=1.e-(k-1)'
      write(6,*) ' k  nlines(k)  nlines_accumulated',
     * ' intensity  accumul.int.  accum.rel.int.'
      naccum = 0
      do 2112 k=55,0,-1
      naccum = naccum + nle(k)
      nlelast = k
      if(naccum .gt. 0) go to 2113
2112  continue
2113  continue
      xisum = 0.d0
      do 2114 k=0,nlelast
2114  xisum = xisum + xinle(k)

      naccum = 0
      xiaccum = 0.d0
      do 2110 k=0,nlelast
      naccum = naccum + nle(k)
      if (naccum .eq. 0) go to 2110
      xiaccum = xiaccum + xinle(k)
      write(6,2111) 
     *k, nle(k), naccum, xinle(k), xiaccum, xiaccum/xisum
2110  continue
2111  format(i3,2i12, 1p2e12.3,0pf12.9)
      return
      end




       subroutine vald_open
C       parameter (ntemp=12,mt1=ntemp+1,nqt=200)
C       parameter (nqt=200)
       implicit real*8 (a-h,o-z)
       include 'atomparameter.inc'

C this program reads the VALD line list and prepares it for 
C a high resolution opacity sampling for spectrum synthesis
C and other applications. It is intended to be a subroutine 
C in the os.f program too.


C      dimension qvin(nqt)
       dimension qvin(ntemp)
       character elm_pick*3,dum*62
       common /cvald/ wb, we, lvald, ion_pick, elm_pick
C      common /cnumatom/ watom(92),namatom(92)
C      character*2 namatom
      common /cwatom/watom(200)
      common /catomname/atomname(200)
      character*2 atomname
      common /cgemat/kgemat(200),kgemion1(200),kgemionm(200)
      common /cprespp/prespp(ndp,nspec)           !partial pressures in dyn/cm^2
      COMMON/CISO/vkmssp,wgtmola,rel_iso_c(3),pc3(3)
     & ,reliso(5),kiso,nhalf,jderiv,kiso_c
     & ,linesc3_direct,lineshcn_direct,linesc2h2_direct
       character elm*2, nameatom*2
       common /cvaldline/ wgtatom, vv0, vexs, vgf,
     &    gam_rad, gam_stark, gam_waal, fac_lande, 
     &    conv_gstel(ntemp),
     &    jatom, jg, ion, lgstel,
     &    elm, nameatom
      common/cvoigt/ pgas(mt1),p_el(mt1),p_h(mt1),hwda(ntemp)
     &              ,pemu(mt1),pp(mt1,nspec)
     &              ,profile,atmos_voigt
      COMMON/COS/wnos(nwl),opacos(ntemp,nwl),wlos(nwl),wlosstep(nwl)
     *   ,tmol(ntemp),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos


c      open(unit=1,file='/ste3/uffegj/vald/atm_lines_cool_obs.dat',
c    &           status='old',readonly)
       open(unit=1,file='/p7/hpg688/hpg688/VALD/vald3.dat',
     &           status='old',readonly)

C if we are to compute a cm2/mol OS for one atom only, then
C compute partition function qvib(it) just once:

             IF(elm_pick.ne.'999') THEN
             do 30 j=1,92
              if(elm_pick(1:2).eq.atomname(j)) then
                 wgtatom = watom(j)
                 wgtmola = wgtatom
                 ion = ion_pick
                 jatom = j
                 call partf(jatom,ion,tmol,ktemp,qvin)
                 if(ion.eq.1) then
                  jg = kgemat(j)
                 else if(ion.eq.2) then
                  jg = kgemion1(j)     !prespp index of atom j (neutral or +)
                 end if
                 write(6,316) j,jg,ion,atomname(j) 
316              format('ID: ',3i4,2x,a2)
                 if (jg.eq.1) goto 33
 
                 do 310 it = 1,ktemp
                 conv_gstel(it) = PP(it,jg)/( PGAS(it)*PEMU(it) )  ![mol_pp/mol*]/[g*/mol*]
310              qvib(it) = qvin(it)

                      write(6,308) (tmol(it),it=1,ktemp)
                      write(6,309) (qvib(it),it=1,ktemp)
                      write(6,315) (pp(it,jg),it=1,ktemp)
                      write(6,312) (pgas(it),it=1,ktemp)
                      write(6,313) (pemu(it),it=1,ktemp)
                      write(6,314) (conv_gstel(it),it=1,ktemp)
308                   format (' T: ',12f8.1)
309                   format ('Qv: ',12f8.2)
315                   format ('Pp: ',1p12e12.3)
312                   format ('Pg: ',1p12e12.3)
313                   format ('emu: ',12f8.2)
314                   format ('conv_gstel: ',1p12e12.3)
                 go to 31
              end if
30           continue

             write(6,32) elm_pick
32           format(' I couldnt identify atom ',a3)
             stop
33           continue
             write(6,34) elm_pick
34           format(' Element ',a3, ' doesnt exist in GEM calc')
             stop

31           continue     !we identified the atom, and PP can be computed
             END IF

C If on the other hand the computation is a cm2/g* OS for a
C combination of several atoms and/or ions, we need to know the 
C partial pressures of all the atoms/ions. We compute them here
C based on the T,Pg structure and call to the gem_routines.

C      call tpgread
C      call model


       read (1,11) dum
       write(6,11) dum
       read (1,11) dum
       write(6,11) dum
11     format(a70)
       ipick = 0
       kentry = 0
       nline_nogem = 0
       nvneutral = 0
       nvion1= 0
       nvion2 = 0
       return

C      do 100 i=1,1000000

       entry vald_line(kadopt)

       kentry = kentry + 1
       kadopt = 0
       read(1,12,end=998,err=995) elm, ion, wl, exc, gflog,
     &            gam_rad, gam_stark, gam_waal, fac_lande
12     format(1x,a2,i2,2x,f10.4,1x,f8.4,1x,f7.3,4(1x,f6.3),1x)
13     format(1x,a2,i2,2x,f11.3,1x,f8.4,1x,f7.3,4(1x,f6.3),1x)
       if(elm.eq.elm_pick .or. elm_pick.eq.'999') then
        if(ion.eq.ion_pick .or. ion_pick.eq.999)  then
         if(wl.ge.wb .or. wb.eq.999.) then
         if(wl.le.we .or. we.eq.999.) then
             if(elm_pick.eq.'999' .or. ion_pick.eq.999) then
             do 20 j=1,92
              if(elm.eq.atomname(j)) then
                 if(ion.le.1) nvneutral = nvneutral + 1
                 if(ion.eq.2) nvion1= nvion1 + 1
                 if(ion.gt.2) nvion2 = nvion2 + 1
                 wgtatom = watom(j)
                 nameatom = atomname(j)
                 jatom = j
                 call partf(jatom,ion,tmol,ktemp,qvin)
                 if(ion.eq.1) then
                  jg = kgemat(j)
                 else if(ion.eq.2) then
                  jg = kgemion1(j)     !prespp index of atom j (neutral or +)
                 end if
                 if(jg.eq.1) then                  !the element/ion is not in the GEM comp.
                   nline_nogem = nline_nogem + 1   !skip this entry in spectrum (i.e.,kadopt=0)
                   go to 999
                 end if
                 do 311 it = 1,ktemp
                 conv_gstel(it) = PP(it,jg)/( PGAS(it)*PEMU(it) )  ![mol_pp/mol*]/[g*/mol*]
311              qvib(it) = qvin(it)
                 go to 21
              end if
20           continue
             write(6,22) elm,i
22           format(' I couldnt identify atom ',a2,' at input line',i4)
             stop
21           continue
             end if
          vgf = 10.**gflog
          vexs= 8067.1 * exc      !conversion from eV to cm^-1
          vv0 = 1.d8 / wl         !conv lambda AA --> wavenumber cm^-1
C           write(66,116)vv0,vexs,vgf,wgtatom
C 116       format(2f11.4,1pe13.5,0pf9.4)
          kadopt = 1
          ipick = ipick + 1
C         if(ipick.le.2. .or. elm.eq.'Cs' .or. elm.eq.'Eu' 
C    &     .or.exc.eq.0.)
          if(ipick.le.5 )
     &    write(6,13) elm, ion, wl, exc, gflog,
     &            gam_rad, gam_stark, gam_waal, fac_lande
         end if
         end if
        end if
       end if
C      ic = i
C 100    continue
       go to 999

995    continue
       backspace(1)
       do 120 i=1,100
       read(1,11,end=999) dum
        if(i.eq.1) then
            if(dum(1:11).eq.' References') go to 998
            write(6,11) dum
            write(6,*) ' Possibly reading problems:'
            write(6,11) dum
        end if  
        write(6,11) dum
120    continue

998    continue

       write(6,*) ' we successfully read ',kentry,' lines from VALD'
       write(6,*) ' we found ',ipick,' lines within specifications'
       write(6,*) nline_nogem,' lines in vald had no gem data (skiped)'
       write(6,*) ' there were ',nvneutral,' neutral lines'
       write(6,*) ' there were ',nvion1,' one time ionized lines'
       write(6,*) ' there were ',nvion2,' two time ionized lines'
       kadopt = 999


999    continue

       return
       end

C____________________________________________________________________

      SUBROUTINE atomw_old
C
C Mean atomic weights are stored into the vector watom, AB2001
C Species 1:92 are included, as in irwin.dat
C

      IMPLICIT real*8 (a-h,o-z)

      common /cnumatom/ watom(92),namatom(92)
      character*2 namatom

      DO i=1,92
         watom(i) = 0.0
      ENDDO

      watom(1)  =   1.0079 !H  Hydrogen
      watom(2)  =   4.0026 !He Helium
      watom(3)  =   6.941  !Li Lithium
      watom(4)  =   9.0121 !Be Beryllium
      watom(5)  =  10.81   !B  Boron
      watom(6)  =  12.011  !C  Carbon
      watom(7)  =  14.0067 !N  Nitrogen
      watom(8)  =  15.9994 !O  Oxygen
      watom(9)  =  18.9984 !F  Fluorine
      watom(10) =  20.179  !Ne Neon
      watom(11) =  22.9897 !Na Sodium
      watom(12) =  24.305  !Mg Magnesium
      watom(13) =  26.9814 !Al Aluminum
      watom(14) =  28.0855 !Si Silicon
      watom(15) =  30.9737 !P  Phosphorus
      watom(16) =  32.06   !S  Sulfur
      watom(17) =  35.453  !Cl Chlorine
      watom(18) =  39.948  !Ar Argon
      watom(19) =  39.0983 !K  Potassium
      watom(20) =  40.08   !Ca Calcium
      watom(21) =  44.9559 !Sc Scandium
      watom(22) =  47.88   !Ti Titanium
      watom(23) =  50.9415 !V  Vanadium
      watom(24) =  51.996  !Cr Chromium
      watom(25) =  54.9380 !Mn Manganese
      watom(26) =  55.847  !Fe Iron
      watom(27) =  58.9332 !Co Cobalt
      watom(28) =  58.96   !Ni Nickel
      watom(29) =  63.546  !Cu Copper
      watom(30) =  65.38   !Zn Zinc
      watom(31) =  69.72   !Ga Gallium
      watom(32) =  72.59   !Ge Germanium
      watom(33) =  74.9216 !As Arsenic
      watom(34) =  78.96   !Se Selenium
      watom(35) =  79.904  !Br Bromine
      watom(36) =  83.80   !Kr Krypton
      watom(37) =  85.4678 !Rb Rubidium
      watom(38) =  87.62   !Sr Strontium
      watom(39) =  88.9059 !Y  Yttrium
      watom(40) =  91.22   !Zr Zirconium
      watom(41) =  92.9064 !Nb Niobium
      watom(42) =  95.94   !Mo Molybdenum
      watom(43) =  97.907  !Tc Technetium
      watom(44) = 101.07   !Ru Ruthenium
      watom(45) = 102.9055 !Rh Rhodium
      watom(46) = 106.42   !Pd Palladium
      watom(47) = 107.868  !Ag Silver
      watom(48) = 112.41   !Cd Cadmium
      watom(49) = 114.82   !In Indium
      watom(50) = 118.69   !Sn Tin
      watom(51) = 121.75   !Sb Antimony
      watom(52) = 127.60   !Te Tellurium
      watom(53) = 126.9045 !I  Iodine
      watom(54) = 131.29   !Xe Xenon
      watom(55) = 132.9054 !Cs Cesium
      watom(56) = 137.33   !Ba Barium
      watom(57) = 138.9055 !La Lanthanum
      watom(58) = 140.12   !Ce Cerium
      watom(59) = 140.9077 !Pr Praseodymium
      watom(60) = 144.24   !Nd Neodymium
      watom(61) = 144.913  !Pm Promethium
      watom(62) = 150.36   !Sm Samarium
      watom(63) = 151.96   !Eu Europium
      watom(64) = 157.25   !Gd Gadolinium
      watom(65) = 158.9254 !Tb Terbium
      watom(66) = 162.50   !Dy Dysprosium
      watom(67) = 164.9304 !Ho Holmium
      watom(68) = 167.26   !Er Erbium
      watom(69) = 168.9342 !Tm Thulium
      watom(70) = 173.04   !Yb Ytterbium
      watom(71) = 174.967  !Lu Lutetium
      watom(72) = 178.49   !Hf Hafnium
      watom(73) = 180.9479 !Ta Tantalum
      watom(74) = 183.85   !W  Tungsten
      watom(75) = 186.207  !Re Rhenium
      watom(76) = 190.2    !Os Osmium
      watom(77) = 192.22   !Ir Iridium
      watom(78) = 195.08   !Pt Platinum
      watom(79) = 196.9665 !Au Gold
      watom(80) = 200.59   !Hg Mercury
      watom(81) = 204.383  !Tl Thallium
      watom(82) = 207.2    !Pb Lead
      watom(83) = 208.9804 !Bi Bismuth
      watom(84) = 208.982  !Po Polonium
      watom(85) = 209.987  !At Astatine
      watom(86) = 222.018  !Rn Radon
      watom(87) = 223.020  !Fr Francium
      watom(88) = 226.0254 !Ra Radium
      watom(89) = 227.0278 !Ac Actinium
      watom(90) = 232.0381 !Th Thorium
      watom(91) = 231.0359 !Pa Protactinium
      watom(92) = 238.051  !U  Uranium


      namatom(1)  = 'H '  ! Hydrogen
      namatom(2)  = 'He'  ! Helium
      namatom(3)  = 'Li'  ! Lithium
      namatom(4)  = 'Be'  ! Beryllium
      namatom(5)  = 'B '  ! Boron
      namatom(6)  = 'C '  ! Carbon
      namatom(7)  = 'N '  ! Nitrogen
      namatom(8)  = 'O '  ! Oxygen
      namatom(9)  = 'F '  ! Fluorine
      namatom(10) = 'Ne'  ! Neon
      namatom(11) = 'Na'  ! Sodium
      namatom(12) = 'Mg'  ! Magnesium
      namatom(13) = 'Al'  ! Aluminum
      namatom(14) = 'Si'  ! Silicon
      namatom(15) = 'P '  ! Phosphorus
      namatom(16) = 'S '  ! Sulfur
      namatom(17) = 'Cl'  ! Chlorine
      namatom(18) = 'Ar'  ! Argon
      namatom(19) = 'K '  ! Potassium
      namatom(20) = 'Ca'  ! Calcium
      namatom(21) = 'Sc'  ! Scandium
      namatom(22) = 'Ti'  ! Titanium
      namatom(23) = 'V '  ! Vanadium
      namatom(24) = 'Cr'  ! Chromium
      namatom(25) = 'Mn'  ! Manganese
      namatom(26) = 'Fe'  ! Iron
      namatom(27) = 'Co'  ! Cobalt
      namatom(28) = 'Ni'  ! Nickel
      namatom(29) = 'Cu'  ! Copper
      namatom(30) = 'Zn'  ! Zinc
      namatom(31) = 'Ga'  ! Gallium
      namatom(32) = 'Ge'  ! Germanium
      namatom(33) = 'As'  ! Arsenic
      namatom(34) = 'Se'  ! Selenium
      namatom(35) = 'Br'  ! Bromine
      namatom(36) = 'Kr'  ! Krypton
      namatom(37) = 'Rb'  ! Rubidium
      namatom(38) = 'Sr'  ! Strontium
      namatom(39) = 'Y '  ! Yttrium
      namatom(40) = 'Zr'  ! Zirconium
      namatom(41) = 'Nb'  ! Niobium
      namatom(42) = 'Mo'  ! Molybdenum
      namatom(43) = 'Tc'  ! Technetium
      namatom(44) = 'Ru'  ! Ruthenium
      namatom(45) = 'Rh'  ! Rhodium
      namatom(46) = 'Pd'  ! Palladium
      namatom(47) = 'Ag'  ! Silver
      namatom(48) = 'Cd'  ! Cadmium
      namatom(49) = 'In'  ! Indium
      namatom(50) = 'Sn'  ! Tin
      namatom(51) = 'Sb'  ! Antimony
      namatom(52) = 'Te'  ! Tellurium
      namatom(53) = 'I '  ! Iodine
      namatom(54) = 'Xe'  ! Xenon
      namatom(55) = 'Cs'  ! Cesium
      namatom(56) = 'Ba'  ! Barium
      namatom(57) = 'La'  ! Lanthanum
      namatom(58) = 'Ce'  ! Cerium
      namatom(59) = 'Pr'  ! Praseodymium
      namatom(60) = 'Nd'  ! Neodymium
      namatom(61) = 'Pm'  ! Promethium
      namatom(62) = 'Sm'  ! Samarium
      namatom(63) = 'Eu'  ! Europium
      namatom(64) = 'Gd'  ! Gadolinium
      namatom(65) = 'Tb'  ! Terbium
      namatom(66) = 'Dy'  ! Dysprosium
      namatom(67) = 'Ho'  ! Holmium
      namatom(68) = 'Er'  ! Erbium
      namatom(69) = 'Tm'  ! Thulium
      namatom(70) = 'Yb'  ! Ytterbium
      namatom(71) = 'Lu'  ! Lutetium
      namatom(72) = 'Hf'  ! Hafnium
      namatom(73) = 'Ta'  ! Tantalum
      namatom(74) = 'W '  ! Tungsten
      namatom(75) = 'Re'  ! Rhenium
      namatom(76) = 'Os'  ! Osmium
      namatom(77) = 'Ir'  ! Iridium
      namatom(78) = 'Pt'  ! Platinum
      namatom(79) = 'Au'  ! Gold
      namatom(80) = 'Hg'  ! Mercury
      namatom(81) = 'Tl'  ! Thallium
      namatom(82) = 'Pb'  ! Lead
      namatom(83) = 'Bi'  ! Bismuth
      namatom(84) = 'Po'  ! Polonium
      namatom(85) = 'At'  ! Astatine
      namatom(86) = 'Rn'  ! Radon
      namatom(87) = 'Fr'  ! Francium
      namatom(88) = 'Ra'  ! Radium
      namatom(89) = 'Ac'  ! Actinium
      namatom(90) = 'Th'  ! Thorium
      namatom(91) = 'Pa'  ! Protactinium
      namatom(92) = 'U '  ! Uranium

      RETURN
      END


      subroutine partf(jatom,ion,temp,ktemp,u)
C      PARAMETER (NQT=200)
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

*
* yields partition functions with polynomial data from
* ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621.
* ln u(temp)=sum(a(i)*(ln(temp))**i) 0<=a<=5
* Input:
* jatom = element number in periodic table
*       NOT implemented for molecules; see ref. in ref.
* ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
* temp  = array with ktemp values of the temperature
* ktemp = number of temperatures for which partf. is calculated
* Output:
* u     = partf. (linear scale) for iat,ion at the ktemp temperatures
*
C      implicit none
*
      integer ktemp,i,j,k,jatom,ion,iread
      real*8 temp(ntemp),u(ntemp),sp,spec
      real*8 a(0:5,1:3,1:92),ulog,t,aa(0:5)
*
      save iread,a
*
      data iread /0/
*
      if(iread.ne.1) then
* read data if first call:

        open(67,file= '/p7/hpg688/hpg688/data/irwin.dat'
     *    ,status='old',readonly)

        read(67,*)
        read(67,*)
        do 20 j=1,92
          do 10 i=1,3
            if(j.eq.1.and.i.eq.3) goto 10
            sp=dfloat(j)+dfloat(i-1)/100.
            read(67,*) spec,aa
            if(sp.ne.spec) then
              print *,'sp.ne.spec:',sp,spec
              stop 'Unexpected data read'
            else
              do 1 k=0,5
                a(k,i,j)=aa(k)
    1         continue
            endif
   10     continue
   20   continue
        close(67)
        iread=1
      endif


      if(ktemp.gt.ntemp) stop 'increase ntemp for partf'
      do 30 i=1,ktemp
        if(temp(i).lt.600.) then
          stop 'partf; temp<600 K'
        else if(temp(i).gt.16000.) then
          stop 'partf; temp>16000 K'
        endif

	t=dlog(dble(temp(i)))
        ulog=   a(0,ion,jatom)+
     &       t*(a(1,ion,jatom)+
     &       t*(a(2,ion,jatom)+
     &       t*(a(3,ion,jatom)+
     &       t*(a(4,ion,jatom)+
     &       t*(a(5,ion,jatom))))))
        u(i)=dexp(ulog)
   30 continue

      return
      end
C
C ________________________________________________________________
      subroutine voigt_init                                        !
C                                                                 |
C ----------------------------------------------------------------|
C                                                                 |
C  Initiation of calculation of Voigt line profile.                                   |
C  The profile is calculated for each temperature T = TMIN to     |
C  TMAX with steps TSTEP. The Voigt profile is also dependent     |
C  on the gas pressure, Pg, and Pg is therefore interpolated from |
C  the model atmosphere values to values corresponding to T.      |
C  As the microturbulence can be different for each molecule,     |
C  through XIturb as well as through kT/M, the Voigt profile will |
C  be a function of molecule/atom/ion as well.                    |
C                                                                 |
C ________________________________________________________________|
C
       implicit real*8 (a-h,o-z)
       PARAMETER (N99=99)
C      PARAMETER (MMOL=12,NTEMP=12,MT1=NTEMP+1,N99=99)
C      PARAMETER (NWL=149000)
       include 'atomparameter.inc'
C
       DIMENSION TLOG(NDP),PGLOG(NDP),PGDER(NDP),W(3,N99)
     &  ,TL(NDP),TLDER(NDP)
     &  ,PeLOG(NDP),PeDER(NDP)
     &  ,PhLOG(NDP),PhDER(NDP)
     &  ,EMULOG(NDP),EMUDER(NDP)
     &  ,PPLOG(NDP),PPDER(NDP)
     &  ,U(N99),H0(N99),H1(N99),H2(N99),H3(N99),H4(N99)
     &  ,H0DER(N99),H1DER(N99),H2DER(N99),H3DER(N99),H4DER(N99)
     &  ,tprofgaus(9999,ntemp),tprofvoigt(9999,ntemp)

      common /chw/ hwst(ntemp), hwvw(ntemp)
       character profile*5,atmos_voigt*45,atmos*45
      common/cvoigt/ pgas(mt1),p_el(mt1),p_h(mt1),hwda(ntemp)
     &              ,pemu(mt1),pp(mt1,nspec)
     &              ,profile,atmos_voigt
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      common /cph/ph(ndp)
      common /cprespp/prespp(ndp,nspec)           !partial pressures in dyn/cm^2
C     COMMON /CMODEL/TEFF,G,ABUND(16)  
C    & ,TAU(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PH(NDP)
C    & ,PTURB(NDP),XKAPR(NDP),NTAU,NLP,ATMOS
      COMMON/COS/WNOS(NWL),OPACOS(NTEMP,NWL),WLOS(NWL),WLOSSTEP(NWL)
     *   ,TMOL(NTEMP),qvib(ntemp)
     *   ,osresl,wnos_first,wnos_last,losresl,ktemp,nwnos
      CHARACTER IDEN*15,ID2*115
      COMMON /CPROFIL/ PROF(MT1,0:100)
      COMMON/CISO/vkmssp,wgtmola,rel_iso_c(3),pc3(3)
     & ,reliso(5),kiso,nhalf,jderiv,kiso_c
     & ,linesc3_direct,lineshcn_direct,linesc2h2_direct
      COMMON /COPAC/HALF,sig_nu(ntemp),hckt(ntemp)
     *    ,sqr_pi_m,sig_sqr,sumop_c(ntemp),a17,n16,n17,npr,j
      COMMON /CCHROMOSPHERE/ nchrom,lchrom,itmin

       
C
C Data for Aller's method of Voigt profile calculation:
      DATA PI,CVEL,ALPHA0,T0,P0 /3.1415926,2.99793d10,0.05,273.
     &                           ,1.01325d6/
      DATA BDA/8.31327d7/
C      DATA BOLTZ,AMU,BDA/1.38042D-16,1.6605D-24,8.31327d7/
C CVEL~velocity of light in cm/s, ALPHA0~collision damping parameter for 
C CO (used for all other molecules, too) measured at T0,P0 = 273K, 1 atm 
C (=1.01325e6 dyn/cm2). BOLTZ~Boltzmanns constant in erg*K-1, AMU~atom
C mass unit on 12C scale, BDA~BOLTZ/AMU, 
C
C
C _____________________________________________________ Voigt profile:
C Calculate the log-log spline interpolated values, Pgas, of the total
C gas pressure in depths corresponding to T(OPAC), from Pg given in the
C NTAU model atmosphere layers.

C
C 1. calculate the derivative of Pg at the NTAU model atmosphere depths:
C

C 1A.Check if T is constant upto a certain value K
C
      N1=1
C      DO 201 I=2,NTAU
C      TDIF=T(I)-T(I-1)
C      IF (TDIF.GT.0) GO TO 202
C      N1=I
C201   CONTINUE
C202   CONTINUE
C
C If we compute a chromospheric OS, find the temperature minimum
C if there is no chromosphere, itmin is in general the top most 
C layer of the atmosphere, but for numerical reasons it can be 
C a few layers down for very flat temperatures, and in that case 
C we skip the very first upper layers
             itmin = 1
             itmax2= 1
             do 161 i = 2,ntau
             if(t(i) .lt. t(i-1)) itmin = i
161          continue
C            tmin = t(itmin)
             n1=itmin
C if there is a secondary maximum, i.e. a cooler region above the chromosphere,
C then for now ignore that in the OS (this is for strongly illuminated
C atmospheres). The maximum temperature (i.e. the top of the chromosphere
C is at itmax2; i.e. we compute the chromospheric OS from itmax2 to itmin-1).
          if(itmin.ge.2) then
             do 162 i = 2,itmin
             if(t(i) .gt. t(i-1)) itmax2 = i
162          continue
          end if
          IF(LCHROM.EQ.1) n1=itmax2
C

C --- in case the atmosphere has a chromosphere (lchrom=1), we reverse
C interpolation, because tb04at demands monotonically increasing function.

      KN=0
      NTAUMAX = NTAU
      IF(LCHROM.EQ.1) NTAUMAX = ITMIN
      DO 200 K=N1,NTAUMAX
      KN=KN+1
      Kj=KN
C       if(t(k)*pg(k)*pe(k)*ph(k) .eq.0.0) 
C     & write(6,251) k,kn,t(k),pg(k),pe(k),ph(k)
C      write(6,251) k,kn,t(k),pg(k),pe(k),ph(k)
251   format('t,pg,pe,ph:',2i3,f7.1,1p3e13.4)
      IF(LCHROM.EQ.1) Kj=NTAUMAX-K+1
      TLOG(Kj)=LOG(T(K))
      TL(Kj)=T(K)
      PGLOG(Kj)=LOG(PG(K))
      PELOG(Kj)=LOG(PE(K))
      PHLOG(Kj)=LOG(PH(K))
      EMULOG(Kj)=LOG(EMU(K))
200   CONTINUE
C
      write(6,*) 'N1,NTAUMAX,ntau=',N1,NTAUMAX,ntau
C      DO 256 K=N1,NTAUMAX
C      write(6,251) k,kn,t(k),pg(k),pe(k),ph(k)
C256   continue

C     CALL TB04AT(NDP,KN,TLOG,PGLOG,PGDER,W)
      CALL TB04AT(NDP,KN,TL,PGLOG,PGDER,W)
      CALL TB04AT(NDP,KN,TL,PELOG,PEDER,W)
      CALL TB04AT(NDP,KN,TL,PHLOG,PHDER,W)
      CALL TB04AT(NDP,KN,TL,EMULOG,EMUDER,W)
C
C 2. calculate the gas pressure, Pgas, in the depths corresponding to
C    the temperatures TMOL(k), k=1,ktemp
C
C     write(6,*) ' interpolated model Pgas, Pe, P(H) values at tmol:'
      write(6,*) ' interpolated model Pgas values at tmol:'
      write(6,*) ' k  temp  Pgas Pe  PP_H:'
      DO 205 kt = 1,ktemp
      temp = tmol(kt)
C     TINT=LOG(TEMP)
      TINT=TEMP
      PGAS(kt)=EXP(TG01BT(-1,NDP,KN,TL,PGLOG,PGDER,TINT))
      P_EL(kt)=EXP(TG01BT(-1,NDP,KN,TL,PeLOG,PeDER,TINT))
      P_H(kt)=EXP(TG01BT(-1,NDP,KN,TL,PhLOG,PhDER,TINT))
      PEMU(kt)=EXP(TG01BT(-1,NDP,KN,TL,EMULOG,EMUDER,TINT))
      write(6,252) kt,temp,pgas(kt),p_el(kt),p_h(kt)
c     write(6,252) kt,temp,pgas(kt)
252   format(i3,f9.1,1p3e13.4)
205   CONTINUE

C now do the same for all the partial pressures. This is necessary to compute
C cm2/g* abs.coef. for atoms, ions, and molecules

      DO 220 kspec = 2,nspec
      KN=0
      DO 221 K=N1,NTAUMAX
      KN=KN+1
      Kj=KN
      IF(LCHROM.EQ.1) Kj=NTAUMAX-K+1
      PPLOG(Kj)=LOG(PRESPP(K,kspec))
221   continue
      CALL TB04AT(NDP,KN,TL,PPLOG,PPDER,W)
      DO 225 kt = 1,ktemp
      temp = tmol(kt)
      TINT=TEMP
      PP(kt,kspec)=EXP(TG01BT(-1,NDP,KN,TL,PPLOG,PPDER,TINT))
225   continue
220   continue


C --------------------------------------------------------------------
C  Voigt profile from L.Aller's method, The Atmosphere of the Sun and 
C  Stars, 1963, p325 and table 7-1 ( file: VOIGT.TABLE, here).
C  New extended table from D.F.Gray 1992 in voigt_gray92.table (Jan2002).
C
C     OPEN (UNIT=9,FILE='voigt.table',STATUS='OLD',readonly)
      OPEN (UNIT=9,FILE='/p7/hpg688/hpg688/data/voigt_gray92.table'
     &     ,STATUS='OLD',readonly)
C
      READ(9,411) DUM
      READ(9,411) DUM
411   format(a25)

      DO 410 I=1,1000
      READ(9,*,END=499) U(I),H0(I),H1(I),H2(I),H3(I),H4(I)
      NHDIM=I
410   CONTINUE
499   CONTINUE
      CLOSE(9)

      write(6,*) ' We read ',nhdim,' values for Voigt profile comp.'
C
C 1. calculate the derivative in the NHDIM Voigt table points:
C
      CALL TB04AT(N99,NHDIM,U,H0,H0DER,W)
      CALL TB04AT(N99,NHDIM,U,H1,H1DER,W)
      CALL TB04AT(N99,NHDIM,U,H2,H2DER,W)
      CALL TB04AT(N99,NHDIM,U,H3,H3DER,W)
      CALL TB04AT(N99,NHDIM,U,H4,H4DER,W)
C
C     CONL=ALPHA0*SQRT(T0)/P0
      XI=vkmssp*1.E5                  ! XI in cm/s; vkmssp in km/s
      ncall = 0
      ws_save = 0.
      it_save = 11
      return


      entry voigt_entry(it,prof1,ws,dw)

C     if(ws.eq.ws_save .and.it.eq.it_save) then
      if(ws.eq.ws_save ) then
      ncall = ncall + 1
      else                      !e.i., a new line
      ncall = 1
      ws_save = ws
      end if
      if(it.eq.it_save) then
        nit = nit + 1
      else
        it_save = it
        nit = 1
      end if
        

C Several pieces from this method is from Bernhard Aringer, January 2002):


C Gauss:  prof = (sqr_pi_m/sig_sqr)*exp(-1.*((dw/sig_sqr)**2))

      voa = hwda(it)/sig_sqr
      vob = dw/sig_sqr

       IF (voa .GE. 0.5) THEN  !Wenn voa>0.5 Voigt-Funktion aus numerischer Integration
                    CALL vofac(voa,vob,vo)
                    PROF1 = (sqr_pi_m/sig_sqr) * vo
                    GOTO 155
       ENDIF

      IF (vob .GE. u(nhdim)) THEN
         voigt0 = 0.0
         voigt1 = 0.56419/(vob**2) + 0.846/(vob**4)
         voigt2 = 0.0
         voigt3 = -0.56/(vob**4)
         voigt4 = 0.0

C         if(v.gt.12.) then
C         prof1 = a/v**2/sqr_pi_m *1.462d0
      ELSE

         VOIGT0=TG01BT(-1,N99,NHDIM,U,H0,H0DER,Vob)
         VOIGT1=TG01BT(-1,N99,NHDIM,U,H1,H1DER,Vob)
         VOIGT2=TG01BT(-1,N99,NHDIM,U,H2,H2DER,Vob)
         VOIGT3=TG01BT(-1,N99,NHDIM,U,H3,H3DER,Vob)
         VOIGT4=TG01BT(-1,N99,NHDIM,U,H4,H4DER,Vob)
      END IF

      PROF1 = (sqr_pi_m/sig_sqr) *
     & (VOIGT0+ voa*(VOIGT1+ voa*(VOIGT2+ voa*(VOIGT3+voa*VOIGT4))) )


155     continue


      if(nit.le.1) then
        profgs = (sqr_pi_m/sig_sqr)*exp(-1.*((dw/sig_sqr)**2))
C        IF (voa .GE. 0.5) THEN
C        write(6,3344) it,profgs,prof1,hwvw(it),hwst(it),voa,vo
C        ELSE
C        write(6,3344) it,profgs,prof1,hwvw(it),hwst(it),voa,
C     &    voigt0,voa*voigt1,voa**2*voigt2
C     &    ,voa**3*voigt3,voa**4*voigt4
C        END IF
3344  format(i2,1p5e8.1,1x,5e8.1)
      end if


      ncall = 999
      if(ncall.gt.1) go to 990

      nlininf = nlininf + 1
      lhalf = 1
      whalf = half
      write(6,*) ' the profile of linr nr ',nlininf,' with half= ',whalf

4041  continue


      do 4040 lt=1,ktemp
      voigtsum = 0.
      gausssum = 0.
      tsig_sqr = ws * sig_nu(lt)
      a = hwda(lt)/tsig_sqr
      vob = dw/tsig_sqr
      voa = a

      DO 420 JQ=0,WHALF
      dv = tsig_sqr * float(jq)
      dvm1 = tsig_sqr * float(jq-1)
      dvp1 = tsig_sqr * float(jq+1)
      delv = (dvp1-dvm1)/2.
      if(jq.eq.0) delv = dvp1/2.
      v = dv/tsig_sqr
      if(v.ge.3.) then
      PROF2 = (sqr_pi_m/tsig_sqr) * a/v**2/sqr_pi_m
      PROF3 = a/v**2/sqr_pi_m *1.462d0
      end if

      vob = v
      IF (vob .GE. u(nhdim)) THEN
         voigt0 = 0.0
         voigt1 = 0.56419/(vob**2) + 0.846/(vob**4)
         voigt2 = 0.0
         voigt3 = -0.56/(vob**4)
         voigt4 = 0.0

C      if(v.gt.12.) then
C      prof1 = a/v**2/sqr_pi_m *1.462d0
C      if(v.gt.12.) then
C      prof1 = prof3

      ELSE
         VOIGT0=TG01BT(-1,N99,NHDIM,U,H0,H0DER,V)
         VOIGT1=TG01BT(-1,N99,NHDIM,U,H1,H1DER,V)
         VOIGT2=TG01BT(-1,N99,NHDIM,U,H2,H2DER,V)
         VOIGT3=TG01BT(-1,N99,NHDIM,U,H3,H3DER,V)
         VOIGT4=TG01BT(-1,N99,NHDIM,U,H4,H4DER,V)
      END IF

      PROF1 = (sqr_pi_m/tsig_sqr) *
     & (VOIGT0+ a*(VOIGT1+ a*(VOIGT2+ a*(VOIGT3+a*VOIGT4))) )
      gaussprof = (sqr_pi_m/tsig_sqr)*exp(-((dv/tsig_sqr)**2))
       IF (voa .GE. 0.5) THEN  !Wenn voa>0.5 Voigt-Funktion aus numerischer Integration
                    CALL vofac(voa,vob,vo)
                    PROF1 = (sqr_pi_m/tsig_sqr) * vo
       ENDIF


       voigtsum = voigtsum + prof1 * delv
       gausssum = gausssum + gaussprof * delv

       tprofgaus(jq,lt) = gaussprof
       tprofvoigt(jq,lt) = prof1

 420   CONTINUE

       write(6,422) lt, tmol(lt), 2.*gausssum, 2.*voigtsum
422    format(' lt, temp, gsusssum, voigtsum: ',i3,f8.1,1p2e12.3)

4040  continue

      if(lhalf.ge.3) go to 990
       lhalf = lhalf+1
       whalf = whalf*10.
       go to 4041

Cc     do 428 jq=0,half
C      do 428 jq=0,5
C      tv = float(jq)
CC     write(6,421) tv, (tprofvoigt(jq,lt),lt=1,ktemp)
C421   format(f6.1,1p9e11.3)
C428   continue


990   CONTINUE   
      RETURN
      END
C
C
C
      SUBROUTINE MODEL_old
C   
C---------------------------------------------------------- 
C                                                         I 
C  subroutine to read T,Pg structure from a model atmosphere
C                                                         I 
C  TEFF  : effective temperature                          I 
C  G     : acceleration of gravity                        I 
C  ABUND : abundance of H,HE,C,N,O,...... (not log() !)   I
C  NTAU  : number of depth points                         I 
C  TAU   : Rosseland optical depth                        I 
C  T       temperatures                                   I 
C  PE    : electron pressures                             I 
C  PG    : gas pressures                                  I 
C  PREAD : radiation pressures                            I 
C  PTURB : turbulent pressures                            I 
C  XKPAPR: Rosseland absorption coefficients              I 
C                                                         I 
C---------------------------------------------------------- 
C
      implicit real*8 (a-h,o-z)
C      PARAMETER (NDP=55)
      include 'atomparameter.inc'
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      common /cph/ph(ndp)
      CHARACTER IDEN*15,ID2*115,ATMOS*45
     & ,name_listmo(16)*8

C     atmos = '/ste1/uffegj/c2/sauval/sun03.dat'

C
      OPEN (UNIT=22,FILE=ATMOS,STATUS='OLD',readonly)
C
C
10      FORMAT(/' ATMOSPHERE model for Voigt profile comp.= ',A45,/
     &     ' (Teff,log(g),C/O) =',F6.0,F5.1,F8.2)
11      FORMAT(A15,A115)
12      FORMAT(29X,F7.0)
13      FORMAT(29X,E11.3)
14      FORMAT(24X,F14.5)
15      FORMAT(16F5.2)
811     FORMAT (A15)
C
C

C
C First count the number of depth points in this atmosphere:
        DO 830 LINE=1,10000
C
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'R R E C')
        if (IDEN.NE.' C O R R E C T ' .and. 
     *     IDEN.NE.'C O R R E C T I') go to 830
C       write(6,811) iden
C       write(6,*) ' ipos = ',ipos
C        IF (IPOS.EQ.0) GO TO 830
        NLIN=0
831     continue
        NLIN = NLIN+1
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'TAUROS')
        IF (IPOS.EQ.0 .and. NLIN.LE.99) GO TO 831
        DO 832 I=1,99                   ! 99=dimension of depth variables
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'M O D E L  A T')
        IF (IPOS.NE.0) GO TO 833
        IF(IDEN(1:6).EQ.'      ') GO TO 833
        NTAU = I
832     CONTINUE
830     CONTINUE
833     CONTINUE
899     CONTINUE
C
        REWIND(22)
C
        IFLAG = 0
        LOGKILL = 0
        DO 4100 MT=1,10000
        READ(22,11,END=4199) IDEN,ID2
C
        IF (IDEN.EQ.'0M O D E L  P A' .OR. IDEN.EQ.' M O D E L  P A') 
     *    THEN
        IFLAG=IFLAG+1            !IFLAG=1
1005    CONTINUE
        READ(22,12) TEFF
        IF (TEFF.EQ.0.0) GO TO 1005
        READ(22,13) FLUX
        READ(22,14) G
        END IF
C
        IF(IDEN.EQ.'  LOG. ABUNDANC') THEN
        IFLAG=IFLAG+1            !IFLAG=2
        READ(22,11) IDEN,ID2
        READ(22,*) (ABUND(I),I=1,16)
        AHA=ABUND(1)
        SUMABUND = 0.D0
        DO 4105 I=1,16
        ABUND(I)=10.**(ABUND(I)-AHA)
        SUMABUND = SUMABUND + ABUND(I)
4105    CONTINUE
        WRITE(6,10) ATMOS,TEFF,LOG10(G),ABUND(3)/ABUND(5)
C       WRITE(6,41051) (ABUND(IE),IE=1,16)
        END IF
41051   FORMAT(1P8E10.3)
C
C
        IF ( IDEN.NE.'1M O D E L  A T' .AND. IDEN.NE.' M O D E L  A T') 
     *    GO TO 4101
        IFLAG=IFLAG+1            !IFLAG=3
1006    READ(22,11,END=4199) IDEN,ID2
        LOOP=LOOP+1
        IF (LOOP.LE.99.AND.IDEN(1:5).NE.'    K'.AND.IDEN(1:5).NE.'   K '
     &       .AND.IDEN(1:5).NE.'  K  '.AND.IDEN(1:5).NE.' K   '
     &       .AND.IDEN(1:5).NE.'0   K'.AND.IDEN(1:5).NE.'0  K '
     &   .AND.IDEN(1:5).NE.'0 K  '.AND.IDEN(1:5).NE.'0K   ') GO TO 1006
        LOOP=0
        WRITE(6,*) ' K1,TAU(I),TAUS,Z,T(I),PE(I),PG(I),XKAPR(I)'
        DO 1410 I=1,NTAU
        READ(22,*) K1,TAUMOD(I),TAUS,Z,T(I),PE(I)
     &  ,PG(I),PRAD(I),PTURB(I),XKAPR(I)
         IF(I.EQ.1.OR.I.EQ.NTAU) WRITE(6,4109)
     &  K1,TAUMOD(I),TAUS,Z,T(I),PE(I),PG(I),XKAPR(I)
4109    format(i2,1p3e10.3,0pf8.1,1p3e10.3)
1410    CONTINUE
4101    CONTINUE

        IPOS = INDEX(IDEN,'L O G A R I T')
        IF ((IPOS.EQ.0).OR.(LOGKILL.EQ.1)) GO TO 4107
        IFLAG=IFLAG+1  !IFLAG=5
C       WRITE(6,*) ' IFLAG5 = ',IFLAG
        LOGKILL = 1    !To prevent another search in next loop (AB 1995-05)

2135    CONTINUE        !read (next) block of partial pressures
        ILOOP = 0
5831    READ(22,811,END=999) IDEN
        IF(IDEN.EQ.' A B S O R P T ' .or.
     &               IDEN.EQ.' P A R T I A L') GO TO 4106  !assumed end of PP blocks
        ILOOP = ILOOP + 1
        IPOS = INDEX(IDEN,' K  ')
        IF (ILOOP.GE.1000) STOP ' ***Error: I found no k in PP'
        IF (IPOS.EQ.0) GO TO 5831
C when here, I identified the line with names of molecules for this block
           backspace(22)
           read(22,2143) (name_listmo(i),i=1,16)
C          write(6,2143) (name_listmo(i),i=1,16)
2143  FORMAT(5x,16a8)


        DO 2130 I=1,NTAU
        READ(22,*) K1,(PRESMP(I,K),K=1,16)
C        if(i.eq.1.or.i.eq.47) write(6,2131) I,(PRESMP(I,K),K=1,10)
C2131    format(i3,10f7.3)
2130    CONTINUE
        if(name_listmo(1).eq.'  P(H)  ') then
        do 2134 i=1,ntau
        ph(i) = 10.**presmp(i,1)
2134    CONTINUE
C        write(6,*) ' I P(H)  P(H)/Pg in each 5th layer in input model'
C        do 2133 i=1,ntau,5
C        write(6,2132) I, PRESMP(i,1),(10.**PRESMP(i,1))/pg(i)
C2133    CONTINUE
2132    format(i3,f7.3,f7.4)
        end if

        go to 2135

4106    CONTINUE
4107    CONTINUE



4100     CONTINUE
4199     CONTINUE
C
      CLOSE (22)

999     CONTINUE

C
      RETURN
      END

C
C----------------------------------------------------------------------------

C
      SUBROUTINE vofac(voa,vob,vo)
C
C This routine produces the Voigt function for voa=a>0.5!
C Do not use it for a<0.01 without decreasing the stepwidth
C for the integration!
C
C AB2002
C

      IMPLICIT REAL*8 (a-h,o-z)

      a = voa
      u = vob
      x1 = 0.0
      x2 = 0.0
      vo = 0.0

C Suchen der Grenzen fuer Integration unter der Annahme, dass Maximum etwa
C bei 0 liegt, stimmt fuer a>0.5

 100  eval = vofunc(a,u,x1)
         IF (eval .GT. 0.0000001) THEN
             x1 = x1 - 0.01
         GOTO 100
         ELSE
             xlo = x1
         ENDIF


 200  eval = vofunc(a,u,x2)
         IF (eval .GT. 0.0000001) THEN
           x2 = x2 + 0.01
         GOTO 200
        ELSE
           xup = x2
        ENDIF

C Integration

      DO i=1,1000
         x1 = xlo + i*(xup-xlo)/1000.0
         x2 = xlo + (i-1)*(xup-xlo)/1000.0
         val = 0.5 * ABS(x1 - x2) * (vofunc(a,u,x1) + vofunc(a,u,x2))
         vo = vo + val
      ENDDO

      RETURN
      END

      FUNCTION vofunc(a,u,x)
      IMPLICIT real*8 (a-h,o-z)
      vofunc = a/3.141592654 * EXP(-1.0*(x**2)) / ((u-x)**2 + a**2)
      RETURN
      END

C
C
C------------------------------------------------------------
C                                                           I
      SUBROUTINE GEM_INIT
C                                                           I
C reads the names and indiexes etc to be used in Gibbs      I
C minimum energy routines                                   I
C------------------------------------------------------------
      implicit real*8 (a-h,o-z)
C      PARAMETER (NDP=55,nspec=892,natms=49)
C      PARAMETER (NDP=55)
      include 'atomparameter.inc'

      dimension abund(natms)
      common /cabundinp/abundatms_inp(natms)
      common /ctotabk/totabk(ndp,natms)
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
      character sunz*1
      NAMELIST /ABUNDANCES/sunz,zscale,abundatms_inp
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
       natms_gem = natms
       nions_gem = 127
       nspec_gem = nspec
        open(unit=2,file='/p7/hpg688/hpg688/data/gfits.data'
     &   ,status='old',readonly)   ! extract molecular etc names from here
        do 2125 i=1,nspec
            read(2,766) name_gem(i)
            read(2,*) ajunk1,ajunk2,ajunk3,ajunk4,ajunk5   ! Fit data
2125    continue
766     format(a8)
        close(unit=2)

       
        open(unit=2,file='/p7/hpg688/hpg688/data/elabund.dat'
     &   ,status='old',readonly)   ! read the elemental abundances
        do 125 i=2,natms
            read(2,*) kelem, abund(i)
125     continue
        write(6,131) abund(2),abund(3),(abund(i),i=7,9),abund(26)
131     format('GEM solar inp. H,He,C,N,O,Fe:',8f7.2)
C        write(6,231) (name_gem(i),i=4,natms)
C        write(6,232) (abund(i),i=4,natms)
C231     format(16a5)
C232     format(16f5.2)

C if sunz = yes, we use solar abundances from unit 2
        do 300 i=1,natms
300     abundatms_inp(i) =  999.9

        READ(12,ABUNDANCES)
        if(sunz.eq.'y'.or.sunz.eq.'Y') then
          write(6,*) 
     &    'We adopted solar abundances from unit=2=elabund.dat'
          go to 309
        end if

        if(zscale.ne.1.) then
          do 310 i=4,natms                     !1,2,3=Pe,H,He not scaled
          abund(i) = log10(zscale) + abund(i)
310       continue
          write(6,311) zscale
311       format
     &    ('We adopted scaled solar abundances; scl.factor:',1pe12.3)
        go to 309
        end if

          nchange_ab = 0
        do 320 i=2,natms
        if(abundatms_inp(i).le.900.0) then
          abund(i) = abundatms_inp(i)
          write(6,*) 'i,abund(i) =', i,abund(i)
          nchange_ab = nchange_ab + 1
        end if
320     continue
        write(6,322) nchange_ab
322     format('We adopted',i3,' individual abundances from input')
        if(natms-nchange_ab-3.gt.0) write(6,323) natms-nchange_ab-3    !1,2,3=e,H,He
323     format('The other',i3,' abundances were solar from inputfile')

309     continue

C       if(sunz.ne.'y'.and.sunz.ne.'Y') then
        write(6,132) abund(2),abund(3),(abund(i),i=7,9),abund(26)
132     format('Adopted abundances of H,He,C,N,O,Fe:',8f7.2)
C       end if


        aha=abund(2)
        sumabund = 0.d0
        do 4105 i=2,natms
        abund(i)=10.**(abund(i)-aha)
        sumabund = sumabund + abund(i)
4105    continue
        do 4205 k=1,ndp
        do 4205 i=2,natms
        totabk(k,i) = abund(i)/sumabund
4205    continue
        close(unit=2)
        write(6,130) totabk(1,2), totabk(1,3), totabk(1,7), totabk(1,8),
     &      totabk(1,9), totabk(1,26)
130     format('Rel.abund,H,He,C,N,O,Fe:',1p6e9.2)


      RETURN
      END
C
C
C
      SUBROUTINE TIMEX
      implicit real*4 (a-h,o-z)
C   This routine prints total accumulated time and time spent since
C  last call.
C
      CHARACTER*20 FORM/'(A,F05.2,A,F05.2,A)'/
      ENTRY TIMEX0
C      TIME_LAST=SECOND(DUMTIM)  
      call etime(time_last)
      time_0 = time_last
C      TIME_LAST=SECOND()
      RETURN
      ENTRY TIMEX1
C      TIME_NOW=SECOND(DUMTIM)
C      TIME_NOW=SECOND()
      call etime(time_now)
      L1=5
      IF(TIME_NOW.GT.99.) L1=6
      IF(TIME_NOW.GT.999.) L1=7
      IF(TIME_NOW.GT.9999.) L1=8
      IF(TIME_NOW.GT.99999.) L1=9
      IF(TIME_NOW.GT.999999.) L1=10
      WRITE(FORM(5:6),'(I2.2)') L1
      DELTA_TIME=TIME_NOW-TIME_LAST
      tot_time = time_now - time_0
      L1=5
      IF(DELTA_TIME.GT.99.) L1=6
      IF(DELTA_TIME.GT.999.) L1=7
      IF(DELTA_TIME.GT.9999.) L1=8
      IF(DELTA_TIME.GT.99999.) L1=9
      IF(DELTA_TIME.GT.999999.) L1=10
      WRITE(FORM(13:14),'(I2.2)') L1
      write(6,FORM) ' Total time spent: ',tot_time,' sec. Added time: ',
     1 DELTA_TIME,' sec.'
      TIME_LAST=TIME_NOW
      RETURN
      END
C
C
C
C
C
C------------------------------------------------------------
C                                                           I
      SUBROUTINE TPGREAD 
C                                                           I
C reads the input T-Pg structure if not a Marcs model       I
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'
      CHARACTER ATMOS*45
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      common /cph/ph(ndp)
      dimension phe(ndp),trpe(ndp),ptot(ndp),trphe(ndp)
      common /ctotabk/totabk(ndp,natms)
C
C     open(unit=22,file=atmos,status='old')
      open(unit=22,file='tpgpe.dat',status='old')

       ntau = 0
       do 110 k=1,ndp   
       read(22,*,end=999) t(k), pg(k), pe(k)
       ntau = ntau+1
110    continue
999    continue

       write(6,*)' Input T Pg Pe Ptot P_He:'
       do 120 i=1,ntau
       phe(i) = pg(i) * totabk(i,3)/(totabk(i,2)+totabk(i,3))
       trpe(i) = 0.1d0 * pe(i)
C050425ptot(i) = pg(i) + pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * phe(i)) / ptot(i)
       if(i.eq.1.or.i.eq.27.or.i.eq.47)
     &          write(6,782) i,t(i),pg(i),pe(i),ptot(i),trphe(i)
120   CONTINUE
782    format(i3,f8.1,1p4e12.3)

      call tstgem(t,ptot,trphe,trpe,ntau)
C
      return
      end

C                                                           I
C------------------------------------------------------------
C
C
      SUBROUTINE MODEL
C     program gemmodel
C   
C---------------------------------------------------------- 
C                                                         I 
C  SUBROUTINE TO READ IN MODEL DATA FROM ARCHIV FILE.     I 
C                                                         I 
C  TEFF  : EFFECTIVE TEMPERATURE                          I 
C  G     : ACCELERATION OF GRAVITY                        I 
C  ABUND : ABUNDANCE OF H,HE,C,N,O,...... (NOT LOG() !)   I
C  NTAU  : NUMBER OF DEPTH POINTS                         I 
C  TAU   : STANDARD OPTICAL DEPTH                         I 
C  T       TEMPERATURES                                   I 
C  PE    : ELECTRON PRESSURES                             I 
C  PG    : GAS PRESSURES                                  I 
C  PREAD : RADIATION PRESSURES                            I 
C  PTURB : TURBULENT PRESSURES                            I 
C  XKPAPR: STANDARD ABSORPTION COEFFICIENTS               I 
C  RO    : DENSITIES                                      I 
C  EMU   : MEAN MOLECULAR WEIGHT (NUCLEONS PER PARTICLE)  I 
C  PRESMP: MOLECULAR PARTIAL PRESSURES (C2H2=15, HCN=16)  I
C        : (IN DYN/CM2) - NOT LOG() !                     I 
C  XL    : WAVELENGTHS FOR ABSKA, SPRIDA (ANGSTROEM)      I 
C  ABSKA : CONTINUOUS ABSORPTION COEFFICIENTS             I 
C  SPRIDA: CONTINUOUS SCATTERING COEFFICIENTS             I 
C                                                         I 
C---------------------------------------------------------- 
C
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'
      dimension SRCH(2)
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),PRESMP(NDP,99)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,ATMOS
      common /cph/ph(ndp)
      COMMON /CMODEL2/P_MOL(NDP), P_NEU_HCNO(NDP), P_ION_HCNO(NDP)
     & ,P_NEU_HE(NDP), P_ION_HE(NDP), P_NON_HHECNO(NDP), P6_JON(NDP)
     & ,PG_JON(NDP), P6(NDP), PHE(NDP)
      COMMON /CTRAN/TAU(NDP),X(NDP),BPLAN(NDP),HSURF,Y1(6),JTAU
      CHARACTER SRCH*15,IDEN*15,ID2*115,ATMOS*45
      DATA SRCH /'O R R E C T I O','M O D E L  P A '/
C common for storage of input partial pressures
      character*8 namesave
      common/cppinp/ pressave(ndp,99),nsave,namesave(99)

C extra for inclusion of GEM routines:
      common/cmxread/pel_ion(ndp,16,4), elabund(ndp,16), 
     &                pprel(ndp,nspec)        !elabund=#atoms_i/atom_in_gas
      common /cgem/pres_gem(ndp,nspec)    !stored input model pressures in gem_indexes
      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character*8 name_listmo(16),name_ex(250)
     &           ,name_comp1,name_comp2,name_comp3
      dimension ptot(ndp),trphe(ndp),trpe(ndp)
      dimension klistmo(99)
      common /cmarcsgem/ kmp(0:99),kidn(99),niden
      common /cabink/abink(ndp,nspec)
      real*8 junk1,junk2,junk3,junk4,junk5! Dummies

       DO 2430 I=1,NTAU
       DO 2430 J=1,NSPEC
2430     pprel(i,j) = 1.d-40
        DO 5113 I=1,NDP
        do 5131 j=1,nspec
        pres_gem(i,j) = -45.d0
5131    CONTINUE
        DO 5113 J=1,16
        DO 5114 JJ=1,4
        pel_ion(i,j,jj) = 1.0d-30
5114    CONTINUE
        DO 5113 K=1,99
        PRESMP(I,K) = -99.d0
5113    CONTINUE

       do 470 i=1,16
470    name_listmo(i) = '        '

C
C     atmos = '/ste1/uffegj/c2/sauval/sun03.dat'


      OPEN (UNIT=22,FILE=ATMOS,STATUS='OLD',readonly)
C
C
10      FORMAT(/' ATMOSPHERE = ',A45,/' (Teff,log(g),C/O) ='
     &         ,F6.0,F5.1,F8.2)
11      FORMAT(A15,A115)
12      FORMAT(29X,F7.0)
13      FORMAT(29X,E11.3)
14      FORMAT(24X,F14.5)
15      FORMAT(16F5.2)
811     FORMAT (A15)
C
C

C
C First count the number of depth points in this atmosphere:
        DO 830 LINE=1,10000
C
        READ(22,811,END=899) IDEN
        IPOS = 0
C        IPOS = INDEX(IDEN,SRCH(1))  !C O R R E C T I O N S  I N  L A S T
        IPOS = INDEX(IDEN,'R R E C')
        if (IDEN.NE.' C O R R E C T ' .and. 
     *     IDEN.NE.'C O R R E C T I') go to 830
C       write(6,811) iden
C       write(6,*) ' ipos = ',ipos
C        IF (IPOS.EQ.0) GO TO 830
        NLIN=0
831     continue
        NLIN = NLIN+1
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'TAUROS')
        IF (IPOS.EQ.0 .and. NLIN.LE.99) GO TO 831
        DO 832 I=1,150                   ! 99=dimension of depth variables
        READ(22,811,END=899) IDEN
        IPOS = 0
        IPOS = INDEX(IDEN,'M O D E L  A T')
        IF (IPOS.NE.0) GO TO 833
        IF(IDEN(1:6).EQ.'      ') GO TO 833
        NTAU = I
832     CONTINUE
830     CONTINUE
833     CONTINUE
899     CONTINUE
C
        REWIND(22)
        write(6,*) ' I found ',ntau,' depth points in model'
C
        IFLAG = 0
        LOGKILL = 0
        DO 4100 MT=1,10000
        READ(22,11,END=4199) IDEN,ID2
C
        IF (IDEN.EQ.'0M O D E L  P A' .OR. IDEN.EQ.' M O D E L  P A') 
     *    THEN
        IFLAG=IFLAG+1            !IFLAG=1
1005    CONTINUE
        READ(22,12) TEFF
        IF (TEFF.EQ.0.0) GO TO 1005
        READ(22,13) FLUX
        READ(22,14) G
        END IF
C
        IF(IDEN.EQ.'  LOG. ABUNDANC') THEN
        IFLAG=IFLAG+1            !IFLAG=2
        READ(22,11) IDEN,ID2
        READ(22,*) (ABUND(I),I=1,16)
        WRITE(6,10) ATMOS,TEFF,LOG10(G),10.**(ABUND(3)-ABUND(5))
        WRITE(6,41053) (ABUND(IE),IE=1,5),ABUND(15)
        AHA=ABUND(1)
        SUMABUND = 0.D0
        DO 4105 I=1,16
        ABUND(I)=10.**(ABUND(I)-AHA)
        SUMABUND = SUMABUND + ABUND(I)
4105    CONTINUE
        DO 4205 k=1,NDP
        DO 4205 I=1,16
        elabund(k,i) = ABUND(I)/SUMABUND
4205    CONTINUE
C       WRITE(6,41051) (ABUND(IE),IE=1,16)
        END IF
C41051   FORMAT(1P8E10.3)
41053   FORMAT('H,He,C,N,O,Fe:'8f8.3)
C
C
        IF ( IDEN.NE.'1M O D E L  A T' .AND. IDEN.NE.' M O D E L  A T') 
     *    GO TO 4101
        IFLAG=IFLAG+1            !IFLAG=3
1006    READ(22,11,END=4199) IDEN,ID2
        LOOP=LOOP+1
        IF (LOOP.LE.99.AND.IDEN(1:5).NE.'    K'.AND.IDEN(1:5).NE.'   K '
     &       .AND.IDEN(1:5).NE.'  K  '.AND.IDEN(1:5).NE.' K   '
     &       .AND.IDEN(1:5).NE.'0   K'.AND.IDEN(1:5).NE.'0  K '
     &   .AND.IDEN(1:5).NE.'0 K  '.AND.IDEN(1:5).NE.'0K   ') GO TO 1006
        LOOP=0
C        READ(22,11,END=4199) IDEN,ID2
C        IF (IDEN(1:1).EQ.' ') READ(22,11,END=4199) IDEN,ID2
C       WRITE(6,*) ' IFLAG3 = ',IFLAG
        DO 1410 I=1,NTAU
        READ(22,*) K1,TAUMOD(I),TAUS,Z,T(I),PE(I)
     &  ,PG(I),PRAD(I),PTURB(I),XKAPR(I)
        TAU(I)=TAUMOD(I)
C        IF(I.EQ.1.OR.I.EQ.NTAU) WRITE(6,*)
C     &  K1,TAU(I),TAUS,Z,T(I),PE(I),PG(I),XKAPR(I)
1410    CONTINUE
4101    CONTINUE
        JTAU=NTAU
        IF ( IDEN.NE.'1T H E R M O D ' .AND. IDEN.NE.' T H E R M O D ') 
     *    GO TO 4102
1007    READ(22,11,END=4199) IDEN,ID2
        LOOP=LOOP+1
        IF (LOOP.LE.99.AND.IDEN(1:5).NE.'    K'.AND.IDEN(1:5).NE.'   K '
     &       .AND.IDEN(1:5).NE.'  K  '.AND.IDEN(1:5).NE.' K   '
     &       .AND.IDEN(1:5).NE.'0   K'.AND.IDEN(1:5).NE.'0  K '
     &   .AND.IDEN(1:5).NE.'0 K  '.AND.IDEN(1:5).NE.'0K   ') GO TO 1007
        LOOP=0
        IFLAG=IFLAG+1            !IFLAG=4
        I = 1
        limit=NTAU
99999   IF (I.LE.limit) THEN
            READ(22,*) IK1,TR,RO(I),EMU(I),CP,CV,AGRAD,Q,U,V,FCF,IK2
            I=I+1
            GOTO 99999

        END IF
        
4112    CONTINUE
4102    CONTINUE

        IPOS = INDEX(IDEN,'L O G A R I T')
        IF ((IPOS.EQ.0).OR.(LOGKILL.EQ.1)) GO TO 4107
        IFLAG=IFLAG+1  !IFLAG=5
C       WRITE(6,*) ' IFLAG5 = ',IFLAG
        LOGKILL = 1    !To prevent another search in next loop (AB 1995-05)

        nmol = 0
        nfail = 0
        nsave = 0
2135    CONTINUE        !read (next) block of partial pressures
        ILOOP = 0
5831    READ(22,811,END=999) IDEN
        IF(IDEN.EQ.' A B S O R P T ' .or.
     &               IDEN.EQ.' P A R T I A L') GO TO 4106  !assumed end of PP blocks
        ILOOP = ILOOP + 1
        IPOS = INDEX(IDEN,' K ')
        IF (ILOOP.GE.1000) STOP ' ***Error: I found no K in PP'
        IF (IPOS.EQ.0) GO TO 5831
C when here, I identified the line with names of molecules for this block
           backspace(22)
           read(22,2143) (name_listmo(i),i=1,16)
2143  FORMAT(5x,16a8)
           innam = 0
           do 2161 i=1,16
2161       if (name_listmo(i).ne.'        ') innam = innam+1
C innam is number of molecule names in this read block
           niden = nmol

           do 2121 j=1,innam
            kn = 0
            name_comp2 = '        '
            name_comp3 = name_listmo(j)
           do 2127 kc = 1,8
           if(name_comp3(kc:kc).ne.' ') then
            kn = kn+1
            name_comp2(kn:kn) = name_comp3(kc:kc)
           end if
2127       continue
           if(name_comp2.eq.'P(H)    ') name_comp2='H       '
           if(name_comp2.eq.'K       ') name_comp2='XXXXXXXX'
C name_comp2 is now the read name_listmo(j), but starting in position 1
           nsave = nsave+1
           namesave(nsave) = name_comp2

           idsucces = 0
           ipos = 0
           do 2122 i=1,nspec
           name_comp1 = name_gem(i)
C just for the comparison with GEM-names:
C (nothing is changed in any name arrays here)
             if(name_comp1.eq.'HN      ') name_comp1 = 'NH      '
             if(name_comp1.eq.'HNa     ') name_comp1 = 'NaH     '
             if(name_comp1.eq.'HMg     ') name_comp1 = 'MgH     '
             if(name_comp1.eq.'HSi     ') name_comp1 = 'SiH     '
             if(name_comp1.eq.'CSi     ') name_comp1 = 'SiC     '
             if(name_comp1.eq.'ClH     ') name_comp1 = 'HCl     '
             if(name_comp1.eq.'CHN     ') name_comp1 = 'HCN     '
             if(name_comp1.eq.'CSi2    ') name_comp1 = 'Si2C    '
           ipos = index( name_comp2,name_comp1 )
           if(ipos.ne.0) then        !we identified the listmo molecule as a GEM name
              niden = niden + 1
              if(niden.gt.99) stop ' increase dimension invl niden !'
              klistmo(niden) = j
              kidn(niden) = i
              idsucces = 1
C             write(6,2128) name_comp1,name_comp2
C             write(6,2129) niden,i,j,kidn(niden),name_gem(i)
2128          format(' name{gem,listmo}=',2(2x,a8))
2129          format(' niden,i,j,kidn(niden), name_gem(i):',4i4,2x,a8)
           end if
2122    CONTINUE
           if(idsucces.eq.0) then 
                 nfail = nfail + 1
                 name_ex(nfail) = name_comp2
           end if
2121    CONTINUE

C At this point namesave(1-nsave) contain the nsave read pp names from input model.
C The corresponding partial pressures (pp) are being saved in pressave(1-nsave).
C name_gem(kidn{1-niden}) contain those of the input model names which has a gem-name,
C while name_ex(1-nfail) contain the input model names not identified in GEM.
C In contrast to pressave, pres_gem(i,kidn(j),j=1,niden) contain only the subset of
C input model pressures of molecules which are in common between GEM and the input model.
C Note that pres_gem is not the GEM-computed pressures, but input pressures.

        DO 2130 I=1,NTAU
        READ(22,*) K1,(PRESMP(I,K),K=1,innam)
        do 2140 k=1,innam
        js=nsave-innam+k
        if (namesave(js).eq.'H       ') then 
            ph(i) = 10.d0**PRESMP(I,K)
C            if(i.eq.1)
C     &        write(6,2240) name_listmo(1),namesave(js),ph(i)
C2240        format(a8,2x,a8,1pe12.3)
        end if

CC 4 lines from model_old:
CC        if(name_listmo(1).eq.'  P(H)  ') then
CC        do 2134 i=1,ntau
CC        ph(i) = 10.**presmp(i,1)
CC2134    CONTINUE

2140    PRESSAVE(I,js) = 10.d0**PRESMP(I,K)
C       if(i.eq.1) then
C         write(6,2142) (namesave(nsave-innam+k),k=1,innam)
C         write(6,2146) 
C    &  (dlog10( max(1.d-40,pressave(i,nsave-innam+k))),k=1,innam)
CC    &  (pressave(i,nsave-innam+k),k=1,innam)
C       end if
2130    CONTINUE
2142    format(2x,8a8/,8a8)
2146    format(8f8.2/,8f8.2)

         if (niden-nmol.eq.0) then
C         write(6,*) 
C    &    ' we found no gem_molecules in this read block of PP'
          go to 2135
         end if
        DO 2131 I=1,NTAU
        do 2131 j=nmol+1,niden
        pres_gem(i,kidn(j)) = presmp(i,klistmo(j))
2131    CONTINUE

2134    FORMAT(i3,16f8.3)
C          write(6,2143)(name_gem(kidn(j)),j=nmol+1,niden)
C          DO 2133 k=1,NTAU,20
C          write(6,2134)k,(pres_gem(k,kidn(j)),j=nmol+1,niden)
C2133    CONTINUE

        NMOL = niden

        go to 2135

4106    CONTINUE

C..        WRITE(6,5110) NMOL,NTAU
C..5110    FORMAT ( ' We identified the following',i3,
C..     &     ' molecular part.pres.in',i3,' depth layers:')
C..        write(6,2144)(name_gem(kidn(j)),j=1,niden)
C..2144    format(15a5,/15a5,/15a5,/15a5)
C..        write(6,2141) nfail,(name_ex(nf),nf=1,nfail)
C..2141    format(' We didnt identify in GEM the following',i3,
C..     &    ' molecules:'/,10a8/,10a8/,10a8/,10a8/,10a8/,10a8)
        DO 2139 I=1,NTAU
        do 2139 j=1,niden
        presmp(i,j) = pres_gem(i,kidn(j))
2139    CONTINUE

C Now, put into the indexes kmp({1-niden}) the GEM index of the molecules which
C we have absorption coefficients for:
        do 2138 j=1,niden
        if(name_gem(kidn(j)).eq.'CO      ') kmp(0) = kidn(j)
        if(name_gem(kidn(j)).eq.'CH      ') kmp(1) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2      ') kmp(2) = kidn(j)
        if(name_gem(kidn(j)).eq.'SiO     ') kmp(3) = kidn(j)
        if(name_gem(kidn(j)).eq.'CN      ') kmp(4) = kidn(j)
        if(name_gem(kidn(j)).eq.'TiO     ') kmp(5) = kidn(j)
        if(name_gem(kidn(j)).eq.'H2O     ') kmp(6) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H2    ') kmp(7) = kidn(j)
        if(name_gem(kidn(j)).eq.'HCN     ') kmp(8) = kidn(j)
        if(name_gem(kidn(j)).eq.'C3      ') kmp(9) = kidn(j)
        if(name_gem(kidn(j)).eq.'H2      ') kmp(10) = kidn(j)
        if(name_gem(kidn(j)).eq.'He      ') kmp(11) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H2    ') kmp(12) = kidn(j)
        if(name_gem(kidn(j)).eq.'CH4     ') kmp(13) = kidn(j)
        if(name_gem(kidn(j)).eq.'CS      ') kmp(14) = kidn(j)
        if(name_gem(kidn(j)).eq.'C2H     ') kmp(15) = kidn(j)
        if(name_gem(kidn(j)).eq.'OH      ') kmp(16) = kidn(j)
2138    continue
        kmp(11) = 3     !He is not in the listmo
        kmp(8) = 372   !HCN as CHN is not in the listmo

C        write(6,*) 
C     &  ' The opacity bearing species were identified as:'
C        do 2136 j=0,16
C        write(6,2137) j,kmp(j),name_gem(kmp(j))
C2137    format(i3,i4,2x,a4)
2136    continue


4107    CONTINUE


        IF(IDEN.NE.' LAMBDA FOR ABS') GO TO 4103
        IFLAG=IFLAG+1            !IFLAG=6
C       WRITE(6,*) ' IFLAG6 = ',IFLAG
        READ(22,*) NLP                      !usually NLP=20 in input, #20 is 8000.AA
C       WRITE(6,811) IDEN
C       WRITE(6,*) ' NLP = ',NLP
        NLP=NLP-1
C       WRITE(6,*) ' NLP = ',NLP
        READ(22,*) (XL(J),J=1,NLP),XLEXTRA   !cont.wavelength scale is 19 freq.
        WRITE(6,*) (XL(J),J=1,NLP),XLEXTRA
        READ(22,11) IDEN,ID2
        READ(22,11) IDEN,ID2
        READ(22,11) IDEN,ID2
C       WRITE(6,*) ' i,ITAU,log(TAUX),absextra,sprextra: '
C       WRITE(6,*) ' ntau = ',ntau
        DO 4202 I=1,NTAU
        READ(22,*,err=4206) ITAU,TAUX
        READ (22,*,err=4206) (ABSKA(K,I),K=1,NLP),ABSEXTRA
        READ (22,*,err=4206) (SPRIDA(K,I),K=1,NLP),SPREXTRA
C       WRITE(6,4207) i,ITAU,log10(TAUX),absextra,sprextra
4207    format(2i4,f8.2,1p2e12.3)
4202    CONTINUE
        go to 4103
4206    write(6,*) 
     & 'i=',i,'log(tau)=',log10(taux), 
     & ' Problems with reading abs/sprid -> skipping'
4103    CONTINUE
C
4100     CONTINUE
4199     CONTINUE
C
      CLOSE (22)
      CLOSE (2)


C at this point pres_gem(i,j) are 10^-49.0 if not identified in the input model,
C or the log10(input partial pressure / (pg+pe)) if identified in the input model.
C  pprel(i,j) = 10.**presrel is then the input relative partial pressure in 
C model layer i for all nspec species, j=1,nspec.


C      open(unit=8,file='/home/ast/uffegj/falkesg/marcs/tpgpe.dat',
C     &                   status='unknown')
C
C      DO 2610 I=1,NTAU
C      write(8,2611) t(i),pg(i),pe(i)
C2610  continue
C2611  format(f10.3,1p2e12.4)
C      close(unit=8)

      DO 2410 I=1,NTAU
        do 166 j=1,nspec
         presrel= pres_gem(i,j) - dlog10( max(1.d-40,(pg(i)+pe(i))) )
C        if(i.eq.1 .or. i.eq.27 .or. i.eq.47) then
C        if (presrel.ge.-9.0d0)
C    & write(6,266) name_gem(j),j,presrel,10.**presrel,pres_gem(i,j)
C        end if
         pprel(i,j) = 10.**presrel
166     continue
266      format(a8,i4,f7.2,f8.5,f8.3)

       phe(i) = dlog10
     &          (pg(i) * elabund(1,2)/(elabund(1,1)+elabund(1,2)))
       trpe(i) = 0.1d0 * pe(i)
C050425ptot(i) = pg(i) + pe(i)
       ptot(i) = pg(i)
       ptot(i) = 0.1d0 * ptot(i)     !1N/m^2 = 10 dynes/cm^2
       trphe(i) = (0.1d0 * 10.**phe(i)) / ptot(i)
C      if(i.eq.1.or.i.eq.27.or.i.eq.47)
C    &          write(6,782) i,t(i),pg(i),pe(i),ptot(i),trphe(i)
2410  CONTINUE
782    format(i3,f8.1,1p4e12.3)


        call tstgem(t,ptot,trphe,trpe,ntau)

C
      IF (IFLAG.LT.5) THEN
        write(6,*)
     &  ' ERROR: problems with reading of model-input. IFLAG=',iflag
        STOP
      ELSE
        write(6,*) 
     &  ' successful reading of input model atmosphere; iflag=',iflag
      END IF
      GO TO 998
999     CONTINUE
      WRITE(6,*) ' MODEL STOP: HOW COULD THIS BE ???????'
998     CONTINUE

C
      RETURN
C   
      END
C
C

      subroutine tstgem (tcal,pcal,phe,pecal,ntau)

!-----------------------------------------------------------------
!
!     This program is designed to test the GEM subroutines.
!
!     27/9-2000 JFF
!
!-----------------------------------------------------------------

C     implicit none
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

C     integer nspec,natms,i,j,un1,un2,ie,k,keq,ipos,ndp,kdp,ntau
      integer i,j,un1,un2,ie,k,keq,ipos,kdp,ntau
      real*8 tiny,small,t,p
      parameter(un1=46,un2=47)            ! Unit numbers for I/O
      parameter(tiny=1.0d-40,small=1.0d-12)
      real*8 totab(natms)                 ! Bulk elemental abundances
      real*8 totabk                       ! ditto, but possibly depth variable
      common /ctotabk/totabk(ndp,natms)
C      real*8 abink(ndp,nspec)             ! Input and output solutions
      real*8 abink
      real*8 abin(nspec),about(nspec)     ! Input and output solutions
      real*8 gibbs(nspec,5)               ! Gibbs' coefficient array
      real*8 compmtx(natms,nspec)         ! Composition matrix
      logical first,firstit,mol_eq_scr            ! Is this the initial run?
      character*8 name(nspec)
      real*8 junk1,junk2,junk3,junk4,junk5! Dummies
      real*8 initsum,finalsum             ! Total fictitious and partial
                                          ! pressures
      common /cabink/abink(ndp,nspec)
      dimension tcal(ndp),pcal(ndp),phe(ndp),pecal(ndp) 
      common /cgeminp/tcall(ndp),pcall(ndp),phecall(ndp),pecall(ndp) 
      common /cmarcsgem/ kmp(0:99),kgem(99),niden
      common /cgem/pres_gem(ndp,nspec) 
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
      character name_gem*8
      common/gemcom/gibbs,compmtx
      common/cmxread/pel_ion(ndp,16,4), elabund(ndp,16), 
     &                pprel(ndp,nspec)        !elabund=#atoms_i/atom_in_gas
      data firstit, mol_eq_scr /.true., .true./
C     mol_eq_scr = .true.  => compute mol.equilibrium from scratch

      do 445 k=1,ntau
      tcall(k) = tcal(k)
      pcall(k) = pcal(k)
      pecall(k) = pecal(k)
      phecall(k) = phe(k)
445   continue


      first=.true.
      initsum=0.0d0
      do 10 i=1,natms
         initsum=initsum+totab(i)              ! Initial abundance sum
 10   continue
      open(unit=un1,file='tstgem.ud',status='unknown')
      open(unit=un2,file='/p7/hpg688/hpg688/data/gfits.data',
     &                                    status='old',readonly)
         do 121 i=1,nspec
            read(un2,166) name(i)
            read(un2,*) junk1,junk2,junk3,junk4,junk5   ! Fit data
 121      continue
166      format(a8)
      do 13 j=1,ntau
         totab(1) = tiny
         do 131 i=2,natms
131      totab(i) = totabk(j,i)
         finalsum=0.0d0
C???         if(firstit .eq. .false.) then
C???         do 3111 i=1,nspec
C???3111     abin(i) = abink(j,i)
C???         end if
         t=tcall(j)                          ! Temperature scale
         p=pcall(j)                          ! Pressure    scale
         call gem(t,p,totab,abin,about,first,j)  ! Actual call to GEM
C        write(un1,*) 'T= ',t,' K; P= ',p,' Pa'
C         if(j.eq.1 .or. j.eq.27
C     &             .or. j.eq.47) then
C            write(6,266) t,p
C266         format(' Now in TSTGEM: T= ',f6.0,' K; P= ',1pe12.3,' Pa')
C            end if
C        write(un1,*) 'name        abundance             log(abundance)'
         do 11 i=1,nspec
            finalsum=finalsum+about(i)         ! Final abundance sum 
 11      continue
         do 12 i=1,nspec
            abin(i)=about(i)                            ! Recycle data
!            about(i)=about(i)
            junk1=about(i)/finalsum                     ! Normalize
!            write(un1,*) name,about(i),dlog10(about(i)),dlog10(junk1)
         if(j.eq.1 .or. j.eq.27
     &             .or. j.eq.47) then
              abink(j,i) = max(abink(j,i),1.d-40)
C           if(junk1.ge.1.d-9) 
C    &        write(6,267) name(i),dlog10(junk1),dlog10(abink(j,i))
            end if
            abink(j,i) = about(i)/finalsum 
 12      continue
267      format(a8,2f8.3)
         if(j.eq.1 .or. j.eq.27
     &             .or. j.eq.47) then
C        write(6,*) 'Initial sum of abundances: ',initsum
C        write(6,*) 'Final sum of abundances  : ',finalsum
            end if
 13   continue
      close(un1)
      close(un2)
      firstit=.false.

       write(6,*)
     & ' pp of the',niden,' species in common betw presm-gem:'
       do 3151 k=1,niden-9,10
       write(6,3152) (kgem(j),name(kgem(j)),j=k,k+9)
       DO 3151 kd=1,NTAU,20
       write(6,3154) 
     &    kd,(dlog10(10.d0*pcall(kd)*abink(kd,kgem(j))),j=k,k+9)
       write(6,3155) 
     &    kd,( pres_gem(kd,kgem(j)),j=k,k+9)
3151   continue
3152   format(4x,10(i3,a4))
3154   format(i2,1x,10f7.2)
3155   format(i2,'M',10f7.2)




      return
      end



!-----------------------------------------------------------------
!
! Inside MARCS main program do the following: Read Gibbs' energy 
! fits (GIBBSREAD), composition matrix (GEMCMSETUP) and abundances 
! (GEMABSET) into right variables.
! After completion/convergence convert abundances back into MARCS
! format. Set NSPEC and NATMS to their respective values {(246,34) 
! for the common subsample and (49,892) for the full JANAF set}.
!              JFF 18/4-2000
!              JFF 19/9-2000 correction of NSPEC and NATMS values.
! Note, NSPEC is including NATMS, the number of non-elemental 
! species is therefore: NSPEC - NATMS = 212 or 843.
!              JFF 21/9-2000 
!
! Note, the common subsample uses another file-set than the full
! sample. These are not yet converted to GEM-readable format.
!              JFF 27/9-2000
!
! The common sample are now converted to GEM-format, in the files
! TSUCOMP and TSUDATA.
!              JFF 6/11-2000
!
! NSPEC=255, NATMS=48 with TSUCOMP and TSUDATA
!
!-----------------------------------------------------------------

      subroutine gem(temp,pres,totabund,abund1,abund2,first,ktau)

!-----------------------------------------------------------------
!
!     This is the control subroutine, the one called from the
!     main program.
!
!     The algorithm is taken from:
!      White, Johnson & Dantzig; 
!      "Chemical Equilibrium in Complex Mixtures";
!      Journal of Chemical Physics, vol. 28, #5, May 1958, p. 751
!
!     Henceforth referred to as "WJD"
!
!   IMPORTANT!
!     The subroutines 'GEMCMSETUP' and 'GIBBSREAD' should be moved
!     outside this subroutine 'GEM' and be called once from the
!     main program if more than one equilibrium composition is 
!     needed.
!
!     22/09-2000: Jens Falkesgaard
!     10/01-2000: Charge conservation-part of MASSCONS updated
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 temp,pres                              ! T and P in WJD
      real*8 totabund(natms)                        ! b in WJD
      real*8 abund1(nspec),abund2(nspec)            ! y and x in WJD
      real*8 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*8 compmtx(natms,nspec)                   ! a in WJD
      logical first

!--internal variables

      integer i,j,count,conprbl,un                  ! iteration params.
      parameter(un=40)                              ! unit number
      logical sofar,newstep                         ! to step or not...
      real*8 absum1,absum2                          ! y-bar and x-bar
      real*8 gef(nspec)                             ! c in WJD
      real*8 pge(nspec)                             ! partial Gibbs' energies
      real*8 eq18a(natms+1,natms+1)                 ! r in WJD (A-matrix)
      real*8 eq18b(natms+1)                         ! B-vector in WJD
      real*8 indx(natms+1)                          ! vector for solvelin
      real*8 delta(nspec)                           ! relative changes
      real*8 oldnorms(6)                            ! delta vector norms
      real*8 norm,avgnorm,lambda

      common/gemcom/gibbs,compmtx

!      print*,'In GEM'
      call gemcmsetup                            ! Make comp. matrix
      call gibbsread                             ! Read polynomials

      if(first) call gemabset(temp,absum1,totabund,abund1,first,ktau) 
C???  call gemabset(temp,absum1,totabund,abund1,first,ktau) 
                                                 ! Sets up abundances
      count=1

      do 10 i=1,6
         oldnorms(i)=1.0d0                       ! Initialize norms
 10   continue
C        if(ktau.eq.1 .or. ktau.eq.27
C    &             .or. ktau.eq.47)
C    & write(6,*)' In GEM: k=',ktau,' T= ',temp,' P= ',pres
!      open(unit=un,file='gem.ud',status='unknown')
!      write(un,*) 'T= ',temp,'P= ',pres
!      print*,'Main iteration no.: ',count  
 11   continue
      if(count.ge.10000) then                     ! Prevent infinite loop
         print *,'Convergence failure!'
C        do 12 i=1,nspec
C           print*, i,abund2(i)                  ! Prints last composition
C12      continue
         stop
      endif
      call gemgibbs(temp,pres,gef,abund1,absum1,pge) ! Computes partial Gibbs' 
                                                     ! energies
      call gemmatrx(abund1,totabund,pge,eq18a,eq18b) ! Set up system matrix

      j=natms+1                                  ! natms is common; j avoids
                                                 ! case-specific code in the
                                                 ! general subrt. SOLVELIN.

      call solvelin(eq18a,j,j,indx,eq18b)            ! Solve the system
      call gemstep(pge,gef,eq18b,absum1,absum2,abund1,abund2) ! Use solution
      call masscons(abund2,totabund)             ! Enforce mass conservation
      sofar=.true.
      newstep=.true.
      conprbl=0                             ! # species that didn't converge
      norm=0.0d0
      avgnorm=0.0d0
      do 13 i=1,6
         avgnorm=avgnorm+oldnorms(i)        ! average of old delta norms
 13   continue
      avgnorm=avgnorm/6.0d0
      do 14 i=1,nspec
         delta(i)=abund2(i)-abund1(i)
         norm=norm+(delta(i)*delta(i))      ! norm-square of delta
         delta(i)=dabs(delta(i)/abund1(i))
         if((delta(i).gt.1.5d-2)) then      ! If rel. change too big (*)
            sofar=.false.
            conprbl=conprbl+1               ! This species didn't converge
!            write(un,*) 'Convergence problems for species: ',i,delta(i)
         endif
 14   continue
C     print*,'Convergence problems for ',conprbl,' species.'
      norm=sqrt(norm)
!      print*, '|delta|= ',norm,', avgnorm= ',avgnorm

      if(norm.gt.avgnorm.and.norm.ne.0.0d0) then ! |delta| was too large
         lambda=avgnorm/norm                      
         do 15 i=1,nspec
            delta(i)=(abund2(i)-abund1(i))*8.5d-2*lambda ! reduce delta
            abund2(i)=abund1(i)+delta(i)
 15      continue
      endif
      do 16 i=1,5
      oldnorms(i)=oldnorms(i+1)                  ! Update 'old' norms
 16   continue
      oldnorms(6)=norm
      if(.not.sofar) then 
         newstep=.true.                          ! (*) do another loop.
       else 
         newstep=.false.                         ! If acceptable stop (+)
      endif
      if(newstep) then
         absum1=0.0d0
         do 17 i=1,nspec
            abund1(i)=abund2(i)                  ! Move values
            absum1=absum1+abund2(i)
 17      continue
         count=count+1                           ! Update counter
         goto 11
      endif                                      ! (+) ...here.
!      close(un)
C     print*,'Iteration steps completed: ',count
!      print*,'Mass checksum  = ',msum
!      print*,'Charge checksum= ',qsum
      end

************

      subroutine gemcmsetup

!-----------------------------------------------------------------
!
!     This subroutine reads and sets up the composition matrix
!     for GEM. 
!
!     To be used inside main program, MARCS or elsewhere, before
!     call to GEMMATRX. There is only need for one run of this
!     subroutine, provided the bulk composition does not change.
!
!     21/9-2000 : Jens Falkesgaard
!     12/10-2000: Tested and found OK, JFF
!
!-----------------------------------------------------------------

C     implicit none
      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables
     
      real*8 gibbs(nspec,5)                         ! not used here
      real*8 compmtx(natms,nspec)                   ! a in WJD

!--internal variables

      integer i,j,un,tmpint
      parameter(un=41)                              ! unit number
      real*8 tmpreal                                ! dummy variable
      character*16 name
      integer listno,nsp                            ! # and # of elements

      common/gemcom/gibbs,compmtx

C     open(unit=un,file='/home/ast/falkesgd/src/comp',
      open(unit=un,file='/p7/hpg688/hpg688/data/comp',
     &                                   status='old',readonly)
      do 12 i=1,nspec
         read(un,166) name                            ! name, discard
         read(un,*) listno                          ! number, discard
         if((listno+1).ne.(i)) print *,'Consistency error!'
         read(un,*) nsp                             ! read # species
         do 11 j=1,nsp                              ! for each
            read(un,*) tmpint, tmpreal              ! read data
            compmtx(tmpint+1,i)=tmpreal             ! & assign
 11         continue
 12      continue
166   format(a16)
      close(un)
      end

***********************

      subroutine gibbsread

!-----------------------------------------------------------------
!
!     This subroutine reads Gibbs' energy fit coefficients from
!     the file 'gfits.data'.
!
!     21/9-2000 : Jens Falkesgaard
!     12/10-2000: Tested OK, JFF
!
!-----------------------------------------------------------------

!--external variables

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'
      real*8 gibbs(nspec,5)                   ! Gibbs' coefficients
      real*8 compmtx(natms,nspec)             ! composition matrix

!--internal variables

      integer i,un
      parameter(un=42)                        ! unit number
      character*8 name
      real*8 coef0,coef1,coef2,coef3,coef4    ! Gibbs' coefficients, tmp

      common/gemcom/gibbs,compmtx

      open(unit=un,file='/p7/hpg688/hpg688/data/gfits.data',
     &                                   status='old',readonly)
      do 11 i=1,nspec
         read(un,166) name
         read(un,*) coef0,coef1,coef2,coef3,coef4
         gibbs(i,1)=coef0
         gibbs(i,2)=coef1                     ! Read and assign for each
         gibbs(i,3)=coef2                     ! species.
         gibbs(i,4)=coef3
         gibbs(i,5)=coef4
 11   continue
166   format(a8)
      close(un)
      end

***********************

      subroutine gemabset(temp,absum1,totabund,abund1,first,ktau)

!-----------------------------------------------------------------
!
!     This subroutine sets the initial abundances.
!     It is only called during the first run.
!     Additionally it computes y-bar from WJD.
!
!     25/9-2000 : Jens Falkesgaard,
!                  explicit values changed to double precision
!
!     12/10-2000: GEMABSET tests nOK! abund1(CO) is not set!!
!                 Fixed, JFF
!     13/10-2000: Re-fixed CO error, JFF
!     30/10-2000: Non-atomic part removed.
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 temp                                   ! T in WJD
      real*8 absum1                                 ! y-bar in WJD
      real*8 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*8 totabund(natms)                        ! b in WJD
      real*8 abund1(nspec)                          ! y in WJD
      real*8 compmtx(natms,nspec)                   ! Composition matrix
      logical first                                 ! First time?

!--internal variables

      common/gemcom/gibbs,compmtx
      common /cabink/abink(ndp,nspec)


      do 11 i=1,natms
C???     abund1(i)=abink(ktau,i)
         abund1(i)=totabund(i)        ! Set elemental abundances
 11   continue
      do 12 i=natms+1,nspec
C        abund1(i)=1.0d-40            ! Set molecular and ionic abs.
         abund1(i)=abink(ktau,i)
 12   continue
      absum1=0.0d0
      do 13 i=1,nspec
         absum1=absum1+abund1(i)      ! find y-bar
 13   continue
      first=.false.
      end

***********************

      subroutine gemgibbs(temp,pres,gef,abund1,absum1,pge)

!-----------------------------------------------------------------
!
!     This subroutine computes the partial Gibbs' energies f_i(Y)
!     used in WJD, equations 2 and 15.
!
!     25/9-2000 : Jens Falkesgaard, 
!                  explicit values changed to double precision
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 absum1
      real*8 temp,pres                              ! T and P in WJD
      real*8 gef(nspec)                             ! c in WJD
      real*8 pge(nspec)                             ! partial Gibbs' energies
      real*8 abund1(nspec)                          ! y in WJD
      real*8 gibbs(nspec,5)                         ! Gibbs' coeffs.
      real*8 compmtx(natms,nspec)                   ! composition matrix

!--internal variables

      integer i,j,un
      parameter(un=43)                              ! unit number
      real*8 tt
      real*8 c(nspec)

      common/gemcom/gibbs,compmtx
      
!      print*,'In GEMGIBBS'
!      open(unit=un,file='gemgibbs.ud',status='unknown')
      tt=temp/1.0d3
      do 11 i=1,nspec
         c(i)=0.0d0                                 ! set to zero
 11   continue   
      do 13 i=1,nspec
         do 12 j=5,1,-1
            c(i)=(c(i)*tt)+gibbs(i,j)               ! using JFF Gibbs' fit
 12      continue
         c(i)=1.0d3*c(i)
         c(i)=c(i)/(8.3145d0*temp)                  ! eq. 2 WJD
C  Was(<May2005): c(i)=c(i)+dlog(pres/1.0d5), but with pres in cgs should be:
         c(i)=c(i)+dlog(pres/1.0d6)
         gef(i)=c(i)
         pge(i)=0.0d0                               ! Partial Gibbs' energy
 13   continue
   
      do 14 i=1,nspec                               ! eq. 15 WJD
         abund1(i) = max(abund1(i),1.0d-40)
C        if(abund1(i).lt.0.0d0) abund1(i)=1.0d-40
         pge(i)=abund1(i)*(c(i)+dlog(dabs(abund1(i))/absum1))
!         write(6,*) 'pge,abund1(',i,')= ',pge(i),abund1(i)
!         print*, 'pge(',i,')= ',pge(i)
 14   continue
!      close(un)
      end

***********************

      subroutine gemmatrx(abund1,totabund,pge,eq18a,eq18b)

!-----------------------------------------------------------------
!
!     This subroutine sets up the matrix equation AX=B (equations
!     17 and 18 in WJD) for use in SOLVELIN.
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 abund1(nspec)                          ! y in WJD
      real*8 totabund(natms)                        ! b in WJD
      real*8 pge(nspec)                             ! partial Gibbs' energies
      real*8 eq18a(natms+1,natms+1)                 ! r in WJD (A-matrix)
      real*8 eq18b(natms+1)                         ! B-vector in WJD
      real*8 gibbs(nspec,5)                         ! Gibbs' coefficients
      real*8 compmtx(natms,nspec)                   ! a in WJD

!--internal variables

      common/gemcom/gibbs,compmtx

      do 12 j=1,natms+1                             ! (natms+1)*(natms+1)
         do 11 k=1,natms+1                          ! matrix...
            eq18a(j,k)=0.0d0                        ! initialize to zero
 11      continue
 12   continue   

      do 13 i=1,natms+1                             ! (natms+1) vector...
         eq18b(i)=0.0d0                             ! initialize to zero
 13   continue

      do 16 j=1,natms
         do 15 k=j,natms
            do 14 i=1,nspec
               eq18a(j,k)=eq18a(j,k)+compmtx(j,i)*compmtx(k,i)*abund1(i)
 14         continue                      ! implementation of eq. 17 WJD
            eq18a(k,j)=eq18a(j,k)         ! upper/lower symmetry
 15      continue
 16   continue

      do 17 i=1,natms
         eq18a(i,natms+1)=totabund(i)     ! setting right and lower edge,
         eq18a(natms+1,i)=eq18a(i,natms+1)! lower right element is zero
 17   continue

      do 19 j=1,natms
         do 18 i=1,nspec                  ! setting B-vector elements
            eq18b(j)=eq18b(j)+compmtx(j,i)*pge(i)
 18      continue
 19   continue
      do 20 i=1,nspec                     ! total Gibbs' energy
         eq18b(natms+1)=eq18b(natms+1)+pge(i)
 20   continue
!      print*,'G= ',eq18b(natms+1)
      end

***********************

      subroutine solvelin(a,n,np,indx,b)

!-----------------------------------------------------------------
!
!     This subroutine solves the matrix equation AX=B for X and
!     returns the solution in place of the original B column.
!
!     The routine is based upon 'LUDCMP' and LUBKSB' from 
!     Numerical Recipies in Fortran, 2nd. Ed. 
!     Since this equation is only to be solved once for each 
!     iteration these two routines have been fused into one.
!
!     ref. Numerical Recipies pp. 36-40
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--internal variables

      integer n,np,nmax
      parameter (nmax=100)
      integer i,ii,imax,j,k,ll
      real*8 d,tiny,aamax,dum,sum
      real*8 a(np,np),vv(nmax),b(n)
      real*8 indx(n)

      tiny=1.0d-40
      d=1.0d0             ! LU-decomposition routine from Numerical Recipies
      do 12 i=1,n         ! modified 7/2-2000, JFF; the two routine fused
         aamax=0.0d0      ! 25/9-2000, JFF; tiny changed to double precision
         do 11 j=1,n
            if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11      continue   
         if (aamax.eq.0.0d0) pause 'Singular matrix!'
         vv(i)=1./aamax
 12   continue
      do 19 j=1,n     
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue
            a(i,j)=sum
 14   continue
      aamax=0.0d0
      do 16 i=j,n
         sum=a(i,j)
         do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
 15      continue
         a(i,j)=sum
         dum=vv(i)*dabs(sum)
         if (dum.ge.aamax) then
            imax=i
            aamax=dum
         endif
 16   continue
      if(j.ne.imax) then
         do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
 17      continue
         d=-d
         vv(imax)=vv(j)
      endif
      indx(j)=imax
      if (a(j,j).eq.0.0d0) a(j,j)=tiny
      if (j.ne.n) then
         dum=1./a(j,j)
         do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18      continue
      endif
 19   continue

      ii=0               ! Original beginning of 'LUBKSB'
      do 21 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0) then
            do 20 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 20         continue
         else if (sum.ne.0.0d0) then
            ii=i
         endif
         b(i)=sum
 21   continue
      do 23 i=n,1,-1
         sum=b(i)
         do 22 j=i+1,n
            sum=sum-a(i,j)*b(j)
 22      continue
         b(i)=sum/a(i,i)
 23   continue
      return
      end

***********************

      subroutine gemstep(pge,gef,eq18b,absum1,absum2,abund1,abund2)

!-----------------------------------------------------------------
!
!     This subroutine computes the revised set of abundance 
!     values according to eqs. 14 and 20 in WJD.
!
!     Furthermore it sets abundances that are too low to trace
!     level: 1.0d-40.
!
!     2000: Jens Falkesgaard
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 absum1,absum2                    ! y-bar and x-bar in WJD
      real*8 abund1(nspec),abund2(nspec)      ! y(i) and x(i) in WJD
      real*8 gef(nspec)                       ! actual Gibbs' energies
      real*8 pge(nspec)                       ! partial Gibbs' energies
      real*8 eq18b(natms+1)
      real*8 compmtx(natms,nspec)             ! composition matrix
      real*8 gibbs(nspec,5)                   ! Gibbs' coefficients

!--internal variables

      integer i,j,un
      parameter (un=44)                       ! unit number
      real*8 eq14sum
      real*8 delta(nspec)

      common/gemcom/gibbs,compmtx

!      open(unit=un,file='gemstep.ud',status='unknown')
      absum2=(eq18b(natms+1)+1)*absum1             ! eq. 19 in WJD
      do 12 i=1,nspec
         eq14sum=0.0d0
         do 11 j=1,natms
            eq14sum=eq14sum+eq18b(j)*compmtx(j,i)
 11      continue
         abund2(i)=-pge(i)+(abund1(i)/absum1)*absum2+(eq14sum*abund1(i))
         delta(i)=abund2(i)-abund1(i)              ! how large change?
!         write(un,*) 'delta(',i,')= ',delta(i)
 12   continue                                     ! eq. 14 in WJD
      do 13 i=1,nspec
         if(abund2(i).lt.1.0d-40) abund2(i)=1.0d-40 
                                                   ! Underflow protection
 13   continue
!      close(un)
      end

***********************

      subroutine masscons(abund2,totabund)

!-----------------------------------------------------------------
!
!     This subroutine enforces mass and charge conservation 
!     onto the solution suggested by GEMSTEP.
!
!     16/10-2000: Jens Falkesgaard
!     26/10-2000: Charge conservation incorporated
!     10/01-2001: Updated charge conservation
!
!-----------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'

!--external variables

      real*8 abund2(nspec)                    ! solution suggested by GEMSTEP
      real*8 totabund(natms)                  ! bulk elemental abundances
      real*8 compmtx(natms,nspec)             ! composition matrix
      real*8 gibbs(nspec,5)                   ! Gibbs' coefficients

!--internal variables

      integer i,j,un
      parameter(un=45)                        ! unit number
      real*8 sum(natms)
      real*8 mcdelta(natms)                   ! total mole change
      real*8 qsum,negions,posions,ratio       ! charge checksums

      common/gemcom/gibbs,compmtx

!      print*, 'In MASSCONS'

!      open(unit=un,file='masscons.ud',status='unknown')
      do 11 i=1,natms                            ! check for each element
         sum(i)=0.0d0
         do 10 j=1,nspec                         ! sum up over all species
            sum(i)=sum(i)+abund2(j)*compmtx(i,j) ! remember stoichiometry
 10      continue
         mcdelta(i)=(sum(i)-totabund(i))/totabund(i)
!         write(un,*) 'mcdelta(',i,')= ',mcdelta(i)
 11   continue

      do 13 i=2,natms                            ! don't mass check electrons
         if(dabs(mcdelta(i)).gt.5.0d-2) then 
            ratio=sum(i)/totabund(i)
            do 12 j=2,nspec                      ! if unbalanced...
               if(compmtx(i,j).ne.0.0d0) then    ! ...then adjust
                  if(ratio.gt.1.0d0) abund2(j)=abund2(j)/ratio   
                  if(ratio.lt.1.0d0) abund2(j)=abund2(j)/ratio
               endif
               if(dabs(abund2(j)).lt.1.0d-40) abund2(j)=1.0d-40
 12         continue                             ! underflow protection
         endif
 13   continue

      qsum=1.0d-50
      negions=1.0d-50
      posions=1.0d-50

      do 14 i=1,nspec                     ! check for all species
         if(compmtx(1,i).eq.1.0d0) then
            negions=negions+abund2(i)
         endif
         if(compmtx(1,i).eq.-1.0d0) then
            posions=posions+abund2(i)
         endif
         qsum=qsum+abund2(i)*compmtx(1,i) ! allways electron contribution
 14   continue

      ratio=posions/negions

      if(dabs(qsum).ge.1.0d-30) then            ! if nOK then adjust
         if(ratio.lt.1.0d0) then                ! if positive imbalance:
            do 15 i=1,nspec
               if(compmtx(1,i).eq.1.0d0) abund2(i)=abund2(i)*ratio
 15         continue                            ! adjust
         endif
         if(ratio.gt.1.0d0) then                ! in negative imbalance:
            do 16 i=1,nspec
               if(compmtx(1,i).eq.-1.0d0) abund2(i)=abund2(i)/ratio
 16         continue                            ! ajust
         endif
      endif
!      close(un)
!      print*, 'Exiting MASSCONS'
      end

!-----------------------------------------------------------------
!
!     End of package
!
!-----------------------------------------------------------------

      subroutine molwgt

C this subroutine computes the molecular weight, wmol(2-nspec), for 
C nspec-1 molecules in the array name_gem(nspec), by
C summing the corresponding atomic weights.
C It also idenfifies the atoms that are included in the 
C GEM computation, and associates a GEM-number kgemat(i)
C = 2-49 for the 48 atoms of the i=1-93 (H-U) in Irwins
C partition function table. Correspondingly for the 
C positive ions, kgemion1(i), and the negative ions,
C kgemionm(i) of the corresponding i=1-93 atoms.
C Actually atom 93 is Deuterium, U=92.
C kgem..(i) = 1 coresponds to the the electron (i.e.,
C name_gem(1) = 'e-', pp(k,1) = Pe, etc, and is therefore not 
C used. It is used here to identify the atoms in the vald list which 
C arenot included in the GEM calculations, and which are therefore
C not included in the OS cm2/g* file.


      implicit real*8 (a-h,o-z)
      include 'atomparameter.inc'
      dimension kat(10),fm(10)
      character atomname*2,atomn*2,atom1*1,nameg*8,name_chr*8
      common /cwatom/watom(200)
      common /catomname/atomname(200)
      common /cwmol/wmol(nspec)
      common /cgemat/kgemat(200),kgemion1(200),kgemionm(200)

      character name_gem*8
      common /cgemnames/natms_gem,nions_gem,nspec_gem,name_gem(nspec)
C atms,ions,spec ~ highest index of neutral atoms, ions, species total
       natms_gem = natms
       nions_gem = 127
       nspec_gem = nspec
        open(unit=2,file='/p7/hpg688/hpg688/data/gfits.data'
     &   ,status='old',readonly)   ! extract molecular etc names from here
        do 2125 i=1,nspec
            read(2,766) name_gem(i)
            read(2,*) ajunk1,ajunk2,ajunk3,ajunk4,ajunk5   ! Fit data
2125    continue
766     format(a8)
        close(unit=2)


        call atomw
        call atomnam


C name_gem(1) = e-, and corresponding pp(k,1) = Pe, etc, so if a line
C in the vald inputlist is from an atom not included in the GEM computation,
C we put the index = 1, and skip the inclusion of this line in the 
C OS computation bacause there is no pp and ionization/neutral/molecular
C equilibrium known for it. We could add the Saha equation for the 
C atoms not included in the GEM to be able to compute spectra for such
C atoms (Eu, U, Y, Sc, La, Nd, Sm, Os, Ir, Pt, Au, ...).
        do 410 i=1,200
        kgemion1(i) = 1
        kgemionm(i) = 1
410     kgemat(i) = 1

C           write(6,*)' mo,at,mol_name,f,atomnam,at_wgt,mol_wgt ='
C split the molecule name up into atoms:
        do 2010 i = 2,nspec
          ipos = 1
          wmol(i) = 0.
          nats = 0
          name_chr = name_gem(i)

2011      continue

          if(name_chr(ipos:ipos) .eq. ' ') then
C            write(6,3021) i,name_chr,wmol(i),
C     &         (atomname(kat(nat)),fm(nat),watom(kat(nat)),nat=1,nats)
C3021        format (i3,1x,a8,f8.3,4x, 6(1x,a2,f3.0,f7.2) )
                go to 2010                         !finished this molecule
          end if

          nats = nats + 1

          IF ( ichar(name_chr(ipos:ipos) ) .le.90
     &     .and. ichar(name_chr(ipos:ipos) ) .ge.65 ) THEN   !A-Z~new atom

C 1 or 2 letter atom ?
          if ( ichar(name_chr(ipos+1:ipos+1) ) .le.121
     &     .and. ichar(name_chr(ipos+1:ipos+1) ) .ge.97 ) then   !a-z~2.let atom
                 do 2015 k = 1,93
                   if(name_chr(ipos:ipos+1) .eq. atomname(k)) then
                   katom = k
                   ipos = ipos + 2
                   go to 2016
                   end if
2015             continue
                   write(6,*) 'I didnt find atom in molecule ',i
2016             continue
          else             !i.e., a one-letter atom
                 do 2017 k = 1,93
                   atomn = atomname(k)
                   if(atomn(2:2).eq.' ' .and.
     &                        name_chr(ipos:ipos) .eq. atomn(1:1)) then
                   katom = k
                   ipos = ipos + 1
                   go to 2018
                   end if
2017             continue
                   write(6,*) 'I didnt find atom in molecule ',i
2018             continue
          end if

C we have now identified the atom as atom number katom

          fmult = 1.
          IF ( ichar(name_chr(ipos:ipos) ) .le.57
     &     .and. ichar(name_chr(ipos:ipos) ) .ge.49 ) THEN   !an integer number 2-9

              if ( ichar(name_chr(ipos+1:ipos+1) ) .le.57
     &          .and. ichar(name_chr(ipos+1:ipos+1) ) .ge.48 ) then   !2 digit number 10-99
              fmult = 10.* float( ichar(name_chr(ipos:ipos)) - 48) +
     &              float( ichar(name_chr(ipos+1:ipos+1)) - 48)
              ipos = ipos + 2
             else
              fmult = float( ichar(name_chr(ipos:ipos)) - 48)
              ipos = ipos + 1
             end if

          END IF

          if (name_chr(ipos:ipos).eq.'+' .or.
     &        name_chr(ipos:ipos).eq.'-')   ipos = ipos + 1

            wmol(i) = wmol(i) + fmult * watom(katom)

            fm(nats) = fmult
            kat(nats) = katom

C      attach gem-index kgem..() to atomic number katom:
C         if(i.le.nions_gem .and. atomname(katom).ne.'D ') then
          if(i.le.nions_gem) then
              if (name_chr(ipos-1:ipos-1).eq.'+') then 
                  kgemion1(katom) = i
C            write(6,3119) katom, kgemion1(katom), name_chr(1:ipos)
3119   format('katom,kgemion1(katom),name_chr(1:ipos)=',2i4,2x,a5)
              else if (name_chr(ipos-1:ipos-1).eq.'-') then 
                  kgemionm(katom) = i
              else
                  kgemat(katom) = i
              end if
          end if


         END IF      !end if for new (or first) atom in the molecule

         go to 2011   !more atoms in this molecule ?

2010     continue



       do 3010 iblock=1,nspec,15
       kmm = 15
       kmtp = nspec-iblock+1
       if (kmtp.lt.15) kmm=kmtp
C      write(6,3066) ( name_gem(iblock-1+km),km=1,kmm)
C      write(6,3065) ( wmol(iblock-1+km),km=1,kmm)
3010   continue
3066    format(4x,15a5)
3065    format(15f5.1)
3064    format(1x,15i5)
3067    format(2x,15a5)

       do 3020 iblock=1,natms,15
       kmm = 15
       kmtp = natms-iblock+1
       if (kmtp.lt.15) kmm=kmtp
C       write(6,3064) ( (iblock-1+km),km=1,kmm)
C       write(6,3066) ( name_gem(iblock-1+km),km=1,kmm)
3020   continue

       do 3022 iblock=natms+1,nions,15
       kmm = 15
       kmtp = nions-iblock+1
       if (kmtp.lt.15) kmm=kmtp
C       write(6,3064) ( (iblock-1+km),km=1,kmm)
C       write(6,3066) ( name_gem(iblock-1+km),km=1,kmm)
3022   continue

C...       do 3024 iblock=1,93,15
C...       kmm = 15
C...       kmtp = 93-iblock+1
C...       if (kmtp.lt.15) kmm=kmtp
C...       write(6,3064) ( (iblock-1+km),km=1,kmm)
C...       write(6,3064) ( kgemat(iblock-1+km),km=1,kmm)
C...       write(6,3064) ( kgemion1(iblock-1+km),km=1,kmm)
C...       write(6,3064) ( kgemionm(iblock-1+km),km=1,kmm)
C...       write(6,3067) ( atomname(iblock-1+km),km=1,kmm)
C...       write(6,3066) ( name_gem(kgemat(iblock-1+km)),km=1,kmm)
C...       write(6,3066) ( name_gem(kgemion1(iblock-1+km)),km=1,kmm)
C...       write(6,3066) ( name_gem(kgemionm(iblock-1+km)),km=1,kmm)
C...3024   continue

      return
      END

C---------------------------------------------------------------


      SUBROUTINE atomw
C
C Mean atomic weights are stored into the vector watom, AB2001
C Species 1:92 are included, as in irwin.dat
C Deuterium is added as number 93, because it is treated as a seperate
C atom in gfits.data of molecules. Actually W(H)=1.007825, while Earth
C mixture of H and D has W(H+D)=1.0079 as given for watom(1).
C

      IMPLICIT real*8 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)


      DO i=1,200
         watom(i) = 0.0
      ENDDO

      watom(1)  =   1.0079 !H  Hydrogen
      watom(2)  =   4.0026 !He Helium
      watom(3)  =   6.941  !Li Lithium
      watom(4)  =   9.0121 !Be Beryllium
      watom(5)  =  10.81   !B  Boron
      watom(6)  =  12.011  !C  Carbon
      watom(7)  =  14.0067 !N  Nitrogen
      watom(8)  =  15.9994 !O  Oxygen
      watom(9)  =  18.9984 !F  Fluorine
      watom(10) =  20.179  !Ne Neon
      watom(11) =  22.9897 !Na Sodium
      watom(12) =  24.305  !Mg Magnesium
      watom(13) =  26.9814 !Al Aluminum
      watom(14) =  28.0855 !Si Silicon
      watom(15) =  30.9737 !P  Phosphorus
      watom(16) =  32.06   !S  Sulfur
      watom(17) =  35.453  !Cl Chlorine
      watom(18) =  39.948  !Ar Argon
      watom(19) =  39.0983 !K  Potassium
      watom(20) =  40.08   !Ca Calcium
      watom(21) =  44.9559 !Sc Scandium
      watom(22) =  47.88   !Ti Titanium
      watom(23) =  50.9415 !V  Vanadium
      watom(24) =  51.996  !Cr Chromium
      watom(25) =  54.9380 !Mn Manganese
      watom(26) =  55.847  !Fe Iron
      watom(27) =  58.9332 !Co Cobalt
      watom(28) =  58.96   !Ni Nickel
      watom(29) =  63.546  !Cu Copper
      watom(30) =  65.38   !Zn Zinc
      watom(31) =  69.72   !Ga Gallium
      watom(32) =  72.59   !Ge Germanium
      watom(33) =  74.9216 !As Arsenic
      watom(34) =  78.96   !Se Selenium
      watom(35) =  79.904  !Br Bromine
      watom(36) =  83.80   !Kr Krypton
      watom(37) =  85.4678 !Rb Rubidium
      watom(38) =  87.62   !Sr Strontium
      watom(39) =  88.9059 !Y  Yttrium
      watom(40) =  91.22   !Zr Zirconium
      watom(41) =  92.9064 !Nb Niobium
      watom(42) =  95.94   !Mo Molybdenum
      watom(43) =  97.907  !Tc Technetium
      watom(44) = 101.07   !Ru Ruthenium
      watom(45) = 102.9055 !Rh Rhodium
      watom(46) = 106.42   !Pd Palladium
      watom(47) = 107.868  !Ag Silver
      watom(48) = 112.41   !Cd Cadmium
      watom(49) = 114.82   !In Indium
      watom(50) = 118.69   !Sn Tin
      watom(51) = 121.75   !Sb Antimony
      watom(52) = 127.60   !Te Tellurium
      watom(53) = 126.9045 !I  Iodine
      watom(54) = 131.29   !Xe Xenon
      watom(55) = 132.9054 !Cs Cesium
      watom(56) = 137.33   !Ba Barium
      watom(57) = 138.9055 !La Lanthanum
      watom(58) = 140.12   !Ce Cerium
      watom(59) = 140.9077 !Pr Praseodymium
      watom(60) = 144.24   !Nd Neodymium
      watom(61) = 144.913  !Pm Promethium
      watom(62) = 150.36   !Sm Samarium
      watom(63) = 151.96   !Eu Europium
      watom(64) = 157.25   !Gd Gadolinium
      watom(65) = 158.9254 !Tb Terbium
      watom(66) = 162.50   !Dy Dysprosium
      watom(67) = 164.9304 !Ho Holmium
      watom(68) = 167.26   !Er Erbium
      watom(69) = 168.9342 !Tm Thulium
      watom(70) = 173.04   !Yb Ytterbium
      watom(71) = 174.967  !Lu Lutetium
      watom(72) = 178.49   !Hf Hafnium
      watom(73) = 180.9479 !Ta Tantalum
      watom(74) = 183.85   !W  Tungsten
      watom(75) = 186.207  !Re Rhenium
      watom(76) = 190.2    !Os Osmium
      watom(77) = 192.22   !Ir Iridium
      watom(78) = 195.08   !Pt Platinum
      watom(79) = 196.9665 !Au Gold
      watom(80) = 200.59   !Hg Mercury
      watom(81) = 204.383  !Tl Thallium
      watom(82) = 207.2    !Pb Lead
      watom(83) = 208.9804 !Bi Bismuth
      watom(84) = 208.982  !Po Polonium
      watom(85) = 209.987  !At Astatine
      watom(86) = 222.018  !Rn Radon
      watom(87) = 223.020  !Fr Francium
      watom(88) = 226.0254 !Ra Radium
      watom(89) = 227.0278 !Ac Actinium
      watom(90) = 232.0381 !Th Thorium
      watom(91) = 231.0359 !Pa Protactinium
      watom(92) = 238.051  !U  Uranium
      watom(93) = 2.01410  !D  Deuterium

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE atomnam
C
C 1 or 2 character names of atoms.
C Species 1:92 are included, as in irwin.dat
C


      implicit REAL*8 (a-h,o-z)
      character atomname*2
      common /cwatom/watom(200)
      common /catomname/atomname(200)

      DO i=1,200
         atomname(i) = '  '
      ENDDO

      atomname(1)  = 'H ' ! Hydrogen
      atomname(2)  = 'He' ! Helium
      atomname(3)  = 'Li' ! Lithium
      atomname(4)  = 'Be' ! Beryllium
      atomname(5)  = 'B ' ! Boron
      atomname(6)  = 'C ' ! Carbon
      atomname(7)  = 'N ' ! Nitrogen
      atomname(8)  = 'O ' ! Oxygen
      atomname(9)  = 'F ' ! Fluorine
      atomname(10) = 'Ne' ! Neon
      atomname(11) = 'Na' ! Sodium
      atomname(12) = 'Mg' ! Magnesium
      atomname(13) = 'Al' ! Aluminum
      atomname(14) = 'Si' ! Silicon
      atomname(15) = 'P ' ! Phosphorus
      atomname(16) = 'S ' ! Sulfur
      atomname(17) = 'Cl' ! Chlorine
      atomname(18) = 'Ar' ! Argon
      atomname(19) = 'K ' ! Potassium
      atomname(20) = 'Ca' ! Calcium
      atomname(21) = 'Sc' ! Scandium
      atomname(22) = 'Ti' ! Titanium
      atomname(23) = 'V ' ! Vanadium
      atomname(24) = 'Cr' ! Chromium
      atomname(25) = 'Mn' ! Manganese
      atomname(26) = 'Fe' ! Iron
      atomname(27) = 'Co' ! Cobalt
      atomname(28) = 'Ni' ! Nickel
      atomname(29) = 'Cu' ! Copper
      atomname(30) = 'Zn' ! Zinc
      atomname(31) = 'Ga' ! Gallium
      atomname(32) = 'Ge' ! Germanium
      atomname(33) = 'As' ! Arsenic
      atomname(34) = 'Se' ! Selenium
      atomname(35) = 'Br' ! Bromine
      atomname(36) = 'Kr' ! Krypton
      atomname(37) = 'Rb' ! Rubidium
      atomname(38) = 'Sr' ! Strontium
      atomname(39) = 'Y ' ! Yttrium
      atomname(40) = 'Zr' ! Zirconium
      atomname(41) = 'Nb' ! Niobium
      atomname(42) = 'Mo' ! Molybdenum
      atomname(43) = 'Tc' ! Technetium
      atomname(44) = 'Ru' ! Ruthenium
      atomname(45) = 'Rh' ! Rhodium
      atomname(46) = 'Pd' ! Palladium
      atomname(47) = 'Ag' ! Silver
      atomname(48) = 'Cd' ! Cadmium
      atomname(49) = 'In' ! Indium
      atomname(50) = 'Sn' ! Tin
      atomname(51) = 'Sb' ! Antimony
      atomname(52) = 'Te' ! Tellurium
      atomname(53) = 'I ' ! Iodine
      atomname(54) = 'Xe' ! Xenon
      atomname(55) = 'Cs' ! Cesium
      atomname(56) = 'Ba' ! Barium
      atomname(57) = 'La' ! Lanthanum
      atomname(58) = 'Ce' ! Cerium
      atomname(59) = 'Pr' ! Praseodymium
      atomname(60) = 'Nd' ! Neodymium
      atomname(61) = 'Pm' ! Promethium
      atomname(62) = 'Sm' ! Samarium
      atomname(63) = 'Eu' ! Europium
      atomname(64) = 'Gd' ! Gadolinium
      atomname(65) = 'Tb' ! Terbium
      atomname(66) = 'Dy' ! Dysprosium
      atomname(67) = 'Ho' ! Holmium
      atomname(68) = 'Er' ! Erbium
      atomname(69) = 'Tm' ! Thulium
      atomname(70) = 'Yb' ! Ytterbium
      atomname(71) = 'Lu' ! Lutetium
      atomname(72) = 'Hf' ! Hafnium
      atomname(73) = 'Ta' ! Tantalum
      atomname(74) = 'W ' ! Tungsten
      atomname(75) = 'Re' ! Rhenium
      atomname(76) = 'Os' ! Osmium
      atomname(77) = 'Ir' ! Iridium
      atomname(78) = 'Pt' ! Platinum
      atomname(79) = 'Au' ! Gold
      atomname(80) = 'Hg' ! Mercury
      atomname(81) = 'Tl' ! Thallium
      atomname(82) = 'Pb' ! Lead
      atomname(83) = 'Bi' ! Bismuth
      atomname(84) = 'Po' ! Polonium
      atomname(85) = 'At' ! Astatine
      atomname(86) = 'Rn' ! Radon
      atomname(87) = 'Fr' ! Francium
      atomname(88) = 'Ra' ! Radium
      atomname(89) = 'Ac' ! Actinium
      atomname(90) = 'Th' ! Thorium
      atomname(91) = 'Pa' ! Protactinium
      atomname(92) = 'U ' ! Uranium
      atomname(93) = 'D ' ! Deuterium

      RETURN
      END

