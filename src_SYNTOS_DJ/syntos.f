!-----------------------------------------------------------------------
! SYNTOS: Calculates the synthetic spectrum, the Planck function, the
! continuum absorption coefficient, the molecular absorption etc. 
!-----------------------------------------------------------------------
      program syntos
       
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      character model_file*45,spectrum_file*60,spect_head*60
      dimension flk(mmol),wlmol(mopav),opimol(mmol,mopav),
     *  ispec1(0:mmol),ispec2(0:mmol),ispec3(0:mmol),
     *  ispec4(0:mmol),ispec5(0:mmol),ispec6(0:mmol),
     *  ispec7(0:mmol),ispec8(0:mmol),ispec9(0:mmol)
      common /cmodel/teff,g,abund(16),
     *  taumod(ndp),t(ndp),pe(ndp),pg(ndp),prad(ndp) ,
     *  pturb(ndp),xkapr(ndp),ro(ndp),emu(ndp),
     *  xl(20),abska(20,ndp),sprida(20,ndp),ntau,nlp,model_file
      common /csplim/vspmin,vspmax,ilist
      common /ctran/tau(ndp),x(ndp),bplan(ndp),hsurf,y1(6),jtau
      common/cfold/flux(0:mopac),fluxinst(mopac),wlinst(0:mopac),
     *  fluxinstd(mopac)
      common /cfilter/flux_bpw(0:mopac),amag(25),nfilt,neso,ifilt(25),
     *  iscan
      common /cmolincl/imol1,imol2,ispec(0:mmol)
      common /tg01ba/i1,in,kki 
      common /cmolabs/ kmol(0:mmol)
      common /cnu/capnu(ndp),capcon(ndp),taunu(ndp)
      common /cwnos/abder(20,ndp),wnos(mopac),nwnossp,nwmin,nwmax,nosav
      common /cwnl/wnl(mopac),nwnl,nwlmin,nwlmax
      common/cos/wlos(mopac),wlosstep(mopac),osresl,ktemp
      common /comdata/ BDA,pi,c_vel
      common/cppinp/ pressave(ndp,99),nsave,namesave(99)
      
      namelist /input/ ilist,vspmin,vspmax,osresl,
     * model_file,spectrum_file,spect_head,
     * lavwr,
     * lcolorsonly,lfilter4all,
     * nfilt,neso,ifilt,iscan,
     * inst_cal,nhalf_inst,resl_inst,nosav
      namelist /input_mol/imol1,imol2,nspect,
     * ispec1,ispec2,ispec3,ispec4,ispec5,ispec6,ispec7,ispec8,ispec9

      data bda,pi,c_vel /8.31327e7, 3.141593, 2.99792e10/
      data i1,in/1,1/  ! extrapolation method

! Read input file
      open(unit=1,file='syntos.input',status='old',readonly)
      read(1,input)
      read(1,input_mol)
      if(imol2 .gt. mmol) stop 'Increase dimension of MMOL'
      if(nspect .gt. mmol) stop 'Increase dimension of FLK'
      
! Read model
      call model
      if(ntau.gt.ndp) stop 'Increase dimension of NDP'

! Colors and/or spectrum
      if(lcolorsonly.eq.1) then
        open(unit=16,file='spectrumfil.add',status='unknown')
      else
        open(unit=16,file=spectrum_file,status='unknown')
      end if

! Calculate continuum abs+ska and derivaties
      call continuum_deriva

! Loop over spectra
      do 143 kspec = 1,nspect+1
        open(unit=12,status='unknown')

        if(kspec.eq.1) then		! Compute continuum only
          ispec(0:mmol) = 0
        else if(kspec.eq.2) then
          ispec(imol1:imol2) = ispec1(imol1:imol2)
        else if(kspec.eq.3) then
          ispec(imol1:imol2) = ispec2(imol1:imol2)
        else if(kspec.eq.4) then
          ispec(imol1:imol2) = ispec3(imol1:imol2)
        else if(kspec.eq.5) then
          ispec(imol1:imol2) = ispec4(imol1:imol2)
        else if(kspec.eq.6) then
          ispec(imol1:imol2) = ispec5(imol1:imol2)
        else if(kspec.eq.7) then
          ispec(imol1:imol2) = ispec6(imol1:imol2)
        else if(kspec.eq.8) then
          ispec(imol1:imol2) = ispec7(imol1:imol2)
        else if(kspec.eq.9) then
          ispec(imol1:imol2) = ispec8(imol1:imol2)
        else if(kspec.eq.10) then
          ispec(imol1:imol2) = ispec9(imol1:imol2)
        end if

! Count the molecules to be included in the spectrum
        nmol = sum(ispec(imol1:imol2))
        molmax = max(nmol,0)
          
        write(6,'(a36,4i3)') 'kspec,imol1,imol2,jm,molmax: ',
     *    kspec,imol1,imol2,nmol
        write(6,'(a8,30i2)') 'ispec = ',
     *    (ispec(imol),imol=imol1,imol2)

! Read absorption coefficients, interpolate, sum up
        call osopacity(kspec)
        close(1)

! Compute the radiative transfer and the spectrum
        do 120 jv=nwmin,nwmax            ! there is NV+1 points
	  if (jv.eq.0) go to 120
          jvm=jv
          call kappaline(jvm)
	  call taunuscale(jvm)
	  call transf(jvm) 

! PLanck function: [B_w] = [erg cm-2 s-1 (cm-1)-1], w: [cm-1]
! The Planck function is intensity; if one wants flux F=pi*intensity
          hsurflam = hsurf*wnos(jv)**2       !flux for wavelength (per cm)
          hsurflam = hsurflam*1.d-8          !flux f.w. in per Angstrom
          hsurfmu  = hsurflam*1.d4           !flux f.w. in per Micron
          flux(jvm)=hsurf            	     !....else F-nu [erg cm-2 s-1 (cm-1)-1]

          flux_bpw(jvm) =
     *      1.191d-5*wnos(jv)**3/(exp(1.439*wnos(jv)/teff)-1.)
          planckmu = flux_bpw(jvm)*wnos(jv)**2*1.d-4/pi


          write(12,'(f14.7,1p30e12.4)') 1.e4/wnos(jvm),
     *      hsurfmu, planckmu
120     continue

! Compute and write average flux to unit 12
        if(lavwr.gt.1) then
          open(unit=88,status='unknown')
          rewind(unit=12)
          avwr = dfloat(lavwr)
          navwr = 0
          do i=1,mopac,lavwr
	    wl = 0.
            fl = 0.
	    fpl = 0.
	    wn = 0.
	    fn = 0.
	    fpn = 0.
	    do ia=1,lavwr
              read(12,*,end=179) wll,fll,fpll
	      fpl = fpl + fpll
	      wl = wl + wll
              fl = fl + fll
            end do
            wl = wl/avwr
            fl = fl/avwr
            fpl = fpl/avwr
            navwr = navwr + 1
            write(88,2021) wl,fl,fpl
          end do
179       continue
          close(unit=12)
       
          open(unit=12,status='unknown')
          rewind(88)
          do i=1,navwr
            read(88,*) wl,fl,fpl
            write(12,2021) wl,fl,fpl
          end do
          close(unit=88)
        end if
          
2021    format(f14.7,1p30e12.4)
2022    format(f14.7,1p2e12.4,0pf17.7,1p2e12.4)


! Unit 12: this kspec, unit 16: previous kspec-1, we now add unit 12 to unit 16
        rewind(12)
        rewind(16)
        nwlfil = 0

        if(kspec.eq.1) then
          write(16,'(a60)') model_file
          write(16,'(a60)') spectrum_file
          write(16,*) ' no atoms included in spectrum'
          nwlsp = nwmax-nwmin+1
          write(16,1692) osresl,nwlsp
1692      format('Resolution, number of wavelengths =',f8.0,i7)
          write(16,*) ' continuum only'
          write(16,1693) molmax
1693      format('Max numb molecules:',i3)
          write(16,*) ' '
          write(16,'(a60)') spect_head
          ncol = kspec+2
          write(16,'(i3,i7)') ncol,nwlsp
          do i = 1, mopac
            read(12,*,end=1441) wl,fl,fpl
            write(16,2021) wl,fl,fpl
            nwlfil = nwlfil + 1
          end do
1441      continue
          nwlsp = nwlfil
        else
          open(unit=13,status='unknown')
          do iad = 1,9
            read(16,*)
          end do
          do i = 1, mopac
            read(16,*,end=1443) wl,(flk(kcol),kcol=1,kspec-1)
            write(13,2021) wl,(flk(kcol),kcol=1,kspec-1)
            if(i.eq.1) wlfirst = wl
            nwlfil = nwlfil + 1
          end do
1443      continue
          wllast  = wl
          close(13)
          open(unit=13,status='unknown')
          close(16)
          if(lcolorsonly.eq.1) then
            open(unit=16,file='spectrumfil.add',status='unknown')
          else
            open(unit=16,file=spectrum_file,status='unknown')
          end if
          nwlfil = 0
          write(16,'(a60)') model_file
          write(16,'(a60)') spectrum_file
            write(16,*) ' no atoms included in spectrum'
          write(16,1692) osresl,nwlsp
          write(16,1694) wlfirst,wllast,1.e4/wllast,1.e4/wlfirst
1694      format('First/last wavelength in mu:',2f10.7,
     *        '=in cm-1:',2f9.2)
          write(16,1693) molmax
          write(16,*) ' '
          write(16,'(a60)') spect_head
          ncol = kspec+2
          write(16,'(i3,i7)') ncol,nwlsp
          do i = 1, mopac
            read(12,*,end=1445) wl,fl,fpl
            read(13,*) wl,(flk(kcol),kcol=1,kspec-1)
            write(16,2021) wl,(flk(kcol),kcol=1,kspec-1),fl,fpl
            nwlfil = nwlfil + 1
          end do
1445      continue
          close(13)
        end if
        close(unit=12)

! For all spectra, or only complete spectra, compute the filter magnitudes:
        if((kspec.eq.nspect+1 .or. lfilter4all.eq.1) 
     *                        .and. nfilt.GT.0) then
          call filter
          write(6,*)'filter magnitudes computed; nfilt=',nfilt
        end if

143   continue

      close(2)
      close(23)
      close(16)

      stop
      end
      
!------------------------------------------------------------
! MODEL: Read in MARCS model atmosphere
!
! TEFF  : EFFECTIVE TEMPERATURE
! G     : ACCELERATION OF GRAVITY
! ABUND : ABUNDANCE OF H,HE,C,N,O,...... (NOT LOG() !)
! NTAU  : NUMBER OF DEPTH POINTS
! TAU   : STANDARD OPTICAL DEPTH
! T       TEMPERATURES
! PE    : ELECTRON PRESSURES
! PG    : GAS PRESSURES
! PREAD : RADIATION PRESSURES
! PTURB : TURBULENT PRESSURES
! XKPAPR: STANDARD ABSORPTION COEFFICIENTS
! RO    : DENSITIES
! EMU   : MEAN MOLECULAR WEIGHT (NUCLEONS PER PARTICLE)
! PRESSAVE: MOLECULAR PARTIAL PRESSURES (C2H2=15, HCN=16)
!       : (IN DYN/CM2) - NOT LOG()
! XL    : WAVELENGTHS FOR ABSKA, SPRIDA (ANGSTROEM)
! ABSKA : CONTINUOUS ABSORPTION COEFFICIENTS
! SPRIDA: CONTINUOUS SCATTERING COEFFICIENTS
!------------------------------------------------------------
      subroutine model

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      real*8 junk1,junk2,junk3,junk4,junk5
      character iden*15,id2*115,model_file*45,namesave*8
      character*8 name_listmo(16),name_ex(250),myline*100,
     *  name_comp1,name_comp2
      dimension pp(ndp,99)
      common /cmolincl/imol1,imol2,ispec(0:mmol)
      common /cmodel/teff,g,abund(16),
     * taumod(ndp),t(ndp),pe(ndp),pg(ndp),prad(ndp),
     * pturb(ndp),xkapr(ndp),ro(ndp),emu(ndp),
     * xl(20),abska(20,ndp),sprida(20,ndp),ntau,nlp,model_file
      common /ctran/tau(ndp),x(ndp),bplan(ndp),hsurf,y1(6),jtau
      common/cppinp/ pressave(ndp,99),nsave,namesave(99)
      common/cmxread/pel_ion(ndp,16,4), elabund(ndp,16),pprel(ndp,nspec)        !elabund=#atoms_i/atom_in_gas
      common /cmarcsgem/ niden
      common /cabink/abink(ndp,nspec)

! Zeroset
      pprel(1:ndp,1:nspec) = 1.d-40
      pel_ion(1:ndp,1:16,1:4) = 1.0d-30
      name_listmo(1:16) = '        '

! Open model file
      open (unit=22,file=model_file,status='old',readonly)



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
        WRITE(6,10) model_file,TEFF,LOG10(G),10.**(ABUND(3)-ABUND(5))
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
        END IF
41053   FORMAT('H,He,C,N,O,Fe:'8f8.3)
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
        DO 1410 I=1,NTAU
        READ(22,*) K1,TAUMOD(I),TAUS,Z,T(I),PE(I)
     &  ,PG(I),PRAD(I),PTURB(I),XKAPR(I)
        TAU(I)=TAUMOD(I)
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
        DO 4112 I=1,NTAU
        READ(22,*) IK1,TR,RO(I),EMU(I),CP,CV,AGRAD,Q,U,V,FCF,IK2
4112    CONTINUE
4102    CONTINUE

        IPOS = INDEX(IDEN,'L O G A R I T')
        IF ((IPOS.EQ.0).OR.(LOGKILL.EQ.1)) GO TO 4107
        IFLAG=IFLAG+1  !IFLAG=5
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
        IF (ILOOP.GE.100) STOP ' ***Error: I found no K in PP'
        IF (IPOS.EQ.0) GO TO 5831
           backspace(22)
           read(22,2143) (name_listmo(i),i=1,16)
2143  FORMAT(5x,16a8)
           innam = 0
           do 2161 i=1,16
2161       if (name_listmo(i).ne.'        ') innam = innam+1
           niden = nmol

! Read mol names
           do 2121 j=1,innam
           name_comp2 = adjustl(name_listmo(j))
           if(name_comp2.eq.'P(H)    ') name_comp2='H       '
           if(name_comp2.eq.'K       ') name_comp2='XXXXXXXX'
           nsave = nsave+1
           namesave(nsave) = name_comp2
2121    CONTINUE

! Read pp
        DO 2130 I=1,NTAU
        READ(22,*) K1,(pp(I,K),K=1,innam)
        do 2140 k=1,innam
        js=nsave-innam+k
2140    PRESSAVE(I,js) = 10.d0**pp(I,K)
2130    CONTINUE
2142    format(2x,8a8/,8a8)
2146    format(8f8.2/,8f8.2)

        if (niden-nmol.eq.0) then
          go to 2135
        end if
        NMOL = niden
        go to 2135

4106    CONTINUE
2139    CONTINUE
2136    continue
4107    CONTINUE

        IF(IDEN.NE.' LAMBDA FOR ABS') GO TO 4103
        IFLAG=IFLAG+1            !IFLAG=6
        READ(22,*) NLP                      !usually NLP=20 in input, #20 is 8000.AA
        NLP=NLP-1
        READ(22,*) (XL(J),J=1,NLP),XLEXTRA   !cont.wavelength scale is 19 freq.
        READ(22,11) IDEN,ID2
        READ(22,11) IDEN,ID2
        DO 4202 I=1,NTAU
        READ(22,*) ITAU,TAUX
        READ (22,*) (ABSKA(K,I),K=1,NLP),ABSEXTRA
        READ (22,*) (SPRIDA(K,I),K=1,NLP),SPREXTRA
312   format(1p8e10.3)
4202    CONTINUE
4103    CONTINUE
C
4100     CONTINUE
4199     CONTINUE
C
      CLOSE (22)
      CLOSE (2)

999     CONTINUE
      WRITE(6,*) ' MODEL STOP: HOW COULD THIS BE ???????'
998     CONTINUE


      RETURN
      END

!-----------------------------------------------------------------------
! CONTINUUM_DERIVA: Calculates the continuum absorption coefficient by 
! use of log-log spline interpolation in the absorption coefficients ABSKA 
! which are read in from the model atmosphere.             
! ABSKA(J,K) = log(absorption+scattering) in wavelength point J, layer K      
! ABDER(J,K) = derivatives of ABSKA as function of log_10(wavelength)                                                           
!-----------------------------------------------------------------------
      subroutine continuum_deriva
      
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      dimension abskk(20),abderk(20),w(3,20)  !dimension of W is (3,NPL+1)
      character model_file*45
      common /csplim/vspmin,vspmax,ilist
      common /cmodel/teff,g,abund(16),
     *  taumod(ndp),t(ndp),pe(ndp),pg(ndp),prad(ndp),
     *  pturb(ndp),xkapr(ndp),ro(ndp),emu(ndp),
     *  xl(20),abska(20,ndp),sprida(20,ndp),ntau,nlp,model_file
      common /ctran/tau(ndp),x(ndp),bplan(ndp),hsurf,y1(6),jtau
      common /cwnos/abder(20,ndp),wnos(mopac),nwnossp,nwmin,nwmax,nosav
      common /cwnl/wnl(mopac),nwnl,nwlmin,nwlmax

      if (nlp.ne.19) then
        write(6,*) ' nlp = ',nlp,' in subr continuum_deriva'
        stop ' nlp not 19'
      end if
      xl(1:nlp)=log10(xl(1:nlp))    ! XL = wavelength in AA from mod.atm.
      do k=1,ntau
        do j=1,nlp
	  abska(j,k)=log(abska(j,k)+sprida(j,k))
	  abskk(j)=abska(j,k)
        end do
        
! Calculate the derivative ABDERK in all NLP nodes when XL and ABSKK are given
	call tb04at(nlp+1,nlp,xl,abskk,abderk,w)
        abder(1:nlp,k)=abderk(1:nlp)
      end do

      return 
      end

!-----------------------------------------------------------------------
! OSOPACITY: Reads the OS-files for all molecules considered. 
! Creates wavenumber array, WNOS(1-NWNOS), for calculation of opacity.
! Molecular opacity stored in OPACITY(MOPAC,NDP).
! Continuum opacity stored in CONTOPAC(MOPAC,NDP)
!-----------------------------------------------------------------------
       subroutine osopacity(kspec)

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      character osfil*55,molid*4,model_file*45
      character osname*4
      character*8 namesave,name_comp
      character nameid*8
      logical first
      dimension abderk(20),abskk(20),osname(0:mmol),sumop(0:mmol,ndp)
      dimension wnb(25),wnstep(25),tmol(mtemp),kmol(0:mmol),sumop2(ndp)
      dimension ac1(mtemp),ac1_pr(mtemp),osfil(0:mmol),reliso(5),
     *   opimol(ndp),ac1ln(mtemp),ac2(ndp),wt(3,mtemp),dadt(mtemp),
     *   xwlg10(mopac),sumop3(ndp)
      dimension pp_sum(ndp),ropp(ndp)
      dimension absk2(200),abder1(ndp,200),abder2(200),w2(3,200)
      dimension pel(ndp), rho(ndp)
      common /cmodel/teff,g,abund(16),
     *  taumod(ndp),tatmos(ndp),pe(ndp),pg(ndp),prad(ndp),
     *  pturb(ndp),xkapr(ndp),ro(ndp),emu(ndp),
     *  xl(20),abska(20,ndp),sprida(20,ndp),ntau,nlp,model_file
      common/cppinp/ pressave(ndp,99),nsave,namesave(99)
      common /tg01ba/i1,in,kki 
      common /cwmol/wmol(nspec)
      common /cabink/abink(ndp,nspec)
      common /cmarcsgem/ niden
      common /copacity/contopac(mopac,ndp),opacity(mopac,ndp)
      common /csplim/vspmin,vspmax,ilist
      common /cmolincl/imol1,imol2,ispec(0:mmol)
      common /cwnos/abder(20,ndp),wnos(mopac),nwnossp,nwmin,nwmax,nosav
      common/cos/wlos(mopac),wlosstep(mopac),osresl,ktemp
      common /kapcon/ absk(ndp,100), sprid(ndp,100)

      data rmol /8.3143d7/    ! Molar gas constant R in erg/(mol*K)
      data first/.true./

      namelist /inos/ osname,osfil
      namelist /inputosmol/ molid, ktemp, tmol, nwnos,
     *  vkms, kiso, reliso, l_per_stellar
     

!  Calculate the wave numbers, WN, for OS opacity and the total number, NWN,
!  of OS frequency points used.
      if(kspec.gt.1) goto 339
      call oswavenumbers

! NWNOSSP is the number of frequencies in the spectrum (NWNOS
! is used for input for individual molecules from OS files)
      nwmin = 1
      nwmax = nwnossp
!      nwmin = 1
!      nwmax = 1
!      print *, nwnossp
!      print *, wnos(1:nwnossp)
!      print *, vspmin, vspmax
!      do iw = 1,nwnossp
!        if(wnos(iw).le.vspmin) then 
!          nwmin = iw
!          cycle
!        end if
!        if(wnos(iw).gt.vspmax) then
!         nwmax = iw
!         exit
!        end if
!      end do
!      print *, nwmin
!      print *, nwmax
!      stop

! Read molecule file names
      read(1,inos)

339   continue

! Zeroset
      sumatdw = 0.
      contopac(1:mopac,1:ndp) = 0.
      opacity(1:mopac,1:ndp) = 0.

! Continuum opacity in [cm^2/g*], interpole MARCS abska/sprida
        if (nlp.ne.19) then
          write(6,*) ' in OSOPACITY:'
          write(6,*) ' nlp = ',nlp
          write(6,*) teff,g,ntau
          stop ' NLP promlem: 19 versus 20?'
        end if
        xwlg10(nwmin:nwmax)=8.-log10(wnos(nwmin:nwmax))

        do k=1,ntau
          abderk(1:nlp)=abder(1:nlp,k)
          abskk(1:nlp)=abska(1:nlp,k)   !at this place abska = log10(abska+sprida)     
          do jv = nwmin,nwmax
            jv1 = min(nwmax,jv+1)
            dw = wnos(jv1)-wnos(jv)
            x1=xwlg10(jv)
C abskk(j) is log10 of the continuums abs.coef. (from MARSC) at the 
C NLP=19 wavelengths it is calculated at. 
C If the wavelength x1(jv) (in aangstrom) is outside interval,
C extrapolation is determined by I1 AND IN in tg01bt, but beware that we
C inter/extrapolate in log10 here and extrapolation values are set to 
C 0 if outside interval if i1=in=0 and hence exp(f(x))=1
            opacity(jv,k) = exp(tg01bt(-1,nlp+1,nlp,xl,abskk,abderk,x1)) !NLP=19
            if(x1.le.xl(1)) then
              if(i1.eq.0) opacity(jv,k) = 0.
            else if(x1.ge.xl(nlp)) then
              if(in.eq.0) opacity(jv,k) = 0.
            end if
            if(opacity(jv,k) .lt.0.d0) then
              write(6,*) ' ***hm:jv,k,wl_mu,contop:',
     *          jv,k,1.e4/wnos(jv),opacity(jv,k)
            end if
          end do
        end do
      contopac(nwmin:nwmax,1:ntau) = opacity(nwmin:nwmax,1:ntau) 
      sumop3(1:ndp) = 0.
      do jv=nwmin+1,nwmax-1
         dw=(wnos(jv+1)-wnos(jv-1))/2.
         do it=1,ntau
              sumop3(it) = sumop3(it) + opacity(jv,it)*dw
         end do
      end do


! Molecules/atoms
      if(kspec.eq.1) go to 199
      kmol(0:mmol) = -1
      do i=0,24
        if(osname(i)(1:4).eq.'at1d' .or. osname(i)(1:4).eq.'dabs' .or.
     *     osname(i)(1:4).eq.'dsca') then 
          kmol(i) = 0
          cycle
        end if
        do j=1,nsave
          if(osname(i)(1:4) .eq. namesave(j)(1:4)) then
             kmol(i) = j
             exit
          end if
        end do
        if(kmol(i).lt.0) then
          stop 'mol not identified'
        end if
      end do

      sumop(0:mmol,1:ndp) = 0.
      sumop2(1:ndp) = 0.
      
      do 123 imol=imol1,imol2


        if(ispec(imol).eq.0) cycle
        
        l_per_stellar = 0
        open(unit=4,file=osfil(imol),status='old',readonly)
        read(4,inputosmol)

        if(ktemp.gt.mtemp) stop 'Increase dimension of MTEMP'
        write(6,2031) molid,imol,osfil(imol)
2031    format (1x,a4,' (imol=',i2,')',a55)

        avdw = 0.
        nav = 0
        nimol = 0
        jv = nwmin

! Initiate reading such that wn<=wnos(jv_first):
        read(4,*,end=499) wn,(ac1(it),it=1,ktemp)
        wn1=wn

        if(wn.gt.wnos(nwmax)) then
          write(6,1278) wn,nwmax,wnos(nwmax)
1278    format('(data not used: first wn:',f9.2,' nwmax,wnos(nwmax) =',
     *    i5,f9.2,')')
          close(unit=4)
          go to 123
        end if

271     continue
        if(wn.gt.wnos(jv)) then
          jv = jv + 1
          go to 271
        end if

        do 400 i=2,nwnos
          wn_pr = wn
          ac1_pr(1:ktemp)= ac1(1:ktemp)
          read(4,*,end=499) wn,(ac1(it),it=1,ktemp)
          wn2 = wn
          if(wn.lt.wnos(1)) go to 400
          if(wn.gt.wnos(nwnossp)) go to 499
          nimol = nimol + 1
          if(nimol.eq.1) then
            wimol_1 = wn
          end if
          wimol_2 = wn
290       continue
  
! begin for wn_pr < wnos(jv) < wn   ...|
          if(wn.ge.wnos(jv) .and. wn_pr.le.wnos(jv)) then
            if(wn-wnos(jv) .le. wnos(jv)-wn_pr) then
              ac1ln(1:ktemp) = log(max(ac1(1:ktemp),1.d-31))
            else
              ac1ln(1:ktemp) = log(max(ac1_pr(1:ktemp),1.d-31))
            end if
          
! The derivative, DADT(it) at each temperature T(it)
            call tb04at(mtemp,ktemp,tmol,ac1ln,dadt,wt)
            
            do 283 it = 1,ntau
              tx = tatmos(it)
              ac2(it) = exp( tg01bt(-5,mtemp,ktemp,tmol,ac1ln,dadt,tx))
              if(tx.lt.tmol(1)) then
                if(i1.eq.0) ac2(it) = 0.
              else if(tx.gt.tmol(ktemp)) then
                if(iN.eq.0) ac2(it) = 0.
              end if

              if(l_per_stellar .eq. 1) then
                 opimol(it)=ac2(it)
              else
                if(kmol(imol) .eq. 0) then
                  partp=1
                  ac2(it)=ac1(it)
                else
                  partp = pressave(it,kmol(imol))/(pg(it)*emu(it))    ![mol_pp/mol*]/[g*/mol*]
                end if

                opimol(it)=partp*ac2(it)      			  ![mol_pp/g*][cm2/mol_pp]=[cm2/g*]
              end if

              if(opimol(it) .lt.0.d0) then
                write(6,*) ' ***hm:jv,k,wl_mu,imol,mol_op:',
     *            jv,it,1.e4/wnos(jv),imol,opimol(it)
              end if
! Add atomic/molecular opacity to continuum
              opacity(jv,it) = opacity(jv,it) + opimol(it)

283         continue


! Total opacity of molecule imol at depth it:
            dw = 0.
            if(jv.ge.2.and.jv.le.mopac-1) dw=(wnos(jv+1)-wnos(jv-1))/2.
            do it=1,ntau
              sumop(imol,it) = sumop(imol,it) + opimol(it)*dw 
              sumop2(it) = sumop2(it) + opimol(it)*dw
            end do
  
            nwimol_os = nwimol_os + 1
            avdw = avdw + abs(wn-wnos(jv))
            jv = jv +1
            go to 290
          end if
400     continue
499     continue
        close(unit=4)

        write(6,2831) (sumop(imol,it),it=1,ntau,20)
2831    format(5x,'opacity in cm/g* at 1,21,41,61,81:',1p8e9.2)
123   continue
199   continue

      do it=1,ntau
        write(*,'(i5,3e15.5)'), it, pg(it), sumop2(it), sumop3(it)
      end do
      print *

 
      write(6,*) ' *** out of OSopacity *** '
      
      return
      end
      
!-----------------------------------------------------------------------
! OSWAVENUMBERS: Computes or reads the wavenumbers for OS
!-----------------------------------------------------------------------
      subroutine oswavenumbers
      
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      common /csplim/vspmin,vspmax,ilist
      common /cwnos/abder(20,ndp),wnos(mopac),nwnossp,nwmin,nwmax,nosav
      common/cos/wlos(mopac),wlosstep(mopac),osresl,ktemp
      character inwnfil*60
      
      if(ilist .ne. 1) then
        wnos(1) = vspmin
        step = 1.d0 + 1.d0/osresl

        do k = 1,mopac-1
          wnkj = wnos(k)
          wnkj = wnkj * step
          nwnossp = k + 1
          wnos(nwnossp) = wnkj
          if(wnos(nwnossp) .gt. vspmax) exit
        end do
      
        if(wnos(nwnossp) .lt. vspmax) then
          wnkj = wnos(1)
          k = 1
          do
            wnkj = wnkj * step
            nwnosneed = k + 1
            if(wnkj .gt. vspmax) exit
            k = k + 1
          end do
          print *, 'Increase MOPAC to ', nwnosneed
          stop
        end if
 
        do i=1,nwnossp
          wlos(nwnossp-i+1)=1.d8/wnos(i)
        end do
      
        do i=2,nwnossp-1
          wlosstep(i) = (wlos(i+1)-wlos(i-1))/2.
        end do
        wlosstep(1) = wlos(2)-wlos(1)
        wlosstep(nwnossp) = wlos(nwnossp)-wlos(nwnossp-1)
      else
        open(unit=15,file='oswn.txt')
        do i=1,15
          read(15,*)
        end do
        
        i = 1
        do
          read(15,*,iostat=io) wl, f_wl
          if(io .ne. 0) exit
          if(f_wl .lt. 0) cycle
          wlos(i) = 1e4*wl
          i = i + 1
        end do
        nwnossp = i-1
        do i=1,nwnossp
          wnos(nwnossp-i+1)=1d8/wlos(i)
        end do
        vspmin = wnos(1)
        vspmax = wnos(nwnossp)
        osresl = 1
      end if

      write(6,'(a16,f6.1,1x,a1,f8.1,a6)') 
     *  'Frequency range:', wnos(1),'-',wnos(nwnossp),'cm^-1'
     

      return
      end

!-----------------------------------------------------------------------
! KAPPALINE: Calculates line absorption coefficient in each spectrum
! point. 
!-----------------------------------------------------------------------
       subroutine kappaline(jvm)

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      DIMENSION ABDERK(20),ABSKK(20),
     * TMOLKAP(MTEMP+1),FT(MTEMP+1),DT(MTEMP+1),Wt(3,mtemp+1)
      CHARACTER model_file*45,molid*4
      COMMON /COPACITY/contopac(mopac,ndp),opacity(mopac,ndp)
      COMMON /CMODEL/TEFF,G,ABUND(16),
     * TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP), 
     * PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP),
     * XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,model_file
      COMMON /CTRAN/TAU(NDP),X(NDP),BPLAN(NDP),HSURF,Y1(6),JTAU
      COMMON /CSPLIM/VSPMIN,VSPMAX
      COMMON /CMOLABS/ KMOL(0:mmol)
      COMMON /CNU/CAPNU(NDP),CAPCON(NDP),TAUNU(NDP)
      COMMON /CWNOS/ABDER(20,NDP),WNOS(MOPAC),NWNOSSP,NWMIN,NWMAX,nosav
      COMMON /CWNL/WNL(MOPAC),NWNL,NWLMIN,NWLMAX

      x(1:ntau)=0.
      capnu(1:ntau)=0.

      do k=1,ntau
        capnu(k)=capnu(k)+opacity(jvm,k)      			!line abs.coef. in cm2/g-stell-mat
        if(xkapr(k).ne.0.)x(k)=x(k)+opacity(jvm,k)/xkapr(k)     !ditto rel. ross mean abs.coef.
      end do
      
      return
      end    

!-----------------------------------------------------------------------
! TAUNUSCALE
!-----------------------------------------------------------------------     
      SUBROUTINE TAUNUSCALE(JVM)

      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
      DIMENSION AI(6),PPL(6),PP(NDP),Y3(3),X3(3),G3(3),GG(NDP)
     *  ,CAPLOG(NDP),DERCAP(NDP),PPLOG(NDP),PGLOG(NDP)
     *  ,W(3,NDP),TAUNUX(NDP),TAUNUY(NDP)
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUDUM(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,model_file
      COMMON /CNU/CAPNU(NDP),CAPCON(NDP),TAUNU(NDP)
      COMMON /CTRAN/TAU(NDP),X(NDP),BPLAN(NDP),HSURF,Y1(6),JTAU
      CHARACTER model_file*45
C
      DO 102 I=1,NTAU
      K=MAX(1,I-1)
      PP(I)=PG(I)+PRAD(I)+0.5*(PTURB(K)+PTURB(K+1))
      GG(I)=G
      CAPLOG(I)=LOG(CAPNU(I))
      PPLOG(I)=LOG(PP(I))
      pGlog(i)=LOG(PG(I))
102   CONTINUE
C       WRITE(6,10) ATMOS,TEFF,log10(G),10.**(ABUND(3)-ABUND(5))
10      FORMAT(' in taunuscale: ATMOSPHERE = ',A45,/
     &         ' (Teff,log(g),C/O) ='
     &         ,F6.0,F6.2,F8.2)
      TAUNU(1)=CAPNU(1)*PP(1)/GG(1)*(PP(2)-PP(1))
C     TAUNU(1)=CAPNU(1)*PP(1)/GG(1)
C..
C THIS ROUTINE CALCULATES THE DERIVATIVE DERCAP IN ALL NTAU NODES WHEN
C CAPNU AND PP ARE GIVEN 
          CALL TB04AT(NDP,NTAU,PPLOG,CAPLOG,DERCAP,W)
C
        DO 100 K=2,NTAU
        TAUNU(K)=TAUNU(K-1)
        DELPP=(PPLOG(K)-PPLOG(K-1))/10.
        DO 100 I=1,10
        PPI=PPLOG(K-1)+(I-0.5)*DELPP
        CAPINT=EXP(TG01BT(-1,NDP,NTAU,PPLOG,CAPLOG,DERCAP,PPI))
        TAUNU(K)=TAUNU(K)+CAPINT*DELPP/GG(K)
100     CONTINUE
C
C
        TAUNUY(1)=TAUNU(1)
        DO 130 K=2,NTAU
        TAUNUY(K)=TAUNUY(K-1)+(CAPNU(K)+CAPNU(K-1))/(GG(K-1)+GG(K))
     *          *(PP(K)-PP(K-1))
130     CONTINUE
C
C        if(jvm/10000*10000 .eq. jvm) then
C        write(6,*) ' jvm= ',jvm
C        write(6,*) 
C     &   'k,log10(taunu(k)),ditto{y},t(k),pg(k),pe(k),emu(k)'
C        do 613 k=1,ntau,20
C        write(6,612) 
C     &   k,log10(taunu(k)),log10(taunuy(k)),t(k),pg(k),pe(k),emu(k)
C612     format(i3,2f7.2,f7.0,1p3e11.3)
C613     continue
C        end if

C
C       DO 140 K=1,25
        DO 140 K=1,NTAU
        TAUNU(K)=TAUNUY(K)
140     CONTINUE
C
        isimple = 0
        if(isimple.eq.1) go to 997
C
C
        TAUNUX(1)=CAPNU(1)/XKAPR(1) * TAU(1)
        DO 120 K=2,NTAU
        TAUNUX(K)=(CAPNU(K)/XKAPR(K)+CAPNU(K-1)/XKAPR(K-1))/2.
     &     *(TAU(K)-TAU(K-1)) + TAUNUX(K-1)
120     CONTINUE
C
        DELTAU=TAUNU(25)-TAUNUX(25)
        DO 141 K=26,NTAU
        TAUNU(K)=TAUNUX(K)+DELTAU
141     CONTINUE
C

997     continue

      RETURN
      END
      
      
!-----------------------------------------------------------------------
C 'TRANSF' calculates the source function (Planck function) for each   
C model layer, and solves the transfer equation, for a given wavenumber
C The flux at the surface (tau=0) is stored in HSURF, and the          
C intensities in Y1.                                                   
C JVM is the number of steps from the spectrum beginning.              
C                                                                      
!-----------------------------------------------------------------------
      SUBROUTINE TRANSF(JVM)   
      
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
c     PARAMETER ( ndp=55, MOPAC=105000 )
      CHARACTER model_file*45
      COMMON /CANGLE/XMU(6),H(6),MMU
      COMMON /CTRAN/TAU(NDP),X(NDP),BPLAN(NDP),HSURF,Y1(6),JTAU
      COMMON P(NDP),SP1(NDP),SP2(NDP),SP3(NDP)
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,model_file
      COMMON /CNU/CAPNU(NDP),CAPCON(NDP),TAUNU(NDP)
      common /cwnos/abder(20,ndp),wnos(mopac),nwnossp,nwmin,nwmax,nosav
      COMMON /CWNL/WNL(MOPAC),NWNL,NWLMIN,NWLMAX
      DATA MMU/1/,XMU(1)/.57735023/,H(1)/1./
C
C PLANCK FUNCTION (ERG/CM**2/SEC/CM-1)
       BPL(TEMP,V)=1.191D-5*V*V*V/(EXP(1.439*V/TEMP)-1.)   
C
        WN = WNOS(JVM)
       DO 150 K=1,NTAU 
        BPLAN(K)=BPL(T(K),WN)
150    CONTINUE
C
C MU LOOP   
      JTAU1=JTAU-1  
      JTAU2=JTAU-2  
      HSURF=0.  
      DO 170 I=1,MMU
C   
C K=1
      DTAUB=.5*(X(1)+X(2))*(TAU(2)-TAU(1))/XMU(I)
C
      IF(DTAUB.EQ.0) WRITE(6,*) 'X(1),X(2),TAU(2,1),I,XMU(I) = '
      IF(DTAUB.EQ.0) WRITE(6,*) X(1),X(2),TAU(2),TAU(1),I,XMU(I)
      A=1./DTAUB
      B=A**2
      SP2(1)=1.+2.*A
      SP3(1)=-2.*B  
      C=2.*A
      TT=TAU(1)*(X(1))/XMU(I)
      EX=TT*(1.-.5*TT*(1.-.3333*TT))
      IF(TT.GT.0.1) EX=1.-EXP(-TT)
      P(1)=BPLAN(1)
      S0=P(1)
      P(1)=P(1)*(1.+C*EX)   
C
C K=2,JTAU-1
      DO 100 K=2,JTAU1  
      DTAUA=DTAUB   
      DTAUB=.5*(X(K)+X(K+1))*(TAU(K+1)-TAU(K))/XMU(I)
      DTAUC=.5*(DTAUA+DTAUB)
      AD=.166667*DTAUA/DTAUC
      CD=.166667*DTAUB/DTAUC
      SP1(K)=-1./(DTAUA*DTAUC)+AD
      SP2(K)=1. 
      SP3(K)=-1./(DTAUB*DTAUC)+CD   
100   P(K)=BPLAN(K)-AD*(BPLAN(K)-BPLAN(K-1))-CD*(BPLAN(K)-BPLAN(K+1))   
C   
C K=JTAU
      P(JTAU)=BPLAN(JTAU)   
      SP2(JTAU)=1.
C   
C ELIMINATE SUBDIAGONAL, SAVE FACTORS IN SP1
      DO 120 K=1,JTAU2  
      SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))  
      P(K+1)=P(K+1)+SP1(K)*P(K) 
      SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)   
120   SP2(K)=SP2(K)-SP3(K)  
121   SP2(JTAU-1)=SP2(JTAU-1)-SP3(JTAU-1)   
C   
C BACKSUBSTITUTE
      DO 160 K=1,JTAU1  
      P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K)*P(JTAU-K+1))/SP2(JTAU-K) 
160   CONTINUE  
C   
C END OF MU LOOP
      R1=P(1)-S0*EX 
      P0=P(1)*(1.-EX)+.5*S0*EX**2   
      HSURF=HSURF+H(I)*XMU(I)*P0
      Y1(I)=2.*P0   
C HSURF AND Y1(6) ARE THE FLUX AND INTENSITIES AT THE SURFACE   
170   CONTINUE  
C
C
      RETURN
      END
      
!-----------------------------------------------------------------------
! BPLAM: Calculate Planck function as function of lambda
! T*X must be greater than 1.7e6 to give BPL > 1.D-37
! BPL therefore limited to be > 1.D-30. UGJ 900510
!-----------------------------------------------------------------------
      function bplam(t,x)

      implicit real*8 (a-h,o-z)
      common /bplc/ex,x5
      data cp/1.191e27/,c2/1.438e8/

      x5=((x**2)**2)*(x/cp)
      ex=exp(-c2/(t*x))
      bplam=ex/((1.-ex)*x5)
      bplam = max(1.0d-30,bplam)
      return

      entry divbp(t,x)
      x6=x5*x
      tex=t*(1.-ex)
      bplam=c2*(ex/tex)/(tex*x6)
      bplam = max(1.0d-30,bplam)
      return
      
      end
      
!-----------------------------------------------------------------------
! TB04AT: Calculates the derivative of X(I) in the N knots. Succesfull
! if W = 0. 
!-----------------------------------------------------------------------
      subroutine tb04at(nd,n,x,f,d,w)

      implicit real*8 (a-h,o-z)
      dimension x(nd),f(nd),d(nd),w(3,nd)   
      
! First point
      cxb=1./(x(2)-x(1))
      cxc=1./(x(3)-x(2))
      dfb=f(2)-f(1) 
      dfc=f(3)-f(2) 
      w(1,1)=cxb*cxb
      w(3,1)=-cxc*cxc   
      w(2,1)=w(1,1)+w(3,1)  
      d(1)=2.*(dfb*cxb*cxb*cxb-dfc*cxc*cxc*cxc) 
 
! Interior points
      n1=n-1
      do k=2,n1 
        cxa=cxb   
        cxb=1./(x(k+1)-x(k))  
        dfa=dfb   
        dfb=f(k+1)-f(k)   
        w(1,k)=cxa
        w(3,k)=cxb
        w(2,k)=2.*(cxa+cxb)   
        d(k)=3.*(dfb*cxb*cxb+dfa*cxa*cxa) 
      end do

! Last point
      w(1,n)=cxa*cxa
      w(3,n)=-cxb*cxb   
      w(2,n)=w(1,n)+w(3,n)  
      d(n)=2.*(dfa*cxa*cxa*cxa-dfb*cxb*cxb*cxb) 

! Eliminate at first point  
      c=-w(3,1)/w(3,2)  
      w(1,1)=w(1,1)+c*w(1,2)
      w(2,1)=w(2,1)+c*w(2,2)
      d(1)=d(1)+c*d(2)  
      w(3,1)=w(2,1) 
      w(2,1)=w(1,1) 

! Eliminate at last point   
      c=-w(1,n)/w(1,n-1)
      w(2,n)=w(2,n)+c*w(2,n-1)  
      w(3,n)=w(3,n)+c*w(3,n-1)  
      d(n)=d(n)+c*d(n-1)
      w(1,n)=w(2,n) 
      w(2,n)=w(3,n) 

! Eliminate subdiagonal 
      do k=2,n  
        c=-w(1,k)/w(2,k-1)
        w(2,k)=w(2,k)+c*w(3,k-1)  
        d(k)=d(k)+c*d(k-1)
      end do

! Backsubstitute
      d(n)=d(n)/w(2,n)  
      do kk=2,n 
        k=(n+1)-kk
        d(k)=(d(k)-w(3,k)*d(k+1))/w(2,k)  
      end do

      return
      end

!-----------------------------------------------------------------------
! TG01BT: Interpolate F(X) at point XX
! Common values i1 and in controls what to do if xx is outside x interval.
!-----------------------------------------------------------------------
!      FUNCTION TG01BT(II,ND,N,X,F,D,XX) 
!      implicit real*8 (a-h,o-z)
!      DIMENSION X(ND),F(ND),D(ND)
!      COMMON /TG01BA/I1,IN,KK
!
!      if(ii.lt.0) kk=2  
!      
!      if(xx .lt. x(1)) then
!        tg01bt=f(1)
!        return
!      elseif(xx .gt. x(n)) then
!        tg01bt=f(n)
!        return
!      else
!        do k=kk,n
!          if(xx.lt.x(k)) exit
!        end do
!        dx=x(kk)-x(kk-1)  
!        df=f(kk)-f(kk-1)  
!        tg01bt = (df/dx)*xx + (f(kk) - (df/dx)*x(kk))
!        return
!      end if
!
!      end 
!      
     
!-----------------------------------------------------------------------
! TG01BT: Interpolate F(X) at point XX
! Common values i1 and in controls what to do if xx is outside x interval.
!----------------------------------------------------------------------- 
      function tg01bt(ii,nd,n,x,f,d,xx) 
      implicit real*8 (a-h,o-z)
      dimension x(nd),f(nd),d(nd)
      common /tg01ba/i1,in,kk

      if(ii.lt.0) kk=2  
      
      if(xx.lt.x(1)) goto 110   
      if(xx.gt.x(n)) goto 120   
      
      do 100 k=kk,n 
        if(xx.lt.x(k)) goto 101   
100   continue  
      kk=n  
      goto 102  
101   kk=k  

102   dx=x(kk)-x(kk-1)  
      df=f(kk)-f(kk-1)  
        tg01bt = (df/dx)*xx + (f(kk) - (df/dx)*x(kk))
      return

110   tg01bt=0.  
      tg01bt=f(1)
      if(i1.eq.1) return
      return

120   tg01bt=0.  
      tg01bt=f(n)
      if(in.eQ.1) RETURN
      END


!-----------------------------------------------------------------------
! FILTER: Computes the magnitude in prespecified filters
!-----------------------------------------------------------------------
      SUBROUTINE FILTER
C
C  This routine computes the magnitude in prespecified filters, by use
C  of the computed spectrum.
C
      implicit real*8 (a-h,o-z)
      include 'parameter.inc'
c     PARAMETER(MOPAC=105000,NDP=55)
      DIMENSION Z(25),FILTER_FILE(25),esofile(25),BPWMAG(25)
     &  ,scan_lam(29),asc(29,4),isc(29),AMAG_lam(25),BPWMAG_lam(25)
     &  ,filter_name(25),namewr(25),amagwr(25)
     &  ,amagnorm(25),amagnorm_lam(25)

      CHARACTER FILTER_FILE*20,model_file*45,filter_name*7
     &         ,namewr*7,eso_filters*7,esofile*20
      LOGICAL first_model
      COMMON/CFOLD/FLUX(0:MOPAC),FLUXINST(MOPAC),WLINST(0:MOPAC)
     &         ,fluxinstd(mopac)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      COMMON /CFILTER/FLUX_BPW(0:MOPAC),AMAG(25),NFILT,NESO,IFILT(25)
     &       ,iscan
      common /CWING8/ amagw8(10), bpwmagw8(10)
      COMMON /CWNOS/ABDER(20,NDP),WNOS(MOPAC),NWNOSSPNWMIN,NWMAX,nosav
      COMMON /CMODEL/TEFF,G,ABUND(16)  
     & ,TAUMOD(NDP),T(NDP),PE(NDP),PG(NDP),PRAD(NDP) 
     & ,PTURB(NDP),XKAPR(NDP),RO(NDP),EMU(NDP)
     & ,XL(20),ABSKA(20,NDP),SPRIDA(20,NDP),NTAU,NLP,model_file
      dimension eso_filters(6)
      dimension wl(1442),tr(1442),dtdw(1442),wtw(3,1442)
      data nfiltdim /1442/
      DATA Z/-20.121,-20.301,-19.598,-17.623,-18.558,-18.571,19*0./
C U,B,V,J,H,K,co,h2o,cont,w8,U,B,B,V,R,I,U,L
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      DATA FILTER_FILE 
     -  / 'trans_u_wing.dat','trans_b_wing.dat','trans_v_wing.dat',
     -  'trans_j_wing.dat', 'trans_h_wing.dat', 'trans_k_wing.dat',
     - 'trans_wing_co.dat', 'trans_wing_h2o.dat', 'trans_wing_cont.dat'
     -  ,'trans_wing8.dat',
     -  'trans_u_bessell.dat','trans_b_bessell.dat',
     -  'trans_bx_bessell.dat','trans_v_bessell.dat',
     -  'trans_r_bessell.dat','trans_i_bessell.dat',
     -  'trans_u_kurucz.dat','trans_l_josef.dat',7*' ' /
C   filter # correspond to:
       data filter_name / ' U_wing', ' B_wing', ' V_wing', 
     &          ' J_wing', ' H_wing', ' K_wing',
     &          ' W_CO  ', ' W_H2O ', ' W_cont', ' Wing8 ', 
     &          ' U_besl', ' B_besl', 'BX_besl', 
     &          ' V_besl', ' R_besl', ' I_besl', 
     &          ' U_kucz', ' L_hron', 7*'       ' /
      data first_model /.true./
      data scan_lam / 7024,  7144,  7564,  7812,  8116,  8140,
     &                8718,  8834,  8880,  9044,
     &  9084,  9164,  9190,  9234,  9284,  9316,  9576,  9704,
     &  9820,  9960, 10104, 10154, 10294, 10334, 10374, 10404,
     & 10564, 10834, 10976/
      DATA esofile
     -  / 'bessel-U632.dat','bessel-B450.dat','bessel-V451.dat',
     &    'bessel-R452.dat','gunn-i425.dat','gunn-zS759.dat',19*' '/
       data eso_filters / 'besU632', 'besB450', 'besV451', 
     &  'besR452', 'guni425', 'gunz759'/ 
 
      if (first_model) then
       OPEN(UNIT=53,FILE='colors.add',STATUS='unknown')
       nfwr = 0
       DO 261 K=1,NFILT
       if (ifilt(k).eq.0) go to 261
       nfwr = nfwr + 1
       namewr(nfwr) = filter_name(k)
       if(neso.eq.1) namewr(nfwr) = eso_filters(k)
261    continue
       write(6,156) nfilt,nfwr,(namewr(k),k=1,nfwr)
156    format( 'nfilt,nfwr,namewr()= '/,i3,i3,9(1x,A7)/,10(1x,A7))
       WRITE(53,262) (namewr(j),j=1,nfwr)
262    FORMAT (' Teff logg logO C/O  ', 8(1x,A7)/,10(1x,A7))
       first_model = .false.
      end if
 
      mwing = 0
 
 
       nfwr = 0
      DO 200 K=1,NFILT
      if (ifilt(k).eq.0) go to 200
       if(neso.eq.1) then
      OPEN(UNIT=52,FILE='../esofilters/'//esofile(k),
     *    STATUS='OLD',readonly)
          write(6,212) k,'../esofilters/'//esofile(k)
       else
      OPEN(UNIT=52,FILE='../filters/'//FILTER_FILE(K),
     *    STATUS='OLD',readonly)
          write(6,212) k,'../filters/'//FILTER_FILE(K)
       end if
212   format(' opened filter #: ',i3,'  = ',a27)
 
      if (filter_file(k).eq.'trans_wing8.dat     ') then
         do 105 i=1,5
105      read(52,*)
      end if
      kwing = 0
197   continue                  !go here if new wing filter in file
      if (filter_file(k).eq.'trans_wing8.dat     ')  then
          read(52,*) nwing,wcw,whdw
          kwing = kwing + 1
      end if
      ifirst = 1
      DO 104 I=1,nfiltdim
      WL(I) = 0.
104   TR(I) = 0.
      i101 = 0
      nos01 = 0
      nos02 = 0
      sumtr = 0.
 
      READ(52,*) nanometers
C            if nanometers = 1 wl(i) is assumed to be in nanometers
C            if nanometers = 0 wl(i) is assumed to be in aangstrom
      DO 102 I=1,10000
      READ(52,*,END=198) WL(I),TR(I)
      if(nanometers.eq.1) wl(i) = wl(i) * 10.
      sumtr = sumtr + tr(i)
      nwfilt = i
      if (tr(i).gt.9.9) go to 198  !end of a filter in the Wing_8color system
      if (i.gt.1442) stop ' increase dimension for filter'
102   continue
198   continue
      wf1 = wl(1)
      wf2 = wl(2)

      if(wl(1).gt.wl(2)) then
       DO 103 i=1,nwfilt
       if(wl(i).le.wl(i+1)) then
          go to 107
       end if
       wl1 = wl(i)
       tr1 = tr(i)
       wl(i) = wl(nwfilt-i+1)
       tr(i) = tr(nwfilt-i+1)
       wl(nwfilt-i+1) = wl1
       tr(nwfilt-i+1) = tr1
103    continue
107    continue
      end if

C      write(6,*)' we flipped the ',nwfilt,' transmission values'
C      write(6,*)
C     & 'first wl{1,2} were ',wf1,wf2,' now they are ',wl(1),wl(2)

       DO 108 i=1,nwfilt
       if(i.lt.nwfilt .and. wl(i).ge.wl(i+1)) then
          write(6,*) 'problems at i, wl(i)= ',i,wl(i)
          write(6,*) 'for filter file k = ',k
          stop ' control filter file'
       end if
      if (tr(i).gt.0. .and.tr(i).lt.9.9) then
         w02 = wl(i)
         if(i101.eq.0) then
         w01 = wl(i)
         i101 = 1
         end if
      end if
108    continue

      CALL TB04AT(nfiltdim,nwfilt,wl,tr,dtdw,wtw)

C nwmin, nwmax are the first and last wnos-index which gives
C nwos inside spectrum interval vspmin-vspmax
      II = -1
      AMAG(K) = 0.
      BPWMAG(K) = 0.
      AMAG_lam(K) = 0.
      BPWMAG_lam(K) = 0.
      ifirstos = 0
      filt_intgr = 0.
      filt_intgr_lam = 0.


      do 110 j=nwmax,nwmin,-1
      wlosf = 1.e8/wnos(j)
      if (wlosf .lt. wl(1)) go to 110
      if (ifirstos.eq.0) then
        jf1 = j
        ifirstos = 1
      end if
      if (wlosf .gt. wl(nwfilt)) go to 111
      jf2 = j

C For all OS-points overlapping with the filter, i.e. 
C from wnos(jf1) to wnos(jf2), compute the 
C contribution of this OS point to the magnitude, by first 
C calculating the transmissioncoefficient, tros, at this OS point

      do 210 iw1 = 1,nwfilt-1
      if (wlosf.ge.wl(iw1) .and. wlosf.lt.wl(iw1+1)) then
         iwf1 = iw1
         iwf2 = iw1+1
         if ( wl(iwf2).lt.w01 ) then
           nos01 = nos01 + 1
           go to 110
         end if
         if ( wl(iw1).gt.w02 ) then
           nos02 = nos02 + 1
           go to 110
         end if
         tros = tr(iw1) +
     &   (wlosf-wl(iw1) ) / (wl(iw1+1)-wl(iw1) ) * (tr(iw1+1)-tr(iw1))
         if(tros.le.1.d-20) write(6,*) 'small_tros at iw1,wl(iw1)=',
     &                  iw1,wl(iw1)
         go to 211
      end if
210   continue
211   continue

C Estimate distance between present OS-points, dwc=delta_cm-1 
C dwc_lam=and delta_AA

C     tros = TG01B(II,nfiltdim,nwfilt,wl,tr,dtdw,wlosf)
      if((j+1).le.nwmax .and. (j-1).ge.nwmin) then
         dwc = (WNOS(J+1)-WNOS(J-1))/2.
         dwc_lam = (1.e8/WNOS(J-1)-1.e8/WNOS(J+1))/2.
      else
         if((j+1).le.nwmax) then 
            dwc = wnos(j+1)-wnos(j)
            dwc_lam = 1.e8/WNOS(J)-1.e8/WNOS(J+1)
         else if((j-1).ge.nwmax) then
            dwc = wnos(j)-wnos(j-1)
            dwc_lam = 1.e8/WNOS(J-1)-1.e8/WNOS(J)
         else
            dwc = 1.
            dwc_lam = 1.
            write(6,*) 'problems with dw for J=',J
         end if
       end if
         if(dwc.le.1.d-20) write(6,*) 'small_dwc at j,wnos(j)=',
     &                  j,wnos(j)
         if(flux(j).le.1.d-20) write(6,*) 'small_flux at j,flux(j)=',
     &                  j,flux(j)

C now sum up the flux*filtertransmission*delta_cm-1 at all OS-points
C compute also the integral of the transmission as
C filtertransmission*delta_cm-1 summed over all OS-points

      AMAG(K) = AMAG(K) + FLUX(J)*TROS*DWC
      AMAG_lam(K) = 
     &  AMAG_lam(K) + FLUX(J)*wnos(j)**2*1.e-8*TROS*DWC_lam
      BPWMAG(K) = BPWMAG(K) + FLUX_BPW(J)*TROS*DWC
      BPWMAG_lam(K) = 
     &  BPWMAG_lam(K) + FLUX_BPW(J)*wnos(j)**2*1.e-8*TROS*DWC_lam
      filt_intgr = filt_intgr + tros*dwc
      filt_intgr_lam = filt_intgr_lam + tros*dwc_lam
      if (flux(j).le.0. .or. flux_bpw(j).le.0. .or.dwc.le.0.
     &    .or.(tros.le.0. .and.tr(iwf1)+tr(iwf2).ne.0.)) write(6,316)
     &    j,1.e8/wnos(j),flux(j),flux_bpw(j),tros,dwc
316   format(' filter error:',i5,f10.2,1p4e12.3)
110   continue
111   continue

      write(6,204)
     &  jf1,jf2,1.e8/wnos(jf1),1.e8/wnos(jf2),wl(1),wl(nwfilt)
204   format(' os wavelength interval: ',2i7,2f8.1/,
     &       ' filter wavelength interval: ',2f8.1)
      write(6,*) ' raw filter fluxes:',amag(k),bpwmag(k),filt_intgr
      write(6,*) ' nwmax,nwmin,jf1,jf2= ', nwmax,nwmin,jf1,jf2

C finally compute the magnitudes as -2.5log_10() and the corresponding
C normalisations (integrals of the filtertransmissions).

      if (amag(k).gt.1.d-20 .and. bpwmag(k).gt.1.d-20
     &   .and. filt_intgr .gt.1.d-20 ) then
C     AMAG(K) = -2.5*LOG10(AMAG(K)/filt_intgr)
C     BPWMAG(K) = -2.5*LOG10(BPWMAG(K)/filt_intgr)
C     AMAG_lam(K) = -2.5*LOG10(AMAG_lam(K)/filt_intgr_lam)
C     BPWMAG_lam(K) = -2.5*LOG10(BPWMAG_lam(K)/filt_intgr_lam)
      AMAG(K) = -2.5*LOG10(AMAG(K))
      BPWMAG(K) = -2.5*LOG10(BPWMAG(K))
      amagnorm(k) = -2.5*log10(filt_intgr)
      write(6,*) 'k,amag(k)= ',k,amag(k) 
      else
      write(6,*) 'k,amag(k)= ',k,amag(k) 
      AMAG(K) = 999.
      BPWMAG(K) = 999.
      amagnorm(k) = 999.
      end if
      if (amag_lam(k).gt.1.d-20 .and.bpwmag_lam(k).gt.1.d-20
     &   .and. filt_intgr_lam .gt.1.d-20 ) then
      AMAG_lam(K) = -2.5*LOG10(AMAG_lam(K))
      BPWMAG_lam(K) = -2.5*LOG10(BPWMAG_lam(K))
      amagnorm_lam(k) = -2.5*log10(filt_intgr_lam)
      else
      AMAG_lam(K) = 999.
      BPWMAG_lam(K) = 999.
      amagnorm_lam(k) = 999.
      end if
C
C
      write(6,313)k,wl(1),wl(nwfilt),nwfilt
     &      ,jf1-jf2+1,nos02-nos01+1
     &      ,amag(k),bpwmag(k), amagnorm(k)
     &      ,amag_lam(k),bpwmag_lam(k), amagnorm_lam(k)
313   format(i3,2f8.1,i4,i5,i3,3f8.3,/47x,3f8.3)
312   format(' Filter ',i3,'; first/last filter wl= ',2f8.1
     &       ,' AA '/
     &       ,' There were ',i4,' filter frequencies in trans-file'/
     &       ,' first/last OS wl= ',2f8.1,'AA'/
     &       ,' There were ',i5,' OS frequencies in filter region'/
     &  ,' Raw *-mag,Planck_mag, filter.norm. (cm-1 and ditto AA):'
     & ,2f8.3)
      if (filter_file(k).eq.'trans_wing8.dat     ') then
          amagw8(kwing) = amag(k)
          bpwmagw8(kwing) = bpwmag(k)
          if (tr(nwfilt+1).gt.9.9) go to 197   !new filter in the Wing_8color system
          mwing = kwing
      end if
      CLOSE(52)
200   CONTINUE
C
      write(6,*) ' We finished filter computations and write first:'
      write(6,*) ' k,filter(k),mag,Planck_mag,mag_norm {cm^-1,AA}:'
       nfwr = 0
      DO 280 K=1,NFILT
      if (ifilt(k).eq.0) go to 280
       nfwr = nfwr + 1
C      write(6,281) k,filter_name(k)
      write(6,281) k,namewr(nfwr)
     &      ,amag(k),bpwmag(k), amagnorm(k)
     &      ,amag_lam(k),bpwmag_lam(k), amagnorm_lam(k)
281   format(i2,1x,a7,3f7.2,1x,3f7.2)
280   continue
      write(6,*) ' Then we second list:'
      write(6,*) ' mag_1-mag_2, magPl_1-magPl_2, magl_1-magl_2'
C234567890 234567890 234567890 234567890 234567890 234567890 234567890 2
      write(6,*) 
     &        ' where mag_x is -2.5(log intregral[flux_x*trans_x*dnu] )'
      write(6,*) 
     &' in cm^-1 units for spectrum & Planck fct, and then in AA-units'
      write(6,*) 
     &         ' Remark: the listed magnitudes are not normalized; i.e.'
      write(6,*) ' no zeropoints and no fluxnormalization (=int(tr*dw))'

      if (nfilt.ge.2 .and. ifilt(1)+ifilt(2).eq.2)
     &WRITE(6,20) 
     & AMAG(1)-AMAG(2),BPWMAG(1)-BPWMAG(2),AMAG_lam(1)-AMAG_lam(2)      !U-Bw
      if (nfilt.ge.3 .and. ifilt(3)+ifilt(2).eq.2)
     &WRITE(6,21)
     &  AMAG(2)-AMAG(3),BPWMAG(2)-BPWMAG(3),AMAG_lam(2)-AMAG_lam(3)      !B-Vw
      if (nfilt.ge.5 .and. ifilt(4)+ifilt(5).eq.2)
     &WRITE(6,23)
     &  AMAG(4)-AMAG(5),BPWMAG(4)-BPWMAG(5),AMAG_lam(4)-AMAG_lam(5)      !J-Hw
      if (nfilt.ge.6 .and. ifilt(5)+ifilt(6).eq.2)
     &WRITE(6,24) 
     &  AMAG(5)-AMAG(6),BPWMAG(5)-BPWMAG(6),AMAG_lam(5)-AMAG_lam(6)      !H-Kw

      if (nfilt.ge.8 .and. ifilt(7)+ifilt(8).eq.2)
     &WRITE(6,20) 
     &  AMAG(7)-AMAG(8),BPWMAG(7)-BPWMAG(8),AMAG_lam(7)-AMAG_lam(8)      !U-Bb
      if (nfilt.ge.9 .and. ifilt(7)+ifilt(9).eq.2)
     &WRITE(6,20) 
     &  AMAG(7)-AMAG(9),BPWMAG(7)-BPWMAG(9),AMAG_lam(7)-AMAG_lam(9)      !U-BXb
      if (nfilt.ge.10 .and. ifilt(8)+ifilt(10).eq.2)
     &WRITE(6,21) 
     & AMAG(8)-AMAG(10),BPWMAG(8)-BPWMAG(10),AMAG_lam(8)-AMAG_lam(10)      !B-Vb
      if (nfilt.ge.10 .and. ifilt(9)+ifilt(10).eq.2)
     &WRITE(6,21) 
     & AMAG(9)-AMAG(10),BPWMAG(9)-BPWMAG(10),AMAG_lam(9)-AMAG_lam(10)      !BX-Vb
      if (nfilt.ge.11 .and. ifilt(10)+ifilt(11).eq.2)
     &WRITE(6,25) AMAG(10)-AMAG(11)
     &             ,BPWMAG(10)-BPWMAG(11),AMAG_lam(10)-AMAG_lam(11)      !V-Rb
      if (nfilt.ge.12 .and. ifilt(11)+ifilt(12).eq.2)
     &WRITE(6,26) AMAG(11)-AMAG(12)
     &             ,BPWMAG(11)-BPWMAG(12),AMAG_lam(11)-AMAG_lam(12)      !R-Ib
      if (nfilt.ge.13 .and. ifilt(9)+ifilt(13).eq.2)
     &WRITE(6,20) 
     & AMAG(9)-AMAG(13),BPWMAG(9)-BPWMAG(13),AMAG_lam(9)-AMAG_lam(13)      !Uk-BXb


      nfwr = 0
      DO 260 K=1,NFILT
      if (ifilt(k).eq.0) go to 260
      nfwr = nfwr + 1
      amagwr(nfwr) = amag(k)
260   continue

      abuo = 12.+log10(abund(5))
      WRITE(53,22) teff, log10(g), abuo, abund(3)/abund(5),
     - (amagwr(j),j=1,nfwr)
22    FORMAT (F6.0,F4.1,2f5.2,1X,8F7.2/,11F7.2)
      if(mwing.gt.0) WRITE(53,28) (amagw8(j),j=1,mwing)
28    FORMAT (' Wing-8 mag: ',10F7.2)

10    FORMAT (' J_MAG    = ',F7.3)
11    FORMAT (' H_MAG    = ',F7.3)
12    FORMAT (' K_MAG    = ',F7.3)
13    FORMAT (' CO_MAG   = ',F7.3)
14    FORMAT (' H2O_MAG  = ',F7.3)
15    FORMAT (' CONT_MAG = ',F7.3)
16    FORMAT (' U_MAG = ',F7.3)
17    FORMAT (' B_MAG = ',F7.3)
18    FORMAT (' V_MAG = ',F7.3)
23    FORMAT (' J - H (direct; no zero correction) = ',3F7.3)
24    FORMAT (' H - K (direct; no zero correction) = ',3F7.3)
20    FORMAT (' U - B (direct; no zero correction) = ',3F7.3)
21    FORMAT (' B - V (direct; no zero correction) = ',3F7.3)
C23    FORMAT (' J - H (+0.071) = ',F7.3)
C24    FORMAT (' H - K (-0.008) = ',F7.3)
C20    FORMAT (' U - B (-0.24) = ',F7.3)
C21    FORMAT (' B - V (+0.596) = ',F7.3)
25    FORMAT (' V - R (direct; no zero correction) = ',3F7.3)
26    FORMAT (' R - I (direct; no zero correction) = ',3F7.3)
C

      if (iscan.eq.0) go to 999
C computation of 29 scanner filter magnitudes normalized to the 
C filter at 1.04mu. Each filter is centered at scan_lam(i), and
C has a squared filterfunction of full width of 30 AA.

      DO 400 K=29,1,-1

      AMAGK = 0.
      BPWMAGK = 0.
      ifirstos = 0
      do 410 j=nwmin,nwmax
      wlosf = 1.e8/wnos(j)
      if (wlosf .gt. scan_lam(k)+15.) go to 410
      if (ifirstos.eq.0) then
        jf1 = j
        ifirstos = 1
      end if
      if (wlosf .lt. scan_lam(k)-15.) go to 411
      jf2 = j
      dwc = (WNOS(J+1)-WNOS(J-1))/2.
      AMAGK = AMAGK + FLUX(J)*DWC
      BPWMAGK = BPWMAGK + FLUX_BPW(J)*DWC
      if (flux(j).le.0. .or. flux_bpw(j).le.0. .or.dwc.le.0.)
     &  write(6,316)
     &    j,1.e8/wnos(j),flux(j),flux_bpw(j),dwc
410   continue
411   continue
C     write(6,*) ' raw filter fluxes:',amag(k),bpwmag(k)

      AMAGK = -2.5*LOG10(AMAGK)
      BPWMAGK = -2.5*LOG10(BPWMAGK)

      asc(k,1) = 1.e8/wnos(jf2)
      asc(k,2) = 1.e8/wnos(jf1)
      isc(k) = jf2-jf1+1
      asc(k,3) = amagk
      asc(k,4) = bpwmagk
400   CONTINUE
C

      write(6,*) ' filter lam_min lam_max wnos_1 wnos_2 nwnos'
     &,' mag mag_Pl mag-mag_1.0404'
      do 470 k=1,29
      write(6,413) k, scan_lam(k)-15.,scan_lam(k)+15.
     &       ,asc(k,1),asc(k,2),isc(k)
     &       ,asc(k,3),asc(k,4),asc(k,3)-asc(26,3)
470   continue
413   format(i3,2f8.1,2f8.1,i5,1x,2f8.3,f10.3)


999   continue

      RETURN
      END





