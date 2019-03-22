************************************************************************
      subroutine smchem ( anHges, tg, eps, anmono, anmol, pel )
************************************************************************
*                                                                      *
*     small chemistry                                                  *
*     ---------------                                                  *
*     Diese Routine berechnet eine GG-Chemie                           *
*                                                                      *
*   - Wird der Parameter "nachit"=FALSE gesetzt, so berechnet          *
*     sie kein exaktes Dissoziationsgleichgewicht, sondern nutzt       *
*     das Wissen aus dem Irsee-paper. Nur die wichtigen Dissozia-      *
*     tionsgleichgewichte werden dann geloest und so die Dichten       *
*     der haeufigsten Molekuele bestimmt. Dabei loest die Routine      *
*     zunaechst die Gleichgewichte der haeufigsten Elemente und        *
*     geht dann weiter zu weniger haeufigen. Die bereits geloesten     *
*     GG werden als von selteneren Elementen unbeeinflusst angenommen. *
*     Hierin liegt natuerlich auch die Schwaeche der Routine: Sie geht *
*     von bestimmten chemischen Haeufigkeiten aus. Weichen diese       *
*     stark von den Vorstellungen der Routine ab, so werden die        *
*     Ergebnisse falsch (z.B. wenn  [C]/[O] ~ 1 ist!!!)                *
*                                                                      *
*     Welche Molekuele in diesem Fall beruecksichtigt werden, steht    *
*     unten im Text.                                                   *
*                                                                      *
*     Die Abweichungen der wichtigen Molekuele liegen in der Regel     *
*     unter einem Faktor 2, oft unter 10%.                             *
*     Die Namen der beteiligten Atome und Molekuele werden ueber den   *
*     commonblock chemnam an die Aussenwelt uebergeben. Mit Hilfe      *
*     der Funktion stindex (siehe unten) kann man sich leicht nach     *
*     dem ersten Aufruf der Routine smchem die Indices der einzelnen   *
*     Molekuele im array anmol verschaffen; z.B.                       *
*     nnco = stindex(cmol,dim,'CO')                                    *
*                                                                      *
*     Die Routine ist extrem schnell und beansprucht wenig Platz.      *
*                                                                      *
*                                                                      *
*   - Wird der Parameter "nachit"=TRUE gesetzt, so berechnet die       *
*     Routine die GG-Chemie korrekt nach Newton-Raphson, wobei         *
*     obige Konzentrationen als Startwerte verwendet werden.           * 
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     e i n g a b e :                                                  *
*     anHges : totale Dichte der Wasserstoffatome (~ nh+2*nh2)         *
*     tg     : Gastemperatur                                           *
*     eps    : Vektor mit linearen chemischen Haeufigkeiten ([H]=1)    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     a u s g a b e :                                                  *
*     anmono : vektor mit den dichten der monomere  falls relevant     *
*     anmol  : vektor mit den dichten der molekuele falls relevant     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     b e n o e t i g t e  u n t e r p r o g r a m m e :               *
*     stindex                                                          *
*                                                                      *
*     b e n o e t i g t e  d a t e i e n :                             *
*     dispol.dat                                                       *
*                                                                      *
************************************************************************
*                                                                      *
*     Die Nummerierung der Elemente ist wie folgt:                     *
*     ____________________________________________________________________ 
*     Element: He  e-  H   C   N   O   Si  Mg  Al  Fe  S   Na  K   Ti  Ca Li Cl Fl
*     Nummer:  1   2   3   4   5   6   7   8   9   10  11  12  13  14  15 16 17 18
*     ____________________________________________________________________
*                                                                      *
************************************************************************
*     (c)       Carsten Dominik                     Do 11. Maerz 1993  *
*     Revision und nachit-Erweiterung Peter Woitke  Do  7. Maerz 1996  *
*     GG-Konstanten fuer TiC von Andreas Gauger gewonnen               *
*     auf der Grundlage spekrtoskopischer Messungen          11.09.97  *
*     Vervollstaendigung der Ionen S+, Si+, Fe+, Ti+         20.10.97  * 
*     Calzium-Chemie Ca+, CaOH, Ca(OH)2, CaO, CaS, Ca2       28.08.00  * 
*     FeH,TiH from Burrows according to Sharp & Huebner ChH  21.01.08  *
*     CaH from Tsuji 1975                               ChH  21.01.08  *
*     LiCl,LiH,LiO, LiOH from Tsuji 1975                ChH  08.06.09  *
*     NaCl, CaCl, CaCl2,KCl, HCl, NaH from Sharp & Huebner 1990  ''    *
************************************************************************
      use drift_data,ONLY: cmol
      implicit none
*-----------------------------------------------------------------------
*  Dimensionierung fuer die Molekuel- und Atom Felder. Hier ist auf
*  Konsistenz mit dem aufrufenden Programm zu achten.
      integer nml,nel
      parameter (nml=167)
      parameter (nel=17)
*-----------------------------------------------------------------------
      real*8 eps(nel),anmono(nel),anmol(nml)
      real*8 anHges,tg
*-----------------------------------------------------------------------
*  Die Variable "alle" entscheidet, ob die nicht unmittelbar 
*  beruecksichtigten Molekuele dennoch inkonsistent mitgerechnet werden. 
      logical alle
      data alle/.true./
*-----------------------------------------------------------------------
*  Die Variable "nachit" entscheidet, ob die GG-Chemie zur Ermittlung
*  der korrekten Loesung nach Newton-Raphson nachiteriert wird.
*  alle muss dafuer alle=TRUE sein.
*  (braucht laenger, manchmal Konvergenzprobleme)
      logical nachit
      data nachit/.true./      
*-----------------------------------------------------------------------
*  Bei merk=.true. merkt sich die Routine die letzte konvergiert Loesung
*  und geht beim naechsten Mal von diesen Startwerten aus.
      logical merk
      data merk/.false./
*-----------------------------------------------------------------------
*  Die Variable "ngestst" entscheidet, ob die Elementerhaltung ueber-
*  prueft werden soll.
      logical ngestst
      data ngestst/.false./
*-----------------------------------------------------------------------
*  Hier kommen die Common-Bloecke fuer die Kommunikation nach aussen.      
      integer   HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,SiII,FeII
      integer   CaII
      common/ionnumm/HII,CII,NII,OII,NaII,MgII,AlII,KII,TiII,SII,SiII,
     &               FeII,CaII
*-----------------------------------------------------------------------
*  Die Variable tdispol bestimmt, welches die niedrigste Temperatur 
*  ist, die in die Dissoziationspolynome eingesetzt werden darf.
      real*8 tdispol
      parameter (tdispol=300.d0) 
*-----------------------------------------------------------------------
*  Um das Program lesbarer zu machen werden fuer jedes Molekuel und 
*  jedes Atom Integer-Variablen definiert, die den Index der Groesse 
*  im entsprechenden Feld beinhalten.
      integer He,el,H,C,N,O,Si,Mg,Al,Fe,S,Na,K,Ti,Ca,Li,Fl,Cl
      parameter(He=1,el=2,H=3,C=4,N=5,O=6,Si=7,Mg=8,Al=9,Fe=10,S=11,
     &Na=12,K=13,Ti=14,Ca=15,Li=16,Cl=17,Fl=18)
*-----------------------------------------------------------------------
      integer m_kind(0:4,nml),m_anz(4,nml)
      integer Al2O,AlH,AlO2H,AlOH,C2,C3,C2H,C3H,C2H2,CH4,CN,CO,CO2,CS
      integer FeS,H2,H2O,H2S,HCN,HS,MgH,MgO,MgOH,MgS,N2,NH3,O2,OH,SO,SO2
      integer Si2C,SiC,SiC2,SiH,SiH4,SiN,SiO,SiO2,SiS,FeO,FeO2H2
      integer TiO,TiO2,TiS,TiC,TiC2,MGO2H2,NAOH,NA2O2H2,CaOH,CaO2H2
      integer KOH,K2O2H2,FeH,CaH,TiH,LiH,LiO,LiOH,LiCl
      integer KCl,CaCl2,CaCl,NaCl,HCl,NaH
      integer stindex
      integer i,j,j1,nmol,jj,kk,ilauf,l,it,m1,m2,pCpOit
      real*8  a(nml,0:4),th1,th2,th3,th4,ppp,qqq
      real*8  pH,pHe,pC,pN,pO,pSi,pMg,pAl,pFe,pS,pNa,pK,pTi,pCa,pel
      real*8  pLi,pCl
      real*8  pHges,pHeges,pCges,pNges,pOges,pSiges,pMgges,pAlges,
     &        pFeges,pSges,pNages,pKges,pTiges,pCages,pCalt,pOalt,
     &        pLiges,pClges
      real*8  g(nml),amerk(nel)
      real*8  bk,kT,kT1,nelek,ng,Sa,Nenner,fak,lth,arg,term
      real*8  aa,bb,cc,dd,ee,dpp(2),func(2),dfunc(2,2),a3
      real*8  f0,f1,f2,f3,ff1,ff2,ff3,f,fs,delta,soll,haben,abw
      real*8  DF(nel,nel),dp(nel),FF(nel),pmol,pat,nges(nel),pmono1(nel)
      real*8  const,eV,me,pi,hplanck,atm,Rcal,RcalT,TT1,TT2,TT3
      real*8  GK,VIETA
      logical abbruch
      character*10 catm(nel)
      character*1  char
      data catm/'He','e-','H','C','N','O','Si','Mg','Al','Fe','S','Na',
     &          'K','Ti','Ca','Li','Cl'/
      data ilauf/0/
      data bk/1.380662D-16/, hplanck/6.626196D-27/ 
      data eV/1.602D-12/, me/9.109558D-28/, pi/3.141592653589793D0/
      data atm/1.013D+6/, Rcal/1.987D+0/
      save
*-----------------------------------------------------------------------
*  Die Formelfunktion fuer die Berechnung der Dissoziationskonstanten
c     gk(i) = DMIN1(1.d+300, DEXP( a(i,0) + a(i,1)*th1 + a(i,2)*th2
c    &                                    + a(i,3)*th3 + a(i,4)*th4 ) )
      GK(i) = DEXP(DMIN1(600.d0, a(i,0) + a(i,1)*th1 + a(i,2)*th2
     &                                  + a(i,3)*th3 + a(i,4)*th4 ) )
*-----------------------------------------------------------------------
*  Die Formelfunktion zur Loesung quadratische Gleichungen mit Vieta
      VIETA(ppp,qqq) = qqq/(-ppp/2.d0-DSQRT(ppp**2/4.d0-qqq))
*-----------------------------------------------------------------------      
*
      ilauf = ilauf + 1
      if ( ilauf .eq. 1 ) then
        
        ! Einlesen der Koeffizienten der Dissoziationspolynome
*       ======================================================
        
        open (unit=12, file='data/DRIFT/dispol_large.dat', status='old')
        write(*,*)
	write(*,*) 'Lese Dissoziationspolynome von dispol_large.dat'
        rewind(12)
        do i=1,nml
 2001     format(a10,10i3)
 2002     format(5e13.5)
          read(12,2001) cmol(i),m_kind(0,i),
     &                  (m_kind(j,i),j=1,m_kind(0,i)),
     &                  (m_anz(j,i),j=1,m_kind(0,i))
          read(12,2002) (a(i,j),j=0,4)
        enddo  
        nmol=nml
        write(*,*) nmol,' Species'
        close(12)
*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*     Heraussuchen der Namen aus cmol mit Hilfe von stindex
*     
        Al2O   = stindex(cmol,nml,'AL2O   ')
        AlH    = stindex(cmol,nml,'ALH    ')
        AlO2H  = stindex(cmol,nml,'ALO2H  ')
        AlOH   = stindex(cmol,nml,'ALOH   ')
        C2     = stindex(cmol,nml,'C2     ')
        C3     = stindex(cmol,nml,'C3     ')
        C2H    = stindex(cmol,nml,'C2H    ')
        C3H    = stindex(cmol,nml,'C3H    ')
        C2H2   = stindex(cmol,nml,'C2H2   ')
        CH4    = stindex(cmol,nml,'CH4    ')
        CN     = stindex(cmol,nml,'CN     ')
        CO     = stindex(cmol,nml,'CO     ')
        CO2    = stindex(cmol,nml,'CO2    ')
        CS     = stindex(cmol,nml,'CS     ')
        FeO    = stindex(cmol,nml,'FEO    ')
        FeO2H2 = stindex(cmol,nml,'FE(OH)2')
        FeS    = stindex(cmol,nml,'FES    ')
        H2     = stindex(cmol,nml,'H2     ')
        H2O    = stindex(cmol,nml,'H2O    ')
        H2S    = stindex(cmol,nml,'H2S    ')
        HCN    = stindex(cmol,nml,'HCN    ')
        HS     = stindex(cmol,nml,'HS     ')
        MgH    = stindex(cmol,nml,'MGH    ')
        MgO    = stindex(cmol,nml,'MGO    ')
        MgOH   = stindex(cmol,nml,'MGOH   ')
        MgO2H2 = stindex(cmol,nml,'MG(OH)2')
        MgS    = stindex(cmol,nml,'MGS    ')
        N2     = stindex(cmol,nml,'N2     ')
        NH3    = stindex(cmol,nml,'NH3    ')
        O2     = stindex(cmol,nml,'O2     ')
        OH     = stindex(cmol,nml,'OH     ')
        SO     = stindex(cmol,nml,'SO     ')
        SO2    = stindex(cmol,nml,'SO2    ')
        Si2C   = stindex(cmol,nml,'SI2C   ')
        SiC    = stindex(cmol,nml,'SIC    ')
        SiC2   = stindex(cmol,nml,'SIC2   ')
        SiH    = stindex(cmol,nml,'SIH    ')
        SiH4   = stindex(cmol,nml,'SIH4   ')
        SiN    = stindex(cmol,nml,'SIN    ')
        SiO    = stindex(cmol,nml,'SIO    ')
        SiO2   = stindex(cmol,nml,'SIO2   ')
        SiS    = stindex(cmol,nml,'SIS    ')
        TiO    = stindex(cmol,nml,'TIO    ')
        TiO2   = stindex(cmol,nml,'TIO2   ')
        TiS    = stindex(cmol,nml,'TIS    ')
        TiC    = stindex(cmol,nml,'TIC    ')
        TiC2   = stindex(cmol,nml,'TIC2   ')
        NaOH   = stindex(cmol,nml,'NAOH   ')
        Na2O2H2= stindex(cmol,nml,'NA2O2H2')
        CaOH   = stindex(cmol,nml,'CAOH   ')
        CaO2H2 = stindex(cmol,nml,'CA(OH)2')
        KOH    = stindex(cmol,nml,'KOH    ')
        K2O2H2 = stindex(cmol,nml,'K2O2H2 ')
        HII    = stindex(cmol,nml,'H+     ')
        CII    = stindex(cmol,nml,'C+     ')
        NII    = stindex(cmol,nml,'N+     ')
        OII    = stindex(cmol,nml,'O+     ')
        NaII   = stindex(cmol,nml,'NA+    ')
        MgII   = stindex(cmol,nml,'MG+    ')
        AlII   = stindex(cmol,nml,'AL+    ')
        KII    = stindex(cmol,nml,'K+     ')
        TiII   = stindex(cmol,nml,'TI+    ')
        SII    = stindex(cmol,nml,'S+     ')
        SiII   = stindex(cmol,nml,'SI+    ')
        FeII   = stindex(cmol,nml,'FE+    ')
        CaII   = stindex(cmol,nml,'CA+    ')
        FeH    = stindex(cmol,nml,'FEH    ')
        CaH    = stindex(cmol,nml,'CAH    ')
        TiH    = stindex(cmol,nml,'TIH    ')
        LiH    = stindex(cmol,nml,'LIH    ')
        LiO    = stindex(cmol,nml,'LIO    ')
        LiOH   = stindex(cmol,nml,'LIOH   ')
        LiCl   = stindex(cmol,nml,'LICL   ')
        KCl    = stindex(cmol,nml,'KCL    ')
        CaCl2  = stindex(cmol,nml,'CACL2  ')
        CaCl   = stindex(cmol,nml,'CACL   ')
        NaCl   = stindex(cmol,nml,'NACL   ')
        HCl    = stindex(cmol,nml,'HCL    ')
      endif
*-----------------------------------------------------------------------
*     ! zu niedrige Temperaturen abfangen und
*     ! Variable fuer die Dissoziationskonstanten berechnen
*     =====================================================
      th1   = 5040.d0 / DMAX1(tdispol,tg)
      th2   = th1*th1
      th3   = th2*th1
      th4   = th3*th1
      kT    = bk * DMAX1(tdispol,tg)
      kT1   = 1.d0/kT
*      
*-----------------------------------------------------------------------
*     ! Vektoren initialisieren
*     =========================
      do j = 1 , nel 
        anmono(j) = 0.d0
      enddo
      do j = 1 , nml 
        anmol(j)  = 0.d0
        g(j)      = 0.d0
      enddo
*
* --------------------------------------------------------------------------
*    TiC Gleichgewichtskonstante von Andreas Gauger ist anders
*        definiert als die Gleichgewichtskonstanten von Gail
*  Gauger: 
*  log Kp = 12.75293-5.4485*th1-1.56672*log(th1)+1.56041*(log(th1))**2
*           - 0.93275(log(th1))**3
*         = log ( p(A)p(B)/p(AB) )
*  Gail:
*   ln Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         =  ln ( p(AB)/p(A)p(B) )
*  Tsuji:
*  log Kp = a0 + a1*th1 + a2*th1**2 + a3*th1**3 + a4*th1**4
*         = log ( p(A)p(B)/p(AB) )
*
*  Umrechnung der Gauger-TiC-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Gauger) = -2.30256*log Kp(Gauger) = ln Kp(Gail)

        lth = DLOG10(th1)
        arg = 12.75293 - 5.44850*th1    - 1.56672*lth
     &                 + 1.56041*lth**2 - 0.93275*lth**3
        g(TiC) = DMIN1(1.D+300, DEXP(-2.30256*arg))

*  Umrechnen der Tsuji-CaH-GG-Konstante in das Gail'sche System
*  -log(10)*log Kp(Tsuji) = -2.30256*log Kp(Tsuji) = ln Kp(Gail)

        arg = 1.13401E+01 
     &       -3.01442E+00*th1 
     &       +4.23487E-01*th2
     &       -6.14674E-02*th3
     &       +3.16392E-03*th4
        g(CaH) = DMIN1(1.D+300, DEXP(-2.30256*arg))
c
        arg = 0.0
        arg =  1.17824E+01
     &        -3.44431E+00*th1  ! this term differs between Tsuji 1073 paper and molecBP data file!  
     &        +3.27412E-01*th2 
     &        -4.84559E-02*th3  
     &        +2.53356E-03*th4
        g(LiH) =  DMIN1(1.D+300, DEXP(-2.30256*arg))
c
        arg = 0.0
        arg =  1.23622E+01
     &        -4.54966E+00*th1
     &        +3.40687E-01*th2
     &        -5.00589E-02*th3
     &        +2.60132E-03*th4
        g(LiO) =  DMIN1(1.D+300, DEXP(-2.30256*arg))
c
c        arg = 0.0
c        arg =  1.24491E+01
c     &        -6.82101E+00*th1
c     &        +2.85775E-01*th2
c     &        -4.16715E-02*th3
c     &        +2.15772E-03*th4
c        g(LiF) =  DMIN1(1.D+300, DEXP(-2.30256*arg))
c
        arg  = 0.0
        arg  = 1.20145E+01
     &        -5.66908E+00*th1
     &        +2.55583E-01*th2
     &        -3.69367E-02*th3
     &        +1.90539E-03*th4
        g(LiCl) =  DMIN1(1.D+300, DEXP(-2.30256*arg))
c   
        arg  = 0.0
        arg  = 2.51219E+01
     &        -1.06248E+01*th1
     &        +5.03575E-01*th2
     &        -7.21409E-02*th3
     &        +3.71830E-03*th4
        g(LiOH) =  DMIN1(1.D+300, DEXP(-2.30256*arg))


*  Umrechnen der Burrow-FeH-deltaG polynome in das Gail'sche System
*  lnKp = - DeltaG_Burrows/RcalT - lnPstd
*  Achtung: [DeltaG_Burrows]= cal/mol ; [Rcal] = 1.987 cal/(mol K)
*  Pstd = 1atm

        TT1 = DMAX1(tdispol,tg)
        TT2 = TT1*TT1
        TT3 = TT2*TT1
        RcalT = Rcal*TT1
        arg = (-3.87740E+04 
     &         +2.47290E+01*TT1
     &         -4.64016E-04*TT2
     &         +6.79410E-08*TT3)/RcalT
        g(FeH) = DMIN1(1.D+300, DEXP(-arg)/atm)
c Sharp & Huebner 1990
        arg = (-3.04561E+05/TT1 
     &         -4.91489E+04 
     &         +2.32539E+01*TT1 
     &         -1.56846E-04*TT2 
     &         +4.53896E-08*TT3)/RcalT
        g(TiH) = DMIN1(1.D+300, DEXP(-arg)/atm)
c
        arg = ( 3.07543E+05/TT1 
     &         -1.00087E+05 
     &         +2.36342E+01*TT1 
     &         +1.24759E-05*TT2 
     &         +1.23522E-08*TT3)/RcalT
        g(NaCl) = DMIN1(1.D+300, DEXP(-arg)/atm)
c
        arg = ( 4.57163e+05/TT1 
     &         -1.04386E+05 
     &         +2.38671E+01*TT1 
     &         -3.71773E-04*TT2 
     &         +6.22590E-08*TT3)/RcalT
        g(KCl) = DMIN1(1.D+300, DEXP(-arg)/atm)
c
        arg = ( 2.89596e+05/TT1 
     &         -9.88707E+04 
     &         +2.13234E+01*TT1 
     &         -1.27956E-04*TT2 
     &         +4.49210E-08*TT3)/RcalT
        g(CaCl) = DMIN1(1.D+300, DEXP(-arg)/atm)
c    
c        print*,'g(CaCl)', CaCl, g(caCl)

        arg = ( 2.75428e+05/TT1 
     &         -2.15484E+05 
     &         +4.91645E+01*TT1 
     &         -4.41153E-04*TT2 
     &         +7.17853E-08*TT3)/RcalT
        g(CaCl2) = DMIN1(1.D+300, DEXP(-arg)/atm)
c
        arg = ( 4.30684e+05/TT1 
     &         -1.06291E+05 
     &         +2.61097e+01*TT1 
     &         +4.66915E-04*TT2 
     &         -2.90088E-08*TT3)/RcalT
        g(HCl) = DMIN1(1.D+300, DEXP(-arg)/atm)


*---------------------------------------------------------------------------
*     ! Elektronendichte
*     ================== 
      ng = anHges
      Sa = gk(HII)*kT1
      nelek = ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(C)
      Sa = gk(CII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(N)
      Sa = gk(NII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(O)
      Sa = gk(OII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(S)
      Sa = gk(SII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Si)
      Sa = gk(SiII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Fe)
      Sa = gk(FeII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Mg)
      Sa = gk(MgII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Al)
      Sa = gk(AlII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Ti)
      Sa = gk(TiII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Na)
      Sa = gk(NaII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(K)
      Sa = gk(KII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      ng = anHges * eps(Ca)
      Sa = gk(CaII)*kT1
      nelek = nelek + ng/(0.5D0 + DSQRT(0.25D0 + ng/Sa))

      anmono(el) = nelek
      pel = nelek*kT
*
*-----------------------------------------------------------------------
*
*     ! Gleichgewicht zwischen H, H+, H2
*     ==================================
      g(HII)      = gk( HII )
      g(H2)       = gk( H2 )
      pHges       = anHges * kT
      ppp         = ( 1.d0 + g(HII)/pel )  / (2.d0*g(H2))
      qqq         = -pHges / (2.d0*g(H2))
      pH          = vieta(ppp,qqq)
      anmono(H)   = pH / kT
      anmol(H2)   = g(H2) * pH**2 / kT
*
*     ! Starterte He: atomar
*     ======================
      anmono(He) = anHges * eps(He) 
*
*-----------------------------------------------------------------------
*
*  Unterscheidung zwischen O-reicher und C-reicher Chemie
*  ======================================================
*
      if (eps(O).gt.eps(C)) then
*
*     ! sauerstoffreicher Fall
*     ========================
*              
*       ! Gleichgewicht  O, O+, C, C+, CO, H2O, CH4
*       ===========================================
        g(OII) = gk( OII )
        g(CII) = gk( CII )
        g(CH4) = gk( CH4 )
        g(H2O) = gk( H2O )
        g(CO)  = gk( CO  )
        pOges  = eps(O) * anHges * kT
        pCges  = eps(C) * anHges * kT
        aa = 1.d0 + g(OII)/pel + pH**2*g(H2O)
        bb = 1.d0 + g(CII)/pel + pH**4*g(CH4)
        ppp = ( (pOges-pCges)*g(CO) + aa*bb) / (g(CO)*bb)
        qqq = -pCges*aa / (g(CO)*bb)
        pC  = VIETA(ppp,qqq)

*       ! Gleichgewicht  O, O+, CO, CO2, H2O
*       ====================================
        g(CO2)  = gk( CO2  )
        g(C2H2) = gk( C2H2 )
        aa  = 2.d0*pC*g(CO2)
        ppp = (1.d0 + g(OII)/pel + pC*g(CO) + pH**2*g(H2O)) /aa
        qqq = -pOges/aa
        pO  = VIETA(ppp,qqq)

*       ! pC/pO-Nachiteration  C,C+,O,O+,CO,CO2,H2O,CH4,C2,C2H,C2H2
*       ===========================================================
        g(C2)  = gk( C2  )
        g(C2H) = gk( C2H )
        aa = 2.d0 * ( g(C2) + pH*g(C2H) + pH**2*g(C2H2) )
        bb = 1.d0 + g(CII)/pel + pH**4*g(CH4)
        cc = g(CO)
        dd = g(CO2)
        ee = 1.d0 + g(OII)/pel + pH**2*g(H2O)
        pCpOit = 0
 900    continue
          func(1) = pC**2*aa + pC*bb + pC*pO*cc + pC*pO**2*dd - pCges
          dfunc(1,1) = 2.d0*pC*aa + bb + pO*cc + pO**2*dd 
          dfunc(1,2) = pC*cc + 2.d0*pC*pO*dd 
          func(2) = pO*ee + pC*pO*cc + 2.d0*pC*pO**2*dd - pOges
          dfunc(2,1) = pO*cc + 2.d0*pO**2*dd 
          dfunc(2,2) = ee + pC*cc + 4.d0*pC*pO*dd 
          call GAUSS(2,dfunc,dpp,func)
          pCalt = pC
          pOalt = pO
          fak = 1.0+exp(-0.2*MAX(0,pCpOit-10))
          pC  = DMAX1(DMIN1(pC-dpp(1),pC*fak),pC/fak)
          pO  = DMAX1(DMIN1(pO-dpp(2),pO*fak),pO/fak)
          delta  = DMAX1(DABS(pCalt/pC-1.d0),DABS(pOalt/pO-1.d0))          
          pCpOit = pCpOit + 1
c         write(*,*) 'pC/pO-Iteration:',pCpOit,pC*kT1,pO*kT1,delta
        if ((delta.gt.1.d-8).and.(pCpOit.lt.10)) goto 900 
        anmono(C)   = pC / kT
        anmono(O)   = pO / kT
        anmol(CO)   = g(CO)   * pC * pO / kT
        anmol(CO2)  = g(CO2)  * pC * pO**2 / kT
        anmol(H2O)  = g(H2O)  * pH**2 * pO / kT
        anmol(CH4)  = g(CH4)  * pC * pH**4 / kT
        anmol(C2)   = g(C2)   * pC**2 / kT
        anmol(C2H)  = g(C2H)  * pC**2 * pH / kT
        anmol(C2H2) = g(C2H2) * pC**2 * pH**2 / kT
*
*       ! Gleichgewicht N, N2, N+, NH3
*       ==============================
        g(NII)      = gk( NII )
        g(N2)       = gk( N2 )
        g(NH3)      = gk( NH3 )
        pNges       = eps(N) * anHges * kT
        ppp         = (1.d0 + g(NII)/pel + pH**3*g(NH3)) / (2.d0*g(N2))
        qqq         = -pNges / (2.d0*g(N2))
        pN          = vieta(ppp,qqq)
        anmono(N)   = pN / kT
        anmol(N2)   = g(N2)  * pN**2 / kT 
        anmol(NH3)  = g(NH3) * pN * pH**3 / kT 
*
*       ! Gleichgewicht  Si,SiH,SiH4,SiO,SiO2,SiN,SiS,S,H2S,HS,SO,SO2,Si+,S+
*       ====================================================================
        g(SiII)     = gk( SiII )
        g(SII)      = gk( SII  )
        g(SiH)      = gk( SiH  )
        g(SiH4)     = gk( SiH4 )
        g(SiO)      = gk( SiO  )
        g(SiO2)     = gk( SiO2 )
        g(SiN)      = gk( SiN  )
        g(SiS)      = gk( SiS  )
        g(H2S)      = gk( H2S  ) 
        g(HS)       = gk( HS   )
        g(SO)       = gk( SO   )
        g(SO2)      = gk( SO2  )
        pSiges      = eps(Si) * anHges * kT
        pSges       = eps( S) * anHges * kT
        aa          = 1.d0 + g(SiH)*pH + g(SiH4)*pH**4 + g(SiO)*pO
     &                + g(SiO2)*pO**2 + g(SiN)*pN + g(SiII)/pel
        bb          = 1.d0 + g(SO)*pO + g(HS)*pH + g(H2S)*pH**2
     &                + g(SO2)*pO**2 + g(SII)/pel
        ppp         = aa/g(SiS) + (pSiges-pSges)/bb
        qqq         = -pSges * aa / bb / g(SiS)
        pS          = vieta(ppp,qqq)
        pSi         = pSiges / (aa + pS * g(SiS))
        anmono(Si)  = pSi / kT
        anmono(S)   = pS  / kT
        anmol(H2S)  = g(H2S)  * pH**2 * pS / kT
        anmol(SiS)  = g(SiS)  * pSi   * pS / kT
        anmol(SO)   = g(SO)   * pS    * pO / kT
        anmol(HS)   = g(HS)   * pH    * pS / kT
        anmol(SiO)  = g(SiO)  * pSi   * pO / kT
        anmol(SiH)  = g(SiH)  * pSi   * pH / kT
        anmol(SiH4) = g(SiH4) * pSi   * pH**4 / kT
        anmol(SiO2) = g(SiO2) * pSi   * pO**2 / kT
        anmol(SiN)  = g(SiN)  * pSi   * pN / kT
        anmol(SO2)  = g(SO2)  * pS    * pO**2 / kT
*        
*       !  Gleichgewicht Mg, MgH, MgO, MgOH, Mg(OH)2, MgS, Mg+
*       ======================================================
        g(MgII)     = gk( MgII )
        g(MgO)      = gk( MgO  )
        g(MgOH)     = gk( MgOH )
        g(MgO2H2)   = gk(MgO2H2)
        g(MgH)      = gk( MgH  )
        g(MgS)      = gk( MgS  )
        pMgges      = eps(8) * anHges * kT
        pMg         = pMgges / ( 1.d0 + g(MgH)*pH + g(MgO)*pO 
     &        + g(MgOH)*pH*pO + g(MgO2H2)*pO**2*pH**2
     &        + g(MgS)*pS + g(MgII)/pel )
        anmono(Mg)  = pMg / kT
        anmol(MgS)  = g(MgS)  * pMg * pS / kT
        anmol(MgO)  = g(MgO)  * pMg * pO / kT
        anmol(MgOH) = g(MgOH) * pMg * pO * pH / kT
        anmol(MgH)  = g(MgH)  * pMg * pH / kT
*
*       ! Gleichgewicht  Fe , FeO , FeS, FeH, Fe(OH)2, Fe+
*       ==================================================
c       g(FeH)        : siehe oben!
        g(FeII)       = gk( FeII   )
        g(FeO)        = gk( FeO    )
        g(FeS)        = gk( FeS    )
        g(FeO2H2)     = gk( FeO2H2 )
        pFeges        = eps(Fe) * anHges * kT
        pFe           = pFeges / ( 1.d0 + pO*g(FeO) + g(FeII)/pel 
     &                  + pH**2*pO**2*g(FeO2H2) + pH*g(FeH) + pS*g(FeS))
        anmono(Fe)    = pFe / kT
        anmol(FeO)    = g(FeO)    * pFe * pO / kT 
        anmol(FeS)    = g(FeS)    * pFe * pS / kT 
        anmol(FeH)    = g(FeH)    * pFe * pH / kT 
        anmol(FeO2H2) = g(FeO2H2) * pFe * pO**2 * pH**2 / kT 
*
*       ! Gleichgewicht Al , AlOH , AlO2H, Al2O, AlH, Al+
*       =================================================
        g(AlII)      = gk( AlII  )
        g(AlH)       = gk( AlH   )
        g(AlOH)      = gk( AlOH  )
        g(AlO2H)     = gk( AlO2H )
        g(Al2O)      = gk( Al2O  )
        pAlges       = eps(Al) * anHges * KT
        ppp          = 1.d0 + pO*pH*g(AlOH) + pO**2*pH*g(AlO2H) 
     &                 + pH*g(AlH) + g(AlII)/pel
        ppp          = ppp / pO / g(Al2O) / 2.d0
        qqq          = -pAlges / pO / g(Al2O) / 2.d0
        pAl          = vieta(ppp,qqq)
        anmono(Al)   = pAl / kT
        anmol(AlOH)  = g(AlOH)  * pAl * pO * pH / kT
        anmol(AlO2H) = g(AlO2H) * pAl * pO**2 * pH / kT
        anmol(Al2O)  = g(Al2O)  * pAl**2 * pO / kT
        anmol(AlH)   = g(AlH)   * pAl * pH / kT
*
*       ! Gleichgewicht Ti, TiH, TiC, TiC2, TiS, Ti+, TiO, TiO2
*       =======================================================
c       g(TiC)  : siehe oben!
c       g(TiH)  : siehe oben!
        g(TiII) = gk( TiII )        
        g(TiC2) = gk( TiC2 )
        g(TiS)  = gk( TiS  )
        g(TiO)  = gk( TiO  )
        g(TiO2) = gk( TiO2 )
	pTiges  = eps(Ti) * anHges * kT    
        pTi     = pTiges / ( 1.D0 + pS*g(TiS) + pH*g(TiH) 
     &                       + pC**2*g(TiC2) + pC*g(TiC) + g(TiII)/pel 
     &                       + pO*g(TiO) + pO**2*g(TiO2) )
        anmono(Ti)  = pTi / kT
        anmol(TiH)  = g(TiH)  * pTi * pH / kT
        anmol(TiC)  = g(TiC)  * pTi * pC / kT
        anmol(TiC2) = g(TiC2) * pTi * pC**2 / kT
        anmol(TiS)  = g(TiS)  * pTi * pS / kT
        anmol(TiO)  = g(TiO)  * pTi * pO / kT 
        anmol(TiO2) = g(TiO2) * pTi * pO**2 / kT 
*
*       ! Gleichgewicht  Na , Na+, NaOH , (NaOH)2, NaH
*       ===============================================
        g(NaII)     = gk( NaII )
c       g(NaH)      : siehe oben!
        g(NaOH)     = gk( NaOH )
        g(Na2O2H2)  = gk( Na2O2H2 )
        pNages      = eps(Na) * anHges * kT
        aa  = 2.d0*pO**2*pH**2*g(Na2O2H2)
        ppp = ( 1.d0 + g(NaII)/pel + pO*pH*g(NaOH) + pH*g(NaH)) /aa 
        qqq = -pNages/aa
        pNa = VIETA(ppp,qqq)
        anmono(Na)  = pNa / kT
*
*       ! Gleichgewicht  K , K+, KOH , (KOH)2
*       =====================================
        g(KII)      = gk( KII )
        g(KOH)      = gk( KOH )
        g(K2O2H2)   = gk( K2O2H2 )
        pKges       = eps(K) * anHges * kT
        aa  = 2.d0*pO**2*pH**2*g(K2O2H2)
        ppp = ( 1.d0 + g(KII)/pel + pO*pH*g(KOH) ) /aa 
        qqq = -pKges/aa
        pK  = VIETA(ppp,qqq)
        anmono(K)  = pK / kT
*
*       ! Gleichgewicht  Ca , Ca+, CaH, CaOH , Ca(OH)2
*       ==============================================
c       g(CaH)      : siehe oben!
        g(CaII)     = gk( CaII )
        g(CaOH)     = gk( CaOH )
        g(CaH)      = gk( CaH )
        g(CaO2H2)   = gk( CaO2H2 )
        pCages      = eps(Ca) * anHges * kT
        pCa         = pCages / ( 1.d0 + g(CaII)/pel + pO*pH*g(CaOH) 
     &                + pH*g(CaH) + pO**2*pH**2*g(CaO2H2) ) 
        anmono(Ca)  = pCa / kT
*
*       ! Gleichgewicht  Cl, KCl, NaCl, CaCl, HCl
*       ================================================
c        g(KCl)     : siehe oben!
c        g(HCl)     : siehe oben!
c        g(NaCl)    : siehe oben!
c        g(CaCl)    : siehe oben!
        pClges     = eps(Cl) * anHges * kT
        pCl        = pClges / ( 1.D0 + pK*g(KCl) + pNa*g(NaCl) 
     &             + pH*g(HCl) + pCa*g(CaCl) )
        anmono(Cl) = pCl/kT 
*
*       ! Gleichgewicht  Li , LiO, LiH, LiOH, LiCl
*       ==========================================
cc        g(LiII)    = gk( LiII )
c        g(LiH)     : siehe oben!
c        g(LiO)     : siehe oben!
c        g(LiCl)    : siehe oben!
c        g(LiOH)    : siehe oben!
        pLiges     = eps(Li) * anHges * kT
        pLi        = pLiges / ( 1.D0 +  pH*g(LiH) + pCl*g(LiCl) 
     &             + pO*g(LiO)  + pO*pH*g(LiOH))
        anmono(Li)  = pLi / kT
*
      else 
*
*     ! kohlenstoffreicher Fall
*     =========================
*              
*       ! Gleichgewicht C, C+, O, O+, CO, H2O, CH4
*       ==========================================
        g(OII) = gk( OII )
        g(CII) = gk( CII )
        g(CO)  = gk( CO  )
        g(CH4) = gk( CH4 )
        g(H2O) = gk( H2O )
        pOges  = eps(O) * anHges * kT
        pCges  = eps(C) * anHges * kT
        aa  = 1.d0 + g(OII)/pel + pH**2*g(H2O)
        bb  = 1.d0 + g(CII)/pel + pH**4*g(CH4)
        ppp = ( (pCges-pOges)*g(CO) + aa*bb ) / (g(CO)*aa)
        qqq = -pOges*bb / (g(CO)*aa)
        pO  = VIETA(ppp,qqq)
*
*       ! Gleichgewicht C, C+, CO, C2 , C2H , C2H2 , CH4
*       ================================================
        g(C2)   = gk( C2   )
        g(C2H)  = gk( C2H  )
        g(C2H2) = gk( C2H2 )
        g(CH4)  = gk( CH4  )
        g(CO2)  = gk( CO2  )
        aa  = 2.d0 * ( pH**2*g(C2H2) + pH*g(C2H) + g(C2) )
        ppp = ( 1.d0 + g(CII)/pel + pO*g(CO) + pH**4*g(CH4) ) / aa
        qqq = -pCges / aa
        pC  = VIETA(ppp,qqq)
        aa  = 2.d0*pC*g(CO2)
        ppp = (1.d0 + g(OII)/pel + pC*g(CO) + pH**2*g(H2O) ) / aa
        qqq = -pOges/aa
        pO  = VIETA(ppp,qqq)
*       
*       ! pC/pO-Nachiteration C,C+,O,O+,CO,CO2,H2O,CH4,C2H2,C2H,C2,C3,C3H
        !================================================================
        g(C3)  = gk( C3  )
        g(C3H) = gk( C3H )
        a3 = 3.d0 * ( g(C3) + pH*g(C3H) )
        aa = 2.d0 * ( g(C2) + pH*g(C2H) + pH**2*g(C2H2) )
        bb = 1.d0 + g(CII)/pel + pH**4*g(CH4)
        cc = g(CO)
        dd = g(CO2)
        ee = 1.d0 + g(OII)/pel + pH**2*g(H2O)
        pCpOit = 0
        do
          func(1)    = pC**3*a3 + pC**2*aa + pC*bb + pC*pO*cc 
     &               + pC*pO**2*dd - pCges
          dfunc(1,1) = 3.d0*pC**2*a3 + 2.d0*pC*aa + bb + pO*cc 
     &               + pO**2*dd 
          dfunc(1,2) = pC*cc + 2.d0*pC*pO*dd 
          func(2)    = pO*ee + pC*pO*cc + 2.d0*pC*pO**2*dd - pOges
          dfunc(2,1) = pO*cc + 2.d0*pO**2*dd 
          dfunc(2,2) = ee + pC*cc + 4.d0*pC*pO*dd 
          call GAUSS(2,dfunc,dpp,func)
          pCalt = pC
          pOalt = pO
          fak = 1.0+exp(-0.2*MAX(0,pCpOit-10))
          pC  = DMAX1(DMIN1(pC-dpp(1),pC*fak),pC/fak)
          pO  = DMAX1(DMIN1(pO-dpp(2),pO*fak),pO/fak)
          delta = DMAX1(DABS(pCalt/pC-1.d0),DABS(pOalt/pO-1.d0))          
          pCpOit = pCpOit + 1
          write(*,'(a16,i3,3(1pE11.4))') 
     &          'pC/pO-Iteration:',pCpOit,pCalt/kT,pOalt/kT,delta
          if ((pCpOit.ge.20).or.(delta.lt.1.d-12)) exit
        enddo  

        anmono(C)   = pC / kT
        anmono(O)   = pO / kT
        anmol(CO)   = g(CO)   * pC * pO / kT
        anmol(CO2)  = g(CO2)  * pC * pO**2 / kT
        anmol(C2)   = g(C2)   * pC**2 / kT
        anmol(C3)   = g(C3)   * pC**3 / kT
        anmol(C2H)  = g(C2H)  * pC**2 * pH / kT
        anmol(C3H)  = g(C3H)  * pC**3 * pH / kT
        anmol(C2H2) = g(C2H2) * pC**2 * pH**2 / kT
        anmol(CH4)  = g(CH4)  * pC * pH**4 / kT
        anmol(H2O)  = g(H2O)  * pH**2 * pO / kT
*
*       ! Gleichgewicht N, N+, N2, CN, HCN, NH3
*       =======================================
        g(NII)     = gk(NII)
        g(N2)      = gk(N2)
        g(CN)      = gk(CN)
        g(HCN)     = gk(HCN)
        g(NH3)     = gk(NH3)
        pNges      = eps(N) * anHges * kT
        ppp        = (1.d0 + g(CN)*pC + g(HCN)*pC*pH 
     &                + g(NH3)*pH**3 + g(NII)/pel)/(2.*g(N2))
        qqq        = -pNges / (2.d0*g(N2))
        pN         = vieta(ppp,qqq)
        anmono(N)  = pN / kT
        anmol(N2)  = g(N2) * pN**2 / kT
        anmol(CN)  = g(CN) * pC * pN / kT
        anmol(HCN) = g(HCN) * pH * pC * pN / kT
        anmol(NH3) = g(NH3) * pN * pH**3 / kT

*       ! Gleichgewicht Si,S,SiS,SiC,SiO,Si2C,SiC2,SiH,SiH4,SiN,CS,HS,H2S
*       =================================================================
        g(SiII)  = gk( SiII )
        g(SII)   = gk( SII )
        g(SiS)   = gk( SiS )
        g(SiC)   = gk( SiC )
        g(SiO)   = gk( SiO )
        g(Si2C)  = gk( Si2C )
        g(SiC2)  = gk( SiC2 )
        g(SiH)   = gk( SiH )
        g(SiH4)  = gk( SiH4 )
        g(SiN)   = gk( SiN )
        g(CS)    = gk( CS )
        g(HS)    = gk( HS )
        g(H2S)   = gk( H2S )
        pSiges   = eps(Si) * anHges * kT
        pSges    = eps( S) * anHges * kT
        aa       = 1.d0 + pC*g(SiC) + pC**2*g(SiC2) + pH*g(SiH) 
     &             + pH**4*g(SiH4) + pN*g(SiN) + pO*g(SiO) + g(SiII)/pel
        bb       = 1.d0 + pC   *g(CS)   + pH*g(HS) + pH**2*g(H2S) 
     &             + g(SII)/pel
        ppp      = aa/g(SiS) + (pSiges-pSges)/bb
        qqq      = -pSges * aa / bb / g(SiS)
        pS       = vieta(ppp,qqq)
        pSi      = pSiges / ( aa + pS * g(SiS) )
*
*       ! Nachiteration wegen Si2C
*       ==========================                
        f3 = 2.d0*g(Si2C)*pC*g(SiS)
        f2 = 2.d0*bb*pc*g(Si2C)+aa*g(SiS)
        f1 = aa*bb+g(SiS)*(pSges-pSiges)
        f0 = -bb*pSiges
 1000   continue
          ff3   = f3*pSi**3
          ff2   = f2*pSi**2
          ff1   = f1*pSi
          f     = ff3+ff2+ff1+f0
          fs    = 3.d0*ff3 + 2.d0*ff2+ff1
          delta = -f/fs
          pSi   = pSi * exp(delta)
        if ( DABS(delta) .gt. 1.d-08 ) goto 1000
        pS = psges / ( bb + g(SiS)*pSi )
        anmono(Si)  = pSi / kT
        anmono(S)   = pS / kT
        anmol(Si2C) = g(Si2C) * pSi**2 * pC    / kT
        anmol(SiC2) = g(SiC2) * pSi    * pC**2 / kT
        anmol(SiH)  = g(SiH)  * pSi    * pH    / kT
        anmol(CS)   = g(CS)   * pS     * pC    / kT
        anmol(HS)   = g(HS)   * pS     * pH    / kT
        anmol(H2S)  = g(H2S)  * pS     * pH**2 / kT
        anmol(SiS)  = g(SiS)  * pSi    * pS    / kT
        anmol(SiC)  = g(SiC)  * pSi    * pC    / kT
        anmol(SiH4) = g(SiH4) * pSi    * pH**4 / kT
        anmol(SiN)  = g(SiN)  * pSi    * pN    / kT
*
*       ! Gleichgewicht Mg, MgH, MgS, Mg+, Mg(OH)2
*       ==========================================
        g(MgII)    = gk(MgII)
        g(MgH)     = gk(MgH)
        g(MgS)     = gk(MgS)
        g(MgO2H2)  = gk(MgO2H2)
	pMgges     = eps(Mg) * anHges * kT
        pMg        = pMgges / (1.d0 + g(MgH)*pH + g(MgS)*pS 
     &             + g(MgII)/pel + g(MgO2H2)*pO**2*pH**2)
        anmono(Mg) = pMg / kT
        anmol(MgH) = g(MgH) * pMg * pH / kT
        anmol(MgS) = g(MgS) * pMg * pS / kT
*
*       ! Gleichgewicht Fe, FeH, FeS, Fe+, Fe(OH)2
*       ==========================================
c       g(FeH)     : siehe oben!
        g(FeII)    = gk(FeII)
        g(FeS)     = gk(FeS)
        g(FeO2H2)  = gk(FeO2H2)
	pFeges     = eps(Fe) * anHges * kT    
        pFe        = pFeges / ( 1.d0 + pS*g(FeS) + g(FeII)/pel 
     &             + g(FeO2H2)*pO**2*pH**2 + g(FeH)*pH)
        anmono(Fe)    = pFe / kT
        anmol(FeS)    = g(FeS) * pFe * pS / kT
        anmol(FeH)    = g(FeH) * pFe * pH / kT
        anmol(FeO2H2) = g(FeO2H2) * pFe * pH**2 * pO**2 / kT
*
*       ! Gleichgewicht Ti, TiC, TiC2, TiS, Ti+, TiO, TiO2
*       ==================================================
c       g(TiC)  : siehe oben!
c       g(TiH)  : siehe oben!
        g(TiII) = gk( TiII )        
        g(TiC2) = gk( TiC2 )
        g(TiS)  = gk( TiS  )
        g(TiO)  = gk( TiO  )
        g(TiO2) = gk( TiO2 )
	pTiges  = eps(Ti) * anHges * kT    
        pTi     = pTiges / ( 1.d0 + pS*g(TiS) + pC**2*g(TiC2)
     &          + pH*g(TiH) + pC*g(TiC) + g(TiII)/pel 
     &          + pO*g(TiO) + pO**2*g(TiO2) )
        anmono(Ti)  = pTi / kT
        anmol(TiC)  = g(TiC)  * pTi * pC / kT
        anmol(TiC2) = g(TiC2) * pTi * pC**2 / kT
        anmol(TiS)  = g(TiS)  * pTi * pS / kT
        anmol(TiO)  = g(TiO)  * pTi * pO / kT 
        anmol(TiO2) = g(TiO2) * pTi * pO**2 / kT 
*
*       ! Gleichgewicht Al , AlOH , AlO2H, Al2O, AlH, Al+
*       =================================================
        g(AlII)      = gk( AlII  )
        g(AlH)       = gk( AlH   )
        g(AlOH)      = gk( AlOH  )
        g(AlO2H)     = gk( AlO2H )
        g(Al2O)      = gk( Al2O  )
        pAlges       = eps(Al) * anHges * kT
        ppp          = 1.d0 + pO*pH*g(AlOH) + pO**2*pH*g(AlO2H) 
     &                      + pH*g(AlH) + g(AlII)/pel
        ppp          = ppp / pO / g(Al2O) / 2.d0
        qqq          = -pAlges / pO / g(Al2O) / 2.d0
        pAl          = vieta(ppp,qqq)
        anmono(Al)   = pAl / kT
        anmol(AlOH)  = g(AlOH)  * pAl    * pO    * pH / kT
        anmol(AlO2H) = g(AlO2H) * pAl    * pO**2 * pH / kT
        anmol(Al2O)  = g(Al2O)  * pAl**2 * pO / kT
        anmol(AlH)   = g(AlH)   * pAl    * pH / kT
*
*       ! Gleichgewicht  Na , Na+, NaOH , (NaOH)2
*       =========================================
        g(NaII)     = gk( NaII )
        g(NaOH)     = gk( NaOH )
        g(Na2O2H2)  = gk( Na2O2H2 )
        pNages      = eps(Na) * anHges * kT
        aa  = 2.d0*pO**2*pH**2*g(Na2O2H2)
        ppp = ( 1.d0 + g(NaII)/pel + pO*pH*g(NaOH) ) /aa 
        qqq = -pNages/aa
        pNa = VIETA(ppp,qqq)
        anmono(Na) = pNa / kT
*
*       ! Gleichgewicht  K , K+, KOH , (KOH)2
*       =====================================
        g(KII)      = gk( KII )
        g(KOH)      = gk( KOH )
        g(K2O2H2)   = gk( K2O2H2 )
        pKges       = eps(K) * anHges * kT
        aa  = 2.d0*pO**2*pH**2*g(K2O2H2)
        ppp = ( 1.d0 + g(KII)/pel + pO*pH*g(KOH) ) /aa 
        qqq = -pKges/aa
        pK  = VIETA(ppp,qqq)
        anmono(K) = pK / kT
*
*       ! Gleichgewicht  Ca , Ca+, CaH, CaOH , Ca(OH)2
*       ==============================================
c       g(CaH)      : siehe oben!
        g(CaII)     = gk( CaII )
        g(CaOH)     = gk( CaOH )
        g(CaO2H2)   = gk( CaO2H2 )
        pCages      = eps(Ca) * anHges * kT
        pCa         = pCages / ( 1.d0 + g(CaII)/pel + pO*pH*g(CaOH) 
     &              + pH*g(CaH) + pO**2*pH**2*g(CaO2H2) ) 
        anmono(Ca) = pCa / kT
*
*       ! Gleichgewicht  Cl, KCl, NaCl, CaCl, HCl
*       ================================================
c        g(KCl)     : siehe oben!
c        g(HCl)     : siehe oben!
c        g(NaCl)    : siehe oben!
c        g(CaCl)    : siehe oben!
        pClges     = eps(Cl) * anHges * kT
        pCl        = pClges / ( 1.D0 + pK*g(KCl) + pNa*g(NaCl) 
     &             + pH*g(HCl) + pCa*g(CaCl))
        anmono(Cl) = pCl/kT 
*
*       ! Gleichgewicht  Li , LiO, LiH, LiOH, LiCl
*       ==========================================
cc        g(LiII)    = gk( LiII )
c        g(LiH)     : siehe oben!
c        g(LiO)     : siehe oben!
c        g(LiCl)    : siehe oben!
c        g(LiOH)    : siehe oben!
        pLiges     = eps(Li) * anHges * kT
        pLi        = pLiges / ( 1.D0 +  pH*g(LiH) + pCl*g(LiCl) 
     &             + pO*g(LiO)  + pO*pH*g(LiOH))
        anmono(Li)  = pLi / kT
*
      endif
*
*     ! nochmal Elektronendichte
*     ==========================
      do i=1,2
        anmol(HII)  = g(HII)  * pH /pel / kT
        anmol(CII)  = g(CII)  * pC /pel / kT
        anmol(OII)  = g(OII)  * pO /pel / kT
        anmol(NII)  = g(NII)  * pN /pel / kT
        anmol(SII)  = g(SII)  * pS /pel / kT
        anmol(SiII) = g(SiII) * pSi/pel / kT
        anmol(MgII) = g(MgII) * pMg/pel / kT
        anmol(FeII) = g(FeII) * pFe/pel / kT
        anmol(AlII) = g(AlII) * pAl/pel / kT
        anmol(TiII) = g(TiII) * pTi/pel / kT
        anmol(NaII) = g(NaII) * pNa/pel / kT
        anmol(KII)  = g(KII)  * pK /pel / kT
        anmol(CaII) = g(CaII) * pCa/pel / kT
        anmono(el)  = anmol(HII) + anmol(CII) + anmol(OII) + anmol(NII)  
     &       + anmol(SII) + anmol(SiII)+ anmol(MgII)+ anmol(FeII) 
     &       + anmol(AlII)+ anmol(TiII)+ anmol(NaII)+ anmol(KII)  
     &       + anmol(CaII)
        pel = anmono(el)*kT
      enddo
*     
*-----------------------------------------------------------------------
*
      if ( alle ) then
        ! alle Molekuele mitrechnen
*       ===========================
        do i=1,nmol
          if ((i.ne.TiC).and.(i.ne.FeH).and.(i.ne.CaH).and.
     &        (i.ne.LiH).and.(i.ne.LiO).and.(i.ne.LiOH).and.
     &        (i.ne.LiCl).and.(i.ne.NaCl).and.(i.ne.KCl).and.
     &        (i.ne.CaCl).and.(i.ne.CaCl2).and.(i.ne.HCl).and.
     &        (i.ne.TiH) )
     &         g(i) = gk(i)
c               print*, 'i, g(i)', i, g(i)
          pmol = g(i)
          do j=1,m_kind(0,i)
c           pmol = pmol * (anmono(m_kind(j,i))*kT) ** m_anz(j,i)
            pat = anmono(m_kind(j,i))*kT
            if (m_anz(j,i).gt.0) then
              do kk=1,m_anz(j,i)
                pmol = pmol*pat
              enddo
            else
              do kk=1,-m_anz(j,i)
                pmol = pmol/pat
              enddo
            endif
          enddo
  	  anmol(i) = pmol*kT1
	enddo

        if (nachit) then
          ! Iteration nach Newton-Raphson
*         ===============================
          if ((ilauf.gt.1).and.merk) then
c           write(*,*) 'benutze Konzentrationen von vorher'
            do i=1,nel
	      anmono(i) = amerk(i) * anhges
            enddo
	  endif
          it = 0
          abbruch = .false.
 1500     continue
*

          do i=1,nel
            FF(i) = anHges*eps(i)*kT - anmono(i)*kT
            do j=1,nel
	      DF(i,j) = 0.d0
            enddo
	    DF(i,i) = -1.d0
            pmono1(i) = 1.d0 / (anmono(i)*kT)
          enddo	
*
          do i=1,nmol
            pmol = g(i)
            do j=1,m_kind(0,i)
c             pmol = pmol * (anmono(m_kind(j,i))*kT) ** m_anz(j,i)
              pat = anmono(m_kind(j,i))*kT
              if (m_anz(j,i).gt.0) then
                do kk=1,m_anz(j,i)
                  pmol = pmol*pat
                enddo
              else
                do kk=1,-m_anz(j,i)
                  pmol = pmol/pat
                enddo
              endif
            enddo
  	    anmol(i) = pmol*kT1
            do j=1,m_kind(0,i)
	      m1     = m_kind(j,i) 
              term   = m_anz(j,i) * pmol
              FF(m1) = FF(m1) - term
              do l=1,m_kind(0,i)
                m2        = m_kind(l,i) 
                DF(m1,m2) = DF(m1,m2) - m_anz(l,i) * term * pmono1(m2)
	      enddo	    
	    enddo
          enddo
*	  
          call GAUSS (nel,DF,dp,FF)
	  delta = 0.d0 
          fak = 1.d0 + 4.d0*DEXP(-DBLE(MAX(0,it-10))/10.d0)
c         if (it.eq.0) then
c           do i=1,nel
c             write(*,2000) catm(i),anHges*eps(i),
c    &                      anmono(i),-dp(i)*pmono1(i)
c           enddo
c         endif
          do i=1,nel
	    delta     = DMAX1(delta, DABS(dp(i)*pmono1(i)))
            anmono(i) = DMAX1( DMIN1( anmono(i)-dp(i)*kT1, 
     &                  anmono(i)*fak), anmono(i)/fak)
            if (eps(i).gt.1.d-30) then
              anmono(i) = DMIN1(anmono(i),eps(i)*anhges)
            endif
            amerk(i)  = anmono(i)/anhges
          enddo
          it = it + 1
          if (it.gt.200) then
            write(*,*) '*** keine Konvergenz in SMCHEM!'
            write(*,*) 'it, delta =',it,delta
            write(*,*) '  n<H>, T =',anhges,Tg
            do i=1,nel
              write(*,*) catm(i),eps(i),dp(i)*pmono1(i)
            enddo  
            abbruch = .true.
            stop
          endif
          if ((.not.(delta.lt.1.d-13)).and.(.not.abbruch)) goto 1500

c         write(*,*) "SMCHEM: ",it,delta
c         write(*,*) tg,it,' Iterationen'
c         do i=1,nel
c           write(*,2000) catm(i),anHges*eps(i),
c    &                    anmono(i),-dp(i)*pmono1(i)
c         enddo
c         read(*,3000) char

          ! final anmol determination
*         ===========================
          do i=1,nmol
            pmol = g(i)
            do j=1,m_kind(0,i)
              pat = anmono(m_kind(j,i))*kT
              if (m_anz(j,i).gt.0) then
                do kk=1,m_anz(j,i)
                  pmol = pmol*pat
                enddo
              else
                do kk=1,-m_anz(j,i)
                  pmol = pmol/pat
                enddo
              endif
            enddo
            anmol(i) = pmol*kT1
          enddo

        endif
      endif      
      
      if ( ngestst ) then
        ! Test auf Elementerhaltung
*       ===========================
        do i = 1 , nel
          nges(i)  = anmono(i)
        enddo
        do i = 1 , nmol
          do j = 1 , m_kind(0,i)
            j1 = m_kind(j,i)
            nges(j1) = nges(j1) + m_anz(j,i)*anmol(i)
          enddo
        enddo
	write(*,*) 'Elementerhaltung ...'
        do i = 1 , nel
          soll  = anHges * eps(i)
          haben = nges(i)
          abw   = DABS(soll-haben)/DMIN1(soll,haben)
	  write(6,2004) catm(i),abw          
          if ( abw .gt. 0.1 ) then
            write(6,2005) catm(i)
          endif
        enddo
 2000   format(A7,99(1pE14.6))
 2004   format(1x,'Element: ',a10,1pe13.5)
 2005   format(1x,'Elementerhaltung(',a10,') staerker als 10% falsch!')
 3000   format(a1)
      endif
      RETURN
      end
