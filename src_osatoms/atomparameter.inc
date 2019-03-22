C  INCLUDE FILE FOR MARCS   version 24-Aug-01  UGJ
C
C Number of depth points in the model.
C
      PARAMETER ( NDP=500 )
C
C Maximal number of rays for the solution of 
C the transfer equation.
C   This should be at least NDP/2 +3.
C
      PARAMETER ( NRAYS=300 ) 
*
*  equation of state parameters   (for tabs and absko)
*
*  mkomp  = max # of absorption coefficient contributors
*  mkompr = max # of temp. indep. absorption coefficient contributors
*
*     integer mkomp,mkompr,ifadim,kfadim
      parameter (mkomp=27,mkompr=16,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp)
C
C
C Number of wavelength points in the  OS calculation
C
      INTEGER NWL
      PARAMETER (NWL=300000)
C
C  for call to INTERP from OSTABLOOK
C
      integer mtemp,mpe
      parameter (mtemp=12,mpe=1)
C
* include file for testing Tsuji's molecular equilibrium
* maxmol is maximum number of molecules to be stored in 'takemolec'
* maxmet is maximum number of atoms
* npoint is the number of copmuting temperatures
      integer maxmol,maxmet
      parameter (maxmol=350, maxmet=50)
C
* maxosmol is the maximum number of molecules for which opacity is included
* (including one number for the atomic lines)
      integer maxosmol
      parameter (maxosmol=15)

C
      integer natms,nspec,nions
      parameter(natms=49,nions=127,nspec=892)       ! # elements and species
* natms = number of atoms considered in Gibbs energy minimalization routines
* nspec = total number of species (atoms, ions, molecules) considered in Gibbs 
*         energy minimalization routines (= 0.5 * number of lines in 
*         /home/ast/uffegj/falkesg/gfits.data

* additional parameters for atoms.f (be careful with ntemp etc, it should 
* be changed consistently in os.f to avoid completely parameter statements in os.f)
      parameter(ntemp=12,mt1=ntemp+1,nso=6000)
*     nso = max number of lines to be sorted after lineintensity
