C © 2022, NTNU
C Author: Mojtaba Mokhtari <mojtaba.mokhtari@ntnu.no>      
C This code is licenced under EUROPEAN UNION PUBLIC LICENCE v. 1.2
C Parameters are in MPa
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElem, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
c
      include 'vaba_param.inc'
c
      dimension jElem(nblock), coordMp(nblock,*), 
     *          direct(nblock,3,3), T(nblock,3,3), 
     *          charLength(nblock), props(nprops), 
     *          stateOld(nblock,nstatev), 
     *          stateNew(nblock,nstatev),
     *          field(nblock,nfieldv)
      character*80 cmname
c
c     Local arrays from vgetvrm are dimensioned to 
c     maximum block size (maxblk)
c
      parameter( nrData=6 )
      character*3 cData(maxblk*nrData)
      dimension rData(maxblk*nrData), jData(maxblk*nrData)
	  real eps0, P2, M, N, Pcut, HP, Failure_Strain, S1, S2, S3, Pl_es
	  parameter(eps0 = 0.01)
	  parameter(P2 = 53.2429898)
	  parameter(M = 1.0)
	  parameter(N = 0.75)
	  parameter(Pcut = -2.0)
C  
      jStatus = 1
      call vgetvrm( 'PEEQ', rData, jData, cData, jStatus )
C
	  do k = 1, nblock
        Pl_es = rData(k)
		stateNew(k,7) = Pl_es
      end do
      RETURN
      END