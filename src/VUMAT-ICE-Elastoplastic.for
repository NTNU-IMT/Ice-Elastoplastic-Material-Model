C © 2022, NTNU
C Author: Mojtaba Mokhtari <mojtaba.mokhtari@ntnu.no>      
C This code is licenced under EUROPEAN UNION PUBLIC LICENCE v. 1.2
C
      subroutine vumat(
C Read only -
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     3 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5 stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      parameter( zero = 0., half = 0.5, one = 1.,
     2 two=2.,threehalf=1.5)
C
      real g,twomu, M, N, eps0
      real p1,p2,prse,p,I1,I2,I1e,diff
      real J2,J2e,sj2f,yieldf,dlambda,eps_eh,denom
      
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1 charLength(nblock), strainInc(nblock,ndir+nshr),
     2 relSpinInc(nblock,nshr), tempOld(nblock),
     3 stretchOld(nblock,ndir+nshr),
     4 defgradOld(nblock,ndir+nshr+nshr),
     5 fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6 stateOld(nblock,nstatev), enerInternOld(nblock),
     7 enerInelasOld(nblock), tempNew(nblock),eigval(nblock,ndir),
     8 stretchNew(nblock,ndir+nshr),
     8 defgradNew(nblock,ndir+nshr+nshr),
     9 fieldNew(nblock,nfieldv),
     1 stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2 enerInternNew(nblock), enerInelasNew(nblock)
      dimension strain_dev(6),epsp_dev(6),stress_dev(6),stress_deve(6),
     2            did(6),dide(6),eps_ed(6),eps_ij(6)
      
C    
      real eps_final_in1, Gf1, eps_final_in2, Gf2
      dimension damage_flag(nblock), damage(nblock,6),
     2 displ_fail(nblock), displ_now(nblock),
     3 eps_init(nblock), eps_integ(nblock,6),
     4 eps_plas(nblock,6), eps_inelas(nblock,6),
     5 eps_damag(nblock,6), stres_undamage(6),
     6 stress_inti_dev(6), stres_init(nblock,6),
     7 stress_inti_equv(nblock)
C
C The state variables are stored as
C STATE 1-6   =   plastic strain components
C STATE 7     =   equivelent plastic strain
C STATE 8     =   hydrostatic pressure
C STATE 9     =   eps_mean=eps_mean-stateOld(4)
C STATE 10    =   erosion flag variable 
C STATE 11    =   stateNew(i,11)=eps_mean-eps_eh
C STATE 12    =   Failure strain
C STATE 13    =   Second invariant of the stress tensor, J2
C STATE 14    =   von Mises stress
C STATE 15    =   elastic limit wrt. J2
C STATE 16    =   elastic limit wrt. von Mises
C STATE 17    =   Criteria between plastic loading and elastic unloading 
C STATE 18-23 =   Initial failure stress components
C STATE 25    =   Damage varible in normal directions
C STATE 26    =   Damage varible in shear directions
C STATE 27    =   Maximum damage varible in time and space
      character*80 cmname
C
C---------------user defined properties-----------------------------------------
      e       =   props(1)        ! E-modulus
      xnu     =   props(2)        ! Poisson ratio
      a0      =   props(3)        ! a0 constant in the Tsai-Wu yield surface
      a1      =   props(4)        ! a1 constant in the Tsai-Wu yield surface
      a2      =   props(5)        ! a2 constant in the Tsai-Wu yield surface
      eps0    =   props(6)        ! Initial failure strain
      pcut    =   props(7)        ! cut off pressure in tension
      M       =   props(8)        ! parameter in failure criteria
      N       =   props(9)        ! parameter in failure criteria
      Gf1      =   props(10)      ! parameter determining the final damage strain in normal directions
      Gf2      =   props(11)      ! parameter determining the final damage strain in shear directions
      alphad    =   props(12)     ! parameter determine the speed of the damage accumulation
C      
      twomu   =   e/(one + xnu)           ! Shear modulusx2
      g       =   twomu/2.0               ! Shear modulus/ 2ndLame parameter
      bulk=e/(3.0*(1.0-2.0*xnu))
      alambda =  twomu*xnu/(one-2*xnu)    ! 1st Lame parameter
      ntens	 = ndir+nshr

C
      p1=(-a1+sqrt(a1**2.0-4.0*a0*a2))/(a2*2.0)
      p2=(-a1-sqrt(a1**2.0-4.0*a0*a2))/(a2*2.0)
C
      do i = 1,nblock
          if (totalTime.eq.0.0) then
              trace=sum(strainInc(i,1:ndir))
              stressNew(i,1:ndir)=stressOld(i,1:ndir)+alambda*trace+twomu*strainInc(i,1:ndir)
              stressNew(i,1+ndir:ntens)=stressOld(i,1+ndir:ntens)+twomu*strainInc(i,1+ndir:ntens)
              stateNew(i,18:ntens+17)=stressNew(i,1:ntens)
          else
C     Stress split
              I1=sum(stressOld(i,1:ndir))
			  stateNew(i,28)=I1
              p=-I1/3.0
              stress_dev(1:ndir)=stressOld(i,1:ndir)+p
              stress_dev(1+ndir:ntens)=stressOld(i,1+ndir:ntens)

              J2=sum(stress_dev(1:ndir)**2.)/2.0+
     2        sum(stress_dev(1+ndir:ntens)**2.)
C
C     	strain split
              eps_mean=-sum(strainInc(i,1:ndir))/3.0
              strain_dev(1:ndir)=strainInc(i,1:ndir)+eps_mean
              strain_dev(1+ndir:ntens)=strainInc(i,1+ndir:ntens)
C
C      Correction due to previous plastic behavior
              strain_dev(1:ntens)=strain_dev(1:ntens)-stateOld(i,1:ntens)
              eps_mean=eps_mean-stateOld(i,11)

C      Elastic trial stress state update
              I1e=I1-bulk*(1-stateOld(i,25))*eps_mean*9.0
              prse=-I1e/3.0
C
              stress_deve(1:ndir)=stress_dev(1:ndir)+
     2        twomu*(1-stateOld(i,25))*strain_dev(1:ndir)
              stress_deve(1+ndir:ntens)=stress_dev(1+ndir:ntens)+
     2        twomu*(1-stateOld(i,26))*strain_dev(1+ndir:ntens)
C
              yieldf=a0+a1*prse+a2*prse**two
C Calculation of J2, the seoncd invariant of the stress tensor
              J2e=sum(stress_deve(1:ndir)**2.0)/2.0+
     2        sum(stress_deve(1+ndir:ntens)**2.0)
              stateNew(i,28)=prse
              stateNew(i,29)=sqrt(3*J2e)
C             
              diff=J2e-yieldf
C
              if (diff.ge.0.0) then
                do j=1,200
                  did(1:ndir)=(a1+2.0*a2*prse)/3.0+stress_deve(1:ndir)
                  did(1+ndir:ntens)=stress_deve(1+ndir:ntens)
                  dide(1:ndir)=bulk*(1-stateOld(i,25))*(a1+2.0*a2*prse)+
     2            twomu*(1-stateOld(i,25))*stress_deve(1:ndir)
                  dide(1+ndir:ntens)=twomu*(1-stateOld(i,26))*stress_deve(1+ndir:ntens)
                  denom=sum(did(1:ndir)*dide(1:ndir))+2.0*sum(did(1+ndir:ntens)*dide(1+ndir:ntens))
C
                  !plastic multiplier
                  dlambda=diff/denom
                  epsp_dev(1:ntens)=dlambda*stress_deve(1:ntens)
                  epspp=dlambda*(a1+2.0*a2*prse)
C     correct the stress state to the yield surface
                  stress_deve(1:ndir)=stress_deve(1:ndir)-
     2            twomu*(1-stateOld(i,25))*epsp_dev(1:ndir)
                  stress_deve(1+ndir:ntens)=stress_deve(1+ndir:ntens)-
     2            twomu*(1-stateOld(i,26))*epsp_dev(1+ndir:ntens)
                  !updating the pressure
                  I1e=I1e-3.0*bulk*(1-stateOld(i,25))*epspp
                  prse=-I1e/3.0
C
                  yieldf=a0+a1*prse+a2*prse**two
C
                  J2e=sum(stress_deve(1:ndir)**2.0)/2.0+
     2            sum(stress_deve(1+ndir:ntens)**2.0)
C
                  diff=J2e-yieldf
                  !cFac = abs(diff/J2e)                      !convergence of return algorithm factor
                  accep = 0.0001                             !accepted level. Use as input parameter?
                  if (abs(diff).lt.accep) go to 100             !return algorithm has converged
                enddo
              endif !end plasticity
C
C Stress Update. Adding the hydrostatic pressure to get the full stress state
  100         stressNew(i,1:ndir)=stress_deve(1:ndir)-prse
              stressNew(i,1+ndir:ntens)=stress_deve(1+ndir:ntens)


              stateNew(i,13) = J2e                      !Second invariant of the strss tensor
              stateNew(i,14) = sqrt(3*J2e)              !von Mises stress
              stateNew(i,15) = yieldf                   !Elastic limit wrt J2
              stateNew(i,16) = sqrt(3.0*yieldf)         !Elastic limit wrt. von mises stress
C
              eps_ed(1:ndir)=(stress_deve(1:ndir)-stress_dev(1:ndir))/
     2        (twomu*(1-stateOld(i,25)))
              eps_ij(1:ndir) = strain_dev(1:ndir)-eps_ed(1:ndir)
              stateNew(i,1:ndir)=eps_ij(1:ndir)
              eps_ed(1+ndir:ntens)=(stress_deve(1+ndir:ntens)-stress_dev(1+ndir:ntens))/
     2        (twomu*(1-stateOld(i,26)))
              eps_ij(1+ndir:ntens) = strain_dev(1+ndir:ntens)-eps_ed(1+ndir:ntens)
              stateNew(i,1+ndir:ntens)=eps_ij(1+ndir:ntens)
C
              eps_eh=(prse-p)/(bulk*(1-stateOld(i,25))*3.0)
C
              stateNew(i,11)=eps_mean-eps_eh
              stateNew(i,8)=prse
              stateNew(i,9)=diff
              
              did(1:ndir)=(a1+2.0*a2*prse)/3.0+stress_deve(1:ndir)
              did(1+ndir:ntens)=stress_deve(1+ndir:ntens)
              
C     Criteria between plastic loading and elastic unloading 
              stateNew(i,17)=sum((bulk*eps_mean*9.0/3.0+twomu*strain_dev(1:ndir))*did(1:ndir))+
     2        sum(twomu*strain_dev(1+ndir:ntens)*did(1+ndir:ntens))
C------------------Calculating the equivalent plastic strain incrementally-----------------------------      
              stateNew(i,7)=stateOld(i,7)+sqrt(2.0/3.0*(sum(stateNew(i,1:ndir)**2.)+
     2        + 2.*sum(stateNew(i,1+ndir:ntens)**2.)))

              eps_fail=eps0+(prse/(M*p2)-(N/M))**2.0

C     if damages starts, initial failure strain is kept              
              if (stateOld(i,27).gt.0.0) then
                  stateNew(i,12) = stateOld(i,12)
              else
                  stateNew(i,12) = eps_fail
              endif
              
C---------------------------------Damage stage judgement-----------------------------                
              if (stateNew(i,7).gt.stateNew(i,12)) then
C     Initial failure stress components
                  stateNew(i,18:ntens+17) = stateOld(i,18:ntens+17)
C     Initial failure equivelent stress
				  I2=sum(stateNew(i,18:ndir+17))/3.
                  stress_inti_equv(i) =
     2            sqrt(3./2.*(sum((stateNew(i,18:ndir+17)-I2)**2.)+
     3            2.*sum(stateNew(i,18+ndir:17+ntens)**2.)))

C     Damage evolution rate in two different parts                  
                  eps_final_in1 = 2.0*Gf1/stress_inti_equv(i)
                  eps_final_in2 = 2.0*Gf2/stress_inti_equv(i)
C                  
                  stressNew(i,1:ndir)=stressNew(i,1:ndir)*(1-stateOld(i,25))
                  stressNew(i,1+ndir:ntens)=stressNew(i,1+ndir:ntens)*(1-stateOld(i,26))
C     Damaged equivelent strain                  
                  stateNew(i,24)=stateOld(i,24)+sqrt(2.0/3.0*(
     2            sum((strain_dev(1:ndir)-(stress_deve(1:ndir)-stress_dev(1:ndir))/twomu)**2.)+
     3            2.*sum((strain_dev(1+ndir:ntens)-(stress_deve(1+ndir:ntens)-stress_dev(1+ndir:ntens))/twomu)**2.)))
C     Damage varible for normal and shear                  
                  stateNew(i,25)=(1-EXP(-alphad*(stateNew(i,24)/eps_final_in1)))/
     2            (1-EXP(-alphad))
                  stateNew(i,26)=(1-EXP(-alphad*(stateNew(i,24)/eps_final_in2)))/
     2            (1-EXP(-alphad))
C     Maximum damage varible in time and space
                  stateNew(i,27)=MAX(stateOld(i,27),0.0,stateNew(i,25),
     2            stateNew(i,26))
              else
                  stateNew(i,18:ntens+17)=stressNew(i,1:ntens)
                  stateNew(i,25)=0.0
                  stateNew(i,26)=0.0
              endif
C     Erosion judgement              
              if (stateNew(i,27).gt.0.95) then
                  stateNew(i,10)=0.0
              elseif (prse.le.pcut) then
                  stateNew(i,10)=0.0
              endif
C Update the specific internal energy -
C
            stressPower = 0.5 * (
     2       sum((stressOld(i,1:ndir) + stressNew(i,1:ndir)) * strainInc(i,1:ndir)) +
     3       2.0*sum((stressOld(i,1+ndir:ntens) + stressNew(i,1+ndir:ntens)) * strainInc(i,1+ndir:ntens)))

            enerInternNew(i) = enerInternOld(i) + stressPower/density(i)

C update the dissipated inelastic specific energy
C     
         plasticWorkInc = dlambda*0.5*(
     2      sum((stressOld(i,1:ndir)+stressNew(i,1:ndir))*stress_deve(1:ndir))
     3      +2.0*sum((stressOld(i,1+ndir:ntens)+stressNew(i,1+ndir:ntens))*stressNew(i,1+ndir:ntens)))      
         
         enerInelasNew(i) = enerInelasOld(i) +
     2   plasticWorkInc/density(i)
C       
      endif! timestep check
      enddo ! end nblock
      return
      end
