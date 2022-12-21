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
      real p1,p2,prse,p,I1,I1e,diff
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
      
      twomu   =   e/(one + xnu)           ! Shear modulusx2
      g       =   twomu/2.0               ! Shear modulus/ 2ndLame parameter
      bulk=e/(3.0*(1.0-2.0*xnu))
      alambda =  twomu*xnu/(one-2*xnu)    ! 1st Lame parameter
C
      p1=(-a1+sqrt(a1**2.0-4.0*a0*a2))/(a2*2.0)
      p2=(-a1-sqrt(a1**2.0-4.0*a0*a2))/(a2*2.0)
C
      do i = 1,nblock
          if (totalTime.eq.0.0) then
              trace=strainInc(i,1)+strainInc(i,2)+strainInc(i,3)
              stressNew(i,1)=stressOld(i,1)+alambda*trace+twomu*strainInc(i,1)
              stressNew(i,2)=stressOld(i,2)+alambda*trace+twomu*strainInc(i,2)
              stressNew(i,3)=stressOld(i,3)+alambda*trace+twomu*strainInc(i,3)
              stressNew(i,4)=stressOld(i,4)+twomu*strainInc(i,4)
              stressNew(i,5)=stressOld(i,5)+twomu*strainInc(i,5)
              stressNew(i,6)=stressOld(i,6)+twomu*strainInc(i,6)
              do it=1,6
                  stateNew(i,it+17)=stressNew(i,it)
              enddo
          else
C     Stress split
              I1=(stressOld(i,1)+stressOld(i,2)+stressOld(i,3))
              p=-I1/3.0
              stress_dev(1)=stressOld(i,1)+p
              stress_dev(2)=stressOld(i,2)+p
              stress_dev(3)=stressOld(i,3)+p
              stress_dev(4)=stressOld(i,4)
              stress_dev(5)=stressOld(i,5)
              stress_dev(6)=stressOld(i,6)

              J2=(stress_dev(1)**2.+stress_dev(2)**2.+stress_dev(3)**2.)/2.0+
     2        stress_dev(4)**2.+stress_dev(5)**2.+stress_dev(6)**2.
C
C     	strain split
              eps_mean=-(strainInc(i,1)+strainInc(i,2)+strainInc(i,3))/3.0
              strain_dev(1)=strainInc(i,1)+eps_mean
              strain_dev(2)=strainInc(i,2)+eps_mean
              strain_dev(3)=strainInc(i,3)+eps_mean
              strain_dev(4)=strainInc(i,4)
              strain_dev(5)=strainInc(i,5)
              strain_dev(6)=strainInc(i,6)
C
C      Correction due to previous plastic behavior
              do it=1,6
                  strain_dev(it)=strain_dev(it)-stateOld(i,it)
              enddo
              eps_mean=eps_mean-stateOld(i,11)

C      Elastic trial stress state update
              I1e=I1-bulk*(1-stateOld(i,25))*eps_mean*9.0
              prse=-I1e/3.0
C
              do it =1,3
                  stress_deve(it)=stress_dev(it)+
     2            twomu*(1-stateOld(i,25))*strain_dev(it)
              enddo
              do it =4,6
                  stress_deve(it)=stress_dev(it)+
     2            twomu*(1-stateOld(i,26))*strain_dev(it)
              enddo
C
              yieldf=a0+a1*prse+a2*prse**two
C Calculation of J2, the seoncd invariant of the stress tensor
              J2e=(stress_deve(1)**2.0+stress_deve(2)**2.0+
     2        stress_deve(3)**2.0)/2.0+
     3        stress_deve(4)**2.0+stress_deve(5)**2.0+stress_deve(6)**2.0
              stateNew(i,28)=prse
              stateNew(i,29)=sqrt(3*J2e)
C             
              diff=J2e-yieldf
C
              if (diff.ge.0.0) then
                  do j=1,200

                      did(1)=(a1+2.0*a2*prse)/3.0+stress_deve(1)
                      did(2)=(a1+2.0*a2*prse)/3.0+stress_deve(2)
                      did(3)=(a1+2.0*a2*prse)/3.0+stress_deve(3)
                      did(4)=stress_deve(4)
                      did(5)=stress_deve(5)
                      did(6)=stress_deve(6)
                      dide(1)=bulk*(1-stateOld(i,25))*(a1+2.0*a2*prse)+
     2                twomu*(1-stateOld(i,25))*stress_deve(1)
                      dide(2)=bulk*(1-stateOld(i,25))*(a1+2.0*a2*prse)+
     2                twomu*(1-stateOld(i,25))*stress_deve(2)
                      dide(3)=bulk*(1-stateOld(i,25))*(a1+2.0*a2*prse)+
     2                twomu*(1-stateOld(i,25))*stress_deve(3)
                      dide(4)=twomu*(1-stateOld(i,26))*stress_deve(4)
                      dide(5)=twomu*(1-stateOld(i,26))*stress_deve(5)
                      dide(6)=twomu*(1-stateOld(i,26))*stress_deve(6)
                      denom=did(1)*dide(1)+did(2)*dide(2)+did(3)*dide(3)+
     2                2.0*(did(4)*dide(4)+did(5)*dide(5)+did(6)*dide(6))
C
                      !plastic multiplier
                      dlambda=diff/denom

                      do it=1,6
                          epsp_dev(it)=dlambda*stress_deve(it)
                      enddo
                      epspp=dlambda*(a1+2.0*a2*prse)
C     correct the stress state to the yield surface
                      do it=1,3
                          stress_deve(it)=stress_deve(it)-
     2                    twomu*(1-stateOld(i,25))*epsp_dev(it)
                      enddo
                      do it=4,6
                          stress_deve(it)=stress_deve(it)-
     2                    twomu*(1-stateOld(i,26))*epsp_dev(it)
                      enddo
                  !updating the pressure
                  I1e=I1e-3.0*bulk*(1-stateOld(i,25))*epspp
                  prse=-I1e/3.0

                  yieldf=a0+a1*prse+a2*prse**two

                  J2e=(stress_deve(1)**2.0+stress_deve(2)**2.0
     2            +stress_deve(3)**2.0)/2.0+
     3            stress_deve(4)**2.0+stress_deve(5)**2.0+stress_deve(6)**2.0

                  diff=J2e-yieldf

                  !cFac = abs(diff/J2e)                      !convergence of return algorithm factor
                  accep = 0.0001                             !accepted level. Use as input parameter?
                  if (abs(diff).lt.accep) go to 100             !return algorithm has converged
                  enddo
              endif !end plasticity
C
C Stress Update. Adding the hydrostatic pressure to get the full stress state
  100         stressNew(i,1)=stress_deve(1)-prse
              stressNew(i,2)=stress_deve(2)-prse
              stressNew(i,3)=stress_deve(3)-prse
              stressNew(i,4)=stress_deve(4)
              stressNew(i,5)=stress_deve(5)
              stressNew(i,6)=stress_deve(6)

              stateNew(i,13) = J2e                      !Second invariant of the strss tensor
              stateNew(i,14) = sqrt(3*J2e)              !von Mises stress
              stateNew(i,15) = yieldf                   !Elastic limit wrt J2
              stateNew(i,16) = sqrt(3.0*yieldf)         !Elastic limit wrt. von mises stress
C
              do it=1,3
                  eps_ed(it)=(stress_deve(it)-stress_dev(it))/
     2            (twomu*(1-stateOld(i,25)))
                  eps_ij(it) = strain_dev(it)-eps_ed(it)
                  stateNew(i,it)=eps_ij(it)
              enddo
              do it=4,6
                  eps_ed(it)=(stress_deve(it)-stress_dev(it))/
     2            (twomu*(1-stateOld(i,26)))
                  eps_ij(it) = strain_dev(it)-eps_ed(it)
                  stateNew(i,it)=eps_ij(it)
              enddo

              eps_eh=(prse-p)/(bulk*(1-stateOld(i,25))*3.0)

              stateNew(i,11)=eps_mean-eps_eh
              stateNew(i,8)=prse
              stateNew(i,9)=diff
              
              did(1)=(a1+2.0*a2*prse)/3.0+stress_deve(1)
              did(2)=(a1+2.0*a2*prse)/3.0+stress_deve(2)
              did(3)=(a1+2.0*a2*prse)/3.0+stress_deve(3)
              did(4)=stress_deve(4)
              did(5)=stress_deve(5)
              did(6)=stress_deve(6)
              
C     Criteria between plastic loading and elastic unloading 
              stateNew(i,17)=(bulk*eps_mean*9.0/3.0+twomu*strain_dev(1))*did(1)+
     2        (bulk*eps_mean*9.0/3.0+twomu*strain_dev(2))*did(2)+
     3        (bulk*eps_mean*9.0/3.0+twomu*strain_dev(3))*did(3)+
     4        twomu*strain_dev(4)*did(4)+twomu*strain_dev(5)*did(5)+
     5        twomu*strain_dev(6)*did(6)        
C------------------Calculating the equivalent plastic strain incrementally-----------------------------      
              stateNew(i,7)=stateOld(i,7)+sqrt(2.0/3.0*(
     2        stateNew(i,1)**2.0+stateNew(i,2)**2.0+stateNew(i,3)**2.0+
     3        stateNew(i,4)**2.0+stateNew(i,5)**2.0+stateNew(i,6)**2.0))

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
                  stateNew(i,18) = stateOld(i,18)
                  stateNew(i,19) = stateOld(i,19)
                  stateNew(i,20) = stateOld(i,20)
                  stateNew(i,21) = stateOld(i,21)
                  stateNew(i,22) = stateOld(i,22)
                  stateNew(i,23) = stateOld(i,23)
                  
C     Initial failure equivelent stress
                  stress_inti_equv(i) =
     2            sqrt(((stateNew(i,18)-stateNew(i,19))**2+
     3            (stateNew(i,19)-stateNew(i,20))**2+
     4            (stateNew(i,20)-stateNew(i,18))**2+
     5            6*(stateNew(i,21)**2+stateNew(i,22)**2+stateNew(i,23)**2))/2)

C     Damage evolution rate in two different parts                  
                  eps_final_in1 = 2.0*Gf1/stress_inti_equv(i)
                  eps_final_in2 = 2.0*Gf2/stress_inti_equv(i)
                  
                  stressNew(i,1)=stressNew(i,1)*(1-stateOld(i,25))
                  stressNew(i,2)=stressNew(i,2)*(1-stateOld(i,25))
                  stressNew(i,3)=stressNew(i,3)*(1-stateOld(i,25))
                  stressNew(i,4)=stressNew(i,4)*(1-stateOld(i,26))
                  stressNew(i,5)=stressNew(i,5)*(1-stateOld(i,26))
                  stressNew(i,6)=stressNew(i,6)*(1-stateOld(i,26))

C     Damaged equivelent strain                  
                  stateNew(i,24)=stateOld(i,24)+sqrt(2.0/3.0*(
     2            (strain_dev(1)-(stress_deve(1)-stress_dev(1))/twomu)**2+
     3            (strain_dev(2)-(stress_deve(2)-stress_dev(2))/twomu)**2+
     4            (strain_dev(3)-(stress_deve(3)-stress_dev(3))/twomu)**2+
     5            2.*(strain_dev(4)-(stress_deve(4)-stress_dev(4))/twomu)**2+
     6            2.*(strain_dev(5)-(stress_deve(5)-stress_dev(5))/twomu)**2+
     7            2.*(strain_dev(6)-(stress_deve(6)-stress_dev(6))/twomu)**2))
C     Damage varible for normal and shear                  
                  stateNew(i,25)=(1-EXP(-alphad*(stateNew(i,24)/eps_final_in1)))/
     2            (1-EXP(-alphad))
                  stateNew(i,26)=(1-EXP(-alphad*(stateNew(i,24)/eps_final_in2)))/
     2            (1-EXP(-alphad))
C     Maximum damage varible in time and space
                  stateNew(i,27)=MAX(stateOld(i,27),0.0,stateNew(i,25),
     2            stateNew(i,26))
              else
                  do it=1,6
                  stateNew(i,it+17)=stressNew(i,it)
                  enddo
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
     2       (stressOld(i,1) + stressNew(i,1)) * strainInc(i,1) +
     3       (stressOld(i,2) + stressNew(i,2)) * strainInc(i,2) +
     4       (stressOld(i,3) + stressNew(i,3)) * strainInc(i,3)+
     5       2.0*(stressOld(i,4) + stressNew(i,4)) * strainInc(i,4) +
     6       2.0*(stressOld(i,5) + stressNew(i,5)) * strainInc(i,5) +
     7       2.0*(stressOld(i,6) + stressNew(i,6)) * strainInc(i,6))
            enerInternNew(i) = enerInternOld(i) + stressPower/density(i)

C update the dissipated inelastic specific energy
C     
         plasticWorkInc = dlambda*0.5*(
     2      (stressOld(i,1)+stressNew(i,1))*stress_deve(1)
     3      +(stressOld(i,2)+stressNew(i,2))*stress_deve(2)
     4      +(stressOld(i,3)+stressNew(i,3))*stress_deve(3)
     5      +2.0*(stressOld(i,4)+stressNew(i,4))*stressNew(i,4)      
     6      +2.0*(stressOld(i,5)+stressNew(i,5))*stressNew(i,5)    
     7      +2.0*(stressOld(i,6)+stressNew(i,6))*stressNew(i,6))
         
         enerInelasNew(i) = enerInelasOld(i) +
     2   plasticWorkInc/density(i)
C       
      endif! timestep check
      enddo ! end nblock
      return
      end
