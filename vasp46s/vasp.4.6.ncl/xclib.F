#include "symbol.inc"
! RCS:  $Id: xclib.F,v 1.2 2000/11/15 08:23:51 kresse Exp $

      FUNCTION ECCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation energy as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),B(2),C(2),D(2),G(2),B1(2),B2(2)
      SAVE A,B,C,D,B1,B2,G
      DATA A/0.0622_q,0.0311_q/,B/-0.0960_q,-0.0538_q/,C/0.0040_q,0.0014_q/, &
     &     D/-0.0232_q,-0.0096_q/,B1/1.0529_q,1.3981_q/,B2/0.3334_q,0.2611_q/, &
     &     G/-0.2846_q,-0.1686_q/
!KRESSE/FURTH---get a continuosly differentiable energy functional
      c(1) = 0.004038664055501747_q
      d(1) =-0.023264632546756681_q
      b2(1)= 0.333390000000000000_q
      c(2) = 0.001395274602717559_q
      d(2) =-0.009602765503781227_q
      b2(2)= 0.261090000000000000_q
!KRESSE/FURTH
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         ECCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         ECCA=G(IFLG)/(1.0_q+B1(IFLG)*RSQ+B2(IFLG)*RS)
      END IF
      RETURN
      END

      FUNCTION VCCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation potential as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),B(2),C(2),D(2),BT1(2),BT2(2)
      PARAMETER(X76=7.0_q/6.0_q, X43=4.0_q/3.0_q, &
     &  AP=0.03110_q*2.0_q, BP=-0.0480_q*2.0_q, CP=0.0020_q*2.0_q, DP=-0.0116_q*2.0_q, &
     &  AF=0.01555_q*2.0_q, BF=-0.0269_q*2.0_q, CF=0.0007_q*2.0_q, DF=-0.0048_q*2.0_q, &
     &  BP1=BP-AP/3.0_q, CP1=2.0_q*CP/3.0_q, DP1=(2.0_q*DP-CP)/3.0_q, &
     &  BF1=BF-AF/3.0_q, CF1=2.0_q*CF/3.0_q, DF1=(2.0_q*DF-CF)/3.0_q)
      SAVE A,B,C,D,BT1,BT2
      DATA A/AP,AF/,B/BP1,BF1/,C/CP1,CF1/,D/DP1,DF1/, &
     &     BT1/1.0529_q,1.3981_q/,BT2/0.3334_q,0.2611_q/
!KRESSE/FURTH---get a continous energy functional
      c(1)  = 0.004038664055501747_q * 2._q/3._q
      d(1)  =-0.023264632546756681_q * 2._q/3._q  -  0.004038664055501747_q / 3._q
      bt2(1)= 0.333390000000000000_q
      c(2)  = 0.001395274602717559_q * 2._q/3._q
      d(2)  =-0.009602765503781227_q * 2._q/3._q  -  0.001395274602717559_q / 3._q
      bt2(2)= 0.261090000000000000_q
!KRESSE/FURTH
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         VCCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         VCCA=ECCA(RS,IFLG)*(1.0_q+X76*BT1(IFLG)*RSQ+X43*BT2(IFLG)*RS) &
     &                     /(1.0_q+    BT1(IFLG)*RSQ+    BT2(IFLG)*RS)
      END IF
      RETURN
      END


      

      FUNCTION ECVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The Ceperley-Alder exchange energy as given by the Pade approximation
! technique of Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with
! the parameters given in table 5, page 1207).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      ECVO=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=-1._q*ECVO*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=ECVO*A(IFLG)
      RETURN
      END

      FUNCTION VCVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ECVO(RS,IFLG)-RS/3.*d(ECVO(RS,IFLG))/d(RS) with function
! ECVO(RS,IFLG) given above (Ceperley-Alder potential derived from the
! approximation for ECVO of Vosko et al. discussed above).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      VCVO=-4._q*B*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      VCVO=VCVO-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      VCVO=ECVO(RS,IFLG)-VCVO*A(IFLG)*X/6._q
      RETURN
      END

      FUNCTION ECGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      X=RS/R(IFLG)
      ECGL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
      END

      FUNCTION VCGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      VCGL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
      END

      FUNCTION ECHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation energy (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      X=RS/R(IFLG)
      ECHL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
      END

      FUNCTION VCHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation potential (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      VCHL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
      END

      FUNCTION ECBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      X=RS/R(IFLG)
      ECBH=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
      END

      FUNCTION VCBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      VCBH=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
      END

      FUNCTION ECWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation energy (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X=CX(IFLG)/R(IFLG)*(1._q+C(IFLG)/RS)
      ECWI=-CX(IFLG)/X/RS
      RETURN
      END

      FUNCTION VCWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation potential (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X1=C(IFLG)/RS
      X2=1._q+X1
      X3=CX(IFLG)/R(IFLG)*X2
      B=1._q+1._q/X3
      F=1._q-X1/X2/(1._q+X3)
      F=1._q+F/3._q
      E=-CX(IFLG)*B/RS
      VCWI=(4._q/3._q-F*B)*CX(IFLG)/RS
      RETURN
      END

      FUNCTION EX(RS,IFLG,TREL)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/0.9163305865663_q,1.1545041946774_q/,CBETA/0.0140_q/
      EX=-CX(IFLG)/RS
      IF (TREL) THEN
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/(B*B)
         F=(SQRT(1+B*B)/B)-F
!jF: the expression given above becomes numerically extremely instable for
!    very small values of B (small difference of two large numbers divided
!    by small number = noise) therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=B*(2._q/3._q-0.7_q*B*B)
         EX=(1._q-1.5_q*F*F)*EX
      END IF
      RETURN
      END

      FUNCTION VX(RS,IFLG,TREL)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/1.2217741154217_q,1.5393389262365_q/,CBETA/0.0140_q/
      VX=-CX(IFLG)/RS
      IF (TREL) THEN
! Warning error in the paper of Bachelet et al. !!
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/B/SQRT(1+(B*B))
!        F=LOG(B+(SQRT(1+B*B)))/B/(1+(B*B))
!jF: potentially the expression given above becomes numerically instable for
!    very small values of B, therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=1._q-B*B*(2._q/3._q-B*B*31._q/30._q)
         VX=(-0.5_q+1.5_q*F)*VX
      END IF
      RETURN
      END

      FUNCTION FZ0(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Interpolation function between paramagnetic and ferromagnetic results
! for exchange energy and exchange potential. The parameter ZETA is:
! ZETA=(RHO[upspin] - RHO[downspin]) / (RHO[upspin] + RHO[downspin]).
      PARAMETER(C43=4._q/3._q,FAC=1.92366105093153632_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ0=0._q
      ELSE IF (Z>=1._q) THEN
         FZ0=1._q
      ELSE
         FZ0=(((1._q+Z)**C43)+((1._q-Z)**C43)-2._q)*FAC
      END IF
      RETURN
      END

      FUNCTION FZ1(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The derivative dFZ0(ZETA)/d(ZETA), FZ0(ZETA) given above.
      PARAMETER(C13=1._q/3._q,FAC=2.56488140124204843_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ1=0._q
      ELSE IF (Z>=1._q) THEN
         FZ1=(2._q**C13)*FAC
      ELSE
         FZ1=(((1._q+Z)**C13)-((1._q-Z))**C13)*FAC
      END IF
      RETURN
      END

      FUNCTION ALPHA0(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The spin stiffness [d(EXC(RS,ZETA))**2/d**2(ZETA) |ZETA=0] as given by
! Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with the parameters
! given on page 1209 [fitting to low-density values using eq.{4.7.}]).
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.019751631321681_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA0=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=-1._q*ALPHA0*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=ALPHA0*A
      RETURN
      END

      FUNCTION ALPHA1(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ALPHA0(RS)-RS/3.*dALPHA0(RS)/dRS, ALPHA0(RS) given above.
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.0032919385536135_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA1=-4._q*B/QQ*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      ALPHA1=ALPHA1-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      ALPHA1=ALPHA0(RS)-ALPHA1*A*X
      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE CORPBE_LDA(RS,ZET,EC,VCUP,VCDN)
!----------------------------------------------------------------------
!  LDA part of the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical*4 lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor_xc(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
      END

!----------------------------------------------------------------------
      FUNCTION PBE_ALPHA(RS)
!----------------------------------------------------------------------
!  Alpha (spin stiffness)
!  from the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!  spin stiffness 
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))

      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units

      PBE_ALPHA = -ALFM/FZZ*2
      RETURN
      END

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE GCOR_XC(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
