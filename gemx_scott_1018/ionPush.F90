!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Ion pre-push
!
      subroutine ppush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1
      INTEGER :: m,i,j,k,l,n
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta
      REAL(8) :: xt,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap
      REAL(8) :: rhox(4),rhoy(4),psp,pzp

      do m=1,mm(1)
         x=x2(m)
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         z = z2(m)
         k = int(z/dzeq)
         k = min(k,nz-1)
         wz0 = ((k+1)*dzeq-z)/dzeq
         wz1 = 1-wz0


         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) &
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1) 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) &
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1) 
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) &
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1) 
         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1) &
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1) 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) &
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1) 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) &
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1) 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 

         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

!  4 pt. avg. done explicitly for vectorization...
         do 200 l=1,lr(1)
!
            xt=x2(m)+rhox(l) !rwx(1,l)*rhog
            zt=z2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!            zeta=modulo(zeta2(m),pi2)
!     
!   particle can go out of bounds during gyroavg...
            if( (xt<0).or.(xt>lx) ) xt=x2(m)
            if( (zt<0).or.(zt>lz) ) zt=z2(m)
            zeta=zeta2(m)
!            xt=modulo(xs,xdim)
!            zt=modulo(zt,zdim)
            include "ppushli.h"
 200     continue
         exp1 = exp1/4.
         ezp = ezp/4.
         ezetap = ezetap/4.
         delbxp = delbxp/4.
         delbzp = delbzp/4.
!
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp         

         vpar = u2(m)
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         dum1 = 1.
         vxdum = (ezp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp

         zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
             +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

         zetadot =  vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)

         pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
         pzdot = pzd0

         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         x3(m) = x2(m) + 0.5*dt*xdot
         z3(m) = z2(m) + 0.5*dt*zdot
         zeta3(m) = zeta2(m) + 0.5*dt*zetadot
         u3(m) = u2(m) + 0.5*dt*pzdot

         dum = 1.0
         vxdum = (ezp/b+vpar/b*delbxp)*dum1
         vzdum = (-exp1/b+vpar/b*delbzp)*dum1
!         vxdum = eyp+vpar/b*delbxp
         w3(m)=w2(m) + 0.5*dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp
         
      if( (x3(m)>0.).and.(x3(m)<lx).and.(z3(m)>0).and.(z3(m)<lz) ) then
      else
          u3(m)=u2(m)
          x3(m)=x2(m)
          z3(m)=z2(m)
          zeta3(m)=zeta2(m)
          w3(m)=0.
      endif

      enddo

      return
      end
!
!-------------- End of subroutine ppush --------------------------------


      subroutine cpush(n)
!-----------------------------------------------------------------------
!              Ion corrector push
!-----------------------------------------------------------------------
      use gem_com
      use equil
      implicit none
      REAL(8) :: exp1,ezp,ezetap,delbxp,delbzp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,vzdum,dum1
      INTEGER :: m,i,j,k,l,n
      REAL(8) :: rhog,vfac,kapxp,kapzp,vpar,pidum,kaptxp,kapnxp,kaptzp,kapnzp,xnp
      REAL(8) :: b,th,r,enerb,qr,ter,x,z,zeta
      REAL(8) :: xt,xs,zt,xdot,zdot,zetadot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdxp,dbdzp,bfldp,bfldxp,bfldzp,bfldzetap
      REAL(8) :: rhox(4),rhoy(4),psp,pzp

      do m=1,mm(1)
         x=x3(m)
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         z = z3(m)
         k = int(z/dzeq)
         k = min(k,nz-1)
         wz0 = ((k+1)*dzeq-z)/dzeq
         wz1 = 1-wz0


         dbdxp = wx0*wz0*dbdx(i,k)+wx0*wz1*dbdx(i,k+1) &
                 +wx1*wz0*dbdx(i+1,k)+wx1*wz1*dbdx(i+1,k+1) 
         dbdzp = wx0*wz0*dbdz(i,k)+wx0*wz1*dbdz(i,k+1) &
                 +wx1*wz0*dbdz(i+1,k)+wx1*wz1*dbdz(i+1,k+1) 
         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         bfldxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
         bfldzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
         bfldzetap = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
         kaptxp = wx0*wz0*captix(i,k)+wx0*wz1*captix(i,k+1) &
                 +wx1*wz0*captix(i+1,k)+wx1*wz1*captix(i+1,k+1) 
         kapnxp = wx0*wz0*capnix(i,k)+wx0*wz1*capnix(i,k+1) &
                 +wx1*wz0*capnix(i+1,k)+wx1*wz1*capnix(i+1,k+1) 

         kaptzp = wx0*wz0*captiz(i,k)+wx0*wz1*captiz(i,k+1) &
                 +wx1*wz0*captiz(i+1,k)+wx1*wz1*captiz(i+1,k+1) 
         kapnzp = wx0*wz0*capniz(i,k)+wx0*wz1*capniz(i,k+1) &
                 +wx1*wz0*capniz(i+1,k)+wx1*wz1*capniz(i+1,k+1) 

         xnp = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 

         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

         rhox(1) = rhog
         rhoy(1) = 0.
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         exp1=0.
         ezp=0.
         ezetap=0.
         delbxp=0.
         delbzp=0.

!  4 pt. avg. done explicitly for vectorization...
         do 200 l=1,lr(1)
!
!SP            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            xt=x3(m)+rhox(l) !rwx(1,l)*rhog
            zt=z3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            if( (xt<0).or.(xt>lx) ) xt=x3(m)
            if( (zt<0).or.(zt>lz) ) zt=z3(m)
            zeta=zeta3(m)
            include "cpushli.h"
 200     continue
         exp1 = exp1/4.
         ezp = ezp/4.
         ezetap = ezetap/4.
         delbxp = delbxp/4.
         delbzp = delbzp/4.
!
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kapxp = kapnxp - (1.5-vfac/ter)*kaptxp
         kapzp = kapnzp - (1.5-vfac/ter)*kaptzp         

         vpar = u3(m)
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         dum1 = 1.
         vxdum = (ezp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin +vpar*bfldxp/b-enerb/bfldp/bfldp*bfldzetap*dbdzp

         zdot = (-exp1/b+vpar/b*delbzp)*dum1*nonlin &
             +vpar*bfldzp/b+enerb/bfldp/bfldp*bfldzetap*dbdxp

         zetadot =  vpar*bfldzetap/(x*b)+enerb/(x*b*b)*(bfldxp*dbdzp-bfldzp*dbdxp)

         pzd0 = -mu(m)/mims(1)/b*(bfldxp*dbdxp+bfldzp*dbdzp)
         pzdot = pzd0

         edot = q(1)*(xdot*exp1+zdot*ezp+zetadot*ezetap)

         x3(m) = x2(m) + dt*xdot
         z3(m) = z2(m) + dt*zdot
         zeta3(m) = zeta2(m) + dt*zetadot
         u3(m) = u2(m) + dt*pzdot

         dum = 1.0
         vxdum = (ezp/b+vpar/b*delbxp)*dum1
         vzdum = (-exp1/b+vpar/b*delbzp)*dum1
!         vxdum = eyp+vpar/b*delbxp
         w3(m)=w2(m) + dt*(vxdum*kapxp + vzdum*kapzp+edot/ter)*dum*xnp
         
         zeta3(m)=modulo(zeta3(m),pi2)

      if( (x3(m)>0.).and.(x3(m)<lx).and.(z3(m)>0).and.(z3(m)<lz) ) then
          u2(m)=u3(m)
          x2(m)=x3(m)
          z2(m)=z3(m)
          zeta2(m)=zeta3(m)
          w2(m)=w3(m)
      else
          u3(m)=u2(m)
          x3(m)=x2(m)
          z3(m)=z2(m)
          zeta3(m)=zeta2(m)
          w2(m)=0.
          w3(m)=0.
      endif
        

      enddo

      return
      end
