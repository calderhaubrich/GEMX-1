      program gemx
#include <petsc/finclude/petscksp.h>
      use gem_com
      use equil

      use petsc
      use petscdmda
      use petscksp
!      use para_com

!      use gem_com
!      use equil
      implicit none

       integer :: status
       integer :: n,i,j,k,ip,m,idx
       real::random
       PetscInt is,js,iw,jw
       PetscInt one,three
       PetscErrorCode petsc_ierr
       PetscScalar, POINTER ::phi_array(:)
       KSP ksp
       DM dm
       Vec petsc_phi
       PetscObject  vec
       PetscViewer viewer
       external ComputeRHS,ComputeMatrix,ComputeInitialGuess
call initialize

       open(935, file='flag_debug',status='unknown',position='append')
       write(935,*)'0001'
       close(935)

  
       one = 1
       three = 3

       PetscCallA(PetscInitialize(petsc_ierr))
!         open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'0002'
!  close(935)
       PetscCallA(KSPCreate(PETSC_COMM_WORLD,ksp,petsc_ierr))
       PetscCallA(DMDACreate2D(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,nx+1,nz+1,PETSC_DECIDE,PETSC_DECIDE,one,one, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, dm, petsc_ierr))
       PetscCallA(DMSetFromOptions(dm,petsc_ierr))
       PetscCallA(DMSetUp(dm,petsc_ierr))
       PetscCallA(KSPSetDM(ksp,dm,petsc_ierr))
 ! open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'0003'
 ! close(935)
       PetscCallA(KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,petsc_ierr))
 !                        open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'0004'
 ! close(935)

       PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,0,petsc_ierr))
       PetscCallA(KSPSetComputeOperators(ksp,ComputeMatrix,0,petsc_ierr))
!                           open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'0005'
 ! close(935)
        	
       PetscCallA(DMDAGetCorners(dm,is,js,PETSC_NULL_INTEGER,iw,jw,PETSC_NULL_INTEGER,petsc_ierr))
!                           open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'0006'
!  close(935)
       PetscCallA(KSPSetFromOptions(ksp,petsc_ierr))
!                           open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'0007'
!  close(935)
       PetscCallA(KSPSetUp(ksp,petsc_ierr))
!                          open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'0008'
!  close(935)

       



         

!      implicit none


!	call particle loader and initialization
       if(iget.eq.0)call loadi

       
        starttm=MPI_WTIME()


 !       call generate_LHS_gkps()

 ! open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'initialize OK'
 ! close(935)
 if (ifield_solver .eq. 1) then         
    ncurr = 1                      
    nm = 1                 
 end if

        do  timestep=ncurr,nm
           tcurr = tcurr+dt

!	   call accumulate(timestep-1,0)
!	   call ezamp
!	   call gkps
!     call field(timestep-1,0)



 !           open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'before solve RHS, time=', timestep
           ! close(935)
           if(ifield_solver .eq. 1) then
           do k=0,1!kmx


              !              if(tcurr==1) then
              do i=0,imx
                 do j=0,jmx
                    call random_number(random)
                    dene(i,j,k)=mask(i,j)*2*(random-0.5)
                 enddo
              enddo
              !              endif
              
              open(unit=11, file = './out/testne',status='unknown',action='write')
              do j=0,jmx

                 write(11,*) dene(:,j,k)
                 enddo
              close(11)
              
                    
   PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))

!              open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'before solve setup, time=', timestep
!  close(935)
!  PetscCallA(KSPSetUp(ksp,petsc_ierr))
!              open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'before solve gk, time=', timestep
!  close(935)
               PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
!                open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'before solve getsolution, time=', timestep
!  close(935)
               PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))
!                open(935, file='flag_debug',status='unknown',position='append')
!  write(935,*)'after getsolution, time=', timestep
!  close(935)

! if (k==0) then
!     PetscCall(VecView(petsc_phi, PETSC_VIEWER_STDOUT_WORLD,ierr))
!endif
               PetscCall(VecGetArrayF90(petsc_phi, phi_array, petsc_ierr))

!       do i=0,imx
    !          do j=0,jmx
               do idx=0,(imx+1)*(jmx+1)-1
                  i=mod(idx,(imx+1))
                  j=idx/(imx+1)
                  phi(i,j,k)=phi_array(idx)!*mask(i,j)
!         enddo
             !      enddo
               enddo
   
               if (k==1) then
                  open(unit=11, file = './out/testphi',status='unknown',action='write')
                  do j=0,jmx

                     write(11,*) phi(:,j,k)
                  enddo
                  close(11)
               endif       

       
        enddo
         

           call get_apar(-1)
           call get_jpar(apars)
           call get_ne(-1)
           

           
  open(935, file='flag_debug',status='unknown',position='append')
  write(935,*)'before push, time=', timestep
  close(935)

  write(*,*)'before push, time=', timestep
  end if
  if (ifield_solver == 1) then
  else
!       if(ision==1)call ppush(timestep)
       if(ifluid==1)call pintef
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif


 ! open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'timestep=', timestep
 ! close(935)

!	   call accumulate(timestep,1)
!	   call ezamp
!	   call gkps
!	   call field(timestep,1)

  if (ifield_solver .eq. 1) then
do k=0,1!kmx
      PetscCallA(KSPSetComputeRHS(ksp,ComputeRHS,k,petsc_ierr))
!      PetscCallA(KSPSetUp(ksp,petsc_ierr))
       PetscCallA(KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,petsc_ierr))
       PetscCallA(KSPGetSolution(ksp,petsc_phi,petsc_ierr))

       !     call VecGetArray(petsc_phi, phi_array, ierr)

        PetscCall(VecGetArrayF90(petsc_phi, phi_array, petsc_ierr))

    do idx=0,(imx+1)*(jmx+1)-1
       i=mod(idx,(imx+1))
       j=idx/(imx+1)
             phi(i,j,k)=phi_array(idx)!*mask(i,j)

          enddo
!      if (k==0) then
 !              open(unit=11, file = 'testphi',status='unknown',action='write')
  !             write(11,*) phi(:,:,0)
   !            close(11)
    !        endif
            
         
    enddo
    

      call get_apar(0)
      call get_jpar(apar)
      call get_ne(0)

   endif
   
           
!      open(unit=11, file = 'flag.dat',status='unknown',action='write')
!      write(11,*) 'OK here'
!      close(11)
    if (ifield_solver == 1) then
    else
        if(ision==1)call cpush(timestep)
        if(ifluid==1)call cintef(timestep)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif

 ! open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'timestep=', timestep
 ! close(935)
!sp    if(myid.eq.master .and. mod(timestep,xnplt)==0)then
!sp       open(9,file='plot',status='unknown',position='append')
!sp       m = 2
!sp       i = 4
!sp       write(9,10)timestep,(x2(m+i*j),z2(m+i*j),j=1,7)
!sp 10    format(1x,i6,16(1x,e10.3))
!sp       close(9)
!sp    endif
          


 !  open(935, file='flag_debug',status='unknown',position='append')
 ! write(935,*)'finial of particle push in a step.'
 ! close(935)

           call outd(timestep)
         end do


!               open(unit=11, file = 'debug.dat',status='unknown',action='write')
!      write(11,*) dR
!      close(11)
	 lasttm=MPI_WTIME()
  tottm=lasttm-starttm


       PetscCallA(KSPDestroy(ksp,petsc_ierr))
       PetscCallA(DMDestroy(dm,petsc_ierr))
       PetscCallA(PetscFinalize(petsc_ierr))

 100     call MPI_FINALIZE(ierr)
         end program gemx
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init
      
      use gem_com
      use equil
      implicit none
      character*(62) dumchar
      INTEGER :: i,j,k,n,ns,idum,i1,k1,m,j1
      INTEGER :: lr1                                      
      REAL(8) :: x,z,dum,zdum
      REAL(8) :: dbdrp,dbdtp,bfldp,upae0p,dnuobdrp,dnuobdtp,btorp,bxp,bzp
      REAL(8) :: gn0ip,gn0ep,gt0ip,gt0ep,capnxp,capnzp
      REAL(8) :: wx0,wx1,wz0,wz1,b

      IU=cmplx(0.,1.)
      pi=4.0*atan(1.0)
      pi2 = pi*2.

      open(115,file='gem.in')
      read(115,*) dumchar
      read(115,*) imx,jmx,kmx,mmx,nmx,nsmx,ntube
      read(115,*) dumchar
      read(115,*) dt,nm,nsm,iez
      read(115,*) dumchar
      read(115,*) iput,iget,ision,peritr
      read(115,*) dumchar
      read(115,*) nplot,xnplt
      read(115,*) dumchar
      read(115,*) cut,amp,tor
      read(115,*) dumchar
      read(115,*) etaohm
      read(115,*) dumchar
      read(115,*) ifluid,amie,rneu
      read(115,*) dumchar
      read(115,*) beta,nonlin,nonline,vcut
      read(115,*) dumchar
      read(115,*) ntracer,ifield_solver                              
      close(115)
      if(myid.eq.master)then
         open(9,file='plot',status='unknown',position='append')
         write(9,*)'dt,beta= ',dt, beta
         write(9,*)'imx,jmx,kmx,ntracer= ',imx,jmx,kmx,ntracer
         close(9)
      end if

      nsm=1
      
      call new_gem_com()
      ns = 1
      tmm(ns)=ntracer                               
      mm(ns)=int(ntracer/numprocs)                   

!     write(*,*)'in init  ',Myid,mm(ns)
      mims(ns)=1.0
      q(ns)=1.0
      lr(ns)=4

      emass = 1./amie
      qel = -1

      call new_equil()
      lx = xdim
      lz = zdim
      
      if(myid.eq.master)then
         open(9,file='plot',status='unknown',position='append')
         write(9,*)'a,rmaj0,lx,lz= ',a,rmaj0,lx,lz
         write(9,*)'xctr,xdim=',xctr,xdim
         close(9)
      end if

      iadi = 0

      if(iget.eq.1) amp=0.

      dx=lx/float(imx)
      dz=lz/float(jmx)
      dzeta=pi2/kmx
!     
      do 10 i=0,imx
         xg(i)=i*dx 
 10   continue
      do 14 k=0,jmx
         zg(k)=k*dz
 14   continue

      do i1 = 0,imx
         x = i1*dx+xctr-xdim/2
         i = int(x/dxeq)
         i = min(i,nx-1)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         do k1 = 0,jmx
            z = k1*dz
            k = int(z/dzeq)
            k = min(k,nz-1)            
            wz0 = ((k+1)*dzeq-z)/dzeq
            wz1 = 1-wz0

            bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
            btorp = wx0*wz0*b0zeta(i,k)+wx0*wz1*b0zeta(i,k+1) &
                 +wx1*wz0*b0zeta(i+1,k)+wx1*wz1*b0zeta(i+1,k+1) 
            bxp = wx0*wz0*b0x(i,k)+wx0*wz1*b0x(i,k+1) &
                 +wx1*wz0*b0x(i+1,k)+wx1*wz1*b0x(i+1,k+1) 
            bzp = wx0*wz0*b0z(i,k)+wx0*wz1*b0z(i,k+1) &
                 +wx1*wz0*b0z(i+1,k)+wx1*wz1*b0z(i+1,k+1) 
            gt0ip = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 
            gt0ep = wx0*wz0*t0e(i,k)+wx0*wz1*t0e(i,k+1) &
                 +wx1*wz0*t0e(i+1,k)+wx1*wz1*t0e(i+1,k+1) 
            gn0ip = wx0*wz0*xn0i(i,k)+wx0*wz1*xn0i(i,k+1) &
                 +wx1*wz0*xn0i(i+1,k)+wx1*wz1*xn0i(i+1,k+1) 
            gn0ep = wx0*wz0*xn0e(i,k)+wx0*wz1*xn0e(i,k+1) &
                 +wx1*wz0*xn0e(i+1,k)+wx1*wz1*xn0e(i+1,k+1) 
            capnxp = wx0*wz0*capnex(i,k)+wx0*wz1*capnex(i,k+1) &
                 +wx1*wz0*capnex(i+1,k)+wx1*wz1*capnex(i+1,k+1) 
            capnzp = wx0*wz0*capnez(i,k)+wx0*wz1*capnez(i,k+1) &
                 +wx1*wz0*capnez(i+1,k)+wx1*wz1*capnez(i+1,k+1) 

            b=1.-tor+tor*bfldp
            bmag(i1,k1) = b
            gbtor(i1,k1) = btorp
            gbx(i1,k1) = bxp
            gbz(i1,k1) = bzp                        

            gt0i(i1,k1) = gt0ip
            gt0e(i1,k1) = gt0ep
            gn0e(i1,k1) = gn0ep
            gn0i(i1,k1) = gn0ip
            gcpnex(i1,k1) = capnxp
            gcpnez(i1,k1) =  capnzp           

            gupae0(i1,k1) = upae0p
!            gnuoby(i1,k1) = (-dydrp*dnuobdtp+r0/q0*qhatp*dnuobdrp)*fp/radiusp*grcgtp
!            gnuobx(i1,k1) = dnuobdtp*fp/radiusp*grcgtp
         end do
      end do


      iseed = -(1777+myid*13)
      idum = ran2(iseed)
      phi = 0.
      apar = 0.
      dene = 0.
      upar = 0.


      do i = 0,imx
         do j = 0,jmx
            do k = 0,kmx
               phi(i,j,k) = amp*(ran2(idum)-0.5)*ifluid*1e-8  !amp*sin(nzcrt*xg(i)*pi/lx) !
               dene(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-8
               apar(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-10 
            end do
         end do
      end do

      if(myid.eq.master)then
         open(9,file='plot',status='unknown',position='append')
         write(9,*)'inner,outer = ',xctr-xdim/2,xctr+xdim/2
         write(9,*)'mi,qi=',mims(1),q(1)
         write(9,*)'mm(1)=',mm(1)
         close(9)
      end if

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grad(ip)
  
!  currently set up for periodic in x,y,z

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,ip
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: tmp(0:imx,0:jmx,0:1),uoverb(0:imx,0:jmx,0:1)
      real(8) :: v(0:imx-1),dum,dum1

      call gradu(phi(:,:,:),ux,uy)
      ex(:,:,:) = -ux(:,:,:)
      ez(:,:,:) = -uy(:,:,:)

      delbx = 0.
      delby = 0.
      if(ifluid.eq.1)then
         call gradu(apar(:,:,:),ux,uy)
         delbx(:,:,:) = uy(:,:,:)
         delby(:,:,:) = -ux(:,:,:)
      end if

      call gradu(tmp(:,:,:),ux,uy)
      dnedx(:,:,:) = ux(:,:,:)
      dnedy(:,:,:) = uy(:,:,:)
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               uoverb(i,j,k) = upar(i,j,k)  !/bfld(i,k)
            end do
         end do
      end do
      call gradu(uoverb(:,:,:),ux,uy)
      dupadx(:,:,:) = ux(:,:,:)
      dupady(:,:,:) = uy(:,:,:)

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid1(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: x,z
      INTEGER :: m,n,i,j,k,l,ns,ip
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,ter
      REAL(8) :: wght,r,b,bfldp,dv
      REAL(8) :: xt,zt,rhog,pidum,vpar,xs,vfac
      real(8) :: myden(0:imx,0:jmx,0:kmx),myjpar(0:imx,0:jmx,0:kmx)
      REAL(8) :: rhox(4),rhoy(4)

      ns=1
      rho=0.
      den=0.
      jpar = 0.
      myden = 0.
      myjpar = 0.

      do m=1,mm(1)
         dv=float(lr(1))*(dx*dz*dzeta)

         x=x3(m)
         i = int(x/dxeq)
         wx0 = ((i+1)*dxeq-x)/dxeq
         wx1 = 1.-wx0

         z = z3(m)
         k = int(z/dzeq)
         wz0 = ((k+1)*dzeq-z)/dzeq
         wz1 = 1-wz0


         bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                 +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
         ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                 +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 

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

         vfac=0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b )
         wght=w3(m)/dv
         if(vfac/ter > vcut)wght=0.
         vpar = u3(m)

!    now do 1,2,4 point average, where lr is the no. of points...
         do 100 l=1,lr(1)
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            zt=z3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=modulo(xs,xdim)
            zt=modulo(zt,zdim)

            include "gridli.h"
 100     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity

      do 110 i=0,imx
         do 120 j=0,jmx
            do 130 k=0,kmx
               den(1,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i)
               jpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i)*ifluid
 130        continue
 120     continue
 110  continue


      do 150 i=0,imx
         do 160 j=0,jmx
            do 170 k=0,kmx
               rho(i,j,k)=rho(i,j,k)+den(1,i,j,k)
 170        continue
 160     continue
 150  continue

 499  continue
      do i = 0,imx
         do j = 0,jmx
            do k = 0,kmx
               rho(i,j,k) = ision*rho(i,j,k) + dene(i,j,k)*qel/ntube
            enddo
         enddo
      enddo      

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!        Normal distribution random no. generator, stand. dev. = 1.
!        Version 2 does it Willy's way...

         subroutine parperp(vpar,vperp2,m,pi,cnt,MyId)

         REAL(8) :: vpar,vperp2,r1,r2,t,pi
         INTEGER :: m,iflag,cnt,MyId
         REAL(8) :: c0,c1,c2
         REAL(8) :: d1,d2,d3
         data c0,c1,c2/2.515517,0.802853,0.010328/
         data d1,d2,d3/1.432788,0.189269,0.001308/


          r1=revers(m+MyId*cnt,7)
          r2=revers(m+MyId*cnt,11)


!.....quiet start---see denavit pf '71(?) & abramowitz hand book
!.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
! warning: we have g1=1 in the x-direction. This surpresses all odd
!          modes in the x-direction!!!

         iflag=1
         if(r1.le.0.5) go to 110
         r1=1.-r1
         iflag=-1
  110    continue
         if(r1.ge.1.e-6) then
           t=sqrt(log(1.0/(r1*r1)))
         else
           t=5.0
           write(*,*)'parperp2 warning  m= ',m
         endif

         vpar=t-(c0+c1*t+c2*t**2)/(1.+d1*t+d2*t**2+d3*t**3)
         vpar=vpar*iflag

          vperp2=-2.0*dlog(r2)

        return
        end

!---------------------------------------------------------------------- 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gkps()
      use gem_com
      use equil
      use petsc
      use petscdmda
      use petscksp
      

      return
      end subroutine gkps
      !real,dimension(nx,nz,nzeta)::nepredict,aparpredic
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_apar(flagnumner)
         use gem_com
         use equil


         integer flagnumner, i, j, k
         if (flagnumner == -1) then
            apars=apar-0.5*dt*(gradpar(phi)) !+Epar)
            else
               apar = apar -dt*(gradpar(phi)) !+Epar)    
            endif
         !   do k=0,kmx-1
         !      apar=apra(0:imx-1,0:jmx-1,k)*mask
         !   end do
            





      CONTAINS

      function  gradpar(matrix)
   !   use gem_com
   !   use equil
      real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradpar
 !     real,dimension(nx,nz,nzeta+4)::ghostmatrix
      

 !     ghostmetric(:,:,0:1)=matrix(:,:,nzeta-2:nzeta-1)
 !     ghostmetric(:,:,2:nzeta+1)=matrix(:,:,0:nzeta-1)
 !     ghostmetric(:,:,nzeta+2:nzeta+3)=matrix(:,:,0:1)
      do k=1,(kmx-2)
         do i=2,(imx-3)
            do j=2,(jmx-3)
               gradpar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
               
            enddo
         enddo
      enddo

         do i=2,(imx-3)
            do j=2,(jmx-3)
               gradpar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,nzeta))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
            enddo
         enddo


         do i=2,(imx-3)
            do j=2,(jmx-3)
               gradpar(i,j,kmx-1)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,kmx-1)-matrix(i-1,j,kmx-1))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,kmx-1)-matrix(i,j-1,kmx-1))*0.5/dz                        &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,kmx-2))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
            enddo
         enddo
      
 !     return
      end function gradpar
      end subroutine get_apar


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine get_jpar(matrix)
         use gem_com
         use equil
         real, dimension(0:imx-1,0:jmx-1,0:kmx-1):: matrix
         integer i,j,k
         do k=0,kmx-1
            do i=2,imx-3
               do j=2,jmx-3
                  jpar(i,j,k)=(-(matrix(i+1,j,k)+matrix(i-1,j,k)-2*matrix(i,j,k))/dx**2       &
                                 -(matrix(i,j+1,k)+matrix(i,j-1,k)-2*matrix(i,j,k))/dz**2   &
                                 -(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/(dx*Rgrid(i)))*mask(i,j)
               enddo
            enddo
         enddo

         !do k=1,kmx-1
         !   jpar=jpar(:,:,k)*mask
        ! end do
         
         end subroutine get_jpar

 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine get_ne(flagnumner)
            use gem_com
            use equil



            integer flagnumner, i, j, k
            if (flagnumner == -1) then
               denes=dene-0.5*dt*(gradpar(jpar)) 
               else
                  dene = dene -dt*(gradpar(jpar))   
               endif

           !    do k=0,kmx-1
           !       dene=dene(:,:,k)*mask
           !    end do
               
   
   
   
   
   
         CONTAINS
   
         function  gradpar(matrix)
      !   use gem_com
      !   use equil
         real,dimension(0:imx,0:jmx,0:kmx)::matrix,gradpar
    !     real,dimension(nx,nz,nzeta+4)::ghostmatrix
    !     integer i,j,k
   
    !     ghostmetric(:,:,0:1)=matrix(:,:,nzeta-2:nzeta-1)
    !     ghostmetric(:,:,2:nzeta+1)=matrix(:,:,0:nzeta-1)
    !     ghostmetric(:,:,nzeta+2:nzeta+3)=matrix(:,:,0:1)
         do k=1,(kmx-2)
            do i=2,(imx-3)
               do j=2,(jmx-3)
                  gradpar(i,j,k)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  &
                  +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                &
                  +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
               enddo
            enddo
         enddo
   
            do i=2,(imx-3)
               do j=2,(jmx-3)
                  gradpar(i,j,0)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  &
                  +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                &
                  +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,nzeta))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
               enddo
            enddo
   


         do i=2,(imx-3)
            do j=2,(jmx-3)
               gradpar(i,j,kmx-1)=(b0x(i,j)/b0(i,j)*(matrix(i+1,j,kmx-1)-matrix(i-1,j,kmx-1))*0.5/dx  &
               +b0z(i,j)/b0(i,j)*(matrix(i,j+1,kmx-1)-matrix(i,j-1,kmx-1))*0.5/dz                        &
               +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,kmx-2))*kmx/(Rgrid(i)*4*pi))*mask(i,j)
            enddo
         enddo


           
    !     return
         end function gradpar
         end subroutine get_ne
!      jpar=-gradper2(apars(:,:,:))
!      denes=dene+0.5*dt*gradpar(jpar(:,:,:))
!      phi=gkpoisson(denes(:,:,:))




!      apar=apar-dt*(gradpar(phi(:,:,:)))!+Epar)
 !     jpar=-gradper2(apar(:,:,:)) 
 !     dene=dene+dt*gradpar(jpar(:,:,:))
  !    phi=gkpoisson(dene(:,:,:))


!       open(unit=11, file = 'flag.dat',status='unknown',action='write')
!       write(11,*) 'OK here!, dx=',dx
!       close(11)
      
! !      implicit none
!      return
      
!CONTAINS

!       function  gradpar(matrix)
!       use gem_com
!       use equil
! !      real,dimension(0:nx,0:nz,0:nzeta)::matrix,gradpar
!  !     real,dimension(nx,nz,nzeta+4)::ghostmatrix
!       integer i,j,k

!  !     ghostmetric(:,:,0:1)=matrix(:,:,nzeta-2:nzeta-1)
!  !     ghostmetric(:,:,2:nzeta+1)=matrix(:,:,0:nzeta-1)
!  !     ghostmetric(:,:,nzeta+2:nzeta+3)=matrix(:,:,0:1)
!       do k=1,(nzeta-2)
!          do i=2,(nx-3)
!             do j=2,(nz-3)
!                gradpar(i,j,k)=b0x(i,j)/b0(i,j)*(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/dx  &
!                +b0z(i,j)/b0(i,j)*(matrix(i,j+1,k)-matrix(i,j-1,k))*0.5/dz                                 &
!                +b0zeta(i,j)/b0(i,j)*(matrix(i,j,k+1)-matrix(i,j,k-1))*kmx/(Rgrid(i)*4*pi)
!             enddo
!          enddo
!       enddo

!          do i=2,(nx-3)
!             do j=2,(nz-3)
!                gradpar(i,j,0)=b0x(i,j)/b0(i,j)*(matrix(i+1,j,0)-matrix(i-1,j,0))*0.5/dx  &
!                +b0z(i,j)/b0(i,j)*(matrix(i,j+1,0)-matrix(i,j-1,0))*0.5/dz                        &
!                +b0zeta(i,j)/b0(i,j)*(matrix(i,j,1)-matrix(i,j,nzeta))*kmx/(Rgrid(i)*4*pi)
!             enddo
!          enddo


!          do i=2,(nx-3)
!             do j=2,(nz-3)
!                gradpar(i,j,nzeta-1)=b0x(i,j)/b0(i,j)*(matrix(i+1,j,nzeta-1)-matrix(i-1,j,nzeta-1))*0.5/dx  &
!                +b0z(i,j)/b0(i,j)*(matrix(i,j+1,nzeta-1)-matrix(i,j-1,nzeta-1))*0.5/dz                        &
!                +b0zeta(i,j)/b0(i,j)*(matrix(i,j,0)-matrix(i,j,nzeta-2))*kmx/(Rgrid(i)*4*pi)
!             enddo
!          enddo
      
!  !     return
!       end function gradpar


! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!       function gradper2(matrix)

!  !     use gem_com
!  !     use equil
!       real,dimension(0:nx,0:nz,0:nzeta)::gradper2,matrix
! !      real(nx,nz,nzeta+4)::ghostmetric
!       integer i,j,k
!       do k=0,nzeta-1
!          do i=2,nx-3
!             do j=2,nz-3
!                gradper2(i,j,k)=(matrix(i+1,j,k)+matrix(i-1,j,k)-2*matrix(i,j,k))/dx**2   &
!                               +(matrix(i,j+1,k)+matrix(i,j-1,k)-2*matrix(i,j,k))/dz**2            &
!                               +(matrix(i+1,j,k)-matrix(i-1,j,k))*0.5/(dx*Rgrid(i))
!             enddo
!          enddo
!       enddo

!   !    return
!       end function gradper2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!      function gkpoisson(matrix)

 !     real,dimension(nx,nz,nzeta)::matrix, gkpoisson

  !    gkpoisson=0!matrix
   !      return
   !   end function gkpoisson


 
  
 !end subroutine gkps
      
      
!      End of gkps....
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine eqmo(ip)
      use gem_com
      use equil
      implicit none
      integer :: i,j,k,ip
      real(8) :: eta

      ez(:,:,:) = 0.
      if(iez==0)return

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine spec(n)
      use gem_com
      use equil
      implicit none
      integer :: i,j,k,l,m,n

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ezamp

      use gem_com
      use equil

      implicit none

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function ran2(idum)
      parameter( IM1=2147483563,  &
                IM2=2147483399, &
                AM=1.0/IM1,&
                IMM1=IM1-1,&
                IA1=40014,&
                IA2=40692,&
                IQ1=53668,&
                IQ2=52774,&
                IR1=12211,&
                IR2=3791,&
                NTAB=32,&
                NDIV=1+IMM1/NTAB,&
                EPS=1.2e-7,&
                RNMX=1.0-EPS &
               )
      integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
      real(8) :: temp

      save idum2, iy,iv
!      write(*,*)'idum2,iy  ',idum2,iy
      if(idum.le.0)then
         if(-idum.lt.1)then
            idum=1
         else
            idum = -idum
         end if
         idum2 = idum
         do j = NTAB+7,0,-1
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0)idum = idum+IM1
            if(j.lt.NTAB)iv(j) = idum
         end do
         iy = iv(0)
      end if

      k = idum/IQ1
      idum = IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0)idum = idum+IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2)-k*IR2
      if(idum2.lt.0)idum2 = idum2+IM2
      j = iy/NDIV
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy<1)iy = iy+IMM1
      temp = AM*iy
      if(temp>RNMX)then
         ran2 = RNMX
      else
         ran2 = temp
      end if
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadi

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1
      REAL(8) :: vpar,vperp2,r,x,z,b,ter,bfldp
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: wx0,wx1,wz0,wz1

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

!      m = 0
      m=1
      do while(m<=mm(1))
!     load a slab of ions...

!         dumx=xdim*ran2(iseed)  !revers(MyId*cnt+j,2) !ran2(iseed)
         dumx=0.85*xdim+0.13*xdim*ran2(iseed)  !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=0.55*zdim+0.1*zdim*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=pi2*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         r = xctr-xdim/2+dumx
         jacp = r/(xctr+xdim/2)
         if(ran2(iseed)<jacp)then
            x2(m)=min(dumx,xdim-1.d-8)
            z2(m)=min(dumy,zdim-1.d-8)
            zeta2(m)=dumz

            call parperp(vpar,vperp2,m,pi,cnt,MyId)

            x=x2(m)
            i = int(x/dxeq)
            wx0 = ((i+1)*dxeq-x)/dxeq
            wx1 = 1.-wx0

            z = z2(m)
            k = int(z/dzeq)
            wz0 = ((k+1)*dzeq-z)/dzeq
            wz1 = 1-wz0

            bfldp = wx0*wz0*b0(i,k)+wx0*wz1*b0(i,k+1) &
                   +wx1*wz0*b0(i+1,k)+wx1*wz1*b0(i+1,k+1) 
            ter = wx0*wz0*t0i(i,k)+wx0*wz1*t0i(i,k+1) &
                   +wx1*wz0*t0i(i+1,k)+wx1*wz1*t0i(i+1,k+1) 

            u2(m)=vpar/sqrt(mims(1)/ter)
            mu(m)=0.5*vperp2/bfldp*ter

            myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
            w2(m)=2.*amp*ran2(iseed)

            myavgw=myavgw+w2(m)
            m = m+1            
         end if
      end do

      myavgw = myavgw/mm(1)

      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(1)
         u2(m)=u2(m)-avgv
         x3(m)=x2(m)
         z3(m)=z2(m)
         zeta3(m)=zeta2(m)
         u3(m)=u2(m)
!         w2(m) = w2(m)-myavgw
         w3(m)=w2(m)
 180  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradu(u,ux,uz)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: ux(0:imx,0:jmx,0:kmx),uz(0:imx,0:jmx,0:kmx)
      integer :: i,j,k,l,m,n,jj,ju,jl
      real(8) :: ydum,wy1,ul

      do j=0,jmx-1
         ju = j+1
         jl = j-1
         if(j.eq.0)jl = jmx-1
         do i=0,imx-1
            do k=0,kmx
               uz(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dz)
            enddo
         enddo
      enddo

      do i=1,imx-1
         do j=0,jmx-1
            do k=0,kmx
               ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
            enddo
         enddo
      enddo

! do boundary i=0
      do j=0,jmx-1
         do k=0,kmx
            ul=u(imx-1,j,k)
            ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
         enddo
      enddo

      return
      end








!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradz(u,uz)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:kmx)
      real(8) :: uz(0:imx,0:jmx,0:kmx)
      integer :: i,j,k,kleft,kright
      real(8) :: wx0,wx1,wz0,wz1,uleft,uright

      uz = 0.
      do k = 0,kmx-1
         kleft = k-1
         if(k==0)kleft=kmx-1
         kright = k+1
         do i = 1,imx-1
            do j=1,jmx-1
               wx0 = ((ileft(i,j)+1)*dx-xg(i))/dx
               wx1 = 1.0-wx0
               wz0 = ((jleft(i,j)+1)*dz-zg(i))/dz
               wz1 = 1.0-wz0
               uleft = wx0*wz0*u(ileft(i,j),jleft(i,j),kleft) &
                      +wx1*wz0*u(ileft(i,j)+1,jleft(i,j),kleft) &
                      +wx0*wz1*u(ileft(i,j),jleft(i,j)+1,kleft) &
                      +wx1*wz1*u(ileft(i,j)+1,jleft(i,j)+1,kleft)
               wx0 = ((iright(i,j)+1)*dx-xg(i))/dx
               wx1 = 1.0-wx0
               wz0 = ((jright(i,j)+1)*dz-zg(i))/dz
               wz1 = 1.0-wz0
               uright = wx0*wz0*u(iright(i,j),jright(i,j),kright) &
                      +wx1*wz0*u(iright(i,j)+1,jright(i,j),kright) &
                      +wx0*wz1*u(iright(i,j),jright(i,j)+1,kright) &
                      +wx1*wz1*u(iright(i,j)+1,jright(i,j)+1,kright)

               uz(i,j,k)=(uright-uleft)/(2.*dzeta)
            enddo
         enddo
      enddo
      uz(:,:,kmx) = uz(:,:,0)

      return
      end






      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initialize
         use gem_com
      use equil

	implicit none
        real(8) :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
!        complex(8),dimension(0:1) :: x,y
        real(8),dimension(0:1) :: x,y
	integer :: n,i,j,k,ip
        
        call ppinit(MyId,numprocs,ntube,TUBE_COMM,GRID_COMM)

!     reset timestep counter.
         Last=numprocs-1
         timestep=0
         tcurr = 0.

         do i=0,Last
            if (MyId.eq.i) call init
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         enddo
         
         dum = 0.
         do i = 0,imx-1
            dum = dum+(jac(i)+jac(i+1))/2
         end do
         call MPI_ALLREDUCE(dum,jacp,1,  &
             MPI_REAL8,MPI_SUM,           &
             tube_comm,ierr)
         totvol = lx*lz*pi2*xctr    
         n0=float(tmm(1))/totvol


         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         ncurr = 1

	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
         use gem_com
         use equil
	implicit none
        integer :: n,i,j,k,ip,i1
        real(8) :: lbfr(0:imx,0:jmx)
        real(8) :: lbfs(0:imx,0:jmx)
        real(8) :: rbfr(0:imx,0:jmx)
        real(8) :: rbfs(0:imx,0:jmx)
        real(8) :: dum
        real(8) :: myrmsphi,rmp(20),myavap(0:imx-1)

	call grad(ip)

        call eqmo(ip)

end subroutine field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pintef
      use gem_com
      use equil
      implicit none

      integer :: i,j,k,ip
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr,ddedt
      real*8 :: ux(0:imx,0:jmx,0:kmx),uz(0:imx,0:jmx,0:kmx)

      phis = phi
      denes = dene
      apars = apar

      call gradz(upar,uz)
      do k = 0,kmx-1
         do i = 1,imx-1
            do j = 1,jmx-1
               ddedt = -uz(i,j,k)*gn0e(i,j)*gbtor(i,j)/((xctr-xdim/2+xg(i))*bmag(i,j))  &
                +gn0e(i,j)*(gcpnex(i,j)*ez(i,j,k)-gcpnez(i,j)*ez(i,j,k))/bmag(i,j)
               dene(i,j,k) = denes(i,j,k)+0.5*dt*ddedt
            end do
         end do
      end do

      call gradz(phi,uz)
      do k = 0,kmx-1
         do i = 1,imx-1
            do j = 1,jmx-1
               ddedt = -uz(i,j,k)*gbtor(i,j)/((xctr-xdim/2+xg(i))*bmag(i,j))
               apar(i,j,k) = apars(i,j,k)+0.5*dt*(ddedt+ezeta(i,j,k))
            end do
         end do
      end do

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cintef(n)
      use gem_com
      use equil
      implicit none

      integer :: i,j,k,n,ip
      real*8 :: tmpa(0:imx,0:jmx,0:1),tmpd(0:imx,0:jmx,0:1),tmp(0:imx,0:jmx,0:1)
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr
      REAL*8 :: myrmsapa
      real*8 :: dmnl1(0:imx,0:jmx,0:1),dmnl2(0:imx,0:jmx,0:1),dmnl3(0:imx,0:jmx,0:1),dmnl4(0:imx,0:jmx,0:1)
      real*8 :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1),exnl(0:imx,0:jmx,0:1)
      real(8) :: mydbr(0:imx-1),v(0:imx-1)


      
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine weight
  
      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m=10,i1,j1
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: x,z,zeta,dzt,rdum,dum1,wz0,wz1,wx0,wx1
      real(8) :: bp,bxp,bzp,radiusp,btorp
      
      dzt = dzeta/float(m)
      do i = 1,imx-1
         x = xg(i)
         do j = 1,jmx-1
            z = zg(j)
            do k = 0,m-1
               i1 = int(x/dxeq)
               j1 = int(z/dzeq)

               i1 = int(x/dxeq)
               wx0 = ((i1+1)*dxeq-x)/dxeq
               wx1 = 1.-wx0

               j1 = int(z/dzeq)
               wz0 = ((j1+1)*dzeq-z)/dzeq
               wz1 = 1-wz0

               bp = wx0*wz0*b0(i1,j1)+wx0*wz1*b0(i1,j1+1) &
                 +wx1*wz0*b0(i1+1,j1)+wx1*wz1*b0(i1+1,j1+1) 
               bxp = wx0*wz0*b0x(i1,j1)+wx0*wz1*b0x(i1,j1+1) &
                 +wx1*wz0*b0x(i1+1,j1)+wx1*wz1*b0x(i1+1,j1+1) 
               bzp = wx0*wz0*b0z(i1,j1)+wx0*wz1*b0z(i1,j1+1) &
                 +wx1*wz0*b0z(i1+1,j1)+wx1*wz1*b0z(i1+1,j1+1) 
               btorp = wx0*wz0*b0zeta(i1,j1)+wx0*wz1*b0zeta(i1,j1+1) &
                 +wx1*wz0*b0zeta(i1+1,j1)+wx1*wz1*b0zeta(i1+1,j1+1) 
               radiusp = xctr-xdim/2+x               
               x = x+dzt*radiusp*bxp/btorp
               x = min(x,lx)
               x = max(x,0.)               
               z = z+dzt*radiusp*bzp/btorp
               z = min(z,lz)
               z = max(z,0.)
            end do
            iright(i,j) = int(x/dx)
            jright(i,j) = int(z/dx)
         end do
      end do

      dzt = -dzeta/float(m)
      do i = 1,imx-1
         x = xg(i)
         do j = 1,jmx-1
            z = zg(j)
            do k = 0,m-1
               i1 = int(x/dxeq)
               j1 = int(z/dzeq)

               i1 = int(x/dxeq)
               wx0 = ((i1+1)*dxeq-x)/dxeq
               wx1 = 1.-wx0

               j1 = int(z/dzeq)
               wz0 = ((j1+1)*dzeq-z)/dzeq
               wz1 = 1-wz0

               bp = wx0*wz0*b0(i1,j1)+wx0*wz1*b0(i1,j1+1) &
                 +wx1*wz0*b0(i1+1,j1)+wx1*wz1*b0(i1+1,j1+1) 
               bxp = wx0*wz0*b0x(i1,j1)+wx0*wz1*b0x(i1,j1+1) &
                 +wx1*wz0*b0x(i1+1,j1)+wx1*wz1*b0x(i1+1,j1+1) 
               bzp = wx0*wz0*b0z(i1,j1)+wx0*wz1*b0z(i1,j1+1) &
                 +wx1*wz0*b0z(i1+1,j1)+wx1*wz1*b0z(i1+1,j1+1) 
               btorp = wx0*wz0*b0zeta(i1,j1)+wx0*wz1*b0zeta(i1,j1+1) &
                 +wx1*wz0*b0zeta(i1+1,j1)+wx1*wz1*b0zeta(i1+1,j1+1) 
               radiusp = xctr-xdim/2+x               
               x = x+dzt*radiusp*bxp/btorp
               x = min(x,lx)
               x = max(x,0.)               
               z = z+dzt*radiusp*bzp/btorp
               z = min(z,lz)
               z = max(z,0.)
            end do
            ileft(i,j) = int(x/dx)
            jleft(i,j) = int(z/dx)
         end do
      end do

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



 
       subroutine ComputeInitialGuess(ksp,b,ctx,ierr)
       use petscksp
       implicit none
       PetscErrorCode  ierr
       KSP ksp
       PetscInt ctx(*)
       Vec b
       PetscScalar  h

       h=0.0
       PetscCall(VecSet(b,h,ierr))
       end subroutine
       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ComputeMatrix(ksp,A,B,dummy,ierr)
      use petscksp
      use gem_com
      use equil
       implicit none
       PetscErrorCode  ierr
       KSP ksp
       Mat A,B
       integer dummy(*)
       DM dm
       integer :: ii,jj

      PetscInt    i,j,mx,my,xm
      PetscInt    ym,xs,ys,i1,i5
      PetscScalar  v(5),Hx,Hy
      PetscScalar  Hx2,Hy2,tmp_r,a_value
      MatStencil   row(4),col(4,5)

      i1 = 1
      i5 = 5
      a_value = 0.5
      PetscCall(KSPGetDM(ksp,dm,ierr))
      PetscCall(DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))

      Hx = (Rgrid(imx)-Rgrid(0)) / real(imx)
      Hy = (Zgrid(jmx)-Zgrid(0)) / real(jmx)
      Hx2 = Hx**2
      Hy2 = Hy**2
      PetscCall(DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr))
      do 10,j=ys,ys+ym-1
         do 20, i=xs,xs+xm-1
           ! i=ii
           ! j=jj
          row(MatStencil_i) = i
          row(MatStencil_j) = j
!          tmp_r = sqrt(X_mat(j+1, i+1)**2 + Y_mat(j+1, i+1)**2)
          if (mask(i,j) <0.99) then
             v(1)=c2_over_vA2(i,j)*(-2.0/Hx2-2.0/Hy2)
             PetscCall(MatSetValuesStencil(B,i1,row,i1,row,v,INSERT_VALUES,ierr))
          else
             if (j > 0) then
                 if (j == jmx) then
                    v(1) = c2_over_vA2(i,j)/Hy2-1.0/(2.0*Hy2)*( c2_over_vA2(i,j)- c2_over_vA2(i,j-1))
                 else
                    v(1) = c2_over_vA2(i,j)/Hy2-1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1))
                 end if
             end if
             col(MatStencil_i, 1) = i
             col(MatStencil_j, 1) = j - 1    

             if (i > 0) then
                 if (i == imx) then
                    v(2) =  c2_over_vA2(i,j)/Hx2-1.0/(2.0*Hx2)*( c2_over_vA2(i,j)- c2_over_vA2(i-1,j))
                 else
                    v(2) =  c2_over_vA2(i,j)/Hx2-1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j))
                 end if
             end if
             col(MatStencil_i, 2) = i - 1
             col(MatStencil_j, 2) = j

             v(3) = -2.0* c2_over_vA2(i,j) / Hx2 - 2.0* c2_over_vA2(i,j) / Hy2
             col(MatStencil_i, 3) = i
             col(MatStencil_j, 3) = j

             if (i < imx) then
                 if (i == 0) then
                    v(4) =  c2_over_vA2(i,j)/Hx2+1.0/(2.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i,j))
                 else
                    v(4) =  c2_over_vA2(i,j)/Hx2+1.0/(4.0*Hx2)*( c2_over_vA2(i+1,j)- c2_over_vA2(i-1,j))
                 end if
             end if
             col(MatStencil_i, 4) = i + 1
             col(MatStencil_j, 4) = j

             if (j < jmx) then
                 if (j == 0) then
                    v(5) =  c2_over_vA2(i,j)/Hy2+1.0/(2.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j))
                 else
                    v(5) =  c2_over_vA2(i,j)/Hy2+1.0/(4.0*Hy2)*( c2_over_vA2(i,j+1)- c2_over_vA2(i,j-1))
                 end if
             end if
             col(MatStencil_i, 5) = i
             col(MatStencil_j, 5) = j + 1
 

!             v(1) = eps(i+1,j+1)/Hy2
!             col(MatStencil_i, 1) = i
!             col(MatStencil_j, 1) = j - 1
!             v(2) = eps(i+1,j+1)/Hx2
!             col(MatStencil_i, 2) = i - 1
!             col(MatStencil_j, 2) = j
!             v(3) = eps(i+1,j+1)*(-2.0/Hx2-2.0/Hy2)
!             col(MatStencil_i, 3) = i
!             col(MatStencil_j, 3) = j
!             v(4) = eps(i+1,j+1)/Hx2
!             col(MatStencil_i, 4) = i + 1
!             col(MatStencil_j, 4) = j
!             v(5) = eps(i+1,j+1)/Hy2
!             col(MatStencil_i, 5) = i
!             col(MatStencil_j, 5) = j + 1
             PetscCall(MatSetValuesStencil(B, i1, row, i5, col, v, INSERT_VALUES, ierr))
          endif

!       enddo
!    enddo
    
20       continue
10    continue
      PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr))
      PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr))
      if (A .ne. B) then
         PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr))
         PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr))
      endif
!      PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD,petsc_ierr))
!      PetscCall(MatView(B,PETSC_VIEWER_STDOUT_WORLD,petsc_ierr))
    end subroutine

    
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


       subroutine ComputeRHS(ksp,b,k,ierr)
       use petscksp
       use gem_com
       use equil
       implicit none
       integer::k,ii,jj
       PetscErrorCode  ierr
       KSP ksp
       Vec b
!       integer dummy(*)
       PetscScalar  h,Hx,Hy
       PetscInt  mx,my,i,j
       DM dm
       PetscInt idx
       PetscScalar tmp_value,a_value,tmp_r

       PetscCall(KSPGetDM(ksp,dm,ierr))
       PetscCall(DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))

!       Hx = Lx / real(mx-1)
!       Hy = Ly / real(my-1)
      ! h = Hx*Hy
      ! print *, 'h=',h
!       a_value = 0.5
       do i = 0,mx 
          do j = 0,my
             ii=i
             jj=j
             idx = j*(imx+1)+i
 !            tmp_r = sqrt(X_mat(j+1, i+1)**2 + Y_mat(j+1, i+1)**2)
             if (mask(i,j) <0.99) then
                tmp_value = 0
             else
                 tmp_value = -dene(ii,jj,k)
             end if    
             PetscCall(VecSetValues(b, 1, idx, tmp_value, INSERT_VALUES, ierr))
          end do
       end do
       !PetscCall(VecView(b,PETSC_VIEWER_STDOUT_WORLD,petsc_ierr))
       
       end subroutine
