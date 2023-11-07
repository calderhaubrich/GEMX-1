MODULE equil
      IMPLICIT NONE
      real :: mimp=2,chgi=1
      real :: beta,rmaj0,a,q0,r0,q0p,q0abs,shat0
      real :: dR,dth
      integer :: nr=200,nr2=100,ntheta=200,isgnf=1,isgnq=-1,isupae0=0
      real :: psi_max=0.32, psi_min=-1.0,R_min=1.0, Z_min=-1.5, Z_internal=-1.2, psi_div=0.305

!     GEM-X
      integer :: nx=449,nz=433,nzeta=32
      real :: xdim,zdim,xctr,zctr,dxeq,dzeq
      real,dimension(:,:),allocatable :: b0,b0x,b0z,b0zeta,dbdx,dbdz,c2_over_vA2
      real,dimension(:,:),allocatable :: t0i,t0e,xn0i,xn0e,captix,captex,capnix,capnex,captiz,captez,capniz,capnez
      
      real,dimension(:),allocatable :: psi,f,psip,sf, &
                                      vpari,vparip, &
                                      zeff,nue0,phinc,phincp,&
                                      er,upari, Rgrid, Zgrid
      real,dimension(:,:),allocatable :: t0s,xn0s,capts,capns,vpars,vparsp,psi_p,mask
      real:: bu,tu,nu,xu,frequ,vu,eru

!     for including bstar effects
      real,dimension(:),allocatable :: psip2
      real,dimension(:,:),allocatable :: curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb

!     for equilibrium current                                                                                                                                         
      real(8),dimension(:,:),allocatable :: upae0,nuob,dnuobdr,dnuobdt
      

contains
      subroutine new_equil()
      use gem_com,only: myid
      implicit none
      real(8) :: pi,pi2,r,th,s

      integer :: i,j,k,m,i1,j1,j2
      real(8) :: dum,x
      real(8) :: e,proton,omegau
 
!     global equilibrium data

      allocate(sf(0:nr), psi(0:nr),&
         psip(0:nr),&
         zeff(0:nr),nue0(0:nr),&
         vpari(0:nr),phinc(0:nr), &
         vparip(0:nr),phincp(0:nr),er(0:nr), &
         upari(0:nr))

      allocate(t0s(1:5,0:nr),xn0s(1:5,0:nr),&
         capts(1:5,0:nr),capns(1:5,0:nr),vpars(1:5,0:nr),&
         vparsp(1:5,0:nr))

      allocate(upae0(0:nr,0:ntheta),nuob(0:nr,0:ntheta),dnuobdr(0:nr,0:ntheta),dnuobdt(0:nr,0:ntheta))      

!     GEM-X
      allocate(b0(0:nx,0:nz),b0x(0:nx,0:nz),b0z(0:nx,0:nz),b0zeta(0:nx,0:nz),dbdx(0:nx,0:nz),dbdz(0:nx,0:nz))
      allocate(t0i(0:nx,0:nz),t0e(0:nx,0:nz),xn0i(0:nx,0:nz),xn0e(0:nx,0:nz), &
               capnix(0:nx,0:nz),capnex(0:nx,0:nz),captix(0:nx,0:nz),captex(0:nx,0:nz), &
               capniz(0:nx,0:nz),capnez(0:nx,0:nz),captiz(0:nx,0:nz),captez(0:nx,0:nz))      

               allocate(Rgrid(0:nx),Zgrid(0:nz))
               allocate(psi_p(0:nx,0:nz),mask(0:nx,0:nz),c2_over_vA2(0:nx,0:nz))

      open(unit=10, file = 'R.dat',status='old',action='read')
	read(10,*) Rgrid
      close(10)    

       open(unit=10, file = 'Z.dat',status='old',action='read')
	read(10,*) Zgrid
        close(10)

        open(unit=10, file='psi_p.dat',status='old',action='read')
        read(10,*) psi_p
        close(10)
        


        open(unit=11, file = 'psi_test',status='unknown',action='write')
 !     do i=0,nx
         write(11,*) psi_p
 !     end do
      close(11)     
      
            
        
        xdim=Rgrid(nx)-Rgrid(0)
        zdim=Zgrid(nz)-Zgrid(0)
        dR=(Rgrid(nx)-Rgrid(0))/nx
!        dx=(Rgrid(nx)-Rgrid(0))/nx
!        dZ=(Zgrid(nz)-Zgrid(0))/nz
      
      open(unit=10, file = 'BR.dat',status='old',action='read')
	read(10,*) b0x
      close(10)    

!      open(unit=11, file = 'debug.dat',status='unknown',action='write')
!      write(11,*) b0x(5,0), 'dR='dR
!      close(11)
      
      open(10, file = 'Bz.dat',status='old',action='read')
	read(10,*) b0z
      close(10)
      
      open(unit=10, file = 'Bt.dat',status='old',action='read')
	read(10,*) b0zeta
      close(10)


      
      do i =0,nx-1
            do j =0,nz-1
		b0(i,j) = sqrt(b0x(i,j)**2+b0z(i,j)**2+b0zeta(i,j)**2)
    	end do
      end do 



          do i=0,nx
           do j=0,nz
              c2_over_vA2(i,j)=2*Rgrid(i)**2/(Rgrid(0)+Rgrid(nx))**2
!              if (psi_p(i,j)<psi_max .and.psi_p(i,j)>psi_min .and. (Zgrid(j)>Z_internal .or. (Zgrid(j)>Z_min .and. psi_p(i,j)>psi_div)) .and. Rgrid(i)>R_min)  then
              if (psi_p(i,j) < 0.3.and.Zgrid(j)>-1.08) then
                  mask(i,j)=1
              else
                  mask(i,j)=0
              endif
        end do
     end do
     open(unit=11, file = './out/mask.dat',status='unknown',action='write')
 !     do i=0,nx
         write(11,*) mask
 !     end do
      close(11)     

!     Normalization
      e = 1.6e-19
      proton = 1.67e-27
      Bu = 1.98562
      Tu = 1000*e
      omegau = e*Bu/proton
      frequ = omegau

      vu = sqrt(Tu/proton)
      xu = proton*vu/(e*Bu)
      nu = 2.5e19
      beta = 4*3.14159*1e-7*nu*Tu/Bu**2 

      b0 = b0/bu
      b0x = b0x/bu
      b0z = b0z/bu
      b0zeta = b0zeta/bu
      
!     write(*,*) 'betaU=', beta

      pi = atan(1.0)*4      
      rmaj0 = 1000.
      a = 360.

!      xctr = a*1.5
!      zctr = 0.


      xctr =  0.5*(Rgrid(nx)+Rgrid(0))/xu
      zctr  = 0.5*(Zgrid(nz)+Zgrid(0))/xu
      xdim=xdim/xu
      zdim=zdim/xu
 !     xdim = a*2
 !     zdim = a*3
      dxeq = xdim/nx
      dzeq = zdim/nz

!     assign T, n profiles 
      do i = 0,nx
         do j = 0,nz
            t0i(i,j) = 1.
            t0e(i,j) = 1.
            xn0i(i,j) = 1.
            xn0e(i,j) = 1.
         end do
      end do
      
      
!     Gradients of T and n profiles for global simulation. Flux-tube parameters done at the end
      captix = 0.
      captex = 0.
      capnix = 0.
      capnex = 0.
      captiz = 0.
      captez = 0.
      capniz = 0.
      capnez = 0.

      dbdx = 0.
      dbdz = 0.
      do i = 1,nx-1
         do j = 1,nz-1
            dbdx(i,j) = (b0(i+1,j)-b0(i-1,j))/(2*dxeq)
            dbdz(i,j) = (b0(i,j+1)-b0(i,j-1))/(2*dzeq)            
         end do
      end do

    if(myid==0)then
       open(11,file='xpp',status='replace')
       write(11,*)'xu,omegau,vu = ', xu, omegau, vu

10     format(1x,i5,10(1x,e16.9))
       do i = 0,nx
          write(11,10)i,b0(i,nz/2),b0x(i,nz/2),b0z(i,nz/2),b0zeta(i,nz/2),dbdx(i,nz/2),dbdz(i,nz/2)
       end do
       close(11)
      end if


 !     open(unit=11, file = 'debug.dat',status='unknown',action='write')
 !     write(11,*) dR
 !     close(11)

      
  end subroutine new_equil

            
END MODULE equil
