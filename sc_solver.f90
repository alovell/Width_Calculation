   real function func(l,r,x)
   implicit none
   real*8, intent(in) :: l,r,x
   real*8 ao,ro,vo,eng,mch
   
   ! energy parameters
   mch = 0.047845
   eng = 1.38
   vo= -61.1
   ro = 1.2*10**0.33333333
   ao = 0.65
   
   func = (mch*(vo/(1.d0 + exp((r-ro)/ao)) - eng) + l*(l+1)/(r**2))*x
   
   end function
   
   real function vr(r)
   implicit none
   real*8, intent(in) :: r
   real*8 ao,ro,vo
   
   vo= -61.1
   ro = 1.2*10**0.33333333
   ao = 0.65
   
   vr = vo/(1.d0 + exp((r-ro)/ao))
   
   end function
   
   subroutine simpson(fx,n,h,ps)
   implicit none
   integer, intent(in) :: n
   real*8, intent(in) :: fx(n),h
   real*8, intent(inout) :: ps
   real*8 fn(3),r(n)
   integer i,j
   
   ! create the radial matrix
   do i=1,n
      r(i) = i*h
   enddo 
   
   ! do the integral with Simpson's Rule
   ps = 0.d0
   do i=1,n-2,2
      do j=1,3
         fn(j) = (fx(i+j-1)*r(i+j-1))**2
      enddo 
      ps = ps + ((r(i+1)-r(i))/3.d0)*(fn(1)+4*fn(2)*fn(3))
   enddo 
   
   end subroutine
   
   subroutine simpson2(fx,n,h,ps)
   implicit none
   integer, intent(in) :: n
   real*8, intent(in) :: fx(n),h
   real*8, intent(inout) :: ps
   real*8 fn(3),mch,eng,r(n),l,vr
   integer i,j,ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   
   mch = 0.047845
   eng = 1.38
   l=2.d0
      
   ! create the radial matrix
   do i=1,n
      r(i) = i*h
   enddo 
   
   ! do the integral with Simpson's Rule
   ps = 0.d0
   do i=1,n-2,2
      do j=1,3
         call coulfg(sqrt(eng*mch)*r(i+j-1),0d0,l,l,fc,gc,fcp,gcp,1,0,ifail,m1)
         fn(j) = fx(i+j-1)*fc(3)*vr(r(i+j-1))
	 !print *, vr(r(i+j-1))
      enddo 
      ps = ps + ((r(i+1)-r(i))/3.d0)*(fn(1)+4*fn(2)*fn(3))
      !print *, ps
   enddo 
   
   end subroutine
   
   subroutine simpsonc(fx,n,h,ps)
   implicit none
   integer, intent(in) :: n
   complex*8, intent(in) :: fx(n)
   real*8, intent(in) :: h
   real*8, intent(inout) :: ps
   real*8 mch,eng,r(n),l,vr
   complex*8 fn(3),psc
   integer i,j,ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   
   mch = 0.047845
   eng = 1.38
   l=2.d0
      
   ! create the radial matrix
   do i=1,n
      r(i) = i*h
   enddo 
   
   ! do the integral with Simpson's Rule
   ps = 0.d0
   do i=1,n-2,2
      do j=1,3
         call coulfg(sqrt(eng*mch)*r(i+j-1),0d0,l,l,fc,gc,fcp,gcp,1,0,ifail,m1)
         fn(j) = fx(i+j-1)*fc(3)*vr(r(i+j-1))
	 !print *, vr(r(i+j-1))
      enddo 
      psc = psc + ((r(i+1)-r(i))/3.d0)*(fn(1)+4*fn(2)*fn(3))
      !print *, ps
   enddo 
   
   ps = abs(psc)
   print *, psc, ps
   
   end subroutine
   
   program calcwf
   implicit none
   real*8 k1,k2,k3,k4,l1,l2,l3,l4,a,func,sum,norm,pi
   real*8 u(5000),r,y(5000),ao,ro,vo,eng,mch,l,uap
   complex*8 wronk,hm,hmp,uc(5000)
   integer i,ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   
   ! energy parameters
   mch = 0.047845
   eng = 1.38
   vo= -61.1
   ro = 1.2*10**0.33333333
   ao = 0.65
   l=2.d0
   
   ! pi
   pi = 4.d0*atan(1.d0)
   
   ! give the potential parameters and the initial wave function
   a = 0.01d0
   u(1) = 0.d0
   !u(2) = 0.5d0
   y(1) = 0.5d0
   !do i=1,500
   !   r(i) = h*i + i
   !enddo 
   
   ! solves a 4th order Runge Kutta to find the scattering wave function
   do i=1,5000-1
      r = i*a
      K1 = y(i)*a
      L1 = func(l, r + 0.5*a, u(i) + 0.5*K1)*a
      K2 = (y(i) + 0.5*L1)*a
      L2 = func(l, r + 0.5*a, u(i) + 0.5*K1)*a
      K3 = (y(i) + 0.5*L2)*a
      L3 = func(l, r + 0.5*a, u(i) + 0.5*K2)*a
      K4 = (y(i) + L3)*a
      L4 = func(l, r + a, u(i) + K3)*a
      u(i + 1) = u(i) + (K1 + 2*K2 + 2*K3 + K4)/6
      y(i + 1) = y(i) + (L1 + 2*L2 + 2*L3 + L4)/6
   enddo 
   
   !do i=1,5000
   !   print *, i*a,u(i)
   !enddo
   
   ! normalize the wave function in a box
   !sum = 0.d0
   !call simpson(u,3500,a,sum)
   !norm = 1.d0/sqrt(sum)
   !print *, (197.32705/846.812)*sqrt(eng*mch)*norm**2
   norm = a**3/(15*u(2))
   
   u = u*norm
   
   !u = u*0.00125
   
   do i=50,5000
      call coulfg(sqrt(eng*mch)*i*a,0d0,l,l,fc,gc,fcp,gcp,1,0,ifail,m1)
   !   print *, i*a,u(i),gc(3)
   !   print *, i*a,(u(i)/gc(3))**2
   enddo
   
   ! then calculate the ratio
   call coulfg(sqrt(eng*mch)*35.0,0d0,l,l,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   print *, (197.32705/846.812)*sqrt(eng*mch)*(u(3500)/gc(3))**2
   !print *, (u(500)/gc(3))**2
   
   print *, (197.32705/846.812)*sqrt(eng*mch)*u(3500)**2/(gc(3)**2+fc(3)**2)
   
   ! calculate the integral relation
   sum = 0.d0
   call simpson2(u,2000,a,sum)
   !print *, sum
   print *, 4*846.812/(197.32705**2*sqrt(eng*mch))*sum**2
   
   ! normalization number two
   ! normalize the Wronskian to 2*sqrt(pi*(2l+1))
   
   uap = (u(5000)-u(4998))/(2*a)
   call coulfg(sqrt(eng*mch)*49.99,0d0,l,l,fc,gc,fcp,gcp,1,0,ifail,m1)
   hm = cmplx(gc(3),-fc(3))
   hmp = cmplx(gcp(3),-fcp(3))
   wronk = u(4999)*hmp - uap*hm
   !print *, wronk
   !uc = 2*sqrt(pi*(2*l+1))*u/wronk
   print *, abs(2*sqrt(3.14159*(2*l+1))/wronk)**2*4*846.812/(197.32705**2*sqrt(eng*mch))*sum**2
   
   !sum = 0.d0
   !call simpsonc(uc,2000,a,sum)
   !print *, 4*846.812/(197.32705**2*sqrt(eng*mch))*sum**2
   !print *, uc(3),u(3)
   
   end program