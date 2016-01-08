   ! subroutines for single channel width calculation
   
   subroutine scpot(wf,pot,a1,a2,hbc,m)
   implicit none
   real*8, intent(in) :: wf(2,500)
   real*8, intent(in) :: hbc,m
   integer, intent(in) :: a1,a2
   real*8, intent(inout) :: pot(500)
   real*8 Vo,ro,ao,l,onepe,mch
   integer i
   
   ! constants in the potential
   Vo = -61.1
   ro = 1.2*a2**0.3333333
   ao = 0.65
   l = 2.d0
   
   mch = 2*m/hbc**2
   mch = 0.0478450
   print *, mch
   
   ! construct the potential for each value of r - wf(1,i)
   open(unit=3,file="scpot.txt")
   do i=1,500
      onepe = 1.d0 + exp((wf(1,i)-ro)/ao)
      pot(i) = Vo/onepe !+ l*(l+1)/(mch*wf(1,i)**2)
      write(3,*) wf(1,i),pot(i)
   enddo 
   
   end subroutine
   
   real function hm(mr,L)
   real*8, intent(in) :: mr,L
   double precision, dimension(101) :: fc,gc,fcp,gcp
   integer ifail,m1
   ifail=0
   m1=1
   
   call coulfg(mr,0d0,L,L,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   hm = sqrt(fc(3)**2 + gc(3)**2)
   return
   end function hm
   
   subroutine scsimpsons(pot,wf,ps,mom)
   implicit none
   real*8, intent(in) :: pot(500),wf(2,500),mom
   real*8, intent(inout) :: ps
   integer i,j,ifail,m1
   real*8 fn(3),l,hm,mr
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   L=2.d0
   
   ! integral using Simpson's Rule
   ! eg mathworld.wolfram.com/SimpsonsRule.html
   ! \int _xo ^x2 f(x)dx = \int _xo ^xo+2h f(x)dx = (h/3)(f(xo)+4f(x1)+f(x2))
   
   ! ifail and m1 must be defined, you can't just put in 0,1
   ! for whatever reason
   ! value is put into fc(L+1)
   
   do i=1,500-2,2 ! shorter when the potential is read in b/c U(r->inf)=0
      do j=1,3
         call coulfg(mom*wf(1,i+j-1),0d0,L,L,fc,gc,fcp,gcp,3,0,ifail,m1)
         fn(j) = fc(3)*pot(i+j-1)*sqrt(wf(2,i+j-1))!/wf(1,i+j-1)
	 !mr = mom*wf(1,i+j-1)
	 !fn(j) = hm(mr,L)*pot(i)*wf(2,i+j-1)/wf(1,i+j-1)
	 !print *, fc(3)
	 !print *, fc(3),fc(2)!,mom,wf(1,i+j-1)
      enddo 
      ps = ps + ((wf(1,i+1)-wf(1,i))/3.d0)*(fn(1)+4*fn(2)+fn(3))
      !print *, ps,i,wf(1,i)*pot(i)
      print *, wf(1,i),ps
   enddo 
   print *, "&"
   
   !print *, fc
   
   end subroutine 
   
   subroutine Gcalc(wf,hbc,mom,mu)
   implicit none
   real*8, intent(in) :: wf(2,500),hbc,mom,mu
   real*8 L,gamma,RG
   integer i,ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   
   ! calculate |u/G|^2 for every value of r
   
   ifail=0
   m1=1
   L=2.0d0
   
   !gamma = (hbc*k/mu)*|u/G|^2
   do i=1,500
      call coulfg(mom*wf(1,i),0d0,L,L,fc,gc,fcp,gcp,1,0,ifail,m1)
      RG = wf(2,i)/gc(3)**2
      print *, wf(1,i),hbc*mom*RG/mu,gc(3)**2
   enddo 
   
   end subroutine