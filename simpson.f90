   subroutine simpsons(ps,cwf,pots,energy,ke,hbc,mu,file,cn,N,nchan)
   implicit none
   integer, intent(in) :: cn,N,nchan
   real*8, intent(in) :: hbc,energy,mu,ke
   real*8, intent(in) :: cwf(N,nchan+1)
   real*8, intent(in) :: pots(400)
   character (len=10), intent(in) :: file
   real*8, intent(inout) :: ps
   integer i,j,ifail,m1
   real*8 fn(3),FL,U,L,temp1,temp2,K
   character (len=10) temp3
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   
   ! open given file to get K
   open(unit=2,file=file)
   rewind(2)
   read(2,*) temp1,temp2,temp3,K
   ! define L
   L = 1.5 + K
   
   ! integral using Simpson's Rule
   ! eg mathworld.wolfram.com/SimpsonsRule.html
   ! \int _xo ^x2 f(x)dx = \int _xo ^xo+2h f(x)dx = (h/3)(f(xo)+4f(x1)+2f(x2))
   
   ! ifail and m1 must be defined, you can't just put in 0,1
   ! for whatever reason
   ! value is put into fc(2)
   
   ! calculate momentum from energy
   !ke=sqrt(mu*energy)/hbc
   
   do i=1,400-2 ! shorter when the potential is read in b/c U(r->inf)=0
      do j=1,3
         call coulfg(ke*cwf(i+j-1,1),0d0,L,L,fc,gc,fcp,gcp,1,0,ifail,m1)
         fn(j) = fc(2)*pots(i)*cwf(i+j-1,cn+1)/cwf(i+j-1,1)
	 !print *, fc(2)
      enddo 
      ps = ps + ((cwf(i+1,1)-cwf(i,1))/3.d0)*(fn(1)+4*fn(2)+2*fn(3))
   enddo 
   
   !ps=0
   !do i=1,400-1
   !   call coulfg(ke*cwf(i,1),0d0,L,L,fc,gc,fcp,gcp,1,0,ifail,m1)
   !   fn(1) = fc(2)*pots(i)*cwf(i,cn+1)/cwf(i,1)
   !   ps = ps + fn(1)
   !enddo 
   !ps = (cwf(i+1,1)-cwf(i,1))*ps/3.d0
   
   end subroutine simpsons