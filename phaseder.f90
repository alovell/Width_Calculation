! calculate the derivative of the phase 
! using two point formula (see main code, below)

subroutine deriv(ph,m,der)
   implicit none
   integer m,n,l
   real*8 ph(19,m),der(19,m)
   
   ! calculate derivative
   do n=2,m+1
      der(1,n-1)=ph(1,n)
      do l=2,19
         der(l,n-1)=(ph(l,n+1)-ph(l,n-1))/(ph(1,n+1)-ph(1,n-1))
      enddo
   enddo 
   
end subroutine deriv

! calculate the width based on the derivative of the phase shift
program phaseder
   implicit none
   integer io,i,j,k,narg
   character(len=50) infile,outfile
   real*8 tempx,tempy
   real*8, dimension(:,:), allocatable :: der,phase,phinterp1,phinterp2,derinterp1
   
   narg = IARGC()
   call getarg(1,infile)
   call getarg(2,outfile)
   
   ! read in file
   j=0
   open(unit=2,file=infile)
   do 
      read(2,*,iostat=io)
      if (io /= 0) exit
      j=j+1
   enddo
   
   allocate(phase(19,j))
   allocate(phinterp1(19,2*(j-1)))

   rewind(2)
   read(2,*) phase
   close(2)
   
   ! do the linear interpolation
   phinterp1(1,1) = phase(1,1)
   do k=2,19
      phinterp1(k,1) = phase(k,1)
   enddo 
   
   do i=1,j-1
      tempx = (phase(1,i) + phase(1,i+1))/2.d0
      phinterp1(1,2*i-1) = tempx
      phinterp1(1,2*i) = phase(1,i+1)
      do k=2,19
         tempy = phase(k,i) + (phase(k,i+1)-phase(k,i))*(tempx-phase(1,i))/(phase(1,i+1)-phase(1,i))
	 phinterp1(k,2*i-1) = tempy
	 phinterp1(k,2*i) = phase(k,i+1)
      enddo 
   enddo 
   
   ! allocate derivative matrix
   allocate (der(19,j-2))
   allocate (derinterp1(19,2*j-3))
   
   ! calculate derivative ((f(x+h) - f(x-h))/2h)
   call deriv(phase,j-2,der)
   call deriv(phinterp1,2*j-3,derinterp1)
   !open(unit=3,file="derivative.txt") 
   open(unit=3, file=outfile)
   open(unit=4, file="interp.txt")
   open(unit=7, file="dorig.txt")
   write(4,*) "# interpolation"
   write(4,20) phinterp1
   !write(4,*) ! put a blank space
   write(7,*) "# derivative of original"
   write(7,20) der
   !write(7,*) ! put a blank space
   write(3,*) "# derivative of interpolation"
   write(3,20) derinterp1   
   close(3)
   close(4)
   close(7)

20    format(19f12.6)    

   ! calculate width (maybe a little trickier)
end program phaseder