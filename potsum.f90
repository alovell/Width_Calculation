   !subroutine potsum(cn,pots,chi,N,nchan)
   subroutine potsum(nchan)
   use channels
   use totwf
   implicit none
   !integer, intent(in) :: cn,N,nchan
   !real*8, intent(in) :: chi(N,nchan+1)
   !real*8, intent(inout) :: pots(400)
   integer, intent(in) :: nchan
   real*8 temppot(400,2),tempwf1,tempwf2
   character (len=20) file,filewf,temp3
   integer i,j
   real*8 temp1,temp2
   
   ! create file name, and read in potentials
   pots = 0
   ! everything after the channel of interest
   do i=channum,nchan
      !print *, i
      write (file, '( I0, "_", I0, ".pot" )' ) channum,i
      open(unit=i+10,file=file)
      do j=1,400
         read(i+10,*) temppot(j,1), temppot(j,2)
	 !print *, temppot(j,2),pots(j),chi(j,i+1)
	 pots(j) = pots(j) + temppot(j,2)*chi(j,i+1)
	 !print *, temppot(j,2),pots(j)
      enddo 
      close(i+10)
   enddo 
   
   ! need everything before the channel of interest
   do i=1,channum-1
      write (file, '( I0, "_", I0, ".pot" )' ) i,channum
      open(unit=i+10,file=file)
      do j=1,400
         read(i+10,*) temppot(j,1), temppot(j,2)
	 !print *, temppot(j,2),pots(j),chi(j,i+1)
	 pots(j) = pots(j) + temppot(j,2)*chi(j,i+1)
	 !print *, temppot(j,2),pots(j)
      enddo 
      close(i+10)
   enddo 
   
   end subroutine potsum