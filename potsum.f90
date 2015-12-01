   subroutine potsum(cn,pots,chi,N,nchan)
   implicit none
   integer, intent(in) :: cn,N,nchan
   real*8, intent(in) :: chi(N,nchan+1)
   real*8, intent(inout) :: pots(400)
   real*8 temppot(400,2),tempwf1,tempwf2
   character (len=20) file,filewf,temp3
   integer i,j
   real*8 temp1,temp2,K
   
   ! create file name, and read in potentials
   pots = 0
   ! everything after the channel of interest
   do i=cn,nchan
      !print *, i
      write (file, '( I0, "_", I0, ".pot" )' ) cn,i
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
   do i=1,cn-1
      write (file, '( I0, "_", I0, ".pot" )' ) i,cn
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