   subroutine readS
   use smat
   implicit none
   integer i,io,j
   
   ! read in the smatrix from the file
   open(unit=13,file="stot.txt")
   
   ! find out how many lines are in the file
   nl = 0
   do 
      read(13,*,iostat=io)
      if (io/=0) exit
      nl = nl + 1
   enddo 
   ! actual number of lines
   nl = nl/(2*schan + 1)
   !print *, "number of energies!",nl
   
   ! allocate the matrices
   allocate (en_ke(2,nl))
   allocate (S(schan+1,2*schan,nl))
   
   rewind(13)
   
   ! read in the smatrix - real and imaginary parts
   do i=1,nl
      do j=1,2*schan+1
         if (j==1) then
            read(13,*) en_ke(1,i),en_ke(2,i)
	    !print *, en_ke(:,i)
	 else
	    read(13,*) S(:,j-1,i)
	    !print *, S(2,j-1,i)
	 end if 
      enddo
   enddo 
   close(13)
   !print *, S(:,1,1)
   
   end subroutine