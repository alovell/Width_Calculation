   subroutine getS(K,Sr,Si,cn,cnp)
   implicit none
   real*8, intent(in) :: K
   integer, intent(in) :: cn,cnp
   real*8, intent(inout) :: Sr,Si
   integer m
   real*8 S(111)
   
   open(unit=12,file="smatrix.txt")
   rewind(12)
   read(12,*) !E and momentum  
   do m=1,220
      read(12,*) S(:)
      !print *, S(:)
      if (m == cn) then
         Sr = S(cnp+1)
	 !print *, Sr, cn,cnp,"real"
      end if 
      if (m == (cn + 110)) then
         Si = S(cnp+1)
	 !print *, Si,cn,cnp,"Imaginary"
      end if 
   enddo 
   
   
   end subroutine getS