   subroutine getS(Sr,Si,cn,cnp,ik)
   use smat!, only: schan
   implicit none
   integer, intent(in) :: cn,cnp,ik
   real*8, intent(inout) :: Sr,Si
   integer m
   !real*8 S(schan+1)
   
   !open(unit=12,file="smatrix.txt")
   !rewind(12)
   !read(12,*) !E and momentum  
   !do m=1,2*schan
      !print *, m
   !   read(12,*) S(:)
      !print *, S(1)
   !   if (m == cn) then
   !      Sr = S(cnp+1)
	 !print *, Sr, cn,cnp,"real"
   !   end if 
   !   if (m == (cn + schan)) then
   !      Si = S(cnp+1)
	 !print *, Si,cn,cnp,"Imaginary"
   !   end if 
   !enddo 
   
   Sr = S(cn+1,cnp,ik)
   Si = S(cn+1,cnp+schan,ik)
   !if (cn==cnp) then
   !print *, Sr,Si,cn
   !end if
   
   end subroutine getS