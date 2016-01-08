   ! integral relation for calculating a width
   ! for a single channel calculation
   ! wave function is already constructed
   
   program scwidth
   implicit none
   real*8 wf(2,500)
   real*8 pot(500)
   real*8 mu,m,hbc,eng,sum,mom,Ca,gamma
   integer a1,a2
   integer i
   
   ! constants
   hbc = 197.32705
   a1 = 1
   a2 = 10
   m = 931.493
   mu = a1*a2*m/(a1+a2)
   eng = 1.39
   mom = sqrt(eng*0.0478450)
   
   ! read in the constructed wave function
   open(unit=2,file="Testwf.txt")
   do i=1,500
      read(2,*) wf(1,i),wf(2,i)
      !print *, wf(1,i), wf(2,i)
   enddo 
   
   ! construct the potential
   call scpot(wf,pot,a1,a2,hbc,mu)
   
   ! perform the integral
   ! using simpson's method
   sum = 0.d0
   call scsimpsons(pot,wf,sum,mom)
   print *, sum
   
   ! calculate Ca
   Ca = 2*mu/(mom*hbc**2)*sum
   
   ! calculate gamma/the width
   gamma = (mom*hbc**2/mu)*Ca**2
   print *, "gamma=",gamma
   
   call Gcalc(wf,hbc,mom,mu)
   
   end program