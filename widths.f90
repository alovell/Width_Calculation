
! calculate the ANC (and widths) from K. Nollett PRC 86, 044330 (2012), eqn. 35

   program widths
   use constants
   use channels
   use gwf
   use totwf
   implicit none
   real*8 partsum,Ca,gamma,temp,rad(3),mass1,mass2
   integer io,nlines,i,j,l,narg,Knum,nchan,npot,schan
   character (len=10) temp2,file(99)
   character (len=20) sfile
   
   namelist/nuclei/mass1,mass2
   namelist/files/sfile,filename
   namelist/chaninfo/nchan,npot,schan
   
   ! default values for namelist objects
   sfile = "smatrix.txt"
   filename = "1K.wf"
   mass1 = 1.d0
   mass2 = 14.d0
   nchan = 1
   npot = 400
   schan = 110
   
   ! read namelists
   !read(5,nml=nuclei)
   !print *, mass1,mass2
   !read(5,nml=files)
   !print *, sfile,filename
   !read(5,nml=chaninfo)
   !print *, nchan,npot,schan
   
   !filename = Kfile
   
   ! read file name from input line (#K.wf)
   narg = IARGC()
   call getarg(1,filename)
   print*, filename
   ! get the number of the channel
   Knum = index(filename,'K')
   read(filename(:Knum-1),*) channum
   print *, Knum, channum
   
   ! constants
   hbc = 197.32705d0
   m = 931.49432d0
   ! this is specifically for n+n+14Be
   mu = 2.d0*14.d0*m/16.d0
   print *, hbc,m,mu
   
   ! read in energy (scattering energy)
   open(unit=7,file="smatrix.txt")
   read(7,*) energy, ke
   rewind(7)
   close(7)
   print *, energy, ke
   
   ! read input file (filename) to get number of radial points
   open(unit=8,file=filename)
   do
      rad(3)=rad(2)
      rad(2)=rad(1)
      read(8,*,iostat=io) rad(1)
      if (io/=0) exit
   enddo 
   print *, rad(1),rad(2),rad(3)
   N = int(rad(1)/(rad(1)-rad(3)))
   print *, N
   rewind(8)
   close(8)
   
   ! read in wave function file
   ! need to figure out how to read this in
   nchan=5 ! number of channels in FaCE - small while testing
   
   ! allocate wf (channel wf) and chi (total wf)
   allocate(wf(N,2))
   allocate(chi(N,nchan+1))

   ! essentially refresh the total wave function file
   ! since it will append in totalwf subroutine
   open(unit=7,file="totalwf.txt")
   write(7,*) "#wave functions for channel gamma"
   close(7)
   
   ! construct all total wave functions
   do i=1,nchan
      write (file(i), '( I0, "K.wf" )') i
      open(unit=3,file=file(i))
      ! get number of poles
      nlines=0
      do 
         read(3,*,iostat=io) temp
	 if (io/=0) exit
	 nlines = nlines + 1
      enddo 
      npoles = nlines/(N+1)
      rewind(3)
      ! allocate pole energies and g_\alpha^p's
      allocate(epole(npoles))
      allocate(channelwf(N,npoles+1))
      do j=1,npoles
         !print *, j
         read(3,*) epole(j), temp, temp2, K
	 !print *, epole(j), K
	 do l=1,N
	    read(3,*) channelwf(l,1), channelwf(l,j+1)
	    !print *, l
	 enddo 
      enddo
      ! construct the total wave function \chi for each channel, i
      !call totalwf(channelwf,npoles,energy,ke,mu,hbc,wf,epole,K,file(i),channum,i,N,nchan)
      call totalwf(file(i),channum,i,nchan)
      ! put each \chi into matrix holding all wave functions
      do j=1,N
         chi(j,1) = wf(j,1)
	 chi(j,i+1) = wf(j,2)
      enddo    
      deallocate(epole)
      deallocate(channelwf)
   enddo 
   
   ! deallocate wf since \chi is constructed
   !deallocate (wf)
   
   ! create the potentials
   ! pseudo potentials
   ! pot = \sum _gamma prime V_{gamma gamma prime} chi_gamma prime
   allocate(pots(400))
   call potsum(nchan)
      
   ! deallocate channelwf
   !deallocate (channelwf)
   
   ! constant for Ca
   mu1 = 0.5
   mu2 = mu/m
   const = 2*mu/(hbc**2*ke*(mu1*mu2)**(1.5))
   
   ! perform the intergration
   ! using Simpson's Rule 
   partsum = 0
   call simpsons(partsum,nchan)
   
   ! deallocate chi
   deallocate (chi)
   
   Ca = const * partsum
   
   ! width
   gamma = ke*(Ca*hbc)**2/mu
   print *, "gamma=",gamma

   end program widths