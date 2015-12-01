   complex function hsh(ke,rho,K,cn)
   ! eventually need channel # for S
   ! calculat hmprime - S*hpprime
   implicit none
   real*8, intent(in) :: ke,rho,K
   integer, intent(in) :: cn
   real*8 S,rhmp,rhpp,ihmp,ihpp,L,Sr,Si
   integer ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   L=K+1.5d0
   S=0.5d0
   
   call coulfg(ke*rho,0d0,1.5d0,1.5d0,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   rhmp = gcp(2)
   ihmp = -1*fcp(2)
   rhpp = gcp(2)
   ihpp = fcp(2)
   
   ! get S matrix element (diagonal)
   call getS(K,Sr,Si,cn,cn)
   
   !hsh = rhmp - S*rhpp
   !hsh = rhmp - Sr*rhpp + Si*ihpp
   hsh = cmplx(rhmp,ihmp) - cmplx(Sr,Si)*cmplx(rhpp,ihpp)
   
   return   
   end function hsh
   
   complex function sh(ke,rho,K,cn,cnp)
   ! eventually need channel # for S
   ! calculat - S*hpprime
   implicit none
   real*8, intent(in) :: ke,rho,K
   integer, intent(in) :: cn,cnp
   real*8 S,rhmp,rhpp,ihmp,ihpp,L,Sr,Si
   integer ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   L=K+1.5d0
   S=1.d0
   
   call coulfg(ke*rho,0d0,1.5d0,1.5d0,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   rhmp = gcp(2)
   ihmp = -1*fcp(2)
   rhpp = gcp(2)
   ihpp = fcp(2)
   
   ! get S matrix element (not diagonal)
   call getS(K,Sr,Si,cn,cnp)
   
   !sh = -2*rhpp
   !sh = -Sr*rhpp + Si*ihpp
   sh = -cmplx(Sr,Si)*cmplx(rhpp,ihpp)
   
   return   
   end function sh
   
   subroutine totalwf(channelwf,npoles,energy,ke,mu,hbc,wf,epole,K,filename,cn,cnprime,N,nchan)
   implicit none
   integer, intent(in) :: npoles,cn,cnprime,N,nchan
   real*8, intent(in) :: energy,mu,hbc,ke
   real*8, intent(in) :: channelwf(N,npoles+1)
   real*8, intent(in) :: epole(npoles)
   real*8, intent(inout) :: wf(N,2)
   integer i,j,io,m,p(npoles),ival
   real*8 S,K,temppole,tempwf(N,2),const
   complex*8 hsh,sh,compwf(N)
   character(len=10) files(nchan),file,filename
   complex*8, allocatable :: Ap(:)
    
   ! calculate total wave function
   ! Nuclear Reactions for Astrophysics IJ Thompson, FM Nunes
   ! Eqn's (6.5.34) and (6.5.35)
   open(unit=7,file="totalwf.txt",access="append")
   write(7,*) "#",filename,K
   print *, filename
   
   open(unit=10,file="wfpoleref.txt",status="replace")
   ! calculate Ap
   p=0
   ! first find the number of files that have the pole
   do i=1,npoles
   !print *, i, epole(i)
   do j=1,nchan
      write (files(j), '( I0, "K.wf" )') j
      open(unit=9,file=files(j))
      !print *, files(j)
      do 
	 read(9,*,iostat=io) temppole
	 if (io/=0) exit
	 if (temppole==epole(i)) then
	    !print *, temppole
	    p(i) = p(i) + 1
	    ! print the file and pole number for easy access
	    write(10,*) i, files(j)
	 end if 
      enddo
      rewind(9)
      close(9)
   enddo 
   enddo
   rewind(10)
   
   allocate (Ap(npoles))
   Ap=cmplx(0.d0,0.d0)
   do 
      read(10,*,iostat=io) ival,file
      if (io/=0) exit
      open(unit=11,file=file)
      do 
	 read(11,*) temppole
	 if (temppole==epole(ival)) then
	    !print *, "ival=",ival,"pole=",temppole
	    do j=1,N
	       read(11,*) tempwf(j,1), tempwf(j,2)
	      ! print *, tempwf(j,1), tempwf(j,2)
	    enddo
	    !print *, temppole, tempwf(N,1),tempwf(N,2)
	    exit
	 end if 
      enddo 
      close(11)
      const = (hbc**2/(2.d0*mu))*(1.d0/(epole(ival)-energy))
      if (file==filename) then
         Ap(ival) = Ap(ival) + const*tempwf(N,2)*hsh(ke,tempwf(N,1),K,cn)
	 !print *, hsh(ke,tempwf(N,1),K,cn)
      else 
         Ap(ival) = Ap(ival) + const*tempwf(N,2)*sh(ke,tempwf(N,1),K,cn,cnprime)
	 !print *, sh(ke,tempwf(N,1),K,cn,cnprime)
      end if 
      !print *, Ap(ival), tempwf(N,2), epole(ival)
   enddo 
   close(10)
   
   ! calculate total wave function
   wf=0.d0
   do i=1,N
      wf(i,1) = channelwf(i,1)
      do j=1,npoles
         !wf(i,2) = wf(i,2) + channelwf(i,j+1)*Ap(j)
	 compwf(i) = compwf(i) + channelwf(i,j+1)*Ap(j)
      enddo 
      wf(i,2) = abs(compwf(i))
      write(7,*) wf(i,1), wf(i,2)
   enddo 
   write(7,*) "&"
   close(7)
   
   end subroutine totalwf
