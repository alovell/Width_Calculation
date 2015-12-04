   ! eventually need channel # for S
   ! calculat hmprime - S*hpprime
   complex function hsh(rho,cn,ik)
   use constants
   use channels
   use smat
   implicit none
   !real*8, intent(in) :: ke,rho,K
   real*8, intent(in) :: rho
   integer, intent(in) :: cn,ik
   real*8 rhmp,rhpp,ihmp,ihpp,L,Sr,Si
   integer ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   L=K+1.5d0
   !S=0.5d0
   
   call coulfg(en_ke(2,ik)*rho,0d0,1.5d0,1.5d0,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   rhmp = gcp(2)
   ihmp = -1*fcp(2)
   rhpp = gcp(2)
   ihpp = fcp(2)
   
   ! get S matrix element (diagonal)
   call getS(Sr,Si,cn,cn,ik)
   
   !hsh = rhmp - S*rhpp
   !hsh = rhmp - Sr*rhpp + Si*ihpp
   hsh = cmplx(rhmp,ihmp) - cmplx(Sr,Si)*cmplx(rhpp,ihpp)
   
   return   
   end function hsh

   ! eventually need channel # for S
   ! calculat - S*hpprime
   complex function sh(rho,cn,cnp,ik)
   use constants
   use channels
   use smat
   implicit none
   real*8, intent(in) :: rho
   integer, intent(in) :: cn,cnp,ik
   real*8 rhmp,rhpp,ihmp,ihpp,L,Sr,Si
   integer ifail,m1
   double precision, dimension(101) :: fc,gc,fcp,gcp
   ifail=0
   m1=1
   L=K+1.5d0
   !S=1.d0
   
   call coulfg(en_ke(2,ik)*rho,0d0,1.5d0,1.5d0,fc,gc,fcp,gcp,1,0,ifail,m1)
   
   rhmp = gcp(2)
   ihmp = -1*fcp(2)
   rhpp = gcp(2)
   ihpp = fcp(2)
   
   ! get S matrix element (not diagonal)
   call getS(Sr,Si,cn,cnp,ik)
   
   !sh = -2*rhpp
   !sh = -Sr*rhpp + Si*ihpp
   sh = -cmplx(Sr,Si)*cmplx(rhpp,ihpp)
   
   return   
   end function sh

   subroutine totalwf(ifile,cn,cnprime,nchan,ik)
   use constants
   use channels
   use gwf
   use smat
   implicit none
   integer, intent(in) :: cn,cnprime,nchan,ik
   integer i,j,io,p(npoles),ival
   real*8 temppole,tempwf(N,2),con
   complex*8 hsh,sh,compwf(N)
   character(len=10) files(nchan),file,ifile !filename
   complex*8, allocatable :: Ap(:)
   
   !print *, hbc,m,energy,ke,mu
   !print *, hbc,m,en_ke(1,1),en_ke(2,1),mu
    
   ! calculate total wave function
   ! Nuclear Reactions for Astrophysics IJ Thompson, FM Nunes
   ! Eqn's (6.5.34) and (6.5.35)
   open(unit=7,file="totalwf.txt",access="append")
   write(7,*) "#",ifile,K
   print *, ifile
   
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

   !print *, en_ke(1,ik)
   
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
      !con = (hbc**2/(2.d0*mu))*(1.d0/(epole(ival)-energy))   
   print *, "here",en_ke(1,2)
      con = (hbc**2/(2.d0*mu))*(1.d0/(epole(ival)-en_ke(1,ik)))
      !print *, con
      !print *, "here"
      if (file==ifile) then
	 Ap(ival) = Ap(ival) + con*tempwf(N,2)*hsh(tempwf(N,1),cn,ik)
	 !print *, hsh(tempwf(N,1),cn,ik)
      else 
	 Ap(ival) = Ap(ival) + con*tempwf(N,2)*sh(tempwf(N,1),cn,cnprime,ik)
	 !print *, sh(tempwf(N,1),cn,cnprime,ik)
      end if 
      !print *, Ap(ival), tempwf(N,2), epole(ival)
   enddo 
   close(10)

   
   ! calculate total wave function
   wf=0.d0
   compwf=cmplx(0.d0,0.d0)
   do i=1,N
      wf(i,1) = channelwf(i,1)
      !print *, channelwf(i,1)
      do j=1,npoles
         !wf(i,2) = wf(i,2) + channelwf(i,j+1)*Ap(j)
	 compwf(i) = compwf(i) + channelwf(i,j+1)*Ap(j)
	 !print  *, Ap(j),channelwf(i,j+1),compwf(i)
	 !print *, Ap(j),channelwf(i,j+1),channelwf(i,j+1)*Ap(j)
	 !print *, channelwf(i,j+1)*Ap(j)
	 !print *, compwf(i)
      enddo 
      wf(i,2) = abs(compwf(i))
      write(7,*) wf(i,1), wf(i,2)
   enddo 
   write(7,*) "&"
   close(7)
   
   deallocate(Ap)
   
   end subroutine totalwf
