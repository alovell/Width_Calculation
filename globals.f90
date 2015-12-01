   ! modules for use in widths code
   
!--------------------------------- Constants throughout the program
   module constants
   parameter(hbc=197.32705d0)
   real*8 ke,mu,const,mu1,mu2,m
   end module constants
!--------------------------------- Channel descriptors
   module channels
   implicit none
   integer channum,N
   real*8 K
   character(len=10) filename
   end module channels
!--------------------------------- Channel wavefunctions
   module gwf
   implicit none
   real*8, allocatable :: channelwf(:,:)
   real*8, allocatable :: wf(:,:)
   real*8, allocatable :: epole(:)
   integer npoles
   end module gwf
!--------------------------------- Total wavefunction
   module totalwf
   implicit none
   real*8, allocatable :: pots(:)
   real*8, allocatable :: chi(:,:)
   end module totalwf