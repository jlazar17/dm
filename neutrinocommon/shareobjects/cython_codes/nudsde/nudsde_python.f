***********************************************************************
*** Program nusigma calculates the neutrino-nucleon cross sections
*** Date: 2005-10-21
*** Joakim Edsjo, edsjo@physto.se
*** Modified by C. Arguelles to calculate 1 point.
*** ModDate : 2011-06-11
***********************************************************************
      function dsde(E1,E2,neutype,targettype,intertype)
      implicit none

      include 'nupar.h'
      
      real*8 dsde
      real*8 NuCrossDiffl
      real*8 E1,E2
      integer neutype
      character targettype
      character*2 intertype

      call nusetup

      dsde = Nucrossdiffl(E1,E2,neutype,targettype,intertype)

      end
      
      function sigma(EE,neutype,targettype,intertype)
      implicit none

      include 'nupar.h'
      
      real*8 sigma
      real*8 NuCross
      real*8 EE
      integer neutype
      character targettype
      character*2 intertype

      call nusetup

      sigma = NuCross(EE,neutype,targettype,intertype,2)

      end


