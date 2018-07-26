***********************************************************************
*** Program nusigma calculates the neutrino-nucleon cross sections
*** Date: 2005-10-21
*** Joakim Edsjo, edsjo@physto.se
*** Modified by C. Arguelles to calculate 1 point.
*** ModDate : 2011-06-11
***********************************************************************

      function dsdxdy(EE,x,y,neutype,targettype,intertype)
      implicit none

      include 'nupar.h'

      real*8 dsdxdy
      real*8 Nudsdxdy
      real*8 EE
      real*8 x
      real*8 y
      integer neutype
      integer targettype
      integer intertype

      character targettype_char
      character*2 intertype_char

      if (targettype == 1) then
        targettype_char = 'p'
      else if (targettype == 0) then
        targettype_char = 'n'
      else if (targettype == 2) then
        targettype_char = 'N'
      end if

      if (intertype == 1) then
        intertype_char = 'CC'
      else
        intertype_char = 'NC'
      end if

      call nusetup

      dsdxdy = Nudsdxdy(EE,x,y,neutype,targettype_char,intertype_char)

      end

      function dsde(E1,E2,neutype,targettype,intertype)
      implicit none

      include 'nupar.h'

      real*8 dsde
      real*8 NuCrossDiffl
      real*8 E1,E2
      integer neutype
      integer targettype
      integer intertype

      character targettype_char
      character*2 intertype_char

      if (targettype == 1) then
        targettype_char = 'p'
      else if (targettype == 0) then
        targettype_char = 'n'
      else if (targettype == 2) then
        targettype_char = 'N'
      end if

      if (intertype == 1) then
        intertype_char = 'CC'
      else
        intertype_char = 'NC'
      end if

      call nusetup

      dsde = Nucrossdiffl(E1,E2,neutype,targettype_char,intertype_char)

      end

      function sigma(EE,neutype,targettype,intertype)
      implicit none

      include 'nupar.h'

      real*8 sigma
      real*8 NuCross
      real*8 EE
      integer neutype
      integer targettype
      integer intertype

      character targettype_char
      character*2 intertype_char

      if (targettype == 1) then
        targettype_char = 'p'
      else if (targettype == 0) then
        targettype_char = 'n'
      else if (targettype == 2) then
        targettype_char = 'N'
      end if

      if (intertype == 1) then
        intertype_char = 'CC'
      else
        intertype_char = 'NC'
      end if

      call nusetup

      sigma = NuCross(EE,neutype,targettype_char,intertype_char,2)

      end


