# defining bodies

cdef class Track:
    cdef Track(self,double,double)
    #cdef Track(self)
    cdef double x
    cdef double xini
    cdef double xend

cdef class Body:
    cdef char* name    
    cdef double density(self,Track)
    cdef double ye(self,Track)