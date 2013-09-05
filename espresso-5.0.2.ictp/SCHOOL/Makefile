# Makefile for school

include ../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG)../PW/src $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../iotk/src/libiotk.a ../flib/flib.a \
         ../clib/clib.a ../flib/ptools.a 

SCHOOLOBJS = \
  stop_pp.o \
  read_export.o \
  openfile_school.o \
  check_wfcs.o \
  check_wfcs_R.o


QEMODS = ../Modules/libqemod.a
PWOBJS = ../PW/src/libpw.a


TLDEPS= bindir libs pw

all : tldeps school.x

school.x : school.o libschool.a  $(SCHOOLOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
		school.o libschool.a $(PWOBJS) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../bin ; ln -fs ../SCHOOL/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libschool.a : $(SCHOOLOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libschool.a

include make.depend
# DO NOT DELETE

