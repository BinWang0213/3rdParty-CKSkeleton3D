CK_SKEL = .
# Version LINUX
  CC = gcc
  CCFLAGS = -g -DUNIXIO
  LIBS = -lm
  ODIR = $(CK_SKEL)/obj
  BDIR = $(CK_SKEL)/bin
  IDIR = $(CK_SKEL)/include
  CDIR = $(CK_SKEL)/src/com
  LDIR = $(CK_SKEL)/src/lib
  OBJ_COMMON = $(ODIR)/mcimage.o

all:	$(BDIR)/raw2pgm $(BDIR)/pgm2raw\
$(BDIR)/skel_MK2 \
$(BDIR)/skel_AK2 \
$(BDIR)/skel_NK2 \
$(BDIR)/skel_MK3 \
$(BDIR)/skel_CK3 \
$(BDIR)/skel_EK3 \
$(BDIR)/skelpar

clean:	
	rm -f $(CK_SKEL)/bin/*
	rm -f $(CK_SKEL)/obj/*

doc:	$(CK_SKEL)/CK.dox
	doxygen $(CK_SKEL)/CK.dox

# ===============================================================
# EXECUTABLES
# ===============================================================

$(BDIR)/raw2pgm:	$(CDIR)/raw2pgm.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(OBJ_COMMON) $(ODIR)/mccodimage.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/raw2pgm.c $(OBJ_COMMON) $(ODIR)/mccodimage.o $(LIBS) -o $(BDIR)/raw2pgm

$(BDIR)/pgm2raw:	$(CDIR)/pgm2raw.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(OBJ_COMMON) $(ODIR)/mccodimage.o
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/pgm2raw.c $(OBJ_COMMON) $(ODIR)/mccodimage.o $(LIBS) -o $(BDIR)/pgm2raw

$(BDIR)/skel_AK2:	$(CDIR)/skel_AK2.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo.o $(ODIR)/lskelpar.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_AK2.c $(OBJ_COMMON) $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_AK2

$(BDIR)/skel_CK3:	$(CDIR)/skel_CK3.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar3d.h $(IDIR)/mctopo3d.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/lskelpar3d.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_CK3.c $(OBJ_COMMON) $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar3d.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_CK3

$(BDIR)/skel_EK3:	$(CDIR)/skel_EK3.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar3d.h $(IDIR)/mctopo3d.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/lskelpar3d.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_EK3.c $(OBJ_COMMON) $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar3d.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_EK3

$(BDIR)/skel_MK2:	$(CDIR)/skel_MK2.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo.o $(ODIR)/lskelpar.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_MK2.c $(OBJ_COMMON) $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_MK2

$(BDIR)/skel_MK3:	$(CDIR)/skel_MK3.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar3d.h $(IDIR)/mctopo3d.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/lskelpar3d.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_MK3.c $(OBJ_COMMON) $(ODIR)/mctopo3d.o $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar3d.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_MK3

$(BDIR)/skel_NK2:	$(CDIR)/skel_NK2.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo.o $(ODIR)/lskelpar.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skel_NK2.c $(OBJ_COMMON) $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skel_NK2

$(BDIR)/skelpar:	$(CDIR)/skelpar.c $(IDIR)/mcimage.h $(IDIR)/mccodimage.h $(IDIR)/lskelpar.h $(OBJ_COMMON) $(ODIR)/mccodimage.o $(ODIR)/mclifo.o $(ODIR)/mctopo.o $(ODIR)/lskelpar.o 
	$(CC) $(CCFLAGS) -I$(IDIR) $(CDIR)/skelpar.c $(OBJ_COMMON) $(ODIR)/mctopo.o $(ODIR)/mccodimage.o $(ODIR)/lskelpar.o $(ODIR)/mclifo.o $(LIBS) -o $(BDIR)/skelpar

# *********************************
# OBJECTS
# *********************************

$(ODIR)/lskelpar.o:	$(LDIR)/lskelpar.c $(IDIR)/mccodimage.h $(IDIR)/mctopo.h $(IDIR)/mctopo3d.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lskelpar.c -o $(ODIR)/lskelpar.o

$(ODIR)/lskelpar3d.o:	$(LDIR)/lskelpar3d.c $(IDIR)/mccodimage.h $(IDIR)/mctopo.h $(IDIR)/mctopo3d.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/lskelpar3d.c -o $(ODIR)/lskelpar3d.o

$(ODIR)/mccodimage.o:	$(LDIR)/mccodimage.c $(IDIR)/mccodimage.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mccodimage.c -o $(ODIR)/mccodimage.o

$(ODIR)/mcimage.o:	$(LDIR)/mcimage.c $(IDIR)/mccodimage.h 
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mcimage.c -o $(ODIR)/mcimage.o

$(ODIR)/mclifo.o:	$(LDIR)/mclifo.c $(IDIR)/mclifo.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mclifo.c -o $(ODIR)/mclifo.o

$(ODIR)/mctopo.o:	$(LDIR)/mctopo.c $(IDIR)/mctopo.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mctopo.c -o $(ODIR)/mctopo.o

$(ODIR)/mctopo3d.o:	$(LDIR)/mctopo3d.c $(IDIR)/mctopo3d.h
	$(CC) -c $(CCFLAGS) -I$(IDIR) $(LDIR)/mctopo3d.c -o $(ODIR)/mctopo3d.o
