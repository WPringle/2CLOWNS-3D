FC = ifort
SRC = 
incdir = -I../include
libdir = -L../lib
FFLAGS = -O3 -c -openmp -fpp -fpe0 -xHost -convert big_endian -r8 $(incdir) -DNORMAL -DOBJECT -DRAN -DSLIDE0 -traceback #-g -check bounds 
LFLAGS =  -O3 -openmp $(libdir) -lNRC
LINK = $(FC)
TARGET =  ../execute/2CLOWNS_3D_Oct6.exe
TARGET2 =  
OBJS =  0mod-type_interface.o \
	0mod-variables.o \
	0mod-arrays.o \
	1sub-2CLOWNS-coupling.o \
	1sub-2CLOWNS-fileout.o \
	1sub-2CLOWNS-local.o \
	1sub-2DH-data_io.o \
	1sub-2DH-initial_reset.o \
	1sub-2DH-momentum.o \
	1sub-2DH-continuity.o \
	1sub-drift.o \
	1sub-drift-Pressure.o \
	1sub-main.o \
	1sub-matrix.o \
	1sub-momentum.o \
	1sub-normal.o \
	1sub-normal-ini.o \
	1sub-normal-PlainAdjust.o \
	1sub-normal-set_surface_normal2.o \
	1sub-normal-set_surface_position_and_info.o \
	1sub-normal-suimen_teimen.o \
	1sub-normal-suimen_teimen_adjust.o \
	1sub-object.o \
	1sub-turbulence.o \
	1sub-turbulence-cal_ke-ex.o \
	1sub-vof.o \
	1sub-wave-2DHvar.o \
	1sub-wave-BCIC.o \
	pcgpme_2DH.o \
	cal_advection_viscosity.o \
	cal_courant.o \
	data_io.o \
	nfd2015.o \
	pcgpme_para.o \
	SegmentShape2014.o \
	SegmentShape2014_mainloop.o \
	set_surface_velocity.o \
	main_2CLOWNS.o \
	
OBJS2 = sample.o \

all: $(TARGET) $(TARGET2)

clean: 
	-rm -rf *.o *.mod

$(TARGET): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LFLAGS) 
#$(TARGET2): $(OBJS2)
#	$(LINK) -o $@ $(OBJS2) -openmp
%.o : $(srcdir)%.F90
	$(FC) $(FFLAGS) $< -o $@
%.o : $(srcdir)%.f90
	$(FC) $(FFLAGS) $< -o $@

#ifort src/pcgpme.F -c -fast -fpp -openmp -parallel -fpe0 -convert big_endian 
#@printf "$(OBJS)\n"
