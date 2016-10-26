FC = ifort
SRC = ../src2CLOWNS/
libdir = -L../lib
FFLAGS = -O3 -c -fpp -fpe0 -xHost -convert big_endian -r8 -traceback #-g -check bounds 
LFLAGS =  -O3 $(libdir) -lNRC
LINK = $(FC)
TARGET  =  ../execute/make_grid.exe
TARGET2 =  ../execute/add_teibou.exe
MODS =  0mod-type_interface.o \
	0mod-variables.o \
	0mod-arrays.o \
	data_io.o \

OBJS =  make_grid_main.o \
	1sub-main.o \
	set_grid.o \
	set_landform.o \
	cal_fluid_area.o \
	
OBJS2 = add_teibou_main.o \

all: $(TARGET) $(TARGET2)

clean: 
	-rm -rf *.o *.mod

$(TARGET): $(MODS) $(OBJS)
	$(LINK) -o $@ $(MODS) $(OBJS) $(LFLAGS) 
$(TARGET2): $(MODS) $(OBJS2)
	$(LINK) -o $@ $(MODS) $(OBJS2) $(LFLAGS)
%.o : $(srcdir)%.F90
	$(FC) $(FFLAGS) $< -o $@
