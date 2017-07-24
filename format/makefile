#----------------------------------------------------------------------
#
#     FILE: Makefile
#
#     DESCRIPTION:
#     
#     A makefile for the gausmake routine that converts
#     data of a netcdf file into a monthfile of ausmap type
#
#     HISTORY:
#
#     WRITTEN:: A. Pepler
#     LAST CHANGE: 23/08/2013
#
#----------------------------------------------------------------------

FFLAGS = -C -frecord-marker=4 #-fbounds-check #-finit-local-zero #-fno-range-check
F90C = gfortran

# Libraries
LIBS= -lnetcdff -lnetcdf -lgfortran
LIBDIRGFORT := $(shell echo $(LIBDIRGFORT))
LIBDIRNETCDF := $(shell echo $(LIBDIRNETCDF))


all: monthly yearly topo getmonthyear

monthly:
	$(F90C) $(FFLAGS) gausmake_monthly.f $(LIBDIGFORT) $(LIBDIRNETCDF) $(LIBS) -o gausmake_monthly

yearly:
	$(F90C) $(FFLAGS) gausmake_yearly.f $(LIBDIGFORT) $(LIBDIRNETCDF) $(LIBS) -o gausmake_yearly

topo:
	$(F90C) $(FFLAGS) gausmake_topo.f $(LIBDIGFORT) $(LIBDIRNETCDF) $(LIBS) -o gausmake_topo

getmonthyear:
	$(F90C) $(FFLAGS) getmonthyear.F $(LIBDIGFORT) $(LIBDIRNETCDF) $(LIBS) -o getmonthyear
	mv getmonthyear ../run/getmonthyear

clean:
	rm -f gausmake_monthly
	rm -f gausmake_yearly
	rm -f gausmake_topo
