FCC = ifort
libroot = /data/users/cdang/CRTM/CRTM_dev_CD/CRTM_dev/src/Build/crtm_v2.3.0
NC4_DIR=/opt/netcdf4/4.6.2-intel-18.0.3
HDF_DIR=/opt/hdf5/1.8.21-intel-18.0.3
HDF4_DIR=/opt/hdf4/4.2.14-intel-18.0.3

FCFLAGS = -I${libroot}/include \
          -I$(NC4_DIR)/include \
          -I$(HDF_DIR)/include 
LDFLAGS = -L${libroot}/lib \
          -L$(NC4_DIR)/lib -lnetcdf -lnetcdff \
          -L$(HDF_DIR)/lib -lhdf5
all:
	${FCC} -O0 -g  ../src/main.f90 ${FCFLAGS} ${LDFLAGS} -lcrtm -o InspCloudCoeff.x
	${FCC} -O0 -g  ../src/main_aerosol.f90 ${FCFLAGS} ${LDFLAGS} -lcrtm -o InspAerosolCoeff.x
