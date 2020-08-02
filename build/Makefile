FCC = gfortran
libroot = /Users/dangch/Documents/CRTM/pyCRTM/pycrtm/CRTM/crtm
FCFLAGS = -I${libroot}/include
LDFLAGS = -L${libroot}/lib
all:
	${FCC} -O0 -g  ../src/main.f90 ${FCFLAGS} ${LDFLAGS} -lcrtm -o InspCloudCoeff.x
	${FCC} -O0 -g  ../src/main_Aerosol.f90 ${FCFLAGS} ${LDFLAGS} -lcrtm -o InspAerosolCoeff.x
