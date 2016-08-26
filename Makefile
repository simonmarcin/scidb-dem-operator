ifeq ($(SCIDB),) 
  X := $(shell which scidb 2>/dev/null)
  ifneq ($(X),)
    X := $(shell dirname ${X})
    SCIDB := $(shell dirname ${X})
  endif
endif

# A way to set the 3rdparty prefix path that is convenient
# for SciDB developers.
ifeq ($(SCIDB_VER),)
  SCIDB_3RDPARTY = $(SCIDB)
else
  SCIDB_3RDPARTY = /opt/scidb/$(SCIDB_VER)
endif

# A better way to set the 3rdparty prefix path that does
# not assume an absolute path. You can still use the above
# method if you prefer.
ifeq ($(SCIDB_THIRDPARTY_PREFIX),) 
  SCIDB_THIRDPARTY_PREFIX := $(SCIDB_3RDPARTY)
endif

$(info Using SciDB path $(SCIDB))
SCIDB_VARIANT=$(shell $(SCIDB)/bin/scidb --ver | head -n 1 | cut -d : -f 2 | cut -d '.' -f 1,2 | sed -e "s/\./ *100 + /" | bc)

INSTALL_DIR = $(SCIDB)/lib/scidb/plugins

# Include the OPTIMIZED flags for non-debug use
OPTIMIZED=-O2 -DNDEBUG
DEBUG=-g -ggdb3
CFLAGS = -pedantic -W -Wextra -Wall -Wno-variadic-macros -Wno-strict-aliasing \
         -Wno-long-long -Wno-unused-parameter -fPIC -D_STDC_FORMAT_MACROS \
         -Wno-system-headers -isystem  $(OPTIMIZED) -D_STDC_LIMIT_MACROS -std=c99
CXXFLAGS = -pedantic -W -Wextra -Wall -Wno-variadic-macros -Wno-strict-aliasing \
         -Wno-long-long -Wno-unused-parameter -fPIC -DSCIDB_VARIANT=$(SCIDB_VARIANT) $(OPTIMIZED)
INC = -I. -DPROJECT_ROOT="\"$(SCIDB)\"" -I"$(SCIDB_THIRDPARTY_PREFIX)/3rdparty/boost/include/" \
      -I"$(SCIDB)/include" -I./extern

LIBS = -shared -Wl,-soname,dem.so -ldl -L. \
       -L"$(SCIDB_THIRDPARTY_PREFIX)/3rdparty/boost/lib" -L"$(SCIDB)/lib" \
       -Wl,-rpath,$(SCIDB)/lib:$(RPATH)

SRCS = LogicalDem.cpp \
       PhysicalDem.cpp \
       PreProcesing.cpp \
       Firdem2.cpp

# Compiler settings for SciDB version >= 15.7
ifneq ("$(wildcard /usr/bin/g++-4.9)","")
  CC := "/usr/bin/gcc-4.9"
  CXX := "/usr/bin/g++-4.9"
  CXXFLAGS+=-std=c++11 -DCPP11
else
  ifneq ("$(wildcard /opt/rh/devtoolset-3/root/usr/bin/gcc)","")
    CC := "/opt/rh/devtoolset-3/root/usr/bin/gcc"
    CXX := "/opt/rh/devtoolset-3/root/usr/bin/g++"
    CXXFLAGS+=-std=c++11 -DCPP11
  endif
endif


all: dem.so

clean:
	rm -rf *.so *.o

dem.so: $(SRCS)
	@if test ! -d "$(SCIDB)"; then echo  "Error. Try:\n\nmake SCIDB=<PATH TO SCIDB INSTALL PATH>"; exit 1; fi
	$(CXX) $(CXXFLAGS) $(INC) -o LogicalDem.o -c LogicalDem.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o PhysicalDem.o -c PhysicalDem.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o PreProcesing.o -c PreProcesing.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o Firdem2.o -c Firdem2.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o dem.so plugin.cpp LogicalDem.o PhysicalDem.o PreProcesing.o Firdem2.o $(LIBS)
	@echo "Now copy *.so to $(INSTALL_DIR) on all your SciDB nodes, and restart SciDB."
