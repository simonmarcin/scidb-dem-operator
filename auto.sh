
export SCIDB_VER=15.12
export SCIDB_INSTALL_PATH=/home/scidb/src/scidb-15.12.1.4cadab5/stage/install
export SCIDB_BUILD_TYPE=Debug
export PATH=$SCIDB_INSTALL_PATH/bin:$PATH
 
cd /home/scidb/Desktop/scidb-15.12.1.4cadab5/examples/DEM/
make
cp dem.so /home/scidb/src/scidb-15.12.1.4cadab5/stage/install/lib/scidb/plugins/libdem.so
scidb.py stopall mydb
scidb.py startall mydb
