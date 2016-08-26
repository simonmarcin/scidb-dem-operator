# SciDB DEM operator
Runs a prototype of Firdem (by Plowman, Joseph, Charles Kankelborg, and Petrus Martens.) as SciDB operator.

## Install
* Compile as SciDB plugin. Use `make` to complie it or adapt the auto.sh script.
* Copy the `dem.so` to your SciDB plugin dir as `libdem.so`

## Usage
Create a 2D or 3D array of AIA images. Each cell should contain a pixel of seven different wavelenght images. Import the DEM operator and use it just like any other SciDB operator.

```
# load AIA images to SciDB
create array aia_256 <a0:int16,a1:int16,a2:int16,a3:int16,a4:int16,a5:int16>[x=0:4095,256,0,y=0:4095,256,0];
insert(redimension(aia_import,aia_256),aia_256);

# load the DEM operator
load_library('dem');

# run the calculation
store(dem(aia_256),'dem');
```
