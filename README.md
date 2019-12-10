# Local build for in-situ post-processing in Catalyst and Nek5000
These instructures were used to compile Catalyst in Nek5000, for Ubuntu 18.04 (09/12/2019).

## Prerequisite

Install packages: 

*sudo apt install build-essential cmake-curses-gui llvm mpich libboost-all-dev*

and Python 3.7.

The source codes of mesa-18.3.3 and ParaView-v5.6.3 are located in *~/InSituPackage*. Binaries file will be placed in *~/InSituPackage/local*.

## Mesa (version 18.3.3):

1) Run in *~/InSituPackage/mesa-18.3.3* the following script:

```bash
export PYTHON=/usr/bin/python3

./configure                                          \
      --enable-opengl --disable-gles1 --disable-gles2   \
      --disable-va --disable-xvmc --disable-vdpau       \
      --enable-shared-glapi                             \
      --with-gallium-drivers=swrast                     \
      --disable-dri --with-dri-drivers=                 \
      --disable-egl --with-egl-platforms= --disable-gbm \
      --disable-glx                                     \
      --disable-osmesa --enable-gallium-osmesa          \
      --enable-shared \
      --enable-texture-float \
      --prefix=/home/marco/InSituPackage/local/
```

**Note:** prefix is the path were the bin will be located.

2) In *~/InSituPackage/mesa-18.3.3*? (***TODO: check if this path is correct***)

*make -j*

*make install*

3) export mesa variable, using the script:

```bash
export OSMESA=/home/marco/InSituPackage/local

export LIBDIR=$OSMESA/lib:$LIBDIR
export LD_LIBRARY_PATH=$OSMESA/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$OSMESA/lib:$LD_RUN_PATH

export OSMESA_INCLUDE_DIR=$OSMESA/include
export OSMESA_LIBRARY=$OSMESA/lib
```

## Python (version 5.6.3):

1) Run the following script:

```bash
cmake \
      -DCMAKE_BUILD_TYPE=Release                                  \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON                     \
      -DPARAVIEW_ENABLE_PYTHON=ON                                 \
      -DVTK_PYTHON_VERSION=3                                      \
      -DPARAVIEW_BUILD_QT_GUI=OFF                                 \
      -DPARAVIEW_USE_MPI=ON                                       \
      -DPARAVIEW_ENABLE_CATALYST=ON                               \
      -DVTK_USE_X=OFF                                             \
      -DVTK_OPENGL_HAS_OSMESA=ON                                  \
      -DVTK_USE_OFFSCREEN=OFF                                     \
      -DOSMESA_INCLUDE_DIR=$OSMESA_INCLUDE_DIR                     \
      -DOSMESA_LIBRARY=$OSMESA_LIBRARY/libOSMesa.so                \
      -DPYTHON_INCLUDE_DIR=/usr/include/python3.6m \
      -DPYTHON_LIBRARY=/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu/libpython3.6m.so   \
      -DCMAKE_INSTALL_PREFIX=/home/marco/InSituPackage/local      \
      /home/marco/InSituPackage/ParaView-v5.6.3
```

**Note:** the last two lines specify the install and the source codes locations, respectively. The forth and third-last lines need to be set according with the system (where Python is installed).

**Note:** in this example, both Mesa and Paraview share the same the path for binaries (*/home/marco/InSituPackage/local*).

2) In *~/InSituPackage/ParaView-v5.6.3*? (***TODO: check if this path is correct***)

*make -j*

*make install*

## Nek5000:

1) Before compiling a case with the InSitu implementation, export Python-path enrivoment variables with the script:

```bash

export PARAVIEW=/home/marco/InSituPackage/local
export PATH=$PARAVIEW/bin:$PATH
export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

export PYTHONPATH=$PARAVIEW/lib/python3.6/site-packages:$PYTHONPATH
export PYTHONPATH=$PARAVIEW/lib/python3.6/site-packages/paraview/:$PYTHONPATH

```

2) In *makenek*, add paraview as optional compiler flag for the C compiler only:

```bash
FFLAGS="-I./inc_src -g"
CFLAGS="-I./inc_src -I$PARAVIEW/include/paraview-5.6"
```

and add Catalyst enrivoment variables:

```bash
CATALYST=1
CATALYST_LIBS=`paraview-config --libs vtkPVPythonCatalyst`
CATALYST_INCS=`paraview-config --include vtkPVPythonCatalyst`
```



