# spQSP-IO
User manual for spatial quantitative systems pharmacology for immuno-oncology (spQSP-IO)

## Ubuntu Operating System Configuration (Only required for Windows User)

This step helps to setup the Ubuntu operating system in for Windows user via virtual machine.
1. The virtual machine host: VirtualBox is available at https://www.virtualbox.org/
2 . The Ubuntu Desktop image (Latest version 20.04.2) is available at   http://www.releases.ubuntu.com/20.04/
3. Enter the “Oracle VM VirtualBox Manager”, press “New” to create the virtual machine with all default settings. (Recommend allocate 20 GB for storage and 2 GB of RAM)

**Notice**: All following operations should be done in the Linux operating systems (the virtual machine), NOT Windows.

## Required Library Installation
Libraries: **SUNDIALS**: version:4.0.1; **Boost**: version 1.70.0

-SUNDIALS <br />
1. Download is available at: https://computing.llnl.gov/projects/sundials/sundials-software <br />
The following files are downloaded: 
 `sundials-4.0.1.tar.gz`

2. Decompress Archieve: <br />
`$ tar xzf sundials-4.0.1.tar.gz`

3. Install cmake if not already available: <br />
`$ sudo apt install cmake-curses-gui`

4.	Create install and build directories: <br />
`$ mkdir -p ~/lib/sundials-4.0.1`
`$ mkdir -p ~/Downloads/sundials-build`
`$ cd ~/Downloads/sundials-build`

5. Configuration <br />
`$ ccmake ~/Downloads/sundials-4.0.1` <br />
Press c key to enter configuration interface
Set install directory: CMAKE_INSTALL_PREFIX set to `~/lib/sundials-4.0.1`
Set example install directory: EXAMPLE_INSTALL_PATH set to `~/lib/sundials-4.0.1/examples`
Press c repeatedly to process configuration; press g to generate Makefile and exit configuration.
6. Build <br />
From `~/Downloads/sundials-build/` <br />
`$ make` <br />
`$ make install`

-Boost Version 1.70.0 <br />

1. Source code available at: https://www.boost.org/users/history/version_1_70_0.html <br />
The following files are downloaded: <br />
`boost_1_70_0.tar.gz`

2.	Decompress the archive: <br />
`$ tar xzf boost_1_70_0.tar.gz` <br />

Official instructions is available at:
https://www.boost.org/doc/libs/1_70_0/more/getting_started/unix-variants.html

3. Building separately-compiled boost libraries <br />
`$ cd ~/Downloads/boost_1_70_0` <br />
`$ ./bootstrap.sh --prefix=$HOME/lib/boost_1_70_0` <br />
`$ ./b2 install` <br />

## Model Simulation 
The Makefile of this model is available at: `~/TNBC_omics/TNBC_single/linux/` <br />
`$ make TNBC_s_sim` <br />
`$ ./TNBC_s_sim -h` <br />
For all options to configure the simulation
