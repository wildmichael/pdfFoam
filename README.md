# Transported JPDF Library and Solver for Reactive Flow Simulation

pdfFoam is a general purpose JPDF algorithm for the simulation of turbulent
reactive flows in OpenFOAM. The algorithms, its implementation and the
validation test cases are described in detail in the authors PhD thesis
[[Wild2013]](doc/references.md#Wild2013).

## Installation

### Documentation
Some documentation is directly available as Markdown:
[doc/README.md](doc/README.md). The code, however, is documented extensively
with Doxygen. Refer to below instructions for details on how to build it.

### Requirements

In the following the prerequisite software packages are listed:

- OpenFOAM: Versions 1.7.x or 2.1.x have been tested. OpenFOAM-extend has not
  been tried at all. Installation instructions are available from the official
  OpenFOAM [download page](http://openfoam.org/download).
- A operating system that is compatible with OpenFOAM. These include popular
  Linux distributions, such as Ubuntu, OpenSUSE and Fedora. Mac OS X works
  with a lot of hackery. Simply forget about MS Windows. Refer to the OpenFOAM
  installation instructions on version requirements. The requirements can
  usually be relaxed when compiling OpenFOAM from source, but that is not for
  the inexperienced.
- A C++ compiler that is compatible with OpenFOAM. For OpenFOAM 1.7.x these
  are g++ versions 4.2 to 4.6. For OpenFOAM 2.1.x this should be g++ 4.3
  through 4.7 and fairly recent versions of Clang.
- GNU Make

On Debian based distributions, the following commands will install OpenFOAM
1.7.1 (the version best tested) and all other dependencies:

```sh
sudo sh -c "echo deb http://www.openfoam.com/download/ubuntu maverick main > /etc/apt/sources.list.d/openfoam171.list"
sudo apt-get update
sudo apt-get install openfoam171 paraviewopenfoam381
```

In case you are using Ubuntu 10.04 LTS, replace `maverick` with `lucid` in
above commands.

If you want to genereate the Doxygen documentation, you also will need:

- Doxygen
- GNU Sed
- GNU Awk
- Graphviz
- LaTeX: Depending on the version of Doxygen, you need a LaTeX distribution
  with `dvipng` installed. This is only necessary if your Doxygen doesn't
  support MathJax and loading of MathJax-extensions for formula rendering.
  Doxygen version 1.7.5 and newer should do.

On Debian based distributions, the following will install these packages
(excluding LaTeX):

```sh
sudo apt-get install doxygen sed gawk graphviz
```

The tutorial cases use [Python Matplotlib](http://matplotlib.org) to create the
graphs. Version 1.2.0 or newer is required. Unfortunately, e.g. the current
Ubuntu 12.10 only provides version 1.1.1. To install a newer version without
affecting the one installed by the system, `virtualenv` can be used to create
an isolated Python installation where the package can be updated:

```sh
sudo apt-get install python-virtualenv
virtualenv --system-site-packages $HOME/Software/python
source $HOME/Software/python/bin/activate
echo "source $HOME/Software/python/bin/activate" >> ~/.bash_profile
pip install --upgrade matplotlib
```

### Building
In order to compile the library `mcParticle` and the `pdfSimpleFoam`
solver, in the top-level source directory issue the following command:

```sh
./Allwmake
```

To build the documentation, use instead

```sh
./Allwmake doc
```

## Acknowledgements
The math renderings in the Markdown files are generously provided by CodeCogs:

[![CodeCogs - An Open Source Scientific Library](http://www.codecogs.com/images/logo.gif)](http://www.codecogs.com)
