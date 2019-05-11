![Sapphire-Ice](https://github.com/geo-fluid-dynamics/sapphire-docs/blob/master/Sapphire-Ice-Logo.png?raw=true)

![](https://travis-ci.com/geo-fluid-dynamics/sapphire-ice.svg?branch=master)

Sapphire-Ice simulates fluid dynamics and heat transfer in phase-changing water and ice systems.

We employ the [Sapphire](https://github.com/geo-fluid-dynamics/sapphire) engine for PDE-based simulations,
which in turn employs sophisticated open-source scientific software libraries, most importantly [Firedrake](https://www.firedrakeproject.org/).


## Portfolio

### Freezing water in a cavity
As a benchmark, we simulated [Kowalewski and Rebow's water freezing experiment](https://www.researchgate.net/publication/243772766_Freezing_of_Water_in_a_Differentially_Heated_Cubic_Cavity).

![Water freezing benchmark](https://github.com/geo-fluid-dynamics/sapphire-docs/blob/master/WaterFreezing.gif?raw=true)


## Setup

First, [set up Sapphire](https://github.com/geo-fluid-dynamics/sapphire). Then...

Download with 

    git clone git@github.com:geo-fluid-dynamics/sapphire-ice.git

Activate the [Firedrake](https://www.firedrakeproject.org/) environment with something like

    . ~/firedrake/bin/activate

The following assumes that the Firedrake virtual environment is already activated.

Test with

    python3 -m pytest sapphire-ice/tests/

Install with

    cd sapphire
    
    python3 setup.py install
    
    
## Development

### Guidelines
Mostly we try to follow PEP proposed guidelines, e.g.
- [The Zen of Python (PEP 20)](https://www.python.org/dev/peps/pep-0020/) 
- [PEP 8](https://www.python.org/dev/peps/pep-0008/)

This package structure mostly follows that suggested by [The Hitchhiker's Guide to Python](http://docs.python-guide.org/en/latest/).
