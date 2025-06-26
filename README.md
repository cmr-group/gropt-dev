# gropt-dev

Staging for the next major update to GrOpt.

## Installation
Clone the respository:

`git clone https://github.com/cmr-group/gropt-dev.git .`

Then install with (modify path to the folder you just cloned):

`pip install path/to/gropt-dev/`

or if you already have it installed:

`pip install --upgrade path/to/gropt-dev/`

## Getting Started

A simple test of operation can be performed in a python console with:
```
import gropt_dev as gropt
gropt.demo()
```
For more, see the jupyter notebooks in `./examples/`

## Critical Next Steps

The highest priority things that are being added back in are:
- [ ] More examples
- [ ] Re-implement the `gropt.gropt(params)` calling interface in the Python code.
- [ ] Remaining constraints (see other versions below)
- [ ] Provide clear feedback if the optimization was feasible or not

## Other Versions

The branch called `mess` contains the more full featured-but messier version that is being cleaned up here.  It is not being merged in with git, it is just there for reference.

It has implementation details for constraints that have not migrated to this version yet.
