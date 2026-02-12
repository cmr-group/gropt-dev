# gropt-dev

[![Windows](https://github.com/cmr-group/gropt-dev/actions/workflows/build-windows.yml/badge.svg)](https://github.com/cmr-group/gropt-dev/actions/workflows/build-windows.yml)
[![Linux](https://github.com/cmr-group/gropt-dev/actions/workflows/build-linux.yml/badge.svg)](https://github.com/cmr-group/gropt-dev/actions/workflows/build-linux.yml)
[![macOS](https://github.com/cmr-group/gropt-dev/actions/workflows/build-mac.yml/badge.svg)](https://github.com/cmr-group/gropt-dev/actions/workflows/build-mac.yml)

Staging for the next major update to GrOpt.

## Installation
Clone the respository:

`git clone https://github.com/cmr-group/gropt-dev.git`

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
