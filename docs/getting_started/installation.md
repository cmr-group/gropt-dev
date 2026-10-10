# Installation

```bash
pip install --pre gropt
```

GrOpt 2.0 is in pre-release, so the `--pre` flag is needed for now.

Wheels are published for Python 3.10 to 3.14 on Linux (x86-64), macOS (Intel and Apple silicon) and Windows
(64-bit). On other platforms pip builds GrOpt from the source distribution, which needs CMake and a C++17 compiler.

The only required dependency is NumPy. The plotting helpers (`gropt.plot`) and the example notebooks also use
Matplotlib, and the notebooks need Jupyter:

```bash
pip install matplotlib jupyterlab
```

## Check the install

```python
import gropt

print(gropt.__version__)
gropt.demo()    # solves a small problem and prints whether it converged
```

Next: the [quick start](../examples/simple_example.ipynb).
