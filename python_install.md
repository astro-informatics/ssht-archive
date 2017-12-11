# Installing the python wrappers

First open setup.py and in the relavant part point it to your fftw3 installation

then run:
```
python setup.py build_ext --inplace
```

For user-wide install of the package use the following:
```
python setup.py install --user
```
