# matlab-units

A MATLAB class that adds units and uncertainty propagation to variables.

## Overview

This repo provides a MATLAB class called `v` (for variable) that allows units and an uncertainty to be associated with a variable.

Try these commands out to start:

```matlab
v.listUnits() % list supported units
x = v(1, 'ft') % make a new v
x = x.convertTo('in') %  convert to inches
y = v(1, 'N') + 10 + x % add some items
x + y % this throws a "not dimensionally equivalent" error
y.convertToFundamentals() % show fundamental units
v(1, 'm/ft').simplifyUnits() % simplify units
```

## Installation

* Check out the repo
* use addpath to automatically include v.m
* run `help v` to learn how to use the class

## Method Reference

* listUnits() - (Static) lists the supported v units
* dimensionallyEquivalent(v, v) - (Static) returns true/1 or false/0
* v(number, 'unit string', uncertainty) - constructs a v
* convertTo('unit string') - returns a new v of equivilant magnitude
* extract('unit string') - returns a v's value after converting it to the given unit string
* convertToFundamentals() - returns a new v in fundamental units
* simplifyUnits() - returns a v with only one unit per dimension
* checkDimension('dimension character')
