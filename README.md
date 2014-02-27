# matlab-units

A MATLAB class that adds units and uncertainty propagation to variables.

## Overview

This repo provides a MATLAB class called `v` (for variable) that allows units and an uncertainty to be associated with a variable.

Try these commands out to start:

```matlab
% list supported units
v.listUnits()

% make a new v
x = v(1, 'ft')

% convert to inches
x = x.convertTo('in')

% add with scalars and other v's
x = x + 10  + v(1, 'cm')

% try something silly (throw a "not dimensionally equivalent" error)
x + v(1, 'N')

% break down into fundamental units
v(1, 'N').convertToFundamentals()

% reduce dimensionally equivalent units to one unit type
v(1, 'm*ft').simplifyUnits()
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
