![Build](https://github.com/SimiPro/DCA/actions/workflows/build.yml/badge.svg?branch=master)

# DCA - Differentiable Collision Avoidance

We implement differentiable collision avoidance using multiple shape primitives.

## Primitives
The library supports the following primitives with the mentioned states and parameters:
- Sphere, where the state is the center of the sphere.
- Capsule, where the state is the two ends of the capsule, stacked.
- Rectangle, where the state is the center of mass and the orientation (expressed in exponential coordinates).
- Box, where the state is the center of mass and the orientation (expressed in exponential coordinates).

Furthermore, the size of the parameters is:
- 0 for a sphere
- 1 for a capsule
- 2 for a rectangle
- 3 for a box

## Example
We refer the reader to the [example](example/example.cpp).

## API
The API can be found inside the [API header](include/DCA/API.h). Generally speaking, one passes two primitives and gets back either the distance or the first or second derivative. If the states of the two primitives do not change, one can also pass the computed parameterization (`t`) using the specialized functions. The parameterization can also be computed with a single API call and the two given primitives.