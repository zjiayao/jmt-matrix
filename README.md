# jmt-matrix

*Johnson's Math ToolBox: Linear Algebra Library*

## Introduction

JMT Linear Algebra Library is a light weighted
fast library for scientific computing with C++.
Though the author does not claim the library can outperform
the celebrated LAPCK library, it is nonetheless very
easy to use and performs pretty well with matrices sized
less than 100.

This library supports common operations for real and complex
matrices including:

- Singular Value Decomposition

Single-sided Jacobi's rotation is used. This operation is
similar to a compact SVD for thin matrices.

- QR/LQ Decomposition

Householder's Reflection is used for numerical robustness.

- Cholesky Decomposition

This directly leads to the ability for matrix inversion.

## Features

- Zero Dependency

The whole library is self-contained, shipped
with a dedicated complex library used for comlex
manipulations.

- Template Programming

This library is written in template classes, that is,
there is no need to compile or link before use.
Directly include the
desired headers (more in the next section) and compile, *voila*!

This also makes the library very easy to extend to
other data types, for examples, 16bit floats, etc.

- Intuitive API

Doing a SVD and print the result is as simple
as if in `MatLAB` or `Python`:

    // a normal 5 * 3 real random matrix
    jmt::dmat mat = jmt::dmat::getNRand(5, 3), u, s, v;
    mat.SVD(u, s, v);
    u.print();

- Portability

This library provides a handful of formats for
printing, hence it is easy to feed the results
to, for example, `Python`.

    // print to file
    jmt::cmat cmat = jmt::cmat::getNRand(6, 3), q, r;
    cmat.QR(q, r);

    q.print(fout); // fout is a FILE pointer

    r.print(stdout, jmt::NUMPY);

this gives:

     [[ 0.671    -    0.741i  ,0.000                 ,0.000                    ],
      [ 0.000                 ,0.975    +    0.223i  ,0.000                    ],
      [ 0.000                 ,0.000                 ,-0.625   -    0.781i     ],
      [ 0.000                 ,0.000                 ,0.000                    ],
      [ 0.000                 ,0.000                 ,0.000                    ],
      [ 0.000                 ,0.000                 ,0.000                    ]]

Of course, directly `std::cout` and `std::cin` (yes!) are supported.


## Installation

First, clone this repo using `git`:

    git clone https://github.com/zjiayao/jmt-matrix

Then `cd` to the repo and build example:

    cd jmt-matrix && make example


## Usage

To use JMT Linear Algebra Library for operating on real numbers,
include `dmatrix.hpp` in the header; analogously, include `cmatrix.hpp`
for complex support.

A sample program for calling common APIs is given in the
`examples` folder.


## Claimer

It is recommended **NOT** to use this library in production environment,
rather, use well-crafted libraries such as LAPACK in lieu.
