#!/bin/bash

set -e

rm -f *RT.out
rm -rf processor*/0.*
rm -rf processor*/[1-9]*
rm -rf processor*/*e-*

rm -rf 0.* *e-* [1-9]*

rm -f *.csv
