#!/bin/bash

set -e

rm -f *.out
rm -rf processor*/0.*
rm -rf processor*/*e-*
rm -rf processor*/[1-9]*
rm -rf 0.* *e-* [1-9]*
rm -rf ../temp

