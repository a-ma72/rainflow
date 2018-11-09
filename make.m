clear rfc
clc
mex -g -v -output rfc -DRFC_USE_INTEGRAL_COUNTS=0 rainflow.c
