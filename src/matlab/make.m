clear rfc
clc
mex -g -v -output rfc -DRFC_USE_INTEGRAL_COUNTS=0 -DRFC_ASTM_SUPPORT=1 -DRFC_TP_SUPPORT=1 -I../../lib ../rfc.c ../../lib/rainflow.c
