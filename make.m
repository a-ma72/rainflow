clear rfc
clc
if 0
    for delegates = 0:1
        for int_counts = 0:1
            for val_type = {'float','double'}
                s = sprintf( ['mex ', ...
                             '-DRFC_USE_DELEGATES=%d ', ...
                             '-DRFC_USE_INTEGRAL_COUNTS=%d ', ...
                             '-DRFC_VALUE_TYPE=%s ', ...
                             '-output rfc ', ...
                             'rainflow.c'], ...
                             delegates, int_counts, val_type{1} );
                system(s, '-echo');
            end
        end
    end
else
%    mex -g -v -c COMPFLAGS='$COMPFLAGS /P' -output rfc rainflow.c
    mex -g -v -output rfc rainflow.c
end
