x            = cumsum( randn( 1e6, 1 ) );
x_max        = max(x);
x_min        = min(x);
class_count  = 100;
class_width  = (x_max - x_min) / (class_count - 1);
class_offset = x_min - class_width / 2;
hysteresis   = class_width;

[pd,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );

close all
plot( re );

figure
surface( rm );

pd
