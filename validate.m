%% Empty series
name         = 'empty';
x            = export_series( name, [] );
x_max        = 1;
x_min        = -1;
class_count  = 100;
class_width  = round( (x_max - x_min) / (class_count - 1), 2 );
class_offset = x_min - class_width / 2;
hysteresis   = class_width;

[~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );

assert( sum( sum( rm ) ) == 0 );

assert( isempty(re) );

save( name, 'rm', 're' );


%% One single cycle (up)
name         = 'one_cycle_up';
x            = export_series( name, [1,3,2,4] );
x_max        = 4;
x_min        = 1;
class_count  = 4;
class_width  = round( (x_max - x_min) / (class_count - 1), 2 );
class_offset = x_min - class_width / 2;
hysteresis   = class_width * 0.99;

[~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );
rm = rm'; % Transpose matrix ("from" in rows, "to" in columns)

assert( sum( sum( rm ) ) == 1 );
assert( rm( 3,2 ) == 1 );

assert( isequal( re, [1;4] ) );

save( name, 'rm', 're' );


%% One single cycle (down)
name         = 'one_cycle_down';
x            = export_series( name, [4,2,3,1] );
x_max        = 4;
x_min        = 1;
class_count  = 4;
class_width  = round( (x_max - x_min) / (class_count - 1), 2 );
class_offset = x_min - class_width / 2;
hysteresis   = class_width * 0.99;

[~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );
rm = rm'; % Transpose matrix ("from" in rows, "to" in columns)

assert( sum( sum( rm ) ) == 1 );
assert( rm( 2,3 ) == 1 );

assert( isequal( re, [4;1] ) );

save( name, 'rm', 're' );


%% Small example (Siemens)
% [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
name         = 'small_example';
x            = export_series( name, [2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2] );
x_max        = max(x);
x_min        = min(x);
class_count  = x_max;
class_width  = round( (x_max - x_min) / (class_count - 1), 2 );
class_offset = x_min - class_width / 2;
hysteresis   = class_width;

[~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );
rm = rm'; % Transpose matrix ("from" in rows, "to" in columns)

assert( sum( sum( rm ) ) == 7 );
assert( rm( 5,3 ) == 2 );
assert( rm( 6,3 ) == 1 );
assert( rm( 1,4 ) == 1 );
assert( rm( 2,4 ) == 1 );
assert( rm( 1,6 ) == 2 );

assert( isequal(re,[2;6;1;5;2] ) );

save( name, 'rm', 're' );


%% Long data series
rng(0);  % Init random seed
name         = 'long_series';
x            = export_series( name, cumsum( randn( 1e4, 1 ) ) );
x_max        = max(x);
x_min        = min(x);
class_count  = 100;
class_width  = round( (x_max - x_min) / (class_count - 1), 2 );
class_offset = x_min - class_width / 2;
hysteresis   = class_width;

[~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis );

save( name, 'rm', 're' );

close all
plot( re );

figure
surface( rm );

pd
