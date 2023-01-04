function validate

  if 1
    build = 'Release';
  else
    build = 'Debug';
  end

  if exist( 'rfc', 'file' ) ~= 3
    if ispc
      addpath( ['../build/', build] );
    else
      addpath( '../build' );
    end
  end

  %% Empty series
  name              =  'empty';
  class_count       =  100;
  x                 =  export_series( name, [], class_count );
  x_max             =  1;
  x_min             = -1;
  [class_width, ...
   class_offset]    =  class_param( x, class_count );
  hysteresis        =  class_width;
  enforce_margin    =  0;  % First and last data point may be excluded in tp
  use_hcm           =  0;  % Use 4 point method, not HCM
  use_astm          =  0;  % Use 4 point method, not ASTM
  residual_method   =  0;
  spread_damage     =  0;
  auto_resize       =  0;

  [~,re,rm] = rfc( 'rfc', x, class_count, class_width, class_offset, hysteresis, ...
                          residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  assert( sum( sum( rm ) ) == 0 );

  assert( isempty(re) );

  save( name, 'rm', 're' );

  %% One single cycle (up)
  name              =  'one_cycle_up';
  class_count       =  4;
  x                 =  export_series( name, [1,3,2,4], class_count );
  x_max             =  4;
  x_min             =  1;
  [class_width, ...
   class_offset]    =  class_param( x, class_count );
  hysteresis        =  class_width * 0.99;
  enforce_margin    =  0;  % First and last data point may be excluded in tp
  use_hcm           =  0;  % Use 4 point method, not HCM
  use_astm          =  0;  % Use 4 point method, not ASTM
  residual_method   =  0;
  spread_damage     =  0;
  auto_resize       =  0;

  [~,re,rm] = rfc( 'rfc', x, class_count, class_width, class_offset, hysteresis, ...
                          residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  assert( sum( sum( rm ) ) == 1 );
  assert( rm( 3,2 ) == 1 );

  assert( isequal( re, [1;4] ) );

  save( name, 'rm', 're' );

  %% One single cycle (down)
  name              =  'one_cycle_down';
  class_count       =  4;
  x                 =  export_series( name, [4,2,3,1], class_count );
  x_max             =  4;
  x_min             =  1;
  [class_width, ...
   class_offset]    =  class_param( x, class_count );
  hysteresis        =  class_width * 0.99;
  enforce_margin    =  0;  % First and last data point may be excluded in tp
  use_hcm           =  0;  % Use 4 point method, not HCM
  use_astm          =  0;  % Use 4 point method, not ASTM
  residual_method   =  0;
  spread_damage     =  0;
  auto_resize       =  0;

  [~,re,rm] = rfc( 'rfc', x, class_count, class_width, class_offset, hysteresis, ...
                          residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  assert( sum( sum( rm ) ) == 1 );
  assert( rm( 2,3 ) == 1 );

  assert( isequal( re, [4;1] ) );

  save( name, 'rm', 're' );

  %% Small example, taken from url:
  % [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
  name              =  'small_example';
  class_count       =  6;
  x                 =  export_series( name, [2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2], class_count );
  x_max             =  max(x);
  x_min             =  min(x);
  [class_width, ...
   class_offset]    =  class_param( x, class_count );
  hysteresis        =  class_width;
  enforce_margin    =  0;  % First and last data point may be excluded in tp
  use_hcm           =  0;  % Use 4 point method, not HCM
  use_astm          =  0;  % Use 4 point method, not ASTM
  residual_method   =  0;
  spread_damage     =  0;
  auto_resize       =  0;

  [~,re,rm] = rfc( 'rfc', x, class_count, class_width, class_offset, hysteresis, ...
                          residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  assert( sum( sum( rm ) ) == 7 );
  assert( rm( 5,3 ) == 2 );
  assert( rm( 6,3 ) == 1 );
  assert( rm( 1,4 ) == 1 );
  assert( rm( 2,4 ) == 1 );
  assert( rm( 1,6 ) == 2 );

  assert( isequal(re,[2;6;1;5;2] ) );

  save( name, 'rm', 're' );

  %% Long data series
  rng(5,'twister')  % Init random seed
  name              =  'long_series';
  class_count       =  100;
  class_offset      = -2025;
  class_width       =  50;
  xx                =  valid_series_generator( class_offset, class_width, class_count, 40.5, 1e4, true );
  xx                =  export_series( name, xx, class_count, '%+ 5.0f' );
  x_int             =  int16(round(xx));
  x                 =  double(x_int);
  hysteresis        =  class_width;
  enforce_margin    =  1;  % Enforce first and last data point included in tp
  use_hcm           =  0;  % Use 4 point method, not HCM
  use_astm          =  0;  % Use 4 point method, not ASTM
  residual_method   =  0;  % 0=RFC_RES_NONE, 7=RFC_RES_REPEATED
  spread_damage     =  1;  % 0=RFC_SD_HALF_23, 1=RFC_SD_RAMP_AMPLITUDE_23
  auto_resize       =  0;

  assert( sum(abs(xx-x)) < 1e-4 );

  [pd,re,rm,rp,lc,tp,dh] = ...
    rfc( 'rfc', double(x_int), class_count, class_width, class_offset, hysteresis, ...
                residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  % With residuum:    pd == 9.8934e-06 (repeated)
  % Without residuum: pd == 1.1486e-07
  assert( abs( sum( tp(:,3) ) / pd - 1 ) < 1e-10 );

  save( name, 'rm', 're' );

  close all
  figure
  [ax, h1, h2] = plotyy( tp(:,1), tp(:,2), tp(:,1), cumsum( tp(:,3) ) );
  set( h1, 'DisplayName', 'time series' );
  set( h2, 'DisplayName', 'cumulative damage (based on TP)' );

  grid

  spread_damage = 8;  % 7=RFC_SD_TRANSIENT_23, 8=RFC_SD_TRANSIENT_23c

  [pd,re,rm,rp,lc,tp,dh] = ...
    rfc( 'rfc', x, class_count, class_width, class_offset, hysteresis, ...
                residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

  assert( abs( sum( dh ) / pd - 1 ) < 1e-10 );
  hold( ax(2), 'all' );
  plot( ax(2), 1:length(dh), cumsum(dh), 'k--', 'DisplayName', 'cumulative damage (based on time series)' )
  plot( ax(2), tp(:,1), cumsum(tp(:,3)), 'g-.', 'DisplayName', 'mapped on turning points' );
  legend( 'show' )

  figure
  plot( re );
  title( 'Residuum' );
  xlabel( 'Index' );
  ylabel( 'Value' );

  figure
  surface( rm );
  axis ij
  xlabel( 'to class' );
  ylabel( 'from class' );
  title( 'Rainflow matrix' );

  disp( pd )

  assert( strcmp( sprintf( '%.4e', pd ), '1.1486e-07' ) )
  assert( sum( rm(:) ) == 640 );
  assert( length(re) == 10 );
  test = re(:)' - [0,142,-609,2950,-2000,2159,1894,2101,1991,2061];
  test = sum( abs( test ) );
  assert( all( test < 1e-3 ))

  if ~isempty( which( 'ftc2' ) )
      y = ftc2( 'rfc', x, 'classcount', 100, 'classwidth', class_width, 'lbound', class_offset, ...
                'dilation', 0, 'hysteresis', class_width, 'residuum', 1 );

      assert( abs( pd / y.bkz - 1 ) < 1e-10 );
      test = sum( abs( y.residuum(:) - re(:) ) );
      assert( all( test < 1e-1 ))
  end

  %% Compare with ASTM E 1049-85 (MATLAB)
  if 0
    enforce_margin    =  1;  % Enforce first and last data point included in tp
    use_hcm           =  0;  % Use 4 point method, not HCM
    use_astm          =  0;  % Use 4 point method, not ASTM
    residual_method   =  4;  % 4=ASTM related
    spread_damage     =  0;  % 0=RFC_SD_HALF_23, 1=RFC_SD_RAMP_AMPLITUDE_23
    auto_resize       =  0;
    class_count       =  1000;
    class_width       =  5;
    class_offset      = -2025;
    hysteresis        =  class_width;

    if 1
      load long_series_csv.mat
    else
      x_int             = [2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2];
      x_max             =  max(x_int)+0.5;
      x_min             =  min(x_int)-0.5;
      class_width       =  1;
      class_offset      =  x_min;
      hysteresis        =  0;
      class_count       =  (x_max - x_min) / class_width;
    end

   residual_method = 7; % 7=repeated
   [pd7,re7,rm7,rp7,lc7,tp7] = ...
      rfc( 'rfc', double(x_int), class_count, class_width, class_offset, hysteresis, ...
                  residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

   residual_method = 4;  % 4=ASTM related
   [pd4,re4,rm4,rp4,lc4,tp4] = ...
      rfc( 'rfc', double(x_int), class_count, class_width, class_offset, hysteresis, ...
                  residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

   use_astm = 1; % Use ASTM algorithm
   residual_method = 4;  % 4=ASTM related
   [pd4a,re4a,rm4a,rp4a,lc4a,tp4a] = ...
      rfc( 'rfc', double(x_int), class_count, class_width, class_offset, hysteresis, ...
                  residual_method, enforce_margin, use_hcm, use_astm, spread_damage, auto_resize );

    % MATLAB - Rainflow counts for fatigue analysis (according to ASTM E 1049)
    c = rainflow( tp4a(:,2), 'ext' );
    edges = (0:class_count) .* class_width;
    [~,bin] = histc( c(:,2), edges );
    N = accumarray( bin, c(:,1) );
    N(class_count+1) = 0;
    N = N(1:class_count);
    Range = (edges(1:end-1)+edges(2:end)) / 2;
    pd_astm = sum( N(:)' ./ (1e7*(Range/2/1e3).^-5) );

    figure
    plot( cumsum(rp7,  'reverse'), edges(1:end-1), 'g-',  'Disp', '4-point method, res=repeated' ), hold all
    plot( cumsum(rp4,  'reverse'), edges(1:end-1), 'k--', 'Disp', '4-point method, res=half cycles' ), hold all
    plot( cumsum(rp4a, 'reverse'), edges(1:end-1), 'b--', 'Disp', '3-point method, ASTM E 1049-85, res=half cycles' ), hold all
    plot( cumsum(N,    'reverse'), edges(1:end-1), 'r-',  'Disp', '3-point method, ASTM E 1049-85 (MATLAB)' ), hold all
    set( gca, 'XScale', 'log' );
    set( gca, 'YScale', 'log' );
    xlim( [0.9 1e3] );
    xlabel( 'Counts' );
    ylabel( 'Range (normalized, prepared)' );
    grid
    legend
    fprintf( 'Damage ratio 4pt,repeated vs. ASTM: %g%%\n', pd7/pd_astm );
  end
  
  %% Long series, turning points only
  y = rfc( 'turningpoints', x, class_width*2, enforce_margin );
  figure
  plot( x, 'k-', 'DisplayName', 'time series' );
  hold all
  plot( y(2,:), y(1,:), 'r--', 'DisplayName', 'turning points' );
  legend( 'show' );

end

function rounded_data = export_series( filename, data, class_count, format )
  if nargin < 4
    format = '%.2f';
  end
  % Round data first
  rounded_data = round( data, 2 );
  % Avoid values near class boundaries, to avoid rounding effects on
  % different machines
  [class_width, ...
   class_offset]    = class_param( rounded_data, class_count );
  % Normalized data
  rounded_data = ( rounded_data - class_offset ) / class_width;
  % Inspect boundaries
  i = mod( rounded_data, 1 );
  i = find( i > 0.999 | i < 0.001 );
  % Avoid them
  data(i) = data(i) + 0.1;
  rounded_data = round( data, 2 );
  [class_width, ...
   class_offset]    = class_param( rounded_data, class_count );
  rounded_data = ( rounded_data - class_offset ) / class_width;
  i = mod( rounded_data, 1 );
  i = find( i > 0.999 | i < 0.001 );
  assert( isempty(i) );
  rounded_data = data;
  write_series( filename, data, format );
end


function data = write_series( filename, data, format )
  fid = fopen( [filename, '.csv'], 'wt' );
  if fid ~= -1
    fprintf( fid, [format,'\n'], data );
    fclose( fid );
  end

  % Build include files
  fid = fopen( [filename, '.h'], 'wt' );
  if fid ~= -1
      fprintf( fid, '#define DATA_LEN %d\n', length(data) );
      fclose( fid );
  end
  fid = fopen( [filename, '.c'], 'wt' );
  if fid ~= -1
      fprintf( fid, 'static double data_export[] = {\n' );
      s = '';
      left = numel(data);
      i = 1;
      while left
          s = [s,sprintf(format, data(i))];
          left = left - 1;
          i = i + 1;
          if left
            s = [s,', '];
          end
          if mod(i-1,16)==0 || ~left
            fprintf( fid, '    %s\n', s );
            s= '';
          end
      end
      fprintf( fid, '};\n' );
      fprintf( fid, 'static size_t data_length = sizeof( data_export ) / sizeof(double);\n' );
  end
end


function [class_width, class_offset] = class_param( data, class_count )
  assert( class_count > 1 );

  if isempty( data )
    class_width  = 1;
    class_offset = 0;
  else
    class_width  = ( max( data ) - min( data ) ) / ( class_count - 1 );
    class_width  = ceil( class_width * 100 ) / 100;
    class_offset = floor( ( min( data ) - class_width / 2 ) * 1000 ) / 1000;
  end
end
