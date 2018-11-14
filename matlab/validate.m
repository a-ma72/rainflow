function validate
  %% Empty series
  name              = 'empty';
  class_count       = 100;
  x                 = export_series( name, [], class_count );
  x_max             = 1;
  x_min             = -1;
  class_width       = round( (x_max - x_min) / (class_count - 1), 2 );
  class_offset      = x_min - class_width / 2;
  hysteresis        = class_width;
  enforce_margin    = 0;
  use_hcm           = 0;

  [~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis, ...
                   enforce_margin, use_hcm );

  assert( sum( sum( rm ) ) == 0 );

  assert( isempty(re) );

  save( name, 'rm', 're' );


  %% One single cycle (up)
  name              = 'one_cycle_up';
  class_count       = 4;
  x                 = export_series( name, [1,3,2,4], class_count );
  x_max             = 4;
  x_min             = 1;
  class_width       = round( (x_max - x_min) / (class_count - 1), 2 );
  class_offset      = x_min - class_width / 2;
  hysteresis        = class_width * 0.99;
  enforce_margin    = 0;
  use_hcm           = 0;

  [~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis, ...
                   enforce_margin, use_hcm );

  assert( sum( sum( rm ) ) == 1 );
  assert( rm( 3,2 ) == 1 );

  assert( isequal( re, [1;4] ) );

  save( name, 'rm', 're' );


  %% One single cycle (down)
  name              = 'one_cycle_down';
  class_count       = 4;
  x                 = export_series( name, [4,2,3,1], class_count );
  x_max             = 4;
  x_min             = 1;
  class_count       = 4;
  class_width       = round( (x_max - x_min) / (class_count - 1), 2 );
  class_offset      = x_min - class_width / 2;
  hysteresis        = class_width * 0.99;
  enforce_margin    = 0;
  use_hcm           = 0;

  [~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis, ...
                   enforce_margin, use_hcm );

  assert( sum( sum( rm ) ) == 1 );
  assert( rm( 2,3 ) == 1 );

  assert( isequal( re, [4;1] ) );

  save( name, 'rm', 're' );


  %% Small example, taken from url:
  % [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
  name              = 'small_example';
  class_count       = 6;
  x                 = export_series( name, [2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2], class_count );
  x_max             = max(x);
  x_min             = min(x);
  class_width       = round( (x_max - x_min) / (class_count - 1), 2 );
  class_offset      = x_min - class_width / 2;
  hysteresis        = class_width;
  enforce_margin    = 0;
  use_hcm           = 0;

  [~,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis, ...
                   enforce_margin, use_hcm );

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
  name              = 'long_series';
  class_count       = 100;
  x                 = export_series( name, cumsum( randn( 1e4, 1 ) ), class_count );
  x_max             = max(x);
  x_min             = min(x);
  class_width       = round( (x_max - x_min) / (class_count - 1), 2 );
  class_offset      = x_min - class_width / 2;
  hysteresis        = class_width;
  enforce_margin    = 0;
  use_hcm           = 0;

  [pd,re,rm] = rfc( x, class_count, class_width, class_offset, hysteresis, ...
                    enforce_margin, use_hcm );

  save( name, 'rm', 're' );

  close all
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
  
  assert( strcmp( sprintf( '%.4e', pd ), '4.8703e-16' ) )
  assert( sum( sum(rm) ) == 601 );
  assert( length(re) == 10 );
  test = re ./ [0.538;2.372;-0.448;17.445;-50.901;114.136;-24.851;31.002;-0.646;16.594];
  test = abs( test - 1 );
  assert( all( test < 1e-3 ))
  
  if ~isempty( which( 'ftc2' ) )
      y = ftc2( 'rfc', x, 'classcount', 100, 'classwidth', class_width, 'lbound', class_offset, ...
                'dilation', 0, 'hysteresis', class_width, 'residuum', 1 );
            
      assert( abs( pd / y.bkz - 1 ) < 1e-10 );
      test = abs( y.residuum ./ re - 1 );
      assert( all( test < 1e-1 ))
  end
end

function rounded_data = export_series( filename, data, class_count )
  % Round data first
  rounded_data = round( data, 2 );
  % Avoid values near class boundaries, to avoid rounding effects on
  % different machines
  class_width = ( max(rounded_data) - min(rounded_data) ) / (class_count - 1);
  class_offset = min(rounded_data) - class_width/2;
  % Normalized data
  rounded_data = ( rounded_data - class_offset) / class_width;
  % Inspect boundaries
  i = mod( rounded_data, 1 );
  i = find( i > 0.999 | i < 0.001 );
  % Avoid them
  data(i) = data(i) + 0.1;
  rounded_data = round( data, 2 );
  class_width = ( max(rounded_data) - min(rounded_data) ) / (class_count - 1);
  class_offset = min(rounded_data) - class_width/2;
  rounded_data = ( rounded_data - class_offset) / class_width;
  i = mod( rounded_data, 1 );
  i = find( i > 0.999 | i < 0.001 );
  assert( isempty(i) );
  rounded_data = data;
  fid = fopen( [filename, '.csv'], 'wt' );
  if fid ~= -1
    fprintf( fid, '%.2f\n', rounded_data );
    fclose( fid );
  end
  
  % Build include file
  fid = fopen( [filename, '.c'], 'wt' );
  if fid ~= -1
      fprintf( fid, 'static double data_export[] = {\n' );
      s = '';
      left = numel(rounded_data);
      i = 1;
      while left
          s = [s,sprintf('%.2f', rounded_data(i))];
          left = left - 1;
          i = i + 1;
          if left
            s = [s,', '];
          end
          if mod(i,16)==0 || ~left
            fprintf( fid, '    %s\n', s );
            s= '';
          end 
      end
      fprintf( fid, '};\n' );
  end
end
