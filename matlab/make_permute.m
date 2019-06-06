addpath( '.' );

clear all
clc
loop = 1;

[~] = rmdir( '../build', 's' );
mkdir( '../build' );
cd( '..' );
projdir = cd( 'build' );

if ispc
  VS_Toolset = 'Visual Studio 14 2015 Win64';
  VS_CommonTools = '%VS140COMNTOOLS%VsDevCmd.bat';
end

for loop = 2^10-1:-1:0
  VALUE_TYPE               = bitconf( loop, 1, 'float', 'double' );
  RFC_USE_INTEGRAL_COUNTS  = bitconf( loop, 2 );
  RFC_USE_HISTOGRAM_FILTER = bitconf( loop, 3 );
  RFC_MINIMAL              = bitconf( loop, 4 );
  RFC_TP_SUPPORT           = bitconf( loop, 5 );
  RFC_HCM_SUPPORT          = bitconf( loop, 6 );
  RFC_GLOBAL_EXTREMA       = bitconf( loop, 7 );
  RFC_DAMAGE_FAST          = bitconf( loop, 8 );
  RFC_DH_SUPPORT           = bitconf( loop, 9 );
  RFC_AT_SUPPORT           = bitconf( loop, 10 );
  RFC_USE_DELEGATES        = bitconf( loop, 11 );
  
  fprintf( 'Configuration %d\n', loop );
  clear defs rfc
  defs(2,:) = who( 'RFC_*' )';
  defs(1,:) = deal( {'-D'} );
  defs(3,:) = deal( {'='} );
  defs(4,:) = sprintfc( '%d', cellfun( @eval, defs(2,:) ) );
  defs(5,:) = deal( {' '} );
  
  if ispc
    fid = fopen( '../build/generate.bat', 'wt' );
    fprintf( fid, '%s\n', 'cd %~dp0' );
    fprintf( fid, '@call "%s"\n', VS_CommonTools );
    fprintf( fid, '%s\n', sprintf( '"%s\\cmake.exe" -DCMAKE_BUILD_TYPE=Debug -G"%s" %s %s %s ', ...
                                   getenv('cmake_root'), VS_Toolset, ...
                                   ['-DRFC_VALUE_TYPE=',VALUE_TYPE], ...
                                   [defs{:}], ...
                                   projdir ) );
    fprintf( fid, '%s\n', 'if %errorlevel% neq 0 exit /b %errorlevel%' );
    fprintf( fid, '%s\n', 'devenv rainflow.sln /build Release' );
    fprintf( fid, '%s\n', 'if %errorlevel% neq 0 exit /b %errorlevel%' );
    fprintf( fid, '%s\n', '.\Release\rfc_test.exe || exit /b 1' );
    fclose( fid );
    i = 0;
    status = 1;
    while status ~= 0
      status = system( '..\build\generate.bat', '-echo' );
      i = i + 1;
      if i > 10, break, end
    end
    if status ~= 0, return, end
  else
    fid = fopen( '../build/generate.sh', 'wt' );
    fprintf( fid, '%s\n', '#!/bin/bash' );
    fprintf( fid, '%s\n', 'set -e' );
    fprintf( fid, '%s\n', sprintf( 'cmake -DCMAKE_BUILD_TYPE=Debug -G"%s" %s %s %s ', ...
                                   'Unix Makefiles', ...
                                   ['-DRFC_VALUE_TYPE=',VALUE_TYPE], ...
                                   [defs{:}], ...
                                   projdir ) );
    fprintf( fid, '%s\n', 'make' );
    fprintf( fid, '%s\n', './rfc_test' );
    fclose( fid );
    system( 'chmod +x generate.sh' );
    i = 0;
    status = 1;
    while status ~= 0
      status = system( '../build/generate.sh', '-echo' );
      i = i + 1;
      if i > 10, break, end
    end
    if status ~= 0, return, end
  end
end
