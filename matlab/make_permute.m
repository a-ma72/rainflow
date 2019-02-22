assert( ispc )
addpath( '.' );

clear all
clc
loop = 1;

[~] = rmdir( '../build', 's' );
mkdir( '../build' );
cd( '..' );
projdir = cd( 'build' );
VS_Toolset = 'Visual Studio 14 2015 Win64';
VS_CommonTools = '%VS140COMNTOOLS%VsDevCmd.bat';

loop = uint32(0);

%for loop = 2^10-1:-1:0
for loop = 364:-1:66
  VALUE_TYPE              = bitconf( loop, 1, 'float', 'double' );
  RFC_USE_INTEGRAL_COUNTS = bitconf( loop, 2 );
  RFC_MINIMAL             = bitconf( loop, 3 );
  RFC_TP_SUPPORT          = bitconf( loop, 4 );
  RFC_HCM_SUPPORT         = bitconf( loop, 5 );
  RFC_GLOBAL_EXTREMA      = bitconf( loop, 6 );
  RFC_DAMAGE_FAST         = bitconf( loop, 7 );
  RFC_DH_SUPPORT          = bitconf( loop, 8 );
  RFC_AT_SUPPORT          = bitconf( loop, 9 );
  RFC_USE_DELEGATES       = bitconf( loop, 10 );
  
  fprintf( 'Configuration %d\n', loop );
  clear defs rfc
  defs(2,:) = who( 'RFC_*' )';
  defs(1,:) = deal( {'-D'} );
  defs(3,:) = deal( {'='} );
  defs(4,:) = sprintfc( '%d', cellfun( @eval, defs(2,:) ) );
  defs(5,:) = deal( {' '} );
  fid = fopen( '../build/generate.bat', 'wt' );
  fprintf( fid, '%s\n', 'cd %~dp0' );
  fprintf( fid, '@call "%s"\n', VS_CommonTools );
  fprintf( fid, '%s\n', sprintf( '"%s\\cmake.exe" -DCMAKE_BUILD_TYPE=Debug -G"%s" %s %s %s" ', ...
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
end
