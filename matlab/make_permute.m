%clear rfc
%clc
%mex -g -v -output rfc -DRFC_USE_INTEGRAL_COUNTS=0 rainflow.c
clear all
for RFC_USE_INTEGRAL_COUNTS = 0:1
  for RFC_MINIMAL = 0:1
    for RFC_TP_SUPPORT = 0:1
      for RFC_HCM_SUPPORT = 0:1
        for RFC_USE_DELEGATES = 0:1
          for RFC_GLOBAL_EXTREMA = 0:1
            for RFC_DAMAGE_FAST = 0:1
              for RFC_DH_SUPPORT = 0:1
                clear defs rfc
                defs(2,:) = who( 'RFC_*' )';
                defs(1,:) = deal( {'-D'} );
                defs(3,:) = deal( {'='} );
                defs(4,:) = sprintfc( '%d', cellfun( @eval, defs(2,:) ) );
                defs(5,:) = deal( {' '} );
                fid = fopen( '../build/generate.bat', 'wt' );
                fprintf( fid, '%s\n', 'cd %~dp0' );
                fprintf( fid, '%s\n', '@call "%VS140COMNTOOLS%VsDevCmd.bat"' );
                fprintf( fid, '%s\n', sprintf( '"%s\\cmake.exe" %s %s" ', ...
                                               getenv('cmake_root'), [defs{:}], ...
                                               'd:\git\gitlab\rainflow\build' ) );
                fprintf( fid, '%s\n', 'devenv rainflow.sln /build Release' );
                fprintf( fid, '%s\n', 'if %errorlevel% neq 0 exit /b %errorlevel%' );
                fprintf( fid, '%s\n', '.\Release\rfc_test.exe || exit /b' );
                fclose( fid );
                status = system( '..\build\generate.bat', '-echo' );
                if status ~= 0, return, end
              end
            end
          end
        end
      end
    end
  end
end