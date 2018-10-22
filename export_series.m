function rounded_data = export_series( filename, data )
  rounded_data = round( data, 2 );  % Auf 2 Nachkommastellen runden
  fid = fopen( [filename, '.csv'], 'wt' );
  if fid ~= -1
    fprintf( fid, '%.2f\n', rounded_data );
    fclose( fid );
  end
end
