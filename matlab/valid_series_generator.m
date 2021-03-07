function result = valid_series_generator( lbound, classwidth, classcount, cls_zero, amount, verbose )

  rfc_param = struct( 'lbound',     0, ...
                      'classwidth', 1, ...
                      'classcount', classcount, ...
                      'hysteresis', 0, ...
                      'dilation',   0, ...
                      'residuum',   1);

  build_result = @(delta) cumsum(delta) + cls_zero;
  delta = 0;
  result = [];

  % Loop until enough data is generated
  while length(result) < amount
    chunk = 1/classwidth * round( randn( 1, 1000 ) / 1.4 * classwidth );
    delta = [delta, chunk];
    while 1
      ok = 1;
      while ok
        % Signal may not exceed class range
        result = build_result( delta );
        ind = find( result >= classcount-1 | result < 0, 1 );
        if isempty(ind)
          break
        else
          delta(ind) = [];
          ok = 0;
        end
      end
      while ok
        % Avoid data near class edges
        pos = rem( result, 1 );
        ind = find( pos < 0.1 | pos > 0.9, 1 );
        if isempty(ind)
          break
        else
          delta(ind) = [];
          result = build_result(delta);
          ok = 0;
        end
      end
      while ok && ~isempty(result)
        % Proof
        r = ftc2( 'rfc', result, rfc_param );
        tp = [r.tp.values,r.tp.idx];
        % Avoid range pair counts near hysteresis threshold
        probe = diff(tp(:,1));
        ind = find( abs(probe) > 0.9 & abs(probe) < 1.1, 1 );
        if isempty(ind)
          break
        else
          delta(tp(ind,2)) = [];
          ok = 0;
        end
      end
      if ok, break, end
    end
  end

  % Ensure outermost classes hold at least one sample
  result = result(1:amount);
  [~,ind] = max( result );
  result(ind) = classcount - 0.5;
  [~,ind] = min( result );
  result(ind) = 0.5;
  
  if verbose
    figure
    plot( result, '.' )
    ylabel( 'Class edges, normalized' );
    xlabel( 'Sample number' );
    title( 'Avoid data near class edges' );
    grid
    
    figure
    plot( probe, '.' );
    ylabel( 'Range pair, normalized' );
    xlabel( 'Pair number' );
    title( 'Avoid range pair counts near hysteresis threshold' );
    grid
    
    % Zahlen im Bereich der Klassengrenzen vermeiden
    pos = rem( result, 1 );
    ind = find( pos<0.1 | pos>0.9 );
    assert( isempty(ind) );
    
    r = ftc2( 'rfc', result, rfc_param );
    tp = [r.tp.values,r.tp.idx];
    delta = diff(tp(:,1));
    ind = find( abs(delta) > 0.9 & abs(delta) < 1.1, 1 );
    assert( isempty(ind) );
    
    % Scale result 
    result = result * classwidth + lbound;
    
    figure
    plot( result, 'k-' )
    ylabel( 'Signal' );
    xlabel( 'Sample number' );
    title( 'Result' );
    grid
  end
  
end
