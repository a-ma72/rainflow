% Create some test series and export to MS Excel(R)

assert( exist( 'ftc2', 'file' ) == 3, 'This script needs the fatigue tool chain (FTC)!' );

classcount =  25;
classwidth =  380;
lbound     = -3610;
ubound     = +5890;   % -3610 + 25 * 380
cls_zero   =  9.5;    % Class number containing zero
amount     =  10000;  % Number of samples
verbose    =  true;

classnr    = @(v) floor((v-lbound)/classwidth)+1;

close all

rng(1,'twister')  % Init random seed
Testseries(1,1:2) = {'Cycle up',            [1,3,2,4] };
Testseries(2,1:2) = {'Cycle down',          [4,2,3,1] };
Testseries(3,1:2) = {'Residue stress test', [2,3,1,4,1,3,2,3, 2,3,1,4,1,3,2,3, 2,3,1,4,1,3,2,3, 1.9] };
Testseries(4,1:2) = {'Small example',       [2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2] };
Testseries(5,1:2) = {'Long series',         valid_series_generator( lbound, classwidth, classcount, cls_zero, amount, verbose ) };

rfc_par = struct( 'lbound', lbound, 'classwidth', classwidth, 'classcount', classcount, ...
                  'dilation', 0, 'residuum', 1 );
rfc_par_no_finalize             = rfc_par;
rfc_par_no_finalize.no_finalize = true;
rfc_par_no_finalize.wl_sd       = classwidth;
rfc_par_no_finalize.wl_nd       = 1;
rfc_par_no_finalize.wl_k        = 4;

for i = 1:5
  name = Testseries{i,1};
  data = Testseries{i,2};
  if i < 5
    % Short test series with class numbers
    data = (2*(data-1)+cls_zero) * classwidth + lbound;
    r = ftc2( 'rfc', data, rfc_par );
  else
    % Long test series with real values
    r = ftc2( 'rfc', data, rfc_par_no_finalize );
  end
  
  % Build rainflow matrix from sparse
  rfc_mat = zeros( classcount + 1 );
  for j = 1:length(r.rfm)
    from = floor((r.rfm(j).from-lbound)/classwidth)+1 +1;
    to   = floor((r.rfm(j).to-lbound)/classwidth)+1 +1;
    rfc_mat( from, to ) = r.rfm(j).counts * 2;
  end
  rfc_mat = num2cell( rfc_mat );
  for j = 1:classcount
    rfc_mat{j+1,1} = sprintf( 'F%d', j );  % From
    rfc_mat{1,j+1} = sprintf( 'T%d', j );  % To
  end
  rfc_mat{1,1} = [];
  wbook = sprintf( 'xl/%s.xlsx', name );
  if i < 5
    data = num2cell( [classnr(data(:)),data(:)] );
    data = [{'ClassNr', 'Value'};data];
  else
    bkz  = ftc2( 'rfc', data, rfc_par_no_finalize );
    bkz  = arrayfun( @(n) getfield( ftc2( 'rfc', data(1:n), rfc_par_no_finalize ), 'bkz' ), 1:length(data) );
    data = num2cell( [classnr(data(:)),data(:),bkz(:)] );
    data = [{'ClassNr', 'Value', 'Pseudo damage'};data];
  end
  residue = num2cell( [classnr(r.residuum(:)),r.residuum(:)] );
  residue = [{'ClassNr', 'Value'};residue];
  % Write to Excel
  writecell( data,    wbook, 'Sheet', 'Data' );
  writecell( rfc_mat, wbook, 'Sheet', 'RF-Matrix' );
  writecell( residue, wbook, 'Sheet', 'Residuum' );
  
  % Pseudo damage:
  % =SUM( B2:Z26/2 *( ABS( ROW(A2:A26)-COL(B1:Z1) ) / 2 )^5 )
end
