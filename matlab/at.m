clearvars all
clc

M = 0.3;

Sa_R_Inf = 1.0 / ( 1.0 - M );                      %/* y = -x && y = Sa(R=-1) - Mx                  */
Sa_R_0   = 1.0 / ( 1.0 + M );                      %/* y =  x && y = Sa(R=-1) - Mx                  */
Sa_R_0p5 = Sa_R_0 * ( 1.0 + M/3 ) / ( 1.0 + M );   %/* Backtrace Sa(R=0) to Sa(R=-1) with M/3, then */
                                                   %/* 3y = x && y = Sa(R=-1) - (M/3)x              */

Sa = 4;
Sm = 1;

%           1         2          3       4             5
Sm_norm = [ -3,       -Sa_R_Inf, Sa_R_0, Sa_R_0p5 * 3, 3        ];
Sa_norm = [ Sa_R_Inf,  Sa_R_Inf, Sa_R_0, Sa_R_0p5,     Sa_R_0p5 ];

M_signed = diff(Sa_norm) ./ diff(Sm_norm);
i = 4;
Sm_norm = Sm_norm(i);
Sa_norm = Sa_norm(i);
M_signed = M_signed(i);
R = (Sm-Sa)/(Sm+Sa);

Sa_norm = ( Sa_norm - M_signed .* Sm_norm ) ./ ( 1 - M_signed .* Sm ./ Sa );
Sa/Sa_norm
ftc2( 'amptransform', Sa, Sm, M, -1, 1 )


%%
Sa = 100;
Sm = 50;
M = 0.3;
target = 400;
target_is_R = 0;

lhs = ftc2( 'amptransform', Sa, Sm, M, target, target_is_R );
rhs = rfc( 'amptransform', Sa, Sm, M, target, target_is_R );

assert( abs( lhs / rhs - 1 ) < 1e-7 );

while 1
  Sa = rand * 50;
  Sm = rand * 100 - 50;
  M  = rand * 0.9 + 0.1;
  target = 400;
  target_is_R = randi(2)-1;
  if target_is_R
    target_m = rand * 100 - 50;
    target_a = rand * 50;
    target = (target_m-target_a) / (target_m+target_a);
  else
    target = rand * 800 - 400;
  end
  lhs = ftc2( 'amptransform', Sa, Sm, M, target, target_is_R );
  rhs = rfc( 'amptransform', Sa, Sm, M, target, target_is_R );

  assert( abs( lhs / rhs - 1 ) < 1e-7 );
end
