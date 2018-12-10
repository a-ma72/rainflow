clearvars all
clc

M = 0.3;

Sa_R_Inf = 1.0 / ( 1.0 - M );                      %/* y = -x && y = Sa(R=-1) - Mx                  */
Sa_R_0   = 1.0 / ( 1.0 + M );                      %/* y =  x && y = Sa(R=-1) - Mx                  */
Sa_R_0p5 = Sa_R_0 * ( 1.0 + M/3 ) / ( 1.0 + M );   %/* Backtrace Sa(R=0) to Sa(R=-1) with M/3, then */
                                                   %/* 3y = x && y = Sa(R=-1) - (M/3)x              */

Sa = 4;
Sm = 1;

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

