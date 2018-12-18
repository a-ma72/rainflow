if 0
PfadName = '\\vwagwofeg204\kundenkollektive$\01_Projekte\05_Streckenvergleiche\2016_IMSBF_137004_Kundenstrecke_WOB\10_Messdaten\03_rpc\EWP';
DispName = 'EWP VW3703-2053 1%K 0p1%BT';
Datei    = 'VW370_3-2053_EWP_80P_To_SF_17M_138070_0177_0180_0182_0170.rsp';
wdh      = [215]; 

x = ftc2( 'acq.load', fullfile( PfadName, Datei ), 'AXKASP' );
x = x.AXKASP.Data;
else
  load datarow
end

x_min  = min(x);
x_max  = max(x);
cc     = 100;
cw     = (x_max-x_min)/(cc-1);
lbound = x_min - cw/2;
ubound = x_min + cw/2;
cm     = (0:cc-1) * cw + lbound;
cn     = @(v) floor( (v-lbound)/cw ) + 1;
c0     = find( cm < 0.0, 1, 'last' );

x = floor( (x-lbound)/cw );
x = (x*1.001+0.5) * cw + lbound;

x_orig = x;
hysteresis = 0;
enforce_margin = 0;
use_hcm = 0;


%%


[bkz,res,rm,rp,lc,tp] = rfc( 'rfc', x_orig, cc, cw, lbound, hysteresis, ...
                             0, enforce_margin, use_hcm ); % 0 = None, 6 = Repeated
rm_none  = rm;
res_none = res;

[bkz,res,rm,rp,lc,tp] = rfc( 'rfc', x_orig, cc, cw, lbound, hysteresis, ...
                             6, enforce_margin, use_hcm ); % 0 = None, 6 = Repeated
rm_rep = rm;
res_rep = res;

[bkz,res,rm,rp,lc,tp] = rfc( 'rfc', res_none, cc, cw, lbound, hysteresis, ...
                             6, enforce_margin, use_hcm ); % 0 = None, 6 = Repeated
rm_res_rep = rm;

[bkz,res,rm,rp,lc,tp] = rfc( 'rfc', [res_none,res_none], cc, cw, lbound, hysteresis, ...
                             0, enforce_margin, use_hcm ); % 0 = None, 6 = Repeated
rm_res_dup = rm;  

clear rfc

Z1 = rm_none + rm_res_rep;
Z2 = rm_none + rm_res_dup;

% rm_rep zaehlt ein Zyklus mehr!
sum( abs( Z1(:) - rm_rep(:) ) )
sum( abs( Z2(:) - rm_rep(:) ) )

y = ftc2( 'rfc', x, 'classcount', cc, 'classwidth', cw, 'lbound', lbound, ...
          'dilation', 0, 'hysteresis', hysteresis, 'residuum', 1 ); % 1 = None, 0/5 = Repeated
y_rfm = zeros(cc);
for i = 1:numel( y.rfm )
  y_rfm( cn(y.rfm(i).from), cn(y.rfm(i).to) ) = y.rfm(i).counts;
end

sum( abs( rm_none(:) - y_rfm(:) ) )
sum( abs( rm_none(:) - y_rfm(:) ) )

y = ftc2( 'rfc', x, 'classcount', cc, 'classwidth', cw, 'lbound', lbound, ...
          'dilation', 0, 'hysteresis', hysteresis, 'residuum', 5 ); % 1 = None, 0/5 = Repeated
y_rfm = zeros(cc);
for i = 1:numel( y.rfm )
  y_rfm( cn(y.rfm(i).from), cn(y.rfm(i).to) ) = y.rfm(i).counts;
end

% Z1 und Z2 sind mit ftc Ergebnis identisch
sum( abs( Z1(:) - y_rfm(:) ) )
sum( abs( Z2(:) - y_rfm(:) ) )
