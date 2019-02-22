function r = bitconf( value, bit, val_1, val_0 )
	r = bitget( value, bit );
	if nargin > 2 && r == 1
		r = val_1;
	elseif nargin > 3 && r == 0
		r = val_0;
	end
end