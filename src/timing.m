fprintf(1,'*****  step %d;  dt = %4.4e;  time = %4.4e [%s]\n\n',step,dt./TimeScale,time./TimeScale,TimeUnits);
TTtime  = tic;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;

if     strcmp(TINT,'be1im') || step==1 || frst % first step / 1st-order backward-Euler implicit scheme
    a1 = 1; a2 = 1; a3 = 0;
    b1 = 1; b2 = 0; b3 = 0;
elseif strcmp(TINT,'bd2im') || step==2         % second step / 2nd-order 3-point backward-difference implicit scheme
    a1 = 3/2; a2 = 4/2; a3 = -1/2;
    b1 = 1;   b2 =  0;  b3 = 0;
elseif strcmp(TINT,'cn2si')                    % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
    a1 = 1;   a2 = 1;   a3 = 0;
    b1 = 1/2; b2 = 1/2; b3 = 0;
elseif strcmp(TINT,'bd2si')                    % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
    a1 = 3/2; a2 = 4/2; a3 = -1/2;
    b1 = 3/4; b2 = 2/4; b3 = -1/4;
end