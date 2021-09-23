vz = W + WBG;
vx = U + UBG;
wp = vz(2:end  ,2:end-1);
wm = vz(1:end-1,2:end-1);
up = vx(2:end-1,2:end  );
um = vx(2:end-1,1:end-1);

a   = 1/f0 - f;
acc = a(2:end-1,2:end-1);
ajp = a(3:end  ,2:end-1);  ajpp = a([4:end,3      ],2:end-1);
ajm = a(1:end-2,2:end-1);  ajmm = a([end-2,1:end-3],2:end-1);
aip = a(2:end-1,3:end  );  aipp = a(2:end-1,[4:end,3      ]);
aim = a(2:end-1,1:end-2);  aimm = a(2:end-1,[end-2,1:end-3]);

Div_fV(2:end-1,2:end-1)   =     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
                          - abs(up).*(-(aipp-aip)./h./8 + (aip - acc)./h./4 - (acc-aim )./h./8) ...
                          -     um .*(-(aip -acc)./h./8 + (acc + aim)./h./2 + (aim-aimm)./h./8) ...
                          + abs(um).*(-(aip -acc)./h./8 + (acc - aim)./h./4 - (aim-aimm)./h./8) ...
                          +     wp .*(-(ajpp-ajp)./h./8 + (ajp + acc)./h./2 + (acc-ajm )./h./8) ...
                          - abs(wp).*(-(ajpp-ajp)./h./8 + (ajp - acc)./h./4 - (acc-ajm )./h./8) ...
                          -     wm .*(-(ajp -acc)./h./8 + (acc + ajm)./h./2 + (ajm-ajmm)./h./8) ...
                          + abs(wm).*(-(ajp -acc)./h./8 + (acc - ajm)./h./4 - (ajm-ajmm)./h./8);

um   =  min(U(jc,ic),0);
up   =  max(U(jc,ic),0);
wm   =  min(W(jc,ic),0);
wp   =  max(W(jc,ic),0);

umo  =  min(uo(jc,ic),0);
upo  =  max(uo(jc,ic),0);
wmo  =  min(wo(jc,ic),0);
wpo  =  max(wo(jc,ic),0);
    
da_dzm =   3/2 * diff(a(2:end-2,ic),1,1)./diff(zz(2:end-2,ic),1,1) ...
    - 1/2 * diff(a(1:end-3,ic),1,1)./diff(zz(1:end-3,ic),1,1);
da_dzp =   3/2 * diff(a(3:end-1,ic),1,1)./diff(zz(3:end-1,ic),1,1) ...
    - 1/2 * diff(a(4:end-0,ic),1,1)./diff(zz(4:end-0,ic),1,1);
da_dxm =   3/2 * diff(a(jc,2:end-2),1,2)./diff(xx(jc,2:end-2),1,2) ...
    - 1/2 * diff(a(jc,1:end-3),1,2)./diff(xx(jc,1:end-3),1,2);
da_dxp =   3/2 * diff(a(jc,3:end-1),1,2)./diff(xx(jc,3:end-1),1,2) ...
    - 1/2 * diff(a(jc,4:end-0),1,2)./diff(xx(jc,4:end-0),1,2);

adv     =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
