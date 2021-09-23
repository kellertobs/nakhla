function    dadx = ddx(a,x)

if length(x) == 1  % constant grid spacing, h = x
    dadx = diff(a,1,2)./x;
else
    dadx = diff(a,1,2)./diff(x,1,2);
end
