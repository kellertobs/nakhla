function    dadz = ddz(a,z)

if length(z) == 1  % constant grid spacing, h = z
    dadz = diff(a,1,1)./z;
else
    dadz = diff(a,1,1)./diff(z,1,1);
end
