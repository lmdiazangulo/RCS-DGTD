function RCS2MAT(fileName)

clear N Nv NODETOL Nsteps time
clear nx ny nz
clear x y z 
clear v
clear ExInc EyInc EzInc
clear Ex Ey Ez  Hx Hy Hz

RCSGlobals;
RCSreader(cat(2,fileName,'.rcs'));
save(cat(2,fileName,'.mat'));

end
