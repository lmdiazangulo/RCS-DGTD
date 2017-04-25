function [] = RCSvisor(x, y, z, fx, fy, fz, time)

for i=1:numel(time)

    quiver3(x,y,z,fx(:,i),fy(:,i),fz(:,i));
    drawnow;

end

end