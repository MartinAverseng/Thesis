function [Y] = sph2cart(theta_x,phi_x,theta_y,phi_y,j)
%sph2cart Convert spherical (angular) coordinates to cartesian
%   theta_y and phi_y are the angular coordinates of a vector y on the
%   sphere, where the angles are expressed in the base (O,x) where x is a
%   point on the sphere given by its angles in the base (O,[0;0;1]) theta_x
%   and phi_x.  The function returns the cartesian coordinates of y.
%   theta_x, phi_x, theta_y, phi_y can be Nx1 arrays (of same size) or
%   theta_x and phi_x can be scalars with theta_y and phi_y arrays.

if ~exist('j','var')
    j = [];
end

if size(theta_y,2)==1
    erx = [sin(theta_x).*cos(phi_x), sin(theta_x).*sin(phi_x), cos(theta_x)];
    ethetax = [cos(theta_x).*cos(phi_x), cos(theta_x).*sin(phi_x), -sin(theta_x)];
    ephix = [-sin(phi_x), cos(phi_x), 0*phi_x];
    
    
    y1 = sin(theta_y).*cos(phi_y).*ethetax(:,1) + sin(theta_y).*sin(phi_y).*ephix(:,1) + cos(theta_y).*erx(:,1);
    y2 = sin(theta_y).*cos(phi_y).*ethetax(:,2) + sin(theta_y).*sin(phi_y).*ephix(:,2) + cos(theta_y).*erx(:,2);
    y3 = sin(theta_y).*cos(phi_y).*ethetax(:,3) + sin(theta_y).*sin(phi_y).*ephix(:,3) + cos(theta_y).*erx(:,3);
    
    
    Y = [y1 y2 y3];
    if ~isempty(j)
        Y = Y(:,j);
    end
    
else
    Y = sph2cart(theta_x,phi_x,theta_y.',phi_y.',j).';
end

end



