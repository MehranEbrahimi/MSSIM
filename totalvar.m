function [u] = totalvar(f, lambda)
   
% Isotropic heat equation
% Input: Initial temperature distribution f
% Output: Solution of the heat equation u

% Parameters
dt = 0.1;   % Time step
T = 20;     % Stopping time

%Initialize u=f.  Convert to double so we can do arithmetic.
u = double(f);
u_init = u;
[m,n] = size(u);

for t = 0:dt:T  % this is actually T/dt+1 iterations
    
   % 2nd derivatives
   u_xx = u(:,[2:n,n]) - 2*u + u(:,[1,1:n-1]);
   u_yy = u([2:m,m],:) - 2*u + u([1,1:m-1],:);
   % 2nd mixed partial derivatives
   u_xy = ( u([2:m,m],[2:n,n]  ) + u([1,1:m-1],[1,1:n-1]) ...
          - u([2:m,m],[1,1:n-1]) - u([1,1:m-1],[2:n,n]  ) ) / 4;
   u_yx = u_xy;
   % 1st partial derivatives (central difference)
   u_x = ( u(:,[2:n,n]) - u(:,[1,1:n-1]) ) / 2;
   u_y = ( u([2:m,m],:) - u([1,1:m-1],:) ) / 2;
%    u_x = u(:,[2:n,n]) - u;  % forward difference -- not good
%    u_y = u([2:m,m],:) - u;  % forward difference -- not good

[mssim, ssim_map, mssim_derivative]=new_ssim(u,u_init);

   u = u + dt * ( (u_xx.*u_y.^2 - u_x.*u_y.*u_xy - u_x.*u_y.*u_yx + u_yy.*u_x.^2) ...
               ./ (0.001 + (u_x.^2 + u_y.^2).^(3/2)) - 2*lambda*mssim_derivative );
   imshow(uint8(u)); 
   title(['t=',num2str(t)]);
   drawnow;  % MatLab will not draw iterates without this command.
end
u = uint8(u);
end