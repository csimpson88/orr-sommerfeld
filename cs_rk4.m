
function [ t, w ] = cs_rk4( fun, t0, tf, y0, h )
% cs_rk4 Runge-Kutta 4th Order Numerical Integration
%   [ t, y ] = cs_rk4( fun, t0, tf, y0, h )
%
%   input   description
%   -----   -----------
%   fun     function handle for ydot function
%   t0      initial time
%   tf      final time
%   y0      starting value(s) (column vector)
%   h       time step (delta time)
%
%   output  description
%   ------  -----------
%   t       row vector containing solution times
%   y       row vector(s) containing solution values
%
% Example:
%   [ t, y ] = cs_rk4( @fun, 1, 3, 0, .2 )
%
%   function yp = fun( t, y )
%       yp = 1 + y/t + (y/t)^2;;
%   end
%
%
% Programmed by:
%   Christopher Simpson
%   27 January 2015
%   www.cdsimpson.net

t = [t0:h:tf];
w = zeros(length(y0),length(t)); % initialize w array
w(:,1) = y0;
% fprintf('\nstp\teqn\tti\t\t\twi\t\t\tk1\t\t\tk2\t\t\tk3\t\t\tk4\t\t\twi+1\n');
% fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:(length(t)-1);
    ti = t(i);
    wi = w(:,i);
    k1 = h*fun(ti,wi);
    k2 = h*fun(ti+h/2,wi+k1/2);
    k3 = h*fun(ti+h/2,wi+k2/2);
    k4 = h*fun(ti+h,wi+k3);
    
    w(:,i+1) = wi + 1/6*(k1+2*k2+2*k3+k4);
%     for j = 1:length(w(:,i))
%         fprintf('%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n',i,j,ti,wi(j),k1(j),k2(j),k3(j),k4(j),w(j,i+1));
%     end
%     fprintf('----------------------------------------------------------------------------------------\n');
end

end
