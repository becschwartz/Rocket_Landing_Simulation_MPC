function xnext = RocketdynNonLinearVxDist(TS,m,l,J,g,x,u,Wind)
% x = [theta, w, h, v, x, vx]
% u = [Fe Fth Ftl]
% Small Angle Aproximation has been made
xnext(1,1) = x(1) + TS*x(2);                                % Theta
xnext(2,1) = x(2) + TS*(l/J)*[u(2)/4 - u(3)/2];             % w
xnext(3,1) = x(3) + TS*x(4);                                % h
xnext(4,1) = x(4) + TS*((u(1)-u(2)*x(1)-u(3)*x(1))/m-g);    % v
xnext(5,1) = x(5) + TS*x(6);                                % x
xnext(6,1) = x(6) + TS*(-u(1)*x(1) - u(2) - u(3) + Wind)/m;        % vx
end