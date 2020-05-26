function [xOpt,uOpt]=Marz_Rocket_Non_Linear(m,l,J,g,TS,umin,umax,xmin,xmax,N,x0,Af,bf)

% NON Linear solver for MARZ Rocket

yalmip('clear')

feas = false;

xN = [0;0;0;0;0;0];                    


x = sdpvar(6,N+1);                  % State variable for Batch
u = sdpvar(3,N);                    % Input variable for Batch 

Q = diag([5000/(xmin(1)^2) 10/(xmin(2)^2) 10/(xmax(3)^2) 1/(xmin(4)^2) 1/(xmin(5)^2) 1/(xmin(6)^2)]);  
R = diag([10/(umax(1)^2) 300/(umax(2)^2) 300/(umax(3)^2)]); 
%Q = diag([100/(thetaMin^2) 10/(wMin^2) 10/(hMax^2) 1/(vMin^2) 1/(xMin^2) 1/(vxMin^2)]);  
%R = diag([1/(FeMax)^2 10/(FthMax^2) 10/(FtlMax^2)]);

% Q = diag([1 1 1 1 1 1]);  
% R = diag([1 1 1]);  


constraints = [x(:,1) == x0, x(:,end) == xN];     % Innital and Terminal constraints 
cost = 0;
for k = 1:N
    cost = cost + x(:,k)'*Q*x(:,k) + u(:,k)'*R*u(:,k) ;   % Cost
    constraints = [constraints x(:,k+1) == RocketdynNonLinearVx(TS,m,l,J,g,x(:,k),u(:,k)), umin<= u(:,k) <= umax];  % Input and Dynamics
    constraints = [constraints xmin <= x(:,k)<= xmax];   % State constraints
    constraints = [constraints umin<= u(:,k) <= umax];
end


options = sdpsettings('verbose',1,'solver','ipopt');
sol = optimize(constraints, cost, options);


xOpt = value(x);
uOpt = value(u);

%------------------
% Analyze error flags
%------------------
% if sol.problem == 0
%     feas = true;
% else
%     xOpt = [];
%     uOpt = [];
%     disp('Batch Problem is infeasible');
%     return;
% end

%
figure;
title('Reference Trajectory')
subplot(4,1,1)
plot(linspace(0,10,N+1), xOpt(1,:))
xlabel('t')
ylabel('\theta (rad)')
subplot(4,1,2)
plot(linspace(0,10,N+1), xOpt(2,:))
xlabel('t')
ylabel('\omega (rad/s)')
subplot(4,1,3)
plot(linspace(0,10,N+1), xOpt(3,:))
xlabel('t (s)')
ylabel('h (m)')
subplot(4,1,4)
plot(linspace(0,10,N+1), xOpt(4,:))
xlabel('t (s)')
ylabel('v (m/s)')

figure()
subplot(2,1,1)
plot(linspace(0,10,N+1), xOpt(5,:))
xlabel('t (s)')
ylabel('x (m)')
subplot(2,1,2)
plot(linspace(0,10,N+1), xOpt(6,:))
xlabel('t (s)')
ylabel('vx (m/s)')

figure;
title('Reference Trajectory')
subplot(3,1,1)
plot(linspace(0,10,N), uOpt(1,:))
xlabel('t (s)')
ylabel('Fe, Engine ?')
subplot(3,1,2)
plot(linspace(0,10,N), uOpt(2,:))
xlabel('t (s)')
ylabel('Fth, Top Thruster ?')
subplot(3,1,3)
plot(linspace(0,10,N), uOpt(3,:))
xlabel('t (s)')
ylabel('Ftl, Lower Thruster ?')


%}

end
