function [feas, xOpt, uOpt] = solve_cftoc_track(A, B, C, P, Q, R, N, x0, xref, uref, xL, xU, uL, uU, bf, Af)

%'Fix dynamics, they are fucked'

yalmip('clear')
nx = size(A,2); %state dimension
nu = size(B,2); %input dimension

%-------------------------------------------------
% Define optimization variables
%------------------
x = sdpvar(nx,N+1);
u = sdpvar(nu,N);

feas = false;

%------------------
% Define constraints
%------------------
%final state polytope constraints
% con_tol=[0;0;1e1;1e1;0;0];
con_tol=[3;pi/180;.2;.1;.2;.1];



%%%%%%%%%%%%%% Innital Constraints
constr = [x(:,1) == x0]; 
%%%%%%%%%%%%%% Innital Constraints



%%%%%%%%%%%%%% Terminal Constraints
if isempty(Af)
    constr = [constr, bf-con_tol <= x(:,end) <= bf+con_tol]; 
else
    constr = [constr -bf<= Af*x(:,N+1) <=bf];
end
%%%%%%%%%%%%% Terminal Constraints



%%%%%%%%%%%%%% Terminal Cost
cost = (x(:,N+1)-xref(:,N+1))'*P*(x(:,N+1)-xref(:,N+1)); 
%%%%%%%%%%%%%% Terminal Cost


for k = 1:N
    constr = [constr, x(:,k+1) == A*x(:,k) + B*u(:,k) + C];                % Dynamics
    constr = [constr, uL <= u(:,k),u(:,k) <= uU];                          % Input Box constraints
    constr = [constr, xL <= x(:,k+1),x(:,k+1)<=xU];                         % State box constraints
    cost = cost + (x(:,k)-xref(:,k))' * Q * (x(:,k)-xref(:,k)) +...
                  (u(:,k)-uref(:,k))' * R * (u(:,k)-uref(:,k));            % Stage Costs
end

%------------------
% Define an objective
%------------------
%Here Objective=cost

%------------------
% Set some options for YALMIP and solver
%------------------

%     options = sdpsettings('verbose',0);
%options = sdpsettings('verbose',0,'solver','quadprog');
%Normal constraint tol is 1e-8
options = sdpsettings('solver','quadprog');

%------------------
% Solve the problem
%------------------
sol = optimize(constr,cost,options);

%------------------
% Analyze error flags
%------------------
if sol.problem == 0
    feas = true;
else
    xOpt = [];
    uOpt = [];
    disp('Problem is infeasible');
    return;
end

xOpt = double(x);
uOpt = double(u);

end