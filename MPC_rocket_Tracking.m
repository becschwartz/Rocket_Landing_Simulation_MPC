
% Todo
% Do limmit constraints outside the trajectory following

%************************************************************
%************************************************************
%Take everything out below this to the stars after debugging and turn back 
%into function

clear all
close all

% Initial conditions
theta0     = 0*pi/180; %10*pi/180;
w0         = 0;
h0         = 1500;        % 1228 m above the ground
v0         = -300;        % moving 205.2 m/s downwards
xo         = 50;
vx0        = 0;

x0 = [theta0; w0; h0; v0; xo; vx0];             % Inital conditions


M=150;
N=30;

%************************************************************
%************************************************************

%M is the problem horizon
%N is the receding horizon
%x0 is the intial state
%P is the final state weight

%-----------------------------
%Define problem variables 
%-----------------------------

%Rocket constants
m = 120000;
l = 50;
%m = 27648;                % 27648 Mass of the Rocket [kg]
%l = 50;                   % Length of the Rocket [m]
J = 1/16*m*l^2;           % Rocket Moment of Inertia MIGHT WANT TO RECALCUATE THIS for thick thing
g = 3.711;                 % Gravity  on mars [m/s]
TS=0.1; %sample time
Wind=4000;

%Linearized model matrices
A=[1  TS   0   0   0   0;
   0   1   0   0   0   0;
   0   0   1   TS  0   0;
   0   0   0   1   0   0;
   0   0   0   0   1   TS;
   0   0   0   0   0   1 ];

B=[   0        0               0;
      0   TS*l/(4*J)    -TS*l/(2*J);
      0        0               0;
    TS/m       0               0;
      0        0               0;
      0      -TS/m           -TS/m];

C=[   0;
      0;
      0;
      -g*TS;
      0;
      0];


nx = size(A,2); %state dimension
nu = size(B,2); %input dimension

% Input Constraints
FeMin       = 0;
FeMax       = 14400*1000; %6 Rapter Engines @2400 kN each
FthMin      = -2*73000;
FthMax      = 2*73000; %Draco engines same as dragon
FtlMin      = -2*73000;
FtlMax      = 2*73000;

% State Constraints
thetaMin    = -2*pi/180;    % rad
thetaMax    =  2*pi/180;    % rad
wMin        = -10;          % rad/s
wMax        =  10;          % rad/s
hMin        = 0;            % m
hMax        = 3000;          % m 
vMin        = -500;          % m/s
vMax        =  100;          % m/s
xMin        = -500;         % m
xMax        =  500;         % m
vxMin       = -100;          % m/s
vxMax       =  100;          % m/s

umin = [FeMin; FthMin; FtlMin];                
umax = [FeMax; FthMax; FtlMax];                
xmin = [thetaMin; wMin; hMin; vMin; xMin; vxMin];   
xmax = [thetaMax; wMax; hMax; vMax; xMax; vxMax];

% constraint sets represented as polyhedra
X = Polyhedron('lb',xmin,'ub',xmax);
U = Polyhedron('lb',umin,'ub',umax);



% umin = [FeMin; FthMin; FtlMin];                
% umax = [FeMax; FthMax; FtlMax];                
% xmin = [thetaMin; wMin; hMin; vMin; xMin; vxMin];   
% xmax = [thetaMax; wMax; hMax; vMax; xMax; vxMax];

% Normalised Cost Matricies Q, R, P
% x = [theta, w, h, v, x, vx]
% u = [Fe Fth Ftl]
Q = diag([3/(xmin(1)^2) 1/(xmin(2)^2) 1/(xmax(3)^2) 1/(xmin(4)^2) 1/10*(xmin(5)^2) 1/(xmin(6)^2)]);  
R = diag([1/(umax(1)^2) 10/(umax(2)^2) 10/(umax(3)^2)]); 
 %Q = diag([100/(thetaMin^2) 10/(wMin^2) 10/(hMax^2) 1/(vMin^2) 1/(xMin^2) 1/(vxMin^2)]);  
 %R = diag([1/(FeMax)^2 10/(FthMax^2) 10/(FtlMax^2)]);  
P = 10*Q; 

%-----------------------------------------------------------
% Initialize the state, input, and prediction error arrays
%-----------------------------------------------------------
xOpt = zeros(nx,M+1); 
uOpt = zeros(nu,M);

xOpt(:,1) = x0; %storing the initial state

xPred = zeros(nx,N+1,M); %Open loop predicted states
predErr = zeros(nx,M-N+1); %Prediction error between open loop and closed loop actual states

feas = false([1,M]);

Afterm = eye(6);
bfterm = [3;pi/180;.2;.1;.2;.1];


%Finding the tracking point
[x_ref, u_ref ]=Marz_Rocket_Non_Linear(m,l,J,g,TS,umin,umax,xmin,xmax,M,x0,Afterm,bfterm); 

if isempty(x_ref)
   fprintf('YOU SUCK PROBLEM INFEASABLE\n')
end

%%

x_end = [0;0;0;0;0;0];
u_end = [m*g;0;0];
x_ref=[x_ref repmat(x_end(:,end),1,N)];
u_ref=[u_ref repmat(u_end(:,end),1,N)];
 % repmat(x_ref(:,end),N))
% figure
%---------------------------
%MPC loop
%---------------------------
for t = 1:M
    
    fprintf('Solving simstep: %i\n',t)
    
    %--------------------------
    %Set the Desired Final set (X_ref N steps ahead)
    %--------------------------
    xN = x_ref(:, t+N);
    if (t == M)
        Af = Afterm;
        bf = bfterm;
    else
        bf=xN;
        Af=[];
    end

    %--------------------------
    %Function below solving the receding horizon problem
    %x
    x_ref(:,t:t+N);
    [feas(t), x, u] = solve_cftoc_track(A,B,C,P,Q,R,N,xOpt(:,t),x_ref(:,t:t+N),...
                                        u_ref(:,t:t+N-1),xmin,xmax,umin,umax,bf,Af);    
    %***********************************
    
    % ~ means "not", or if feas(t)==false
    if ~feas(t)
        warning('MPC problem infeasible--exiting simulation')
%         xOpt = [];
%         uOpt = [];
%         predErr = [];
        return;
    end
    
    % Save open loop predictions
    %xPred(:,:,t) = x;
    
    % Save closed loop trajectory
%     xOpt(:,t+1) = x(:,2); %-The initial state is in x(:,1)
%                           %-The propagated forward step is in x(:,2) from
%                           %dynamics
%                           %-The predicted states are in x(:,3:N)
%     xNext = RocketdynNonLinearVx(TS,m,L,J,g,x(:,1),u(:,1));
%             RocketdynNonLinearVx(TS,m,l,J,g,x,u)

    %xOpt(:,t+1) = RocketdynNonLinearVx(TS,m,l,J,g,xOpt(:,t),u(:,1));   % Now we use non linear dynamics
    xOpt(:,t+1) = RocketdynNonLinearVxDist(TS,m,l,J,g,xOpt(:,t),u(:,1),Wind);
    uOpt(:,t) = u(:,1); %The input we use at time t
    
    % Plot Open Loop    
    %---------------------------------------------
    %State Plots
    figure(4)
    times=TS*((t-1):1:(t-1+N));
    
    subplot(4,1,1);
    OL_1=plot(times,x(1,:)*180/pi,'r--'); hold on
    ylabel('\theta [degree]')
    set(gca,'fontsize', 12)

    subplot(4,1,2);
    OL_2=plot(times,x(2,:),'r--'); hold on
    ylabel('\omega [rad/s]')
    set(gca,'fontsize', 12)

    subplot(4,1,3);
    OL_3=plot(times,x(3,:),'r--' ); hold on
    ylabel('h [m]')
    set(gca,'fontsize', 12)

    subplot(4,1,4);
    OL_4=plot(times,x(4,:),'r--' ); hold on
    ylabel('v [m/s]')
    set(gca,'fontsize', 12)
    xlabel('time [s]');
    set(gca,'fontsize', 12)

    figure(5)
    subplot(2,1,1);
    OL_5=plot(times,x(5,:),'r--' ); hold on
    ylabel('x [m]')
    set(gca,'fontsize', 12)

    subplot(2,1,2);
    OL_6=plot(times,x(6,:),'r--' ); hold on
    ylabel('v_x [m/s]')
    xlabel('time [s]');
    set(gca,'fontsize', 12)
        %---------------------------------------------    
        %Input Plots
        times=TS*((t-1):1:(t-2+N));
    figure(6)
    subplot(3,1,1);
    OLU_1=plot(times,u(1,:),'r--'); hold on
    ylabel('F_e [N]')
    set(gca,'fontsize', 12)

    subplot(3,1,2);
    OLU_2=plot(times,u(2,:),'r--'); hold on
    ylabel('F_{th} [N]')
    set(gca,'fontsize', 12)

    subplot(3,1,3);
    OLU_3=plot(times,u(3,:),'r--' ); hold on
    ylabel('F_{tL} [N]')
    set(gca,'fontsize', 12)
    pause(0.01)
    
    %N=N-1; %Receding horizon MPC
end

% Plot Closed Loop
times_2=TS*(0:1:M);
figure(4)
subplot(4,1,1);
CL_1=plot(times_2,xOpt(1,:)*180/pi,'-bo'); hold on
plot(times_2, xmin(1,1)*ones(size(times_2))*180/pi ); hold on
plot(times_2, xmax(1,1)*ones(size(times_2))*180/pi ); hold on
ylabel('\theta [degree]')
% legend([OL_1 CL_1],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)

subplot(4,1,2);
CL_2=plot(times_2,xOpt(2,:),'-bo'); hold on
plot(times_2, xmin(2,1)*ones(size(times_2)) ); hold on
plot(times_2, xmax(2,1)*ones(size(times_2)) ); hold on
ylabel('\omega [rad/s]')
% legend([OL_2 CL_2],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)

subplot(4,1,3);
CL_3=plot(times_2,xOpt(3,:),'-bo' ); hold on
plot(times_2, xmin(3,1)*ones(size(times_2)) ); hold on
plot(times_2, xmax(3,1)*ones(size(times_2)) ); hold on
ylabel('h [m]')
% legend([OL_3 CL_3],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)

subplot(4,1,4);
CL_4=plot(times_2,xOpt(4,:),'-bo' ); hold on
plot(times_2, xmin(4,1)*ones(size(times_2)) ); hold on
plot(times_2, xmax(4,1)*ones(size(times_2)) ); hold on
ylabel('v [m/s]')
% legend([OL_4 CL_4],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)
title('Closed Loop vs. Open Loop States')

figure(5)
subplot(2,1,1);
CL_5=plot(times_2,xOpt(5,:),'-bo' ); hold on
plot(times_2, xmin(5,1)*ones(size(times_2)) ); hold on
plot(times_2, xmax(5,1)*ones(size(times_2)) ); hold on
ylabel('x [m]')
% legend([OL_5 CL_5],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)
title('Closed Loop vs. Open Loop States')

subplot(2,1,2);
CL_6=plot(times_2,xOpt(6,:),'-bo' ); hold on
plot(times_2, xmin(6,1)*ones(size(times_2)) ); hold on
plot(times_2, xmax(6,1)*ones(size(times_2)) ); hold on
ylabel('v_x [m/s]')
xlabel('time [s]');
set(gca,'fontsize', 12)
% legend([OL_6 CL_6],'OL Predictions','CL Actual');

    %---------------------------------------------    
    %Input Plots
    times_2=TS*(0:1:M-1);
    
figure(6)
subplot(3,1,1);
CLU_1=plot(times_2,uOpt(1,:),'-bo'); hold on
plot(times_2, umin(1,1)*ones(size(times_2))); hold on
plot(times_2, umax(1,1)*ones(size(times_2))); hold on
ylabel('F_e [N]')
% legend([OLU_1 CLU_1],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)
title('Closed Loop vs. Open Loop Inputs')

subplot(3,1,2);
CLU_2=plot(times_2,uOpt(2,:),'-bo'); hold on
plot(times_2, umin(2,1)*ones(size(times_2)) ); hold on
plot(times_2, umax(2,1)*ones(size(times_2)) ); hold on
ylabel('F_{th} [N]')
% legend([OLU_2 CLU_2],'OL Predictions','CL Actual');
set(gca,'fontsize', 12)

subplot(3,1,3);
CLU_3=plot(times_2,uOpt(3,:),'-bo' ); hold on
plot(times_2, umin(3,1)*ones(size(times_2)) ); hold on
plot(times_2, umax(3,1)*ones(size(times_2)) ); hold on
ylabel('F_{tL} [N]')
set(gca,'fontsize', 12)
% legend([OLU_3 CLU_3],'OL Predictions','CL Actual');

figure
hold on
plot(1:M,u_ref(1,1:M))
plot(1:M,uOpt(1,:))
title('Engine Force')
legend('Reference','Actual')

figure
hold on
plot(1:M,u_ref(2,1:M))
plot(1:M,uOpt(2,:))
title('Upper Thruster Force')
legend('Reference','Actual')

figure
hold on
plot(1:M,u_ref(3,1:M))
plot(1:M,uOpt(3,:))
title('Lower Thruster Force')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(1,1:M+1)*180/pi)
plot(1:M+1,xOpt(1,:)*180/pi)
title('Theta')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(2,1:M+1))
plot(1:M+1,xOpt(2,:))
title('Omega w')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(3,1:M+1))
plot(1:M+1,xOpt(3,:))
title('Height')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(4,1:M+1))
plot(1:M+1,xOpt(4,:))
title('Virtical Velocity')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(5,1:M+1))
plot(1:M+1,xOpt(5,:))
title('Horizontal Position')
legend('Reference','Actual')

figure
hold on
plot(1:M+1,x_ref(6,1:M+1))
plot(1:M+1,xOpt(6,1:M+1))
title('Horizontal Velocity')
legend('Reference','Actual')

figure
hold on
plot(1:M,atand(uOpt(2,:)./uOpt(1,:)))
title('Calculated Delta')


save('rocket_data.mat')
 
%end