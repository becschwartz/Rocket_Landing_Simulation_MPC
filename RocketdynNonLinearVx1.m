%% Demo of recording a moving-line Movie
% This is a simple script to demonstrate some of the ideas to capturing a
% sequence of "frames", collecting them into an array, and finally viewing
% them sequentially as a movie, at a fixed, specified frame rate.

%%
close all
clf
clc
load('rocket_data.mat')

    %states
% x = [theta, w, h, v, x, vx]
%xOpt(1,:)(1 theta = angle to virtical   [rad]
%xOpt(2,:) w = thetaDot                [rad/s]
%xOpt(3,:) h = vertical height         [m]
%xOpt(4,:) v = hDot                    [m/s]   
%xOpt(5,:) x = horizontal position     [m]
%xOpt(6,:) vx = horizontal velocity    [m/s]

    %inputs
% u = [Fe Fth Ftl]
%uOpt(1,:) Fe = Main Engine force      [N]
%uOpt(2,:) Fth = Thruster 1/4 from CM  [N]
%uOpt(3,:) Ftl = Thruster bottom       [N]

% Initial center and orientation of line (uncaptured - see below)
cx = xOpt(5,1); cy = xOpt(3,1)+l/2; theta = xOpt(1,1)+pi/2; 

% Create initial line object (xdata, and ydata)
Lh1 = line([cx-l/2*cos(theta) cx+l/2*cos(theta)],...
    [cy-l/2*sin(theta) cy+l/2*sin(theta)],'color','blue');
nFrames = length(xOpt(1,:));
% Allocate a 1-by-nFrames STRUCT array with 2 fields
F(nFrames) = struct('cdata',[],'colormap',[]);
disp('Creating and recording frames...')
v = VideoWriter('new.avi')
open(v)
for j = 1:nFrames
    % Change center and angle, in a sensible manner, based on Frame#
    cx = xOpt(5,j);
    cy = xOpt(3,j)+l/2;
    theta = xOpt(1,j)+pi/2;
    % Move the line to new location/orientation
    set(Lh1,'xdata',[cx-l/2*cos(theta) cx+l/2*cos(theta)],...
        'ydata', [cy-l/2*sin(theta) cy+l/2*sin(theta)]);
    
    % Make sure the axis stays fixed (and square)
    % axis([-750 750 0 1500]); axis square
    axis([-100 100 0 1500]); axis square
    % Flush the graphics buffer to ensure the line is moved on screen
    drawnow
    % Capture frame
    F(j) = getframe;
    writeVideo(v,F(j));
    writeVideo(v,F(j));
    writeVideo(v,F(j));
end
disp('Playing movie...')
Fps = 10;
nPlay = 1;
movie(F,nPlay,Fps)
close(v)

%% Attribution
% ME C231A and EECS C220B, Fall 2017
