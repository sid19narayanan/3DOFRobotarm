% Inverse for Project.
% Takes the DH table and final homogeneous transformation matrix and
% gives the joint values.
close all; clear all; clc;


% DH table for the Project
% Robot dimensions (sample values)

d1=1;
a2=1;
a3=1;

% Joint variables (these will be calculated)
syms theta1 theta2 theta3 theta4 theta5 theta6;
dhtable = [theta1,  d1, 0,  pi/2;
           theta2,  0, a2,  0;
           theta3,  0, a3,  0];
       
% Final homogeneous transformation matrix (sample values)
% These sample values were found using a forward kinematics calculation

Tfinal =[   -0.8660    0.5000         0   -1.5731;
         0         0   -1.0000         0;
   -0.5000   -0.8660         0    1.2071;
         0         0         0    1.0000];
        
% Get joint values
jointValues = AdeptJV(dhtable, Tfinal);
disp('Calculated joint values for the Viper robot manipulator.')
disp(' ')
disp(vpa(jointValues,3))
disp('Units: [radians, radians, radians, radians, radians, radians]')

jointValues = AdeptJV1(dhtable, Tfinal);
disp('Calculated joint values for the Viper robot manipulator.')
disp(' ')
disp(vpa(jointValues,3))
disp('Units: [radians, radians, radians, radians, radians, radians]')

jointValues = AdeptJV2(dhtable, Tfinal);
disp('Calculated joint values for the Viper robot manipulator.')
disp(' ')
disp(vpa(jointValues,3))
disp('Units: [radians, radians, radians, radians, radians, radians]')

jointValues = AdeptJV3(dhtable, Tfinal);
disp('Calculated joint values for the Viper robot manipulator.')
disp(' ')
disp(vpa(jointValues,3))
disp('Units: [radians, radians, radians, radians, radians, radians]')



function [jointVals] = AdeptJV(DHtbl, Tf)

    % Define the elementary rotation and transformation matrices needed
    syms theta d;
    % Rotation homogeneous transformations
    Rx(theta) = [1 0 0 0;0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1];
    Rz(theta) = [cos(theta) -sin(theta) 0 0;sin(theta) cos(theta) 0 0;0 0 1 0;0 0 0 1];
    % Translation homogeneous transformations
    Tx(d) = [1 0 0 d;0 1 0 0;0 0 1 0;0 0 0 1];
    Tz(d) = [1 0 0 0;0 1 0 0;0 0 1 d;0 0 0 1];
    
    % Define the DH-convention composite transformation matrix
    syms thetaz dz ax alphax
    T(thetaz,dz,ax,alphax) = Rz(thetaz)*Tz(dz)*Tx(ax)*Rx(alphax);
    
    % Solve for first three joint values
    
    % Extract some known values from the DH table
    d1 = double(DHtbl(1,2));
    a2 = double(DHtbl(2,3));
    a3 = double(DHtbl(3,3));
    
    % Solve for theta2&3
    
    xyz_bar = double(Tf*[0;0;0;1]);
    x_bar = xyz_bar(1);
    y_bar = xyz_bar(2);
    z_bar = xyz_bar(3);
    L = sqrt(x_bar^2 + y_bar^2);
    
    
    % Solve for theta2
    
    beta = atan2((z_bar-d1),(L));
    gamma = acos((L^2+(z_bar-d1)^2+a2^2-a3^2)/(2*a2*sqrt((L)^2+(z_bar-d1)^2)));
    theta2 = gamma+beta;
   
    
    % Solve for theta3
    theta3 = fsolve(@(theta)[cos(theta); sin(theta)]-(1/(a3^2))*[a3*cos(theta2), a3*sin(theta2); -a3*sin(theta2), a3*cos(theta2)]*[L-a2*cos(theta2);z_bar-a2*sin(theta2)-d1],0);
    
    % Solve for theta1
    theta1 = fsolve(@(theta) [cos(theta); sin(theta)]-(1/(x_bar^2+y_bar^2))*[x_bar, y_bar; y_bar, -x_bar]*[a3*cos(theta3)*cos(theta2)+a2*cos(theta2)-a3*sin(theta2)*sin(theta3); 0], 0);
    
    jointVals = [theta1, theta2, theta3];
end

function [jointVals] = AdeptJV1(DHtbl, Tf)

    % Define the elementary rotation and transformation matrices needed
    syms theta d;
    % Rotation homogeneous transformations
    Rx(theta) = [1 0 0 0;0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1];
    Rz(theta) = [cos(theta) -sin(theta) 0 0;sin(theta) cos(theta) 0 0;0 0 1 0;0 0 0 1];
    % Translation homogeneous transformations
    Tx(d) = [1 0 0 d;0 1 0 0;0 0 1 0;0 0 0 1];
    Tz(d) = [1 0 0 0;0 1 0 0;0 0 1 d;0 0 0 1];
    
    % Define the DH-convention composite transformation matrix
    syms thetaz dz ax alphax
    T(thetaz,dz,ax,alphax) = Rz(thetaz)*Tz(dz)*Tx(ax)*Rx(alphax);
    
    % Solve for first three joint values
    
    % Extract some known values from the DH table
    d1 = double(DHtbl(1,2));
    a2 = double(DHtbl(2,3));
    a3 = double(DHtbl(3,3));
    
    % Solve for theta2&3
    
    xyz_bar = double(Tf*[0;0;0;1]);
    x_bar = xyz_bar(1);
    y_bar = xyz_bar(2);
    z_bar = xyz_bar(3);
    L = sqrt(x_bar^2 + y_bar^2);
    
    
    % Solve for theta2
    
    beta = atan2((z_bar-d1),(L));
    gamma = acos((L^2+(z_bar-d1)^2+a2^2-a3^2)/(2*a2*sqrt((L)^2+(z_bar-d1)^2)));
    theta2 = -gamma+beta;
   
    
    % Solve for theta3
    theta3 = fsolve(@(theta)[cos(theta); sin(theta)]-(1/(a3^2))*[a3*cos(theta2), a3*sin(theta2); -a3*sin(theta2), a3*cos(theta2)]*[L-a2*cos(theta2);z_bar-a2*sin(theta2)-d1],0);
    
    % Solve for theta1
    theta1 = fsolve(@(theta) [cos(theta); sin(theta)]-(1/(x_bar^2+y_bar^2))*[x_bar, y_bar; y_bar, -x_bar]*[a3*cos(theta3)*cos(theta2)+a2*cos(theta2)-a3*sin(theta2)*sin(theta3); 0], 0);
    
    jointVals = [theta1, theta2, theta3];
end

function [jointVals] = AdeptJV2(DHtbl, Tf)

    % Define the elementary rotation and transformation matrices needed
    syms theta d;
    % Rotation homogeneous transformations
    Rx(theta) = [1 0 0 0;0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1];
    Rz(theta) = [cos(theta) -sin(theta) 0 0;sin(theta) cos(theta) 0 0;0 0 1 0;0 0 0 1];
    % Translation homogeneous transformations
    Tx(d) = [1 0 0 d;0 1 0 0;0 0 1 0;0 0 0 1];
    Tz(d) = [1 0 0 0;0 1 0 0;0 0 1 d;0 0 0 1];
    
    % Define the DH-convention composite transformation matrix
    syms thetaz dz ax alphax
    T(thetaz,dz,ax,alphax) = Rz(thetaz)*Tz(dz)*Tx(ax)*Rx(alphax);
    
    % Solve for first three joint values
    
    % Extract some known values from the DH table
    d1 = double(DHtbl(1,2));
    a2 = double(DHtbl(2,3));
    a3 = double(DHtbl(3,3));
    
    % Solve for theta2&3
    
    xyz_bar = double(Tf*[0;0;0;1]);
    x_bar = xyz_bar(1);
    y_bar = xyz_bar(2);
    z_bar = xyz_bar(3);
    L = -sqrt(x_bar^2 + y_bar^2);
    
    
    % Solve for theta2
    
    beta = atan2((z_bar-d1),(L));
    gamma = acos((L^2+(z_bar-d1)^2+a2^2-a3^2)/(2*a2*sqrt((L)^2+(z_bar-d1)^2)));
    theta2 = gamma+beta;
   
    
    % Solve for theta3
    theta3 = fsolve(@(theta)[cos(theta); sin(theta)]-(1/(a3^2))*[a3*cos(theta2), a3*sin(theta2); -a3*sin(theta2), a3*cos(theta2)]*[L-a2*cos(theta2);z_bar-a2*sin(theta2)-d1],0);
    
    % Solve for theta1
    theta1 = fsolve(@(theta) [cos(theta); sin(theta)]-(1/(x_bar^2+y_bar^2))*[x_bar, y_bar; y_bar, -x_bar]*[a3*cos(theta3)*cos(theta2)+a2*cos(theta2)-a3*sin(theta2)*sin(theta3); 0], 0);
    
    jointVals = [theta1, theta2, theta3];
end

function [jointVals] = AdeptJV3(DHtbl, Tf)

    % Define the elementary rotation and transformation matrices needed
    syms theta d;
    % Rotation homogeneous transformations
    Rx(theta) = [1 0 0 0;0 cos(theta) -sin(theta) 0;0 sin(theta) cos(theta) 0;0 0 0 1];
    Rz(theta) = [cos(theta) -sin(theta) 0 0;sin(theta) cos(theta) 0 0;0 0 1 0;0 0 0 1];
    % Translation homogeneous transformations
    Tx(d) = [1 0 0 d;0 1 0 0;0 0 1 0;0 0 0 1];
    Tz(d) = [1 0 0 0;0 1 0 0;0 0 1 d;0 0 0 1];
    
    % Define the DH-convention composite transformation matrix
    syms thetaz dz ax alphax
    T(thetaz,dz,ax,alphax) = Rz(thetaz)*Tz(dz)*Tx(ax)*Rx(alphax);
    
    % Solve for first three joint values
    
    % Extract some known values from the DH table
    d1 = double(DHtbl(1,2));
    a2 = double(DHtbl(2,3));
    a3 = double(DHtbl(3,3));
    
    % Solve for theta2&3
    
    xyz_bar = double(Tf*[0;0;0;1]);
    x_bar = xyz_bar(1);
    y_bar = xyz_bar(2);
    z_bar = xyz_bar(3);
    L = -sqrt(x_bar^2 + y_bar^2);
    
    
    % Solve for theta2
    
    beta = atan2((z_bar-d1),(L));
    gamma = acos((L^2+(z_bar-d1)^2+a2^2-a3^2)/(2*a2*sqrt((L)^2+(z_bar-d1)^2)));
    theta2 = -gamma+beta;
   
    
    % Solve for theta3
    theta3 = fsolve(@(theta)[cos(theta); sin(theta)]-(1/(a3^2))*[a3*cos(theta2), a3*sin(theta2); -a3*sin(theta2), a3*cos(theta2)]*[L-a2*cos(theta2);z_bar-a2*sin(theta2)-d1],0);
    
    % Solve for theta1
    theta1 = fsolve(@(theta) [cos(theta); sin(theta)]-(1/(x_bar^2+y_bar^2))*[x_bar, y_bar; y_bar, -x_bar]*[a3*cos(theta3)*cos(theta2)+a2*cos(theta2)-a3*sin(theta2)*sin(theta3); 0], 0);
    
    jointVals = [theta1, theta2, theta3];
end


