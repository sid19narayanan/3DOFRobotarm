t_space = [0 30];

Mt = [1 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 1 1 1 1; 1 2 4 8 16 32 64 128 256 512; 1 3 9 27 81 243 729 2187 6561 19683;...
        1 4 16 64 256 1024 4096 16384 65536 262144; 1 5 25 125 625 3125 15625 78125 390625 1953125;...
        0 1 0 0 0 0 0 0 0 0; 0 1 2 3 4 5 6 7 8 9; 0 1 4 12 32 80 192 448 1024 2304; 0 1 6 27 108 405 1458 5103 17496 59049;...
        0 1 8 48 256 1280 6144 28672 131072 589824; 0 1 10 75 500 3125 18750 109375 625000 3515625;]; %Matrix for trajec finding 3sec
    c1 = [0; pi/4; pi/2; pi/2; pi/4; 0; 0; 1; 0; 0; -1; 0]; %theta_posand_velo
    c2 = [3*pi/4; 2*pi/3; 3*pi/4; 3*pi/4; 2*pi/3 ; 3*pi/4; 0; 0.1; 0; 0; -0.1; -0.01]; %alpha_posand velo
    c3 = [5*pi/12; pi/3; 5*pi/12; 5*pi/12;pi/3 ; 5*pi/12; 0; 0.1; 0; 0 ; -0.1; -0.01]; %beta_posandvelo
    
    
t1 = Mt\c1;
t2 = Mt\c2;
t3 = Mt\c3;

t0 = 0;
theta_des = t1(1); 
theta_dot_des = t1(2);
alpha_des = t2(1); 
alpha_dot_des = t2(2);
beta_des = t3(1); 
beta_dot_des = t3(2);

q_0 = [theta_des; alpha_des; beta_des; theta_dot_des; alpha_dot_des; beta_dot_des]+[0;0.1;1;0.2;0.2;0.2];
[t,xact] = ode45(@EOM2,t_space,q_0);

theta_d = zeros(size(t,1),1);
alpha_d = zeros(size(t,1),1);
beta_d  = zeros(size(t,1),1);
theta_dot_d = zeros(size(t,1),1);
alpha_dot_d = zeros(size(t,1),1);
beta_dot_d = zeros(size(t,1),1);


for i = 1:size(t,1)

    theta_d(i) = t1(1)+ t1(2)*mod(t(i),5) + t1(3)*(mod(t(i),5)^2) + t1(4)*(mod(t(i),5)^3) + t1(5)*(mod(t(i),5)^4)+ t1(6)*(mod(t(i),5)^5) + t1(7)*(mod(t(i),5)^6)+ t1(8)*(mod(t(i),5)^7)+ t1(9)*(mod(t(i),5)^8)+t1(10)*(mod(t(i),5)^9);
    theta_dot_d(i) = t1(2)+ 2*t1(3)*(mod(t(i),5)) + 3*t1(4)*(mod(t(i),5))^2 + 4*t1(5)*(mod(t(i),5))^3 + 5*t1(6)*(mod(t(i),5))^4 + 6*t1(7)*(mod(t(i),5))^5 + 7*t1(8)*(mod(t(i),5))^6+ 8*t1(9)*(mod(t(i),5))^7+ 9*t1(10)*(mod(t(i),5))^8;
    alpha_d(i) = t2(1)+ t2(2)*(mod(t(i),5)) + t2(3)*(mod(t(i),5))^2 + t2(4)*(mod(t(i),5))^3 + t2(5)*(mod(t(i),5))^4+ t2(6)*(mod(t(i),5))^5 + t2(7)*(mod(t(i),5))^6 + t2(8)*(mod(t(i),5))^7+ t2(9)*(mod(t(i),5)^8)+t2(10)*(mod(t(i),5)^9);
    alpha_dot_d(i) = t2(2)+ 2*t2(3)*(mod(t(i),5)) + 3*t2(4)*(mod(t(i),5))^2 + 4*t2(5)*(mod(t(i),5))^3 + 5*t2(6)*(mod(t(i),5))^4 + 6*t2(7)*(mod(t(i),5))^5 + 7*t2(8)*(mod(t(i),5))^6+ 8*t2(9)*(mod(t(i),5))^7+ 9*t2(10)*(mod(t(i),5))^8;
    beta_d(i)= t3(1)+ t3(2)*(mod(t(i),5)) + t3(3)*(mod(t(i),5))^2 + t3(4)*(mod(t(i),5))^3 + t3(5)*(mod(t(i),5))^4+ t3(6)*(mod(t(i),5))^5 + t3(7)*(mod(t(i),5))^6 + t3(8)*(mod(t(i),5))^7+ t3(9)*(mod(t(i),5)^8)+t3(10)*(mod(t(i),5)^9);
    beta_dot_d(i)= t3(2)+ 2*t3(3)*(mod(t(i),5)) + 3*t3(4)*(mod(t(i),5))^2 + 4*t3(5)*(mod(t(i),5))^3 + 5*t3(6)*(mod(t(i),5))^4 + 6*t3(7)*(mod(t(i),5))^5 + 7*t3(8)*(mod(t(i),5))^6+ 8*t3(9)*(mod(t(i),5))^7+ 9*t3(10)*(mod(t(i),5))^8;
           
end
 
figure(1)
subplot(3,1,1);
plot(t,(xact(:,1)),'LineWidth', 2);
hold on
plot(t,theta_d,'--r');
ylim([0 2]);
xlabel('t(s)');
ylabel('theta(rad)');
title('Computed torque');
legend('Actual OP','Desired')


subplot(3,1,2)
plot(t,(xact(:,2)),'LineWidth', 2);
hold on
plot(t,alpha_d,'--r');
ylim([1.8 3]);
xlabel('t(s)');
ylabel('alpha(rad)');
title('Computed torque');
legend('Actual OP','Desired')


subplot(3,1,3)
plot(t,(xact(:,3)),'LineWidth', 2);
hold on
plot(t,beta_d,'--r');
ylim([0 2]);
xlabel('t(s)');
ylabel('beta(rad)');
title('Computed torque');
legend('Actual OP','Desired')

figure(2)
subplot(3,1,1)
plot(t,(xact(:,4)), 'LineWidth', 2);
hold on
plot(t,theta_dot_d,'--r');
ylim([-2 2]);
xlabel('t(s)');
ylabel('theta dot(rad/s)');
title('Computed torque');
legend('Actual OP','Desired')

subplot(3,1,2)
plot(t,(xact(:,5)), 'LineWidth', 2);
hold on
plot(t,alpha_dot_d,'--r');
xlabel('t(s)');
ylabel('alpha dot(rad/s)');
title('Computed torque');
legend('Actual OP','Desired')

subplot(3,1,3)
plot(t,(xact(:,6)), 'LineWidth', 2);
hold on
plot(t,beta_dot_d,'--r');
ylim([-2 2]);
xlabel('t(s)');
ylabel('beta dot(rad/s)');
title('Computed torque');
legend('Actual OP','Desired')