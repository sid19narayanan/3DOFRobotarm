function dxdt = EOM3(t,x)
%mass of links
    m1 = 110; %kg
    m2 = 100; %kg
    m3 = 120; %kg
%link para    
    d1 = 1; %link1
    b  = 0.5; %cglink2
    a2 = 1;  %link2
    a3 = 1;  %link3
    g  = 0.5;  %cglink3
%moment of inertia   
    I1 = 113; %kgm^2
    I2 = 112; %kgm^2
    I3 = 100; %kgm^2
    I_yy1 = 35; 
    I_xx2 = 40; 
    I_yy2 = 40;
    I_zz2 = 42;
    I_xx3 = 30;
    I_yy3 = 20;
    I_zz3 = 50;
    lamda = 130;
    gamma1 = 5250;
    gamma2 = 30;
    gamma3 =20;
    gamma4 =50;
%State variables    
    theta = x(1);
    theta_dot = x(4);
    alpha = x(2);
    alpha_dot = x(5);
    beta = x(3);
    beta_dot = x(6);
    m3_hat = x(7);
    I_xx3hat = x(8);
    I_yy3hat = x(9);
    I_zz3hat = x(10);
%Mass matrix
    M11 = m2*(b^2)*(cos(alpha))^2 + m3*(a2*cos(alpha)-g*cos(alpha+beta))^2 ...
        + I_yy1 + I_xx2*(sin(alpha))^2 ...
        + I_yy2*(cos(alpha))^2 + I_xx3*(sin(alpha+beta))^2+I_yy3*(cos(alpha+beta))^2;
    M12 = 0;
    M13 = 0;
    M21 = 0;
    M22 = m2*b^2 + I_zz2 + m3*a2*a2 + m3*g*g + I_zz3 +2*m3*a2*g*cos(beta);
    M23 = m3*(a2*g*cos(beta)+g^2);
    M31 = 0;
    M32= m3*(a2*g*cos(beta)+g^2);
    M33= m3*g*g+I_zz3;
    M = [M11 M12 M13;M21 M22 M23; M31 M32 M33];
%Coriolis matrix    
    C11 = alpha_dot*cos(alpha)*sin(alpha)*(I_xx2-I_yy2-m2*b^2-m3*a2^2) + ...
            (alpha_dot+beta_dot)*cos(alpha+beta)*(sin(alpha+beta))*(I_xx3-I_yy3-m3*g*g)...
            +a2*g*m3*sin(alpha+beta)*cos(alpha)*(alpha_dot+beta_dot)+a2*g*m3*sin(alpha)*cos(alpha+beta)*alpha_dot;
    C12 = theta_dot*cos(alpha)*sin(alpha)*(I_xx2-I_yy2-m2*b^2-m3*a2^2)...
          + sin(alpha+beta)*cos(alpha+beta)*theta_dot*(I_xx3-I_yy3-m3*g*g) ...
          + theta_dot*(m3*a2*g*sin(alpha)*(cos(alpha+beta))+ m3*a2*g*cos(alpha)*sin(alpha+beta));
    C13 = sin(alpha+beta)*cos(alpha+beta)*theta_dot*(I_xx3-I_yy3-m3*g*g)+ theta_dot*(a2*g*m3*sin(alpha+beta)*cos(alpha));
      
    C21 = -1*C12;
    C22 = -m3*a2*g*sin(beta)*beta_dot;
    C23 = -m3*a2*g*sin(beta)*alpha_dot - a2*g*sin(beta)*beta_dot*m3;
    
    C31 = -1*C13;
    C32 = m3*a2*g*sin(beta)*alpha_dot;
    C33 = 0;
    C = [C11 C12 C13;C21 C22 C23; C31 C32 C33];
%Potential energy matrix        
    G11 = 0;
    G21 = m2*9.81*g*a2*cos(alpha)+m3*g*a2*cos(alpha)+m3*g*a3*cos(alpha+beta);
    G31 = m3*g*a3*cos(alpha+beta);
    G = [G11; G21; G31];
%Trajectories
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
    

    
   theta_des = t1(1)+ t1(2)*mod(t,5) + t1(3)*(mod(t,5)^2) + t1(4)*(mod(t,5)^3) +...
        t1(5)*(mod(t,5)^4)+ t1(6)*(mod(t,5)^5) + t1(7)*(mod(t,5)^6)+ t1(8)*(mod(t,5)^7) + t1(9)*(mod(t,5)^8) +...
        t1(10)*(mod(t,5)^9); 
    theta_dot_des = t1(2)+ 2*t1(3)*(mod(t,5)^1) + 3*t1(4)*(mod(t,5)^2) + 4*t1(5)*(mod(t,5)^3) +...
        5*t1(6)*(mod(t,5)^4) + 6*t1(7)*(mod(t,5)^5) + 7*t1(8)*(mod(t,5)^6)+8*t1(9)*(mod(t,5)^7) +...
        9*t1(10)*(mod(t,5)^8) ;
    theta_dot_dot_des = 2*t1(3)+6*t1(4)*(mod(t,5)^1) + 12*t1(5)*(mod(t,5)^2) + 20*t1(6)*(mod(t,5)^3)...
        + 30*t1(7)*(mod(t,5)^4) + 42*t1(8)*(mod(t,5)^5) + 56*t1(9)*(mod(t,5)^6) + 72*t1(10)*(mod(t,5)^7);
    
    theta_dot_ref = theta_dot_des - lamda*(theta-theta_des);
    theta_dot_dot_ref = theta_dot_dot_des - lamda*(theta_dot - theta_dot_des);
    
    alpha_des = t2(1)+ t2(2)*mod(t,5) + t2(3)*(mod(t,5)^2) + t2(4)*(mod(t,5)^3) +...
        t2(5)*(mod(t,5)^4)+ t2(6)*(mod(t,5)^5) + t2(7)*(mod(t,5)^6)+ t2(8)*(mod(t,5)^7) + t2(9)*(mod(t,5)^8) +...
        t2(10)*(mod(t,5)^9);
    alpha_dot_des = t2(2)+ 2*t2(3)*(mod(t,5)^1) + 3*t2(4)*(mod(t,5)^2) +...
        4*t2(5)*(mod(t,5)^3) + 5*t2(6)*(mod(t,5)^4) + 6*t2(7)*(mod(t,5)^5) + 7*t2(8)*(mod(t,5)^6)+ 8*t2(9)*(mod(t,5)^7) +...
        9*t2(10)*(mod(t,5)^8) ;
    alpha_dot_dot_des = 2*t2(3)+6*t2(4)*(mod(t,5)^1) + 12*t2(5)*(mod(t,5)^2) +...
        20*t2(6)*(mod(t,5)^3) + 30*t2(7)*(mod(t,5)^4) + 42*t2(8)*(mod(t,5)^5) +...
        56*t2(9)*(mod(t,5)^6) + 72*t2(10)*(mod(t,5)^7);
    alpha_dot_ref = alpha_dot_des - lamda*(alpha-alpha_des);
    alpha_dot_dot_ref = alpha_dot_dot_des - lamda*(alpha_dot - alpha_dot_des);
    
    beta_des = t3(1)+ t3(2)*mod(t,5) + t3(3)*(mod(t,5)^2) + t3(4)*(mod(t,5)^3) + t3(5)*(mod(t,5)^4)+ ...
        t3(6)*(mod(t,5)^5) + t3(7)*(mod(t,5)^6)+ t3(8)*(mod(t,5)^7)+ t3(9)*(mod(t,5)^8) +...
        t3(10)*(mod(t,5)^9);
    beta_dot_des = t3(2)+ 2*t3(3)*(mod(t,5)^1) + 3*t3(4)*(mod(t,5)^2) + 4*t3(5)*(mod(t,5)^3) +...
        5*t3(6)*(mod(t,5)^4) + 6*t3(7)*(mod(t,5)^5) + 7*t3(8)*(mod(t,5)^6)+ 8*t3(9)*(mod(t,5)^7) +...
        9*t3(10)*(mod(t,5)^8) ;
    beta_dot_dot_des = 2*t3(3)+6*t3(4)*(mod(t,5)^1) + 12*t3(5)*(mod(t,5)^2) +...
        20*t3(6)*(mod(t,5)^3) + 30*t3(7)*(mod(t,5)^4) + 42*t3(8)*(mod(t,5)^5)+...
        56*t3(9)*(mod(t,5)^6) + 72*t3(10)*(mod(t,5)^7);
    
    beta_dot_ref = beta_dot_des - lamda*(beta-beta_des);
    beta_dot_dot_ref = beta_dot_dot_des - lamda*(beta_dot - beta_dot_des);
    
    %Adaptive controlpart
    
    Y_0 = [m2*b^2*(cos(alpha))^2+I_yy1+ I_xx2*(sin(alpha))^2+ I_yy2*(cos(alpha))^2 0 0; 0 m2*b*b+I_zz2 0; 0 0 0]*...
        [theta_dot_dot_ref; alpha_dot_dot_ref; beta_dot_dot_ref] + ...
        [alpha_dot*cos(alpha)*sin(alpha)*(I_xx2-I_yy2-m2*b^2) theta_dot*cos(alpha)*sin(alpha)*(I_xx2-I_yy2-m2*b^2) 0;...
        -theta_dot*cos(alpha)*sin(alpha)*(I_xx2-I_yy2-m2*b^2) 0 0; 0 0 0]*[theta_dot_ref; alpha_dot_ref; beta_dot_ref]...
        +[0; m2*9.81*a2*cos(alpha); 0];
    
    Y_1 = [(a2*cos(alpha)-g*cos(alpha+beta))^2 0 0; 0 a2^2+g^2+2*a2*g*cos(beta) a2*g*cos(beta)+g^2;...
        0 a2*g*cos(beta)+g^2 g^2]*[theta_dot_dot_ref; alpha_dot_dot_ref; beta_dot_dot_ref] +...
        [-alpha_dot*cos(alpha)*sin(alpha)*a2^2-(alpha_dot+beta_dot)*cos(alpha+beta)*...
        (sin(alpha+beta))+a2*g*sin(alpha+beta)*cos(alpha)*(alpha_dot+beta_dot)+a2*g*sin(alpha)*cos(alpha+beta)*alpha_dot,...
        sin(alpha)*cos(alpha)*theta_dot*(-a2^2)+sin(alpha+beta)*cos(alpha+beta)*theta_dot*(-g^2)+(a2*g*sin(alpha)*...
        cos(alpha+beta)+a2*g*cos(alpha)*sin(alpha+beta))*theta_dot, -sin(alpha+beta)*cos(alpha+beta)*theta_dot*g^2+theta_dot*a2*g*sin(alpha+beta)*cos(alpha);
        sin(alpha)*cos(alpha)*theta_dot*(a2^2)+sin(alpha+beta)*cos(alpha+beta)*theta_dot*(g^2)-(a2*g*sin(alpha)*...
        cos(alpha+beta)-a2*g*cos(alpha)*sin(alpha+beta))*theta_dot, a2*g*sin(beta)*beta_dot, -a2*g*sin(beta)*alpha_dot-a2*g*sin(beta)*beta_dot;...
        sin(alpha+beta)*cos(alpha+beta)*theta_dot*g^2-theta_dot*a2*g*sin(alpha+beta)*cos(alpha), a2*g*sin(beta)*alpha_dot, 0]*...
        [theta_dot_ref; alpha_dot_ref; beta_dot_ref] + [0; 9.81*a2*cos(alpha)+9.81*a3*cos(alpha+beta); 9.81*a3*cos(alpha+beta)];
    
    Y_2 = [(sin(alpha+beta))^2 0 0;0 0 0;0 0 0]*[theta_dot_dot_ref; alpha_dot_dot_ref; beta_dot_dot_ref] + ...
        [cos(alpha+beta)*sin(alpha+beta)*(alpha_dot+beta_dot), sin(alpha+beta)*cos(alpha+beta)*theta_dot,...
        sin(alpha+beta)*cos(alpha+beta)*theta_dot; -sin(alpha+beta)*cos(alpha+beta)*theta_dot 0 0;...
        -sin(alpha+beta)*cos(alpha+beta)*theta_dot 0 0]*[theta_dot_ref; alpha_dot_ref; beta_dot_ref];
    
    Y_3 = [(cos(alpha+beta))^2 0 0;0 0 0;0 0 0]*[theta_dot_dot_ref; alpha_dot_dot_ref; beta_dot_dot_ref] + ...
        [-cos(alpha+beta)*sin(alpha+beta)*(alpha_dot+beta_dot), -sin(alpha+beta)*cos(alpha+beta)*theta_dot,...
        -sin(alpha+beta)*cos(alpha+beta)*theta_dot; sin(alpha+beta)*cos(alpha+beta)*theta_dot 0 0;...
        sin(alpha+beta)*cos(alpha+beta)*theta_dot 0 0]*[theta_dot_ref; alpha_dot_ref; beta_dot_ref];
    
    Y_4 = [0 0 0;0 1 0;0 0 1]*[theta_dot_dot_ref; alpha_dot_dot_ref; beta_dot_dot_ref];
    
    S= [theta_dot-theta_dot_ref;alpha_dot-alpha_dot_ref; beta_dot-beta_dot_ref];
    
    m3_hatdot = -(transpose(S)*Y_1)/gamma1;
    I_xx3hatdot = -(transpose(S)*Y_2)/gamma2;
    I_yy3hatdot = -(transpose(S)*Y_3)/gamma3;
    I_zz3hatdot = -(transpose(S)*Y_4)/gamma4;   
    
    
    
   % Kp = 2;   
    Kd = 120;

    T_in = Y_0 + Y_1*m3_hat + Y_2*I_xx3hat + Y_3*I_yy3hat + Y_4*I_zz3hat - Kd.*S;

%    T_in = [10; 10; 10];
%    acc = inv(M) * (T_in - C * [theta_dot ; alpha_dot; beta_dot] - G);

    acc = inv(M) * (T_in - C * [theta_dot ; alpha_dot; beta_dot] - G);
    dxdt = [theta_dot ; alpha_dot; beta_dot; acc(1); acc(2);acc(3); m3_hatdot; I_xx3hatdot; I_yy3hatdot; I_zz3hatdot];
end