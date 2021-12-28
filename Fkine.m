%Forward Kinematics

%Link definitions
clc;

F = input("Enter a DH table of form a,alpha,d, theta:\n");

%linkl


T10=DH1(F(1,1),F(1,2),F(1,3),F(1,4))

%link2 

T21=DH1(F(2,1),F(2,2),F(2,3),F(2,4));

T20=T10*T21

%link3 

T32=DH1(F(3,1),F(3,2),F(3,3),F(3,4));

T30 = T10*T21*T32