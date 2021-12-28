function [A] = DH1(a,alpha,d,theta)
A= [cosd(theta), -sind(theta)*round(cosd(alpha)), sind(theta)*round(sind(alpha)), a*cosd(theta);
    sind(theta), cosd(theta)*round(cosd(alpha)), -cosd(theta)*round(sind(alpha)), a*sind(theta);
    0, round(sind(alpha)), round(cosd(alpha)), d;
    0, 0, 0, 1];
end