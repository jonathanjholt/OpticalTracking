% input x y z coordinates (circLocs) as a nx3 matrix
% e.g. circLocs = [-9 15 0;28 6 0;7 -20 0]
% CircFit3D function will then give you centerLoc, circleNormal & radius

clear
clc

p1 = [-110 -55 -60];
p2 = [56 54 0];
p3 = [192 6 20];

circLocs = [p1;p2;p3];

[centerLoc, circleNormal, radius] = CircFit3D(circLocs);

centerLoc