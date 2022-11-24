function Pint=CircleIntersectionFinder(x1,y1,r1,x2,y2,r2)
% Find the intersections of two circles
% Output: 2-dimensional array with points in rows, x-coordinate in first column, y-coordinate in second column
% P1 and P2 are column vectors r1 and r2 are their respective radius. % P1 = [x1;y1]; P2 = [x2;y2];


    P1=[x1; y1];
    P2=[x2; y2];
    d2 = sum((P2-P1).^2);
    P0 = (P1+P2)/2+(r1^2-r2^2)/d2/2*(P2-P1);
    t = ((r1+r2)^2-d2)*(d2-(r2-r1)^2);
    if t <= 0
        fprintf('The circles don''t intersect.\n')
        Pa=[NaN; NaN];
        Pb=[NaN; NaN];
    else
        T = sqrt(t)/d2/2*[0 -1;1 0]*(P2-P1);
        Pa = P0 + T; % Pa and Pb are circles' intersection points
        Pb = P0 - T;
    end
    Pint = [Pa, Pb]';

end