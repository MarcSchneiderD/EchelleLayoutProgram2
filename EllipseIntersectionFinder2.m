
function EllipseIntersections = EllipseIntersectionFinder(app, x1, y1, x2, y2, x3, y3, a1, a2, IntersectionFinderCenterAngle, IntersectionFinderOpeningAngle, color)
% Calculates the intersections of two ellipses with a common focus point at
% (x1,y1).
% The second focus point of the second ellipse with large half
% axis a1 is at (x2,y2). 
% The second focus point of the third ellipse with large half
% axis a2 is at (x2,y2).
% IntersectionFinderCenterAngle, IntersectionFinderOpeningAngle: defines the wedge in which an intersection
%      is searched for output, angles in degrees
% color: color of intersection marker in diagram
% 
% (Output: Two dimensional Array with intersection points:
% column 1: x-position, column 2: y-position)
% rows: intersections from 0° counterclockwise -> 90° -> 180° -> 270° -> 0°
% at 0° is a calculation gap
% Output: Intersection point [x, y, r, phi]   x,y absolute coordinates, r,phi polar coordinates relative to x1,y1
% If there are several intersections in the search wedge, the one with the
% largest angle is choosen
% Ellipse formulas in polar coordinates from Wikipedia
%

%    clf
    axis(app,'equal');
%    hold on


%    x1=10;
%    y1=0;
%    x2=20;
%    y2=0;
%    x3=21;
%    y3=2;
%    a1=10;
%    a2=9;
    discretization=120;  %discretization for intersection finding before refinement
    intersectionAngleAccuracy=1e-15;    % angle accuracy for calculating the intersections
    
    
    k=(linspace(0,2*pi,discretization+1))';

    focusconnectionvector1=[x2, y2]-[x1, y1];
    focusconnectionvector2=[x3, y3]-[x1, y1];
    
    %line([0, focusconnectionvector(1)], [0, focusconnectionvector(2)]);
    if norm(focusconnectionvector1)>=2*a1
        warning('on');
        warning('Degenerate Ellipse: major axis too small');
    end
    
    if norm(focusconnectionvector2)>=2*a2
        warning('on');
        warning('Degenerate Ellipse: major axis too small');
    end
    
    b1=sqrt(a1^2-(norm(focusconnectionvector1)/2)^2);
%    p1=b1^2/a1;    %Halbparameter
%    e1=sqrt(a1^2-b1^2);
    epsilon1=sqrt(1-(b1/a1)^2);    %Eccentricity
    
    b2=sqrt(a2^2-(norm(focusconnectionvector2)/2)^2);
%    p2=b2^2/a2;    %Halbparameter
%    e2=sqrt(a2^2-b2^2);
    epsilon2=sqrt(1-(b2/a2)^2);    %Eccentricity
    
    
    
    
    % calculate theta according to
    % https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
    theta1=atan2(focusconnectionvector1(2),focusconnectionvector1(1));
    theta2=atan2(focusconnectionvector2(2),focusconnectionvector2(1));

    epointspolar=zeros(discretization+1,3); % Angle, Radius 1. ellipse, Radius 2. Ellipse
    epoints1=zeros(discretization+1,2);
    epoints2=zeros(discretization+1,2);
    
%    testpoints=zeros(discretization+1,2);
    
    
    epointspolar(:,1)=k;
    
    %calculate some points on the ellipse for visualization
    for i=1:(discretization+1)
        phi=epointspolar(i,1);
        epointspolar(i,2)=(a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-phi));
        epointspolar(i,3)=(a2*(1-epsilon2^2))/(1-epsilon2*cos(theta2-phi));
    end
    %epointspolar

    %find all intersections first
    lowPhi=0;%-pi/2;
    lowPhiR1=(a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-lowPhi));
    lowPhiR2=(a2*(1-epsilon2^2))/(1-epsilon2*cos(theta2-lowPhi));
    RdiffLowPhi=lowPhiR2-lowPhiR1;
    intersectionList=[];
    if sign(RdiffLowPhi)==0 %first intersection found!
        intersectionList=[lowPhi];
    end
    for i=1:(discretization+1-1)
        highPhi=i*2*pi/(discretization);%-pi/2;   %break at 0°, shouldn't matter 

        highPhiR1=(a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-highPhi));
        highPhiR2=(a2*(1-epsilon2^2))/(1-epsilon2*cos(theta2-highPhi));
        
        RdiffHighPhi=highPhiR2-highPhiR1;
        
        if (sign(RdiffHighPhi)==0)   %Volltreffer: Intersection at highPhi
            intersectionList=[intersectionList; highPhi];
        elseif (sign(RdiffLowPhi)==0)
            %do nothing, because case already handled before
        elseif (sign(RdiffHighPhi)~=sign(RdiffLowPhi))
            %intersection somewhere between lowPhi and highPhi ==> find it
            %intersectionList=[intersectionList; (lowPhi+highPhi)/2];    %test with rough estimation
            
            hPhi=highPhi;
            lPhi=lowPhi;
            RdiffHPhi=RdiffHighPhi;
            RdiffLPhi=RdiffLowPhi;
            %iPhi=(highPhi+lowPhi)/2;
            %hPhiR1=highPhiR1;
            %hPhiR2=highPhiR2;
            %lPhiR1=lowPhiR1;
            %lPhiR2=lowPhiR2;
            while((hPhi-lPhi>intersectionAngleAccuracy))
                iPhi=(hPhi+lPhi)/2;
                iPhiR1=(a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-iPhi));
                iPhiR2=(a2*(1-epsilon2^2))/(1-epsilon2*cos(theta2-iPhi));
                RdiffIPhi=iPhiR2-iPhiR1;
                if sign(RdiffIPhi)==0
                    % we found the intersection and are finished are finished
                    break;
                elseif (sign(RdiffIPhi)~=sign(RdiffLPhi))
                    % radius difference of intermediate phi has different sign than for low phi
                    % ==> intersection between low phi and intermediate phi
                    % ==> make intermediate phi the new high phi
                    hPhi=iPhi;
                    %hPhiR1=iPhiR1;
                    %hPhiR2=iPhiR2;
                    RdiffHPhi=RdiffIPhi;
                else
                    % intersection between high phi and intermediate phi
                    % ==> make intermediate phi the new low phi
                    lPhi=iPhi;
                    %lPhiR1=iPhiR1;
                    %lPhiR2=iPhiR2;
                    RdiffLPhi=RdiffIPhi;
                end
                
            end
            intersectionList=[intersectionList; (lPhi+hPhi)/2];
            
        end
        
%        testpoints(i,1)=10*cos(highPhi);
%        testpoints(i,2)=10*sin(highPhi);
%        plot(10*cos(highPhi),10*sin(highPhi),'g+');

        lowPhi=highPhi; %the new lowPhi becomes the current highPhi
        lowPhiR1=highPhiR1;
        lowPhiR2=highPhiR2;
        RdiffLowPhi=RdiffHighPhi;
    end
    
    %intersectionList
    
%    plot(testpoints(:,1),testpoints(:,2),'g--');
    
    
%{    
    for i=1:(discretization+1)
        epoints1(i,1)=epointspolar(i,2)*cos(epointspolar(i,1))+x1;
        epoints1(i,2)=epointspolar(i,2)*sin(epointspolar(i,1))+y1;
        epoints2(i,1)=epointspolar(i,3)*cos(epointspolar(i,1))+x1;
        epoints2(i,2)=epointspolar(i,3)*sin(epointspolar(i,1))+y1;
    end
%}
    %epoints1
    %epoints2
    
%    plot(app,epoints1(:,1),epoints1(:,2),'b--');
%    plot(app,epoints2(:,1),epoints2(:,2),'c--');
%    plot(app,x1,y1,"ro");
%    plot(app,x2,y2,"o","MarkerEdgeColor",[0 0.0 1.0]);
%    plot(app,x3,y3,"o","MarkerEdgeColor",[0 1.0 1.0]);

    % now find the last intersection within the search wedge
    IntersectionFinderCenterAngle=deg2rad(IntersectionFinderCenterAngle);
    IntersectionFinderOpeningAngle=deg2rad(IntersectionFinderOpeningAngle);

    FoundIntersectionAngle=NaN;
    if isempty(intersectionList)==false

        for i=1:size(intersectionList,1)
            if (intersectionList(i)>=IntersectionFinderCenterAngle-IntersectionFinderOpeningAngle/2)&&(intersectionList(i)<=IntersectionFinderCenterAngle+IntersectionFinderOpeningAngle/2)
                FoundIntersectionAngle=intersectionList(i);
            end;
        end
    end
    %FoundIntersectionAngle
    IntersectionPhi=rad2deg(FoundIntersectionAngle);
    IntersectionR=(a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-FoundIntersectionAngle));
    IntersectionX=IntersectionR*cos(FoundIntersectionAngle)+x1;
    IntersectionY=IntersectionR*sin(FoundIntersectionAngle)+y1;
    plot(app, IntersectionX , IntersectionY, "d","MarkerEdgeColor",color);
    %{
    intersectionListPlot=[];
    if isempty(intersectionList)==false

        for i=1:size(intersectionList,1)
            intersectionListPlot=[intersectionListPlot; (a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-intersectionList(i)))*cos(intersectionList(i))+x1 , (a1*(1-epsilon1^2))/(1-epsilon1*cos(theta1-intersectionList(i)))*sin(intersectionList(i))+y1];
        end
        plot(app, intersectionListPlot(:,1) , intersectionListPlot(:,2),"d","MarkerEdgeColor",color);
    end
    EllipseIntersections=intersectionListPlot;
    %}
    EllipseIntersections=[IntersectionX, IntersectionY, IntersectionR, IntersectionPhi];
%    hold off
end    
    
    
    
    