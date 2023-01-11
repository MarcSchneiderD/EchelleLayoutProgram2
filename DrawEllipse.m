function DrawEllipse(app,x1,y1,x2,y2,a,discretization,color)
% x1,y1: 1. focal point
% x2,y2: 2. focal point
% a: major axis radius (axis through focal points)
% discretization: number of points to use for drawing the ellipse
%      (in fact discretization+1, as first and last points are the same)


    %first make an array with angles for full circle with the required
    %discretization
    k=(linspace(0,2*pi,discretization+1))';
    
    focusconnectionvector=[x2, y2]-[x1, y1];
    %line([0, focusconnectionvector(1)], [0, focusconnectionvector(2)]);
    if norm(focusconnectionvector)>=2*a
        warning('on');
        warning('Degenerate Ellipse: major axis too small');
    end
    
    b=sqrt(a^2-(norm(focusconnectionvector)/2)^2);
    
    %alpha=acos( (focusconnectionvector *[1;0]) / (norm(focusconnectionvector)*1) ) %wrong results because no angles <0
    
    % calculate alpha according to
    % https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
    alpha=atan2(focusconnectionvector(2),focusconnectionvector(1));
    
    rotationmatrix=[ cos(alpha) , -sin(alpha) ; sin(alpha) , cos(alpha) ];
    
    %calculate the center of the ellipse
    ecenter=([x1,y1]+[x2,y2])./2;
    
    %create an ellipse with major axis horizontal
    epoints=[a*cos(k) , b*sin(k)];
    
    %now rotate the ellipse to its correct angle
    for i=1:(discretization+1)
        epoints(i,:)=(rotationmatrix*epoints(i,:)')';
    end
    %and shift it to the correct center point
    epoints=epoints+ecenter;
    
    plot(app,epoints(:,1),epoints(:,2),'Color',color,'LineStyle','-');
    axis(app,'equal');
%    hold on
%    plot(app,x1,y1,"ro");
%    plot(app,x2,y2,"*","MarkerEdgeColor",[0 0.5 0]);
%    hold off
end
