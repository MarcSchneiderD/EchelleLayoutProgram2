% COMSOL Multiphysics project generation for Echelle layout program 2
% requires COMSOL Multiphysics and MATLAB started through "COMSOL Multiphysics with MATLAB"


function BuildComsolProject(uifig, pathname, filename, BraggReflectors, PortTapers, PortTrenches, PortTrenchesCon, SiBoundary, EffectiveRefractiveIndex, ComsolMinWavelength, ComsolMaxWavelength, ComsolWavelengthStep, ShiftRotate, CommentString, NoOfPartitions, ComsolSBends, WaveguideWidth, TrenchWidth, Trench2Width, WaveguideLength, WaveguidePitch, PortEndPositions, PortDirections, Discretization, CladdingBlockMode, BlockModeTrenchwidth)
    import com.comsol.model.*;
    import com.comsol.model.util.*;
    
%    NoOfPartitions=96;  % no of partitions for more efficient meshing of SiBoundary

    % Variable preconditioning
    AbsoluteMinimumRadius=1e100;

    function BorderCoordinates = BezierBend(StartPos, StartDir, EndPos, EndDir, ControlPointDistance, Discretization, Width1, Width2)
    % Calculates cubic bezier curve between start and end point with directions (start point: direction to curve, end point
    % direction out of curve), Control point construction through distance parameter.

        DirS=[cos(StartDir), sin(StartDir)];    %direction vector Start
        DirSP=[-DirS(2), DirS(1)];              %perpendicular direction vector, rotated counterclockwise

        DirE=[cos(EndDir), sin(EndDir)];        %direction vector End
        DirEP=[-DirE(2), DirE(1)];              %perpendicular direction vector, rotated counterclockwise

        ConVec=EndPos-StartPos;   % connection vector between Start and End
        ConVecP=[-ConVec(2), ConVec(1)]/norm(ConVec);   % direction vector perpendicular to ConVec
        ConVecBS=StartPos+ControlPointDistance*ConVec;  % Point between Start and End to calculate new Bezier point nearer to Start
        ConVecBE=EndPos-ControlPointDistance*ConVec;    % Point between Start and End to calculate new Bezier point nearer to End

        BSPos=LineLineIntersection2(StartPos, DirS, ConVecBS, ConVecP);
        BEPos=LineLineIntersection2(EndPos, DirE, ConVecBE, ConVecP);

        L1=zeros(3,2);
        L2=zeros(2,2);
        Path1=zeros(Discretization+1,2);
        Path2=zeros(Discretization+1,2);

        %Calculate cubic Bezier Curve using the De-Castteljau algorithm
        for BBi=0:Discretization
            BBx=BBi*1/Discretization;
            %calculate first level lines for De-Castteljau algorithm...
            L1(1,:)=StartPos+BBx*(BSPos-StartPos);
            L1(2,:)=BSPos+BBx*(BEPos-BSPos);
            L1(3,:)=BEPos+BBx*(EndPos-BEPos);
            %calculate second level lines for De-Castteljau algorithm...
            L2(1,:)=L1(1,:)+BBx*(L1(2,:)-L1(1,:));
            L2(2,:)=L1(2,:)+BBx*(L1(3,:)-L1(2,:));
            %calculate the point of the Bezier spline
            BPos=L2(1,:)+BBx*(L2(2,:)-L2(1,:));
            BDir=L2(2,:)-L2(1,:);   %Direction vector...
            BDirP=[-BDir(2), BDir(1)]/norm(BDir);  %...and direction vector perpendicular to that, normalized

            %Width=Width1+(Width2-Width1)*BBx;
            Width=Width1+(Width2-Width1)*(BBx * BBx * (3 - 2 * BBx));    % https://en.wikipedia.org/wiki/Smoothstep
            %Width=Width1+(Width2-Width1)*(BBx * BBx * BBx * (BBx * (BBx * 6 - 15) + 10)); %https://en.wikipedia.org/wiki/Smoothstep
            Path1(BBi+1,:)=BPos+Width/2*BDirP;    %Path on one side of Bezier curve
            Path2(BBi+1,:)=BPos-Width/2*BDirP;    % Path on other side of Bezier curve
        end
        
        %calculate minimum bending radius
        minRad=1e100;
        for BBi=1:Discretization
            iPos=LineLineIntersectionB2( Path1(BBi,:) , Path2(BBi,:) , Path1(BBi+1,:) , Path2(BBi+1,:) );
            pPos=(Path1(BBi,:)+Path2(BBi,:)+Path1(BBi+1,:)+Path2(BBi+1,:))/4;
            pRad=norm(pPos-iPos);
            if pRad<minRad
                minRad=pRad;
            end
        end
        %fprintf('Minimum wavguide bending radius: %g µm\n', minRad);
        if (minRad<AbsoluteMinimumRadius)
            AbsoluteMinimumRadius=minRad;
        end
        BorderCoordinates=[Path1;flip(Path2,1)];    % border coordinates for GDS structure with Path1 forth and Path2 back
    end

    function ipoint = LineLineIntersection2(p1, p1dir, p3, p3dir)
        ipoint=NaN(1,2);
        %finds intersection point of two lines
        %aline=p1+s*p1dir
        %bline=p3+t*p3dir
        %first convert the vector form of the lines into the coordinate form
        p2=p1+p1dir;
        p4=p3+p3dir;
        % then use the formula from https://de.wikipedia.org/wiki/Schnittpunkt
        %denom=(p4(2)-p3(2))*(p2(1)-p1(1))-(p2(2)-p1(2))*(p4(1)-p3(1));
        %ps(1)=(p4(1)-p3(1))*(p2(1)*p1(2)-p1(1)*p2(2))-(p2(1)-p1(1))*(p4(1)*p3(2)-p3(1)*p4(2))/denom;
        %ps(2)=(p1(2)-p2(2))*(p4(1)*p3(2)-p3(1)*p4(2))-(p3(2)-p4(2))*(p2(1)*p1(2)-p1(1)*p2(2))/denom;
        %ipoint=ps;
        %and that's wrong. SHIT! Änderung eingereicht...

        % a*x+b*y=c
        % with p1=(x1,y1) and p2=(x2,y2)
        % ==> (x,y)=(x1,y1)+t*(x2-x1, y2-y1)
        % ==> a=y2-y1   b=x1-x2   c=x1*y2-x2*y1
        %  => a1=y2-y1  b1=x1-x2  c1=x1*y2-x2*y1
        % and a2=y4-y3  b2=x3-x4  c2=x3*y4-x4*y3
        % Cramersche Regel (https://de.wikipedia.org/wiki/Schnittpunkt ,
        %                   https://de.wikipedia.org/wiki/Cramersche_Regel)
        % ==> xs=(b2*c1-b1*c2)/(a1*b2-a2*b1)   ys=(a1*c2-a2*c1)/(a1*b2-a2*b1)
        denom=((p2(2)-p1(2))*(p3(1)-p4(1))-(p4(2)-p3(2))*(p1(1)-p2(1)));
        ipoint(1)=((p3(1)-p4(1))*(p1(1)*p2(2)-p2(1)*p1(2))-(p1(1)-p2(1))*(p3(1)*p4(2)-p4(1)*p3(2)))/denom;
        ipoint(2)=((p2(2)-p1(2))*(p3(1)*p4(2)-p4(1)*p3(2))-(p4(2)-p3(2))*(p1(1)*p2(2)-p2(1)*p1(2)))/denom;
    end


    function ipoint = LineLineIntersectionB2(p1, p2, p3, p4)
        ipoint=NaN(1,2);
        %finds intersection point of two lines
        %aline: line through p1 and p2
        %bline: line through p3 and p4

        % a*x+b*y=c
        % with p1=(x1,y1) and p2=(x2,y2)
        % ==> (x,y)=(x1,y1)+t*(x2-x1, y2-y1)
        % ==> a=y2-y1   b=x1-x2   c=x1*y2-x2*y1
        %  => a1=y2-y1  b1=x1-x2  c1=x1*y2-x2*y1
        % and a2=y4-y3  b2=x3-x4  c2=x3*y4-x4*y3
        % Cramersche Regel (https://de.wikipedia.org/wiki/Schnittpunkt ,
        %                   https://de.wikipedia.org/wiki/Cramersche_Regel)
        % ==> xs=(b2*c1-b1*c2)/(a1*b2-a2*b1)   ys=(a1*c2-a2*c1)/(a1*b2-a2*b1)
        denom=((p2(2)-p1(2))*(p3(1)-p4(1))-(p4(2)-p3(2))*(p1(1)-p2(1)));
        ipoint(1)=((p3(1)-p4(1))*(p1(1)*p2(2)-p2(1)*p1(2))-(p1(1)-p2(1))*(p3(1)*p4(2)-p4(1)*p3(2)))/denom;
        ipoint(2)=((p2(2)-p1(2))*(p3(1)*p4(2)-p4(1)*p3(2))-(p4(2)-p3(2))*(p1(1)*p2(2)-p2(1)*p1(2)))/denom;
    end









    %first shift the geometry to origin and then rotate 
    %new origin is given in x=ShiftRotate(1), y=ShiftRotate(2), angle=ShiftRotate(3)
    
    ShiftX=ShiftRotate(1);
    ShiftY=ShiftRotate(2);
    RotationAngle=ShiftRotate(3);
    
    %first shift everything...
    for i=1:size(BraggReflectors,1)
        for j=1:size(BraggReflectors,2)
            BraggReflectors{i,j}(1,:)=BraggReflectors{i,j}(1,:)-ShiftX;
            BraggReflectors{i,j}(2,:)=BraggReflectors{i,j}(2,:)-ShiftY;
        end
    end
    
    for i=1:size(PortTrenches,1)
        PortTrenches{i}(1,:)=PortTrenches{i}(1,:)-ShiftX;
        PortTrenches{i}(2,:)=PortTrenches{i}(2,:)-ShiftY;
    end
    
    PortTrenchesCon(1,:)=PortTrenchesCon(1,:)-ShiftX;
    PortTrenchesCon(2,:)=PortTrenchesCon(2,:)-ShiftY;

    for i=1:size(PortTapers,1)
        PortTapers{i}(1,:)=PortTapers{i}(1,:)-ShiftX;
        PortTapers{i}(2,:)=PortTapers{i}(2,:)-ShiftY;
    end
    
    SiBoundary(1,:)=SiBoundary(1,:)-ShiftX;
    SiBoundary(2,:)=SiBoundary(2,:)-ShiftY;
        
    PortEndPositions(:,1)=PortEndPositions(:,1)-ShiftX;
    PortEndPositions(:,2)=PortEndPositions(:,2)-ShiftY;
    
    %...now rotate...
    RotMat = [cos(RotationAngle) -sin(RotationAngle); sin(RotationAngle) cos(RotationAngle)];   %Rotation matrix
    % rotated vector vR is then vR = v*RotMat;


    for i=1:size(BraggReflectors,1)
        for j=1:size(BraggReflectors,2)
            for k=1:size(BraggReflectors{i,j},2)
                help=[BraggReflectors{i,j}(1,k) BraggReflectors{i,j}(2,k)]*RotMat;
                BraggReflectors{i,j}(1,k)=help(1);
                BraggReflectors{i,j}(2,k)=help(2);
            end
        end
    end

    for i=1:size(PortTrenches,1)
        for k=1:size(PortTrenches{i},2)
            help=[PortTrenches{i}(1,k) PortTrenches{i}(2,k)]*RotMat;
            PortTrenches{i}(1,k)=help(1);
            PortTrenches{i}(2,k)=help(2);
        end
    end
    
    for k=1:size(PortTrenchesCon,2)
        help=[PortTrenchesCon(1,k) PortTrenchesCon(2,k)]*RotMat;
        PortTrenchesCon(1,k)=help(1);
        PortTrenchesCon(2,k)=help(2);
    end
    
    for i=1:size(PortTapers,1)
        for k=1:size(PortTapers{i},2)
            help=[PortTapers{i}(1,k) PortTapers{i}(2,k)]*RotMat;
            PortTapers{i}(1,k)=help(1);
            PortTapers{i}(2,k)=help(2);
        end
    end
    
    for k=1:size(SiBoundary,2)
        help=[SiBoundary(1,k) SiBoundary(2,k)]*RotMat;
        SiBoundary(1,k)=help(1);
        SiBoundary(2,k)=help(2);
    end

    for k=1:size(PortEndPositions,1)
        help=[PortEndPositions(k,1) PortEndPositions(k,2)]*RotMat;
        PortEndPositions(k,1)=help(1);
        PortEndPositions(k,2)=help(2);
        
    end
    
    PortDirections=PortDirections-RotationAngle;
    
    
    
    %find leftmost, rightmost, topmost, and bottommost coordinate of silicon boundary for mesh partitioning
    SiB_top=SiBoundary(2,1);
    SiB_left=SiBoundary(1,1);
    SiB_bottom=SiBoundary(2,1);
    SiB_right=SiBoundary(1,1);
    for k=1:size(SiBoundary,2)
        if SiBoundary(1,k)<SiB_left
            SiB_left=SiBoundary(1,k);
        end
        if SiBoundary(1,k)>SiB_right
            SiB_right=SiBoundary(1,k);
        end
        if SiBoundary(2,k)>SiB_top
            SiB_top=SiBoundary(2,k);
        end
        if SiBoundary(2,k)<SiB_bottom
            SiB_bottom=SiBoundary(2,k);
        end
    end
    

    model = ModelUtil.create('Model');
    comp=model.component.create('comp1');
    geom=model.component('comp1').geom.create('geom', 2);   %create 2-dimensional geometry

    comp.comments(sprintf(CommentString));
    comp.author('Echelle Layout Program v2 by Dr.-Ing. Marc Schneider, Institute for Data Processing and Electronics (IPE), Karlsruhe Institute of Technology (KIT), Germany');
    comp.version('1.0');

    % define global basic parameters
    model.param.set('lambda0', sprintf('%u[nm]',(ComsolMinWavelength+ComsolMaxWavelength)/2));
    model.param.set('lambda_min', sprintf('%u[nm]',ComsolMinWavelength));
    model.param.set('lambda_max', sprintf('%u[nm]',ComsolMaxWavelength));
    model.param.set('lambda_step', sprintf('%u[nm]',ComsolWavelengthStep));
    model.param.set('freq0', 'c_const/lambda0');
    model.param.set('w0', '2*pi*freq0');
    
    model.geom('geom').lengthUnit('µm'); % µm ist Standardeinheit
    model.geom('geom').label('EchelleGeometry');

    d = uiprogressdlg(uifig,'Title','COMSOL project generation','Message','Adding physics...');
    
    %==============
    %adding Physics / this part also checks for licenses
    %==============
    % frequency domain
    comp.physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom');
    comp.physics('ewfd').prop('components').set('components', 'inplane');
    
    
    %==============
    %adding Geometry
    %==============
    d.Value=0.1;
    d.Message = 'Adding free space area...';
    % now add the large free space area, the SiBoundary
    
    mpsitag = sprintf('meshpolygon_si');
    meshpolysi_geom_name = geom.feature.create(mpsitag, 'Polygon');
    meshpolysi_geom_name.set('type', 'solid');
    meshpolysi_geom_name.set('x', SiBoundary(1,:));
    meshpolysi_geom_name.set('y', SiBoundary(2,:));
    meshpolysi_geom_name.label(mpsitag);
    meshpolysi_geom_name.set('selresult', 'on');

    
    if (NoOfPartitions>1)
        %we want to have several partitions, so generate Lines to cut there
        polytags=cell(NoOfPartitions-1,1);
        for i=1:NoOfPartitions-1
            polytag = sprintf('poly_%u',i);
            polytags{i}=polytag;
            poly_geom_name = geom.feature.create(polytag, 'Polygon');
            poly_geom_name.set('type', 'open');
            poly_geom_name.set('source', 'table');
            poly_y=(SiB_top-SiB_bottom)/NoOfPartitions*i+SiB_bottom;
            poly_geom_name.set('table', [SiB_left-10 poly_y; SiB_right+10 poly_y]);
    %        poly_x=(SiB_right-SiB_left)/NoOfPartitions*i+SiB_left;
    %        poly_geom_name.set('table', [poly_x SiB_bottom-10; poly_x SiB_top+10]);

        end


        partition_geom_name = geom.feature.create('partition1', 'Partition');
        partition_geom_name.selection('input').set(mpsitag);
        partition_geom_name.selection('tool').set(polytags);
    end;    

    

    geom.selection.create('csel5', 'CumulativeSelection');
    geom.selection('csel5').label('FreeSpace');
    geom.create('free_space', 'Union');
    if (NoOfPartitions>1)
    %    geom.feature('free_space').selection('input').set('meshpolygon_si');
        geom.feature('free_space').selection('input').set('partition1');
    else
       geom.feature('free_space').selection('input').set('meshpolygon_si');
    end;
    geom.feature('free_space').set('contributeto', 'csel5');
    
    
    d.Value=0.2;
    d.Message = 'Adding reflectors...';
    
    %Build reflector geometries
    rtags=cell(size(BraggReflectors,1)*size(BraggReflectors,2),1);  % full list required for union and selection
    for i=1:size(BraggReflectors,1)
        for j=1:size(BraggReflectors,2)
            rtag = sprintf('refl_%u_%u',i,j);
            rtags{(i-1)*size(BraggReflectors,2)+j}=rtag;
            refl_geom_name = geom.feature.create(rtag, 'Polygon');
            refl_geom_name.set('type', 'solid');
            refl_geom_name.set('x', BraggReflectors{i,j}(1,:));
            refl_geom_name.set('y', BraggReflectors{i,j}(2,:));
            refl_geom_name.label(rtag);
        end
    end
    
    d.Value=0.4;
    d.Message = 'Adding ports...';

    
    ttag = sprintf('WGtrench');
    trench_geom_name = geom.feature.create(ttag, 'Polygon');
    trench_geom_name.set('type', 'solid');
    trench_geom_name.set('x', PortTrenchesCon(1,:));
    trench_geom_name.set('y', PortTrenchesCon(2,:));
    trench_geom_name.label(ttag);

    
    dtags=cell(size(PortTapers,1),1);  % full list required for union and selection
    ptags=cell(size(PortTapers,1),1);  % full list required for union and selection
    for i=1:size(PortTapers,1)
        if i==size(PortTapers,1)
            ptag = sprintf('common_port');
            %ttag = sprintf('common_trench');
            %dtag = sprintf('common_difference');
        else
            ptag = sprintf('port_%u',i);
            %ttag = sprintf('trench_%u',i);
            %dtag = sprintf('difference_%u',i);
        end
        %dtags{i}=dtag;
        ptags{i}=ptag;
        
%        trench_geom_name = geom.feature.create(ttag, 'Polygon');
%        trench_geom_name.set('type', 'solid');
%        trench_geom_name.set('x', PortTrenches{i}(1,:));
%        trench_geom_name.set('y', PortTrenches{i}(2,:));
%        trench_geom_name.label(ttag);
        
        port_geom_name = geom.feature.create(ptag, 'Polygon');
        port_geom_name.set('type', 'solid');
        port_geom_name.set('x', PortTapers{i}(1,:));
        port_geom_name.set('y', PortTapers{i}(2,:));
        port_geom_name.label(ptag);

%        difference_geom_name = geom.feature.create(dtag, 'Difference');
%        difference_geom_name.selection('input').set(ttag);
%        difference_geom_name.selection('input2').set(ptag);
%        difference_geom_name.label(dtag);
%        difference_geom_name.set('keep', 'on');
    end
    
    dtag = sprintf('WGdifference');
    difference_geom_name = geom.feature.create(dtag, 'Difference');
    difference_geom_name.selection('input').set(ttag);
    difference_geom_name.selection('input2').set(ptags);
    difference_geom_name.label(dtag);
    difference_geom_name.set('keep', 'on');
    
    
    
    if (ComsolSBends==1)
        % Build waveguide bundle

        sdtags=cell(size(PortTapers,1),1);  % full list required for union and selection
        swgtags=cell(size(PortTapers,1),1);  % full list required for union and selection
        for i=1:size(PortTapers,1)
            if i==size(PortTapers,1)
                swgtag = sprintf('SBend_common_port');
                sttag = sprintf('SBend_common_trench');
                sdtag = sprintf('SBend_common_difference');
            else
                swgtag = sprintf('SBend_port_%u',i);
                sttag = sprintf('SBend_trench_%u',i);
                sdtag = sprintf('SBend_difference_%u',i);
            end
            sdtags{i}=sdtag;
            swgtags{i}=swgtag;

            swgtrench_geom_name = geom.feature.create(sttag, 'Polygon');
            swgtrench_geom_name.set('type', 'solid');
            if i==size(PortTapers,1)
                SBendTrench=BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3);
            else
                SBendTrench=BezierBend(PortEndPositions(i,:), PortDirections(i), PortEndPositions(end,:)+[i*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3);
            end
            swgtrench_geom_name.set('x', SBendTrench(:,1));
            swgtrench_geom_name.set('y', SBendTrench(:,2));
            swgtrench_geom_name.label(sttag);

            swg_geom_name = geom.feature.create(swgtag, 'Polygon');
            swg_geom_name.set('type', 'solid');
            if i==size(PortTapers,1)
                SBendWG=BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3 , WaveguideWidth/1e3 );
            else
                SBendWG=BezierBend(PortEndPositions(i,:), PortDirections(i), PortEndPositions(end,:)+[i*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3 , WaveguideWidth/1e3 );    
            end
            swg_geom_name.set('x', SBendWG(:,1));
            swg_geom_name.set('y', SBendWG(:,2));
            swg_geom_name.label(swgtag);

            sdifference_geom_name = geom.feature.create(sdtag, 'Difference');
            sdifference_geom_name.selection('input').set(sttag);
            sdifference_geom_name.selection('input2').set(swgtag);
            sdifference_geom_name.label(sdtag);
            sdifference_geom_name.set('keep', 'on');
        end
        fprintf('     Minimum wavguide bending radius: %g µm\n', AbsoluteMinimumRadius);
        
    
    end
    
    
    
    
    
    
    
    
    
    
    
    %Adding Elements to list for easier Comsol selection

    geom.selection.create('csel3', 'CumulativeSelection');
    geom.selection('csel3').label('Reflector_union');
    geom.create('uni_refl', 'Union');
    geom.feature('uni_refl').selection('input').set(rtags);
    geom.feature('uni_refl').set('contributeto', 'csel3');
    
    geom.selection.create('csel4', 'CumulativeSelection');
    geom.selection('csel4').label('Trench_union');
    geom.create('uni_trenches', 'Union');
    if (ComsolSBends==0)
%        geom.feature('uni_trenches').selection('input').set(dtags);
        geom.feature('uni_trenches').selection('input').set(dtag);
        geom.feature('uni_trenches').set('contributeto', 'csel4');
    else        
%        geom.feature('uni_trenches').selection('input').set( cat(1,dtags,sdtags) );
        geom.feature('uni_trenches').selection('input').set( cat(1,dtag,sdtags) );
        geom.feature('uni_trenches').set('contributeto', 'csel4');
    end
   
    geom.selection.create('csel6', 'CumulativeSelection');
    geom.selection('csel6').label('Port_union');
    geom.create('uni_ports', 'Union');
    if (ComsolSBends==0)
        geom.feature('uni_ports').selection('input').set(ptags);
        geom.feature('uni_ports').set('contributeto', 'csel6');
    else
        geom.feature('uni_ports').selection('input').set( cat(1,ptags,swgtags) );
        geom.feature('uni_ports').set('contributeto', 'csel6');
    end

    geom.selection.create('csel7', 'CumulativeSelection');
    geom.selection('csel7').label('PortTrench_union');
    geom.create('uni_portstrenches', 'Union');
    geom.feature('uni_portstrenches').selection('input').set({'uni_trenches', 'uni_ports'});
    geom.feature('uni_portstrenches').set('contributeto', 'csel7');
    
    geom.run
    
    %Boxselection for all Objects
    geom.create('boxsel1', 'BoxSelection');
    geom.feature('boxsel1').set('entitydim', '-1');
    geom.feature('boxsel1').label('AllObjects');

    
    
    %==============
    %adding Material
    %==============
    d.Value=0.50;
    d.Message = 'Adding materials...';

    % silicon material with effective refractive index
    ncore_table=cell(size(EffectiveRefractiveIndex,1),2);
    for i=1:size(EffectiveRefractiveIndex,1)
        ncore_table{i,1} = sprintf('%fE-9',EffectiveRefractiveIndex(i,1));
        ncore_table{i,2} = sprintf('%f',EffectiveRefractiveIndex(i,2));
    end
    comp_mat = comp.material.create('mat1', 'Common');
    comp_mat.label('Si_neff');
    comp_mat_refi = comp_mat.propertyGroup.create('RefractiveIndex', 'Refractive index');
    comp_mat_refi.set('n', 'n_interp(1[1/m]*c_const/freq)');
    comp_mat_refi.func.create('n_interp', 'Interpolation');
    comp_mat_refi.func('n_interp').set('sourcetype', 'user');
    comp_mat_refi.func('n_interp').set('source', 'table');
    comp_mat_refi.func('n_interp').set('funcname', 'n_interp');
    comp_mat_refi.func('n_interp').set('table', ncore_table);
    comp_mat_refi.func('n_interp').set('interp', 'piecewisecubic');
    comp_mat_refi.func('n_interp').set('extrap', 'linear');
    comp_mat_refi.addInput('frequency');
    
    % from COMSOL material library "SiO2 (Malitson)":
    ncladd_table = {'2.0999999999999997E-7' '1.5383576204905378';
                    '2.45E-7' '1.510272436589456';
                    '2.8E-7' '1.4941636611187716';
                    '3.15E-7' '1.4839008951422648';
                    '3.5E-7' '1.476891413495998';
                    '3.8499999999999997E-7' '1.4718556531995413';
                    '4.1999999999999995E-7' '1.4680936900401065';
                    '4.5499999999999993E-7' '1.4651930999599743';
                    '4.9E-7' '1.4628966820387057';
                    '5.25E-7' '1.4610366660573582';
                    '5.599999999999999E-7' '1.4594995356592282';
                    '5.949999999999999E-7' '1.4582061049260293';
                    '6.3E-7' '1.4570996888768784';
                    '6.65E-7' '1.4561387969802702';
                    '7.0E-7' '1.4552924662622837';
                    '7.35E-7' '1.45453719287602';
                    '7.699999999999999E-7' '1.4538548630588606';
                    '8.049999999999999E-7' '1.4532313266004242';
                    '8.399999999999999E-7' '1.4526553936728075';
                    '8.75E-7' '1.4521181167939423';
                    '9.099999999999999E-7' '1.4516122686289965';
                    '9.45E-7' '1.4511319566976737';
                    '9.8E-7' '1.4506723353352597';
                    '1.015E-6' '1.4502293877558508';
                    '1.05E-6' '1.4497997593262841';
                    '1.085E-6' '1.4493806287125588';
                    '1.12E-6' '1.4489696073536897';
                    '1.155E-6' '1.4485646603469264';
                    '1.1899999999999998E-6' '1.4481640436751153';
                    '1.2249999999999997E-6' '1.4477662540206953';
                    '1.26E-6' '1.4473699883562041';
                    '1.2949999999999999E-6' '1.446974111188911';
                    '1.33E-6' '1.4465776278425762';
                    '1.3649999999999998E-6' '1.4461796625342722';
                    '1.4E-6' '1.4457794402848239';
                    '1.435E-6' '1.4453762719132304';
                    '1.47E-6' '1.4449695415265835';
                    '1.5049999999999998E-6' '1.4445586960404744';
                    '1.5399999999999999E-6' '1.4441432363602276';
                    '1.575E-6' '1.4437227099273546';
                    '1.6099999999999998E-6' '1.4432967043935492';
                    '1.645E-6' '1.4428648422301078';
                    '1.6799999999999998E-6' '1.4424267761167082';
                    '1.7149999999999999E-6' '1.4419821849821588';
                    '1.75E-6' '1.4415307705926665';
                    '1.7849999999999999E-6' '1.4410722546016008';
                    '1.82E-6' '1.4406063759896086';
                    '1.8549999999999998E-6' '1.4401328888360165';
                    '1.89E-6' '1.439651560372282';
                    '1.9249999999999998E-6' '1.4391621692763015';
                    '1.96E-6' '1.4386645041729973';
                    '1.995E-6' '1.4381583623120469';
                    '2.03E-6' '1.4376435483981267';
                    '2.0649999999999997E-6' '1.4371198735527937';
                    '2.1E-6' '1.4365871543902402';
                    '2.1350000000000003E-6' '1.4360452121917664';
                    '2.17E-6' '1.435493872166012';
                    '2.205E-6' '1.4349329627838137';
                    '2.2399999999999997E-6' '1.434362315178124';
                    '2.275E-6' '1.4337817626007234';
                    '2.31E-6' '1.4331911399285864';
                    '2.3449999999999996E-6' '1.4325902832137039';
                    '2.3799999999999997E-6' '1.4319790292709682';
                    '2.4149999999999997E-6' '1.4313572152994267';
                    '2.4500000000000003E-6' '1.4307246785327923';
                    '2.485E-6' '1.4300812559156133';
                    '2.52E-6' '1.4294267838019374';
                    '2.555E-6' '1.4287610976736893';
                    '2.5899999999999998E-6' '1.4280840318762884';
                    '2.625E-6' '1.4273954193693414';
                    '2.66E-6' '1.4266950914904597';
                    '2.6949999999999996E-6' '1.4259828777304846';
                    '2.7299999999999997E-6' '1.4252586055185756';
                    '2.765E-6' '1.4245221000157793';
                    '2.8E-6' '1.4237731839158354';
                    '2.835E-6' '1.423011677252103';
                    '2.87E-6' '1.4222373972095803';
                    '2.9049999999999997E-6' '1.421450157941106';
                    '2.94E-6' '1.4206497703868863';
                    '2.975E-6' '1.419836042096581';
                    '3.0099999999999996E-6' '1.419008777053232';
                    '3.0449999999999996E-6' '1.4181677754983741';
                    '3.0799999999999997E-6' '1.4173128337577083';
                    '3.115E-6' '1.4164437440667712';
                    '3.15E-6' '1.41556029439605';
                    '3.185E-6' '1.4146622682750325';
                    '3.2199999999999997E-6' '1.4137494446147074';
                    '3.2549999999999998E-6' '1.4128215975280431';
                    '3.29E-6' '1.4118784961479969';
                    '3.325E-6' '1.410919904442614';
                    '3.3599999999999996E-6' '1.4099455810267978';
                    '3.3949999999999997E-6' '1.4089552789703252';
                    '3.43E-6' '1.4079487456017001';
                    '3.465E-6' '1.4069257223074338';
                    '3.5E-6' '1.4058859443263447';
                    '3.535E-6' '1.404829140538468';
                    '3.5699999999999997E-6' '1.4037550332481645';
                    '3.605E-6' '1.4026633379610098';
                    '3.64E-6' '1.4015537631540418';
                    '3.6749999999999995E-6' '1.4004260100389325';
                    '3.7099999999999996E-6' '1.3992797723176442'};

    comp_mat2 = comp.material.create('mat2', 'Common');
    comp_mat2.label('SiO2');
    comp_mat2_refi = comp_mat2.propertyGroup.create('RefractiveIndex', 'Refractive index');
    comp_mat2_refi.set('n', 'n_interp(1[1/m]*c_const/freq)');
    comp_mat2_refi.func.create('n_interp', 'Interpolation');
    comp_mat2_refi.func('n_interp').set('sourcetype', 'user');
    comp_mat2_refi.func('n_interp').set('source', 'table');
    comp_mat2_refi.func('n_interp').set('funcname', 'n_interp');
    comp_mat2_refi.func('n_interp').set('table', ncladd_table);
    comp_mat2_refi.func('n_interp').set('interp', 'piecewisecubic');
    comp_mat2_refi.func('n_interp').set('extrap', 'linear');
    comp_mat2_refi.addInput('frequency');
    
    geom.run

    % I have NO idea, why you get the entities of the domain in the selections csel3 and csel4 in this way
    % (mphgetselection(model.selection('geom_csel3_dom')).entities) and who has found it out. It's an absolute mystery.
    % Maybe, because the model geometry is named 'geom' and the selections 'csel3' and 'csel4' and the domain
    % abbrevated with 'dom'...
    model_sel = [mphgetselection(model.selection('geom_csel3_dom')).entities, mphgetselection(model.selection('geom_csel4_dom')).entities];
    comp_mat2.selection.set( model_sel );


    
    %==============
    %adding Physics (additional)
    %==============
    d.Value=0.65;
    d.Message = 'Adding additional physics...';

    
    comp.physics('ewfd').create('sctr1', 'Scattering', 1);
    comp.physics('ewfd').feature('sctr1').selection.all;
    comp.physics('ewfd').feature('sctr1').set('Order', 'SecondOrder');

    %comp.physics('ewfd').feature('sctr1').selection.entities(1)
    
    geom.run;

    %adding ports to physics

    for i=1:size(PortTrenches,1)
        if i==size(PortTrenches,1)
            ptag = sprintf('PortCommon');
        else
            ptag = sprintf('Port%u',i);
        end
        
        %PortTrenches{i}(1,:));
        % [<x0> <x1>;<y0> <y1>]
        if (ComsolSBends==0)
            objsel(1,1) = PortTrenches{i}(1,3); % 1. point x
            objsel(2,1) = PortTrenches{i}(2,3); % 1. point y
            objsel(1,2) = PortTrenches{i}(1,4); % 2. point x
            objsel(2,2) = PortTrenches{i}(2,4); % 2. point y
        else
            if i==size(PortTrenches,1)
                objsel(1,1) = PortEndPositions(end,1)-WaveguideWidth/2/1e3-Trench2Width/1e3; % 1. point x
                objsel(2,1) = PortEndPositions(end,2)-WaveguideLength/1e3; % 1. point y
                objsel(1,2) = PortEndPositions(end,1)+WaveguideWidth/2/1e3+Trench2Width/1e3; % 2. point x
                objsel(2,2) = PortEndPositions(end,2)-WaveguideLength/1e3; % 2. point y
            else
                objsel(1,1) = PortEndPositions(end,1)+i*WaveguidePitch/1e3-WaveguideWidth/2/1e3-Trench2Width/1e3; % 1. point x
                objsel(2,1) = PortEndPositions(end,2)-WaveguideLength/1e3; % 1. point y
                objsel(1,2) = PortEndPositions(end,1)+i*WaveguidePitch/1e3+WaveguideWidth/2/1e3+Trench2Width/1e3; % 2. point x
                objsel(2,2) = PortEndPositions(end,2)-WaveguideLength/1e3; % 2. point y
            end
            
            
        end

        SelectionBoxEnhancement=1e-3;  % make selection box bigger in all directions by this in µm
        if objsel(1,1) > objsel(1,2)
            objsel(1,1) = objsel(1,1)+SelectionBoxEnhancement;
            objsel(1,2) = objsel(1,2)-SelectionBoxEnhancement;
        else
            objsel(1,1) = objsel(1,1)-SelectionBoxEnhancement;
            objsel(1,2) = objsel(1,2)+SelectionBoxEnhancement;
        end
        if objsel(2,1) > objsel(2,2)
            objsel(2,1) = objsel(2,1)+SelectionBoxEnhancement;
            objsel(2,2) = objsel(2,2)-SelectionBoxEnhancement;
        else
            objsel(2,1) = objsel(2,1)-SelectionBoxEnhancement;
            objsel(2,2) = objsel(2,2)+SelectionBoxEnhancement;
        end
        
        boundrylabel = mphselectbox(model,'geom',objsel,'boundary');

        model.physics('ewfd').create(ptag, 'Port', 1);
        model.physics('ewfd').feature(ptag).set('PortType', 'Numeric');
        model.physics('ewfd').feature(ptag).label(ptag);
        
        if i==size(PortTrenches,1)
            model.physics('ewfd').feature(ptag).set('PortExcitation', 'on');
        else
            model.physics('ewfd').feature(ptag).set('PortExcitation', 'off');
        end

        if isempty(boundrylabel)==1 % first selection attempt not successful, try it in a different way...
            boundarylabel_1 = mphselectcoords(model,'geom',objsel(:,1)','boundary');
            boundarylabel_2 = mphselectcoords(model,'geom',objsel(:,2)','boundary');
            boundrylabel = intersect(boundarylabel_1,boundarylabel_2);
        end
        if isempty(boundrylabel)~=1 % now make port with first or second selection
            model.physics('ewfd').feature(ptag).selection.set(boundrylabel);
            comp.physics('ewfd').feature('sctr1').selection.remove(model.physics('ewfd').feature(ptag).selection.entities(1));
        end
    end
    
    %==============
    %adding Mesh
    %==============
    d.Value=0.8;
    d.Message = 'Adding mesh...';

    
    comp_mesh = comp.mesh.create('mesh1', 'geom');
    comp_mesh.create('size1', 'Size');

    comp_mesh.feature('size1').set('custom', 'on');
    comp_mesh.feature('size1').set('hmaxactive', 'on');
    comp_mesh.feature('size1').set('hmax', 'lambda_min/3/5');
    comp_mesh.feature('size1').set('hmin', 'lambda_min/3/8');

    comp_mesh.feature('size1').selection.geom('geom', 2);
    comp_mesh.feature('size1').selection.all;
    comp_mesh.create('ftri1', 'FreeTri');
    comp_mesh.feature('ftri1').selection.remaining;


    comp_mesh2 = comp.mesh.create('mesh2', 'geom');
    
    comp_mesh2.create('copy1', 'Copy');
    comp_mesh2.feature('copy1').selection('source').geom(2);
    comp_mesh2.feature('copy1').selection('destination').geom(2);
    %comp_mesh2.feature.move('copy1', 1);
    comp_mesh2.feature('copy1').set('mesh', 'mesh1');
    comp_mesh2.feature('copy1').set('buildsource', true);
    comp_mesh2.feature('copy1').set('dimension', 1);
    comp_mesh2.feature('copy1').selection('source').named('geom_csel7_bnd');
    comp_mesh2.feature('copy1').selection('destination').named('geom_csel7_bnd');

    comp_mesh2.create('size1', 'Size');
    comp_mesh2.feature('size1').set('custom', 'on');
    comp_mesh2.feature('size1').set('hmaxactive', 'on');
    comp_mesh2.feature('size1').set('hmax', 'lambda_min/3/1');
    comp_mesh2.feature('size1').set('hmin', 'lambda_min/3/8');

    comp_mesh2.feature('size1').selection.geom('geom', 2);
    comp_mesh2.feature('size1').selection.all;
    
    comp_mesh2.create('size2', 'Size');

    comp_mesh2.feature('size2').set('custom', 'on');
    comp_mesh2.feature('size2').set('hmaxactive', 'on');
    comp_mesh2.feature('size2').set('hmax', (SiB_top-SiB_bottom)/NoOfPartitions*2);
    comp_mesh2.feature('size2').set('hmin', (SiB_top-SiB_bottom)/NoOfPartitions/10);
    comp_mesh2.feature('size2').set('hgrad', '2');
    comp_mesh2.feature('size2').set('hnarrowactive', true);
    comp_mesh2.feature('size2').set('hnarrow', 0.1);

    comp_mesh2.feature('size2').selection.named('geom_csel3_dom');
    comp_mesh2.feature('size2').selection.named('geom_csel5_dom');
    
    comp_mesh2.create('ftri1', 'FreeTri');
    comp_mesh2.feature('ftri1').selection.remaining;


    
    %==============
    %adding Study
    %==============
    d.Value=0.9;
    d.Message = 'Adding study...';

    stdy=model.study.create('std1');
    
    stdy.create('param', 'Parametric');
    stdy.feature('param').setIndex('pname', 'lambda0', 0);
    stdy.feature('param').setIndex('plistarr', 'range(lambda_min,lambda_step,lambda_max)', 0);
    stdy.feature('param').setIndex('punit', 'nm', 0);
    stdy.feature('param').set('plot', true);
    stdy.feature('param').set('plotgroup', 'Default');


    %Boundary mode analysis nodes here
    for i=1:size(PortTapers,1)
        if i==size(PortTapers,1)
            bma=stdy.create(sprintf('bma'), 'BoundaryModeAnalysis');
        else
            bma=stdy.create(sprintf('bma%u', i), 'BoundaryModeAnalysis');
        end
        
        bma.set('PortName', num2str(i));
        bma.set('appnreigs', '1');

        bma.set('modeFreq', 'freq0');
        bma.set('neigsactive', true);
        bma.set('shiftactive', true);
        bma.set('shift', '3.5');
        bma.set('modeFreq', 'freq0');
        bma.set('physselection', 'ewfd');
        bma.setIndex('mesh', 'mesh2', 1);
        
        
    end    
    
    %Frequency domain study here
    stdy.create('freq', 'Frequency');
    stdy.feature('freq').set('punit', 'Hz');
    stdy.feature('freq').set('plist', 'freq0');
    stdy.feature('freq').setIndex('mesh', 'mesh1', 1);
    


    %ComsolMinWavelength
    %ComsolMaxWavelength
    %ComsolWavelengthStep
    
    
    
    
    
    
    d.Value=1.0;
    d.Message = 'Saving project...';

    model.save([pathname filename]);    %finally save the model
    %mphsave(model,[pathname filename]);
    
    close(d);
    ModelUtil.remove('Model');
    
end