function BuildGDSFile(uifig, pathname, filename, BraggReflectors, PortTapers, PortTrenches, PortTrenchesCon, SiBoundary, SiBoundaryT, SiO2Boundary, AttenuatorLeft, AttenuatorRight, ShiftRotate, EchelleID, WaveguideWidth, TrenchWidth, Trench2Width, WaveguideLength, WaveguidePitch, PortEndPositions, PortDirections, Discretization, GDS2DBUnit, GDS2UserUnit, BevelAngle, BevelDistance, CladdingBlockMode, BlockModeTrenchwidth)
% GDSII Layout generation for Echelle layout program 2
% uses GDSII Toolbox v1.41 from Ulf Griesmann
% https://sites.google.com/site/ulfgri/numerical/gdsii-toolbox

% SiBoundary not used anymore

    % Layer number definitions:
    Layer_BG = 37;     % Bragg gratings (layer)
    Layer_BG_dt = 6;    % (data type)
    
    Layer_WG = 37;       % Waveguides
    Layer_WG_dt = 4;

    Layer_Taper_Obsolete = 200; %Layer_WG;    % Waveguide tapers (generally the same as Layer_WG)
    Layer_Taper_Obsolete_dt = 4; %Layer_WG_dt;
    
    Layer_Trench = 37; % Trenches near waveguides
    Layer_Trench_dt = 5;

    Layer_Trench_Obsolete = 200; % Trenches near waveguides
    Layer_Trench_Obsolete_dt = 5;
    
    Layer_Slab = Layer_WG;   % Free space region
    Layer_Slab_dt = Layer_WG_dt;
    
    Layer_SlabClad = 37;   % Cladding of free space region
    Layer_SlabClad_dt = 5;
    
    Layer_NoFill = 1158;    % Free space region and Ports: indicate that no dummy structures of any sort should be placed
    Layer_NoFill_dt = 0;
    
    
    Layer_ID = 100;      % ID text for identification
    Layer_ID_dt = 0;
    
    Layer_Doping = 25;  % Doping for attenuation
    Layer_Doping_dt = 0;
    
    Layer_Marker = 100;  % Markers (if there are any)
    Layer_Marker_dt = 0;

    % start to write a small spt-script-file for Synopsys OptoDesigner to load and show the grating
    SPTstring='';
    SPTstring=[SPTstring '// Clear Info Window every script execution.\n'];
    SPTstring=[SPTstring 'dsp::clearInfoWin();\n'];
    SPTstring=[SPTstring 'mask::clearFiles();\n'];
    SPTstring=[SPTstring '#include @layout;\n'];
    SPTstring=[SPTstring '#include @mask/tech;\n'];
    SPTstring=[SPTstring '\n'];
    SPTstring=[SPTstring 'int gdsLayer_BG = mask::AddLayer("WGTre", LD(' num2str(Layer_BG) ',' num2str(Layer_BG_dt) '), RGB(0,   128, 255), true, false);\n'];
    SPTstring=[SPTstring 'int gdsLayer_WG = mask::AddLayer("WGCor", LD(' num2str(Layer_WG) ',' num2str(Layer_WG_dt) '), RGB(128, 128, 255), true, false);\n'];
%    SPTstring=[SPTstring 'int gdsLayer_WG = mask::AddLayer("WG", LD(' num2str(Layer_Taper) ',' num2str(Layer_Taper_dt) '), RGB(128, 128, 255), true, true);\n'];
%    SPTstring=[SPTstring 'int gdsLayer_Trench = mask::AddLayer("Trench", LD(' num2str(Layer_Trench) ',' num2str(Layer_Trench_dt) '), RGB(0 ,  200, 255), true, false);\n'];
    SPTstring=[SPTstring 'int gdsLayer_Trench = mask::AddLayer("WGClad", LD(' num2str(Layer_Trench) ',' num2str(Layer_Trench_dt) '), RGB(0 ,  200, 255), true, false);\n'];
%    SPTstring=[SPTstring 'int gdsLayer_Slab = mask::AddLayer("Slab", LD(' num2str(Layer_Slab) ',' num2str(Layer_Slab_dt) '), RGB(192, 192, 192), true, false);\n'];
%    SPTstring=[SPTstring 'int gdsLayer_SlabClad = mask::AddLayer("SlabClad", LD(' num2str(Layer_SlabClad) ',' num2str(Layer_SlabClad_dt) '), RGB(128, 128, 128), true, false);\n'];
    SPTstring=[SPTstring 'int gdsLayer_ID  = mask::AddLayer("LogoTxt", LD(' num2str(Layer_ID) ',' num2str(Layer_ID_dt) '), RGB(128, 128, 128  ), true, false);\n'];
    SPTstring=[SPTstring 'int gdsLayer_Doping = mask::AddLayer("NBody", LD(' num2str(Layer_Doping) ',' num2str(Layer_Doping_dt) '), RGB(255, 128, 0  ), true, false);\n'];
%    SPTstring=[SPTstring 'int gdsLayer_Marker = mask::AddLayer("Marker", LD(' num2str(Layer_Marker) ',' num2str(Layer_Marker_dt) '), RGB(255, 128, 0  ), true, false);\n'];
    SPTstring=[SPTstring '\n\n'];
    

    % Unit definitions
%    UserUnits = 1e-6;        % use Micrometers as standard unit
%    DatabaseUnits = 1e-10;    % GDS calculations are performed using 1/10 of a Nanometer accuracy
    UserUnits = GDS2UserUnit;      % unit for coordinates, e.g. 1e-6 are micrometers
    DatabaseUnits = GDS2DBUnit;    % GDS calculations are performed on this grid, e.g. 1e-9 uses 1 nm accuracy

    % Variable preconditioning
    AbsoluteMinimumRadius=1e100;


%{    
    function WGCoords=MakeStraightWaveguide(Start, End, Width)
    % Generates the coordinates [x1 y1; x2 y2; x3 y3; x4 y4] of a straight waveguide between Start [xs ys]
    % and End [xe ye] with Width (µm)
        dir1=(End-Start)/norm(End-Start);
        perpend=[-dir1(2) dir1(1)];
        WGCoords= [ Start+perpend*Width/2;
                    Start-perpend*Width/2;
                    End-perpend*Width/2;
                    End+perpend*Width/2];
        
    end
%}    

    function BorderCoordinates = BezierBend(StartPos, StartDir, EndPos, EndDir, ControlPointDistance, Discretization, Width1, Width2)
    % Calculates cubic bezier curve between start and end point with directions (start point: direction to curve, end point
    % direction out of curve), Control point construction through distance parameter.

%        hold off;

%        StartPos=[40, 10]
%        StartDir=-1.3
%        EndPos=[2.5, -50]
%        EndDir=-pi/2

    %    Width=0.4;
%        Width=2;

%        ControlPointDistance=0.3;

%        Discretization=100;


        %dir1=(End-Start)/norm(End-Start);

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
        for BC_i=0:Discretization
            x=BC_i*1/Discretization;
            %calculate first level lines for De-Castteljau algorithm...
            L1(1,:)=StartPos+x*(BSPos-StartPos);
            L1(2,:)=BSPos+x*(BEPos-BSPos);
            L1(3,:)=BEPos+x*(EndPos-BEPos);
            %calculate second level lines for De-Castteljau algorithm...
            L2(1,:)=L1(1,:)+x*(L1(2,:)-L1(1,:));
            L2(2,:)=L1(2,:)+x*(L1(3,:)-L1(2,:));
            %calculate the point of the Bezier spline
            BPos=L2(1,:)+x*(L2(2,:)-L2(1,:));
            BDir=L2(2,:)-L2(1,:);   %Direction vector...
            BDirP=[-BDir(2), BDir(1)]/norm(BDir);  %...and direction vector perpendicular to that, normalized

            %Width=Width1+(Width2-Width1)*x;
            Width=Width1+(Width2-Width1)*(x * x * (3 - 2 * x));    % https://en.wikipedia.org/wiki/Smoothstep
            %Width=Width1+(Width2-Width1)*(x * x * x * (x * (x * 6 - 15) + 10)); %https://en.wikipedia.org/wiki/Smoothstep
            Path1(BC_i+1,:)=BPos+Width/2*BDirP;    %Path on one side of Bezier curve
            Path2(BC_i+1,:)=BPos-Width/2*BDirP;    % Path on other side of Bezier curve
        end
        
        minRad=1e100;
        for BC_i=1:Discretization
            iPos=LineLineIntersectionB2( Path1(BC_i,:) , Path2(BC_i,:) , Path1(BC_i+1,:) , Path2(BC_i+1,:) );
            pPos=(Path1(BC_i,:)+Path2(BC_i,:)+Path1(BC_i+1,:)+Path2(BC_i+1,:))/4;
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




    
    d = uiprogressdlg(uifig,'Title','GDSII layout generation','Message','Rotating and shifting geometry...');
    
    
    gdsii_units(UserUnits, DatabaseUnits);   % set the units for the GDSII library

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

    SiBoundaryT(1,:)=SiBoundaryT(1,:)-ShiftX;
    SiBoundaryT(2,:)=SiBoundaryT(2,:)-ShiftY;

    SiO2Boundary(1,:)=SiO2Boundary(1,:)-ShiftX;
    SiO2Boundary(2,:)=SiO2Boundary(2,:)-ShiftY;

    AttenuatorLeft(1,:)=AttenuatorLeft(1,:)-ShiftX;
    AttenuatorLeft(2,:)=AttenuatorLeft(2,:)-ShiftY;
    
    AttenuatorRight(1,:)=AttenuatorRight(1,:)-ShiftX;
    AttenuatorRight(2,:)=AttenuatorRight(2,:)-ShiftY;
    
    PortEndPositions(:,1)=PortEndPositions(:,1)-ShiftX;
    PortEndPositions(:,2)=PortEndPositions(:,2)-ShiftY;
    
%    [ShiftX ShiftY]
%    PortEndPositions
    
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
    
    for k=1:size(SiBoundaryT,2)
        help=[SiBoundaryT(1,k) SiBoundaryT(2,k)]*RotMat;
        SiBoundaryT(1,k)=help(1);
        SiBoundaryT(2,k)=help(2);
    end
    
    for k=1:size(SiO2Boundary,2)
        help=[SiO2Boundary(1,k) SiO2Boundary(2,k)]*RotMat;
        SiO2Boundary(1,k)=help(1);
        SiO2Boundary(2,k)=help(2);
    end
    
    for k=1:size(AttenuatorLeft,2)
        help=[AttenuatorLeft(1,k) AttenuatorLeft(2,k)]*RotMat;
        AttenuatorLeft(1,k)=help(1);
        AttenuatorLeft(2,k)=help(2);
    end
    
    for k=1:size(AttenuatorRight,2)
        help=[AttenuatorRight(1,k) AttenuatorRight(2,k)]*RotMat;
        AttenuatorRight(1,k)=help(1);
        AttenuatorRight(2,k)=help(2);
    end
    
    for k=1:size(PortEndPositions,1)
        help=[PortEndPositions(k,1) PortEndPositions(k,2)]*RotMat;
        PortEndPositions(k,1)=help(1);
        PortEndPositions(k,2)=help(2);
        
    end


    PortDirections=PortDirections-RotationAngle;
    
    
    
    %find leftmost and topmost coordinate of silicon boundary for positioning the ID text
    SiB_topmost=-1e100;
    SiB_leftmost=1e100;
    for k=1:size(SiBoundary,2)
        if SiBoundary(1,k)<SiB_leftmost
            SiB_leftmost=SiBoundary(1,k);
        end;
        if SiBoundary(2,k)>SiB_topmost
            SiB_topmost=SiBoundary(2,k);
        end;
    end


    
    Structure_NoFill = gds_structure('Echelle_NoFill');
    
    %Build reflector geometries
    d.Value=0.1;
    d.Message = 'Adding reflectors...';
    fprintf('-- GDSII reflectors...\n');

    Structure_Grating = gds_structure('Echelle_Grating');
    
    for i=1:size(BraggReflectors,1)
        for j=1:size(BraggReflectors,2)
            Structure_Grating(end+1) = gds_element('boundary', 'xy', BraggReflectors{i,j}' , 'layer',Layer_BG, 'dtype',Layer_BG_dt); 
        end
    end


    %Build port geometry
    d.Value=0.2;
    d.Message = 'Adding ports...';
    fprintf('-- GDSII port geometry...\n');
    
    Structure_TaperWGObsolete = gds_structure('Echelle_Tapers_Obsolete');
    Structure_TaperTrench = gds_structure('Echelle_Tapers_Trenches');
    Structure_TaperTrenchObsolete = gds_structure('Echelle_Tapers_Trenches_Obsolete');

    for i=1:size(PortTapers,1)
        Structure_TaperWGObsolete(end+1) = gds_element('boundary', 'xy', PortTapers{i}' , 'layer',Layer_Taper_Obsolete, 'dtype',Layer_Taper_Obsolete_dt);
        Structure_TaperTrenchObsolete(end+1) = gds_element('boundary', 'xy', PortTrenches{i}' , 'layer',Layer_Trench_Obsolete, 'dtype',Layer_Trench_Obsolete_dt);
    end
    Structure_TaperTrench(end+1) = gds_element('boundary', 'xy', PortTrenchesCon' , 'layer',Layer_Trench, 'dtype',Layer_Trench_dt);

    Structure_NoFill(end+1) = gds_element('boundary', 'xy', PortTrenchesCon' , 'layer',Layer_NoFill, 'dtype',Layer_NoFill_dt);

    % Build large free space area, the SiBoundary
    d.Value=0.3;
    d.Message = 'Adding free space area...';
    fprintf('-- GDSII slab...\n');
    
    Structure_Slab = gds_structure('Echelle_Slab');
    
%    Structure_Slab(end+1) = gds_element('boundary', 'xy', SiBoundary' , 'layer',Layer_Slab, 'dtype',Layer_Slab_dt);
    Structure_Slab(end+1) = gds_element('boundary', 'xy', SiBoundaryT' , 'layer',Layer_Slab, 'dtype',Layer_Slab_dt);


    Structure_SlabClad = gds_structure('Echelle_SlabClad');
    Structure_SlabClad(end+1) = gds_element('boundary', 'xy', SiO2Boundary' , 'layer',Layer_SlabClad, 'dtype',Layer_SlabClad_dt);

    Structure_NoFill(end+1) = gds_element('boundary', 'xy', SiO2Boundary' , 'layer',Layer_NoFill, 'dtype',Layer_NoFill_dt);



    % Build attenuation regions
    d.Value=0.4;
    d.Message = 'Adding attenuators...';
    fprintf('-- GDSII attenuators...\n');
    
    Structure_AttenuatorLeft = gds_structure('Echelle_Attenuator_Left');
    Structure_AttenuatorRight = gds_structure('Echelle_Attenuator_Right');
    
    Structure_AttenuatorLeft(end+1) = gds_element('boundary', 'xy', AttenuatorLeft' , 'layer', Layer_Doping, 'dtype',Layer_Doping_dt);
    Structure_AttenuatorRight(end+1) = gds_element('boundary', 'xy', AttenuatorRight' , 'layer', Layer_Doping, 'dtype',Layer_Doping_dt);



    % Build ID text
    d.Value=0.5;
    d.Message = 'Adding ID text...';
    fprintf('-- GDSII ID text...\n');

    Structure_ID = gds_structure('Echelle_ID');
%    Structure_ID(end+1) = gdsii_ptext(EchelleID, [SiB_leftmost, SiB_topmost+1], 50, Layer_ID);
%    Structure_ID(end+1) = gdsii_boundarytext(EchelleID, [SiB_leftmost, SiB_topmost+1], 50, 0, Layer_ID, Layer_ID_dt);
    Structure_ID(end+1) = gdsii_boundarytext_Bevel2(EchelleID, [SiB_leftmost, SiB_topmost+1], 50, 0, Layer_ID, Layer_ID_dt, BevelAngle, BevelDistance/1000);




    % Build waveguide bundle
    d.Value=0.6;
    d.Message = 'Adding waveguide bundle...';
    fprintf('-- GDSII waveguide bundle...\n');
    
    Structure_WGBundle = gds_structure('Waveguide_Bundle');
    Structure_TrenchBundle = gds_structure('Waveguide_Bundle_Trenches');
%{
    % first make some 1nm long 'markers', migth be deleted later...
    WGEndMarker = [ -WaveguideWidth/2/1e3 0;
                    -WaveguideWidth/2/1e3 1e-3;
                     WaveguideWidth/2/1e3 1e-3;
                     WaveguideWidth/2/1e3 0];
    for k=1:size(PortEndPositions,1)
        if (k==size(PortEndPositions,1))
            i=0;
        else
            i=k;
        end
        Marker=[WGEndMarker(:,1)+PortEndPositions(end,1)-0+i*WaveguidePitch/1e3 , WGEndMarker(:,2)+PortEndPositions(end,2)-WaveguideLength/1e3+0];
        Structure_WGBundle(end+1) = gds_element('boundary', 'xy', Marker , 'layer', Layer_Marker);
    end
%}
    % first (the common) waveguide is straight, so a long straight waveguide is suitable
%    Structure_WGBundle(end+1) = gds_element('boundary', 'xy', MakeStraightWaveguide( PortEndPositions(end,:), PortEndPositions(end,:)+[0 -WaveguideLength/1e3],WaveguideWidth/1e3) , 'layer', Layer_WG);
%    Structure_TrenchBundle(end+1) =  gds_element('boundary', 'xy', MakeStraightWaveguide( PortEndPositions(end,:), PortEndPositions(end,:)+[0 -WaveguideLength/1e3],WaveguideWidth/1e3+2*TrenchWidth/1e3) , 'layer', Layer_Trench);

    % first (the common) waveguide is built, because it's the first waveguide, but comes last in the position arrays
    Structure_WGBundle(end+1) =  gds_element('boundary', 'xy',...
        BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3 , WaveguideWidth/1e3),...
        'layer', Layer_WG, 'dtype',Layer_WG_dt);
    if (CladdingBlockMode~=true)
        Structure_TrenchBundle(end+1) =  gds_element('boundary', 'xy',...
            BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3),...
            'layer', Layer_Trench, 'dtype',Layer_Trench_dt);
        Structure_NoFill(end+1) =  gds_element('boundary', 'xy',...
            BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3),...
            'layer', Layer_NoFill, 'dtype',Layer_NoFill_dt);
    end

    % then also all other waveguides are built
    for k=1:size(PortEndPositions,1)-1
        Structure_WGBundle(end+1) =  gds_element('boundary', 'xy',...
            BezierBend(PortEndPositions(k,:), PortDirections(k), PortEndPositions(end,:)+[k*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3 , WaveguideWidth/1e3),...
            'layer', Layer_WG, 'dtype',Layer_WG_dt);
        if (CladdingBlockMode~=true)
            Structure_TrenchBundle(end+1) =  gds_element('boundary', 'xy',...
                BezierBend(PortEndPositions(k,:), PortDirections(k), PortEndPositions(end,:)+[k*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3),...
                'layer', Layer_Trench, 'dtype',Layer_Trench_dt);
            Structure_NoFill(end+1) =  gds_element('boundary', 'xy',...
                BezierBend(PortEndPositions(k,:), PortDirections(k), PortEndPositions(end,:)+[k*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3),...
                'layer', Layer_NoFill, 'dtype',Layer_NoFill_dt);
        end
    end

    if (CladdingBlockMode==true)
        if (TrenchWidth<BlockModeTrenchwidth) TW1=BlockModeTrenchwidth;
        else                                  TW1=TrenchWidth;
        end
        
        if (Trench2Width<BlockModeTrenchwidth) TW2=BlockModeTrenchwidth;
        else                                   TW2=Trench2Width;
        end
        
        DummyCoordinates = BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3);
        DummyCoordinatesB = BezierBend(PortEndPositions(end,:), PortDirections(end), PortEndPositions(end,:)+[0, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TW1/1e3 , WaveguideWidth/1e3+2*TW2/1e3);
        
        
        BuildBorder = [DummyCoordinates((Discretization+1),:)];
        BuildBorder = [BuildBorder ; DummyCoordinates((Discretization+2),:)];
        BuildBorder = [BuildBorder ; DummyCoordinatesB((Discretization+2):end,:)];
        BuildBorder = [BuildBorder ; DummyCoordinates(1,:)];
        
        BuildBorderBottom = [];
        
        for k=1:size(PortEndPositions,1)-2
            DummyCoordinates = BezierBend(PortEndPositions(k,:), PortDirections(k), PortEndPositions(end,:)+[k*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3);
            BuildBorder = [ DummyCoordinates((Discretization+1),:) ; DummyCoordinates((Discretization+2),:) ; BuildBorder ; DummyCoordinates(end,:) ; DummyCoordinates(1,:)];  % first add bottom border points before rest, then already existing points and then border points on top to keep clockwise border definition
        end
        
        DummyCoordinates = BezierBend(PortEndPositions(size(PortEndPositions,1)-1,:), PortDirections(size(PortEndPositions,1)-1), PortEndPositions(end,:)+[(size(PortEndPositions,1)-1)*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TrenchWidth/1e3 , WaveguideWidth/1e3+2*Trench2Width/1e3);        
        DummyCoordinatesB = BezierBend(PortEndPositions(size(PortEndPositions,1)-1,:), PortDirections(size(PortEndPositions,1)-1), PortEndPositions(end,:)+[(size(PortEndPositions,1)-1)*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , 0.4, Discretization, WaveguideWidth/1e3+2*TW1/1e3 , WaveguideWidth/1e3+2*TW2/1e3);        
        
        BuildBorder = [BuildBorder ; DummyCoordinates(end,:)];
        BuildBorder = [BuildBorder ; DummyCoordinatesB(1:(Discretization+1),:)];
        BuildBorder = [BuildBorder ; DummyCoordinates((Discretization+1),:)];
        BuildBorder = [BuildBorder ; DummyCoordinates((Discretization+2),:)];

        Structure_TrenchBundle(end+1) =  gds_element('boundary', 'xy', BuildBorder, 'layer', Layer_Trench, 'dtype',Layer_Trench_dt);
        Structure_NoFill(end+1)       =  gds_element('boundary', 'xy', BuildBorder, 'layer', Layer_NoFill, 'dtype',Layer_NoFill_dt);

    end


    fprintf('     Minimum waveguide bending radius: %g µm\n', AbsoluteMinimumRadius);
%    Structure_WGBundle=MakeCurve(Structure_WGBundle, PortEndPositions(1,:), PortDirections(1), PortEndPositions(end,:)+[1*WaveguidePitch/1e3, -WaveguideLength/1e3+0], -pi/2 , WaveguideWidth/1e3, CurveRadius/1e3, Layer_WG);
    
    
    
    fprintf('-- GDSII finishing...\n');

%{
    marker_size = [-500 -500; -500 500; 500 500; 500 -500];
    marker_pos = [2000 2000];
    
    marker1 = [marker_size(:,1)+marker_pos(1) marker_size(:,2)+marker_pos(2)];
    marker2 = [marker_size(:,1)+marker_pos(1) marker_size(:,2)-marker_pos(2)];
    marker3 = [marker_size(:,1)-marker_pos(1) marker_size(:,2)-marker_pos(2)];
    marker4 = [marker_size(:,1)-marker_pos(1) marker_size(:,2)+marker_pos(2)];

    marker = gds_element('boundary', 'xy', marker1 , 'layer',Layer_Marker, 'dtype',Layer_Marker_dt)+...
    gds_element('boundary', 'xy', marker2 , 'layer',Layer_Marker, 'dtype',Layer_Marker_dt)+...
    gds_element('boundary', 'xy', marker3 , 'layer',Layer_Marker, 'dtype',Layer_Marker_dt)+...
    gds_element('boundary', 'xy', marker4 , 'layer',Layer_Marker, 'dtype',Layer_Marker_dt);
%}
    
    Structure_Marker = gds_structure('Marker');
%    Structure_Marker(end+1) = marker;
    
    Structure_Full = gds_structure('Echelle');
    Structure_Full = add_ref(Structure_Full, [Structure_ID Structure_Marker Structure_Grating Structure_TaperWGObsolete Structure_TaperTrench Structure_TaperTrenchObsolete Structure_Slab Structure_SlabClad Structure_AttenuatorLeft Structure_AttenuatorRight Structure_WGBundle Structure_TrenchBundle Structure_NoFill]);
    
%    GDSLib = gds_library('Echelle', 'uunit', 1e-6, 'dbunit', 1e-10, Structure_ID, Structure_Marker);
    GDSLib = gds_library('Echelle', 'uunit', UserUnits, 'dbunit', DatabaseUnits, Structure_Full, Structure_ID, Structure_Marker, Structure_Grating, Structure_TaperWGObsolete, Structure_TaperTrench, Structure_TaperTrenchObsolete, Structure_Slab, Structure_SlabClad, Structure_AttenuatorLeft, Structure_AttenuatorRight, Structure_WGBundle, Structure_TrenchBundle, Structure_NoFill);

    fprintf('-- GDSII saving...\n');

   
    d.Value=1.0;
    d.Message = 'Saving layout...';
    
    write_gds_library(GDSLib, ['!', pathname, filename]);


    % continue to write the small spt-script-file for Synopsys OptoDesigner to load and show the grating
    SPTstring=[SPTstring 'layout myEchelle(devMode="DEMUX" KamusDropDown "MUX\tDEMUX")\n'];
    SPTstring=[SPTstring '    dlgname "Echelle grating"\n'];
    SPTstring=[SPTstring '    Domain_Optics\n'];
    SPTstring=[SPTstring '    AuthorInfo     "Dr.-Ing. Marc Schneider, KIT, Germany"\n'];
    SPTstring=[SPTstring '    VersionInfo    "Generated by Echelle Layout Program 2 v1.3 from 2022-09-05 by Dr.-Ing. Marc Schneider (KIT) "\n'];
    SPTstring=[SPTstring '    LicenseInfo    "TBD"\n'];
    SPTstring=[SPTstring '    Maturity       Research\n'];
    SPTstring=[SPTstring '    MaskLayers     "based on IMEC ISIPP50G"\n'];
    SPTstring=[SPTstring '    Disclosure     "For Internal Use Only"\n'];
    SPTstring=[SPTstring '    TexDoc         "Just the title"\n'];
    SPTstring=[SPTstring '{\n'];
    SPTstring=[SPTstring '    string mcs=mask::CSget();\n'];
    SPTstring=[SPTstring '    // ----------------------------------------------------\n'];
    SPTstring=[SPTstring '    // Change this when moving or/and renaming the GDS file\n'];
    SPTstring=[SPTstring '    gdsfile("' strrep(pathname,'\','\\') strrep(filename,'\','\\') '", 0) myGDS;\n'];
    SPTstring=[SPTstring '    // ----------------------------------------------------\n'];
    SPTstring=[SPTstring '    \n'];
    SPTstring=[SPTstring '    myGDS.place( wher->this@origin : "Echelle", 0);\n'];

    PE=PortEndPositions(end,:)+[0*WaveguidePitch/1e3, -WaveguideLength/1e3+0];

    SPTstring=[SPTstring '    if (devMode=="MUX") { ml::setPort(this:out0->this@origin+[' num2str(PE(1)) ', ' num2str(PE(2)) ', 270]); }\n'];
    SPTstring=[SPTstring '    else                { ml::setPort(this:in0 ->this@origin+[' num2str(PE(1)) ', ' num2str(PE(2)) ',  90]); }\n'];
    SPTstring=[SPTstring '    \n'];
    SPTstring=[SPTstring '    for (int ii=0; ii<=' num2str(size(PortEndPositions,1)-1-1) '; ii++) {\n'];
    SPTstring=[SPTstring '        if (devMode=="MUX") {ml::setPort(this:"in"+ii->this@origin+[' num2str(PE(1)) '+(ii+1)*' num2str(WaveguidePitch/1e3) ', ' num2str(PE(2)) ',  90]);}\n'];
    SPTstring=[SPTstring '        else                {ml::setPort(this:"out"+ii->this@origin+[' num2str(PE(1)) '+(ii+1)*' num2str(WaveguidePitch/1e3) ',' num2str(PE(2)) ', 270]);}\n'];
    SPTstring=[SPTstring '    }\n'];
    SPTstring=[SPTstring '    \n'];
    SPTstring=[SPTstring '    if (devMode=="MUX") { this{"nin"}=' num2str(size(PortEndPositions,1))-1 '; this{"nout"}=1; } // number of input and output ports\n'];
    SPTstring=[SPTstring '    else                { this{"nout"}=' num2str(size(PortEndPositions,1))-1 '; this{"nin"}=1; } // number of input and output ports\n'];
    SPTstring=[SPTstring '    \n'];
    SPTstring=[SPTstring '    mask::CSselect(mcs);\n'];
    SPTstring=[SPTstring '}\n'];
    SPTstring=[SPTstring '   \n'];
    SPTstring=[SPTstring 'var BB1 = ml::myEchelle( in0 -> [0,0,90] : "DEMUX");\n'];
    SPTstring=[SPTstring '//var BB2 = ml::myEchelle( out0 -> [100,0,-90] : "MUX");\n'];
    SPTstring=[SPTstring '\n'];


    % save the small spt-script-file for Synopsys OptoDesigner to load and show the grating
    filename4OD=[filename,'_ports.spt'];    
    ODfileID=fopen([pathname, filename4OD],'w');
    fprintf(ODfileID,SPTstring);
    fclose(ODfileID);

    
    
    
    close(d);


end
