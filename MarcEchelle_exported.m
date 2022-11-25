classdef MarcEchelle_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        NumberofChannelsLabel           matlab.ui.control.Label
        NumberofChannelsSpinner         matlab.ui.control.Spinner
        ChannelPropsUITable             matlab.ui.control.Table
        CommonPropsUITable              matlab.ui.control.Table
        CommonPortPositionLabel         matlab.ui.control.Label
        DistinctPortsLabel              matlab.ui.control.Label
        ChannelNumberSpinnerLabel       matlab.ui.control.Label
        TSPChannel1                     matlab.ui.control.Spinner
        ChannelNumberSpinner_2Label     matlab.ui.control.Label
        TSPChannel2                     matlab.ui.control.Spinner
        TSPChannel1Label                matlab.ui.control.Label
        TSPChannel2Label                matlab.ui.control.Label
        ProgramTitleLabel               matlab.ui.control.Label
        ProgramSubTitle1Label           matlab.ui.control.Label
        MaterialChooserDropDownLabel    matlab.ui.control.Label
        MaterialChooserDropDown         matlab.ui.control.DropDown
        ITUChannelChooserButton         matlab.ui.control.Button
        Footnote1                       matlab.ui.control.Label
        GratingorderEditFieldLabel      matlab.ui.control.Label
        GratingorderEditField           matlab.ui.control.NumericEditField
        MinradiusmEditFieldLabel        matlab.ui.control.Label
        MinradiusmEditField             matlab.ui.control.NumericEditField
        MaxradiusmEditFieldLabel        matlab.ui.control.Label
        MaxradiusmEditField             matlab.ui.control.NumericEditField
        TSPradiusoffsetmEditFieldLabel  matlab.ui.control.Label
        TSPradiusoffsetmEditField       matlab.ui.control.NumericEditField
        GratingopeningangleLabel        matlab.ui.control.Label
        GratingopeningangleEditField    matlab.ui.control.NumericEditField
        CalculateButton                 matlab.ui.control.Button
        ImageKIT                        matlab.ui.control.Image
        Image                           matlab.ui.control.Image
        ResultsPanel                    matlab.ui.container.Panel
        toeasecalculationsasseenfromcommonportpositionLabel  matlab.ui.control.Label
        RealminradiusmEditFieldLabel    matlab.ui.control.Label
        RealminradiusmField             matlab.ui.control.NumericEditField
        RealmaxradiusmEditFieldLabel    matlab.ui.control.Label
        RealmaxradiusmField             matlab.ui.control.NumericEditField
        RealopeningangleEditFieldLabel  matlab.ui.control.Label
        RealopeningangleField           matlab.ui.control.NumericEditField
        GratingpointsfoundEditFieldLabel  matlab.ui.control.Label
        GratingpointsfoundField         matlab.ui.control.NumericEditField
        reducedLabel                    matlab.ui.control.Label
        AreaReducedField                matlab.ui.control.NumericEditField
        fullLabel                       matlab.ui.control.Label
        AreaFullField                   matlab.ui.control.NumericEditField
        AreaoffreespaceareamLabel       matlab.ui.control.Label
        SiO2Label                       matlab.ui.control.Label
        AreaSiO2Field                   matlab.ui.control.NumericEditField
        PortSearchBoundingBoxUITable    matlab.ui.control.Table
        BoundingBoxforPortSearchLabel   matlab.ui.control.Label
        ShowEllipsesCheckBox            matlab.ui.control.CheckBox
        ShowCirclesCheckBox             matlab.ui.control.CheckBox
        ShowPortSearchIntersectionsCheckBox  matlab.ui.control.CheckBox
        ShowGratingAdditionsCheckBox    matlab.ui.control.CheckBox
        BragggratingdefinitionPanel     matlab.ui.container.Panel
        PeriodsEditFieldLabel           matlab.ui.control.Label
        PeriodsEditField                matlab.ui.control.NumericEditField
        PeriodlengthnmEditFieldLabel    matlab.ui.control.Label
        PeriodlengthEditField           matlab.ui.control.NumericEditField
        SiO2lengthnmEditFieldLabel      matlab.ui.control.Label
        SiO2lengthEditField             matlab.ui.control.NumericEditField
        BoundarydistancenmLabel         matlab.ui.control.Label
        BoundarydistanceEditField       matlab.ui.control.NumericEditField
        ReflectorspacingnmLabel         matlab.ui.control.Label
        ReflectorSpacingEditField       matlab.ui.control.NumericEditField
        PortlayoutdefinitionPanel       matlab.ui.container.Panel
        WaveguidewidthnmEditFieldLabel  matlab.ui.control.Label
        WaveguidewidthEditField         matlab.ui.control.NumericEditField
        TaperendwidthnmLabel            matlab.ui.control.Label
        TaperendwidthEditField          matlab.ui.control.NumericEditField
        TaperlengthnmLabel              matlab.ui.control.Label
        TaperlengthEditField            matlab.ui.control.NumericEditField
        WGstublengthnmLabel             matlab.ui.control.Label
        WaveguidestublengthEditField    matlab.ui.control.NumericEditField
        TapershiftnmLabel               matlab.ui.control.Label
        TapershiftEditField             matlab.ui.control.NumericEditField
        TrenchwidthnmLabel              matlab.ui.control.Label
        TrenchwidthEditField            matlab.ui.control.NumericEditField
        Trench2widthnmLabel             matlab.ui.control.Label
        Trench2widthEditField           matlab.ui.control.NumericEditField
        COMSOLModelButton               matlab.ui.control.Button
        GDS2LayoutButton                matlab.ui.control.Button
        COMSOLSimulationPanel           matlab.ui.container.Panel
        MinWavelengthnmLabel            matlab.ui.control.Label
        ComsolMinWavelengthEditField    matlab.ui.control.NumericEditField
        MaxWavelengthnmLabel            matlab.ui.control.Label
        ComsolMaxWavelengthEditField    matlab.ui.control.NumericEditField
        WavelengthstepnmLabel           matlab.ui.control.Label
        ComsolWavelengthStepEditField   matlab.ui.control.NumericEditField
        ComsolBorderButton              matlab.ui.control.StateButton
        ComsolSBendsButton              matlab.ui.control.StateButton
        FreespaceareaLabel              matlab.ui.control.Label
        NoofPartitionsEditFieldLabel    matlab.ui.control.Label
        NoOfPartitionsEditField         matlab.ui.control.NumericEditField
        GDSIIDefinitionsPanel           matlab.ui.container.Panel
        EchelleNameIDEditFieldLabel     matlab.ui.control.Label
        EchelleNameIDEditField          matlab.ui.control.EditField
        GDSdatabaseunitLabel            matlab.ui.control.Label
        GDS2DBUnitEditField             matlab.ui.control.NumericEditField
        GDSuserunitLabel                matlab.ui.control.Label
        GDS2UserUnitEditField           matlab.ui.control.NumericEditField
        SbendWaveguidesPanel            matlab.ui.container.Panel
        WaveguidelengthnmLabel          matlab.ui.control.Label
        WaveguideLengthEditField        matlab.ui.control.NumericEditField
        EndpointpitchnmLabel            matlab.ui.control.Label
        WaveguidePitchEditField         matlab.ui.control.NumericEditField
        CurvediscretizationLabel        matlab.ui.control.Label
        WGCurveDiscretizationEditField  matlab.ui.control.NumericEditField
        ProgramSubTitle2Label           matlab.ui.control.Label
        BevelingPanel                   matlab.ui.container.Panel
        AttenuatorBevelnmLabel          matlab.ui.control.Label
        GDS2AttBevelEditField           matlab.ui.control.NumericEditField
        minAngleLabel                   matlab.ui.control.Label
        GDS2AttBevelAngleEditField      matlab.ui.control.NumericEditField
        BevelnmLabel                    matlab.ui.control.Label
        GDS2BevelEditField              matlab.ui.control.NumericEditField
        BlazeshiftmEditFieldLabel       matlab.ui.control.Label
        BlazeshiftmEditField            matlab.ui.control.NumericEditField
        ParallelComputingCheckBox       matlab.ui.control.CheckBox
        CladdingBlockModeCheckBox       matlab.ui.control.CheckBox
        BlockTrenchwidthnmLabel         matlab.ui.control.Label
        BlockModeTrenchwidthEditField   matlab.ui.control.NumericEditField
    end

    
    properties (Access = private)
        % Echelle Layout Program v2
        % by Dr.-Ing. Marc Schneider (KIT)
        % programming started at 2021-03-30
        % parallelized channel search added 22.02.2022
        % Block mode for increased cladding width added for GDS output 23.08.2022
        % ReflectorSpacing added to avoid DRC errors 05.09.2022
        %
        % required toolboxes/add-ons/programs:
        % - GDSII Toolbox v1.41 from Ulf Griesmann
        % - Comsol
        %
        % prefered additional toolboxes:
        % - Parallel Computing Toolbox (Matlab)
        %
        
        ProgramVersionString1='Echelle Layout Program 2'
        ProgramVersionString2='by Dr.-Ing. Marc Schneider (KIT)'
        ProgramVersionString3='v1.3 from 2022-09-05'
        
        %convenience-'#defines'
        PushFigureWindowToSpecificPosition=true;
        
        % Table of effective refractive indices for different Si-thicknesses calculated by COMSOL with "E:\Dokumente\Simulationen\COMSOL\EffectiveRefractiveIndex\EffectiveRefractiveIndex_ModalAnalysis_li293.mph"
        % Si thicknesses: 215-225nm + 245-255nm, SiO2 thickness bottom 2µm, top 1µm, data from COMSOL built in: Si: Li-293K, SiO2: Malitson
%        EffectiveRefractiveIndicesFileName='E:\Dokumente\Simulationen\COMSOL\EffectiveRefractiveIndex\EffectiveRefractiveIndices_Wavelength1480-1630nm_SiliconThickness215-225nm+245-255nm.xlsx';
%        EffectiveRefractiveIndicesFileName='D:\EigeneProgramme\Matlab\Test2\EffectiveRefractiveIndices_Wavelength1480-1630nm_SiliconThickness215-225nm+245-255nm.xlsx';
        EffectiveRefractiveIndicesFileName='.\EffectiveRefractiveIndices_Wavelength1480-1630nm_SiliconThickness205-225nm+245-255nm.xlsx';
        EffectiveRefractiveIndices  % Table with the effective refractive indices and the wavelength in the first column
        EffectiveRefractiveIndexDefault='SiThickness220nm'
        EffectiveRefractiveIndex % Current matrix of wavelength in first column and effective refractive index of choosen material in second column
        
        
        CommonProps  % Table with properties of common port (X,Y-Position in µm)
        ChannelProps % Table with properties of the channels (Channel Number, X,Y-Position in µm, Wavelength in nm, refractive index)
        NoOfChannels % Number of distinct channels, excluding the common channel
        TSPChannels  % Channel numbers of the two initial focus points of the 'Two Stigmatic Points' method.
        GratingOrder % Diffraction order of the Echelle grating
        RadiusMin    % Radius of large half axis of the smallest ellipse defining the grating points with the common port and the 1. stigmatic point as foci
        RadiusMax    % Maximum Radius of large half axis of the largest ellipse defining the grating points with the common port and the 1. stigmatic point as foci
        RadiusOffset2TSP % Offset for the radius of the large half axis of the ellipse with the common port and the 2. stigmatic point as foci to shift the grating point trajectory
        GratingOpeningAngle % Opening angle of the grating defining the maximum radius of the ellipses besids RadiusMax
        BlazeShift   % Offset to PortsCenterOfGravity to calculate blaze angle of grating in the direction of the connection between the common port and the center of gravity of the distinct ports
        BlazeShiftXY % Absolute position offset for BlazeShift
        
        PortSearchBoundingBox   %Table with coordinates of bounding box in which the circle intersections for searching the position of non-TSP-ports are averaged; 1. column x-coordinates, 2. column y-coordinates
        
        fastFigure  % handle to additional figure window
        fastFigureNumber % Number of that figure window
        fastAxes    % handle to axes in additional figure windows
        fastFigurePosition  % Position and size of figure window
        
        ITUChannelChooserApp    % Sub-App to choose and set ITU channels
        ITUChannelChooserParameters % Parameters to ITUChannelChooserApp
        
        DialogBoxApp    % Dialog box app
        DialogBoxParameters %Parameters to test dialog box window


        EllipseIntersectionsSearchWedgeOpeningAngle=20; % Opening angle to search consecutive ellipse intersections
        EllipseIntersections    %Array of intersection points of the two stigmatic point ellipses with [AbsoluteX, AbsoluteY, RelativeR, RelativePhi], AbsoluteX and-Y are the positions of the intersection points, RelativeR and -Phi are the positions in polar coordinates, relative to the common port position

        PortsCenterOfGravity    %Average position of all ports, used to calculate reflector tilt angle [x, y]
        GratingCenter   %Center of grating, used to calculate tilt angle of ports [x, y]
        GratingPointsReflectorDirections %angle in radians of each reflector on the grating, calculated with respect to port center of gravity
        GratingPointsReflectorHalfAngles    %angle in radians of line through each pair of neighboring grating points, calculated with respect to port center of gravity, defining the borders of the bragg gratings; first and last are outside the grating trajectory
        PortDirections  %angle in radians of each port, calculated with respect to grating center; common port is LAST port (app.PortDirections(app.NoOfChannels+1)) !!!
        PortWGPositions %center point of wavveguide end of each port [x, y; x, y; ...]; common port is LAST port (app.PortDirections(app.NoOfChannels+1)) !!!
        GratingTrajectoryReversed=0;    % flag which indicates, if the grating trajectory goes from left to right (=0) or from right to left (=1)
        
        BraggPeriods    % number of periods of Bragg reflectors as reflecting elements of Echelle grating
        BraggPeriodLength   % length of one Bragg reflector period consisting of SiO2-part and Si-part
        BraggSiO2Length     % length of the SiO2-part of one Bragg reflector
        BoundaryDistance    % distance between bragg reflectors and Si-boundary, used for COMSOL model and layout (in layout also distance between boundary and ports)
%        BraggReflectors % Cell array with all corner coordinates of all Bragg reflectors: rows: reflector (one for each grating point), columns: each period of the reflectors, cell: coordinates of the corners [x1 y1; x2 y2; x3 y3; x4 y4]
        BraggReflectors % Cell array with all corner coordinates of all Bragg reflectors: rows: reflector (one for each grating point), columns: each period of the reflectors, cell: coordinates of the corners [x1 x2 x3 x4; y1 y2 y3 y4]
        ReflectorSpacing    % distance between two bragg reflectors to avoid sharp angles and small feature errors from design rule checks
        
        TaperLength     % length of taper for channel port
        TaperWidth      % width of taper end at free space end
        TaperWGWidth    % width of taper end at waveguide end
        WaveguideStubLength % length of waveguide stub at the waveguide end of the taper with constant width
        TaperShift      % axial shift of tapers; negative values towards grating trajectory, positive values away from grating trajectory
        TrenchWidth     % width of SiO2 trench around the waveguide; the width is the additional width of the trench to one side, the resulting full width is 2*TrenchWidth+TaperWGWidth
        Trench2Width    % width of SiO2 trench around the waveguide at S-bend port end; the width is the additional width of the trench to one side, the resulting full width is 2*TrenchWidth+TaperWGWidth
%        PortTapers      % Cell array with all corner coordinates of all port waveguide tapers: rows: taper (common port is LAST port (app.PortDirections(app.NoOfChannels+1)) ), column: just one containing the taper, cell: coordinates of the corners [x1 y1; x2 y2; x3 y3; x4 y4]
        PortTapers      % Cell array with all corner coordinates of all port waveguide tapers: rows: taper (common port is LAST port (app.PortDirections(app.NoOfChannels+1)) ), column: just one containing the taper, cell: coordinates of the corners [x1 x2 x3 x4 x5 x6; y1 y2 y3 y4 y5 y6]
        PortTrenches    % Cell array with all corner coordinates of all port trenches: rows: trench(common port is LAST port (app.PortDirections(app.NoOfChannels+1)) ), column: just one containing the trench, cell: coordinates of the corners [x1 x2 x3 x4 x5 x6; y1 y2 y3 y4 y5 y6]
        PortTrenchesCon % Array with boundary coordinates for cladding (Trenches) of port waveguides. Consolidated from port trenches. 1. row x-values, 2. row y-values
        
        SiBoundary      % Array with coordinates of SI boundary, used for COMSOL model (and layout); 1. row x-values, 2. row y-values
        SiBoundaryGDS2  % Array with coordinates of SI boundary for GDS2 layout, used for COMSOL model (and layout); 1. row x-values, 2. row y-values
        SiBoundaryGDS2T % Array with coordinates of SI boundary, merged with port tapers for GDS2 layout; 1. row x-values, 2. row y-values
        %SiO2SlabEnhancement % Value in nanometers to shift the borders of the SiBoundaryGDS2 out for its cladding
        SiO2BoundaryGDS2 % Array with coordinates of SiO2 boundary for GDS2 layout; 1. row x-values, 2. row y-values
        
        AreaSiBoundary  % Area of SiBoundary polygon in µm²
        AreaSiBoundaryGDS2  % % Area of SiBoundaryGDS2 polygon in µm²
        AreaSiO2BoundaryGDS2  % % Area of SiO2BoundaryGDS2 polygon in µm²

        AttenuatorLeft  % Array with coordinates of left attenuator boundary, used for GDSII layout; 1. row x-values, 2. row y-values
        AttenuatorRight % Array with coordinates of right attenuator boundary, used for GDSII layout; 1. row x-values, 2. row y-values
        
        EchelleID       % Text: Name or ID of Echelle grating for GDS2 layout
        GDSWaveguideLength  % length of additional waveguide for GDSII layout to assemble all channels, also used for Comsol model if s-bends are activated
        GDSWaveguidePitch   % pitch of additional waveguides for GDSII layout to assemble all channels, also used for Comsol model if s-bends are activated
        GDSWGCurveDiscretization    % Number of linear sections to approximate the S-bends between port tapers and well-ordered connection points, also used for Comsol model if s-bends are activated
        GDS2DBUnit      % database unit of the GDS file (coordinate discretization); for 1e-9 all coordinates are snapped to a nm grid
        GDS2UserUnit    % user unit of the GDS file     (unit for any coordinate descriptions); for 1e-6 all coordinates are in µm
        GDS2AttBevelAngle   % minimum corner angle below which beveling for attenuators is done
        GDS2AttBevel        % side length of bevels for attenuators in nm
        GDS2Bevel        % side length of bevels for everything else in nm
        
        ComsolMinWavelength % Minimum wavelength for COMSOL simulations, defines mesh size and sweep start
        ComsolMaxWavelength % Maximum wavelength for COMSOL simulations, defines sweep end
        ComsolWavelengthStep    % Step size for COMSOL wavelength sweep
        ComsolNoOfPartitions    % The number of domains into which the free space area is divided for (much) faster meshing in Comsol
        ComsolBorder        % defines, which border should be used for Comsol model: 1: Full border as in GDS file, 0: reduced border (reduced width at waveguide ports)
        ComsolSBends        % defines, if S-Bend-waveguides are attached to the Comsol model or not
        
        BlockModeTrenchwidth    % Trench or cladding width if Block mode is active; special for processes, which need the cladding a certain amount larger than the core
    end
    
    methods (Access = private)
        
        function NewTable = CreateChannelPropsTable(app,NoOfChannels)
            NewTable=table([1:NoOfChannels]',zeros(NoOfChannels,1),zeros(NoOfChannels,1),ones(NoOfChannels,1)*1550,zeros(NoOfChannels,1));
            NewTable.Properties.VariableNames={'ChannelNo','X','Y','Wavelength','ERI'};
            NewTable.Properties.VariableDescriptions={'Channel No.','X (µm)','Y (µm)','Wavelength (nm)','Effective Refractive Index'};
        end
        
        function UpdateTableColoring(app)
            addStyle(app.ChannelPropsUITable,uistyle('BackgroundColor',[1 1 1]));
            addStyle(app.ChannelPropsUITable,uistyle('BackgroundColor',[0.8 1.0 0.8]),'cell',[app.TSPChannels(1), 1]);
            addStyle(app.ChannelPropsUITable,uistyle('BackgroundColor',[0.8 1.0 1.0]),'cell',[app.TSPChannels(2), 1]);
            if app.TSPChannels(1)==app.TSPChannels(2)
                addStyle(app.ChannelPropsUITable,uistyle('BackgroundColor',[1.0 0.2 0.2]),'cell',[app.TSPChannels(1), 1]);
            end
        end
        
        function ERIMatrix = UpdateERIMatrix(app,VariableName)
            % build new matrix/array with wavelengths in first column and
            % effective refractive index of currently choosen material
            % system in second column
            ERIMatrix=app.EffectiveRefractiveIndices{:,1};
            ERIMatrix=[ERIMatrix, app.EffectiveRefractiveIndices.(VariableName)];
        end
        
        function ERI=CalculateEffectiveRefractiveIndex(app, CurrentWavelength)
            % find indices from wavelengths in refractive index table adjacent to
            % desired wavelength
            x1 = find(app.EffectiveRefractiveIndex(:,1)<=CurrentWavelength,1,'last');
            x2 = find(app.EffectiveRefractiveIndex(:,1)>CurrentWavelength,1,'first');
            if isempty(x1)
                x1=x2;
                x2=x2+1;
                uialert(app.UIFigure,sprintf('Wavelength of %d nm too small, refractive index extrapolated and highly inaccurate',CurrentWavelength),'DANGER!',"Icon","error","Modal",true);
            end
            if isempty(x2)
                x2=x1;
                x1=x1-1;
                uialert(app.UIFigure,sprintf('Wavelength of %d nm too large, refractive index extrapolated and highly inaccurate',CurrentWavelength),'DANGER!',"Icon","error","Modal",true);
            end
            % get wavelengths to these indices
            lambda1 = app.EffectiveRefractiveIndex(x1,1);
            lambda2 = app.EffectiveRefractiveIndex(x2,1);
            % get effective refractive indices to these wavelengths
            n1 = app.EffectiveRefractiveIndex(x1,2);
            n2 = app.EffectiveRefractiveIndex(x2,2);
            % calculate the intermediate position of the desired wavelength
            d=(CurrentWavelength-lambda1)/(lambda2-lambda1);
            % linear interpolation of effective refractive index for the
            % intermediate wavelength; as the graph for the effective
            % refractive index is rather linear and the number of
            % calculated points is high, linear interpolation is sufficient
            ERI=n1+d*(n2-n1);
            if isempty(ERI)
                ERI=1;
            end
        end
        
        function TheTable = UpdateChannelERI(app,TheTable)
            % update the effective refractive index in the ChannelProps
            % table
            for (i=1:app.NoOfChannels)
                TheTable.ERI(i)=CalculateEffectiveRefractiveIndex(app, TheTable.Wavelength(i));
            end
        end
        
        function CalculateEllipseIntersections(app,myAxes)
            firstpoint=true;
            app.EllipseIntersections=[];
            firstAngle=NaN;
            MaxAngleReached=false;
            i=0;
            plot(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),"ro");
            plot(myAxes,app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)),"*","MarkerEdgeColor",[0 0.5 0]);
            plot(myAxes,app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)),"*","MarkerEdgeColor",[0 0.5 0]);
            while (MaxAngleReached==false)
%            for i=0:40
                Radius1=app.RadiusMin+i*(app.ChannelProps.Wavelength(app.TSPChannels(1))/app.ChannelProps.ERI(app.TSPChannels(1))*app.GratingOrder/2*1e-3);
                Radius2=app.RadiusMin+app.RadiusOffset2TSP+i*(app.ChannelProps.Wavelength(app.TSPChannels(2))/app.ChannelProps.ERI(app.TSPChannels(2))*app.GratingOrder/2*1e-3);
%                DrawEllipse(app.UIAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.RadiusMin+i*(app.ChannelProps.Wavelength(app.TSPChannels(1))/app.ChannelProps.ERI(app.TSPChannels(1))*app.GratingOrder*1e-3),120);
%                DrawEllipse(app.UIAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), app.RadiusMin+app.RadiusOffset2TSP+i*(app.ChannelProps.Wavelength(app.TSPChannels(2))/app.ChannelProps.ERI(app.TSPChannels(2))*app.GratingOrder*1e-3),120);
                if app.ShowEllipsesCheckBox.Value==true
                    DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), Radius1,120, [ 0.75, 0.75, 0.75]);
                    DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius2,120, [ 0.75, 0.75, 1.0]);
                end
                if ~firstpoint
                    % all subsequent intersetions must be within a 20°-wedge of the previous intersection
                    app.EllipseIntersections=[app.EllipseIntersections; EllipseIntersectionFinder2(myAxes, app.CommonProps.X(1),app.CommonProps.Y(1), app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius1, Radius2, app.EllipseIntersections(end, 4), 20, 'red')];
                    
                else
                    % First intersection must be in upper half and is the leftmost intersection, if there are more than one
                    app.EllipseIntersections=EllipseIntersectionFinder2(myAxes, app.CommonProps.X(1),app.CommonProps.Y(1), app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius1, Radius2, 90, 180, 'red');
                    firstAngle=app.EllipseIntersections(1,4);
                    if (isnan(firstAngle)==false)
                        firstpoint=false;
                    end
                end
                if (abs(firstAngle-app.EllipseIntersections(end,4))>=app.GratingOpeningAngle) || (Radius1>=app.RadiusMax) || (Radius2>=app.RadiusMax)
                    MaxAngleReached=true;
                end
                    
                i=i+1;
            end
            %app.EllipseIntersections
            %Strip NaN-values

            %if the intersection list is not empty, look for and strip NaN values
            if isempty(app.EllipseIntersections)==false
                lastrow=0;
                i=1;
                while (i<=size(app.EllipseIntersections,1)) && (isnan(app.EllipseIntersections(i,4))==false)
                    lastrow=i;
                    i=i+1;
                end
                if lastrow==0   %List doesn't start with valid number; shouldn't occure, but...
                    app.EllipseIntersections=[];
                else
                    app.EllipseIntersections=app.EllipseIntersections(1:lastrow,:); %cut out all intersection until the first NaN
                end
            end
            %app.EllipseIntersections
            
        end


        
        function CalculateChannelPositions(app,myAxes)
            app.TSPChannels(1);
            app.TSPChannels(2);
            app.NoOfChannels;
            app.RadiusOffset2TSP;
            
            d2 = uiprogressdlg(app.UIFigure,'Title','Calculating Channel Positions','Message','Channel No.');
            d2.Value = 0; 
            %d2.Message = 'Calculating channel positions';

            
            %calculate the radius offset for each single channel
            RadiusOffset=zeros(app.NoOfChannels,1);
            for i=1:app.NoOfChannels
                RadiusOffset(i)=(app.RadiusOffset2TSP/(app.TSPChannels(2)-app.TSPChannels(1)))*(i-app.TSPChannels(1));
            end
            %RadiusOffset
            
            %create array with the number channels colums and the number of grating points rows to store all radii for
            %calculating the remaining ports
            Radius=zeros( size(app.EllipseIntersections,1), app.NoOfChannels);
            %Radius for each single channel is
            for chan=1:app.NoOfChannels
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                    for intersection=1:size(app.EllipseIntersections,1)
                        Radius(intersection,chan)=app.RadiusMin+RadiusOffset(chan)+ (intersection-1)*(app.ChannelProps.Wavelength(chan)/app.ChannelProps.ERI(chan)*app.GratingOrder/2*1e-3);
                    end
                end
            end
            %Radius
            
            %Calculate the remaining pathlengths from the intersection points to the respective ports now
            Pathlength=zeros( size(app.EllipseIntersections,1), app.NoOfChannels);
            for chan=1:app.NoOfChannels
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                    for intersection=1:size(app.EllipseIntersections,1)
                        Pathlength(intersection,chan)=Radius(intersection,chan)*2-sqrt((app.EllipseIntersections(intersection,1)-app.CommonProps.X(1))^2+(app.EllipseIntersections(intersection,2)-app.CommonProps.Y(1))^2);
                    end
                end
            end
            %Pathlength            
            
            if app.ShowCirclesCheckBox.Value==true
                for chan=1:app.NoOfChannels
                    if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                        for intersection=1:size(app.EllipseIntersections,1)
                            app.plotcircle(myAxes, app.EllipseIntersections(intersection,1), app.EllipseIntersections(intersection,2), Pathlength(intersection,chan),'black');
                        end
                    end
                end
            end

            
            
            plot(myAxes,[app.PortSearchBoundingBox.X(1) app.PortSearchBoundingBox.X(2) app.PortSearchBoundingBox.X(2) app.PortSearchBoundingBox.X(1) app.PortSearchBoundingBox.X(1)],[app.PortSearchBoundingBox.Y(1) app.PortSearchBoundingBox.Y(1) app.PortSearchBoundingBox.Y(2) app.PortSearchBoundingBox.Y(2) app.PortSearchBoundingBox.Y(1)],'Color',[0.5 , 0.5, 0.5],'LineStyle',':');

%            fprintf('\nFinding channel positions...');
            %find for each channel find the intersections of circles around grating points with radius in 'Pathlength'
            for chan=1:app.NoOfChannels
%                fprintf('\n%i',chan);
                d2.Value = chan/app.NoOfChannels;
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))

                    finalPortIntersection=[0, 0];
                    ValidIntersectionCounter=0;
                    for intersection1=1:(size(app.EllipseIntersections,1)-1)
%                        fprintf('.');
                        for intersection2=(intersection1+1):size(app.EllipseIntersections,1)
%                            fprintf('.');
                            CurrentIntersections=CircleIntersectionFinder(app.EllipseIntersections(intersection1,1), app.EllipseIntersections(intersection1,2), Pathlength(intersection1,chan),  app.EllipseIntersections(intersection2,1), app.EllipseIntersections(intersection2,2), Pathlength(intersection2,chan));
                            if app.ShowPortSearchIntersectionsCheckBox.Value==true
                                plot(myAxes,CurrentIntersections(:,1),CurrentIntersections(:,2),"d","MarkerEdgeColor",[0 0.5 0]);
                            end
                            
                            %Check, if and which of the IntersectionsAre inside the port search bounding box
                            if app.isInPortSearchBoundingBox(CurrentIntersections(1,1),CurrentIntersections(1,2))==true
                                finalPortIntersection=finalPortIntersection+CurrentIntersections(1,:);
                                ValidIntersectionCounter=ValidIntersectionCounter+1;
                            end
                            if app.isInPortSearchBoundingBox(CurrentIntersections(2,1),CurrentIntersections(2,2))==true
                                finalPortIntersection=finalPortIntersection+CurrentIntersections(2,:);
                                ValidIntersectionCounter=ValidIntersectionCounter+1;
                            end
                        end
                    end
                    if ValidIntersectionCounter>0
                        finalPortIntersection=finalPortIntersection/ValidIntersectionCounter;
                        %finalPortIntersection
                        plot(myAxes,finalPortIntersection(1),finalPortIntersection(2),'*',"MarkerEdgeColor",[1 0 0]);
                        app.ChannelProps.X(chan)=finalPortIntersection(1);
                        app.ChannelProps.Y(chan)=finalPortIntersection(2);
                        app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(1:app.NoOfChannels,:));
                    else
                        uialert(app.UIFigure,sprintf('No valid port position found for channel no. %d',chan),'Warning!',"Icon","warning","Modal",true);
                    end
                end
            end
            fprintf('\n');
            
            close(d2);
%            CircleIntersections
%            circinter=CircleIntersectionFinder(app.EllipseIntersections(intersection1,1), app.EllipseIntersections(intersection1,2), Pathlength(intersection1,chan),  app.EllipseIntersections(intersection2,1), app.EllipseIntersections(intersection2,2), Pathlength(intersection2,chan));
%            plot(myAxes,circinter(:,1),circinter(:,2),"d","MarkerEdgeColor",[0 0.5 0]);
        end
       
        
        
        function CalculateChannelPositionsPar(app,myAxes)
            app.TSPChannels(1);
            app.TSPChannels(2);
            app.NoOfChannels;
            app.RadiusOffset2TSP;
            
            d2 = uiprogressdlg(app.UIFigure,'Title','Calculating Channel Positions','Message','Channel No.');
            d2.Value = 0; 
            %d2.Message = 'Calculating channel positions';

            
            %calculate the radius offset for each single channel
            RadiusOffset=zeros(app.NoOfChannels,1);
            for i=1:app.NoOfChannels
                RadiusOffset(i)=(app.RadiusOffset2TSP/(app.TSPChannels(2)-app.TSPChannels(1)))*(i-app.TSPChannels(1));
            end
            %RadiusOffset
            
            %create array with the number channels colums and the number of grating points rows to store all radii for
            %calculating the remaining ports
            Radius=zeros( size(app.EllipseIntersections,1), app.NoOfChannels);
            %Radius for each single channel is
            for chan=1:app.NoOfChannels
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                    for intersection=1:size(app.EllipseIntersections,1)
                        Radius(intersection,chan)=app.RadiusMin+RadiusOffset(chan)+ (intersection-1)*(app.ChannelProps.Wavelength(chan)/app.ChannelProps.ERI(chan)*app.GratingOrder/2*1e-3);
                    end
                end
            end
            %Radius
            
            %Calculate the remaining pathlengths from the intersection points to the respective ports now
            Pathlength=zeros( size(app.EllipseIntersections,1), app.NoOfChannels);
            for chan=1:app.NoOfChannels
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                    for intersection=1:size(app.EllipseIntersections,1)
                        Pathlength(intersection,chan)=Radius(intersection,chan)*2-sqrt((app.EllipseIntersections(intersection,1)-app.CommonProps.X(1))^2+(app.EllipseIntersections(intersection,2)-app.CommonProps.Y(1))^2);
                    end
                end
            end
            %Pathlength            
            
            if app.ShowCirclesCheckBox.Value==true
                for chan=1:app.NoOfChannels
                    if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))
                        for intersection=1:size(app.EllipseIntersections,1)
                            app.plotcircle(myAxes, app.EllipseIntersections(intersection,1), app.EllipseIntersections(intersection,2), Pathlength(intersection,chan),'black');
                        end
                    end
                end
            end

            
            
            plot(myAxes,[app.PortSearchBoundingBox.X(1) app.PortSearchBoundingBox.X(2) app.PortSearchBoundingBox.X(2) app.PortSearchBoundingBox.X(1) app.PortSearchBoundingBox.X(1)],[app.PortSearchBoundingBox.Y(1) app.PortSearchBoundingBox.Y(1) app.PortSearchBoundingBox.Y(2) app.PortSearchBoundingBox.Y(2) app.PortSearchBoundingBox.Y(1)],'Color',[0.5 , 0.5, 0.5],'LineStyle',':');

%            fprintf('\nFinding channel positions...');
            %find for each channel find the intersections of circles around grating points with radius in 'Pathlength'
            for chan=1:app.NoOfChannels
%                fprintf('\n%i',chan);
                d2.Value = chan/app.NoOfChannels;
                if (chan~=app.TSPChannels(1))&&(chan~=app.TSPChannels(2))

                    finalPortIntersection=[0, 0];
                    ValidIntersectionCounter=0;
                    
                    % Vorbereitungen für parfor-Loop
                    appShowPortSearchIntersectionsCheckBoxValue=app.ShowPortSearchIntersectionsCheckBox.Value;
                    appEllipseIntersections=app.EllipseIntersections;
                    appPortSearchBoundingBox=app.PortSearchBoundingBox;
                    BBxmin=app.PortSearchBoundingBox.X(1);
                    BBxmax=app.PortSearchBoundingBox.X(2);
                    BBymin=app.PortSearchBoundingBox.Y(1);
                    BBymax=app.PortSearchBoundingBox.Y(2);
                    if (BBxmin>BBxmax)
                        help=BBxmax;
                        BBxmax=BBxmin;
                        BBxmin=help;
                    end
                    if (BBymin>BBymax)
                        help=BBymax;
                        BBymax=BBymin;
                        BBymin=help;
                    end

                    
                    parfor intersection1=1:(size(appEllipseIntersections,1)-1)
%                        fprintf('.');
                        for intersection2=(intersection1+1):size(appEllipseIntersections,1)
%                            fprintf('.');
                            CurrentIntersections=CircleIntersectionFinder(appEllipseIntersections(intersection1,1), appEllipseIntersections(intersection1,2), Pathlength(intersection1,chan),  appEllipseIntersections(intersection2,1), appEllipseIntersections(intersection2,2), Pathlength(intersection2,chan));
%%                            if appShowPortSearchIntersectionsCheckBoxValue==true
%%                                plot(myAxes,CurrentIntersections(:,1),CurrentIntersections(:,2),"d","MarkerEdgeColor",[0 0.5 0]);
%%                            end
                            
                            %Check, if and which of the IntersectionsAre inside the port search bounding box
                            %if isInPortSearchBoundingBoxPar(CurrentIntersections(1,1),CurrentIntersections(1,2),appPortSearchBoundingBox)==true
                            if ((CurrentIntersections(1,1)>=BBxmin)&&(CurrentIntersections(1,1)<=BBxmax)&&(CurrentIntersections(1,2)>=BBymin)&&(CurrentIntersections(1,2)<=BBymax)==true)
                                finalPortIntersection=finalPortIntersection+CurrentIntersections(1,:);
                                ValidIntersectionCounter=ValidIntersectionCounter+1;
                            end
                            %if isInPortSearchBoundingBoxPar(CurrentIntersections(2,1),CurrentIntersections(2,2),appPortSearchBoundingBox)==true
                            if ((CurrentIntersections(2,1)>=BBxmin)&&(CurrentIntersections(2,1)<=BBxmax)&&(CurrentIntersections(2,2)>=BBymin)&&(CurrentIntersections(2,2)<=BBymax)==true)
                                finalPortIntersection=finalPortIntersection+CurrentIntersections(2,:);
                                ValidIntersectionCounter=ValidIntersectionCounter+1;
                            end
                        end
                    end
                    
                    
                    if ValidIntersectionCounter>0
                        finalPortIntersection=finalPortIntersection/ValidIntersectionCounter;
                        %finalPortIntersection
                        plot(myAxes,finalPortIntersection(1),finalPortIntersection(2),'*',"MarkerEdgeColor",[1 0 0]);
                        app.ChannelProps.X(chan)=finalPortIntersection(1);
                        app.ChannelProps.Y(chan)=finalPortIntersection(2);
                        app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(1:app.NoOfChannels,:));
                    else
                        uialert(app.UIFigure,sprintf('No valid port position found for channel no. %d',chan),'Warning!',"Icon","warning","Modal",true);
                    end
                end
            end
            fprintf('\n');
            
            close(d2);
%            CircleIntersections
%            circinter=CircleIntersectionFinder(app.EllipseIntersections(intersection1,1), app.EllipseIntersections(intersection1,2), Pathlength(intersection1,chan),  app.EllipseIntersections(intersection2,1), app.EllipseIntersections(intersection2,2), Pathlength(intersection2,chan));
%            plot(myAxes,circinter(:,1),circinter(:,2),"d","MarkerEdgeColor",[0 0.5 0]);
        end
       
        
        
        function plotcircle(app, myAxes, x,y,r,color)
            th = 0:pi/10000:2*pi;
            f = r * exp(j*th) + x+j*y;
            plot(myAxes, real(f), imag(f),'Color', color, 'Linestyle','-');
        end
        
        
        function result = isInPortSearchBoundingBox(app,x,y)
            xmin=app.PortSearchBoundingBox.X(1);
            xmax=app.PortSearchBoundingBox.X(2);
            ymin=app.PortSearchBoundingBox.Y(1);
            ymax=app.PortSearchBoundingBox.Y(2);
            if (xmin>xmax)
                help=xmax;
                xmax=xmin;
                xmin=help;
            end
            if (ymin>ymax)
                help=ymax;
                ymax=ymin;
                ymin=help;
            end
            if (x>=xmin)&&(x<=xmax)&&(y>=ymin)&&(y<=ymax)
                result=true;
            else
                result=false;
            end
        end
        
        
        

        function CalculatePortsCenterOfGravity(app,myAxes)
%            Pos=[app.CommonProps.X(1), app.CommonProps.Y(1)];
            Pos=[0, 0];
            for chan=1:app.NoOfChannels
                Pos=Pos+[app.ChannelProps.X(chan), app.ChannelProps.Y(chan)];
            end
            Pos=Pos./app.NoOfChannels;  %Center of gravity of distinct ports
            
            shiftdir=(Pos-[app.CommonProps.X(1), app.CommonProps.Y(1)])/norm( Pos-[app.CommonProps.X(1), app.CommonProps.Y(1)] );
            app.BlazeShiftXY=app.BlazeShift*shiftdir;
            
            Pos=(Pos+[app.CommonProps.X(1), app.CommonProps.Y(1)])./2;  % average of center of gravity of distinct ports and of common port
            Pos=Pos+app.BlazeShiftXY;   % shift the ports center of gravity for BlazeShift
            
            if app.ShowGratingAdditionsCheckBox.Value==true
                plot(myAxes,Pos(1),Pos(2),'x',"MarkerEdgeColor",[0.5 1 0], 'MarkerSize', 12);
            end
            app.PortsCenterOfGravity=Pos;
        end



        function CalculateGratingCenter(app, myAxes)
            app.GratingCenter=[];
            if isempty(app.EllipseIntersections)==false
                %calculate the line through the center of the grating with respect to the center of gravity of the ports
                %first shift the end points of the grating so that the ports center of gravity is the new origin
                x1=app.EllipseIntersections(1,1)-app.PortsCenterOfGravity(1);
                y1=app.EllipseIntersections(1,2)-app.PortsCenterOfGravity(2);
                x2=app.EllipseIntersections(end,1)-app.PortsCenterOfGravity(1);
                y2=app.EllipseIntersections(end,2)-app.PortsCenterOfGravity(2);
                %convert to polar coordinates
                phi1=atan2(y1,x1);
                phi2=atan2(y2,x2);
                %r1=sqrt(x1^2+y1^2);
                %r2=sqrt(x2^2+y2^2);
                phicenter=(phi1+phi2)/2;
                cosphicenter=cos(phicenter);
                sinphicenter=sin(phicenter);
                %calculate distance from each grating point to line through port center of gravity with angle phicenter
                %first calculate line
                %two points on the line:
                %centerline_x1=app.PortsCenterOfGravity(1);
                %centerline_y1=app.PortsCenterOfGravity(2);
                %centerline_x2=cos(phicenter)+app.PortsCenterOfGravity(1);
                %centerline_y2=sin(phicenter)+app.PortsCenterOfGravity(2);
                %line formula
                %gvec=posvec+t*dirvec
                %gvec=[centerline_x1; centerline_y1]+t*([centerline_x2; centerline_y2]-[centerline_x1; centerline_y1])
                %gvec=[centerline_x1; centerline_y1]+t*([cos(phicenter); sin(phicenter)]
                
                %distance to qvec=[app.EllipseIntersections(i,1); app.EllipseIntersections(i,2)]?
                %d=abs( dirvec cross (qvec-posvec) )/abs(dirvec)
                %abs(dirvec)=1 ==> d=abs( dirvec cross (qvec-posvec) )
                %d=abs([cos(phicenter); sin(phicenter)] * ([app.EllipseIntersections(i,1); app.EllipseIntersections(i,2)]-[centerline_x1; centerline_y1])   );
                %in 2D vector '*' is the determinant of the two vectors put into a matrix
                %d=abs( det( [dirvec , (qvec-posvec)] ) )
                d=1e6;  %set initial d to large value
                gratingpointnumber=0;
                for i=1:size(app.EllipseIntersections,1)
                    dnew=abs(det([[cosphicenter; sinphicenter],[app.EllipseIntersections(i,1); app.EllipseIntersections(i,2)]-[app.PortsCenterOfGravity(1); app.PortsCenterOfGravity(2)]])   );    
                    if dnew<d
                        d=dnew; %set new lowest distance
                        gratingpointnumber=i;   %and save the respective grating point for that lowest distance
                    end
                end
                %get the distance of the nearest grating point to the ports center of gravity
                dx=app.EllipseIntersections(gratingpointnumber,1)-app.PortsCenterOfGravity(1);
                dy=app.EllipseIntersections(gratingpointnumber,2)-app.PortsCenterOfGravity(2);
                r=sqrt(dx^2+dy^2);
                %now get the point on the line with distance r to ports center of gravity
                
                %x=r*cosphicenter+app.PortsCenterOfGravity(1);
                %y=r*sinphicenter+app.PortsCenterOfGravity(2);
                %app.GratingCenter=[x, y];
                app.GratingCenter=[r*cosphicenter+app.PortsCenterOfGravity(1), r*sinphicenter+app.PortsCenterOfGravity(2)];
            end
            if app.ShowGratingAdditionsCheckBox.Value==true
                plot(myAxes,app.GratingCenter(1),app.GratingCenter(2),'x',"MarkerEdgeColor",[0.5 1 0], 'MarkerSize', 12);
                
                plot(myAxes,[app.PortsCenterOfGravity(1); app.GratingCenter(1)], [app.PortsCenterOfGravity(2); app.GratingCenter(2)], "r-");
                plot(myAxes,[app.PortsCenterOfGravity(1); app.EllipseIntersections(1,1)], [app.PortsCenterOfGravity(2); app.EllipseIntersections(1,2)], "r-");
                plot(myAxes,[app.PortsCenterOfGravity(1); app.EllipseIntersections(end,1)], [app.PortsCenterOfGravity(2); app.EllipseIntersections(end,2)], "r-");

            end
        end
        
        
        
        function CalculateGratingReflectorDirections(app, myAxes)
            if isempty(app.EllipseIntersections)==false
                app.GratingPointsReflectorDirections=zeros(size(app.EllipseIntersections,1),1);
                for i=1:size(app.EllipseIntersections,1)
                    app.GratingPointsReflectorDirections(i)=atan2(app.EllipseIntersections(i,2)-app.PortsCenterOfGravity(2) , app.EllipseIntersections(i,1)-app.PortsCenterOfGravity(1));
                    if app.ShowGratingAdditionsCheckBox.Value==true
                        plot(myAxes,[app.EllipseIntersections(i,1); app.EllipseIntersections(i,1)+10*cos(app.GratingPointsReflectorDirections(i))], [app.EllipseIntersections(i,2); app.EllipseIntersections(i,2)+10*sin(app.GratingPointsReflectorDirections(i))], 'g-');
                    end
                end
            end
            if (app.GratingPointsReflectorDirections(end)-app.GratingPointsReflectorDirections(1))<0
                app.GratingTrajectoryReversed=0;
            else
                app.GratingTrajectoryReversed=1;
            end
        end
        
        
        
        function CalculateGratingReflectors(app, myAxes)
            if isempty(app.EllipseIntersections)==false
                
                app.BraggReflectors=cell(size(app.EllipseIntersections,1),app.BraggPeriods);
                app.SiBoundary=zeros(2,size(app.EllipseIntersections,1));
                app.SiO2BoundaryGDS2=app.SiBoundary;

                % first calculate the angles of the lines between grating points, which define the borders of the bragg
                % gratings
                app.GratingPointsReflectorHalfAngles=zeros(size(app.EllipseIntersections,1)+1,1);
                                
                for i=2:size(app.EllipseIntersections,1)
                    app.GratingPointsReflectorHalfAngles(i)=(atan2(app.EllipseIntersections(i-1,2)-app.PortsCenterOfGravity(2) , app.EllipseIntersections(i-1,1)-app.PortsCenterOfGravity(1)) + atan2(app.EllipseIntersections(i,2)-app.PortsCenterOfGravity(2) , app.EllipseIntersections(i,1)-app.PortsCenterOfGravity(1)) )/2;
                end
                app.GratingPointsReflectorHalfAngles(1)= 2*atan2(app.EllipseIntersections(1,2)-app.PortsCenterOfGravity(2) , app.EllipseIntersections(1,1)-app.PortsCenterOfGravity(1))-app.GratingPointsReflectorHalfAngles(2);
                app.GratingPointsReflectorHalfAngles(end)= 2*atan2(app.EllipseIntersections(size(app.EllipseIntersections,1),2)-app.PortsCenterOfGravity(2) , app.EllipseIntersections(size(app.EllipseIntersections,1),1)-app.PortsCenterOfGravity(1)) - app.GratingPointsReflectorHalfAngles(end-1);
                %app.GratingPointsReflectorHalfAngles
%{                
                if app.ShowGratingAdditionsCheckBox.Value==true
                    for i=1:size(app.EllipseIntersections,1)+1
                        plot(myAxes,[app.PortsCenterOfGravity(1)+app.RadiusMin*cos(app.GratingPointsReflectorHalfAngles(i)); app.PortsCenterOfGravity(1)+app.RadiusMax*cos(app.GratingPointsReflectorHalfAngles(i))], [app.PortsCenterOfGravity(2)+app.RadiusMin*sin(app.GratingPointsReflectorHalfAngles(i)); app.PortsCenterOfGravity(2)+app.RadiusMax*sin(app.GratingPointsReflectorHalfAngles(i))], 'b-');
                    end
                end
%}              
                
                % calculate the corner points of all Bragg elements and store them into app.BraggReflectors
                
                % vector representation of border lines
                % bvec=[app.PortsCenterOfGravity(1);app.PortsCenterOfGravity(2)] + s*[cos(app.GratingPointsReflectorHalfAngles(i));sin(app.GratingPointsReflectorHalfAngles(i))]
                
                % points, which define the bragg grating, beginning with grating points
                % pvec=[app.EllipseIntersections(j,1);app.EllipseIntersections(j,2)] + t*[cos(app.GratingPointsReflectorDirections(j)); sin(app.GratingPointsReflectorDirections(j))]
                % t=0, app.BraggSiO2Length*1e-3, app.BraggPeriodLength*1e-3, (app.BraggPeriodLength+app.BraggSiO2Length)*1e-3,
                % 2*app.BraggPeriodLength*1e-3, ...
                
                % vector representation of first bragg grating line through grating points
                % qvec=pvec + t*[sin(app.GratingPointsReflectorDirections(j));-cos(app.GratingPointsReflectorDirections(j))]
                

                % part of line definitions for left and right boundaries
                pos1=[app.PortsCenterOfGravity(1);app.PortsCenterOfGravity(2)];
                for i=1:size(app.EllipseIntersections,1)
                    % second part of line definitions for left and right boundaries
                    dir1a=[cos(app.GratingPointsReflectorHalfAngles(i));sin(app.GratingPointsReflectorHalfAngles(i))];
                    dir1b=[cos(app.GratingPointsReflectorHalfAngles(i+1));sin(app.GratingPointsReflectorHalfAngles(i+1))];
                    % part of line definitions for bottom boundary
                    dir2=[sin(app.GratingPointsReflectorDirections(i));-cos(app.GratingPointsReflectorDirections(i))];
                    for j=1:app.BraggPeriods
                        % second part of line definitions for bottom boundary
                        pos2=[app.EllipseIntersections(i,1)+((j-1)*app.BraggPeriodLength)*1e-3*cos(app.GratingPointsReflectorDirections(i));
                              app.EllipseIntersections(i,2)+((j-1)*app.BraggPeriodLength)*1e-3*sin(app.GratingPointsReflectorDirections(i))];
                        % additional line definition for top boundary (direction vector is the same as before: dir2
                        pos3=[app.EllipseIntersections(i,1)+((j-1)*app.BraggPeriodLength +app.BraggSiO2Length)*1e-3*cos(app.GratingPointsReflectorDirections(i));
                              app.EllipseIntersections(i,2)+((j-1)*app.BraggPeriodLength +app.BraggSiO2Length)*1e-3*sin(app.GratingPointsReflectorDirections(i))];
                        % find the corner points of the Bragg element of the current period...
                        %   bottom left
                        bp1=app.LineLineIntersection(pos1 , dir1a , pos2 , dir2 );
                        %   bottom right
                        bp2=app.LineLineIntersection(pos1 , dir1b , pos2 , dir2 );
                        %   top left
                        bp3=app.LineLineIntersection(pos1 , dir1a , pos3 , dir2 );
                        %   top right
                        bp4=app.LineLineIntersection(pos1 , dir1b , pos3 , dir2 );

                        %to avoid sharp angles and small features error in certain design rule checks, the calculated points are shifted a little to make the reflectors smaller and give them some space in between 
                        newbp1=app.MovePoint(bp1, bp2, app.ReflectorSpacing*1E-3/2);
                        newbp2=app.MovePoint(bp2, bp1, app.ReflectorSpacing*1E-3/2);
                        newbp3=app.MovePoint(bp3, bp4, app.ReflectorSpacing*1E-3/2);
                        newbp4=app.MovePoint(bp4, bp3, app.ReflectorSpacing*1E-3/2);
                        
                        % ...and put the coordinates into the respective cell
%                        app.BraggReflectors{i,j}=[ bp1(1), bp1(2); bp2(1), bp2(2); bp4(1), bp4(2); bp3(1), bp3(2)];
                        app.BraggReflectors{i,j}=[ newbp1(1), newbp2(1), newbp4(1), newbp3(1); newbp1(2), newbp2(2), newbp4(2), newbp3(2)];

                        % let's show the stuff, if the user wants
                        if app.ShowGratingAdditionsCheckBox.Value==true
                            plot(myAxes,[newbp1(1);newbp2(1);newbp4(1);newbp3(1);newbp1(1)],[newbp1(2);newbp2(2);newbp4(2);newbp3(2);newbp1(2)], 'b-');
                            %plot(myAxes,[app.BraggReflectors{i,j}(:,1); app.BraggReflectors{i,j}(1,1)], [app.BraggReflectors{i,j}(:,2); app.BraggReflectors{i,j}(1,2)], 'b-');
                        end
                        %plot(myAxes, bp3(1), bp3(2), "rx");
                    end
                    % and calculate some border points
                    app.SiBoundary(:,i)=[app.EllipseIntersections(i,1)+(app.BraggPeriodLength*app.BraggPeriods+app.BoundaryDistance)*1e-3*cos(app.GratingPointsReflectorDirections(i));
                                         app.EllipseIntersections(i,2)+(app.BraggPeriodLength*app.BraggPeriods+app.BoundaryDistance)*1e-3*sin(app.GratingPointsReflectorDirections(i))];
%                    app.SiO2BoundaryGDS2(:,i)=app.SiBoundary(:,i)+[app.SiO2SlabEnhancement*1e-3*cos(app.GratingPointsReflectorDirections(i)) ; app.SiO2SlabEnhancement*1e-3*sin(app.GratingPointsReflectorDirections(i))]; %shift the border for the cladding even further out
                    app.SiO2BoundaryGDS2(:,i)=app.SiBoundary(:,i);
                end
                %celldisp(app.BraggReflectors)
                
                
                if app.GratingTrajectoryReversed==0
                    LeftDirectionIndex=1;
                    LeftHalfAngleIndex=1;
                    LeftReflectorIndex=1;
                    RightDirectionIndex=size(app.GratingPointsReflectorDirections,1);
                    RightHalfAngleIndex=size(app.GratingPointsReflectorHalfAngles,1);
                    RightReflectorIndex=size(app.BraggReflectors,1);
                    
                else
                    LeftDirectionIndex=size(app.GratingPointsReflectorDirections,1);
                    LeftHalfAngleIndex=size(app.GratingPointsReflectorHalfAngles,1);
                    LeftReflectorIndex=size(app.BraggReflectors,1);
                    RightDirectionIndex=1;
                    RightHalfAngleIndex=1;
                    RightReflectorIndex=1;
                    app.SiBoundary=flip(app.SiBoundary,2); % reverse order of path of border points
                    app.SiO2BoundaryGDS2=flip(app.SiO2BoundaryGDS2,2);
                end
                
                % additional most 'top left' point of boundary around bragg grating
                dir1a=[sin(app.GratingPointsReflectorDirections(LeftDirectionIndex));-cos(app.GratingPointsReflectorDirections(LeftDirectionIndex))];
                dir1b=[cos(app.GratingPointsReflectorHalfAngles(LeftHalfAngleIndex));sin(app.GratingPointsReflectorHalfAngles(LeftHalfAngleIndex))];
                pos1=[app.BraggReflectors{LeftReflectorIndex,end}(1,4); app.BraggReflectors{LeftReflectorIndex,end}(2,4)]-dir1a*app.BoundaryDistance*1e-3+dir1b*app.BoundaryDistance*1e-3;
                %pos1SiO2=[app.BraggReflectors{LeftReflectorIndex,end}(1,4); app.BraggReflectors{LeftReflectorIndex,end}(2,4)]-dir1a*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3+dir1b*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3;
                pos1SiO2=[app.BraggReflectors{LeftReflectorIndex,end}(1,4); app.BraggReflectors{LeftReflectorIndex,end}(2,4)]-dir1a*(app.BoundaryDistance)*1e-3+dir1b*(app.BoundaryDistance)*1e-3;
                
                % additional most 'bottom left' point of boundary around bragg grating
                dir2=[sin(app.GratingPointsReflectorDirections(LeftDirectionIndex));-cos(app.GratingPointsReflectorDirections(LeftDirectionIndex))];
                pos2=[app.BraggReflectors{LeftReflectorIndex,1}(1,1); app.BraggReflectors{LeftReflectorIndex,1}(2,1)]-dir2*app.BoundaryDistance*1e-3;
                %pos2SiO2=[app.BraggReflectors{LeftReflectorIndex,1}(1,1); app.BraggReflectors{LeftReflectorIndex,1}(2,1)]-dir2*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3;
                pos2SiO2=[app.BraggReflectors{LeftReflectorIndex,1}(1,1); app.BraggReflectors{LeftReflectorIndex,1}(2,1)]-dir2*(app.BoundaryDistance)*1e-3;

                % additional most 'top right' point of boundary around bragg grating
                dir3a=[sin(app.GratingPointsReflectorDirections(RightDirectionIndex));-cos(app.GratingPointsReflectorDirections(RightDirectionIndex))];
                dir3b=[cos(app.GratingPointsReflectorHalfAngles(RightHalfAngleIndex));sin(app.GratingPointsReflectorHalfAngles(RightHalfAngleIndex))];
                pos3=[app.BraggReflectors{RightReflectorIndex,end}(1,3); app.BraggReflectors{RightReflectorIndex,end}(2,3)]+dir3a*app.BoundaryDistance*1e-3+dir3b*app.BoundaryDistance*1e-3;
                %pos3SiO2=[app.BraggReflectors{RightReflectorIndex,end}(1,3); app.BraggReflectors{RightReflectorIndex,end}(2,3)]+dir3a*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3+dir3b*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3;
                pos3SiO2=[app.BraggReflectors{RightReflectorIndex,end}(1,3); app.BraggReflectors{RightReflectorIndex,end}(2,3)]+dir3a*(app.BoundaryDistance)*1e-3+dir3b*(app.BoundaryDistance)*1e-3;
                
                % additional most 'bottom right' point of boundary around bragg grating
                dir4=[sin(app.GratingPointsReflectorDirections(RightDirectionIndex));-cos(app.GratingPointsReflectorDirections(RightDirectionIndex))];
                pos4=[app.BraggReflectors{RightReflectorIndex,1}(1,2); app.BraggReflectors{RightReflectorIndex,1}(2,2)]+dir4*app.BoundaryDistance*1e-3;
                %pos4SiO2=[app.BraggReflectors{RightReflectorIndex,1}(1,2); app.BraggReflectors{RightReflectorIndex,1}(2,2)]+dir4*(app.BoundaryDistance+app.SiO2SlabEnhancement)*1e-3;
                pos4SiO2=[app.BraggReflectors{RightReflectorIndex,1}(1,2); app.BraggReflectors{RightReflectorIndex,1}(2,2)]+dir4*(app.BoundaryDistance)*1e-3;

                
%{                
                % additional most 'top left' point of boundary around bragg grating
                dir1a=[sin(app.GratingPointsReflectorDirections(1));-cos(app.GratingPointsReflectorDirections(1))];
                dir1b=[cos(app.GratingPointsReflectorHalfAngles(1));sin(app.GratingPointsReflectorHalfAngles(1))];
                pos1=[app.BraggReflectors{1,end}(1,4); app.BraggReflectors{1,end}(2,4)]-dir1a*app.BoundaryDistance*1e-3+dir1b*app.BoundaryDistance*1e-3;
                
                % additional most 'bottom left' point of boundary around bragg grating
                dir2=[sin(app.GratingPointsReflectorDirections(1));-cos(app.GratingPointsReflectorDirections(1))];
                pos2=[app.BraggReflectors{1,1}(1,1); app.BraggReflectors{1,1}(2,1)]-dir2*app.BoundaryDistance*1e-3;

                % additional most 'top right' point of boundary around bragg grating
                dir3a=[sin(app.GratingPointsReflectorDirections(end));-cos(app.GratingPointsReflectorDirections(end))];
                dir3b=[cos(app.GratingPointsReflectorHalfAngles(end));sin(app.GratingPointsReflectorHalfAngles(end))];
                pos3=[app.BraggReflectors{end,end}(1,3); app.BraggReflectors{end,end}(2,3)]+dir3a*app.BoundaryDistance*1e-3+dir3b*app.BoundaryDistance*1e-3;
                
                % additional most 'bottom right' point of boundary around bragg grating
                dir4=[sin(app.GratingPointsReflectorDirections(end));-cos(app.GratingPointsReflectorDirections(end))];
                pos4=[app.BraggReflectors{end,1}(1,2); app.BraggReflectors{end,1}(2,2)]+dir4*app.BoundaryDistance*1e-3;
%}                
                app.SiBoundary=[pos2 pos1 app.SiBoundary(:,:) pos3 pos4];
                app.SiO2BoundaryGDS2=[pos2SiO2 pos1SiO2 app.SiO2BoundaryGDS2(:,:) pos3SiO2 pos4SiO2];
                
                % use some points of the Si-boundary for the top parts of the attenuators in the GDS2 layout...
                app.AttenuatorLeft = [ pos2 pos1 app.SiBoundary(:,3)*0.8+pos1*0.2 ];        % fixed to 20% of the way between the respective points
                app.AttenuatorRight = [ app.SiBoundary(:,end-2)*0.8+pos3*0.2 pos3 pos4 ];   % fixed to 20% of the way between the respective points
                
                %plot(myAxes,app.BraggReflectors{1,end}(1,4),app.BraggReflectors{1,end}(2,4),'rd');
                
%                if app.ShowGratingAdditionsCheckBox.Value==true
%                    plot(myAxes,app.SiBoundary(1,:),app.SiBoundary(2,:), 'b-');
%                end
            end
        end
        
        
        function ipoint = LineLineIntersection(app, p1, p1dir, p3, p3dir)
            ipoint=NaN(2,1);
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

        
        
        
        function ipoint = MovePoint(app, p1, p2, dist)
            % finds a new point on straight between p1 and p2 with distance dist to p1
            ipoint=NaN(2,1);
            dir=(p2-p1)/norm(p2-p1);
            ipoint=p1+dir*dist;
        end

        
        
        
        function CalculatePortDirections(app, myAxes)
            if isempty(app.GratingCenter)==false
                app.PortDirections=zeros(app.NoOfChannels+1,1);
                for i=1:app.NoOfChannels
                    app.PortDirections(i)=atan2(app.ChannelProps.Y(i)-app.GratingCenter(2), app.ChannelProps.X(i)-app.GratingCenter(1));
                    if app.ShowGratingAdditionsCheckBox.Value==true
                        plot(myAxes,[app.ChannelProps.X(i); app.ChannelProps.X(i)+10*cos(app.PortDirections(i))], [app.ChannelProps.Y(i); app.ChannelProps.Y(i)+10*sin(app.PortDirections(i))], 'g-');
                    end
                end
                app.PortDirections(app.NoOfChannels+1)=atan2(app.CommonProps.Y(1)-app.GratingCenter(2), app.CommonProps.X(1)-app.GratingCenter(1));
                if app.ShowGratingAdditionsCheckBox.Value==true
                    plot(myAxes,[app.CommonProps.X(1); app.CommonProps.X(1)+10*cos(app.PortDirections(app.NoOfChannels+1))], [app.CommonProps.Y(1); app.CommonProps.Y(1)+10*sin(app.PortDirections(app.NoOfChannels+1))], 'g-');
                end
                
            end
        end
        
        
        % calculate angle between lines connecting p1 to p2 to p3
        function angle=CalculateAngle(app, p1, p2, p3)
            CAdir1 = p1-p2;
            CAdir2 = p3-p2;
            %Formula from https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
            angle = atan2d(CAdir1(1)*CAdir2(2)-CAdir1(2)*CAdir2(1),CAdir1(1)*CAdir2(1)+CAdir1(2)*CAdir2(2));
        end
        
        
        
        % This function shifts the new bevel points a certain length from the corner
        % calculate bevel at p2, newPoints is matrix [p4x p5x; p4y p5y] with p4 between p1 and p2, and p5 between p2 and
        % p3 with 'distance' to p2
        function newPoints=CalculateBevel(app, p1, p2, p3, distance)
            newPoints=zeros(2,2);
            CBdir1 = p1-p2;
            CBdir2 = p3-p2;
            CBdir1 = CBdir1/norm(CBdir1);
            CBdir2 = CBdir2/norm(CBdir2);
            newPoints(:,1)=p2+CBdir1*distance/1e3;
            newPoints(:,2)=p2+CBdir2*distance/1e3;
            %newPoints
        end
        
            
        
        % This function makes the chord length of the bevel a certain length
        % calculate bevel at p2, newPoints is matrix [p4x p5x; p4y p5y] with p4 between p1 and p2, and p5 between p2 and
        % p3 with 'distance' to p2
        function newPoints=CalculateBevel2(app, p1, p2, p3, distance, angle)
            newPoints=zeros(2,2);
            CBdir1 = p1-p2;
            CBdir2 = p3-p2;
            CBdir1 = CBdir1/norm(CBdir1);
            CBdir2 = CBdir2/norm(CBdir2);
            movedistance=distance/(2*sind(angle/2));
            newPoints(:,1)=p2+CBdir1*movedistance/1e3;
            newPoints(:,2)=p2+CBdir2*movedistance/1e3;
            %newPoints
        end
        
            
        
        function newPath = BevelEdges(app, oldPath, minAngle, Bevel)
            newPath=[];
            %oldPath
            for i=1:size(oldPath,2)
                %calculate angle for i-th corner
                if (i==1) p1=oldPath(:,end);
                else      p1=oldPath(:,i-1);
                end
                
                p2=oldPath(:,i);
                
                if (i==size(oldPath,2)) p3=oldPath(:,1);
                else                    p3=oldPath(:,i+1);
                end
                
                angle=abs(app.CalculateAngle(p1,p2,p3));
                
                if (angle<minAngle)
                    % bevel edge here
                    %BevelPoints=app.CalculateBevel(p1, p2, p3, Bevel);
                    BevelPoints=app.CalculateBevel2(p1, p2, p3, Bevel, angle);
                    newPath(:,end+1)=BevelPoints(:,1);
                    newPath(:,end+1)=BevelPoints(:,2);
                else
                    newPath(:,end+1)=oldPath(:,i);
                end
            end
            %newPath
        end
        
        


        function CalculatePortTapers(app, myAxes)
            %fprintf('Trench width: %f\n', app.TrenchWidth);
            app.PortTapers=cell(app.NoOfChannels+1,1);
            app.PortTrenches=cell(app.NoOfChannels+1,1);
            app.PortTrenchesCon=[];
            app.PortWGPositions=zeros(app.NoOfChannels+1,2);
            
            portCoordinates=zeros(app.NoOfChannels+1,2);
            portCoordinates(1:end-1,1)=app.ChannelProps.X(1:app.NoOfChannels);
            portCoordinates(1:end-1,2)=app.ChannelProps.Y(1:app.NoOfChannels);
            portCoordinates(end,1)=app.CommonProps.X;
            portCoordinates(end,2)=app.CommonProps.Y;
            for i=1:app.NoOfChannels+1
                %calculate top end (free space end) of taper
                pos1=[portCoordinates(i,1) ; portCoordinates(i,2)] +app.TaperShift*1e-3*[cos(app.PortDirections(i));sin(app.PortDirections(i))];
                dir1=[sin(app.PortDirections(i)) ; -cos(app.PortDirections(i))];
                %calculate bottom end (waveguide end) of taper
                pos2=pos1+app.TaperLength*1e-3*[cos(app.PortDirections(i));sin(app.PortDirections(i))];
                %direction is the same as above: dir1

                %calculate bottom end of waveguide stub
                pos3=pos2+app.WaveguideStubLength*1e-3*[cos(app.PortDirections(i));sin(app.PortDirections(i))];
                %direction is the same as above: dir1
                
                %for the ports...
                %top left
                bp1=pos1+app.TaperWidth*1e-3/2*dir1;
                %top right
                bp2=pos1-app.TaperWidth*1e-3/2*dir1;
                
                %mid left
                bp3=pos2+app.TaperWGWidth*1e-3/2*dir1;
                %mid right
                bp4=pos2-app.TaperWGWidth*1e-3/2*dir1;
                
                %bottom left
                bp5=pos3+app.TaperWGWidth*1e-3/2*dir1;
                %bottom right
                bp6=pos3-app.TaperWGWidth*1e-3/2*dir1;
                
                
%                app.PortTapers{i}=[ bp1(1), bp1(2); bp3(1), bp3(2); bp4(1), bp4(2); bp2(1), bp2(2)];
                app.PortTapers{i}=[ bp1(1), bp3(1), bp5(1), bp6(1), bp4(1), bp2(1); bp1(2), bp3(2), bp5(2), bp6(2), bp4(2), bp2(2)];
                app.PortWGPositions(i,1)=pos3(1);
                app.PortWGPositions(i,2)=pos3(2);
                
                %... and for the Trenches
                %top left
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels+1)))
                    tp1=pos1+(app.TaperWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp1=pos1+(app.TaperWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
                %top right
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels)))
                    tp2=pos1-(app.TaperWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp2=pos1-(app.TaperWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
                
                %mid left
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels+1)))
                    tp3=pos2+(app.TaperWGWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp3=pos2+(app.TaperWGWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
                %mid right
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels)))
                    tp4=pos2-(app.TaperWGWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp4=pos2-(app.TaperWGWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
                
                %bottom left
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels+1)))
                    tp5=pos3+(app.TaperWGWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp5=pos3+(app.TaperWGWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
                %bottom right
                if ((app.CladdingBlockModeCheckBox.Value==true) && (i==(app.NoOfChannels)))
                    tp6=pos3-(app.TaperWGWidth+2*app.BlockModeTrenchwidth)*1e-3/2*dir1;
                else
                    tp6=pos3-(app.TaperWGWidth+2*app.TrenchWidth)*1e-3/2*dir1;
                end
%                app.PortTrenches{i}=[ tp1(1), tp3(1), tp5(1), tp6(1), tp4(1), tp2(1); tp1(2), tp3(2), tp5(2), tp6(2), tp4(2), tp2(2)];
                app.PortTrenches{i}=[ tp1(1), tp3(1), tp5(1), tp6(1), tp4(1), tp2(1), bp2(1), bp1(1); tp1(2), tp3(2), tp5(2), tp6(2), tp4(2), tp2(2), bp2(2), bp1(2)];
                
                
                
%                if app.ShowGratingAdditionsCheckBox.Value==true
%                    plot(myAxes,[bp1(1); bp3(1); bp5(1); bp6(1); bp4(1); bp2(1); bp1(1)], [bp1(2); bp3(2); bp5(2); bp6(2); bp4(2); bp2(2); bp1(2)], 'b-');
%                    plot(myAxes,pos3(1), pos3(2), 'rx');
%                end
            end
            
%{
            %modify trenches between taper ends
            %The common port stays as it is, there must be enough distance to the first port, which has to be on the
            %right (larger x-coordinates). The second port has to be on the right of the first port and so on
            % modify the 1. and the 6. coordinate
            for i=1:app.NoOfChannels-1
                %calculate the mid-point of the line between the waveguide taper corners of adjacent channels, which are the 7. and 8. coordinates
                midpoint=(app.PortTrenches{i}(:,7)+app.PortTrenches{i+1}(:,8))./2;
                % this mid-point becomes the new corner coordinates of the trenches
                app.PortTrenches{i}(:,6)=midpoint;
                app.PortTrenches{i+1}(:,1)=midpoint;

            end
            %and now modify the top left corner of channel 1 and put it to the top right corner of the common channel
            app.PortTrenches{1}(:,1)=app.PortTrenches{end}(:,6);
%}
            %modify trenches between taper ends
            %The first port has to be on the right (larger x-coordinates) of the common port.
            %The second port has to be on the right of the first port and so on.
            %Modify the 1. and the 6. coordinate to become the Si wavguide corner.
            for i=1:app.NoOfChannels-1
                app.PortTrenches{i}(:,6)=app.PortTrenches{i+1}(:,8);
                app.PortTrenches{i+1}(:,1)=app.PortTrenches{i}(:,7);
            end
            %The common port (and the 1. port) have to be treated seperately, due to the odd numbering and geometric
            %placing (common port is the last number, but the first port to the left)
            app.PortTrenches{end}(:,6)=app.PortTrenches{1}(:,8);
            app.PortTrenches{1}(:,1)=app.PortTrenches{end}(:,7);
            
            
            %Now calculate a single 'trench' boundary for all tapers; required for angle beveling for certain design rules
            %Start with common port
            if (app.CladdingBlockModeCheckBox.Value~=true)
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,1);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,2);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,3);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,4);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,5);
                for i=1:app.NoOfChannels
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,2);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,3);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,4);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,5);
                end
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels}(:,6);
                for i=1:app.NoOfChannels
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1-i}(:,7);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1-i}(:,8);
                end
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,7);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,8);
                app.PortTrenchesCon=app.BevelEdges(app.PortTrenchesCon, app.GDS2AttBevelAngle, app.GDS2Bevel);
                %plot(myAxes,[app.PortTrenchesCon(1,:) app.PortTrenchesCon(1,1)], [app.PortTrenchesCon(2,:) app.PortTrenchesCon(2,1)], 'k-');
            else
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,1);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,2);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,3);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,4);
                %app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,5);
                for i=1:app.NoOfChannels
                    %app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,2);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,3);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,4);
                    %app.PortTrenchesCon(:,end+1)=app.PortTrenches{i}(:,5);
                end
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels}(:,5);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels}(:,6);
                for i=1:app.NoOfChannels
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1-i}(:,7);
                    app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1-i}(:,8);
                end
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,7);
                app.PortTrenchesCon(:,end+1)=app.PortTrenches{app.NoOfChannels+1}(:,8);
                app.PortTrenchesCon=app.BevelEdges(app.PortTrenchesCon, app.GDS2AttBevelAngle, app.GDS2Bevel);
                %plot(myAxes,[app.PortTrenchesCon(1,:) app.PortTrenchesCon(1,1)], [app.PortTrenchesCon(2,:) app.PortTrenchesCon(2,1)], 'k-');
            end
            
            % make a copy of the current SiBoundary, as the closing of the polygon will be different for the Comsol and
            % the GDS2 versions, it could be modified in the final boundary, but it's easier to get the respective
            % points this way...
            app.SiBoundaryGDS2 = app.SiBoundary;

            % close the Si-boundary now
%            HelpBoundary=zeros(2,4*(app.NoOfChannels+1));
            HelpBoundary=zeros(2,2*(app.NoOfChannels+1));
            HelpBoundaryT=zeros(2,6*(app.NoOfChannels+1));  % lower boundary of free space region, including waveguide core tapers
            for i=1:app.NoOfChannels
                %start with most right point of upper part of Trenches to add the port section to the grating section
                %boundary clockwise
%                HelpBoundary(:,(i-1)*4+1)=app.PortTrenches{app.NoOfChannels-i+1}(:,6);
                HelpBoundary(:,(i-1)*2+1)=app.PortTrenches{app.NoOfChannels-i+1}(:,7);
                HelpBoundary(:,(i-1)*2+2)=app.PortTrenches{app.NoOfChannels-i+1}(:,8);
%                HelpBoundary(:,(i-1)*4+4)=app.PortTrenches{app.NoOfChannels-i+1}(:,1);

                for k=1:6
                    HelpBoundaryT(:,(i-1)*6+k)=app.PortTapers{app.NoOfChannels-i+1}(:,(6+1)-k);
                end
%                HelpBoundaryT(:,(i-1)*6+1)=app.PortTapers{app.NoOfChannels-i+1}(:,6);
%                HelpBoundaryT(:,(i-1)*6+2)=app.PortTapers{app.NoOfChannels-i+1}(:,5);
%                HelpBoundaryT(:,(i-1)*6+3)=app.PortTapers{app.NoOfChannels-i+1}(:,4);
%                HelpBoundaryT(:,(i-1)*6+4)=app.PortTapers{app.NoOfChannels-i+1}(:,3);
%                HelpBoundaryT(:,(i-1)*6+5)=app.PortTapers{app.NoOfChannels-i+1}(:,2);
%                HelpBoundaryT(:,(i-1)*6+6)=app.PortTapers{app.NoOfChannels-i+1}(:,1);
            end
            %and add the common channel
%            HelpBoundary(:,app.NoOfChannels*4+1)=app.PortTrenches{app.NoOfChannels+1}(:,6);
            HelpBoundary(:,app.NoOfChannels*2+1)=app.PortTrenches{app.NoOfChannels+1}(:,7);
            HelpBoundary(:,app.NoOfChannels*2+2)=app.PortTrenches{app.NoOfChannels+1}(:,8);
%            HelpBoundary(:,app.NoOfChannels*4+4)=app.PortTrenches{app.NoOfChannels+1}(:,1);
            for k=1:6
                HelpBoundaryT(:,app.NoOfChannels*6+k)=app.PortTapers{app.NoOfChannels+1}(:,(6+1)-k);
            end


            HelpBoundary=[app.PortTrenches{app.NoOfChannels}(:,6) HelpBoundary app.PortTrenches{app.NoOfChannels+1}(:,1)];
            HelpBoundaryT=[app.PortTrenches{app.NoOfChannels}(:,6) HelpBoundaryT app.PortTrenches{app.NoOfChannels+1}(:,1)];

            app.SiBoundary=[app.SiBoundary HelpBoundary];
            
            % this was for the comsol version
            % the GDS2 version gets a wider base

%            HelpBoundary2=HelpBoundary; % most points are identical, so just copy them all and modify then...
            HelpBoundary2=HelpBoundary; % most points are identical, so just copy them all and add modified ones (which is better for boolean operations on the GDS file...)
            dir1=[sin(app.PortDirections(app.NoOfChannels)) ; -cos(app.PortDirections(app.NoOfChannels))];  % most right distinct port
            dir2=[sin(app.PortDirections(app.NoOfChannels+1)) ; -cos(app.PortDirections(app.NoOfChannels+1))];  % common port, most left
%            HelpBoundary2(:,1)=HelpBoundary2(:,1)-app.BoundaryDistance*1e-3*dir1;   % shift the most right point a little further to the right
%            HelpBoundary2(:,app.NoOfChannels*4+4)=HelpBoundary2(:,app.NoOfChannels*4+4)+app.BoundaryDistance*1e-3*dir2; % and the most left point a little further to the left
            HelpBoundary2=[HelpBoundary2(:,1)-app.BoundaryDistance*1e-3*dir1 HelpBoundary2];   % shift the most right point a little further to the right and add as first point
            HelpBoundary2(:,end+1)=HelpBoundary2(:,end)+app.BoundaryDistance*1e-3*dir2; % shift the most left point a little further to the left and add as last point
            HelpBoundaryT=[HelpBoundaryT(:,1)-app.BoundaryDistance*1e-3*dir1 HelpBoundaryT];   % shift the most right point a little further to the right and add as first point
            HelpBoundaryT(:,end+1)=HelpBoundaryT(:,end)+app.BoundaryDistance*1e-3*dir2; % shift the most left point a little further to the left and add as last point
            
            HelpBoundary3=HelpBoundary2;    % make copy of GDS-HelpBoundary2 to add points for SiO2-Layer of slab region
            %HelpBoundary3=[HelpBoundary3(:,1)-app.SiO2SlabEnhancement*1e-3*dir1  HelpBoundary3]; % shift the most right point a little further to the right and add as first point
            %HelpBoundary3(:,end+1)=HelpBoundary3(:,end)+app.SiO2SlabEnhancement*1e-3*dir2; % shift the most left point a little further to the left and add as last point
            
            UpperSiBoundaryGDS2=app.SiBoundaryGDS2;
            app.SiBoundaryGDS2=[app.SiBoundaryGDS2 HelpBoundary2];  % complete the Si-boundary for GDS2 layout
            app.SiBoundaryGDS2T=[UpperSiBoundaryGDS2 HelpBoundaryT];  % complete the Si-boundary for GDS2 layout with necessary chamfers
            app.SiBoundaryGDS2T=app.BevelEdges(app.SiBoundaryGDS2T, app.GDS2AttBevelAngle, app.GDS2Bevel);

%            plot(myAxes,HelpBoundary3(1,:),HelpBoundary3(2,:), 'k-');
%            plot(myAxes,app.SiO2BoundaryGDS2(1,:),app.SiO2BoundaryGDS2(2,:), 'r-');
            
            %app.SiO2BoundaryGDS2
            %HelpBoundary3
            app.SiO2BoundaryGDS2=[app.SiO2BoundaryGDS2 HelpBoundary3];   % complete the SiO2-boundary for the GDS2 layout

%            plot(myAxes,app.SiO2BoundaryGDS2(1,:),app.SiO2BoundaryGDS2(2,:), 'r-');
%            plot(myAxes,HelpBoundary(1,:),HelpBoundary(2,:), 'k-');
%            plot(myAxes,HelpBoundary(1,:),HelpBoundary(2,:), 'k-');

            % and add the lower parts of the attenuators...
            app.AttenuatorLeft=[app.AttenuatorLeft HelpBoundary(:,end) HelpBoundary2(:,end)];
            app.AttenuatorRight=[app.AttenuatorRight HelpBoundary2(:,1) HelpBoundary(:,1)];
            
            app.AttenuatorLeft=app.BevelEdges(app.AttenuatorLeft, app.GDS2AttBevelAngle, app.GDS2AttBevel);
            app.AttenuatorRight=app.BevelEdges(app.AttenuatorRight, app.GDS2AttBevelAngle, app.GDS2AttBevel);

            % now expand SiO2BoundaryGDS2 by Block Trench Width
            if (app.CladdingBlockModeCheckBox.Value==true)
                %app.SiO2BoundaryGDS2
                SiO2BoundaryGDS2help = app.SiO2BoundaryGDS2;
                
                for i=2:size(app.SiO2BoundaryGDS2,2)-1
                    SiO2BoundaryGDS2help(:,i) = app.NewPointOnLeft( app.SiO2BoundaryGDS2(:,i-1), app.SiO2BoundaryGDS2(:,i), app.SiO2BoundaryGDS2(:,i+1), app.BlockModeTrenchwidth);
                end
                SiO2BoundaryGDS2help(:,1) = app.NewPointOnLeft( app.SiO2BoundaryGDS2(:,end), app.SiO2BoundaryGDS2(:,1), app.SiO2BoundaryGDS2(:,2), app.BlockModeTrenchwidth);
                SiO2BoundaryGDS2help(:,end) = app.NewPointOnLeft( app.SiO2BoundaryGDS2(:,end-1), app.SiO2BoundaryGDS2(:,end), app.SiO2BoundaryGDS2(:,1), app.BlockModeTrenchwidth);
                %SiO2BoundaryGDS2help
                app.SiO2BoundaryGDS2 = SiO2BoundaryGDS2help;
            end
            
            %show to the user if required
            if app.ShowGratingAdditionsCheckBox.Value==true
                for i=1:app.NoOfChannels+1
                    plot(myAxes,[app.PortTapers{i}(1,:) app.PortTapers{i}(1,1)], [app.PortTapers{i}(2,:) app.PortTapers{i}(2,1)], 'b-');
                    plot(myAxes,[app.PortTrenches{i}(1,:) app.PortTrenches{i}(1,1)], [app.PortTrenches{i}(2,:) app.PortTrenches{i}(2,1)], 'b-');
                    plot(myAxes,[app.PortTrenchesCon(1,:) app.PortTrenchesCon(1,1)], [app.PortTrenchesCon(2,:) app.PortTrenchesCon(2,1)], 'b:', 'LineWidth',2);
                    plot(myAxes,app.PortWGPositions(i,1), app.PortWGPositions(i,2), 'rx');
                end
                plot(myAxes,[app.SiBoundary(1,:) app.SiBoundary(1,1)],[app.SiBoundary(2,:) app.SiBoundary(2,1)], 'r:');
                plot(myAxes,[app.SiBoundaryGDS2(1,:) app.SiBoundaryGDS2(1,1)],[app.SiBoundaryGDS2(2,:) app.SiBoundaryGDS2(2,1)], 'g:');
                plot(myAxes,[app.SiO2BoundaryGDS2(1,:) app.SiO2BoundaryGDS2(1,1)],[app.SiO2BoundaryGDS2(2,:) app.SiO2BoundaryGDS2(2,1)], 'k--');
            end
        end

        
        
        
        function NewMidPoint = NewPointOnLeft(app, p1, p2, p3, Dist)
            p1p2 = p2-p1;
            p2p3 = p3-p2;
            p1p2P = [-p1p2(2); p1p2(1)]/norm(p1p2);
            p2p3P = [-p2p3(2); p2p3(1)]/norm(p2p3);
            
            %as the following method becomes numerically instable, we have to choose another algorithm for points
            %(almost) on a single line
                %%LineLineIntersection(app, p1, p1dir, p3, p3dir)
                %NewMidPoint = app.LineLineIntersection( p1+p1p2P*Dist*1e-3, p1p2, p2+p2p3P*Dist*1e-3, p2p3);
            % So time for a new method:
            
            % Displacement direction is on the line from p2 with direction (p1p2P+p2p3P)
            dispdir=p1p2P+p2p3P;
            
            % Now we have to find, how much to displace. This is the intersection of the described line and a line
            % parallel to p1 -> p2 with an offset of Dist.
            
            NewMidPoint = app.LineLineIntersection( p1+p1p2P*Dist*1e-3, p1p2, p2, dispdir);
        end

        
        
        function Area=CalculateArea(app, Boundary)
            % calculates the area of a polygon, given by Boundary
            PolyArea=0;
            for i=1:size(Boundary,2)
                inext=i+1;
                if inext>size(Boundary,2)
                    inext=1;
                end
                %Method from https://www.inf.hs-flensburg.de/lang/algorithmen/geo/polygon.htm
                PolyArea=PolyArea+ (Boundary(1,inext)-Boundary(1,i))*(Boundary(2,inext)+Boundary(2,i))/2;
            end
            Area=PolyArea;
        end
        
        
        
        function results = BuildCommentString(app)
            CommStr='';
%            CommStr=sprintf('%s\n',CommStr);
%            CommStr=sprintf('Material system: %s\nNo. of channels: %u\n1. stigmatic point: %u\n2. stigmatic point:%u\n',app.EffectiveRefractiveIndex, app.NoOfChannels, app.TSPChannels(1), app.TSPChannels(2));
%            CommStr=sprintf('%s',CommStr);
%            CommStr=['Material system: ' app.EffectiveRefractiveIndex '\nNo. of channels: ' num2str(app.NoOfChannels) '\n1. stigmatic point: ' num2str(app.TSPChannels(1)) '2. stigmatic point: ' num2str(app.TSPChannels(2))];
            ChanStr=[];
            for i=1:app.NoOfChannels
                ChanStr=[ChanStr num2str(i) '  ' num2str(app.ChannelProps.X(i)) '  ' num2str(app.ChannelProps.Y(i)) '  ' num2str(app.ChannelProps.Wavelength(i)) '  ' num2str(app.ChannelProps.ERI(i),'%.15f') '\n'];
            end
            
            CommStr=['Material system: ' app.MaterialChooserDropDown.Value...
                '\nNo. of channels: ' num2str(app.NoOfChannels)...
                '\n1. stigmatic point: ' num2str(app.TSPChannels(1))...
                '\n2. stigmatic point: ' num2str(app.TSPChannels(2))...
                '\nGrating order: ' num2str(app.GratingOrder)...
                '\nMin. radius: ' num2str(app.RadiusMin)...
                ' µm\nMax. radius: ' num2str(app.RadiusMax)...
                ' µm\n2. TSP radius offset: ' num2str(app.RadiusOffset2TSP)...
                ' µm\nGrating opening angle: ' num2str(app.GratingOpeningAngle)...
                '°\nBlaze shift: ' num2str(app.BlazeShift)...
                ' µm\nCommon port: X=' num2str(app.CommonProps.X(1)) ' Y=' num2str(app.CommonProps.Y(1))...
                '\nDistinct ports:'...
                '\nNo.  X(µm)  Y(µm)  Wavelength(nm)  Eff.RefractiveIndex'...
                '\n'...
                ChanStr...
                '\n--BraggGratings--'...
                '\nPeriods: ' num2str(app.BraggPeriods)...
                '\nPeriodLength: ' num2str(app.BraggPeriodLength)...
                ' nm\nSiO2 length: ' num2str(app.BraggSiO2Length)...
                ' nm\nBoundary distance: ' num2str(app.BoundaryDistance)...
                ' nm\nReflector spacing: ' num2str(app.ReflectorSpacing)...
                ' nm\n\n--Port layout--'...
                '\nWaveguide width: ' num2str(app.TaperWGWidth)...
                ' nm\nTaper end width: ' num2str(app.TaperWidth)...
                ' nm\nTaper length: ' num2str(app.TaperLength)...
                ' nm\nWaveguide stub length: ' num2str(app.WaveguideStubLength)...
                ' nm\nTaper shift: ' num2str(app.TaperShift)...
                ' nm\nTrench width: ' num2str(app.TrenchWidth)...
                ' nm\nTrench2 width: ' num2str(app.Trench2Width)...
                ' nm\n\nComsol min. wavelength: ' num2str(app.ComsolMinWavelength)...
                ' nm\nComsol max. wavelength: ' num2str(app.ComsolMaxWavelength)...
                ' nm\nComsol wavelength step: ' num2str(app.ComsolWavelengthStep)...
                ' nm\nComsol free space area partitions: ' num2str(app.ComsolNoOfPartitions)...
                '\nComsol border (0=reduced, 1=full (GDS)): ' num2str(app.ComsolBorder)...
                '\nComsol S-bends (0=no, 1=yes): ' num2str(app.ComsolSBends)...
                '\nGDS Echelle ID: ' sprintf('%s',app.EchelleID)...
                '\n\nGDS (and Comsol) waveguide length for S-bends: ' num2str(app.GDSWaveguideLength)...
                ' nm\nGDS (and Comsol) waveguide pitch: ' num2str(app.GDSWaveguidePitch)...
                ' nm\nGDS (and Comsol) curve discretization: ' num2str(app.GDSWGCurveDiscretization)...
                '\n\nGDS Beveling min. Angle: ' num2str(app.GDS2AttBevelAngle)...
                '°\nGDS Beveling Bevel size: ' num2str(app.GDS2Bevel)...
                'nm\nGDS Beveling Bevel size for attenuators: ' num2str(app.GDS2AttBevel)...
                'nm\nCladding Block Mode: ' num2str(app.CladdingBlockModeCheckBox.Value)...
                '\nBlock Mode Trench Width: ' num2str(app.BlockModeTrenchwidth)...
                'nm\n\nGrating points found: ' num2str(app.GratingpointsfoundField.Value)...
                '\nReal min. radius: ' num2str(app.RealminradiusmField.Value)...
                ' µm\nReal max. radius: ' num2str(app.RealmaxradiusmField.Value)...
                ' µm\nReal opening angle: ' num2str(app.RealopeningangleField.Value)...
                '°\nArea of reduced free space area: ' num2str(app.AreaSiBoundary)...
                ' µm²\nArea of full free space area: ' num2str(app.AreaSiBoundaryGDS2)...
                ' µm²\nArea of full free space area cladding: ' num2str(app.AreaSiO2BoundaryGDS2)...
                ' µm²\n'...
                '\nProject generated: ' datestr(now,31)...
                '\nwith: ' app.ProgramVersionString1...
                '\n      ' app.ProgramVersionString2...
                '\n      ' app.ProgramVersionString3...
                '\n'...
                %'\n ' num2str(app.)...
                ];
            
            %CommStr
            results=CommStr;
            
        end
    end

    
    
    
    methods (Access = public)
        
        function UpdateDialogParameters(app,myParameters)
            % Store inputs as properties
            app.DialogBoxParameters=myParameters;
            % do any updates, which might be required for changed dialog
            % box parameters...
            
            % Re-enable the OpenDialogBoxTest-button
            app.OpenDialogBoxTestButton.Enable = 'on';
        end
        
        function UpdateChannelWavelengths(app, ChannelWavelengths)
            %put new wavelengths into ChannelProperties table and update
            %uitable here
            
            app.ChannelProps.Wavelength(1:app.NoOfChannels)=ChannelWavelengths(:);
            app.ChannelProps=UpdateChannelERI(app,app.ChannelProps);
            %app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(:,:));
            app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(1:app.NoOfChannels,:));
            
            % Re-enable the ITUChannelChooser-button
            app.ITUChannelChooserButton.Enable = 'on';
        end
        
        function DontUpdateChannelWavelengths(app)
            % Re-enable the ITUChannelChooser-button
            app.ITUChannelChooserButton.Enable = 'on';
        end
        
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %set nice number format for tables
            format longG
            app.ChannelPropsUITable.ColumnFormat={'shortG','longG','longG','longG','longG'};
            app.CommonPropsUITable.ColumnFormat={'longG','longG'};
            app.PortSearchBoundingBox.ColumnFormat={'longG','longG'};
            
            app.ProgramTitleLabel.Text=app.ProgramVersionString1;
            app.ProgramSubTitle1Label.Text=app.ProgramVersionString2;
            app.ProgramSubTitle2Label.Text=app.ProgramVersionString3;
                                    
            % --- Effective Refractive Index chooser ---
            app.EffectiveRefractiveIndices=readtable(app.EffectiveRefractiveIndicesFileName);
            ERI_VNames=app.EffectiveRefractiveIndices.Properties.VariableNames;
            if (strcmp(ERI_VNames{1},'Wavelength_nm')==0)
                %error('Wrong Table for effective refractive inidices: first column has wrong name (not "Wavelength_nm")');
                uialert(app.UIFigure,'Wrong Table for effective refractive inidices: first column has wrong name (not "Wavelength_nm")','ERROR!',"Icon","error","Modal",true);
            end
            ERI_Thicknesses=ERI_VNames(2:end);
            app.MaterialChooserDropDown.Items=ERI_Thicknesses;
            app.MaterialChooserDropDown.Value=app.EffectiveRefractiveIndexDefault;
            app.EffectiveRefractiveIndex=UpdateERIMatrix(app,app.MaterialChooserDropDown.Value);
            % --- End of Effective Refractive Index chooser ---
            
            X=0;    % X-position of common port in µm
            Y=0;    % Y-position of common port in µm
            app.CommonProps=table(X,Y);
            app.NoOfChannels=app.NumberofChannelsSpinner.Value;
%            app.ChannelProps=table([1:app.NumberofChannelsSpinner.Value]',zeros(app.NumberofChannelsSpinner.Value,1),zeros(app.NumberofChannelsSpinner.Value,1),zeros(app.NumberofChannelsSpinner.Value,1));
%            app.ChannelProps.Properties.VariableNames={'ChannelNo','X','Y','Wavelength'};
%            app.ChannelProps.Properties.VariableDescriptions={'Channel No.','X (µm)','Y (µm)','Wavelength (nm)'};
            app.ChannelProps=CreateChannelPropsTable(app,app.NumberofChannelsSpinner.Value);
            app.ChannelProps.X(1)=20;
            app.ChannelProps.X(2)=25;
            app.ChannelProps.X(3)=30;
            app.ChannelProps.Wavelength(1)=1540;
            app.ChannelProps.Wavelength(2)=1550;
            app.ChannelProps.Wavelength(3)=1560;
            app.ChannelProps=UpdateChannelERI(app,app.ChannelProps);
            app.CommonPropsUITable.Data=table2cell(app.CommonProps);
            app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(:,:));
            app.TSPChannels=[app.TSPChannel1.Value ; app.TSPChannel2.Value];
            UpdateTableColoring(app);
            
            app.GratingOrder=9;
            app.GratingorderEditField.Value=app.GratingOrder;
            app.RadiusMin=100;
            app.MinradiusmEditField.Value=app.RadiusMin;
            app.RadiusMax=1000;
            app.MaxradiusmEditField.Value=app.RadiusMax;
            app.RadiusOffset2TSP=0;
            app.TSPradiusoffsetmEditField.Value=app.RadiusOffset2TSP;
            app.GratingOpeningAngle=30;
            app.GratingopeningangleEditField.Value=app.GratingOpeningAngle;
            app.BlazeShift=0;
            app.BlazeshiftmEditField.Value=app.BlazeShift;
            
            app.RealminradiusmField.Value=0;
            app.RealmaxradiusmField.Value=0;
            app.RealopeningangleField.Value=0;
            
            X=[0; 50];
            Y=[-15; 15];
            app.PortSearchBoundingBox=table(X,Y);
            app.PortSearchBoundingBoxUITable.Data=table2cell(app.PortSearchBoundingBox);
            
            app.BraggPeriods=4;
            app.PeriodsEditField.Value=app.BraggPeriods;
            app.BraggPeriodLength=361.0;
            app.PeriodlengthEditField.Value=app.BraggPeriodLength;
            app.BraggSiO2Length=180.5;
            app.SiO2lengthEditField.Value=app.BraggSiO2Length;
            app.BoundaryDistance=10000;
            app.BoundarydistanceEditField.Value=app.BoundaryDistance;
            app.ReflectorSpacing=185;
            app.ReflectorSpacingEditField.Value=app.ReflectorSpacing;

            app.TaperLength=30000;
            app.TaperlengthEditField.Value=app.TaperLength;
            app.TaperWidth=2000;
            app.TaperendwidthEditField.Value=app.TaperWidth;
            app.TaperWGWidth=450;
            app.WaveguidewidthEditField.Value=app.TaperWGWidth;
            app.WaveguideStubLength=2000;
            app.WaveguidestublengthEditField.Value=app.WaveguideStubLength;
            app.TaperShift=0.0;
            app.TapershiftEditField.Value=app.TaperShift;
            app.TrenchWidth=1000;
            app.TrenchwidthEditField.Value=app.TrenchWidth;
            app.Trench2Width=1000;
            app.Trench2widthEditField.Value=app.Trench2Width;
            
            app.ComsolMinWavelength=1500;
            app.ComsolMinWavelengthEditField.Value=app.ComsolMinWavelength;
            app.ComsolMaxWavelength=1600;
            app.ComsolMaxWavelengthEditField.Value=app.ComsolMaxWavelength;
            app.ComsolWavelengthStep=2;
            app.ComsolWavelengthStepEditField.Value=app.ComsolWavelengthStep;
            app.ComsolNoOfPartitions=48;
            app.NoOfPartitionsEditField.Value=app.ComsolNoOfPartitions;
            app.ComsolBorder=0;
            app.ComsolBorderButton.Value=app.ComsolBorder;
            app.ComsolBorderButton.Text={'Red.','Border'};
            app.ComsolSBends=0;
            app.ComsolSBendsButton.Value=app.ComsolSBends;
            

            
            app.EchelleID = 'EG_01';
            app.EchelleNameIDEditField.Value=app.EchelleID;
            
            app.GDSWaveguideLength = 200000;
            app.WaveguideLengthEditField.Value=app.GDSWaveguideLength;
            
            app.GDSWaveguidePitch = 2500;
            app.WaveguidePitchEditField.Value=app.GDSWaveguidePitch;
            
            app.GDSWGCurveDiscretization = 200;
            app.WGCurveDiscretizationEditField.Value=app.GDSWGCurveDiscretization;
            
            app.GDS2DBUnit = 1e-9;
            app.GDS2DBUnitEditField.Value=app.GDS2DBUnit;
            
            app.GDS2UserUnit = 1e-6;
            app.GDS2UserUnitEditField.Value=app.GDS2UserUnit;
            
            app.GDS2AttBevelAngle=85;
            app.GDS2AttBevelAngleEditField.Value=app.GDS2AttBevelAngle;
            app.GDS2Bevel=151;
            app.GDS2BevelEditField.Value=app.GDS2Bevel;
            app.GDS2AttBevel=310;
            app.GDS2AttBevelEditField.Value=app.GDS2AttBevel;
            %app.SiO2SlabEnhancement=0;
            %app.SiO2SlabEnhancementEditField.Value=app.SiO2SlabEnhancement;

            app.CladdingBlockModeCheckBox.Value=true;
            app.BlockModeTrenchwidth=1600;
            app.BlockModeTrenchwidthEditField.Value=app.BlockModeTrenchwidth;

            
%            app.fastFigurePosition=[1020, 60, 1040, 900];
%            app.fastFigurePosition=[1020, 60, 940, 900];
            app.fastFigurePosition=[1075, 60, 940, 900];
            app.fastFigure=figure('Name','Grating Ellipses');
            app.fastFigureNumber=app.fastFigure.Number;
            if app.PushFigureWindowToSpecificPosition
                app.fastFigure.Position=app.fastFigurePosition;
            end
            app.fastAxes=axes;
            
            
            app.DialogBoxParameters=[1 2 3];    % dummy parameters to test dialog box window
        end

        % Value changed function: NumberofChannelsSpinner
        function NumberofChannelsSpinnerValueChanged(app, event)
            value = app.NumberofChannelsSpinner.Value;
            Dummy=CreateChannelPropsTable(app,value);   % create new table
            %Dummy=table([1:value]',zeros(value,1),zeros(value,1),zeros(value,1));
            %Dummy.Properties.VariableNames={'Channel No.','X (µm)','Y (µm)','Wavelength (nm)'};
            app.NoOfChannels=value;
            if value>size(app.ChannelProps,1)   % new no. of channels larger than old largest number of channels --> copy all old data to new table
                Dummy(1:size(app.ChannelProps,1),:)=app.ChannelProps(1:size(app.ChannelProps,1),:);
                app.ChannelProps=Dummy; %make the new table the new ChannelProps-table
            end
            app.ChannelProps=UpdateChannelERI(app,app.ChannelProps);
            app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(1:value,:));
            if app.TSPChannels(1)>value
                app.TSPChannels(1)=value;
                app.TSPChannel1.Value=app.TSPChannels(1);
            end
            if app.TSPChannels(2)>value
                app.TSPChannels(2)=value;
                app.TSPChannel2.Value=app.TSPChannels(2);
            end
            UpdateTableColoring(app);
            
            %app.ChannelPropsUITable.Data=app.ChannelProps{:,:};
            %ChannelProps
        end

        % Cell edit callback: ChannelPropsUITable
        function ChannelPropsUITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            %msgbox(num2str( indices(1) ));
            app.ChannelProps{indices(1),indices(2)}=newData;
            app.ChannelProps=UpdateChannelERI(app,app.ChannelProps);
            app.ChannelPropsUITable.Data=table2cell(app.ChannelProps(:,:));
        end

        % Value changed function: TSPChannel1
        function TSPChannel1ValueChanged(app, event)
            value = app.TSPChannel1.Value;
            if value<=app.NoOfChannels
                app.TSPChannels(1)=value;
                UpdateTableColoring(app);
            else
                app.TSPChannel1.Value=app.NoOfChannels;
            end

        end

        % Value changed function: TSPChannel2
        function TSPChannel2ValueChanged(app, event)
            value = app.TSPChannel2.Value;
            if value<=app.NoOfChannels
                app.TSPChannels(2)=value;                
                UpdateTableColoring(app);
            else
                app.TSPChannel2.Value=app.NoOfChannels;
            end
            
        end

        % Cell edit callback: CommonPropsUITable
        function CommonPropsUITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            app.CommonProps{indices(1),indices(2)}=newData;
        end

        % Callback function
        function OpenDialogBoxTestButtonPushed(app, event)
            % disable button to prevent users pressing it multiple times
            % opening multiple windows
            app.OpenDialogBoxTestButton.Enable = 'off';
            
            % Open the options dialog and pass inputs
            app.DialogBoxApp = MarcEchelle_DialogBox(app, app.DialogBoxParameters);

        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            % clean up...
            
            %close all apps
            if isgraphics(app.fastFigure)~=0
                close(app.fastFigure);
            end
            delete(app.DialogBoxApp);
            delete(app.ITUChannelChooserApp);
            delete(app);
        end

        % Value changed function: MaterialChooserDropDown
        function MaterialChooserDropDownValueChanged(app, event)
            value = app.MaterialChooserDropDown.Value;
            %msgbox(num2str( indices(1) ));
            %msgbox(value);
            app.EffectiveRefractiveIndex=UpdateERIMatrix(app,app.MaterialChooserDropDown.Value);
            %msgbox(num2str( app.EffectiveRefractiveIndex(1,2) ));
            app.ChannelProps=UpdateChannelERI(app,app.ChannelProps);
            app.ChannelPropsUITable.Data=app.ChannelProps{:,:};
        end

        % Button pushed function: ITUChannelChooserButton
        function ITUChannelChooserButtonPushed(app, event)
            % disable button to prevent users pressing it multiple times
            % opening multiple windows
            app.ITUChannelChooserButton.Enable = 'off';
            
            app.ITUChannelChooserParameters=[app.NoOfChannels];
            
            % Open the options dialog and pass inputs
            app.ITUChannelChooserApp = MarcEchelle_ITUChooser(app, app.ITUChannelChooserParameters);

            
        end

        % Value changed function: GratingorderEditField
        function GratingorderEditFieldValueChanged(app, event)
            value = app.GratingorderEditField.Value;
            value=round(value);
            app.GratingOrder=value;
            app.GratingorderEditField.Value=app.GratingOrder;
        end

        % Value changed function: MinradiusmEditField
        function MinradiusmEditFieldValueChanged(app, event)
            value = app.MinradiusmEditField.Value;
            app.RadiusMin=value;
        end

        % Value changed function: MaxradiusmEditField
        function MaxradiusmEditFieldValueChanged(app, event)
            value = app.MaxradiusmEditField.Value;
            app.RadiusMax=value;
        end

        % Value changed function: TSPradiusoffsetmEditField
        function TSPradiusoffsetmEditFieldValueChanged(app, event)
            value = app.TSPradiusoffsetmEditField.Value;
            app.RadiusOffset2TSP=value;
        end

        % Value changed function: GratingopeningangleEditField
        function GratingopeningangleEditFieldValueChanged(app, event)
            value = app.GratingopeningangleEditField.Value;
            app.GratingOpeningAngle=value;
        end

        % Value changed function: BlazeshiftmEditField
        function BlazeshiftmEditFieldValueChanged(app, event)
            value = app.BlazeshiftmEditField.Value;
            app.BlazeShift=value;
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            % draw ellipses to seperate window, because it's waaaaay faster
            % than using a figure in the App Designer window
%            figure(app.fastFigureNumber);
            if isgraphics(app.fastAxes)==0
%                app.fastAxes=axes;
                app.fastFigure=figure('Name','Grating Ellipses');
                app.fastFigureNumber=app.fastFigure.Number;
                if app.PushFigureWindowToSpecificPosition
                    app.fastFigure.Position=app.fastFigurePosition;
                end
                app.fastAxes=axes;
            end
            

            
            myAxes=app.fastAxes;
            axes(myAxes);
            axis('equal');
%            cla(app.UIAxes);    % clear axes: clears the diagram field
%            hold(app.UIAxes,'on');  % keep all successive draws to the diagram field
            cla(myAxes);    % clear axes: clears the diagram field
            hold(myAxes,'on');  % keep all successive draws to the diagram field
            drawnow limitrate
            
            
            d = uiprogressdlg(app.UIFigure,'Title','Calculating','Message','Calculating ellipse intersections');
            
            app.CalculateEllipseIntersections(myAxes);

%{            
            firstpoint=true;
            
            app.EllipseIntersections=[];
            firstAngle=NaN;
            MaxAngleReached=false;
            i=0;
            while (MaxAngleReached==false)
%            for i=0:40
                Radius1=app.RadiusMin+i*(app.ChannelProps.Wavelength(app.TSPChannels(1))/app.ChannelProps.ERI(app.TSPChannels(1))*app.GratingOrder*1e-3);
                Radius2=app.RadiusMin+app.RadiusOffset2TSP+i*(app.ChannelProps.Wavelength(app.TSPChannels(2))/app.ChannelProps.ERI(app.TSPChannels(2))*app.GratingOrder*1e-3);
%                DrawEllipse(app.UIAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.RadiusMin+i*(app.ChannelProps.Wavelength(app.TSPChannels(1))/app.ChannelProps.ERI(app.TSPChannels(1))*app.GratingOrder*1e-3),120);
%                DrawEllipse(app.UIAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), app.RadiusMin+app.RadiusOffset2TSP+i*(app.ChannelProps.Wavelength(app.TSPChannels(2))/app.ChannelProps.ERI(app.TSPChannels(2))*app.GratingOrder*1e-3),120);
                DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), Radius1,120, [ 0.75, 0.75, 0.75]);
                DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius2,120, [ 0.75, 0.75, 1.0]);
                if ~firstpoint
                    % all subsequent intersetions must be within a 20°-wedge of the previous intersection
                    app.EllipseIntersections=[app.EllipseIntersections; EllipseIntersectionFinder2(myAxes, app.CommonProps.X(1),app.CommonProps.Y(1), app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius1, Radius2, app.EllipseIntersections(end, 4), 20, 'red')];
                    
                else
                    % First intersection must be in upper half and is the leftmost intersection, if there are more than one
                    app.EllipseIntersections=EllipseIntersectionFinder2(myAxes, app.CommonProps.X(1),app.CommonProps.Y(1), app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), Radius1, Radius2, 90, 180, 'red');
                    firstAngle=app.EllipseIntersections(1,4);
                    if (isnan(firstAngle)==false)
                        firstpoint=false;
                    end
                end
                if (firstAngle-app.EllipseIntersections(end,4)>=app.GratingOpeningAngle) || (Radius1>=app.RadiusMax) || (Radius2>=app.RadiusMax)
                    MaxAngleReached=true;
                end
                    
                i=i+1;
            end
            %app.EllipseIntersections
            %Strip NaN-values

            %if the intersection list is not empty, look for and strip NaN values
            if isempty(app.EllipseIntersections)==false
                lastrow=0;
                i=1;
                while (i<=size(app.EllipseIntersections,1)) && (isnan(app.EllipseIntersections(i,4))==false)
                    lastrow=i;
                    i=i+1;
                end
                if lastrow==0   %List doesn't start with valid number; shouldn't occure, but...
                    app.EllipseIntersections=[];
                else
                    app.EllipseIntersections=app.EllipseIntersections(1:lastrow,:); %cut out all intersection until the first NaN
                end
            end
            app.EllipseIntersections;
%}            
            app.GratingpointsfoundField.Value=size(app.EllipseIntersections,1);
            if isempty(app.EllipseIntersections)==false
                app.RealminradiusmField.Value=app.EllipseIntersections(1,3);
                app.RealmaxradiusmField.Value=app.EllipseIntersections(end,3);
                app.RealopeningangleField.Value=abs(app.EllipseIntersections(end,4)-app.EllipseIntersections(1,4));
            else
                app.RealminradiusmField.Value=0;
                app.RealmaxradiusmField.Value=0;
                app.RealopeningangleField.Value=0;
            end

                d.Value = .3; 
                d.Message = 'Calculating channel positions';
            
            if ((matlab.addons.isAddonEnabled('Parallel Computing Toolbox')==1) && (app.ParallelComputingCheckBox.Value==true))
                app.CalculateChannelPositionsPar(myAxes);    
            else
                app.CalculateChannelPositions(myAxes);    
            end
            %app.CalculateChannelPositions(myAxes);    
                d.Value = .4; 
                d.Message = 'Calculating ports center of gravity';
            app.CalculatePortsCenterOfGravity(myAxes);
                d.Value = .5; 
                d.Message = 'Calculating grating center';
            app.CalculateGratingCenter(myAxes); %must have a valid PortsCenterOfGravity as input!!!
                d.Value = .6; 
                d.Message = 'Calculating grating reflector directions';
            app.CalculateGratingReflectorDirections(myAxes);
                d.Value = .7; 
                d.Message = 'Calculating port directions';
            app.CalculatePortDirections(myAxes);
                d.Value = .8; 
                d.Message = 'Calculating grating reflectors';
            app.CalculateGratingReflectors(myAxes);
                d.Value = .9; 
                d.Message = 'Calculating port tapers';
            app.CalculatePortTapers(myAxes);
            
            app.AreaSiBoundary=app.CalculateArea(app.SiBoundary);
            app.AreaSiBoundaryGDS2=app.CalculateArea(app.SiBoundaryGDS2);
            app.AreaSiO2BoundaryGDS2=app.CalculateArea(app.SiO2BoundaryGDS2);
            app.AreaReducedField.Value=app.AreaSiBoundary;
            app.AreaFullField.Value=app.AreaSiBoundaryGDS2;
            app.AreaSiO2Field.Value=app.AreaSiO2BoundaryGDS2;
            
                close(d);
            
%            DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.RadiusMin,120,'black');
%            DrawEllipse(myAxes,app.CommonProps.X(1),app.CommonProps.Y(1),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), app.RadiusMin+app.RadiusOffset2TSP,120,'blue');

%            DrawEllipse(myAxes,app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), app.RadiusMin,120,'red');
            
%            DrawEllipse(myAxes, app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)),app.ChannelProps.X(app.TSPChannels(1)),app.ChannelProps.Y(app.TSPChannels(1)), app.RadiusMin*2, 120  )
%            DrawEllipse(myAxes, app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)),app.ChannelProps.X(app.TSPChannels(2)),app.ChannelProps.Y(app.TSPChannels(2)), (app.RadiusMin+app.RadiusOffset2TSP)*2, 120  )

%            hold(app.UIAxes,'off');  % an upcoming draw to the diagram field would clear all previous drawings now
            hold(myAxes,'off');  % an upcoming draw to the diagram field would clear all previous drawings now

        end

        % Cell edit callback: PortSearchBoundingBoxUITable
        function PortSearchBoundingBoxUITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            app.PortSearchBoundingBox{indices(1),indices(2)}=newData;
        end

        % Value changed function: PeriodsEditField
        function PeriodsEditFieldValueChanged(app, event)
            value = app.PeriodsEditField.Value;
            app.BraggPeriods=value;
        end

        % Value changed function: PeriodlengthEditField
        function PeriodlengthEditFieldValueChanged(app, event)
            value = app.PeriodlengthEditField.Value;
            app.BraggPeriodLength=value;
        end

        % Value changed function: SiO2lengthEditField
        function SiO2lengthEditFieldValueChanged(app, event)
            value = app.SiO2lengthEditField.Value;
            app.BraggSiO2Length=value;
        end

        % Value changed function: BoundarydistanceEditField
        function BoundarydistanceEditFieldValueChanged(app, event)
            value = app.BoundarydistanceEditField.Value;
            app.BoundaryDistance=value;
        end

        % Value changed function: ReflectorSpacingEditField
        function ReflectorSpacingEditFieldValueChanged(app, event)
            value = app.ReflectorSpacingEditField.Value;
            app.ReflectorSpacing=value;
        end

        % Value changed function: WaveguidewidthEditField
        function WaveguidewidthEditFieldValueChanged(app, event)
            value = app.WaveguidewidthEditField.Value;
            app.TaperWGWidth=value;            
        end

        % Value changed function: TaperendwidthEditField
        function TaperendwidthEditFieldValueChanged(app, event)
            value = app.TaperendwidthEditField.Value;
            app.TaperWidth=value;         
        end

        % Value changed function: TaperlengthEditField
        function TaperlengthEditFieldValueChanged(app, event)
            value = app.TaperlengthEditField.Value;
            app.TaperLength=value;
        end

        % Value changed function: WaveguidestublengthEditField
        function WaveguidestublengthEditFieldValueChanged(app, event)
            value = app.WaveguidestublengthEditField.Value;
            app.WaveguideStubLength=value;
        end

        % Value changed function: TapershiftEditField
        function TapershiftEditFieldValueChanged(app, event)
            value = app.TapershiftEditField.Value;
            app.TaperShift=value;
        end

        % Value changed function: TrenchwidthEditField
        function TrenchwidthEditFieldValueChanged(app, event)
            value = app.TrenchwidthEditField.Value;
            app.TrenchWidth=value;
        end

        % Value changed function: Trench2widthEditField
        function Trench2widthEditFieldValueChanged(app, event)
            value = app.Trench2widthEditField.Value;
            app.Trench2Width=value;
        end

        % Value changed function: ComsolMinWavelengthEditField
        function ComsolMinWavelengthEditFieldValueChanged(app, event)
            value = app.ComsolMinWavelengthEditField.Value;
            app.ComsolMinWavelength=value;
        end

        % Value changed function: ComsolMaxWavelengthEditField
        function ComsolMaxWavelengthEditFieldValueChanged(app, event)
            value = app.ComsolMaxWavelengthEditField.Value;
            app.ComsolMaxWavelength=value;
        end

        % Value changed function: ComsolWavelengthStepEditField
        function ComsolWavelengthStepEditFieldValueChanged(app, event)
            value = app.ComsolWavelengthStepEditField.Value;
            app.ComsolWavelengthStep=value;
        end

        % Value changed function: ComsolBorderButton
        function ComsolBorderButtonValueChanged(app, event)
            value = app.ComsolBorderButton.Value;
            app.ComsolBorder=value;
            if value==0
                app.ComsolBorderButton.Text={'Red.','Border'};
            else
                app.ComsolBorderButton.Text={'Full','Border'};
            end
        end

        % Value changed function: ComsolSBendsButton
        function ComsolSBendsButtonValueChanged(app, event)
            value = app.ComsolSBendsButton.Value;
            app.ComsolSBends=value;
        end

        % Value changed function: NoOfPartitionsEditField
        function NoOfPartitionsEditFieldValueChanged(app, event)
            value = app.NoOfPartitionsEditField.Value;
            app.ComsolNoOfPartitions=value;
        end

        % Value changed function: EchelleNameIDEditField
        function EchelleNameIDEditFieldValueChanged(app, event)
            value = app.EchelleNameIDEditField.Value;
            app.EchelleID = value;
        end

        % Value changed function: WaveguideLengthEditField
        function WaveguideLengthEditFieldValueChanged(app, event)
            value = app.WaveguideLengthEditField.Value;
            app.GDSWaveguideLength = value;
        end

        % Value changed function: WaveguidePitchEditField
        function WaveguidePitchEditFieldValueChanged(app, event)
            value = app.WaveguidePitchEditField.Value;
            app.GDSWaveguidePitch = value;
        end

        % Value changed function: WGCurveDiscretizationEditField
        function WGCurveDiscretizationEditFieldValueChanged(app, event)
            value = app.WGCurveDiscretizationEditField.Value;
            app.GDSWGCurveDiscretization = value;
        end

        % Button pushed function: COMSOLModelButton
        function COMSOLModelButtonPushed(app, event)
            filetime = datestr(now,30);
            filename=sprintf("Echelle_%s.mph",filetime);
            %filename=sprintf("Echelle.mph");
            [filename, pathname] = uiputfile({'*.mph','COMSOL project';'*.*','All files (*.*)'},'Filename for COMSOL Model File',filename);
            if filename~=0
                % Now generate a COMSOL simulation project
                fprintf('\n\nStarting Comsol simulation transfer...\n');
                %check if comsolserver is running
                [~,result] = system('tasklist /FI "imagename eq comsolmphserver.exe" /fo table /nh');
                
                if strcmp(result(2:16), 'comsolmphserver')
                    %do nothing, just lets skript run further
                    fprintf('Comsolserver seems to be running\nBuilding COMSOL project %s%s ...\n',pathname,filename);
                    CommentString=app.BuildCommentString;
                    if app.ComsolBorder==0
                        BuildComsolProject(app.UIFigure,pathname,filename,app.BraggReflectors, app.PortTapers, app.PortTrenches, app.PortTrenchesCon, app.SiBoundary, app.EffectiveRefractiveIndex, app.ComsolMinWavelength, app.ComsolMaxWavelength, app.ComsolWavelengthStep,[app.PortWGPositions(end,1), app.PortWGPositions(end,2), (app.PortDirections(end)-3*pi/2)], CommentString, app.ComsolNoOfPartitions, app.ComsolSBends, app.TaperWGWidth, app.TrenchWidth, app.Trench2Width, app.GDSWaveguideLength, app.GDSWaveguidePitch, app.PortWGPositions, app.PortDirections, app.GDSWGCurveDiscretization, app.CladdingBlockModeCheckBox.Value, app.BlockModeTrenchwidth);
                    else
                        BuildComsolProject(app.UIFigure,pathname,filename,app.BraggReflectors, app.PortTapers, app.PortTrenches, app.PortTrenchesCon, app.SiBoundaryGDS2, app.EffectiveRefractiveIndex, app.ComsolMinWavelength, app.ComsolMaxWavelength, app.ComsolWavelengthStep,[app.PortWGPositions(end,1), app.PortWGPositions(end,2), (app.PortDirections(end)-3*pi/2)], CommentString, app.ComsolNoOfPartitions, app.ComsolSBends, app.TaperWGWidth, app.TrenchWidth, app.Trench2Width, app.GDSWaveguideLength, app.GDSWaveguidePitch, app.PortWGPositions, app.PortDirections, app.GDSWGCurveDiscretization, app.CladdingBlockModeCheckBox.Value, app.BlockModeTrenchwidth);
                    end
%                    test=[filename '.params']
                else
                    fprintf('\n\nComsolserver does not seem to be running. Start "COMSOL Multiphysics X.x with MATLAB.exe" to use this Funktion.\n\n');
                end
                
            end
        end

        % Value changed function: GDS2DBUnitEditField
        function GDS2DBUnitEditFieldValueChanged(app, event)
            value = app.GDS2DBUnitEditField.Value;
            app.GDS2DBUnit = value;
        end

        % Value changed function: GDS2UserUnitEditField
        function GDS2UserUnitEditFieldValueChanged(app, event)
            value = app.GDS2UserUnitEditField.Value;
            app.GDS2UserUnit = value;
        end

        % Value changed function: GDS2AttBevelAngleEditField
        function GDS2AttBevelAngleEditFieldValueChanged(app, event)
            value = app.GDS2AttBevelAngleEditField.Value;
            app.GDS2AttBevelAngle=value;
        end

        % Value changed function: GDS2BevelEditField
        function GDS2BevelEditFieldValueChanged(app, event)
            value = app.GDS2BevelEditField.Value;
            app.GDS2Bevel=value;
        end

        % Value changed function: GDS2AttBevelEditField
        function GDS2AttBevelEditFieldValueChanged(app, event)
            value = app.GDS2AttBevelEditField.Value;
            app.GDS2AttBevel=value;
        end

        % Button pushed function: GDS2LayoutButton
        function GDS2LayoutButtonPushed(app, event)
            filetime = datestr(now,30);
            filename=sprintf("Echelle_%s.gds",filetime);
%            filename=sprintf("Echelle.gds");
            [filename, pathname] = uiputfile({'*.gds','GDSII file';'*.*','All files (*.*)'},'Filename for GDSII Layout File',filename);
            filename4comments=[filename,'_comments.txt'];
            
            if filename~=0
                fprintf('\n\nWriting parameters to _comments.txt file...\n');
                commentFileID=fopen([pathname, filename4comments],'w');
                fprintf(commentFileID,app.BuildCommentString);
                fclose(commentFileID);

                % Now generate a GDS2 layout
                fprintf('\n\nStarting GDSII output...\n');
                
%                CommentString=app.BuildCommentString;
%                BuildGDSFile(app.UIFigure,pathname,filename,app.BraggReflectors, app.PortTapers, app.PortTrenches, app.SiBoundaryGDS2, app.AttenuatorLeft, app.AttenuatorRight, [app.PortWGPositions(end,1), app.PortWGPositions(end,2), (app.PortDirections(end)-3*pi/2)], app.EchelleID, app.TaperWGWidth, app.TrenchWidth, app.GDSWaveguideLength, app.GDSWaveguidePitch, app.PortWGPositions, app.PortDirections);
                BuildGDSFile(app.UIFigure,pathname,filename,app.BraggReflectors, app.PortTapers, app.PortTrenches, app.PortTrenchesCon, app.SiBoundaryGDS2, app.SiBoundaryGDS2T, app.SiO2BoundaryGDS2, app.AttenuatorLeft, app.AttenuatorRight, [app.PortWGPositions(end,1), app.PortWGPositions(end,2), (pi/2+app.PortDirections(end))], app.EchelleID, app.TaperWGWidth, app.TrenchWidth, app.Trench2Width , app.GDSWaveguideLength, app.GDSWaveguidePitch, app.PortWGPositions, app.PortDirections, app.GDSWGCurveDiscretization, app.GDS2DBUnit, app.GDS2UserUnit, app.GDS2AttBevelAngle, app.GDS2Bevel, app.CladdingBlockModeCheckBox.Value, app.BlockModeTrenchwidth);
                
            end
            
        end

        % Value changed function: BlockModeTrenchwidthEditField
        function BlockModeTrenchwidthEditFieldValueChanged(app, event)
            value = app.BlockModeTrenchwidthEditField.Value;
            app.BlockModeTrenchwidth=value;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 950 1029];
            app.UIFigure.Name = 'UI Figure';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create NumberofChannelsLabel
            app.NumberofChannelsLabel = uilabel(app.UIFigure);
            app.NumberofChannelsLabel.HorizontalAlignment = 'right';
            app.NumberofChannelsLabel.Position = [289 856 119 22];
            app.NumberofChannelsLabel.Text = 'Number of Channels:';

            % Create NumberofChannelsSpinner
            app.NumberofChannelsSpinner = uispinner(app.UIFigure);
            app.NumberofChannelsSpinner.Limits = [2 Inf];
            app.NumberofChannelsSpinner.ValueChangedFcn = createCallbackFcn(app, @NumberofChannelsSpinnerValueChanged, true);
            app.NumberofChannelsSpinner.Tooltip = {'Number of wavelength channels'};
            app.NumberofChannelsSpinner.Position = [339 829 71 22];
            app.NumberofChannelsSpinner.Value = 3;

            % Create ChannelPropsUITable
            app.ChannelPropsUITable = uitable(app.UIFigure);
            app.ChannelPropsUITable.ColumnName = {'Channel No.'; 'X (µm)'; 'Y (µm)'; 'Wavelength (nm)'; 'Eff. ref. index'};
            app.ChannelPropsUITable.ColumnWidth = {70, 70, 70, 90, 'auto'};
            app.ChannelPropsUITable.RowName = {};
            app.ChannelPropsUITable.ColumnEditable = [false true true true false];
            app.ChannelPropsUITable.CellEditCallback = createCallbackFcn(app, @ChannelPropsUITableCellEdit, true);
            app.ChannelPropsUITable.Tooltip = {'Data of single wavelength ports'};
            app.ChannelPropsUITable.Position = [22 451 388 350];

            % Create CommonPropsUITable
            app.CommonPropsUITable = uitable(app.UIFigure);
            app.CommonPropsUITable.ColumnName = {'X (µm)'; 'Y (µm)'};
            app.CommonPropsUITable.RowName = {};
            app.CommonPropsUITable.ColumnEditable = true;
            app.CommonPropsUITable.CellEditCallback = createCallbackFcn(app, @CommonPropsUITableCellEdit, true);
            app.CommonPropsUITable.Tooltip = {'Data of common port for all wavelengths'};
            app.CommonPropsUITable.Position = [22 829 222 49];

            % Create CommonPortPositionLabel
            app.CommonPortPositionLabel = uilabel(app.UIFigure);
            app.CommonPortPositionLabel.Position = [22 877 126 22];
            app.CommonPortPositionLabel.Text = 'Common Port Position';

            % Create DistinctPortsLabel
            app.DistinctPortsLabel = uilabel(app.UIFigure);
            app.DistinctPortsLabel.Position = [22 800 76 22];
            app.DistinctPortsLabel.Text = 'Distinct Ports';

            % Create ChannelNumberSpinnerLabel
            app.ChannelNumberSpinnerLabel = uilabel(app.UIFigure);
            app.ChannelNumberSpinnerLabel.HorizontalAlignment = 'right';
            app.ChannelNumberSpinnerLabel.Position = [128 407 96 22];
            app.ChannelNumberSpinnerLabel.Text = 'Channel Number';

            % Create TSPChannel1
            app.TSPChannel1 = uispinner(app.UIFigure);
            app.TSPChannel1.Limits = [1 Inf];
            app.TSPChannel1.ValueChangedFcn = createCallbackFcn(app, @TSPChannel1ValueChanged, true);
            app.TSPChannel1.BackgroundColor = [0.8 1 0.8];
            app.TSPChannel1.Tooltip = {'Channel number of first stigmatic point of the two stigmatic points method'};
            app.TSPChannel1.Position = [243 407 55 22];
            app.TSPChannel1.Value = 1;

            % Create ChannelNumberSpinner_2Label
            app.ChannelNumberSpinner_2Label = uilabel(app.UIFigure);
            app.ChannelNumberSpinner_2Label.HorizontalAlignment = 'right';
            app.ChannelNumberSpinner_2Label.Position = [128 378 96 22];
            app.ChannelNumberSpinner_2Label.Text = 'Channel Number';

            % Create TSPChannel2
            app.TSPChannel2 = uispinner(app.UIFigure);
            app.TSPChannel2.Limits = [1 Inf];
            app.TSPChannel2.ValueChangedFcn = createCallbackFcn(app, @TSPChannel2ValueChanged, true);
            app.TSPChannel2.BackgroundColor = [0.8 1 1];
            app.TSPChannel2.Tooltip = {'Channel number of second stigmatic point of the two stigmatic points method'};
            app.TSPChannel2.Position = [243 378 55 22];
            app.TSPChannel2.Value = 3;

            % Create TSPChannel1Label
            app.TSPChannel1Label = uilabel(app.UIFigure);
            app.TSPChannel1Label.Position = [24 407 99 22];
            app.TSPChannel1Label.Text = '1. Stigmatic Point';

            % Create TSPChannel2Label
            app.TSPChannel2Label = uilabel(app.UIFigure);
            app.TSPChannel2Label.Position = [24 378 99 22];
            app.TSPChannel2Label.Text = '2. Stigmatic Point';

            % Create ProgramTitleLabel
            app.ProgramTitleLabel = uilabel(app.UIFigure);
            app.ProgramTitleLabel.FontSize = 18;
            app.ProgramTitleLabel.FontWeight = 'bold';
            app.ProgramTitleLabel.Position = [22 990 239 22];
            app.ProgramTitleLabel.Text = 'Echelle Layout Program v2';

            % Create ProgramSubTitle1Label
            app.ProgramSubTitle1Label = uilabel(app.UIFigure);
            app.ProgramSubTitle1Label.Position = [22 960 187 22];
            app.ProgramSubTitle1Label.Text = 'by Dr.-Ing. Marc Schneider';

            % Create MaterialChooserDropDownLabel
            app.MaterialChooserDropDownLabel = uilabel(app.UIFigure);
            app.MaterialChooserDropDownLabel.HorizontalAlignment = 'right';
            app.MaterialChooserDropDownLabel.Position = [16 909 97 22];
            app.MaterialChooserDropDownLabel.Text = 'Material Chooser';

            % Create MaterialChooserDropDown
            app.MaterialChooserDropDown = uidropdown(app.UIFigure);
            app.MaterialChooserDropDown.ValueChangedFcn = createCallbackFcn(app, @MaterialChooserDropDownValueChanged, true);
            app.MaterialChooserDropDown.Position = [128 909 164 22];

            % Create ITUChannelChooserButton
            app.ITUChannelChooserButton = uibutton(app.UIFigure, 'push');
            app.ITUChannelChooserButton.ButtonPushedFcn = createCallbackFcn(app, @ITUChannelChooserButtonPushed, true);
            app.ITUChannelChooserButton.BackgroundColor = [0.8392 0.9686 0.8392];
            app.ITUChannelChooserButton.Position = [311 907 153 26];
            app.ITUChannelChooserButton.Text = 'ITU Channel Chooser';

            % Create Footnote1
            app.Footnote1 = uilabel(app.UIFigure);
            app.Footnote1.FontSize = 10;
            app.Footnote1.Position = [22 434 352 22];
            app.Footnote1.Text = '(Only coordinates for stigmatic points are important, others will be calculated)';

            % Create GratingorderEditFieldLabel
            app.GratingorderEditFieldLabel = uilabel(app.UIFigure);
            app.GratingorderEditFieldLabel.HorizontalAlignment = 'right';
            app.GratingorderEditFieldLabel.Position = [545 898 76 22];
            app.GratingorderEditFieldLabel.Text = 'Grating order';

            % Create GratingorderEditField
            app.GratingorderEditField = uieditfield(app.UIFigure, 'numeric');
            app.GratingorderEditField.Limits = [1 Inf];
            app.GratingorderEditField.ValueChangedFcn = createCallbackFcn(app, @GratingorderEditFieldValueChanged, true);
            app.GratingorderEditField.Position = [636 898 100 22];
            app.GratingorderEditField.Value = 9;

            % Create MinradiusmEditFieldLabel
            app.MinradiusmEditFieldLabel = uilabel(app.UIFigure);
            app.MinradiusmEditFieldLabel.HorizontalAlignment = 'right';
            app.MinradiusmEditFieldLabel.Position = [529 869 92 22];
            app.MinradiusmEditFieldLabel.Text = 'Min. radius (µm)';

            % Create MinradiusmEditField
            app.MinradiusmEditField = uieditfield(app.UIFigure, 'numeric');
            app.MinradiusmEditField.Limits = [0 Inf];
            app.MinradiusmEditField.ValueDisplayFormat = '%12.12g';
            app.MinradiusmEditField.ValueChangedFcn = createCallbackFcn(app, @MinradiusmEditFieldValueChanged, true);
            app.MinradiusmEditField.Position = [636 869 100 22];
            app.MinradiusmEditField.Value = 100;

            % Create MaxradiusmEditFieldLabel
            app.MaxradiusmEditFieldLabel = uilabel(app.UIFigure);
            app.MaxradiusmEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxradiusmEditFieldLabel.Position = [525 840 96 22];
            app.MaxradiusmEditFieldLabel.Text = 'Max. radius (µm)';

            % Create MaxradiusmEditField
            app.MaxradiusmEditField = uieditfield(app.UIFigure, 'numeric');
            app.MaxradiusmEditField.Limits = [0 Inf];
            app.MaxradiusmEditField.ValueDisplayFormat = '%12.12g';
            app.MaxradiusmEditField.ValueChangedFcn = createCallbackFcn(app, @MaxradiusmEditFieldValueChanged, true);
            app.MaxradiusmEditField.Position = [636 840 100 22];
            app.MaxradiusmEditField.Value = 200;

            % Create TSPradiusoffsetmEditFieldLabel
            app.TSPradiusoffsetmEditFieldLabel = uilabel(app.UIFigure);
            app.TSPradiusoffsetmEditFieldLabel.HorizontalAlignment = 'right';
            app.TSPradiusoffsetmEditFieldLabel.Position = [482 811 139 22];
            app.TSPradiusoffsetmEditFieldLabel.Text = '2. TSP radius offset (µm)';

            % Create TSPradiusoffsetmEditField
            app.TSPradiusoffsetmEditField = uieditfield(app.UIFigure, 'numeric');
            app.TSPradiusoffsetmEditField.ValueDisplayFormat = '%12.12g';
            app.TSPradiusoffsetmEditField.ValueChangedFcn = createCallbackFcn(app, @TSPradiusoffsetmEditFieldValueChanged, true);
            app.TSPradiusoffsetmEditField.Position = [636 811 100 22];

            % Create GratingopeningangleLabel
            app.GratingopeningangleLabel = uilabel(app.UIFigure);
            app.GratingopeningangleLabel.HorizontalAlignment = 'right';
            app.GratingopeningangleLabel.Position = [481 782 140 22];
            app.GratingopeningangleLabel.Text = 'Grating opening angle (°)';

            % Create GratingopeningangleEditField
            app.GratingopeningangleEditField = uieditfield(app.UIFigure, 'numeric');
            app.GratingopeningangleEditField.Limits = [0 180];
            app.GratingopeningangleEditField.ValueChangedFcn = createCallbackFcn(app, @GratingopeningangleEditFieldValueChanged, true);
            app.GratingopeningangleEditField.Position = [636 782 100 22];
            app.GratingopeningangleEditField.Value = 30;

            % Create CalculateButton
            app.CalculateButton = uibutton(app.UIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.BackgroundColor = [0.8392 0.9686 0.8392];
            app.CalculateButton.Position = [482 576 268 32];
            app.CalculateButton.Text = {'Calculate'; ''};

            % Create ImageKIT
            app.ImageKIT = uiimage(app.UIFigure);
            app.ImageKIT.Position = [769 937 160 73];
            app.ImageKIT.ImageSource = 'KIT_Logo_Englisch_verysmall_Alpha.png';

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.Position = [480 948 254 62];
            app.Image.ImageSource = 'IPE-Logo_Offiziell_en_Tiny_Alpha.png';

            % Create ResultsPanel
            app.ResultsPanel = uipanel(app.UIFigure);
            app.ResultsPanel.Title = 'Results';
            app.ResultsPanel.Position = [482 377 454 185];

            % Create toeasecalculationsasseenfromcommonportpositionLabel
            app.toeasecalculationsasseenfromcommonportpositionLabel = uilabel(app.ResultsPanel);
            app.toeasecalculationsasseenfromcommonportpositionLabel.FontSize = 10;
            app.toeasecalculationsasseenfromcommonportpositionLabel.Position = [7 96 254 22];
            app.toeasecalculationsasseenfromcommonportpositionLabel.Text = 'to ease calculations as seen from common port position';

            % Create RealminradiusmEditFieldLabel
            app.RealminradiusmEditFieldLabel = uilabel(app.ResultsPanel);
            app.RealminradiusmEditFieldLabel.HorizontalAlignment = 'right';
            app.RealminradiusmEditFieldLabel.Position = [19 70 120 22];
            app.RealminradiusmEditFieldLabel.Text = 'Real min. radius (µm)';

            % Create RealminradiusmField
            app.RealminradiusmField = uieditfield(app.ResultsPanel, 'numeric');
            app.RealminradiusmField.Limits = [0 Inf];
            app.RealminradiusmField.ValueDisplayFormat = '%12.12g';
            app.RealminradiusmField.Editable = 'off';
            app.RealminradiusmField.Position = [154 70 100 22];
            app.RealminradiusmField.Value = 100;

            % Create RealmaxradiusmEditFieldLabel
            app.RealmaxradiusmEditFieldLabel = uilabel(app.ResultsPanel);
            app.RealmaxradiusmEditFieldLabel.HorizontalAlignment = 'right';
            app.RealmaxradiusmEditFieldLabel.Position = [15 41 124 22];
            app.RealmaxradiusmEditFieldLabel.Text = 'Real max. radius (µm)';

            % Create RealmaxradiusmField
            app.RealmaxradiusmField = uieditfield(app.ResultsPanel, 'numeric');
            app.RealmaxradiusmField.Limits = [0 Inf];
            app.RealmaxradiusmField.ValueDisplayFormat = '%12.12g';
            app.RealmaxradiusmField.Editable = 'off';
            app.RealmaxradiusmField.Position = [154 41 100 22];
            app.RealmaxradiusmField.Value = 2000;

            % Create RealopeningangleEditFieldLabel
            app.RealopeningangleEditFieldLabel = uilabel(app.ResultsPanel);
            app.RealopeningangleEditFieldLabel.HorizontalAlignment = 'right';
            app.RealopeningangleEditFieldLabel.Position = [14 12 125 22];
            app.RealopeningangleEditFieldLabel.Text = 'Real opening angle (°)';

            % Create RealopeningangleField
            app.RealopeningangleField = uieditfield(app.ResultsPanel, 'numeric');
            app.RealopeningangleField.Limits = [0 180];
            app.RealopeningangleField.ValueDisplayFormat = '%12.12g';
            app.RealopeningangleField.Editable = 'off';
            app.RealopeningangleField.Position = [154 12 100 22];
            app.RealopeningangleField.Value = 30;

            % Create GratingpointsfoundEditFieldLabel
            app.GratingpointsfoundEditFieldLabel = uilabel(app.ResultsPanel);
            app.GratingpointsfoundEditFieldLabel.HorizontalAlignment = 'right';
            app.GratingpointsfoundEditFieldLabel.Position = [25 127 114 22];
            app.GratingpointsfoundEditFieldLabel.Text = 'Grating points found';

            % Create GratingpointsfoundField
            app.GratingpointsfoundField = uieditfield(app.ResultsPanel, 'numeric');
            app.GratingpointsfoundField.Limits = [0 Inf];
            app.GratingpointsfoundField.Editable = 'off';
            app.GratingpointsfoundField.Position = [154 127 100 22];

            % Create reducedLabel
            app.reducedLabel = uilabel(app.ResultsPanel);
            app.reducedLabel.HorizontalAlignment = 'right';
            app.reducedLabel.Position = [278 99 49 22];
            app.reducedLabel.Text = 'reduced';

            % Create AreaReducedField
            app.AreaReducedField = uieditfield(app.ResultsPanel, 'numeric');
            app.AreaReducedField.Limits = [0 Inf];
            app.AreaReducedField.ValueDisplayFormat = '%12.12g';
            app.AreaReducedField.Editable = 'off';
            app.AreaReducedField.Position = [342 99 100 22];

            % Create fullLabel
            app.fullLabel = uilabel(app.ResultsPanel);
            app.fullLabel.HorizontalAlignment = 'right';
            app.fullLabel.Position = [302 70 25 22];
            app.fullLabel.Text = 'full';

            % Create AreaFullField
            app.AreaFullField = uieditfield(app.ResultsPanel, 'numeric');
            app.AreaFullField.Limits = [0 Inf];
            app.AreaFullField.ValueDisplayFormat = '%12.12g';
            app.AreaFullField.Editable = 'off';
            app.AreaFullField.Position = [342 70 100 22];

            % Create AreaoffreespaceareamLabel
            app.AreaoffreespaceareamLabel = uilabel(app.ResultsPanel);
            app.AreaoffreespaceareamLabel.Position = [278 128 166 22];
            app.AreaoffreespaceareamLabel.Text = 'Area of free space area (µm²):';

            % Create SiO2Label
            app.SiO2Label = uilabel(app.ResultsPanel);
            app.SiO2Label.HorizontalAlignment = 'right';
            app.SiO2Label.Position = [295 41 32 22];
            app.SiO2Label.Text = 'SiO2';

            % Create AreaSiO2Field
            app.AreaSiO2Field = uieditfield(app.ResultsPanel, 'numeric');
            app.AreaSiO2Field.ValueDisplayFormat = '%12.12g';
            app.AreaSiO2Field.Editable = 'off';
            app.AreaSiO2Field.Position = [342 41 100 22];

            % Create PortSearchBoundingBoxUITable
            app.PortSearchBoundingBoxUITable = uitable(app.UIFigure);
            app.PortSearchBoundingBoxUITable.ColumnName = {'X (µm)'; 'Y (µm)'};
            app.PortSearchBoundingBoxUITable.RowName = {};
            app.PortSearchBoundingBoxUITable.ColumnEditable = true;
            app.PortSearchBoundingBoxUITable.CellEditCallback = createCallbackFcn(app, @PortSearchBoundingBoxUITableCellEdit, true);
            app.PortSearchBoundingBoxUITable.Tooltip = {'Data of common port for all wavelengths'};
            app.PortSearchBoundingBoxUITable.Position = [482 636 164 80];

            % Create BoundingBoxforPortSearchLabel
            app.BoundingBoxforPortSearchLabel = uilabel(app.UIFigure);
            app.BoundingBoxforPortSearchLabel.Position = [483 715 164 22];
            app.BoundingBoxforPortSearchLabel.Text = 'Bounding Box for Port Search';

            % Create ShowEllipsesCheckBox
            app.ShowEllipsesCheckBox = uicheckbox(app.UIFigure);
            app.ShowEllipsesCheckBox.Text = 'Show Ellipses';
            app.ShowEllipsesCheckBox.Position = [681 696 97 22];

            % Create ShowCirclesCheckBox
            app.ShowCirclesCheckBox = uicheckbox(app.UIFigure);
            app.ShowCirclesCheckBox.Text = 'Show Circles';
            app.ShowCirclesCheckBox.Position = [681 675 92 22];

            % Create ShowPortSearchIntersectionsCheckBox
            app.ShowPortSearchIntersectionsCheckBox = uicheckbox(app.UIFigure);
            app.ShowPortSearchIntersectionsCheckBox.Text = 'Show Port Search Intersections';
            app.ShowPortSearchIntersectionsCheckBox.Position = [681 654 191 22];

            % Create ShowGratingAdditionsCheckBox
            app.ShowGratingAdditionsCheckBox = uicheckbox(app.UIFigure);
            app.ShowGratingAdditionsCheckBox.Text = 'Show Grating Additions';
            app.ShowGratingAdditionsCheckBox.Position = [681 633 147 22];
            app.ShowGratingAdditionsCheckBox.Value = true;

            % Create BragggratingdefinitionPanel
            app.BragggratingdefinitionPanel = uipanel(app.UIFigure);
            app.BragggratingdefinitionPanel.Title = 'Bragg grating definition';
            app.BragggratingdefinitionPanel.Position = [16 189 218 174];

            % Create PeriodsEditFieldLabel
            app.PeriodsEditFieldLabel = uilabel(app.BragggratingdefinitionPanel);
            app.PeriodsEditFieldLabel.HorizontalAlignment = 'right';
            app.PeriodsEditFieldLabel.Position = [97 123 46 22];
            app.PeriodsEditFieldLabel.Text = 'Periods';

            % Create PeriodsEditField
            app.PeriodsEditField = uieditfield(app.BragggratingdefinitionPanel, 'numeric');
            app.PeriodsEditField.Limits = [1 Inf];
            app.PeriodsEditField.ValueChangedFcn = createCallbackFcn(app, @PeriodsEditFieldValueChanged, true);
            app.PeriodsEditField.Position = [155 123 57 22];
            app.PeriodsEditField.Value = 4;

            % Create PeriodlengthnmEditFieldLabel
            app.PeriodlengthnmEditFieldLabel = uilabel(app.BragggratingdefinitionPanel);
            app.PeriodlengthnmEditFieldLabel.HorizontalAlignment = 'right';
            app.PeriodlengthnmEditFieldLabel.Position = [39 95 104 22];
            app.PeriodlengthnmEditFieldLabel.Text = 'Period length (nm)';

            % Create PeriodlengthEditField
            app.PeriodlengthEditField = uieditfield(app.BragggratingdefinitionPanel, 'numeric');
            app.PeriodlengthEditField.Limits = [0 Inf];
            app.PeriodlengthEditField.ValueDisplayFormat = '%12.12g';
            app.PeriodlengthEditField.ValueChangedFcn = createCallbackFcn(app, @PeriodlengthEditFieldValueChanged, true);
            app.PeriodlengthEditField.Position = [155 95 57 22];
            app.PeriodlengthEditField.Value = 361;

            % Create SiO2lengthnmEditFieldLabel
            app.SiO2lengthnmEditFieldLabel = uilabel(app.BragggratingdefinitionPanel);
            app.SiO2lengthnmEditFieldLabel.HorizontalAlignment = 'right';
            app.SiO2lengthnmEditFieldLabel.Position = [39 66 104 22];
            app.SiO2lengthnmEditFieldLabel.Text = 'SiO2 length (nm)';

            % Create SiO2lengthEditField
            app.SiO2lengthEditField = uieditfield(app.BragggratingdefinitionPanel, 'numeric');
            app.SiO2lengthEditField.Limits = [0 Inf];
            app.SiO2lengthEditField.ValueDisplayFormat = '%12.12g';
            app.SiO2lengthEditField.ValueChangedFcn = createCallbackFcn(app, @SiO2lengthEditFieldValueChanged, true);
            app.SiO2lengthEditField.Position = [155 66 57 22];
            app.SiO2lengthEditField.Value = 180.5;

            % Create BoundarydistancenmLabel
            app.BoundarydistancenmLabel = uilabel(app.BragggratingdefinitionPanel);
            app.BoundarydistancenmLabel.HorizontalAlignment = 'right';
            app.BoundarydistancenmLabel.Position = [10 37 133 22];
            app.BoundarydistancenmLabel.Text = 'Boundary distance (nm)';

            % Create BoundarydistanceEditField
            app.BoundarydistanceEditField = uieditfield(app.BragggratingdefinitionPanel, 'numeric');
            app.BoundarydistanceEditField.Limits = [0 Inf];
            app.BoundarydistanceEditField.ValueDisplayFormat = '%12.12g';
            app.BoundarydistanceEditField.ValueChangedFcn = createCallbackFcn(app, @BoundarydistanceEditFieldValueChanged, true);
            app.BoundarydistanceEditField.Position = [155 37 57 22];
            app.BoundarydistanceEditField.Value = 10000;

            % Create ReflectorspacingnmLabel
            app.ReflectorspacingnmLabel = uilabel(app.BragggratingdefinitionPanel);
            app.ReflectorspacingnmLabel.HorizontalAlignment = 'right';
            app.ReflectorspacingnmLabel.Position = [17 8 126 22];
            app.ReflectorspacingnmLabel.Text = 'Reflector spacing (nm)';

            % Create ReflectorSpacingEditField
            app.ReflectorSpacingEditField = uieditfield(app.BragggratingdefinitionPanel, 'numeric');
            app.ReflectorSpacingEditField.ValueDisplayFormat = '%12.12g';
            app.ReflectorSpacingEditField.ValueChangedFcn = createCallbackFcn(app, @ReflectorSpacingEditFieldValueChanged, true);
            app.ReflectorSpacingEditField.Position = [155 8 57 22];
            app.ReflectorSpacingEditField.Value = 185;

            % Create PortlayoutdefinitionPanel
            app.PortlayoutdefinitionPanel = uipanel(app.UIFigure);
            app.PortlayoutdefinitionPanel.Title = 'Port layout definition';
            app.PortlayoutdefinitionPanel.Position = [244 131 208 232];

            % Create WaveguidewidthnmEditFieldLabel
            app.WaveguidewidthnmEditFieldLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.WaveguidewidthnmEditFieldLabel.HorizontalAlignment = 'right';
            app.WaveguidewidthnmEditFieldLabel.Position = [6 181 124 22];
            app.WaveguidewidthnmEditFieldLabel.Text = 'Waveguide width (nm)';

            % Create WaveguidewidthEditField
            app.WaveguidewidthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.WaveguidewidthEditField.Limits = [1 Inf];
            app.WaveguidewidthEditField.ValueDisplayFormat = '%12.12g';
            app.WaveguidewidthEditField.ValueChangedFcn = createCallbackFcn(app, @WaveguidewidthEditFieldValueChanged, true);
            app.WaveguidewidthEditField.Position = [143 181 57 22];
            app.WaveguidewidthEditField.Value = 400;

            % Create TaperendwidthnmLabel
            app.TaperendwidthnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.TaperendwidthnmLabel.HorizontalAlignment = 'right';
            app.TaperendwidthnmLabel.Position = [12 153 118 22];
            app.TaperendwidthnmLabel.Text = 'Taper end width (nm)';

            % Create TaperendwidthEditField
            app.TaperendwidthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.TaperendwidthEditField.Limits = [1 Inf];
            app.TaperendwidthEditField.ValueDisplayFormat = '%12.12g';
            app.TaperendwidthEditField.ValueChangedFcn = createCallbackFcn(app, @TaperendwidthEditFieldValueChanged, true);
            app.TaperendwidthEditField.Position = [143 153 57 22];
            app.TaperendwidthEditField.Value = 2000;

            % Create TaperlengthnmLabel
            app.TaperlengthnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.TaperlengthnmLabel.HorizontalAlignment = 'right';
            app.TaperlengthnmLabel.Position = [26 124 104 22];
            app.TaperlengthnmLabel.Text = 'Taper length (nm)';

            % Create TaperlengthEditField
            app.TaperlengthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.TaperlengthEditField.Limits = [0 Inf];
            app.TaperlengthEditField.ValueDisplayFormat = '%12.12g';
            app.TaperlengthEditField.ValueChangedFcn = createCallbackFcn(app, @TaperlengthEditFieldValueChanged, true);
            app.TaperlengthEditField.Position = [143 124 57 22];
            app.TaperlengthEditField.Value = 30000;

            % Create WGstublengthnmLabel
            app.WGstublengthnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.WGstublengthnmLabel.HorizontalAlignment = 'right';
            app.WGstublengthnmLabel.Position = [14 95 116 22];
            app.WGstublengthnmLabel.Text = 'WG stub length (nm)';

            % Create WaveguidestublengthEditField
            app.WaveguidestublengthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.WaveguidestublengthEditField.Limits = [0 Inf];
            app.WaveguidestublengthEditField.ValueDisplayFormat = '%12.12g';
            app.WaveguidestublengthEditField.ValueChangedFcn = createCallbackFcn(app, @WaveguidestublengthEditFieldValueChanged, true);
            app.WaveguidestublengthEditField.Position = [143 95 57 22];
            app.WaveguidestublengthEditField.Value = 2000;

            % Create TapershiftnmLabel
            app.TapershiftnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.TapershiftnmLabel.HorizontalAlignment = 'right';
            app.TapershiftnmLabel.Position = [41 66 89 22];
            app.TapershiftnmLabel.Text = 'Taper shift (nm)';

            % Create TapershiftEditField
            app.TapershiftEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.TapershiftEditField.ValueDisplayFormat = '%12.12g';
            app.TapershiftEditField.ValueChangedFcn = createCallbackFcn(app, @TapershiftEditFieldValueChanged, true);
            app.TapershiftEditField.Position = [143 66 57 22];

            % Create TrenchwidthnmLabel
            app.TrenchwidthnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.TrenchwidthnmLabel.HorizontalAlignment = 'right';
            app.TrenchwidthnmLabel.Position = [28 37 102 22];
            app.TrenchwidthnmLabel.Text = 'Trench width (nm)';

            % Create TrenchwidthEditField
            app.TrenchwidthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.TrenchwidthEditField.Limits = [1 Inf];
            app.TrenchwidthEditField.ValueDisplayFormat = '%12.12g';
            app.TrenchwidthEditField.ValueChangedFcn = createCallbackFcn(app, @TrenchwidthEditFieldValueChanged, true);
            app.TrenchwidthEditField.Position = [143 37 57 22];
            app.TrenchwidthEditField.Value = 1000;

            % Create Trench2widthnmLabel
            app.Trench2widthnmLabel = uilabel(app.PortlayoutdefinitionPanel);
            app.Trench2widthnmLabel.HorizontalAlignment = 'right';
            app.Trench2widthnmLabel.Position = [22 8 108 22];
            app.Trench2widthnmLabel.Text = 'Trench2 width (nm)';

            % Create Trench2widthEditField
            app.Trench2widthEditField = uieditfield(app.PortlayoutdefinitionPanel, 'numeric');
            app.Trench2widthEditField.Limits = [1 Inf];
            app.Trench2widthEditField.ValueDisplayFormat = '%12.12g';
            app.Trench2widthEditField.ValueChangedFcn = createCallbackFcn(app, @Trench2widthEditFieldValueChanged, true);
            app.Trench2widthEditField.Position = [143 8 57 22];
            app.Trench2widthEditField.Value = 1000;

            % Create COMSOLModelButton
            app.COMSOLModelButton = uibutton(app.UIFigure, 'push');
            app.COMSOLModelButton.ButtonPushedFcn = createCallbackFcn(app, @COMSOLModelButtonPushed, true);
            app.COMSOLModelButton.BackgroundColor = [0.8392 0.9686 0.8392];
            app.COMSOLModelButton.Position = [760 214 176 32];
            app.COMSOLModelButton.Text = 'COMSOL Model';

            % Create GDS2LayoutButton
            app.GDS2LayoutButton = uibutton(app.UIFigure, 'push');
            app.GDS2LayoutButton.ButtonPushedFcn = createCallbackFcn(app, @GDS2LayoutButtonPushed, true);
            app.GDS2LayoutButton.BackgroundColor = [0.8392 0.9686 0.8392];
            app.GDS2LayoutButton.Position = [760 69 176 32];
            app.GDS2LayoutButton.Text = 'GDS2 Layout';

            % Create COMSOLSimulationPanel
            app.COMSOLSimulationPanel = uipanel(app.UIFigure);
            app.COMSOLSimulationPanel.Title = 'COMSOL Simulation';
            app.COMSOLSimulationPanel.Position = [482 251 454 115];

            % Create MinWavelengthnmLabel
            app.MinWavelengthnmLabel = uilabel(app.COMSOLSimulationPanel);
            app.MinWavelengthnmLabel.HorizontalAlignment = 'right';
            app.MinWavelengthnmLabel.Position = [10 64 122 22];
            app.MinWavelengthnmLabel.Text = 'Min. Wavelength (nm)';

            % Create ComsolMinWavelengthEditField
            app.ComsolMinWavelengthEditField = uieditfield(app.COMSOLSimulationPanel, 'numeric');
            app.ComsolMinWavelengthEditField.Limits = [1 Inf];
            app.ComsolMinWavelengthEditField.ValueDisplayFormat = '%12.12g';
            app.ComsolMinWavelengthEditField.ValueChangedFcn = createCallbackFcn(app, @ComsolMinWavelengthEditFieldValueChanged, true);
            app.ComsolMinWavelengthEditField.Position = [145 64 57 22];
            app.ComsolMinWavelengthEditField.Value = 1500;

            % Create MaxWavelengthnmLabel
            app.MaxWavelengthnmLabel = uilabel(app.COMSOLSimulationPanel);
            app.MaxWavelengthnmLabel.HorizontalAlignment = 'right';
            app.MaxWavelengthnmLabel.Position = [6 35 126 22];
            app.MaxWavelengthnmLabel.Text = 'Max. Wavelength (nm)';

            % Create ComsolMaxWavelengthEditField
            app.ComsolMaxWavelengthEditField = uieditfield(app.COMSOLSimulationPanel, 'numeric');
            app.ComsolMaxWavelengthEditField.Limits = [1 Inf];
            app.ComsolMaxWavelengthEditField.ValueDisplayFormat = '%12.12g';
            app.ComsolMaxWavelengthEditField.ValueChangedFcn = createCallbackFcn(app, @ComsolMaxWavelengthEditFieldValueChanged, true);
            app.ComsolMaxWavelengthEditField.Position = [145 35 57 22];
            app.ComsolMaxWavelengthEditField.Value = 1600;

            % Create WavelengthstepnmLabel
            app.WavelengthstepnmLabel = uilabel(app.COMSOLSimulationPanel);
            app.WavelengthstepnmLabel.HorizontalAlignment = 'right';
            app.WavelengthstepnmLabel.Position = [10 7 122 22];
            app.WavelengthstepnmLabel.Text = 'Wavelength step (nm)';

            % Create ComsolWavelengthStepEditField
            app.ComsolWavelengthStepEditField = uieditfield(app.COMSOLSimulationPanel, 'numeric');
            app.ComsolWavelengthStepEditField.Limits = [0 Inf];
            app.ComsolWavelengthStepEditField.ValueDisplayFormat = '%12.12g';
            app.ComsolWavelengthStepEditField.ValueChangedFcn = createCallbackFcn(app, @ComsolWavelengthStepEditFieldValueChanged, true);
            app.ComsolWavelengthStepEditField.Position = [145 7 57 22];
            app.ComsolWavelengthStepEditField.Value = 2;

            % Create ComsolBorderButton
            app.ComsolBorderButton = uibutton(app.COMSOLSimulationPanel, 'state');
            app.ComsolBorderButton.ValueChangedFcn = createCallbackFcn(app, @ComsolBorderButtonValueChanged, true);
            app.ComsolBorderButton.Text = {'Full'; 'Border'};
            app.ComsolBorderButton.Position = [389 50 53 36];
            app.ComsolBorderButton.Value = true;

            % Create ComsolSBendsButton
            app.ComsolSBendsButton = uibutton(app.COMSOLSimulationPanel, 'state');
            app.ComsolSBendsButton.ValueChangedFcn = createCallbackFcn(app, @ComsolSBendsButtonValueChanged, true);
            app.ComsolSBendsButton.Text = 'S-Bends';
            app.ComsolSBendsButton.FontSize = 11;
            app.ComsolSBendsButton.Position = [215 50 53 36];
            app.ComsolSBendsButton.Value = true;

            % Create FreespaceareaLabel
            app.FreespaceareaLabel = uilabel(app.COMSOLSimulationPanel);
            app.FreespaceareaLabel.Position = [284 64 96 22];
            app.FreespaceareaLabel.Text = 'Free space area:';

            % Create NoofPartitionsEditFieldLabel
            app.NoofPartitionsEditFieldLabel = uilabel(app.COMSOLSimulationPanel);
            app.NoofPartitionsEditFieldLabel.Position = [319 9 91 22];
            app.NoofPartitionsEditFieldLabel.Text = 'No. of Partitions';

            % Create NoOfPartitionsEditField
            app.NoOfPartitionsEditField = uieditfield(app.COMSOLSimulationPanel, 'numeric');
            app.NoOfPartitionsEditField.Limits = [1 Inf];
            app.NoOfPartitionsEditField.ValueChangedFcn = createCallbackFcn(app, @NoOfPartitionsEditFieldValueChanged, true);
            app.NoOfPartitionsEditField.Position = [412 9 30 22];
            app.NoOfPartitionsEditField.Value = 48;

            % Create GDSIIDefinitionsPanel
            app.GDSIIDefinitionsPanel = uipanel(app.UIFigure);
            app.GDSIIDefinitionsPanel.Title = 'GDSII Definitions';
            app.GDSIIDefinitionsPanel.Position = [482 69 268 138];

            % Create EchelleNameIDEditFieldLabel
            app.EchelleNameIDEditFieldLabel = uilabel(app.GDSIIDefinitionsPanel);
            app.EchelleNameIDEditFieldLabel.HorizontalAlignment = 'right';
            app.EchelleNameIDEditFieldLabel.Position = [8 82 96 22];
            app.EchelleNameIDEditFieldLabel.Text = 'Echelle Name/ID';

            % Create EchelleNameIDEditField
            app.EchelleNameIDEditField = uieditfield(app.GDSIIDefinitionsPanel, 'text');
            app.EchelleNameIDEditField.ValueChangedFcn = createCallbackFcn(app, @EchelleNameIDEditFieldValueChanged, true);
            app.EchelleNameIDEditField.Position = [117 82 137 22];
            app.EchelleNameIDEditField.Value = 'EG_01';

            % Create GDSdatabaseunitLabel
            app.GDSdatabaseunitLabel = uilabel(app.GDSIIDefinitionsPanel);
            app.GDSdatabaseunitLabel.HorizontalAlignment = 'right';
            app.GDSdatabaseunitLabel.Position = [77 53 107 22];
            app.GDSdatabaseunitLabel.Text = 'GDS database unit';

            % Create GDS2DBUnitEditField
            app.GDS2DBUnitEditField = uieditfield(app.GDSIIDefinitionsPanel, 'numeric');
            app.GDS2DBUnitEditField.Limits = [0 Inf];
            app.GDS2DBUnitEditField.ValueDisplayFormat = '%12.12g';
            app.GDS2DBUnitEditField.ValueChangedFcn = createCallbackFcn(app, @GDS2DBUnitEditFieldValueChanged, true);
            app.GDS2DBUnitEditField.Position = [197 53 57 22];
            app.GDS2DBUnitEditField.Value = 1e-09;

            % Create GDSuserunitLabel
            app.GDSuserunitLabel = uilabel(app.GDSIIDefinitionsPanel);
            app.GDSuserunitLabel.HorizontalAlignment = 'right';
            app.GDSuserunitLabel.Position = [103 24 81 22];
            app.GDSuserunitLabel.Text = 'GDS user unit';

            % Create GDS2UserUnitEditField
            app.GDS2UserUnitEditField = uieditfield(app.GDSIIDefinitionsPanel, 'numeric');
            app.GDS2UserUnitEditField.Limits = [0 Inf];
            app.GDS2UserUnitEditField.ValueDisplayFormat = '%12.12g';
            app.GDS2UserUnitEditField.ValueChangedFcn = createCallbackFcn(app, @GDS2UserUnitEditFieldValueChanged, true);
            app.GDS2UserUnitEditField.Position = [197 24 57 22];
            app.GDS2UserUnitEditField.Value = 1e-06;

            % Create SbendWaveguidesPanel
            app.SbendWaveguidesPanel = uipanel(app.UIFigure);
            app.SbendWaveguidesPanel.Title = 'S-bend Waveguides';
            app.SbendWaveguidesPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.SbendWaveguidesPanel.Position = [16 65 218 118];

            % Create WaveguidelengthnmLabel
            app.WaveguidelengthnmLabel = uilabel(app.SbendWaveguidesPanel);
            app.WaveguidelengthnmLabel.HorizontalAlignment = 'right';
            app.WaveguidelengthnmLabel.Position = [13 67 129 22];
            app.WaveguidelengthnmLabel.Text = 'Waveguide length (nm)';

            % Create WaveguideLengthEditField
            app.WaveguideLengthEditField = uieditfield(app.SbendWaveguidesPanel, 'numeric');
            app.WaveguideLengthEditField.Limits = [0 Inf];
            app.WaveguideLengthEditField.ValueDisplayFormat = '%12.12g';
            app.WaveguideLengthEditField.ValueChangedFcn = createCallbackFcn(app, @WaveguideLengthEditFieldValueChanged, true);
            app.WaveguideLengthEditField.Position = [155 67 57 22];
            app.WaveguideLengthEditField.Value = 200000;

            % Create EndpointpitchnmLabel
            app.EndpointpitchnmLabel = uilabel(app.SbendWaveguidesPanel);
            app.EndpointpitchnmLabel.HorizontalAlignment = 'right';
            app.EndpointpitchnmLabel.Position = [32 39 110 22];
            app.EndpointpitchnmLabel.Text = 'Endpoint pitch (nm)';

            % Create WaveguidePitchEditField
            app.WaveguidePitchEditField = uieditfield(app.SbendWaveguidesPanel, 'numeric');
            app.WaveguidePitchEditField.Limits = [0 Inf];
            app.WaveguidePitchEditField.ValueDisplayFormat = '%12.12g';
            app.WaveguidePitchEditField.ValueChangedFcn = createCallbackFcn(app, @WaveguidePitchEditFieldValueChanged, true);
            app.WaveguidePitchEditField.Position = [155 38 57 22];
            app.WaveguidePitchEditField.Value = 2500;

            % Create CurvediscretizationLabel
            app.CurvediscretizationLabel = uilabel(app.SbendWaveguidesPanel);
            app.CurvediscretizationLabel.HorizontalAlignment = 'right';
            app.CurvediscretizationLabel.Position = [31 10 111 22];
            app.CurvediscretizationLabel.Text = 'Curve discretization';

            % Create WGCurveDiscretizationEditField
            app.WGCurveDiscretizationEditField = uieditfield(app.SbendWaveguidesPanel, 'numeric');
            app.WGCurveDiscretizationEditField.Limits = [0 Inf];
            app.WGCurveDiscretizationEditField.ValueDisplayFormat = '%12.12g';
            app.WGCurveDiscretizationEditField.ValueChangedFcn = createCallbackFcn(app, @WGCurveDiscretizationEditFieldValueChanged, true);
            app.WGCurveDiscretizationEditField.Position = [155 9 57 22];
            app.WGCurveDiscretizationEditField.Value = 200;

            % Create ProgramSubTitle2Label
            app.ProgramSubTitle2Label = uilabel(app.UIFigure);
            app.ProgramSubTitle2Label.Position = [222 960 118 22];
            app.ProgramSubTitle2Label.Text = 'v0.1';

            % Create BevelingPanel
            app.BevelingPanel = uipanel(app.UIFigure);
            app.BevelingPanel.Title = 'Beveling';
            app.BevelingPanel.Position = [244 10 208 115];

            % Create AttenuatorBevelnmLabel
            app.AttenuatorBevelnmLabel = uilabel(app.BevelingPanel);
            app.AttenuatorBevelnmLabel.HorizontalAlignment = 'right';
            app.AttenuatorBevelnmLabel.Position = [8 9 122 22];
            app.AttenuatorBevelnmLabel.Text = 'Attenuator Bevel (nm)';

            % Create GDS2AttBevelEditField
            app.GDS2AttBevelEditField = uieditfield(app.BevelingPanel, 'numeric');
            app.GDS2AttBevelEditField.Limits = [0 Inf];
            app.GDS2AttBevelEditField.ValueDisplayFormat = '%12.12g';
            app.GDS2AttBevelEditField.ValueChangedFcn = createCallbackFcn(app, @GDS2AttBevelEditFieldValueChanged, true);
            app.GDS2AttBevelEditField.Position = [143 9 57 22];
            app.GDS2AttBevelEditField.Value = 310;

            % Create minAngleLabel
            app.minAngleLabel = uilabel(app.BevelingPanel);
            app.minAngleLabel.HorizontalAlignment = 'right';
            app.minAngleLabel.Position = [68 65 62 22];
            app.minAngleLabel.Text = 'min. Angle';

            % Create GDS2AttBevelAngleEditField
            app.GDS2AttBevelAngleEditField = uieditfield(app.BevelingPanel, 'numeric');
            app.GDS2AttBevelAngleEditField.Limits = [0 180];
            app.GDS2AttBevelAngleEditField.ValueDisplayFormat = '%12.12g';
            app.GDS2AttBevelAngleEditField.ValueChangedFcn = createCallbackFcn(app, @GDS2AttBevelAngleEditFieldValueChanged, true);
            app.GDS2AttBevelAngleEditField.Position = [143 65 57 22];
            app.GDS2AttBevelAngleEditField.Value = 85;

            % Create BevelnmLabel
            app.BevelnmLabel = uilabel(app.BevelingPanel);
            app.BevelnmLabel.HorizontalAlignment = 'right';
            app.BevelnmLabel.Position = [66 37 64 22];
            app.BevelnmLabel.Text = 'Bevel (nm)';

            % Create GDS2BevelEditField
            app.GDS2BevelEditField = uieditfield(app.BevelingPanel, 'numeric');
            app.GDS2BevelEditField.Limits = [0 Inf];
            app.GDS2BevelEditField.ValueDisplayFormat = '%12.12g';
            app.GDS2BevelEditField.ValueChangedFcn = createCallbackFcn(app, @GDS2BevelEditFieldValueChanged, true);
            app.GDS2BevelEditField.Position = [143 37 57 22];
            app.GDS2BevelEditField.Value = 140;

            % Create BlazeshiftmEditFieldLabel
            app.BlazeshiftmEditFieldLabel = uilabel(app.UIFigure);
            app.BlazeshiftmEditFieldLabel.HorizontalAlignment = 'right';
            app.BlazeshiftmEditFieldLabel.Position = [532 753 89 22];
            app.BlazeshiftmEditFieldLabel.Text = 'Blaze shift (µm)';

            % Create BlazeshiftmEditField
            app.BlazeshiftmEditField = uieditfield(app.UIFigure, 'numeric');
            app.BlazeshiftmEditField.ValueChangedFcn = createCallbackFcn(app, @BlazeshiftmEditFieldValueChanged, true);
            app.BlazeshiftmEditField.Position = [636 753 100 22];

            % Create ParallelComputingCheckBox
            app.ParallelComputingCheckBox = uicheckbox(app.UIFigure);
            app.ParallelComputingCheckBox.Text = 'parallel computing';
            app.ParallelComputingCheckBox.Position = [772 581 119 22];
            app.ParallelComputingCheckBox.Value = true;

            % Create CladdingBlockModeCheckBox
            app.CladdingBlockModeCheckBox = uicheckbox(app.UIFigure);
            app.CladdingBlockModeCheckBox.Text = 'Cladding Block Mode';
            app.CladdingBlockModeCheckBox.Position = [97 36 135 22];
            app.CladdingBlockModeCheckBox.Value = true;

            % Create BlockTrenchwidthnmLabel
            app.BlockTrenchwidthnmLabel = uilabel(app.UIFigure);
            app.BlockTrenchwidthnmLabel.HorizontalAlignment = 'right';
            app.BlockTrenchwidthnmLabel.Position = [24 12 134 22];
            app.BlockTrenchwidthnmLabel.Text = 'Block Trench width (nm)';

            % Create BlockModeTrenchwidthEditField
            app.BlockModeTrenchwidthEditField = uieditfield(app.UIFigure, 'numeric');
            app.BlockModeTrenchwidthEditField.Limits = [1 Inf];
            app.BlockModeTrenchwidthEditField.ValueDisplayFormat = '%12.12g';
            app.BlockModeTrenchwidthEditField.ValueChangedFcn = createCallbackFcn(app, @BlockModeTrenchwidthEditFieldValueChanged, true);
            app.BlockModeTrenchwidthEditField.Position = [171 12 57 22];
            app.BlockModeTrenchwidthEditField.Value = 1600;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MarcEchelle_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end