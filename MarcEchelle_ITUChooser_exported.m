classdef MarcEchelle_ITUChooser_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        SubmitButton                  matlab.ui.control.Button
        ReferenceChannelNoofEchellegratingLabel  matlab.ui.control.Label
        ReferenceChannelNo            matlab.ui.control.Spinner
        ITUChannelDropDownLabel       matlab.ui.control.Label
        ITUChannelDropDown            matlab.ui.control.DropDown
        WavelengthDropDownLabel       matlab.ui.control.Label
        WavelengthDropDown            matlab.ui.control.DropDown
        ChannelchooserLabel           matlab.ui.control.Label
        accordingtoITUTG6941102020Label  matlab.ui.control.Label
        ReferenceITUchannelin100GHzgridLabel  matlab.ui.control.Label
        FrequencyDropDownLabel        matlab.ui.control.Label
        FrequencyDropDown             matlab.ui.control.DropDown
        ITUgridDropDownLabel          matlab.ui.control.Label
        ITUgridDropDown               matlab.ui.control.DropDown
        WavelengthTableUITable        matlab.ui.control.Table
        thischannelshouldbecomeLabel  matlab.ui.control.Label
        CancelButton                  matlab.ui.control.Button
    end

    
    properties (Access = private)
        CallingApp   % Main app object
        StartParameters % Parameters from calling function
        NoOfChannels    % available number of channels
        ReferenceChannel % Channel for which you choose an ITU channel and from which the other channels are calculated
        ITUWavelengths  % Wavelengths in ITU-grid for calling app
        ReferenceITUchannel % selected ITU channel number for reference channel
        GridSpacing % selected spacing for channel grid
        ChannelTable    % table with echelle grating channel number as first column and wavelength as second column
    end
    
    methods (Access = private)
        
        function CalcWavelengths(app,RefChannel,RefITUchannel,GridSpacing)
            app.ChannelTable.Wavelength(RefChannel)= (2.99792458e8/( (190.1+(RefITUchannel-1)*0.1)*1e12))*1e9;
            FrequencyList=zeros(app.NoOfChannels,1);
            %FrequencyList(app.ReferenceChannel)=(190.1+(RefITUchannel-1)*0.1)
            %GridSpacing
            for f=1:app.NoOfChannels
                FrequencyList(f)=(190.1+(RefITUchannel-1)*0.1+(f-app.ReferenceChannel)*GridSpacing*1e-3);
                app.ChannelTable.Wavelength(f)=(2.99792458e8/(FrequencyList(f)*1e12))*1e9;
            end
            %FrequencyList
            %app.ChannelTable.Wavelength
        end
        
        function UpdateTableColoring(app)
            addStyle(app.WavelengthTableUITable,uistyle('BackgroundColor',[1 1 1]));
            addStyle(app.WavelengthTableUITable,uistyle('BackgroundColor',[0.8 1.0 0.8]),'cell',[app.ReferenceChannel, 1]);
        end

    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp, myParameters)
            % Store main app object
            app.CallingApp = mainapp;
            app.StartParameters = myParameters;
            app.NoOfChannels=app.StartParameters;
            
            app.ReferenceChannel=1;
            app.ReferenceITUchannel=34;    % default channel near 1550 nm
            app.GridSpacing=800;           % default grid spacing
            
            app.ReferenceChannelNo.Value=app.ReferenceChannel;
            
            %populate the ITU channel chooser drop down menues
            %myarray=cell(1,73)
            for i=1:73
                myChannelArray{1,i}=sprintf('ch %u',i); % Channel numbers
                f=190.1+(i-1)*0.1;
                myFrequencyArray{1,i}=sprintf('%.1f THz',f);    % frequencies
                lambda=(2.99792458e8/(f*1e12))*1e9; %calculate Wavelength in nanometers
                myWavelengthArray{1,i}=sprintf('%.4f nm',lambda);   % wavelengths
            end
            app.ITUChannelDropDown.Items=myChannelArray;
            app.ITUChannelDropDown.ItemsData=(1:73);    % this is to get the item number as value instead of the item's content (5 instead of 'ch 5')
            app.FrequencyDropDown.Items=myFrequencyArray;
            app.FrequencyDropDown.ItemsData=(1:73);
            app.WavelengthDropDown.Items=myWavelengthArray;
            app.WavelengthDropDown.ItemsData=(1:73);
            
            app.ITUChannelDropDown.Value=app.ReferenceITUchannel;
            app.WavelengthDropDown.Value=app.ITUChannelDropDown.Value;
            app.FrequencyDropDown.Value=app.ITUChannelDropDown.Value;
            
            
            app.ITUgridDropDown.Items={'-1600 GHz','-800 GHz','-400 GHz','-200 GHz','-100 GHz','-50 GHz','-25 GHz','-12.5 GHz', '12.5 GHz', '25 GHz', '50 GHz', '100 GHz', '200 GHz', '400 GHz', '800 GHz', '1600 GHz'};
            app.ITUgridDropDown.ItemsData=[-1600, -800, -400, -200, -100, -50, -25, -12.5, 12.5, 25, 50, 100, 200, 400, 800, 1600];
            app.ITUgridDropDown.Value=app.GridSpacing;
            
            app.ChannelTable=table([1:app.NoOfChannels]',zeros(app.NoOfChannels,1));
            app.ChannelTable.Properties.VariableNames={'ChannelNo','Wavelength'};
            app.ChannelTable.Properties.VariableDescriptions={'Channel No.','Wavelength (nm)'};
            %app.ChannelTable
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.ColumnFormat={'shortG','longG'};
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));
            UpdateTableColoring(app);
        end

        % Button pushed function: SubmitButton
        function SubmitButtonPushed(app, event)
            % Call main app's public function
            app.ITUWavelengths=app.ChannelTable.Wavelength(:);
            %app.ITUWavelengths
            UpdateChannelWavelengths(app.CallingApp, app.ITUWavelengths);

            % Delete the dialog box
            delete(app);
        end

        % Value changed function: ReferenceChannelNo
        function ReferenceChannelNoValueChanged(app, event)
            value = app.ReferenceChannelNo.Value;
            if value<=app.NoOfChannels
                app.ReferenceChannel=value;
            else
                app.ReferenceChannelNo.Value=app.NoOfChannels;
            end
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));
            UpdateTableColoring(app);
        end

        % Value changed function: ITUChannelDropDown
        function ITUChannelDropDownValueChanged(app, event)
            value = app.ITUChannelDropDown.Value;
            app.WavelengthDropDown.Value=value;
            app.FrequencyDropDown.Value=value;
            app.ReferenceITUchannel=value;
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));
            
        end

        % Value changed function: WavelengthDropDown
        function WavelengthDropDownValueChanged(app, event)
            value = app.WavelengthDropDown.Value;
            app.ITUChannelDropDown.Value=value;
            app.FrequencyDropDown.Value=value;
            app.ReferenceITUchannel=value;
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));

        end

        % Value changed function: FrequencyDropDown
        function FrequencyDropDownValueChanged(app, event)
            value = app.FrequencyDropDown.Value;
            app.ITUChannelDropDown.Value=value;
            app.WavelengthDropDown.Value=value;
            app.ReferenceITUchannel=value;
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));

        end

        % Value changed function: ITUgridDropDown
        function ITUgridDropDownValueChanged(app, event)
            value = app.ITUgridDropDown.Value;
            app.GridSpacing=value;
            CalcWavelengths(app,app.ReferenceChannel,app.ReferenceITUchannel,app.GridSpacing);
            app.WavelengthTableUITable.Data=table2cell(app.ChannelTable(:,:));

        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            % Call main app's public function
            DontUpdateChannelWavelengths(app.CallingApp);
            % Delete the dialog box
            delete(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 639 721];
            app.UIFigure.Name = 'UI Figure';

            % Create SubmitButton
            app.SubmitButton = uibutton(app.UIFigure, 'push');
            app.SubmitButton.ButtonPushedFcn = createCallbackFcn(app, @SubmitButtonPushed, true);
            app.SubmitButton.BackgroundColor = [0.902 0.9608 0.7961];
            app.SubmitButton.Position = [475 36 141 37];
            app.SubmitButton.Text = 'Submit';

            % Create ReferenceChannelNoofEchellegratingLabel
            app.ReferenceChannelNoofEchellegratingLabel = uilabel(app.UIFigure);
            app.ReferenceChannelNoofEchellegratingLabel.HorizontalAlignment = 'right';
            app.ReferenceChannelNoofEchellegratingLabel.Position = [30 564 229 22];
            app.ReferenceChannelNoofEchellegratingLabel.Text = 'Reference Channel No. of Echelle grating';

            % Create ReferenceChannelNo
            app.ReferenceChannelNo = uispinner(app.UIFigure);
            app.ReferenceChannelNo.Limits = [1 Inf];
            app.ReferenceChannelNo.ValueChangedFcn = createCallbackFcn(app, @ReferenceChannelNoValueChanged, true);
            app.ReferenceChannelNo.Position = [274 564 100 22];
            app.ReferenceChannelNo.Value = 1;

            % Create ITUChannelDropDownLabel
            app.ITUChannelDropDownLabel = uilabel(app.UIFigure);
            app.ITUChannelDropDownLabel.HorizontalAlignment = 'right';
            app.ITUChannelDropDownLabel.Position = [30 498 73 22];
            app.ITUChannelDropDownLabel.Text = 'ITU Channel';

            % Create ITUChannelDropDown
            app.ITUChannelDropDown = uidropdown(app.UIFigure);
            app.ITUChannelDropDown.Items = {'ch 1', 'ch 2', 'ch 3', 'ch 4'};
            app.ITUChannelDropDown.ValueChangedFcn = createCallbackFcn(app, @ITUChannelDropDownValueChanged, true);
            app.ITUChannelDropDown.Position = [118 498 100 22];
            app.ITUChannelDropDown.Value = 'ch 1';

            % Create WavelengthDropDownLabel
            app.WavelengthDropDownLabel = uilabel(app.UIFigure);
            app.WavelengthDropDownLabel.HorizontalAlignment = 'right';
            app.WavelengthDropDownLabel.Position = [417 498 68 22];
            app.WavelengthDropDownLabel.Text = 'Wavelength';

            % Create WavelengthDropDown
            app.WavelengthDropDown = uidropdown(app.UIFigure);
            app.WavelengthDropDown.ValueChangedFcn = createCallbackFcn(app, @WavelengthDropDownValueChanged, true);
            app.WavelengthDropDown.Position = [500 498 110 22];

            % Create ChannelchooserLabel
            app.ChannelchooserLabel = uilabel(app.UIFigure);
            app.ChannelchooserLabel.FontSize = 20;
            app.ChannelchooserLabel.Position = [30 673 157 24];
            app.ChannelchooserLabel.Text = 'Channel chooser';

            % Create accordingtoITUTG6941102020Label
            app.accordingtoITUTG6941102020Label = uilabel(app.UIFigure);
            app.accordingtoITUTG6941102020Label.Position = [30 643 514 22];
            app.accordingtoITUTG6941102020Label.Text = 'according to ITU-T G.694.1 (10/2020) from https://www.itu.int/rec/T-REC-G.694.1-202010-I/en';

            % Create ReferenceITUchannelin100GHzgridLabel
            app.ReferenceITUchannelin100GHzgridLabel = uilabel(app.UIFigure);
            app.ReferenceITUchannelin100GHzgridLabel.Position = [30 622 217 22];
            app.ReferenceITUchannelin100GHzgridLabel.Text = 'Reference ITU channel in 100 GHz grid';

            % Create FrequencyDropDownLabel
            app.FrequencyDropDownLabel = uilabel(app.UIFigure);
            app.FrequencyDropDownLabel.HorizontalAlignment = 'right';
            app.FrequencyDropDownLabel.Position = [233 498 62 22];
            app.FrequencyDropDownLabel.Text = 'Frequency';

            % Create FrequencyDropDown
            app.FrequencyDropDown = uidropdown(app.UIFigure);
            app.FrequencyDropDown.ValueChangedFcn = createCallbackFcn(app, @FrequencyDropDownValueChanged, true);
            app.FrequencyDropDown.Position = [310 498 100 22];

            % Create ITUgridDropDownLabel
            app.ITUgridDropDownLabel = uilabel(app.UIFigure);
            app.ITUgridDropDownLabel.HorizontalAlignment = 'right';
            app.ITUgridDropDownLabel.Position = [30 457 48 22];
            app.ITUgridDropDownLabel.Text = 'ITU grid';

            % Create ITUgridDropDown
            app.ITUgridDropDown = uidropdown(app.UIFigure);
            app.ITUgridDropDown.ValueChangedFcn = createCallbackFcn(app, @ITUgridDropDownValueChanged, true);
            app.ITUgridDropDown.Position = [93 457 100 22];

            % Create WavelengthTableUITable
            app.WavelengthTableUITable = uitable(app.UIFigure);
            app.WavelengthTableUITable.ColumnName = {'Channel #'; 'Wavelength (nm)'};
            app.WavelengthTableUITable.ColumnWidth = {70, 'auto'};
            app.WavelengthTableUITable.RowName = {};
            app.WavelengthTableUITable.Position = [31 36 302 389];

            % Create thischannelshouldbecomeLabel
            app.thischannelshouldbecomeLabel = uilabel(app.UIFigure);
            app.thischannelshouldbecomeLabel.Position = [35 529 158 22];
            app.thischannelshouldbecomeLabel.Text = 'this channel should become:';

            % Create CancelButton
            app.CancelButton = uibutton(app.UIFigure, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.BackgroundColor = [0.9608 0.8118 0.8118];
            app.CancelButton.Position = [475 88 141 37];
            app.CancelButton.Text = 'Cancel';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MarcEchelle_ITUChooser_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

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