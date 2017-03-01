function NARWVocalPitchEnergyV0_4()
% Right whale vocalization quantification

% V0.3: add the frequency range for energy calculation; May 23, 2016

% V0.2: move the PitchContourEnergy to NextSound & Terminate; May 18, 2016


% Version 1.4:
% Add UI:
% Sound file path
% Slider
% save button

% Version 1.3:
% Check Pv; design the way to get the boundaries: t1, t2
% How to make the continuation stronger?
% 

% Revision: use sound files only but not label.csv
% Simply run pitch estimation on every single sound file in a folder

% Version 1: plus formant
% Oct 15, 2015

% Version 0: pitch contour
% Oct 1, 2015
% Author: Yu Shiu

    %clear;
    clear global;
    if ishandle(1)
        close(1)
    end

    HOW_SOUND = 1;
    switch HOW_SOUND
        case 1
            WorkingPlace0 = 'C:\Users\ys587\Google Drive\Right whale calls for Yu Shiu\Focals';
            %WorkingPlace0 = 'C:\ASE_Data\__NARWHolly';
            WorkingPlace = uigetdir(WorkingPlace0,'Select the folder that has Right Whale calls');
            %WorkingPlace = uigetdir(matlabroot,'Select the folder that has Right Whale calls');
        case 2
            CurrFolder = pwd;
            WorkingPlace = uigetdir(CurrFolder,'Select the folder that has Right Whale calls');
    end

    SoundFolder = fullfile(WorkingPlace);
    SoundList = dir(fullfile(SoundFolder,'*.wav'));
    if(size(SoundList,1)==0)
        SoundList = dir(fullfile(SoundFolder,'*.aif'));
    end

    p.fstep = 5;              % frequency resolution of initial spectrogram (Hz)
    %p.fstep = 2; % May 31, 2016
    p.fmax = 4000;            % maximum frequency of initial spectrogram (Hz)
    %p.fmax = 2000;
    p.fres = 10;            % bandwidth of initial spectrogram (Hz)
    %p.fres = 5; % May 31, 2016
    
    p.tinc = 0.01;
    p.flim0 = [20, 500];
    p.flim = [20, 500];
    p.HarmNumMax = round(p.fmax/20.); % MaxF0 = 10Hz

    p.TimeStep = 0.01;
    %TimeStep = 0.05;

    p2.PvThre = 0.4;
    %BufferSize = 0.1;
    %BoundLower = 0.05;
    %BoundUpper = 1.0 - BoundLower;
    p2.MedianHalfWin = 10;

    % Programmatic UI
    %Data = struct('p');

    %global Weight0;
    %global PitchRange
    PitchRange = p.flim;
    Weight0 = [0.6500, 20.0000, 1.000, 10000.0, 500.0, 0.0];
    p2.Weight0 = Weight0;
    %p2.Weight0 = [0.6500, 20.0000, 1.000, 10000.0, 500.0, 0.0];

    FreqRange = 2000; % freq range for spectrogram; display only
    %DEBUG_FLAG = 1;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%  <<<<<<<<<<<<<<<<<<<  DEBUG FLAG

    %% Calculate the feature matrix
    if exist('ax1', 'var')
        cla reset
    else
        f1 = figure('Position',[1, 1, 900, 640]); 
        %f1 = figure('Position',[1, 1, 900, 480]); 
        %f1 = figure(1);
        ax1 = axes('Parent', f1, 'position', [0.05, 0.30, 0.90, 0.65]);
    end

    % Filename display
    UiFilename = uicontrol('Style','text','Position',[20 5 400 20]);
    
    CurrSound =1;
    %global Fx Tx;
    %global IndL IndU;
    %IndLower = 0;
    %IndUpper = 0;

    % "TargetCallCalc": Call TargetCallCalc to estimate pitch contour for the first time
    Pv = [];
    Fv = [];

    MIX0 = [];
    txSpec = [];
    IndL = 0;
    IndU = 0;
    
    EneFreqL = 0;
    EneFreqU = p.fmax;
    
    %MIX0, Fx, IndL, IndU
    %h1 = [];
    
    %global Samples Fs;
    Samples = [];
    Fs = 0;
    %%[Fx, Tx, h1, h2, h3] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);

    % "Score" popupmenu
    %global CallScore;
    CallScore = 0;
    %CallScore = 1;
    b12 = uicontrol('Style','popup','Position',[440,5+25*2,100,20], 'String', {'0: NA', '1: Superb','2: Comme Si Comme Sa','3: Mauvais'} ,'Callback', @StoreScore); 

    % "Suivant": Push button next sound
    %global PitchContour;
    PitchContour = cell(0,13); % Output measurement variables!
    uicontrol('Style','pushbutton','Position',[440+50+2,5+25*1,48,20],'String','Suivant','Callback',{@NextSound, SoundList, SoundFolder, p, p2, f1, ax1, UiFilename, FreqRange});

    % "Re-initialization"
    uicontrol('Style','pushbutton','Position',[440,5+25*1,48,20],'String','Reinitialiser','Callback',{@ResetSound, SoundList, SoundFolder, p, p2, f1, ax1, UiFilename, FreqRange});

    % "Terminer": Write out the pitch contours to a file and finish the work
    uicontrol('Style','pushbutton','Position',[440,5+25*0,100,20],'String','Terminer','Callback',{@TermiTask, SoundList, p, p2});

    % set the energy results global
    %global EneFxTot EneTot;
    EneFxTot = 0; EneTot = 0;

    % "Jouer" button: play sound button
    uicontrol('Style','pushbutton','Position',[440,5+25*3,100,20],'String','Jouer', 'Callback', {@PlayCallback});

    %% ReDrawWeight1
    %global Weight0;
    %PitchRange = p.flim;
    % Slider for weights of pitch estimation algorithm
    uicontrol('Style','text','Position',[20 5+25*3 50 20],'String','Weight 1');
    b1 = uicontrol('Parent',f1,'Style','slider','Position',[80,5+25*3,350,20],'min', 0.0, 'max',30.0, 'Value',15.0);
    b1.Callback = @(es,ed) ReDrawWeight1(1, es.Value);

    uicontrol('Style','text','Position',[20 5+25*2 50 20],'String','Weight 2');
    b2 = uicontrol('Parent',f1,'Style','slider','Position',[80,5+25*2,350,20],'min', 0.0, 'max',20.0, 'Value',10.0);
    b2.Callback = @(es,ed) ReDrawWeight1(2, es.Value);

    uicontrol('Style','text','Position',[20 5+25*1 50 20],'String','Weight 3');
    b3 = uicontrol('Parent',f1,'Style','slider','Position',[80,5+25*1,350,20],'min', 0.0, 'max',20.0, 'Value',10.0);
    b3.Callback = @(es,ed) ReDrawWeight1(3, es.Value);

    %% Currnet pitch range
    b10 = uicontrol('Style','text','Position',[225 5+25*4 30 20],'String',num2str(p.flim(1,1)));
    b11 = uicontrol('Style','text','Position',[260 5+25*4 30 20],'String',num2str(p.flim(1,2)));
    
    %% Redraw Pitch Range
    % Slider for frequency range: low and high
    uicontrol('Style','text','Position',[20 5+25*4 80 20],'String','Pitch Range', 'HorizontalAlignment','left');
    b4 = uicontrol('Parent',f1,'Style','slider','Position',[80,5+25*4,140,20],'min', 0, 'max',500.0, 'Value', p.flim(1,1), 'SliderStep', [1/500, 20/500]);
    FlagPitch = 0;
    b4.Callback = @(es,ed) ReDrawPitchRange(es.Value, FlagPitch);

    b5 = uicontrol('Parent',f1,'Style','slider','Position',[290,5+25*4,140,20],'min', 0, 'max',1000.0, 'Value', p.flim(1,2), 'SliderStep', [1/1000, 20/1000]);
    FlagPitch = 1;
    b5.Callback = @(es,ed) ReDrawPitchRange(es.Value, FlagPitch);

    %% Currnet time boundary
    b8 = uicontrol('Style','text','Position',[225 5+25*5 30 20],'String', num2str(0.0,'%2.3f'));
    b9 = uicontrol('Style','text','Position',[260 5+25*5 30 20],'String', num2str(1.0,'%2.3f'));
    
    %% Adjust Time Boundary: start time and end time
    SoundLength = 100;
    uicontrol('Style','text','Position',[20 5+25*5 80 20],'String','Time Range', 'HorizontalAlignment','left');
    b6 = uicontrol('Parent',f1,'Style','slider','Position',[80,5+25*5,140,20],'min', 0, 'max', SoundLength, 'Value', IndL, 'SliderStep', [1/SoundLength, 10/SoundLength]);
    %FlagTime = 0; % Start Time
    b6.Callback = @(es,ed) RefineTimeBoundary(es.Value, 0);
    
    b7 = uicontrol('Parent',f1,'Style','slider','Position',[290,5+25*5,140,20],'min', 0, 'max', SoundLength, 'Value', IndU, 'SliderStep', [1/SoundLength, 10/SoundLength]);
    %FlagTime = 1; % End time
    b7.Callback = @(es,ed) RefineTimeBoundary(es.Value, 1);
    
    %% Currnet energy freq range
    b13 = uicontrol('Style','text','Position',[660 5+25*5 30 20],'String', num2str(EneFreqL,'%4.0f'));
    b14 = uicontrol('Style','text','Position',[690 5+25*5 30 20],'String', num2str(EneFreqU,'%4.0f'));
    
    %% Freq range for energy measurement
    uicontrol('Style','text','Position',[440 5+25*5 80 20],'String','FreqEne Range', 'HorizontalAlignment','left');
    b15 = uicontrol('Parent',f1,'Style','slider','Position',[520,5+25*5,140,20],'min', 0, 'max', p.fmax, 'Value', EneFreqL, 'SliderStep', [10/p.fmax, 100/p.fmax]);
    b15.Callback = @(es,ed) RefineEneFreqRange(es.Value, 0);
    b16 = uicontrol('Parent',f1,'Style','slider','Position',[720,5+25*5,140,20],'min', 0, 'max', p.fmax, 'Value', EneFreqU, 'SliderStep', [10/p.fmax, 100/p.fmax]);
    b16.Callback = @(es,ed) RefineEneFreqRange(es.Value, 1);
    
    % Run for the 1st time
    %[Fx, Tx, h1, h2, h3] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
    [Fx, Tx, h1, h2, h3, h4, h5] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
    
    %% Callback functions
    function StoreScore(src, evt)
        %global CallScore;
        CallScore = evt.Source.Value - 1;
        %print CallScore;
    end

    function ResetSound(src, evt, SoundList, SoundFolder, p, p2, f1, ax1, UiFilename, FreqRange)
        %[Fx, Tx, h1, h2, h3] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
        [Fx, Tx, h1, h2, h3, h4, h5] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
    end

    % Finish the current sound and move to the next one
    function NextSound(src, evt, SoundList, SoundFolder, p, p2, f1, ax1, UiFilename, FreqRange)
        CallMeasurement;
        
        % Run the next one
        if(CurrSound+1 <= size(SoundList,1))
            CurrSound = CurrSound + 1;
            %[Fx, Tx, h1, h2, h3] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
            [Fx, Tx, h1, h2, h3, h4, h5] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange);
        else
            % Terminate the session
            TermiTaskCore;
        end
    end

    function TermiTask(src, evt, SoundList, p, p2)
        CallMeasurement;
        
        TermiTaskCore;
        if 0
        CurrFolder = pwd;
        MatPath = uigetdir(CurrFolder, 'Select the Output Folder');
        OutputFileName = strcat('NARWPitch-',char(datestr(now,'yyyymmdd-HHMM')), '.mat');
        save(fullfile(MatPath, OutputFileName), 'PitchContour');
        end
    end

    function TermiTaskCore
        CurrFolder = pwd;
        MatPath = uigetdir(CurrFolder, 'Terminer: select the output folder');
        OutputFileName = strcat('NARWPitch-',char(datestr(now,'yyyymmdd-HHMM')), '.mat');
        save(fullfile(MatPath, OutputFileName), 'PitchContour');
    end

    function CallMeasurement
        % Calculate PitchContourEnergy
        %%[IndL, IndU] = SegIntervalEst(Pv, p2.PvThre, p2.MedianHalfWin);
        %[EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, IndL, IndU, p);
        %[EneFxTot, EneTot] = PitchContourEnergy;
        %[EneFxTot, EneTot, FxPeak, FxPeakPoint] = PitchContourEnergy;
        [EneFxTot, EneTot, FxPeak, FxPeakPoint, FxPeakHarmon, FxPeakPointHarmon, F0Bandwidth] = PitchContourEnergy

        % Peak freq measurement: 
        %FxPeak = PeakFreqContour;
        
        IndLower = txSpec(IndL,1);
        IndUpper = txSpec(IndU,1);

        PitchContour{CurrSound,1} = SoundList(CurrSound,1).name;
        PitchContour{CurrSound,2} = CallScore;
        PitchContour{CurrSound,3} = Fx;
        PitchContour{CurrSound,4} = Tx;
        PitchContour{CurrSound,5} = [IndLower, IndUpper];
        PitchContour{CurrSound,6} = [IndL, IndU];
        PitchContour{CurrSound,7} = [EneFxTot, EneTot];
        PitchContour{CurrSound,8} = Fx(IndL:IndU,1); % Fx between time bars
        PitchContour{CurrSound,9} = FxPeak; % Peak frequency contour
        PitchContour{CurrSound,10} = FxPeakPoint; % Peak single frequency
        PitchContour{CurrSound,11} = FxPeakHarmon; % Peak frequency contour from harmonics
        PitchContour{CurrSound,12} = FxPeakPointHarmon; % Peak single frequency from harmonics
        PitchContour{CurrSound,13} = F0Bandwidth; % Bandwidth around F0
        
        PitchContour{CurrSound,14} = EneFxTot; %         
        PitchContour{CurrSound,15} = EneTot; %
        PitchContour{CurrSound,16} = EneFxTot/EneTot; % 
    end
    
    function ReDrawWeight1(VarIndex, SliderVal1)
        % Re-calculate the new pitch contour
        %global Weight0
        switch VarIndex
            case 1
                Weight0(VarIndex) = SliderVal1/1.0;
            case 2
                Weight0(VarIndex) = SliderVal1/1.0;
            case 3
                Weight0(VarIndex) = SliderVal1/1.0;
        end

        %global Fx Tx;
        %H.XDataSource = 'Tx';
        %H.YDataSource = 'Fx';
        h1.XDataSource = 'Tx';
        h1.YDataSource = 'Fx';
        %global PitchRange
        p.flim = PitchRange;

        if PitchRange(1,1) < PitchRange(1,2)
            [Fx, Tx, Pv, Fv] = fxpefacYuV3(Samples, Fs, p.TimeStep,'',p, Weight0);
            refreshdata(h1, 'caller')

            %PitchContourEnergy
            %[EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, Pv, p, p2);
            %%[IndL, IndU] = SegIntervalEst(Pv, p2.PvThre, p2.MedianHalfWin);
            %%[EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, IndL, IndU, p, p2);
        else
            disp('Pitch range low needs to be smaller than pitch range high!')
        end    
    end

    function ReDrawPitchRange(PitchNew, FlagPitch)
        % FlagPitch:
        % 0: low 
        % 1: high
        MinFreqSep = 40; % 40Hz 
        if FlagPitch == 0
            if(PitchNew + MinFreqSep < PitchRange(1,2))
                PitchRange(1,1) = PitchNew;
                
                h1.XDataSource = 'Tx';
                h1.YDataSource = 'Fx';
                p.flim = PitchRange;
                [Fx, Tx, Pv, Fv] = fxpefacYuV3(Samples, Fs, p.TimeStep,'',p, Weight0);
                refreshdata(h1, 'caller')
            else
                set(b4, 'Value', PitchRange(1,1)); % reset to previous value
            end
            set(b10, 'String', PitchRange(1,1));
        else
            if(PitchRange(1,1) + MinFreqSep < PitchNew)
                PitchRange(1,2) = PitchNew;
                
                h1.XDataSource = 'Tx';
                h1.YDataSource = 'Fx';
                p.flim = PitchRange;
                [Fx, Tx, Pv, Fv] = fxpefacYuV3(Samples, Fs, p.TimeStep,'',p, Weight0);
                refreshdata(h1, 'caller')
            else
                set(b5, 'Value', PitchRange(1,2)); % reset to previous value
            end
            set(b11, 'String', PitchRange(1,2));
        end
    end

    function RefineTimeBoundary(TimeNew, FlagTime)
        TimeNew = round(TimeNew);
        if FlagTime == 0 % Start time
            if(TimeNew < IndU)
                IndL = TimeNew;
                %uicontrol('Style','text','Position',[225 5+25*5 30 20],'String',num2str(Tx(IndL, 1),'%2.3f'));
                
                h2.XDataSource = 'XDataTemp';
                h2.YDataSource = 'YDataTemp';
                XDataTemp = Tx(IndL,1)*ones(1,p.fmax);
                YDataTemp = (1:p.fmax);
                refreshdata(h2, 'caller')
            else
                set(b6, 'Value', IndL); % reset to previous value
                %uicontrol('Style','text','Position',[225 5+25*5 30 20],'String',num2str(Tx(IndL, 1),'%2.3f'));
            end
            set(b8, 'String', num2str(Tx(IndL, 1),'%2.2f'));
        else % FlagTime == 1 % end time
            if( IndL < TimeNew)
                IndU = TimeNew;
                %uicontrol('Style','text','Position',[260 5+25*5 30 20],'String',num2str(Tx(IndU, 1),'%2.3f'));
                
                h3.XDataSource = 'XDataTemp';
                h3.YDataSource = 'YDataTemp';
                XDataTemp = Tx(IndU,1)*ones(1,p.fmax);
                YDataTemp = (1:p.fmax);
                refreshdata(h3, 'caller')
            else
                set(b7, 'Value', IndU);
                %uicontrol('Style','text','Position',[260 5+25*5 30 20],'String',num2str(Tx(IndU, 1),'%2.3f'));
            end
            set(b9, 'String', num2str(Tx(IndU, 1),'%2.2f'));
        end
    end

    function RefineEneFreqRange(EneFreqNew, FlagEneFreq)
        EneFreqNew = round(EneFreqNew);
        if FlagEneFreq == 0 % Low freq
            if(EneFreqNew < EneFreqU)
                EneFreqL = EneFreqNew;
                
                h4.XDataSource = 'XDataTemp';
                h4.YDataSource = 'YDataTemp';
                XDataTemp = Tx;
                YDataTemp = EneFreqL*ones(length(Tx),1);
                refreshdata(h4, 'caller')
            else
                set(b15, 'Value', EneFreqL); % reset to previous value
            end
            set(b13, 'String', num2str(EneFreqL,'%4.0f'));
        else % FlagEneFreq == 1 % high freq
            if( IndL < EneFreqNew)
                EneFreqU = EneFreqNew;
                
                h5.XDataSource = 'XDataTemp';
                h5.YDataSource = 'YDataTemp';
                XDataTemp = Tx;
                YDataTemp = EneFreqU*ones(length(Tx),1);
                refreshdata(h5, 'caller')
            else
                set(b16, 'Value', EneFreqU);
            end
            set(b14, 'String', num2str(EneFreqU,'%4.0f'));
        end
    end

    function PlayCallback(src, evt)
        soundsc(Samples, Fs);
    end

    %function [IndLower, IndUpper] = SegIntervalEst(Pv, txSpec, PvThre, MedianHalfWin)
    function [IndL, IndU] = SegIntervalEst(Pv, PvThre, MedianHalfWin)
        % Use Pv to get the boundaries t1 & t2
        %PvThre = 0.4;
        %PvInd = find(Pv > PvThre);

        %MedianWin = MedianHalfWin*2+1;
        Pv2 = zeros(length(Pv),1);
        for pp = MedianHalfWin+1:length(Pv)-MedianHalfWin
           Pv2(pp,1) =  median(Pv(pp-MedianHalfWin:pp+MedianHalfWin,1));
        end
        Pv2Ind = find(Pv2 > PvThre);

        % Pitch estimation method
        SegEnd = find(diff(Pv2Ind)-1>0);
        if(isempty(Pv2Ind))
            IndL = 1;
            IndU = 2;
            disp('No voice in this clip!!')
        elseif(isempty(SegEnd))
            %IndL = txSpec(min(Pv2Ind),1);
            %IndU = txSpec(max(Pv2Ind),1);
            IndL = min(Pv2Ind);
            IndU = max(Pv2Ind);
        elseif(length(SegEnd)>=1)
            SegTime = [Pv2Ind(1,1),Pv2Ind(SegEnd(1,1),1)]; % first segment
            tt1 = 1;
            while(tt1 < length(SegEnd))
                SegTime = [SegTime; [Pv2Ind(SegEnd(tt1,1)+1,1), Pv2Ind(SegEnd(tt1+1,1),1)]];
                tt1 = tt1 + 1;
            end
            SegTime = [SegTime; [Pv2Ind(SegEnd(tt1,1)+1,1), max(Pv2Ind)]];
            SegDur = diff(SegTime,[],2);
            [~, SegDurInd] = max(SegDur);
            %IndL = txSpec(SegTime(SegDurInd,1),1);
            %IndU = txSpec(SegTime(SegDurInd,2),1);
            IndL = SegTime(SegDurInd,1);
            IndU = SegTime(SegDurInd,2);
        end
    end

    %function [Fx, Tx, h1, h2, h3] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange)
    function [Fx, Tx, h1, h2, h3, h4, h5] = TargetCallCalc(SoundList, SoundFolder, CurrSound, p, p2, f1, UiFilename, FreqRange)
        Filename0 = SoundList(CurrSound,1).name;
        Filename = char(fullfile(SoundFolder,Filename0 ));
        set(UiFilename, 'HorizontalAlignment', 'left')
        set(UiFilename, 'String', strcat('Call ', num2str(CurrSound),': ',Filename0));

        fprintf('%d: %s\n', CurrSound, Filename0);

        [Samples, Fs] = audioread(Filename);
        Samples = Samples(:,1); % Use only the first channel if stereo recordings ====>>>> need change!

        fmin = 0; fstep = p.fstep; fmax = p.fmax;
        fres = p.fres;

        % Spectrogram
        [txSpec,f,MIX0]=spgrambw(Samples,Fs,fres,[fmin fstep fmax],[],p.tinc); % frm VoiceBox
        % Pitch estimation
        [Fx, Tx, Pv, Fv] = fxpefacYuV3(Samples, Fs, p.TimeStep,'',p, p2.Weight0);
        % Interval estimation

        %global IndL IndU;
        %[IndLower, IndUpper] = SegIntervalEst(Pv, txSpec, p2.PvThre, p2.MedianHalfWin);
        [IndL, IndU] = SegIntervalEst(Pv, p2.PvThre, p2.MedianHalfWin);
        IndLower = txSpec(IndL,1);
        IndUpper = txSpec(IndU,1);

        %PitchContourEnergy
        %global EneFxTot EneTot;
        %[EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, Pv, p, p2);
        %%[EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, IndL, IndU, p, p2);

        imagesc(txSpec, f, MIX0'.^.1); axis xy; 
        hold on;
        % draw pitch contour
        h1 = plot(Tx, Fx,'r');
        %Fx2 = [Fx, Fx*2, Fx*3, Fx*4, Fx*5];
        %h1 = plot(Tx, Fx2,'r');
        % draw lower time bar
        h2 = plot(IndLower*ones(1,p.fmax),1:p.fmax,'k-');
        % draw higher time bar
        h3 = plot(IndUpper*ones(1,p.fmax),1:p.fmax,'k-');
        
        % draw lower energy freq bar
        h4 = plot(Tx, 1.0*ones(length(Tx),1), 'k--');
        % draw upper energy freq bar
        h5 = plot(Tx, p.fmax*ones(length(Tx),1),'k--');
        hold off;
        %ylim([0, p.fmax]);
        ylim([0, FreqRange]);
        
        %% Initalize values for new sound clip
        CallScore = 0;
        %set(b12, 'String', {'1: Superb','2: Comme Si Comme Sa','3: Mauvais'});
        set(b12, 'String', {'Choose Score:','1: Superb','2: Comme Si Comme Sa','3: Mauvais'});
        set(b12, 'Value', 1);
        
        % Weights
        % Reset ReDrawWeight1 & ReDrawPitchRange to initial values
        set(b1, 'Value', 15.0);
        set(b2, 'Value', 10.0);
        set(b3, 'Value', 10.0);

        %set(b4, 'Value', p.flim0(1,1));
        %set(b5, 'Value', p.flim0(1,2));
        
        % Pitch range
        set(b10, 'String',num2str(p.flim(1,1)));
        set(b11, 'String',num2str(p.flim(1,2)));
        set(b4, 'Value', p.flim(1,1));
        set(b5, 'Value', p.flim(1,2));
        
        % time boundaries
        set(b8, 'String', num2str(Tx(IndL, 1),'%2.3f'));
        set(b9, 'String', num2str(Tx(IndU, 1),'%2.3f'));
        
        SoundLength0 = length(Fx);
        set(b6, 'max', SoundLength0);
        set(b6, 'SliderStep', [1/SoundLength0, 10/SoundLength0]);
        set(b6, 'Value', IndL);
        set(b7, 'max', SoundLength0);
        set(b7, 'SliderStep', [1/SoundLength0, 10/SoundLength0]);
        set(b7, 'Value', IndU);
        
        % Freq range of energy measurement
        set(b13, 'String', num2str(0,'%4.0f'));
        set(b14, 'String', num2str(p.fmax,'%4.0f'));
        
        set(b15, 'max', p.fmax);
        set(b15, 'SliderStep', [10/p.fmax, 100/p.fmax]);
        set(b15, 'Value', 0);
        set(b16, 'max', p.fmax);
        set(b16, 'SliderStep', [10/p.fmax, 100/p.fmax]);
        set(b16, 'Value', p.fmax);
    end

%     %function [EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, Pv, p, p2)
%     function [EneFxTot, EneTot] = PitchContourEnergy(MIX0, Fx, IndL, IndU, p, p2)
%         %[IndL, IndU] = SegIntervalEst(Pv, p2.PvThre, p2.MedianHalfWin);
%         MIX1 = MIX0(IndL:IndU, :);
%         %figure; imagesc(txSpec(IndL:IndU,1), f, MIX1'.^.1);axis xy;
% 
%         %% Energy calculation, given spectrogram MIX1 and pitch contour F0 sequence
%         %imagesc(txSpec, f, MIX1'.^.1); axis xy; 
%         %hold on;
%         HarmNum = 1;
%         FxMax = max(Fx(IndL:IndU,:));
%         EneFxTot = 0.0;
%         EneTot = sum(MIX1(:));
%         %FreqWinHalf = round((round(Fx(IndL:IndU,:)/p.fstep)+1)*0.1); % window size: 20% of a harmonics; 10% on one side
%         FreqWinHalf = ceil((ceil(Fx(IndL:IndU,:)/p.fstep)+1)*0.1); % window size: 20% of a harmonics; 10% on one side
%         %%while(FxMax*HarmNum < p.fmax) % try all harmonics until it hits the freq_max_range
%         while(FxMax*HarmNum < p.fmax && HarmNum <= p.HarmNumMax)
%             FxInd = round(Fx(IndL:IndU,:)*HarmNum/p.fstep)+1;
%             %h1 = plot(Tx(IndL:IndU,1), (FxInd-1)*fstep,'r');
%             EneFx = 0.0;
% 
%             for tt = 1:length(FxInd)
%                 %EneFx = EneFx + MIX1(tt, FxInd(tt,1));
%                 if (FxInd(tt,1) - FreqWinHalf(tt,1) <= 0)
%                     FxIndLow = 1;
%                 else
%                     FxIndLow = FxInd(tt,1) - FreqWinHalf(tt,1);
%                 end
%                 if (FxInd(tt,1)+FreqWinHalf(tt,1) > size(MIX1,2))
%                     FxIndHigh = size(MIX1,2);
%                 else
%                     FxIndHigh = FxInd(tt,1)+FreqWinHalf(tt,1);
%                 end
%                 %MIX1Temp = MIX1(tt, FxInd(tt,1)-FreqWinHalf:FxInd(tt,1)+FreqWinHalf);
%                 MIX1Temp = MIX1(tt, FxIndLow:FxIndHigh);
%                 EneFx = EneFx + sum(MIX1Temp(:));
%             end
%             EneFxTot = EneFxTot + EneFx;
%             %fprintf('Harmonic energy percentage: %.3f as the harmonic number %d\n', EneFxTot/EneTot, HarmNum);
%             HarmNum = HarmNum + 1;
%         end
%         fprintf('Harmonic energy percentage: %.3f as the harmonic number %d\n', EneFxTot/EneTot, HarmNum-1);
%         %fprintf('\n');
%     end

    %function [EneFxTot, EneTot] = PitchContourEnergy
    function [EneFxTot, EneTot, FxPeak, FxPeakPoint, FxPeakHarmon, FxPeakPointHarmon, F0Bandwidth] = PitchContourEnergy
        % EneFxTot: total energy from fundamental frequency F0 and its harmonics
        % EneTot: total energy
        % FxPeak: peak frequency contour
        % FxPeakPoint: peak single frequency
        
        MIX1 = MIX0(IndL:IndU, 1:round(EneFreqU/p.fstep)+1);

        % Peak freq contour: The contour of freqs with the largest amplitude; doing arg max on spectrum
        [ValFxPeak, IndFxPeak] = max(MIX1, [],2); 
        FxPeak = IndFxPeak*p.fstep;
        
        % Paek single freq: The only one freq with the largest amplitude
        [~, IndT1] = max(ValFxPeak);
        FxPeakPoint = IndFxPeak(IndT1,1)*p.fstep;
        
        % The same with FxPeak but the freqs come from the pitch or its harmonics
        %FreqWinHalf = ceil((ceil(Fx(IndL:IndU,:)/p.fstep)+1)*0.1);
        %MIXMask = zeros(size(MIX1));
        
        %% Energy calculation, given spectrogram MIX1 and pitch contour F0 sequence
        if 0
        HarmNum = 1;
        FxMax = max(Fx(IndL:IndU,:));
        
        %EneFxTot = 0.0;
        %EneTot = sum(MIX1(:));

        %FreqWinHalf = ceil((ceil(Fx(IndL:IndU,:)/p.fstep)+1)*0.1); % window size: 20% of a harmonics; 10% on one side
        while(FxMax*HarmNum < EneFreqU && HarmNum <= p.HarmNumMax)
            FxInd = round(Fx(IndL:IndU,:)*HarmNum/p.fstep)+1;
            %EneFx = 0.0;

            for tt = 1:length(FxInd)
                %EneFx = EneFx + MIX1(tt, FxInd(tt,1));
                if (FxInd(tt,1) - FreqWinHalf(tt,1) <= 0)
                    FxIndLow = 1;
                else
                    FxIndLow = FxInd(tt,1) - FreqWinHalf(tt,1);
                end
                if (FxInd(tt,1)+FreqWinHalf(tt,1) > size(MIX1,2))
                    FxIndHigh = size(MIX1,2);
                else
                    FxIndHigh = FxInd(tt,1)+FreqWinHalf(tt,1);
                end
                %MIX1Temp = MIX1(tt, FxIndLow:FxIndHigh);
                %EneFx = EneFx + sum(MIX1Temp(:));
                MIXMask(tt, FxIndLow:FxIndHigh) = 1;
            end
            %EneFxTot = EneFxTot + EneFx;
            %fprintf('Harmonic energy percentage: %.3f as the harmonic number %d\n', EneFxTot/EneTot, HarmNum);
            HarmNum = HarmNum + 1;
        end
        %figure; imagesc(MIXMask'); axis xy; grid;
        %disp('');
        
        MIX2 = MIX1.*MIXMask;
        EneFxTot = sum(sum(MIX2,1),2);
        EneTot = sum(MIX1(:));
        
        % Peak freq contour from harmonics:
        [ValFxPeak2, IndFxPeak2] = max(MIX2, [],2); 
        FxPeakHarmon = IndFxPeak2*p.fstep;
        
        % Paek single freq: The only one freq with the largest amplitude
        [~, IndT2] = max(ValFxPeak2);
        FxPeakPointHarmon = IndFxPeak2(IndT2,1)*p.fstep;
        end
        
        if 1
            EneTot = sum(MIX1(:));
            
            WinDim = max(ceil((ceil(Fx(IndL:IndU,:)/p.fstep)+1)*0.2)); % 20%
            
            EneRatio = zeros(WinDim+1, 1);
            EneFxTotCand = zeros(WinDim+1, 1);
            FxMax = max(Fx(IndL:IndU,:));
            FxPeakHarmonCell = cell(WinDim+1,1);
            FxPeakPointHarmonCell = cell(WinDim+1,1);
            for FreqWinHalf = 0:WinDim
                MIXMask = zeros(size(MIX1));
                HarmNum = 1;
                while(FxMax*HarmNum < EneFreqU && HarmNum <= p.HarmNumMax)
                    FxInd = round(Fx(IndL:IndU,:)*HarmNum/p.fstep)+1;
                    %EneFx = 0.0;

                    for tt = 1:length(FxInd)
                        %EneFx = EneFx + MIX1(tt, FxInd(tt,1));
                        if (FxInd(tt,1) - FreqWinHalf <= 0)
                            FxIndLow = 1;
                        else
                            FxIndLow = FxInd(tt,1) - FreqWinHalf;
                        end
                        if (FxInd(tt,1)+FreqWinHalf > size(MIX1,2))
                            FxIndHigh = size(MIX1,2);
                        else
                            FxIndHigh = FxInd(tt,1)+FreqWinHalf;
                        end
                        %MIX1Temp = MIX1(tt, FxIndLow:FxIndHigh);
                        %EneFx = EneFx + sum(MIX1Temp(:));
                        MIXMask(tt, FxIndLow:FxIndHigh) = 1;
                    end
                    %EneFxTot = EneFxTot + EneFx;
                    %fprintf('Harmonic energy percentage: %.3f as the harmonic number %d\n', EneFxTot/EneTot, HarmNum);
                    HarmNum = HarmNum + 1;
                end
                MIX3 = MIX1.*MIXMask;
                EneFxTotCand(FreqWinHalf+1,1) = sum(sum(MIX3,1),2);
                EneRatio(FreqWinHalf+1,1) = EneFxTotCand(FreqWinHalf+1,1)/EneTot;
                %EneTot = sum(MIX1(:));
                
                % Peak freq contour from harmonics:
                [ValFxPeak2, IndFxPeak2] = max(MIX3, [],2); 
                FxPeakHarmonCell{FreqWinHalf+1,1} = IndFxPeak2*p.fstep;

                % Paek single freq: The only one freq with the largest amplitude
                [~, IndT2] = max(ValFxPeak2);
                FxPeakPointHarmonCell{FreqWinHalf+1,1} = IndFxPeak2(IndT2,1)*p.fstep;
            end
            EneRatioDiff = diff(EneRatio,1)/2; % divided by 2 since the window is symmetric
            EneF0Thre = EneRatio(1,1)*.50; % 50% of the energy of F0
            bb = 1;
            while(bb < WinDim && EneRatioDiff(bb,1)>EneF0Thre )
               bb = bb + 1; 
            end
            bb = bb - 1; % get bb back lower than the threshold EneF0Thre
            EneFxTot = EneFxTotCand(bb+1,1);
            F0Bandwidth = bb*2*p.fstep;
            FxPeakHarmon = FxPeakHarmonCell{bb+1,1};
            FxPeakPointHarmon = FxPeakPointHarmonCell{bb+1,1}; 
        end
        
        fprintf('EneFxTot: %.3f; EneTot: %.3f\n', EneFxTot, EneTot);
        %fprintf('Harmonic energy percentage: %.3f as the harmonic number %d\n', EneFxTot/EneTot, HarmNum-1);
        %fprintf('\n');
    end

    function [FxPeak, FxPeakPoint] = PeakFreqContour
        %MIX1 = MIX0(IndL:IndU, :);
        MIX1 = MIX0(IndL:IndU, 1:round(EneFreqU/p.fstep)+1);
        
        % Peak freq contour: The contour of freqs with the largest amplitude; doing arg max on spectrum
        [ValFxPeak, IndFxPeak] = max(MIX1, [],2); 
        FxPeak = IndFxPeak*p.fstep;
        
        % Paek single freq: The only one freq with the largest amplitude
        [ValT1, IndT1] = max(ValFxPeak);
        FxPeakPoint = IndFxPeak(IndT1,1)*p.fstep;
        
        % The same with FxPeak but the freqs come from the pitch or its harmonics
        FreqWinHalf = ceil((ceil(Fx(IndL:IndU,:)/p.fstep)+1)*0.1);
        MIXMask = zeros(size(MIX1));
        
        HarmNum = 1;
        FxMax = max(Fx(IndL:IndU,:));
        while(FxMax*HarmNum < EneFreqU && HarmNum <= p.HarmNumMax)
            FxInd = round(Fx(IndL:IndU,:)*HarmNum/p.fstep)+1;
            
            for tt = 1:length(FxInd)
                %EneFx = EneFx + MIX1(tt, FxInd(tt,1));
                if (FxInd(tt,1) - FreqWinHalf(tt,1) <= 0)
                    FxIndLow = 1;
                else
                    FxIndLow = FxInd(tt,1) - FreqWinHalf(tt,1);
                end
                if (FxInd(tt,1)+FreqWinHalf(tt,1) > size(MIX1,2))
                    FxIndHigh = size(MIX1,2);
                else
                    FxIndHigh = FxInd(tt,1)+FreqWinHalf(tt,1);
                end
                MIXMask(tt, FxIndLow:FxIndHigh) = 1;
            end
            HarmNum = HarmNum + 1;
        end
        
        %figure; imagesc(MIXMask'); axis xy; grid;
        disp('');
    end 
end













