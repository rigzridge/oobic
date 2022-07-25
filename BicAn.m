%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%    ____                 ______             
%   /\  _`\    __        /\  _  \                       Version 4.00
%   \ \ \L\ \ /\_\    ___\ \ \L\ \    ___    
%    \ \  _ <'\/\ \  /'___\ \  __ \ /' _ `\          (c) 2022 -- G.Riggs
%     \ \ \L\ \\ \ \/\ \__/\ \ \/\ \/\ \/\ \ 
%      \ \____/ \ \_\ \____\\ \_\ \_\ \_\ \_\      WVU Physics & Astronomy
%       \/___/   \/_/\/____/ \/_/\/_/\/_/\/_/      
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Object-oriented bicoherence analysis and signal processing toolbox
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% The Bispectrum
% B_xyz(f1,f2) = < X(f1)Y(f2)Z(f1+f2)* >, where x,y,z are time series with 
% corresponding Fourier transforms X,Y,Z, and <...> denotes averaging.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% The (squared) Bicoherence spectrum
% b^2_xyz(f1,f2) =           |B_xyz(f1,f2)|^2
%                          --------------------
%                ( <|X(f1)Y(f2)|^2> <|Z(f1+f2)|^2> + eps ),
% where eps is a small number meant to prevent 0/0 = NaN catastrophe
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Inputs
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% inData    -> time-series {or structure}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% additional options... (see below for instructions)
% - - - - - - - - - - - - - - - - - - - - 
% autoscale -> autoscaling in figures                  [default :: False]
% bispectro -> computes bispectrogram                  x[default :: False]
% cbarnorth -> control bolorbar location               [default :: True]
% cmap      -> adjust colormap                         [default :: 'viridis']
% dealias   -> applies antialiasing (LP) filter        x[default :: False]
% detrend   -> remove linear trend from data           [default :: False]
% errlim    -> mean(fft) condition                     [default :: 1e15] 
% filter    -> xxxxxxxxxxxxxxx                         x[default :: 'none']
% freqres   -> desired frequency resolution [Hz]       [default :: 0]
% fscale    -> scale for plotting frequencies          [default :: 0]
% justspec  -> true for just spectrogram               [default :: False]
% lilguy    -> set epsilon                             [default :: 1e-6]
% note      -> optional string for documentation       [default :: ' '] 
% plotit    -> start plotting tool when done           [default :: False]
% plottype  -> set desired plottable                   [default :: 'bicoh']
% samprate  -> sampling rate in Hz                     [default :: 1]
% sigma     -> parameter for wavelet spectrum          [default :: 0]
% spectype  -> set desired time-freq. method           [default :: 'stft']
% step      -> step size for Welch method in samples   [default :: 512]
% subint    -> subinterval size in samples             [default :: 128]
% sizewarn  -> warning for matrix size                 x[default :: True]
% smooth    -> smooths FFT by n samples                [default :: 1]
% tscale    -> scale for plotting time                 [default :: 0]
% tzero     -> initial time                            [default :: 0]
% verbose   -> allow printing of info structure        [default :: True]
% window    -> select window function                  [default :: 'hann']
% zpad      -> add zero-padding to end of time-series  [default :: False]
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Version History
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/18/2022 -> Added loading bar in CalcMean(), line now wider in b2 space.
% Working on support for "BicOfTime()" i.e., time slices of bicoherence.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/14/2022 -> Added a bunch of toolbar icons as PNGs (16x16x3).  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/13/2022 -> Small tweak here and there; cross-coherence in main loop.
% Fixed error in TestSignal('cross_circle'), and fiddled with power-spec
% axes limits. Benchmarked against old code (v3.1) on DIII-D data. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/12/2022 -> Been messing with a bunch of things! Getting the point-outs
% for cross-bicoherence stuff squared away, added warning for attempts to
% plot empty data [e.g., before CalcMean() is run], tracking down a bug or
% two. [...] Fixed "v" floating around... now BicVec property. CalcMean()
% finally works with cross-b^2, and "x"s are plotted where you click!
% Finaly got to the issue with inputs to cross-wavelet-b^2 analy-sizauce.  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/11/2022 -> Tried like hell to dynamically adjust colorbar height/width
% to save some space in GUI... absolute catastrophe. Also, I realize now
% that the callbacks for the GUI figure will either have to be updated as
% the object's data is changed, or I'll need a new f^*%ing idea. [...]
% Okay, so I've calmed down and thought about it! The real question is why
% are we worried about having the callbacks be methods (static or not)? Why
% not treat it like a regular GUI and just get/set our way there? So that's
% what I've done! PlotPointOut() has been completely rewritten to either
% plot probability density [PlotType='bicoh'] or line-outs of time [else].
% Condensed the guts of many set.Methods into a single subfunction. =^] 
% Debugging time/freq-scaling issues with GUI and sitch. Added auto sigma.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/10/2022 -> Working on GUI... Real slog, might change to handle class.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/09/2022 -> Mostly messed with PyBic... Small clean-up here and there.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/08/2022 -> Made get.Samples check branchless, and added get.Raw check;
% adjusted "ProcessData" method to stop plotting, now PlotIt = true goes
% straight to GUI. Added "WhichPlot" helper method to identify the plot of 
% interest, "PlotPointOut" now reports location in bifrequency-space.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/07/2022 -> Had a good idea about making cross-bicoherence branchless.
% [Worked quite a bit with PyBic in the evening] Made general function to
% handle time/freq scaling "set" functions. Further progress on GUI. Clicks
% on bicoherence spectra allow grabbing of multiple points now; still bugs.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/06/2022 -> Working on GUI support. Want to essentially do the work of 
% GUIDE with "axes()" calls, and use ginput for point-outs time-resolved 
% biphases, etc. Learned about toolbar functions ("uitoolbar", ...), so 
% I'm thinking about implementing the functions as toolbar calls. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/05/2022 -> Added to "pybic" last night. Weekend was fun! Trying to 
% clean up what I have (i.e., a couple to-dos) before I take some vacation
% time. Should get to line-out and instantaneous freq. stuff today...
% Added support for more literate "SpecType" inputs ('fft,'wavelet',etc.),
% implemented cross-bicoherence, and messed with KeyPressFcn of PlotBispec
% figure so that SHIFT+[B,A,R,I,P] plots b2, abs(B), real/imag(B), biphase.
% Trying out support for plotting mean and std dev (if calculated); done
% with cross-spectrum, cross-coherence, etc. Want to do GUIs tomorrow. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/02/2022 -> Adjusted "PlotBispec" method to include mean and std dev
% plots (have not bug-tested this!); wedding this weekend so there might
% not be too much progress on the BicAn front.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/01/2022 -> Messed with signal generator ("SignalGen") function; wrapped
% it with new "TestSignal" method that uses a (more) convenient switchyard.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 6/30/2022 -> Lost time today due to a fiasco with my DPP poster. In any
% case, I've cleaned up a couple things, added support for auto-labeling of
% time and frequency axes, and started moving <b2> routine over. Somewhat
% worried that I haven't made a cross-bispectrum yet... Debugged the input
% parsing for multiple time-series. Mean and std dev of b2 now implemented.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 6/29/2022 -> Added wavelet stat. method and checks in set.Window. Started
% to ease into plotting support -> spectrograms implemented, also brought 
% "PlotLabels" over from pplk_bispec. Fixed FreqRes issue when nargin==1.
% Added 2 static methods for bispectral stuff: 1) for individual points,
% where random phases can be used to find PDFs, and 2) for fast production 
% of B and b2 maps from spectrograms.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 6/28/2022 -> Fixed a couple things! Moved initial (default) values to 
% "properties" block, added some access-protected props, and implemented
% static "ApplySTFT" method. Thus, "SpectroSTFT" is just a wrapper! This is
% pretty much the approach I want to take. (BicAn should be an analyzer,
% but also a nice package of signal processing tools!)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 6/27/2022 -> First code. Inspired to do O.O. b/c of Phase Space Synth!
% Honestly think that it will be glove in hand with a project like this.
% Anyway: I'm trying to write a nice, clean constructor. [...] Boom! =^O
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% STUFF TO DO!
% **Move "WinVec" to dependent properties, that way you can rid yourself of
%   the dumb double-check for window stuff.
% __Fix multi-time-series input to wavelet spec
% __Fix input "v" vector! Should be some kind of variable!
% *_Implement sizewarn
% __Line-out
% **Inst. freq/ interpolation
% __Cross-bicoherence  
% __Add CalcMean support for cross stuff (just look!)
% __Fix issue with cross-bicoh and grabbing points
% **Filter stuff!
% **Movies
% **Can colorbars be resized without weirdness?
% **Why worry about class=BicAn inputs? Isn't that the point of OOP?
% __Line-drawing issue with CWTs
% **Domain coloring for complex plots? 
% **Add time-selection if time vector given
% **Phasor diagram of biphase
% **Print default value when parsing!

% KNOWN BUGS!
% -> If data has not been created (bic.mb, say), then using PlotBispec()
%    only throws error dialog and bails. This leaves bic.PlotType as whatever
%    it was, so using a ClickPlot() call will try to produce data from an
%    empty array, despite the GUI figure displaying "good" data.
% __ FIX? Well, this is trivial in principle, but it isn't obvious what the
%    best practice is. Personally, I think WhichPlot() should throw a flag
%    to the other plotters, but it's something that will be annoying to
%    rewrite!


classdef BicAn
% Bicoherence analysis class for DSP

    % Globals 
    properties (Constant)
        WarnSize  = 1024;
        Date      = datestr(now);
    end

    % Dependents
    properties (Dependent)
        MaxRes
        Samples
        NFreq
        LineColor
    end
    
    % Private
    properties (Access=public)
        InGUI     = false;
        RunBicAn  = false;
        IsPlaying = false;
        NormToNyq = false;
        Nseries   = [];
        WinVec    = [];
        Slider
    end

    % Editables
    properties
        Note      = ' ';
        Raw       = [];
        Processed = [];
        History   = ' ';
        
        SampRate  = 1;
        FreqRes   = 0;
        SubInt    = 512;
        Step      = 128;
        Window    = 'hann';       
        Sigma     = 0;
        JustSpec  = false;
        SpecType  = 'stft';
        Bispectro = false;
        
        ErrLim    = inf;
        FScale    = 0;
        TScale    = 0;
        Filter    = 'none';
        Smooth    = 1;
        LilGuy    = 1e-6;
        SizeWarn  = true;
        BicVec    = [1 1 1];
        
        PlotIt    = true;
        CMap      = 'viridis';
        CbarNorth = true;
        PlotType  = 'bicoh';
        ScaleAxes = 'manual';
        LineWidth = 2;
        FontSize  = 20;
        PlotSlice = 0; 
        BicOfTime = false;
        
        Verbose   = false;
        Detrend   = false;
        ZPad      = false;
        Cross     = false;
        Vector    = false;
        TZero     = 0;
        
        Figure    
        AxHands   = [];
        TBHands   = [];
        
        tv = []; % Time vector
        fv = []; % Frequency vector
        ff = []; % Full frequency vector
        ft = []; % Fourier amplitudes
        sg = []; % Spectrogram (complex)
        cs = []; % Cross-spectrum
        cc = []; % Cross-coherence
        cg = []; % Coherence spectrum
        bs = []; % Bispectrum
        bc = []; % Bicoherence spectrum
        bp = []; % Biphase proxy
        bg = []; % Bispectrogram
        er = []; % Mean & std dev of FFT
        mb = []; % Mean b^2
        sb = []; % Std dev of b^2
    end % properties

    %% Class methods
    methods 

        function bic = BicAn(varargin)
        % ------------------
        % Constructor
        % ------------------  
            if nargin~=0
                bic = bic.ParseInput(varargin);
                if bic.RunBicAn
                    bic = bic.ProcessData;
                end
            end
        end % BicAn

        % ------------------  
        % "Get" functions
        % ------------------  
        function val = get.MaxRes(bic)      % Maximum resolution
            val = bic.SampRate / bic.SubInt;
        end
        function val = get.NFreq(bic)       % Number of Fourier bins
            val = floor( bic.SampRate / bic.FreqRes );
        end
        function val = get.Nseries(bic)     % Number of time-series
            [val,~] = size(bic.Raw);
        end
        function val = get.Samples(bic)     % Samples in data
            val = (~isempty(bic.Processed))*length(bic.Processed) + isempty(bic.Processed)*length(bic.Raw);
        end
        function val = get.LineColor(bic)   % For coloring plots
            val = eval(sprintf('%s(256)',bic.CMap)); 
        end

        % ------------------  
        % "Set" functions
        % ------------------  
        function bic = set.Raw(bic,val)     % Checks for raw data
            if isnumeric(val)
                dum = size(val);                % Get data dimensions
                if dum(1)>dum(2)                % Check if column vector
                    val = val.';                % Transpose
                    dum = fliplr(dum);          % Flip dimensions
                end
                if dum(1)>3
                    error('BicAn:Input','Invalid input!\nMaximum time-series is 3...\nSee "help BicAn"'); 
                end
                bic.Raw = val;
            else
                error('BicAn:Input','Invalid input!\ninData must be numeric!\nSee "help BicAn"'); 
            end
        end
        function bic = set.Window(bic,val)  % Checks for window
            bic.Window = lower(val);
            try
                window(val,10);            % Try MATLAB's windowing function
            catch winbeef 
                warning('BicAn:wrongWindow','\n"%s" window unknown... Using Hann.',bic.Window)
                bic.Window = 'hann';
            end
        end
        function bic = set.SpecType(bic,val)% Checks for spectrogram choice
            opts = {'fft ','stft ','fourier ','wave ','wavelet ','cwt '};
            bic.SpecType = CheckInString(bic.SpecType,val,opts);
            
        end
        function bic = set.PlotType(bic,val)% Checks for plotter choice
            opts = {'bicoh ','abs ','real ','imag ','angle ','mean ','std '};
            bic.PlotType = CheckInString(bic.PlotType,val,opts);
        end
        function bic = set.TScale(bic,val)  % Checks for time-scaling
            bic.TScale = CheckScales(bic.TScale,val);
        end
        function bic = set.FScale(bic,val)  % Checks for freq-scaling
            bic.FScale = CheckScales(bic.FScale,val);               
        end
        
        %% Inputs, etc.
        function bic = ParseInput(bic,vars)
        % ------------------
        % Handle inputs
        % ------------------
            Ninputs = length(vars);
            bic.RunBicAn = true;
            
            if Ninputs==1
                if isobject(vars{1}) && isequal(class(vars{1}),'BicAn')
                    % If BicAn object, output or go to GUI
                    bic = vars{1};
                    bic.RunBicAn = false;
                    try
                        bic = bic.PlotGUI;
                    catch
                    end
                elseif isnumeric(vars{1})
                    % If array input, use normalized frequencies
                    bic.Raw       = vars{1}; 
                    bic.FreqRes   = 1/bic.SubInt;     
                    bic.NormToNyq = true;
                elseif ischar(vars{1})
                    % Check string inputs
                    bic.RunBicAn = false;
                    dum = lower(vars{1});

                    %%%% Should this be global?
                    siglist = {'classic','tone','noisy','2tone','3tone','line','circle','fast_circle',...
                                'cross_2tone','cross_3tone','cross_circle','demo'};
                    switch dum
                        case 'input'
                            % Start getfile prompt
                            str = uigetfile();

                            % {PARSE FILE...}
                        case siglist
                            % If explicit test signal (or demo), confirm with user, then recursively call ParseInputs
                            str = 'Test signal!';
                            if isequal(dum,'demo')
                                dum = 'circle';  
                            end
                            dumstr = sprintf('Run the "%s" demo?',dum);

                            qwer = questdlg(dumstr,'RunDemo','Yes!','Rather not...','Yes!');
                            if isequal(qwer,'Yes!') 
                                [sig,~,fS] = bic.TestSignal(dum);
                                bic = bic.ParseInput({sig,'SampRate',fS});  
                            end
                        otherwise
                            str = 'Hmmm. That string isn`t supported yet... Try ''demo''.';
                    end
                    disp(str)
                else
                    error('BicAn:improperInput','\nInput must be BicAn object or array. "%s" class is not supported.',...
                        class(vars{1}))
                end
            else
                fprintf('Checking inputs...') 
                % Parse user-defined inputs
                bic.Raw = vars{1};
                options = fieldnames(bic); % Should I use properties(...) instead?
                %options = properties(bic);
                try
                    for i=2:2:Ninputs
                        switch lower(vars{i})            
                            case lower(options)
                                for j=1:length(options)                       
                                    if isequal(lower(options{j}),lower(vars{i}))
                                        cl = eval(sprintf('class(bic.%s)',options{j}));                   
                                        if isequal(cl,class(vars{i+1}))
                                            eval(sprintf('bic.%s = vars{i+1};',options{j}));
                                        else
                                            warning('BicAn:errOption','\n"%s" must be a %s... Using default value ().',...
                                                        vars{i},cl)
                                        end
                                    end
                                end             
                            otherwise
                                warning('BicAn:unknownOption','\n"%s" is not a known option...',vars{i})
                        end
                    end 
                catch inbeef
                    inbeef.message
                end

                % These input checks must be done in this order! Can't use set.Property functions =^\
                bic.SubInt = floor(abs(bic.SubInt));       % Remove sign and decimals
                if bic.SubInt==0 || bic.SubInt>bic.Samples % Check subinterval <= total samples
                    bic.SubInt = min(512,bic.Samples);     % Choose 512 as long as data isn't too short
                    warning('BicAn:subintWarn','Subinterval too large for time-series... Using %d.',bic.SubInt)
                end

                bic.FreqRes = abs(bic.FreqRes);            % Remove sign
                if bic.FreqRes==0                          % Check max res option
                   bic.FreqRes = bic.MaxRes;               % Maximum resolution  
                elseif bic.FreqRes<bic.MaxRes || bic.FreqRes>bic.SampRate/2
                    warning('BicAn:resWarn','Requested resolution not possible... Using maximum.')
                    bic.FreqRes = bic.MaxRes;
                end
                if bic.NFreq>bic.SubInt                    % Check if Fourier bins exceed subinterval
                    warning('BicAn:subintSmall','Subinterval too small for requested resolution... Using required.')
                    bic.FreqRes = bic.MaxRes;              % Really hate repeating code, but...
                end     

                bic.Step = floor(abs(bic.Step));           % Remove sign and decimals
                if bic.Step==0 || bic.Step>bic.SubInt      % Check step <= subinterval
                    bic.Step = floor(bic.SubInt/4);        % This seems fine?
                    warning('BicAn:stepError','Step must be nonzero and less than subint... Using %d.',bic.Step)     
                end
                fprintf('done.\n')
            end
        end % ParseInput


        function bic = ApplyFilter(bic)
        % ------------------
        % LP/BP/HP filter
        % ------------------
            %%%%%%%%%%%%%%%%
        end % ApplyFilter
        
        
        function bic = ApplyZPad(bic)
        % ------------------
        % Zero-padding
        % ------------------
            if bic.ZPad
                tail_error = mod(bic.Samples,bic.SubInt);
                if tail_error~=0
                    % Add enough zeros to make subint evenly divide samples
                    bic.Processed = [bic.Raw zeros(bic.Nseries,bic.SubInt-tail_error)];  
                else
                    bic.Processed = bic.Raw;
                end
            else
                % Truncate time series to fit integer number of stepped subintervals
                samplim = bic.Step*floor((bic.Samples - bic.SubInt)/bic.Step) + bic.SubInt;
                bic.Processed = bic.Raw(:,1:samplim);
            end
        end % ApplyZPad

       
        function bic = ProcessData(bic)
        % ------------------
        % Main processing loop
        % ------------------
            tic
            lb = LoadingBar(4,'hsv');
            
            bic = bic.ApplyZPad;
            lb = lb.TickBar(1);
            
            %bic.Processed = bic.Raw;
            switch lower(bic.SpecType)
                case {'fft','stft','fourier'}
                    bic = bic.SpectroSTFT;
                    bic.SpecType = 'stft';
                case {'wave','wavelet','cwt'}
                    
                    bic = bic.SpectroWavelet;
                    bic.SpecType = 'wave';
            end
            lb = lb.TickBar(2);
            
            if bic.Cross
                bic = bic.Coherence;
            end     
            lb = lb.TickBar(3);
            if ~bic.JustSpec
                bic = bic.Bicoherence;
            end  
            lb = lb.TickBar(4);
            fprintf('Complete! Process required %.5f seconds.\n',toc) % Done!
            
            if bic.Verbose
                disp(bic)
            end    
            lb.KillBar;

            if bic.PlotIt       
                bic = bic.PlotGUI;
            end

        end % ProcessData

        
        function [inData,t,fS] = TestSignal(bic,in)
        % ------------------
        % Provides FM test signal
        % ------------------
            fS   = 200;
            tend = 100;
            noisy = 2;
            f1 = 19;
            f2 = 45;
            switch lower(in)
                case 'classic'
                    [inData,t] = bic.SignalGen(fS,tend,1,f2,6,1,f1,10,1,1/20,noisy);
                case 'tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                case 'noisy'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,5*noisy);
                case '2tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,1,f2,0,0,0,noisy);
                case '3tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,1,f2,0,1,0,noisy);
                case 'line'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,1,f2,20,1,1/20,noisy);
                case 'circle'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,10,1,f2,10,1,1/20,noisy);
                case 'fast_circle'
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,5,1,f2,5,1,5/20,noisy);
                case 'quad_couple'
                    [x,t]  = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,f2,0,0,0,0,0,0,noisy);
                    inData = x + y + x.*y; 
                case 'coherence'
                    [x,t]  = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,f2,0,0,0,0,0,0,noisy);
                    z      = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                    inData = [x; y+z]; 
                case 'cross_2tone'
                    [x,t]  = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,f2,0,0,0,0,0,0,noisy);
                    inData = [x; x+y]; 
                case 'cross_3tone'
                    [x,t]  = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,f2,0,0,0,0,0,0,noisy);
                    z      = bic.SignalGen(fS,tend,1,f1+f2,0,0,0,0,0,0,noisy);
                    inData = [x; y; z]; 
                case 'cross_circle'
                    [x,t]  = bic.SignalGen(fS,tend,1,f1,10,0,0, 0, 0,1/20,noisy);
                    y      = bic.SignalGen(fS,tend,0,0 ,0 ,1,f2,10,0,1/20,noisy);
                    z      = bic.SignalGen(fS,tend,0,f1,10,0,f2,10,1,1/20,noisy);
                    inData = [x; y; z]; 
                otherwise
                    warning('BicAn:unknownTest','\n"%s" test signal unknown... Using single tone.',in)
                    [inData,t] = bic.SignalGen(fS,tend,1,f1,0,0,0,0,0,0,0);
            end        
        end % TestSignal   


        %% Analysis
        function bic = SpectroSTFT(bic)
        % ------------------
        % STFT method
        % ------------------     
            [spec,f,t,err,Ntoss] = bic.ApplySTFT(bic.Processed,bic.SampRate,bic.SubInt,...
                bic.Step,bic.Window,bic.NFreq,bic.TZero,bic.Detrend,bic.ErrLim,bic.Smooth);

            M = length(t);
            fprintf('%d of %d intervals (%.0f%%) exceed given tolerance\n',Ntoss,bic.Nseries*M,Ntoss/M*100);
            
            bic.tv = t;
            bic.fv = f;
            for k=1:bic.Nseries
                bic.ft(k,:) = mean(abs(spec(:,:,k)'));
            end
            bic.sg = spec;
            bic.er = err;        
        end % SpectroSTFT

        
        function bic = SpectroWavelet(bic)
        % ------------------
        % Wavelet method
        % ------------------
            if bic.Sigma==0  % Check auto
               bic.Sigma = 5*bic.Samples/bic.SampRate; 
            end

            if bic.Detrend
                for k=1:bic.Nseries
                    bic.Processed(k,:) = bic.ApplyDetrend(bic.Processed(k,:));
                end
            end
            
            for k=1:bic.Nseries
                bic.Processed(k,:) = bic.Processed(k,:) - mean(bic.Processed(k,:));
            end
            
            if bic.Samples>bic.WarnSize && bic.SizeWarn
                bic.SizeWarnPrompt(bic.Samples);
            end 

            [CWT,f,t] = bic.ApplyCWT(bic.Processed,bic.SampRate,bic.Sigma);

            bic.tv = t + bic.TZero;
            bic.fv = f;
            for k=1:bic.Nseries
                bic.ft(k,:) = mean(abs(CWT(:,:,k)'));
            end
            bic.sg = CWT;
        end % SpectroWavelet
        

        function bic = Coherence(bic)
        % ------------------
        % Cross-spectrum/coh
        % ------------------
            if bic.Nseries~=2
                warning('BicAn:crossSpec','\nCross-coherence requires exactly 2 signals!');
            else
                [cspec,crosscoh,coh] = bic.SpecToCoherence(bic.sg,bic.LilGuy);
                bic.cs = cspec;
                bic.cc = crosscoh;
                bic.cg = coh;           
            end
        end % Coherence


        function bic = Bicoherence(bic)
        % ------------------
        % Calculate bicoherence
        % ------------------       
            dum = bic.sg; 
            if isequal(bic.SpecType,'wave')
                WTrim = 50*2;
                dum = bic.sg(:,WTrim:end-WTrim,:); 
            end
            if bic.Nseries==1
                bic.BicVec = [1 1 1];
                [b2,B] = bic.SpecToBispec(dum,bic.BicVec,bic.LilGuy);
            else
                if bic.Nseries==2
                    bic.BicVec = [1 2 2];
                else
                    bic.BicVec = [1 2 3];
                end
                [b2,B] = bic.SpecToCrossBispec(dum,bic.BicVec,bic.LilGuy);
                bic.ff = [-bic.fv(end:-1:2) bic.fv];
            end

            bic.bs = B;
            bic.bc = b2;
        end % Bicoherence

        
        function bic = CalcMean(bic,Ntrials)
        % ------------------
        % Calculate mean of b^2
        % ------------------
            [n,m,r] = size(bic.sg);
            A = abs(bic.sg);
            
            if nargin==0; Ntrials = 100; end
           
            bic.mb = zeros(size(bic.bc));
            bic.sb = zeros(size(bic.bc));
            
            lb = LoadingBar(Ntrials,'hsv');
            for k=1:Ntrials
                lb = lb.TickBar(k);
                
                P = exp( 2i*pi*(2*rand(n,m,r) - 1) );

                if bic.Nseries==1
                    dum = bic.SpecToBispec(A.*P,bic.BicVec,bic.LilGuy);
                else
                    dum = bic.SpecToCrossBispec(A.*P,bic.BicVec,bic.LilGuy);
                end
                old_est = bic.mb/(k - 1 + eps);
                
                bic.mb = bic.mb + dum;
                % "Online" algorithm for variance 
                bic.sb = bic.sb + (dum - old_est).*(dum - bic.mb/k);
            end
            lb.KillBar;

            bic.mb = bic.mb/Ntrials;
            bic.sb = sqrt(bic.sb/(Ntrials-1));
        end % CalcMean
        
        
        %% Plot methods
        function PlotPowerSpec(bic)
        % ------------------
        % Plot power spectrum
        % ------------------
            if bic.PlotSlice~=0
                dum = bic.sg(:,bic.PlotSlice,:);
                [n,~,m] = size(dum);
                dum = abs(reshape(dum,n,m).').^2;
            else
                dum = (bic.ft).^2;
            end
        
            for k=1:bic.Nseries
                semilogy(bic.fv/10^bic.FScale,dum(k,:),'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:))
                if k==1; hold on; end
            end
            hold off
            
            %%% test this!
            maxspec = ceil( log10( max(abs(bic.sg(:))) ) );
            ylim([1e-5 10^maxspec]) 
            xlim([bic.fv(1)/10^bic.FScale bic.fv(end)/10^bic.FScale])
            grid on
            
            set(gca,'ytick',exp(-(5:-1:maxspec)*log(10)))
            set(gca,'xtick',linspace(0,bic.SampRate/2/10^bic.FScale,6));
            fstr = sprintf('f_n [%sHz]',bic.ScaleToString(bic.FScale));
            bic.PlotLabels({fstr,'|P|^2 [arb.]'},bic.FontSize,bic.CbarNorth);
        end % PlotPowerSpec


        function PlotSpectro(bic)
        % ------------------
        % Plot spectrograms
        % ------------------
            tstr = sprintf('Time [%ss]',bic.ScaleToString(bic.TScale));
            fstr = sprintf('f [%sHz]',bic.ScaleToString(bic.FScale));
            if isequal(bic.SpecType,'stft')
                sstr = 'log_{10}|P(t,f)|^2';
            else
                sstr = 'log_{10}|W(t,f)|^2';
            end
            for k=1:bic.Nseries
                %figure
                imagesc(bic.tv/10^bic.TScale,bic.fv/10^bic.FScale,2*log10(abs(bic.sg(:,:,k))));
                %imagesc(bic.tv/10^bic.TScale,bic.fv/10^bic.FScale,abs(bic.sg));
                bic.PlotLabels({tstr,fstr,sstr},bic.FontSize,bic.CbarNorth);
                colormap(bic.CMap);
            end
        end % PlotSpectro

                
        function PlotBispec(bic)
        % ------------------
        % Plot bispectrum
        % ------------------
        
            [~,cbarstr] = bic.WhichPlot;
            if bic.PlotSlice~=0 && bic.BicOfTime
                dumsg = bic.sg(:,bic.PlotSlice,:);
                
                if bic.Nseries==1
                    dum = bic.SpecToBispec(dumsg,bic.BicVec,bic.LilGuy);
                else
                    dum = bic.SpecToCrossBispec(dumsg,bic.BicVec,bic.LilGuy);
                end
                
            else
                [dum,~] = bic.WhichPlot;
            end
           
            %[dum,cbarstr] = bic.WhichPlot;
            if isempty(dum)
                warnstr = sprintf('"%s" data has not been created!',bic.PlotType);
                h = warndlg(warnstr,'!! Warning !!');
                uiwait(h)
                return
            end

            if bic.Nseries==1
                f = bic.fv/10^bic.FScale;
                imagesc(f,f/2,dum);
                %ylim([0 f(end)/2]); 
                line([0 f(end)/2],[0 f(end)/2],     'linewidth',2.5,'color',0.5*[1 1 1])
                line([f(end)/2 f(end)],[f(end)/2 0],'linewidth',2.5,'color',0.5*[1 1 1])
            else 
                f = bic.ff/10^bic.FScale;
                imagesc(f,f,dum) 
                line([0 f(end)],[f(end) 0],     'linewidth',2.5,'color',0.5*[1 1 1])
                line([f(1) 0],[0 f(1)],     'linewidth',2.5,'color',0.5*[1 1 1])
            end
            
            fstr1 = sprintf('f_1 [%sHz]',bic.ScaleToString(bic.FScale));
            fstr2 = sprintf('f_2 [%sHz]',bic.ScaleToString(bic.FScale));
            bic.PlotLabels({fstr1,fstr2,cbarstr},bic.FontSize,bic.CbarNorth);
            
            colormap(bic.CMap);          
        end % PlotBispec


        function [dum,cbarstr] = WhichPlot(bic)
        % ------------------
        % Helper method for plots
        % ------------------
            guy = bic.PlotType;
            switch guy
                case 'bicoh'
                    dum = bic.bc;
                    cbarstr = 'b^2(f_1,f_2)';
                case {'abs','real','imag','angle'}
                    dum = eval(sprintf('%s(bic.bs)',guy));
                    cbarstr = sprintf('%s%s B(f_1,f_2)',upper(guy(1)),guy(2:end));
                    if isequal(guy,'angle'); cbarstr = '\beta(f_1,f_2)'; end
                case 'mean'
                    dum = bic.mb;
                    cbarstr = '\langleb^2(f_1,f_2)\rangle';
                case 'std'
                    dum = bic.sb;
                    cbarstr = '\sigma_{b^2}(f_1,f_2)';
            end
        end


        function PlotConfidence(bic)
        % ------------------
        % Plot confidence interval
        % ------------------
            old_plot = bic.PlotType;
            old_dats = bic.bc;
            bic.PlotType = 'bicoh';
            noise_floor  = -bic.mb*log(1-0.999);
            bic.bc       = bic.bc .* (bic.bc>noise_floor);
            %%%%%%%%%%%%%%%%%%
            bic.bc = noise_floor;
            bic.PlotBispec;
            bic.bc       = old_dats;
            bic.PlotType = old_plot;
        end % PlotConfidence


        function PlotPointOut(bic,X,Y)
        % ------------------
        % Plot value of b^2 over time
        % ------------------
            figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fLocX = X;
            fLocY = Y;
            
            [~,ystr] = bic.WhichPlot; 
            
            dum = bic.fv/10^bic.FScale;
            if bic.Nseries>1
                dum = bic.ff/10^bic.FScale;
                X = X - (length(bic.fv)-1);
                Y = Y - (length(bic.fv)-1);
            end

            if isequal(bic.PlotType,'bicoh')
                
                Ntrials = 2000; 
                g = zeros(1,Ntrials);
                xstr = sprintf('(%3.1f,%3.1f) %sHz',dum( fLocX(1) ),dum( fLocY(1) ),bic.ScaleToString(bic.FScale));
                
                fprintf('Calculating distribution for %s...      ',xstr)  
                for k=1:Ntrials
                    LoadBar(k,Ntrials)
                    [g(k),~,~] = bic.GetBispec(bic.sg,bic.BicVec,bic.LilGuy,Y(1),X(1),true);
                end
                fprintf('\b\b^]\n')  
                
                % Limit b^2, create vector, and produce histogram 
                b2lim = 0.5;
                b2vec = linspace(0,b2lim,1000);
                cnt = hist(g,b2vec);

                % Integrate count
                intcnt = sum(cnt)*(b2vec(2)-b2vec(1));
                % exp dist -> (1/m)exp(-x/m)
                m = mean(g);
                plot(b2vec,smooth(cnt/intcnt,10),'linewidth',bic.LineWidth,'color',bic.LineColor(110,:),'marker','x','linestyle','none')
                hold on
                % More accurate distibution... Just more complicated! (Get to it later...)
                %semilogy(b2vec,(1/m)*exp(-b2vec/m).*(1-b2vec),'linewidth',bic.LineWidth,'color','red'); 
                plot(b2vec,(1/m)*exp(-b2vec/m),'linewidth',bic.LineWidth,'color','red'); 

                bic.PlotLabels({['b^2' xstr],'Probability density'},bic.FontSize,bic.CbarNorth);
                grid on
                axis tight
                legend('randomized','(1/\mu)e^{-b^2/\mu}')
                
            else

                pntstr = cell(1,length(X));
                dumt = bic.tv/10^bic.TScale;
                for k=1:length(X)
                    
                    % Calculate "point-out"
                    [~,~,Bi] = bic.GetBispec(bic.sg,bic.BicVec,bic.LilGuy,Y(k),X(k),false);
                    if isempty(Bi)
                        return
                    end
                    
                    switch bic.PlotType
                        case {'abs','imag','real'}
                            umm = eval(sprintf('%s(Bi)',bic.PlotType));
                            if isequal(bic.PlotType,'abs')
                                semilogy(dumt,umm,'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:))
                            else
                                plot(dumt,umm,'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:))
                            end
                        case 'angle'
                            plot(dumt,unwrap(angle(Bi))/pi,'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:),...
                                'linestyle','-.','marker','x')
                    end
                    if k==1; hold on; end
                    pntstr{k} = sprintf('(%3.2f,%3.2f) %sHz',dum( fLocX(k) ),dum( fLocY(k) ),bic.ScaleToString(bic.FScale));
                end            
                hold off
                xlim([dumt(1) dumt(end)])
                grid on    

                if isequal(bic.PlotType,'angle'); ystr = [ystr '/\pi']; end
                tstr = sprintf('Time [%ss]',bic.ScaleToString(bic.TScale));
                bic.PlotLabels({tstr,ystr},bic.FontSize,bic.CbarNorth);

                legend(pntstr,'fontsize',14,'fontweight','bold')                
            end
        end % PlotPointOut
        
        
        function bic = PlotGUI(bic)
        % ------------------
        % Interactive figure!
        % ------------------
            
            bic.Figure = figure('Position',[100 0 800 600],...
                'Resize','on',...
                'WindowKeyPressFcn',@SwitchPlot,...
                'WindowButtonDownFcn',@ClickPlot,...
                'DeleteFcn',@KillPlot);
            
            ax1 = axes('parent',bic.Figure,'DrawMode','fast','outerposition',[ 0   0   0.5 0.95 ]);
            ax2 = axes('parent',bic.Figure,'DrawMode','fast','outerposition',[ 0.5 0.5 0.5 0.45 ]);
            ax3 = axes('parent',bic.Figure,'DrawMode','fast','outerposition',[ 0.5 0   0.5 0.5 ]);
            
            bic.Slider = uicontrol('Style','slider',...
                'Parent',bic.Figure,...
                'Max',1,...
                'Min',0,...
                'Value',0,...
                'TooltipString','0',...
                'Units','normalized',...
                'Position',[.95 .1 .02 .2],...
                'SliderStep',[.01 .10]); %,...               
                %'Callback',@(src,event)slider_cb(pss));
                
            %%% Toolbar
            % Play/Pause
            % CalcMean()
            % Save to workspace
            % New input
            % Plot data
            % Interpolate
            
            
            % Bispectrogram
            
            
            th = uitoolbar(bic.Figure);
            % Push buttons on toolbar
            
            % Thug Life
            im = load('assets/test_pic1.mat');
            pth1a = uipushtool(th,'CData',im.M,...
                       'TooltipString','Boom!',...
                       'HandleVisibility','off');
            % Credits
            im = load('assets/test_pic3.mat');
            pth1b = uipushtool(th,'CData',im.M,...
                       'TooltipString','Credits?',...
                       'HandleVisibility','off'); 
            % CalcMean()
            im = importdata('assets/CalcMean.png');
            pth2 = uipushtool(th,'CData',im.cdata,'Separator','on',...
                       'TooltipString','Calculate <b^2>',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@CalcMeanButton);
            % SaveToWorkspace
            im = importdata('assets/SaveToWorkspace.png');
            pth3 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Save to workspace',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@CalcMeanButton);
            % NewInput
            im = importdata('assets/NewInput.png');
            pth4 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','New inputs',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@InputButton);
            % PlotData
            im = importdata('assets/PlotData.png');
            pth5 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Plot input data',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@PlotButton);
            % Interpolate
            im = importdata('assets/Interpolate.png');
            pth6 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Interpolate subinterval',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@InterpButton);
            % Bullseye
            im = importdata('assets/Bullseye.png');
            pth7 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Inspect point',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@InspectButton);
            % ConfidenceInterval
            im = importdata('assets/ConfidenceInterval.png');
            pth8 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Show noise floor',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@ConfIntButton);
            % BicOfTime
            im = importdata('assets/BicOfTime.png');
            pth9 = uitoggletool(th,'CData',im.cdata,...
                       'TooltipString','Show bispectral evolution',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@BicOfTimeButton); 
            % Play/Pause
            im = load('assets/test_pic2.mat');
            tth = uitoggletool(th,'CData',im.M,'Separator','on',...          
                                'TooltipString','Play',...
                                'HandleVisibility','off',...
                                'OnCallback',@PlayButton,...
                                'OffCallback',@PauseButton);
            % MakeMovie
            im = importdata('assets/MakeMovie1.png');
            pth10 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Output video',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@CalcMeanButton);
            % SelectInterval
            im = importdata('assets/SelectInterval2.png');
            pth11 = uipushtool(th,'CData',im.cdata,...
                       'TooltipString','Select time limits',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@CalcMeanButton);
                   
                            
            bic.AxHands = [ax1, ax2, ax3];
            bic.TBHands = [th, pth1a, pth1b, pth2, pth3, pth4, pth5, pth6,...
                            pth7, pth8, pth9, pth10, pth11 tth];
            
            bic.RefreshGUI;       
            bic.InGUI = true;
            
            set(bic.Figure,'UserData',bic);      
        end

        
        function RefreshGUI(bic)
        % ------------------
        % One stop shop for refreshment!
        % ------------------            
            set(bic.Figure,'CurrentAxes',bic.AxHands(1))
                bic.PlotBispec;
                
                if bic.PlotSlice~=0
                    m = bic.PlotSlice;
                    dum = bic.tv/10^bic.TScale;
                    tstr = sprintf('t = %6.4g %ss',dum(m),bic.ScaleToString(bic.TScale));
                    text(0.05,0.95,tstr,'units','normalized','fontsize',14,...
                                        'fontweight','bold','color','white')
                end
                
            set(bic.Figure,'CurrentAxes',bic.AxHands(2))
                
                bic.PlotSpectro;           
                if bic.PlotSlice~=0
                    dumf = bic.fv/10^bic.FScale;
                    dt   = bic.SubInt/bic.SampRate/10^bic.TScale;
                    line([dum(m),dum(m)],[0,dumf(end)],'color','white','linewidth',2)
                    line([dum(m)+dt,dum(m)+dt],[0,dumf(end)],'color','white','linewidth',2)
                end
                
            set(bic.Figure,'CurrentAxes',bic.AxHands(3))
                bic.PlotPowerSpec;    
                
            tags = {'bispec','spectro','fft'};
            for k=1:3
                set(bic.AxHands(k),'Tag',tags{k});
            end
        end % Refresh GUI
        
        
        function bic = MakeMovie(bic)
        % ------------------
        % Output movie
        % ------------------
        end % MakeMovie

        
        function bic = SizeWarnPrompt(bic,n)
        % ------------------
        % Prompt for CPU health
        % ------------------
            str = sprintf('nf = %d',n);
            qwer = questdlg({'FFT elements exceed 1000...';str;'Continue?'},...
                'FFTWarning','Absolutely!','No way...','No way...');
            pause(eps)
            if ~isequal(qwer,'Absolutely!') 
                fprintf('\nOperation terminated by user.\n')
                return         % Bail if that seems scary! 
            end
            pause(eps);
        end % SizeWarn
        
    end % methods

    %% Static methods
    methods (Static)

        function y = ApplyDetrend(y)
        % ------------------
        % Remove linear trend
        % ------------------
           fprintf('Applying detrend...') 
           n = length(y);
           s = (6/(n*(n^2-1)))*(2*sum((1:n).*y) - sum(y)*(n+1));
           % Convenient form assuming x = 1:n
           y = y - s*(1:n);
           fprintf('done.\n')
        end % ApplyDetrend
        
        
        function [b2,B] = SpecToBispec(spec,v,lilguy)
        % ------------------
        % Turns spectrogram to b^2
        % ------------------
            [nfreq,slices] = size(spec); 

            lim = nfreq;

            B  = zeros(floor(lim/2),lim);
            b2 = zeros(floor(lim/2),lim);
            
            fprintf('Calculating bicoherence...      ')     
            for j=1:floor(lim/2)+1
                LoadBar(j,floor(lim/2)+1);
                
                for k=j:lim-j+1
                    
                    p1 = spec(k,:,v(1));
                    p2 = spec(j,:,v(2));
                    s  = spec(j+k-1,:,v(3));

                    Bi  = p1.*p2.*conj(s);
                    e12 = abs(p1.*p2).^2;
                    e3  = abs(s).^2;  

                    Bjk = sum(Bi);                    
                    E12 = sum(e12);             
                    E3  = sum(e3);                      

                    b2(j,k) = (abs(Bjk).^2)./(E12.*E3+lilguy); 
                    B(j,k)  = Bjk;

                end
            end
            B = B/slices;
            fprintf('\b\b^]\n')                   
        end % SpecToBispec
        
        
        function [b2,B] = SpecToCrossBispec(spec,v,lilguy)
        % ------------------
        % Turns 2 or 3 spectrograms to b^2
        % ------------------
            [nfreq,slices] = size(spec); 
            vec = -(nfreq-1):(nfreq-1);
            lim = 2*nfreq-1;

            B  = zeros(lim);
            b2 = zeros(lim);
            
            fprintf('Calculating cross-bicoherence...      ')     
            for j=vec
                LoadBar(j+nfreq,lim);             
                for k=vec
                    if abs(j+k) < nfreq
                        p1 = real( spec(abs(k)+1,:,v(1)) ) + 1i*sign(k)*imag( spec(abs(k)+1,:,v(1)) );
                        p2 = real( spec(abs(j)+1,:,v(2)) ) + 1i*sign(j)*imag( spec(abs(j)+1,:,v(2)) );
                        s  = real( spec(abs(j+k)+1,:,v(3)) ) + 1i*sign(j+k)*imag( spec(abs(j+k)+1,:,v(3)) );

                        Bi  = p1.*p2.*conj(s);
                        e12 = abs(p1.*p2).^2;   
                        e3  = abs(s).^2;  

                        Bjk = sum(Bi);                    
                        E12 = sum(e12);             
                        E3  = sum(e3);                      

                        b2(j+nfreq,k+nfreq) = (abs(Bjk).^2)./(E12.*E3+lilguy); 

                        B(j+nfreq,k+nfreq)  = Bjk;
                    end
                end
            end
            B = B/slices;
            fprintf('\b\b^]\n')                    
        end % SpecToCrossBispec
        
        
        function [w,B,Bi] = GetBispec(spec,v,lilguy,j,k,rando)
        % ------------------
        % Calculates the bicoherence of a single (f1,f2) value
        % ------------------            
            [nfreq,slices] = size(spec);
            if abs(j+k) < nfreq
            
                p1 = real( spec(abs(k)+1,:,v(1)) ) + 1i*sign(k)*imag( spec(abs(k)+1,:,v(1)) );
                p2 = real( spec(abs(j)+1,:,v(2)) ) + 1i*sign(j)*imag( spec(abs(j)+1,:,v(2)) );
                s  = real( spec(abs(j+k)+1,:,v(3)) ) + 1i*sign(j+k)*imag( spec(abs(j+k)+1,:,v(3)) );

                if rando
                    p1 = abs(p1).*exp( 2i*pi*(2*rand(size(p1)) - 1) );
                    p2 = abs(p2).*exp( 2i*pi*(2*rand(size(p2)) - 1) );
                    s  = abs(s).* exp( 2i*pi*(2*rand(size(s)) - 1) );
                end

                Bi  = p1.*p2.*conj(s);
                e12 = abs(p1.*p2).^2;
                e3  = abs(s).^2;  

                B   = sum(Bi);                    
                E12 = sum(e12);             
                E3  = sum(e3);                      

                w = (abs(B).^2)./(E12.*E3+lilguy); 

                B = B/slices;
            else
                w  = 0;
                B  = 0;
                Bi = [];
                warning('BicAn:outofBounds','\nSelection is out of bounds... Returning null.')
            end
        end % GetBispec

        
        function [cs,cc,xx] = SpecToCoherence(spec,lilguy)
        % ------------------
        % Cross-spectrum, cross-coherence, coherogram
        % ------------------
            fprintf('Calculating cross-coherence...\n')     
            C = conj(spec(:,:,1)).*spec(:,:,2);
            N1 = mean( (abs(spec(:,:,1)).^2)' );
            N2 = mean( (abs(spec(:,:,2)).^2)' );
            
            cs = C;
            cc = abs( mean(C.') ).^2;
            cc = cc./(N1.*N2);

            xx = (abs(C).^2) ./ ( (abs(spec(:,:,1)).^2) .* (abs(spec(:,:,2)).^2) + lilguy );
        end % SpecToCoherence


        function win = HannWindow(N)
        % ------------------
        % Hann window
        % ------------------
           win = sin(pi*(0:N-1)/(N-1)).^2;
        end % HannWindow

                
        function [sig,t] = SignalGen(fS,tend,Ax,fx,Afx,Ay,fy,Afy,Az,Ff,noisy)
        % ------------------
        % Provides 3-osc FM test signal
        % ------------------
        % [sig,t] = SignalGen(fS,tend,Ax,fx,Afx,Ay,fy,Afy,Az,Ff,noisy)
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        % INPUTS:
        % fS    --> Sampling frequency in Hz
        % tend  --> End time [t = 0:1/fS:tend]
        % Ax    --> Amplitude of oscillation #1
        % fx    --> Frequency "       "      #1
        % Afx   --> Amplitude of frequency sweep
        % Ay    --> Amplitude of oscillation #2
        % fy    --> Frequency "       "      #2
        % Afy   --> Amplitude of frequency sweep
        % Az    --> Amplitude of oscillation #3
        % Ff    --> Frequency of frequency mod.
        % noisy --> Noise floor
        % OUTPUTS:
        % sig   --> Signal 
        % t     --> Time vector
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
            t = 0:1/fS:tend;  % Time-vector sampled at "fS" Hz

            dfx = Afx*sin(2*pi*t*Ff);                     % delta f1
            dfy = Afy*cos(2*pi*t*Ff);                     % delta f2
            x = Ax*sin( 2*pi*(fx*t + dfx) );              % f1
            y = Ay*sin( 2*pi*(fy*t + dfy) );              % f2
            z = Az*sin( 2*pi*(fx*t + fy*t + dfx + dfy) ); % f1 + f2

            sig = x + y + z + noisy*(0.5*rand(1,length(t)) - 1);
        end % SignalGen

             
        function cbar = PlotLabels(strings,fsize,cbarNorth)
        % ------------------
        % Convenience function
        % ------------------
            cbar = [];
            n = length(strings);
            fweight = 'bold';
            
            if n>2
                if cbarNorth
                    cbar = colorbar('location','NorthOutside');                      
                    xlabel(cbar,strings{3},'fontsize',fsize,'fontweight',fweight)
                else
                    cbar = colorbar;
                    ylabel(cbar,strings{3},'fontsize',fsize,'fontweight',fweight)
                end
            end
            
            xlabel(strings{1},'fontsize',fsize,'fontweight',fweight)
            if n>1; ylabel(strings{2},'fontsize',fsize,'fontweight',fweight); end;
            
            set(gca,'YDir','normal',...  YDir is crucial w/ @imagesc!
                'fontsize',fsize,...
                'xminortick','on',...
                'yminortick','on'); 
        end % PlotLabels

        
        function s = ScaleToString(scale)
        % ------------------
        % Frequency scaling
        % ------------------
            tags = {'n',[],[],'\mu',[],[],'m','c','d','',...
                    [],[],'k',[],[],'M',[],[],'G',[],[],'T'};   
            s = tags{10+scale};
        end % ScaleToString        

        
        function [spec,freq_vec,time_vec,err,Ntoss] = ApplySTFT(sig,samprate,subint,step,windoe,nfreq,t0,detrend,errlim,smoo)
        % ------------------
        % STFT static method
        % ------------------
            % FFT OPTIMIZATION
            fftw('planner', 'patient');
            fft(1:nfreq);
            %%%%
            [N,~] = size(sig);          % Number of time-series
            M = 1 + floor( (length(sig) - subint)/step );
            lim  = floor(nfreq/2);      % lim = |_ Nyquist/res _|
            time_vec = zeros(1,M);      % Time vector
            err  = zeros(N,M);          % Mean information
            spec = zeros(lim,M,N);      % Spectrogram
            fft_coeffs = zeros(N,lim);  % Coeffs for slice
            Ntoss = 0;                  % Number of removed slices
            
            try
                win = window(windoe,nfreq)';  % Try MATLAB's windowing function
            catch winbeef 
                warning('BicAn:wrongWindow','\n"%s" window unknown... Using Hann.',windoe)
                win = sin(pi*(0:nfreq-1)/(nfreq-1)).^2; % Apply Hann window
            end
            
            fprintf('Applying STFT...      ')
            
            for m=1:M
                LoadBar(m,M);
               
                time_vec(m) = t0 + (m-1)*step/samprate;
                for k=1:N
                    Ym = sig(k,(1:subint) + (m-1)*step); % Select subinterval     
                    Ym = Ym(1:nfreq);         % Take only what is needed for res
                    if detrend                % Remove linear least-squares fit
                        n = length(Ym);
                        s = (6/(n*(n^2-1)))*(2*sum((1:n).*Ym) - sum(Ym)*(n+1));
                        Ym = Ym - s*(1:n);
                    end
                    Ym = win.*(Ym-mean(Ym));  % Remove DC offset, multiply by window

                    DFT = fft(Ym)/nfreq;      % Limit and normalize by vector length
                    if smoo>1                 % Smooth DFT
                        DFT = smooth(DFT,floor(smoo));
                    end

                    fft_coeffs(k,:) = DFT(1:lim);            % Get interested parties
                    err(k,m) = mean(abs(fft_coeffs(k,:)));   % Find mean distribution

                    if err(k,m)>=errlim 
                        fft_coeffs(k,:) = 0*fft_coeffs(k,:); % Blank if mean excessive
                        Ntoss = Ntoss + 1; 
                    end
                    spec(:,m,k) = fft_coeffs(k,:);           % Build spectrogram

                end
            end
            fprintf('\b\b^]\n')
            
            freq_vec = (0:nfreq-1)*samprate/nfreq;
            freq_vec = freq_vec(1:lim); 
        end % ApplySTFT


        function [CWT,freq_vec,time_vec] = ApplyCWT(sig,samprate,sigma)
        % ------------------
        % Wavelet static method
        % ------------------

            [N,Nsig] = size(sig);
            nyq  = floor(Nsig/2);

            f0 = samprate/Nsig;
            freq_vec = (0:nyq-1)*f0;
            
            CWT = zeros(nyq,nyq,N);
            
            for k=1:N

                fft_sig = fft(sig(k,:));
                fft_sig = fft_sig(1:nyq);

                % Morlet wavelet in frequency space
                Psi = @(a) (pi^0.25)*sqrt(2*sigma/a) .*...
                    exp( -2 * pi^2 * sigma^2 * ( freq_vec/a - f0).^2 );

                fprintf('Applying CWT...      ')
                for a=1:nyq  
                    LoadBar(a,nyq);
                    % Apply for each scale (read: frequency)
                    CWT(a,:,k) = ifft(fft_sig.*Psi(a)); 
                end
                fprintf('\b\b^]\n')
            end

            time_vec = (0:2:Nsig-1)/samprate;
        end % ApplyCWT
        

        function bic = RunDemo
        % ------------------
        % Demonstration
        % ------------------
            bic = BicAn;
            [x,t,fS] = bic.TestSignal('circle');
            %N = 512*1;
            N = length(t);
            x = x(:,1:N);
            dT = t(N)-t(1);
            bic = BicAn(x,'sigma',dT*5,...
                            'spectype','stft',...
                            'sizewarn',false,...
                            'verbose' ,false,...
                            'samprate',fS,...
                            'justspec',~true,...
                            'plottype','bicoh');
        end % RunDemo
            

    end % static methods

end % BicAn

%% Subfunctions
function LoadBar(m,M)
% ------------------
% Help the user out!
% ------------------
    ch1 = '||\-/|||';
    ch2 = '_.:"^":.';
    fprintf('\b\b\b\b\b\b%3.0f%%%s%s',100*m/M,ch1(mod(m,8)+1),ch2(mod(m,8)+1))
end % LoadBar

function out = CheckScales(inscale,val)
    if sum(val==[-9,-6,-3,-2,-1,0,3,6,9,12])
        out = val;
    else
        warning('BicAn:scaleFail','\nIncompatable scale requested. Using ''''')
        out = inscale; 
    end   
end

function out = CheckInString(intype,val,opts)
    switch deblank(lower(val))
        case deblank(lower(opts))
            out = val;
        otherwise
            warning('BicAn:stringInput','Invalid input!\n%s does not correspond to an option.',val); 
            disp('Valid are:')
            disp(char(opts))
            out = intype;
    end
end

function SwitchPlot(obj,event)
% For changing desired plottable
% - - - - - - - - 
    bic = get(obj,'UserData');
    if ~isempty(event.Modifier)
        if isequal(event.Modifier{1},'shift')
            ind = [];
            if length(event.Character)==1
                ind = find(event.Character == 'BARIPMS',1);
            end
            if ~isempty(ind)
                figs = {'bicoh','abs','real','imag','angle','mean','std'};
                bic.PlotType = figs{ind};
            elseif isequal(event.Key,'rightarrow')
                bic.PlotSlice = mod(bic.PlotSlice, length(bic.tv)) + 1;
            elseif isequal(event.Key,'leftarrow')
                bic.PlotSlice = mod(bic.PlotSlice - 1, length(bic.tv));
            else
                return
            end
            % Activate!
            bic.RefreshGUI;
            set(obj,'UserData',bic);
        end
    end
end

function ClickPlot(obj,~)
% ------------------
% Callback for clicks
% ------------------
    bic = get(obj,'UserData');
    switch get(gca,'Tag')
        case 'spectro'
            P = get(gca,'CurrentPoint');     % Get point of mouse click  
            dum = bic.tv/10^bic.TScale;
            if P(1,1)>=dum(1) && P(1,1)<=dum(end)
                [~,m] = min(abs(dum-P(1,1))); % Find closest point in time
                bic.PlotSlice = m;
                bic.RefreshGUI;
            else
                return                       % Somewhere bad on figure
            end           
        case 'fft'
            [val,m] = max(abs(bic.sg(:,1)))  
            return
        case 'bispec' 
            button = 0;
            X = []; Y = [];
            ClickLim = 5;
            dum = bic.fv/10^bic.FScale;
            if bic.Nseries>1
                dum = bic.ff/10^bic.FScale;
                % Need to subtract something from index now!!!!
            end
            while button~=3 && length(X)<ClickLim
                [x,y,button] = ginput(1);
                if button==1 && x>=dum(1) && x<=dum(end) && y>=dum(1) && y<=dum(end)/2
                    text(x,y,'x','color',0.5*[1 1 1],'fontsize',14,'fontweight','bold')
                    [~,Ix] = min(abs(dum-x));
                    [~,Iy] = min(abs(dum-y));
                       X = [X Ix]; %#ok<*AGROW>
                       Y = [Y Iy];
                end
            end
            if ~isempty(X)
                bic.PlotPointOut(X,Y);
            end
    end
    set(obj,'UserData',bic);
end

function KillPlot(obj,~)
% ------------------
% Callback for exiting GUI
% ------------------
    bic = get(obj,'UserData');
    bic.InGUI = false;
    % Set up dialog
    options.Default = 'Yes!';
    options.WindowStyle = 'modal';
    options.Interpreter = 'TeX';
    qwer = questdlg('Save object to workspace?',...
        'GUISaveData?','Yes!','Rather not...',options);
    if isequal(qwer,'Yes!') 
        % Get variables from workspace
        lstvars = evalin('base','who');
        switch 'b'
            case lstvars
                % Output "b" for convenience
                outstr = 'b';
            otherwise
                % Otherwise do generic
                outstr = 'bicOut';
        end
        % Send out to workspace!
        assignin('base',outstr,bic);
        fprintf('\nData outputted successfully as "%s"!\n',outstr)
    end 
end

function PlayButton(~,~)
% ------------------
% Callback for play
% ------------------
    bic = get(gcbf,'UserData');
    set(gcbo,'TooltipString','Playing... Click to pause');
    
    fprintf('Note during callback = "%s"\n',bic.Note)
    %[~,M,~] = size(bic.sg);
    M = length(bic.tv);
    
    bic.IsPlaying = true;
    set(gcbf,'UserData',bic);
    
    while bic.IsPlaying
        bic = get(gcbf,'UserData');
        
        bic.PlotSlice = mod(bic.PlotSlice+10,M)+1;
        bic.RefreshGUI;
    
        set(gcbf,'UserData',bic);
        pause(eps)
    end
    
end 

function PauseButton(~,~)
% ------------------
% Callback for pause
% ------------------
    set(gcbo,'TooltipString','Play');
    bic = get(gcbf,'UserData');
    bic.IsPlaying = false;
    set(gcbf,'UserData',bic);
end

function InputButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    
    prompt = {'SubInt','Step','LilGuy','Sigma','Window',...
                'JustSpec','Detrend','ZPad','Cross','BicVec'};
    name = 'NewInput';
    for k=1:4
        def{k} = num2str(bic.(prompt{k}));
    end
    for k=6:9
        def{k} = 'false';
        if bic.(prompt{k})
            def{k} = 'true';
        end
    end
    def{5} = bic.Window;
    def{10} = sprintf('%d',bic.BicVec);
    numlines = 1;
    answer = inputdlg(prompt,name,numlines,def);
    
    if ~isempty(answer)
        
        fig = bic.Figure;
        ax  = bic.AxHands;
        
        str = ['BicAn(bic.Raw,''SampRate'',' num2str(bic.SampRate)];
        for k=[1:4 6:9]
            str = [str sprintf(',''%s'',%s',prompt{k},answer{k})];
        end
        dum = answer{10};
        bvstr = sprintf('[%d,%d,%d]',str2double(dum(1)),str2double(dum(2)),str2double(dum(3)));
        str = [str ',''Window'',''' answer{5} ''',''BicVec'',' bvstr ',''PlotIt'',false)'];

        bic = eval(str);
        
        bic.Figure  = fig;
        bic.AxHands = ax;
        bic.RefreshGUI; 
        
        set(gcbf,'UserData',bic);
    end
end

function CalcMeanButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    bic = bic.CalcMean(5);
    set(gcbf,'UserData',bic);
end

function PlotButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    figure;
    
    for k=1:bic.Nseries
        plot((0:bic.Samples-1)/10^bic.FScale/bic.SampRate,bic.Processed,'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:))
        if k==1; hold on; end
    end
    hold off
    grid on

    fstr = sprintf('t [%ss]',bic.ScaleToString(bic.TScale));
    bic.PlotLabels({fstr,'Signal [arb.]'},bic.FontSize,bic.CbarNorth);
    set(gcbf,'UserData',bic);
end

function InterpButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    
    set(gcbf,'UserData',bic);
end

function BicOfTimeButton(~,~)
% ------------------
% Callback for BicOfTime()
% ------------------
    bic = get(gcbf,'UserData');
    bic.BicOfTime = ~bic.BicOfTime;
    set(gcbf,'UserData',bic);
end

function InspectButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    
    set(gcbf,'UserData',bic);
end

function ConfIntButton(~,~)
% ------------------
% Callback for CalcMean()
% ------------------
    bic = get(gcbf,'UserData');
    
    set(gcbf,'UserData',bic);
end
