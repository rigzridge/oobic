%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%/\\\\\\\\\\\\\___________________________/\\\\\\\\\__________________        
%\/\\\/////////\\\_______________________/\\\\\\\\\\\\\________________       
%_\/\\\_______\/\\\__/\\\________________/\\\/////////\\\_______________      
% _\/\\\\\\\\\\\\\\__\///______/\\\\\\\\_\/\\\_______\/\\\__/\\/\\\\\\___     
%  _\/\\\/////////\\\__/\\\___/\\\//////__\/\\\\\\\\\\\\\\\_\/\\\////\\\__    
%   _\/\\\_______\/\\\_\/\\\__/\\\_________\/\\\/////////\\\_\/\\\__\//\\\_   
%    _\/\\\_______\/\\\_\/\\\_\//\\\________\/\\\_______\/\\\_\/\\\___\/\\\_  
%     _\/\\\\\\\\\\\\\/__\/\\\__\///\\\\\\\\_\/\\\_______\/\\\_\/\\\___\/\\\ 
%      _\/////////////____\///_____\////////__\///________\///__\///____\///
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Object-oriented bicoherence analysis and signal processing toolbox
% v4.0 (c) 2022 G.A. Riggs, WVU Dept. of Physics & Astronomy 
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
% autoscale -> autoscaling in figures [default :: false]
% bispectro -> computes bispectrogram [default :: false]
% cbarnorth -> control bolorbar location [default :: true]
% cmap      -> adjust colormap [default :: viridis]
% dealias   -> applies antialiasing (LP) filter [default :: false]
% detrend   -> remove linear trend from data [default :: false]
% errlim    -> mean(fft) condition [default :: inf] 
% filter    -> xxxxxxxxxxxxxxx [default :: 'none']
% freqres   -> desired frequency resolution [Hz]
% fscale    -> scale for plotting frequencies [default :: 0]
% justspec  -> true for just spectrogram [default :: false]
% lilguy    -> set epsilon [default :: 1e-6]
% note      -> optional string for documentation [default :: ' '] 
% plotit    -> start plotting tool when done [default :: false]
% plottype  -> set desired plottable [default :: 'bicoh']
% samprate  -> sampling rate in Hz [default :: 1]
% sigma     -> parameter for wavelet spectrum [default :: 1]
% spectype  -> set desired time-freq. method [default :: 'stft']
% step      -> step size for Welch method in samples [default :: 512]
% subint    -> subinterval size in samples [default :: 128]
% sizewarn  -> warning for matrix size [default :: true]
% smooth    -> smooths FFT by n samples [default :: 1]
% tscale    -> scale for plotting time [default :: 0]
% verbose   -> allow printing of info structure [default :: true]
% window    -> select window function [default :: 'hann' (see @window)]
% zpad      -> add zero-padding to end of time-series [default :: true]
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Version History
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 7/11/2022 -> Tried lke hell to dynamically adjust colorbar's height/width
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
% 7/10/2022 -> Working on GUI stuff... Real slog
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
% on bicoherence spectra allow grabbing of multiple points now, ...
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
% Anyway: I'm trying to write a nice, clean constructor. [...]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% STUFF TO DO!
% **Move "WinVec" to dependent properties, that way you can rid yourself of
%   the dumb double-check for window stuff.
% **Fix multi-time-series input to wavelet spec
% *_Fix input "v" vector! Should be some kind of variable!
% *_Implement sizewarn
% *_Line-out
% **Inst. freq/ interpolation
% *_Cross-bicoherence
% **set(gcf,'WindowKeyPressFcn',@(src,event)SwitchPlot(bic,event)) is
%   necessary to keep keypresses up-to-date with object...  
% **Add CalcMean support for cross stuff (just look!)
% **Fix issue with cross-bicoh and grabbing points


% **WTF IS HAPPENING WITH INGUI NOT SETTING?


classdef BicAn
% Bicoherence analysis class for DSP

    % Globals 
    properties (Constant)
        LineWidth = 2;
        FontSize  = 20;
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
    properties (Access=private)
        InGUI     = false;
        RunBicAn  = false;
        IsPlaying = false;
        NormToNyq = false;
        Nseries   = [];
        WinVec    = [];
        Figure
        Axes
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
        ErrLim    = inf;
        FScale    = 0;
        TScale    = 0;
        Filter    = 'none';
        Bispectro = false;
        Smooth    = 1;
        PlotIt    = true;
        LilGuy    = 1e-6;
        SizeWarn  = true;
        CMap      = 'viridis';
        CbarNorth = true;
        PlotType  = 'bicoh';
        ScaleAxes = 'manual';
        Verbose   = false;
        Detrend   = false;
        ZPad      = false;
        Cross     = false;
        Vector    = false;
        TZero     = 0;
        PlotSlice = 0; 
        
        TBHands   = [];
        
        tv = []; % Time vector
        fv = []; % Frequency vector
        ff = []; % Full frequency vector
        ft = []; % Fourier amplitudes
        sg = []; % Spectrogram (complex)
        xs = []; % Cross-spectrum
        xc = []; % Cross-coherence
        cs = []; % Coherence spectrum
        bs = []; % Bispectrum
        bc = []; % Bicoherence spectrum
        bp = []; % Biphase proxy
        bg = []; % Bispectrogram
        er = []; % Mean & std dev of FFT
        mb = []; % Mean b^2
        sb = []; % Std dev of b^2
    end % properties

    % Functions
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
        
        
        function bic = ParseInput(bic,vars)
        % ------------------
        % Handle inputs
        % ------------------
            Ninputs = length(vars);
            bic.RunBicAn = true;
            
            if Ninputs==1
                if isobject(vars{1})
                    % If object, output or go to GUI
                    bic = vars{1};
                    bic.RunBicAn = false;
                elseif isnumeric(vars{1})
                    % If array input, use normalized frequencies
                    bic.Raw       = vars{1};
                    bic.SampRate  = 1;  
                    bic.FreqRes   = 1/bic.SubInt;     
                    bic.NormToNyq = true;
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
                                            warning('BicAn:errOption','\n"%s" must be a %s... Using default value.',...
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

       
        function bic = ProcessData(bic)
        % ------------------
        % Main processing loop
        % ------------------
            tic
            bic = bic.ApplyZPad;
            %bic.Processed = bic.Raw;
            switch lower(bic.SpecType)
                case {'fft','stft','fourier'}
                    bic = bic.SpectroSTFT;
                    bic.SpecType = 'stft';
                case {'wave','wavelet','cwt'}
                    if bic.Sigma==0
                       bic.Sigma = 5*bic.Samples/bic.SampRate; 
                    end
                    bic = bic.SpectroWavelet;
                    bic.SpecType = 'wave';
            end        
            if ~bic.JustSpec
                bic = bic.Bicoherence;
            end            
            toc
            
            if bic.Verbose
                disp(bic)
            end       

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
            switch lower(in)
                case 'classic'
                    [inData,t] = bic.SignalGen(fS,tend,1,45,6,1,22,10,1,1/20,noisy);
                case 'tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,noisy);
                case 'noisy'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,5*noisy);
                case '2tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,1,45,0,0,0,noisy);
                case '3tone'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,1,45,0,1,0,noisy);
                case 'line'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,1,45,10,1,1/20,noisy);
                case 'circle'
                    [inData,t] = bic.SignalGen(fS,tend,1,22,10,1,45,10,1,1/20,noisy);
                case 'cross_2tone'
                    [x,t]  = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,45,0,0,0,0,0,0,noisy);
                    inData = [x; x+y]; 
                case 'cross_3tone'
                    [x,t]  = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,noisy);
                    y      = bic.SignalGen(fS,tend,1,45,0,0,0,0,0,0,noisy);
                    z      = bic.SignalGen(fS,tend,1,67,0,0,0,0,0,0,noisy);
                    inData = [x; y; z]; 
                case 'cross_circle'
                    [x,t]  = bic.SignalGen(fS,tend,1,22,10,0,0,0,0,1/20,noisy);
                    y      = bic.SignalGen(fS,tend,1,45,10,0,0,0,0,1/20,noisy);
                    z      = bic.SignalGen(fS,tend,0,22,10,0,45,10,1,1/20,noisy);
                    inData = [x; y; z]; 
                otherwise
                    warning('BicAn:unknownTest','\n"%s" test signal unknown... Using single tone.',in)
                    [inData,t] = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,0);
            end        
        end % TestSignal   


        function bic = ApplyFilter(bic)
        % ------------------
        % LP/BP/HP filter
        % ------------------
            %%%%%%%%%%%%%%%%
        end % ApplyFilter

            
        function bic = SpectroSTFT(bic)
        % ------------------
        % STFT method
        % ------------------     
            [spec,f,t,err,Ntoss] = bic.ApplySTFT(bic.Processed,bic.SampRate,bic.SubInt,...
                bic.Step,bic.Window,bic.NFreq,bic.TZero,bic.Detrend,bic.ErrLim,bic.Smooth);
            
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
            if bic.Detrend
                bic.Processed = bic.ApplyDetrend(bic.Processed);
            end
            bic.Processed = bic.Processed - mean(bic.Processed);
            
            if length(bic.Processed)>bic.WarnSize && bic.SizeWarn
                bic.SizeWarnPrompt(length(bic.Processed));
            end 

            %[CWT,f,t] = bic.ApplyCWT(bic.Processed,bic.SampRate,bic.Sigma);
            for k=1:bic.Nseries
                [CWT(:,:,k),f,t] = bic.ApplyCWT(bic.Processed,bic.SampRate,bic.Sigma);
            end

            bic.tv = t + bic.TZero;
            bic.fv = f;
            bic.ft = mean(abs(CWT'));    
            bic.sg = CWT;
        end % SpectroWavelet
        

        function bic = Coherence(bic)
        % ------------------
        % Cross-spectrum/coh
        % ------------------
            if bic.Nseries~=2
                error('BicAn:crossSpec','\nCross-coherence requires exactly 2 signals!');
            else
                [cspec,crosscoh,coh] = bic.SpecToCoherence(bic.sg,bic.LilGuy);
                bic.cs = cspec;
                bic.xc = crosscoh;
                bic.xs = coh;           
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
                v = [1 1 1];
                [b2,B] = bic.SpecToBispec(dum,v,bic.LilGuy);
            else
                if bic.Nseries==2
                    v = [1 2 2];
                else
                    v = [1 2 3];
                end
                [b2,B] = bic.SpecToCrossBispec(dum,v,bic.LilGuy);
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
            v = [1 1 1];
            A = abs(bic.sg);
            
            if nargin==0; Ntrials = 100; end
            bic.mb = zeros(floor(n/2),n);

            bic.sb = bic.mb;
            for k=1:Ntrials

                P = exp( 2i*pi*(2*rand(n,m,r) - 1) );

                dum = bic.SpecToBispec(A.*P,v,bic.LilGuy);
                old_est = bic.mb/(k - 1 + eps);
                
                bic.mb = bic.mb + dum;
                % "Online" algorithm for variance 
                bic.sb = bic.sb + (dum - old_est).*(dum - bic.mb/k);
            end

            bic.mb = bic.mb/Ntrials;
            bic.sb = sqrt(bic.sb/(Ntrials-1));
        end % CalcMean
        
        
        function PlotPowerSpec(bic)
        % ------------------
        % Plot power spectrum
        % ------------------
            if bic.PlotSlice~=0
                dum = abs(bic.sg(:,bic.PlotSlice,:).').^2;
            else
                dum = bic.ft;
            end
        
            for k=1:bic.Nseries
                semilogy(bic.fv/10^bic.FScale,dum(k,:),'linewidth',bic.LineWidth,'color',bic.LineColor(50+40*k,:))
                if k==1; hold on; end
            end
            hold off
            %ylim([1e-10 1e1]) %%%%%%%%%
            xlim([bic.fv(1)/10^bic.FScale bic.fv(end)/10^bic.FScale])
            grid on
            
            %set(gca,'xticklabels',linspace(bic.fv(1),bic.fv(end)+bic.fv(1),6));
            fstr = sprintf('f_n [%sHz]',bic.ScaleToString(bic.FScale));
            bic.PlotLabels({fstr,'|P| [arb.]'},bic.FontSize,bic.CbarNorth);
        end % PlotPowerSpec


        function PlotSpectro(bic)
        % ------------------
        % Plot spectrograms
        % ------------------
            tstr = sprintf('Time [%ss]',bic.ScaleToString(bic.TScale));
            fstr = sprintf('f [%sHz]',bic.ScaleToString(bic.FScale));
            if isequal(bic.SpecType,'stft')
                sstr = 'log_{10}|P(t,f)|';
            else
                sstr = 'log_{10}|W(t,f)|';
            end
            for k=1:bic.Nseries
                %figure
                imagesc(bic.tv/10^bic.TScale,bic.fv/10^bic.FScale,log10(abs(bic.sg(:,:,k))));
                bic.PlotLabels({tstr,fstr,sstr},bic.FontSize,bic.CbarNorth);
                colormap(bic.CMap);
            end
        end % PlotSpectro

                
        function PlotBispec(bic)
        % ------------------
        % Plot bispectrum
        % ------------------
           
            [dum,cbarstr] = bic.WhichPlot;

            if bic.Nseries==1
                f = bic.fv/10^bic.FScale;
                h = imagesc(f,f/2,dum);
                %ylim([0 f(end)/2]); 
                line([0 f(end)/2],[0 f(end)/2],'linewidth',2.5,'color',0.5*[1 1 1])
                line([f(end)/2 f(end)],[f(end)/2 0],'linewidth',2.5,'color',0.5*[1 1 1])
            else 
                f = bic.ff/10^bic.FScale;
                imagesc(f,f,dum) 
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
            
            [~,ystr] = bic.WhichPlot; 
            
            v = [1 1 1];
            dum = bic.fv/10^bic.FScale;
            if bic.Nseries>1
                dum = bic.ff/10^bic.FScale;
                crossplot = true;
            end
            
            if isequal(bic.PlotType,'bicoh')
                
                Ntrials = 1000; 
                g = zeros(1,Ntrials);
                
                for k=1:Ntrials
                    [g(k),~,~] = bic.GetBispec(bic.sg,v,bic.LilGuy,Y(1),X(1),true);
                end
                
                % Limit b^2, create vector, and produce histogram 
                b2lim = 0.2;
                b2vec = linspace(0,b2lim,1000);
                cnt = hist(g,b2vec);

                % Integrate count
                intcnt = sum(cnt)*(b2vec(2)-b2vec(1));
                % exp dist -> (1/m)exp(-x/m)
                m = mean(g);
                plot(b2vec,cnt/intcnt,'linewidth',bic.LineWidth,'color',bic.LineColor(110,:))
                hold on
                plot(b2vec,(1/m)*exp(-b2vec/m).*(1-b2vec),'linewidth',bic.LineWidth,'color','red'); 

                xstr = sprintf('(%3.1f,%3.1f) %sHz',dum(X(1)),dum(Y(1)),bic.ScaleToString(bic.FScale));
                bic.PlotLabels({['b^2' xstr],'Probability density'},bic.FontSize,bic.CbarNorth);
                grid on
                legend('randomized','(1/\mu)e^{-b^2/\mu}')
                
            else

                pntstr = cell(1,length(X));
                dumt = bic.tv/10^bic.FScale;
                for k=1:length(X)
                    
                    % Calculate "point-out"
                    [~,~,Bi] = bic.GetBispec(bic.sg,v,bic.LilGuy,Y(k),X(k),false);
                    
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
                    pntstr{k} = sprintf('(%3.2f,%3.2f) %sHz',dum(X(k)),dum(Y(k)),bic.ScaleToString(bic.FScale));
                end            
                hold off
                xlim([dumt(1) dumt(end)])
                grid on
                

                tstr = sprintf('Time [%ss]',bic.ScaleToString(bic.TScale));
                bic.PlotLabels({tstr,ystr},bic.FontSize,bic.CbarNorth);

                legend(pntstr,'fontsize',14,'fontweight','bold')                
            end
        end % PlotPointOut
        
        
        function bic = PlotGUI(bic)
            
            h = figure('Position',[100 0 800 600],...
                'Resize','on',...
                'WindowKeyPressFcn',@SwitchPlot,...
                'WindowButtonDownFcn',@ClickBispec);
            
            bic.Axes = axes('DataAspectRatio',[1 1 1],...
                'XLimMode','manual','YLimMode','manual',...
                'DrawMode','fast',...
                'Parent',h);
            bic.Slider = uicontrol('Style','slider',...
                'Parent',h,...
                'Max',1,...
                'Min',0,...
                'Value',0,...
                'TooltipString','0',...
                'Units','normalized',...
                'Position',[.95 .1 .02 .2],...
                'SliderStep',[.01 .10]); %,...               
                %'Callback',@(src,event)slider_cb(pss));
            
            th = uitoolbar(h);
            % Push buttons on toolbar
            img1 = rand(16,16,3);
            pth1 = uipushtool(th,'CData',img1,...
                       'TooltipString','My push tool',...
                       'HandleVisibility','off');
                   
            img2 = rand(16,16,3);
            pth2 = uipushtool(th,'CData',img2,...
                       'TooltipString','CalcMean()',...
                       'HandleVisibility','off',...
                       'ClickedCallback',@CalcMeanButton);
                   
            % Add a toggle tool to the toolbar
            img3 = rand(16,16,3);
            tth = uitoggletool(th,'CData',img3,'Separator','on',...          
                                'TooltipString','Play',...
                                'HandleVisibility','off',...
                                'OnCallback',@PlayButton,...
                                'OffCallback',@PauseButton);
                            
            bic.TBHands = [th, pth1, pth2, tth];
            
            bic.RefreshGUI;       
            bic.InGUI = true;
            
            set(h,'UserData',bic);
            
        end

        function RefreshGUI(bic)
        % ------------------
        % Callback for clicks
        % ------------------
            %%%%%%%%   clf
            ax(1) = subplot(2,2,[1 3]);
                bic.PlotBispec;
                set(ax(1),'outerposition',[ 0 0 0.5 1 ]);
                
                if bic.PlotSlice~=0
                    m = bic.PlotSlice;
                    dum = bic.tv/10^bic.TScale;
                    tstr = sprintf('t = %6.4g %ss',dum(m),bic.ScaleToString(bic.TScale));
                    text(0.05,0.95,tstr,'units','normalized','fontsize',14,...
                                        'fontweight','bold','color','white')
                end
                
            ax(2) = subplot(2,2,2);
                
                bic.PlotSpectro;           
                if bic.PlotSlice~=0
                    %m = bic.PlotSlice;
                    dt = bic.SubInt/bic.SampRate/10^bic.TScale;
                    line([dum(m),dum(m)],[0,dum(end)],'color','white','linewidth',2)
                    line([dum(m)+dt,dum(m)+dt],[0,dum(end)],'color','white','linewidth',2)
                end
                set(ax(2),'outerposition',[ 0.5 0.5 0.5 0.5 ]);
                
            ax(3) = subplot(2,2,4);
                bic.PlotPowerSpec;    
                set(ax(3),'outerposition',[ 0.5 0 0.5 0.5 ]);
                
            tags = {'bispec','spectro','fft'};
            for k=1:3
                set(ax(k),'Tag',tags{k});
            end
        end % Refresh GUI
        
        
        function bic = Dum(bic,event)
        % ------------------
        % Callback for clicks
        % ------------------
            fprintf('Note during callback = "%s"\n',bic.Note)
            for m=1:2:100
                bic.PlotSlice = m;
                bic.RefreshGUI;
                pause(eps)
            end
        end 
        
        function bic = MakeMovie(bic)
        % ------------------
        % Output movie
        % ------------------

        end % MakeMovie

       
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

        
        function bic = SizeWarnPrompt(bic,n)
        % ------------------
        % Prompt for 
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

            %p1 = spec(k,:,v(1));
            %p2 = spec(j,:,v(2));
            %s  = spec(j+k-1,:,v(3));

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
            
            B = B/length(Bi);
        end % GetBispec

        
        function [cs,cc,xx] = SpecToCoherence(spec,lilguy)
        % ------------------
        % Cross-spectrum, cross-coherence, coherogram
        % ------------------
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
                    xlabel(cbar,strings{3},'fontsize',fsize,'fontweight','bold')
                else
                    cbar = colorbar;
                    ylabel(cbar,strings{3},'fontsize',fsize,'fontweight','bold')
                end
            end
            
            xlabel(strings{1},'fontsize',fsize,'fontweight','bold')
            if n>1; ylabel(strings{2},'fontsize',fsize,'fontweight','bold'); end;
            
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
                warning('BicAn:wrongWindow','\n"%s" window unknown... Using Hanning.',windoe)
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

            Nsig = length(sig);
            nyq  = floor(Nsig/2);

            f0 = samprate/Nsig;
            freq_vec = (0:nyq-1)*f0;
            
            CWT = zeros(nyq);

            fft_sig = fft(sig);
            fft_sig = fft_sig(1:nyq);

            % Morlet wavelet in frequency space
            Psi = @(a) (pi^0.25)*sqrt(2*sigma/a) .*...
                exp( -2 * pi^2 * sigma^2 * ( freq_vec/a - f0).^2 );

            fprintf('Applying CWT...      ')
            for a=1:nyq  
                LoadBar(a,nyq);
                % Apply for each scale (read: frequency)
                CWT(a,:) = ifft(fft_sig.*Psi(a)); 
            end
            fprintf('\b\b^]\n')

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
        warning('BicAn:scaleFail','\nIncompatable time scale requested.')
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
    press = event.Character;
    switch press
        case {'B','A','R','I','P','M','S'}
            if isequal(press,'B') 
                bic.PlotType = 'bicoh';
            elseif isequal(press,'A') 
                bic.PlotType = 'abs';
            elseif isequal(press,'R') 
                bic.PlotType = 'real';
            elseif isequal(press,'I') 
                bic.PlotType = 'imag';
            elseif isequal(press,'P') 
                bic.PlotType = 'angle';
            elseif isequal(press,'M') 
                bic.PlotType = 'mean';
            elseif isequal(press,'S') 
                bic.PlotType = 'std';
            end
            
            % Activate!
            bic.RefreshGUI;
            set(obj,'UserData',bic);
            
        otherwise
            return
    end
end


function ClickBispec(obj,~)
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
                if button==1 && x>=dum(1) && x<=dum(end) && y>=dum(1) && y<=dum(end)
                    [~,Ix] = min(abs(dum-x));
                    [~,Iy] = min(abs(dum-y));
                       X = [X Ix];
                       Y = [Y Iy];
                end
            end
            if ~isempty(X)
                bic.PlotPointOut(X,Y);
            end
    end
    set(obj,'UserData',bic);
end

function PlayButton(~,~)
% ------------------
% Callback for play
% ------------------
    bic = get(gcbf,'UserData');
    set(gcbo,'TooltipString','Playing... Click to pause');
    
    fprintf('Note during callback = "%s"\n',bic.Note)
    bic.CMap = 'cmr';
    [~,M,~] = size(bic.sg);
    
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
% Callback for play
% ------------------
    bic = get(gcbf,'UserData');
    bic.IsPlaying = false;
    set(gcbf,'UserData',bic);
end


function CalcMeanButton(~,~)
% ------------------
% Callback for play
% ------------------
    bic = get(gcbf,'UserData');
    bic = bic.CalcMean(5);
    set(gcbf,'UserData',bic);
end



