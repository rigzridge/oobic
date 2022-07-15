%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%             Object-oriented bicoherence analyzer (OOBic) 
%               o\/\/\/\/\/\/\/\/ .__. \/\/\/\/\/\/\/\/\o
%          by G. Riggs for Plasma Physics Laboratory - Koepke
%                 - - - - - - - - - - - - - - - - - - -
%                   (c) 2022 West Virginia University
%                 - - - - - - - - - - - - - - - - - - -
%                   contact me :: gariggs@mix.wvu.edu
%               o\/\/\/\/\/\/\/\/ ^--^ \/\/\/\/\/\/\/\/\o
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% The Bispectrum
% B_xyz(f1,f2) = < X(f1)Y(f2)Z(f1+f2)* >, where x,y,z are time series with 
% corresponding Fourier transforms X,Y,Z, and <...> denotes averaging.
% - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% The (squared) Bicoherence spectrum
% b^2_xyz(f1,f2) =           |B_xyz(f1,f2)|^2
%                          --------------------
%                ( <|X(f1)Y(f2)|^2> <|Z(f1+f2)|^2> + eps ),
% where eps is a small number meant to prevent 0/0 = NaN catastrophe
% - - - - - - - - - - - - - - - - - - - -
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Inputs
% - - - - - - - - - - - - - - - - - - - - 
% inData    -> time-series {or structure}
% - - - - - - - - - - - - - - - - - - - - 
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
% note      -> optional string for documentation [default :: {DATE & TIME}] 
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


% STUFF TO DO!
% **Move "WinVec" to dependent properties, that way you can rid yourself of
%   the dumb double-check for window stuff.
% **Try to make anonymous function for t/fscale set methods
% **Fix multi-time-series input to wavelet spec
% *_Fix input "v" vector! Should be some kind of variable!
% **Implement sizewarn
% **Line-out
% **Inst. freq/ interpolation
% *_Cross-bicoherence
% **set(gcf,'WindowKeyPressFcn',@(src,event)SwitchPlot(bic,event)) is
%   necessary to keep keypresses up-to-date with object...  
% **Add CalcMean support for cross stuff (just look!)


classdef BicAn
% Bicoherence analysis class for DSP

    % Globals 
    properties (Constant)
        FontSize = 20;
        WarnSize = 1000;
        Date     = datestr(now);
    end

    % Dependents
    properties (Dependent)
        MaxRes
        Samples
        NFreq
    end
    
    % Private
    properties (Access=private)
        RunBicAn  = false;
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
        Sigma     = 1;
        JustSpec  = false;
        SpecType  = 'stft';
        ErrLim    = inf;
        FScale    = 0;
        TScale    = 0;
        Filter    = 'none';
        Bispectro = false;
        Smooth    = 1;
        PlotIt    = false;
        LilGuy    = 1e-6;
        SizeWarn  = true;
        CMap      = 'viridis';
        CbarNorth = true;
        PlotType  = 'bicoh';
        ScaleAxes = 'manual';
        Verbose   = true;
        Detrend   = false;
        ZPad      = false;
        Cross     = false;
        Vector    = false;
        TZero     = 0;
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
        function val = get.MaxRes(bic)
            val = bic.SampRate / bic.SubInt;
        end
        function val = get.NFreq(bic)
            val = floor( bic.SampRate / bic.FreqRes );
        end
        function val = get.Nseries(bic)
            [val,~] = size(bic.Raw);
        end
        function val = get.Samples(bic)
            if ~isempty(bic.Processed)
                val = length(bic.Processed);
            else 
                val = length(bic.Raw);
            end
        end

        % ------------------  
        % "Set" functions
        % ------------------  
        function bic = set.Window(bic,val)
            bic.Window = lower(val);
            try
                window(val,10);  % Try MATLAB's windowing function
            catch winbeef 
                warning('BicAn:wrongWindow','\n"%s" window unknown... Using Hann.',bic.Window)
                bic.Window = 'hann';
            end
        end
        function bic = set.SpecType(bic,val)
            switch lower(val)
                case {'fft','stft','fourier','wave','wavelet','cwt'}
                    bic.SpecType = val;
                otherwise
                    warning('BicAn:specInput','Invalid input!\n%s does not correspond to a type of spectrogram.\nSee "help BicAn"',val); 
                    bic.SpecType = 'stft';
            end
        end
        function bic = set.TScale(bic,val)
            if sum(val==[-9,-6,-3,-2,-1,0,3,6,9,12])
                bic.TScale = val;
            else
                warning('BicAn:scaleFail','\nIncompatable time scale requested.')
            end                
        end
        function bic = set.FScale(bic,val)
            if sum(val==[-9,-6,-3,-2,-1,0,3,6,9,12])
                bic.FScale = val;
            else
                warning('BicAn:scaleFail','\nIncompatable frequency scale requested.')
            end                
        end

        function bic = ParseInput(bic,vars)
        % ------------------
        % Handle inputs
        % ------------------
            Ninputs = length(vars);
            bic.RunBicAn = true;
            
            if isnumeric(vars{1})
                [nrows,ncols] = size(vars{1});   % Get data dimensions
                if nrows>ncols                   % Check if column vector
                    vars{1} = vars{1}.';         % Transpose
                end
                if bic.Nseries>3
                    error('BicAn:Input','Invalid input!\nMaximum time-series is 3...\nSee "help BicAn"'); 
                end
            end
            
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
                    inbeef
                end

                % Do checks on inputs
                bic.SubInt = floor(abs(bic.SubInt));       % Remove sign and decimals
                if bic.SubInt==0 || bic.SubInt>bic.Samples % Check subinterval <= total samples
                    bic.SubInt = min(512,bic.Samples);     % Choose 512 as long as data isn't short
                    warning('BicAn:subintWarn','Subinterval too large for time-series... Using %d.',bic.SubInt)
                end

                bic.FreqRes = floor(abs(bic.FreqRes));     % Remove sign and decimals
                if bic.FreqRes==0                          % Check max res option
                   bic.FreqRes = bic.MaxRes;               % Maximum resolution  
                elseif bic.FreqRes<bic.MaxRes || bic.FreqRes>bic.SampRate/2
                    warning('Bispec:resWarn','Requested resolution not possible... Using maximum.')
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
                    bic = bic.SpectroWavelet;
                    bic.SpecType = 'wave';
            end        
            if ~bic.JustSpec
                bic = bic.Bicoherence;
                %bic = bic.CalcMean(15);
                h = figure('WindowKeyPressFcn',@(src,event)SwitchPlot(bic,event));
                bic.PlotBispec;
            end            
            toc
            
            bic.PlotSpectro;
            
            if bic.Verbose
                disp(bic)
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
                    [x,t] = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,noisy);
                    y     = bic.SignalGen(fS,tend,1,45,0,0,0,0,0,0,noisy);
                    inData = [x; y]; 
                case 'cross_3tone'
                    [x,t] = bic.SignalGen(fS,tend,1,22,0,0,0,0,0,0,noisy);
                    y     = bic.SignalGen(fS,tend,1,45,0,0,0,0,0,0,noisy);
                    z     = bic.SignalGen(fS,tend,1,67,0,0,0,0,0,0,noisy);
                    inData = [x; y; z]; 
                case 'cross_circle'
                    [x,t] = bic.SignalGen(fS,tend,1,22,10,0,0,0,0,1/20,noisy);
                    y     = bic.SignalGen(fS,tend,1,45,10,0,0,0,0,1/20,noisy);
                    z     = x.*y;
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

            [CWT,f,t] = bic.ApplyCWT(bic.Processed,bic.SampRate,bic.Sigma);
            %for k=1:bic.Nseries
            %    [CWT(:,:,k),f,t] = bic.ApplyCWT(bic.Processed,bic.SampRate,bic.Sigma);
            %end

            bic.tv = t;
            bic.fv = f;
            bic.ft = mean(abs(CWT'));    
            bic.sg = CWT;
        end % SpectroWavelet
        
        function [b2,B] = SpecToBispec_slow(bic,spec,v,lilguy)
        % ------------------
        % Turns spectrogram to b^2
        % ------------------
            limx = floor(bic.NFreq/2);
            limy = floor(limx/2)+1;
            
            B    = zeros(limy,limx);
            b2   = zeros(limy,limx);
            
            fprintf('Working...      ')     
            for j=1:limy
                LoadBar(j,limy);
                for k=j:limx-j+1
                    [b2jk,Bjk] = bic.GetBispec(spec,v,lilguy,j,k,false);
                    B(j,k)  = Bjk;
                    b2(j,k) = b2jk;
                end
            end
            fprintf('\b\b^]\n')                   
        end % SpecToBispec
          
        function bic = PlotSpectro(bic)
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
                
        function bic = PlotBispec(bic)
        % ------------------
        % Plot bispectrum
        % ------------------
           
            guy = bic.PlotType;
            switch guy
                case 'bicoh'
                    dum = bic.bc;
                    cbarstr = 'b^2(f_1,f_2)';
                case {'abs','real','imag','angle'}
                    dum = eval(sprintf('%s(bic.bs)',guy));
                    cbarstr = sprintf('%s%s B(f_1,f_2)',upper(guy(1)),guy(2:end));
                case 'mean'
                    dum = bic.mb;
                    cbarstr = '\langleb^2(f_1,f_2)\rangle';
                case 'std'
                    dum = sqrt(bic.sb);
                    cbarstr = '\sigma_{b^2}(f_1,f_2)';
            end

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
            
            get(h)
            bic.History = get(h,'CData');
        end % PlotBispec
        

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
            WTrim = 50*2;
            dum = bic.sg; 
            if isequal(bic.SpecType,'wave')
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
            bic.sb = bic.sb/(Ntrials-1);
        end % CalcMean

        function bic = PlotConfidence(bic)
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
        
        
        
        
        
        
        
        
        
        

        function bic = PlotGUI(bic)
        % ------------------
        % Convenience
        % ------------------
            bic.Figure = figure('Position',[100 100 600 600],...
                'Resize','off'   ,...
                'WindowKeyPressFcn',@(src,event)SwitchPlot(bic,event),...
                'WindowButtonDownFcn',@(src,event)ClickPlot(bic,event));   
            bic.Axes = axes('DataAspectRatio',[1 1 1],...
                'XLimMode','manual','YLimMode','manual',...
                'DrawMode','fast',...
                'Parent',bic.Figure);
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
                
            % Create the toolbar
            th = uitoolbar(bic.Figure);

            % Add a push tool to the toolbar
            a = [.20:.05:0.95]
            img1(:,:,1) = repmat(a,16,1)'
            img1(:,:,2) = repmat(a,16,1);
            img1(:,:,3) = repmat(flipdim(a,2),16,1);
            pth = uipushtool(th,'CData',img1,...
                       'TooltipString','My push tool',...
                       'HandleVisibility','off')
            % Add a toggle tool to the toolbar
            img2 = rand(16,16,3);
            tth = uitoggletool(th,'CData',img2,'Separator','on',...          
           'TooltipString','Your toggle tool',...
           'HandleVisibility','off')
                
            h = figure('position',[100 100 600 600])
            a1 = axes('position',[0.1 0.1 0.4 0.4])
            a2 = axes('position',[0.5 0.5 0.4 0.4])
                    
        end % PlotGUI
        
        
        
        function ClickPlot(bic,event)
        % For changing desired plottable
        % - - - - - - - - 
            [x,y,button] = ginput(1)
            if button==1
                [~,Ix] = min(abs(bic.fv-x))
                [~,Iy] = min(abs(bic.fv-y))

                v = [1 1 1];
                [w,B,Bi] = bic.GetBispec(bic.sg,v,bic.LilGuy,Iy,Ix,false);
                figure
                plot(unwrap(angle(Bi))/pi)
                %plot(log10(abs(Bi)))
                %ylim([-10 -1])
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
                end
                bic.Processed = bic.Raw;
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
        
                
        function SwitchPlot(bic,event)
        % For changing desired plottable
        % - - - - - - - - 
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
                    bic.PlotBispec;
                otherwise
                    return
            end
        end

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

            B = zeros(floor(lim/2),lim);
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

                    B(j,k) = Bjk;

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
            vec = [-(nfreq-1:-1:1) 0:nfreq-1];
            lim = 2*nfreq-1;

            B = zeros(lim);
            b2 = zeros(lim);
            
            fprintf('Calculating cross-bicoherence...      ')     
            for j=vec
                LoadBar(j+nfreq,lim);             
                for k=vec
                    if abs(j+k) < nfreq
                        p1 = (k>=0)*spec(abs(k)+1,:,v(1)) + (k<0)*conj( spec(abs(k)+1,:,v(1)) );
                        p2 = (j>=0)*spec(abs(j)+1,:,v(2)) + (j<0)*conj( spec(abs(j)+1,:,v(2)) );
                        s  = (j+k>=0)*spec(abs(j+k)+1,:,v(3)) +...
                             (j+k<0)*conj( spec(abs(j+k)+1,:,v(3)) );

                        Bi  = p1.*p2.*conj(s);
                        e12 = abs(p1.*p2).^2;   
                        e3  = abs(s).^2;  

                        Bjk = sum(Bi);                    
                        E12 = sum(e12);             
                        E3  = sum(e3);                      

                        b2(j+nfreq,k+nfreq) = (abs(Bjk).^2)./(E12.*E3+lilguy); 

                        B(j+nfreq,k+nfreq) = Bjk;
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

            p1 = spec(k,:,v(1));
            p2 = spec(j,:,v(2));
            s  = spec(j+k-1,:,v(3));

            % Negatives
            %p2 = conj(spec(j,:));
            %s  = conj(spec(j+k-1,:));

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
        end
        
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
        % sig --> Signal 
        % t   --> Time vector
        % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
            t = 0:1/fS:tend;  % Time-vector sampled at "fS" Hz

            dfx = Afx*sin(2*pi*t*Ff);                     % delta f1
            dfy = Afy*cos(2*pi*t*Ff);                     % delta f2
            x = Ax*sin( 2*pi*(fx*t + dfx) );              % f1
            y = Ay*sin( 2*pi*(fy*t + dfy) );              % f2
            z = Az*sin( 2*pi*(fx*t + fy*t + dfx + dfy) ); % f1 + f2

            sig = x + y + z + noisy*(0.5*rand(1,length(t)) - 1);

        end % SignalGen
             
        function PlotLabels(strings,fsize,cbarNorth)
        % ------------------
        % Convenience function
        % ------------------
            n = length(strings);
            fweight = 'bold';
            xlabel(strings{1},'fontsize',fsize,'fontweight','bold')
            if n>1; ylabel(strings{2},'fontsize',fsize,'fontweight','bold'); end;
            if n>2
                if cbarNorth
                    cbar = colorbar('location','NorthOutside'); 
                    xlabel(cbar,strings{3},'fontsize',fsize,'fontweight','bold')
                else
                    cbar = colorbar;
                    ylabel(cbar,strings{3},'fontsize',fsize,'fontweight','bold')
                end
            end
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
            fftw('planner', 'patient')
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
            [x,t,fS] = bic.TestSignal('3tone');
            %N = 512*1;
            N = length(t);
            x = x(:,1:N);
            dT = t(N)-t(1);
            bic = BicAn(x,'sigma',dT*5,...
                            'spectype','stft',...
                            'sizewarn',false,...
                            'verbose',false,...
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
    ch1 = {'|','|','/','-','\','|','|','|'};  
    ch2 = {'_','.',':','''','^','''',':','.'};
    %fprintf(' Working...      ')
    fprintf('\b\b\b\b\b\b%3.0f%%%s%s',100*m/M,ch1{mod(m,8)+1},ch2{mod(m,8)+1})
end % LoadBar




