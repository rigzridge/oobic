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
% Y         -> time-series {or structure}
% samprate  -> sampling rate [Hz]
% res       -> desired frequency resolution* [Hz]
% subint    -> subinterval size** [samples]
% step      -> step size for Welch method [samples]
% - - - - - - - - - - - - - - - - - - - - 
% additional options... (see below for instructions)
% - - - - - - - - - - - - - - - - - - - - 
% window    -> select window function [default :: "hann" (see @window)]
% justfft   -> true for just spectrogram [default :: false]
% errlim    -> mean(fft) condition [default :: inf] 
% fscale    -> scale for plotting frequencies [default :: 1]
% dealias   -> applies antialiasing (LP) filter [default :: false]
% bispectro -> computes bispectrogram [default :: false]
% smooth    -> smooths FFT by n samples [default :: 1]
% plot      -> start plotting tool when done [default :: false]
% lilguy    -> set epsilon [default :: 1e-6]
% sizewarn  -> warning for matrix size [default :: true]
% cmap      -> adjust colormap [default :: parula]
% plottype  -> set desired plottable [default :: 'bicoh']
% autoscale -> autoscaling in figures [default :: false]
% verbose   -> allow printing of info structure [default :: true]
% detrend   -> remove linear trend from data [default :: false]
% zpad      -> add zero-padding to end of time-series [default :: true]
% note      -> optional string for notes [default :: ' '] 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Version History
% - - - - - - - - - - - - - - - - - - - - 
% 6/27/2022 -> First code. Inspired to do O.O. b/c of Phase Space Synth!
% Honestly think that it will be glove in hand with a project like this.
% Anyway: I'm trying to write a nice, clean constructor. [...]


classdef BicAn < handle
% Bicoherence analysis class for DSP

    % Globals 
    properties (Constant)
        FontSize = 20;
    end

    % Dependents
    properties (Dependent)
        MaxRes
        Samples
    end

    % Editables
    properties
        Raw
        Processed
        History
        Samples
        SampRate
        FreqRes
        SubInt
        Step
        Window
        JustFFT
        ErrLim
        FScale
        TScale
        Filter
        Bispectro
        Smooth
        PlotIt
        LilGuy
        SizeWarn
        CMap
        PlotType
        ScaleAxes
        Verbose
        Detrend
        ZPad
        Note
        Cross
        Vector
        StartTime
        tv
        fv
        ff
        ft
        sg
        xs
        xc
        cs
        bs
        bc
        bp
        bg
        er
    end % properties

    % Functions
    methods 

        function bic = BicAn(varargin)
        % ------------------
        % Constructor
        % ------------------
            if nargin==0
                bic.Raw       = [];
                bic.Processed = [];
                bic.History   = ' ';
                bic.Samples   = [];
                bic.SampRate  = [];
                bic.FreqRes   = [];
                bic.SubInt    = 512;
                bic.Step      = 128;
                bic.Window    = 'hann';
                bic.JustFFT   = false;
                bic.ErrLim    = inf;
                bic.FScale    = 1;
                bic.TScale    = 1;
                bic.Filter    = 'none';
                bic.Bispectro = false;
                bic.Smooth    = 1;
                bic.PlotIt    = false;
                bic.LilGuy    = 1e-6;
                bic.SizeWarn  = true;
                bic.CMap      = 'viridis';
                bic.PlotType  = 'bicoh';
                bic.ScaleAxes = 'manual';
                bic.Verbose   = true;
                bic.Detrend   = false;
                bic.ZPad      = true;
                bic.Note      = ' ';
                bic.Cross     = false;
                bic.Vector    = false;
                bic.StartTime = 0;
                bic.tv = []; % Time vector
                bic.fv = []; % Frequency vector
                bic.ff = []; % Full frequency vector
                bic.ft = []; % Fourier amplitudes
                bic.sg = []; % Spectrogram (complex)
                bic.xs = []; % Cross-spectrum
                bic.xc = []; % Cross-coherence
                bic.cs = []; % Coherence spectrum
                bic.bs = []; % Bispectrum
                bic.bc = []; % Bicoherence spectrum
                bic.bp = []; % Biphase proxy
                bic.bg = []; % Bispectrogram
                bic.er = []; % Mean & std dev of FFT
            else
                bic = ParseInput(varargin);
            end


        end % BicAn

        % "Get properties" functions
        function val = get.MaxRes(bic)
            val = bic.SampRate / bic.SubInt;
        end
        function val = get.Samples(bic)
            val = length(bic.Processed);
        end

        % "Set properties" functions
        function bic = set.FreqRes(bic,val)
            if isequal(val,'max')
                bic.FreqRes = bic.MaxRes;
            elseif val>0 && val<bic.MaxRes
                bic.FreqRes = val;
            else
                warning('Bispec:resWarn','Requested resolution not possible... Using maximum.')
                bic.FreqRes = bic.MaxRes;
            end
        end


        function bic = ParseInput(bic,varargin)
        % ------------------
        % Handle inputs
        % ------------------
        % 
            N = length(varargin);
            if N>5
                for i=1:2:(N-5)
                    switch lower(varargin{i})            
                        case inputs
                            for j=1:length(options)                       
                                if isequal(inputs{j},lower(varargin{i}))
                                    cl = eval(sprintf('class(%s)',options{j}));                    
                                    if isequal(cl,class(varargin{i+1}))
                                        eval(sprintf('bic.%s = varargin{i+1};',options{j}));
                                    else
                                        warning('BicAn:errOption','\n"%s" must be a %s... Using default value.',...
                                                    varargin{i},cl)
                                    end
                                end
                            end             
                        otherwise
                            warning('BicAn:unknownOption','\n"%s" is not a known option...',varargin{i})
                    end
                end
            end

        end % ParseInput


        function bic = Filter(bic)
        % ------------------
        % LP/BP/HP filter
        % ------------------

        end % Filter


        function bic = SpectroSTFT(bic)
        % ------------------
        % STFT method
        % ------------------

        end % SpectroSTFT


        function bic = SpectroWavelet(bic)
        % ------------------
        % Wavelet method
        % ------------------

        end % SpectroWavelet


        function bic = Coherence(bic)
        % ------------------
        % Cross-spectrum/coh
        % ------------------

        end % Coherence

        
        function bic = Bicoherence(bic)
        % ------------------
        % STFT method
        % ------------------

        end % Bicoherence


        function bic = Confidence(bic)
        % ------------------
        % Find confidence interval
        % ------------------

        end % Confidence


        function bic = Plot(bic)
        % ------------------
        % Plot stuff
        % ------------------

        end % Plot
        

        function bic = PlotGUI(bic)
        % ------------------
        % Convenience
        % ------------------

        end % PlotGUI
        

        function bic = RunDemo(bic)
        % ------------------
        % Demonstration
        % ------------------

        end % RunDemo
        

        function bic = MakeMovie(bic)
        % ------------------
        % Output movie
        % ------------------

        end % MakeMovie

    end % methods


    methods (Static)

        function y = Detrend(y)
        % ------------------
        % Remove linear trend
        % ------------------
           n = length(y);
           s = (6/(n*(n^2-1)))*(2*sum((1:n).*y) - sum(y)*(n+1));
           % Convenient form assuming x = 1:n
           y = y - s*(1:n);
        end % Detrend

    end % static methods

end % BicAn




