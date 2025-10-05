%% INITIALIZE
clear
close all
clc

%% LOAD DATAFILE

Particiante = 20 % ===============================================> SET THIS

Run         = 2;
Carpeta     = 'Datos';

% P0xx_UTB00y
subcarpeta  = ['P' num2str(Particiante) '_UTB001'];

filename = [ Carpeta '/' subcarpeta '/P' num2str(Particiante) '_UTBS001R0' num2str(Run) '.dat'];
[ signal, states, parameters ] = load_bcidat(filename);

%filename = [ Carpeta '/P0' num2str(Particiante) '_UTBS001R01.dat'];
%[ signal, states, parameters ] = load_bcidat(filename);


%% CREATE TIME VECTOR

% Number of samples
Nsamples = size(signal,1);

% Number of de channels
Nchannels = size(signal,2);

% OSampling frequency
fs = parameters.SamplingRate.NumericValue;

% Time vector
t  = (0:1:Nsamples-1)'/fs;

% % Number of epochs
% Nepochs = 80*4; % ===============================================> SET THIS


%% CALCULAR Y VERIFICR EL NUMERO DE ETIMULOS Y POR ENDE DE EPOCAS

% Define a vector with ones where se presenta el estimulo (i.e., donde the epoch begins)
states.StimulusBegin2 = diff(states.StimulusBegin);
states.StimulusBegin2 = [0 ; states.StimulusBegin2];
% Esto es para confirmar el "StimulusBegin"

% Obtener el number of epochs
Nepochs = sum(states.StimulusBegin2);
fprintf('El numero de estimulos (o epocas) es: %d \n',Nepochs)

% Confirmar que el numero de epocas sea 320, si no, error
if Nepochs~=320
    error('PILAS PERRITO: el numero de estimulos en los datos esta mal')
end
% Si hay error: determinar porque se da ese error


%% PLOT DATA

% figure, hold on
% plot(states.StimulusBegin)
% plot(states.StimulusBegin2)
% ylabel('Stimulus begin')
% xlabel('Tiempo (s)')
% set(gca,'ylim',[0  1.1])
% box on, grid on


%% SAMPLE INITATION OF EACH EPOCH

% Sample of the initiation of each epoch
Sini = find(states.StimulusBegin2==1);

% Verify
if length(Sini)~=Nepochs
    error('PILAS: the length of Sini is not Nepochs')
end


%% GET EPOCHS/TRIALS

% Epoch time information
tini       = -1.0;% =============================================> SET THIS
tfin       = +5.0;% =============================================> SET THIS

duration   = tfin - tini;
numsamples = duration*fs + 1;

% Esto que es?
INC4KeyDown = 0; % Muestras adicional despues de termianr la epoca

% Initialize 3D matrix for the epochs --> Nepochs ; Nsamples ; Nchannles
epochs        = zeros( Nepochs, numsamples, Nchannels );
times         = zeros( Nepochs, numsamples );
stimcode      = zeros( Nepochs, 1 );
KeyDown       = zeros( Nepochs, numsamples+INC4KeyDown);
responsetime  = zeros( Nepochs, 1 );
epochtype     = zeros( Nepochs, 1 ); % 1: epochtype; 0: inepochtype

% Get epoch
for i = 1:1:Nepochs
    %fprintf("Numero de eepoca")
    disp(i)
    
    % Get sampleini and samplefin
    sampini = Sini(i)-abs(tini)*fs;
    sampfin = Sini(i)+abs(tfin)*fs+1;
    if (sampfin-sampini)~=numsamples
        error('PILAS: el numero de muestras de la epoca actual no esta bien')
    end
    
    % Get epoch
    epochs(i,:,:)   = signal(sampini:sampfin-1,:);
    % Get time
    times(i,:,:)    = t(sampini:sampfin-1);
    
    % Get stimcode of the epoch
    stimcode(i)     = states.StimulusCode( Sini(i) );
    
    % Get KeyDown signal in the entire epoch
    KeyDown(i,:)    = states.KeyDown(sampini:sampfin-1+INC4KeyDown);
    % Compute response time
    if isempty(find(KeyDown(i,:)~=0))
        responsetime(i) = inf;
    else
        SampleKeyDown   = find(KeyDown(i,:)~=0);
        SampleKeyDown   = SampleKeyDown(1);
        responsetime(i) = (SampleKeyDown)/fs - 0.5; % tini ****************
    end
    
    if (stimcode(i)<=4)
        epochtype(i) = 1; % congruente
    end
    
    % Plot for debugging
    if (0)        
        figure(2), clf
        
        subplot(3,1,1), hold on
        yyaxis left
        plot(t,signal,'-','color',[0.75,0.75,0.75])
        plot(times(i,:,:), squeeze(epochs(i,:,:)),'-' )
        ylabel('Amplitud (uV)')
        yyaxis right
        plot(t,states.StimulusBegin2)
        ylabel('Stimulus begin')
        set(gca,'ylim',[-0.1  1.1])
        
        subplot(3,1,2), hold on
        yyaxis left
        plot(t,signal,'-','color',[0.75,0.75,0.75])
        plot(times(i,:,:), squeeze(epochs(i,:,:)),'-' )
        ylabel('Amplitud (uV)')
        set(gca,'xlim',[ min(times(i,:,:)) max(times(i,:,:)) ] + [ -0.5 +0.5])
        yyaxis right
        plot(t,states.StimulusBegin2)
        ylabel('Stimulus begin')
        set(gca,'ylim',[-0.1  1.1])
        set(gca,'xlim',[ min(times(i,:,:)) max(times(i,:,:)) ] + [ -0.5 +0.5])        
        xlabel('Tiempo (s)')
        title('Señales EEG')
        box on, grid on
        
        subplot(3,1,3), hold on
        plot(times(i,:,:),KeyDown(i,:))
        title(['Response time = ' num2str(responsetime(i)) 's'])
        xlabel([ 'StimCode=' num2str(stimcode(i)) '    |    epochtype=' num2str(epochtype(i))])
        set(gca,'xlim',[ min(times(i,:,:)) max(times(i,:,:)) ] + [ -0.5 +0.5]) 
        box on, grid on
        
        pause
    end % if (1)
    
end % for i = 1:1:Nepochs

% Create a commun time vector across all epochs
t = squeeze(times(1,:,:));
t = t-t(1)-1; % PILAS PERRITO *******************************************

% Clear garbage
clear ans duration i Nchannels Nsamples numsamples parameters
clear sampfin sampini signal Sini states tfin times tini SampleKeyDown
clear KeyDown


%% OBTENER VARIABLE DE FIELDTRIP

% Agregar fieldtrip a matlab
%path(path,'C:\Users\L01046417\Documents\fieldtrip-20240211')


% Inicializar estructura de fieldtrip
EEG             = [];
EEG.fsample     = fs;
EEG.label       = {'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'C3';'C4';'T5';'P3';'Pz';'P4';'T6';'O1';'O2'};
EEG.trial       = cell(1,Nepochs);
EEG.time        = cell(1,Nepochs);


% Meter en cada celda, la matriz de datos de la epoca/trial correspondiente
for i=1:1:Nepochs
    EEG.trial{i} = squeeze(epochs(i,:,:))';
    EEG.time{i} = t;
end

% Dummy processing in fieldtrip para asegurarse de que todo quedo bien
%cfg             = [];
%cfg.feedback    = '';
%EEG             = ft_preprocessing(cfg,EEG);

% Crear SubjectInfo
SubjectInfo     = [];
SubjectInfo.responsetime = responsetime;
SubjectInfo.epochtype    = epochtype;
SubjectInfo.stimcode     = stimcode;


%% SAVE DATA

filename2save = [ Carpeta '/' subcarpeta '/P' num2str(Particiante) '_DataEpochs'];

save(filename2save,'epochs','t','epochtype','stimcode','responsetime') %,'EEG','SubjectInfo'
