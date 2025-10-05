%% INITIALIZE
clear
close all
clc


%% LOAD DATAFILE


Particiante = 20;
Run         = 2;

Carpeta     = 'Datos';
subcarpeta  = ['P' num2str(Particiante) '_UTB001'];

filename = [ Carpeta '/' subcarpeta '/P' num2str(Particiante) '_UTBS001R0' num2str(Run) '.dat'];
[ signal, states, parameters ] = load_bcidat(filename);


%% CREATE TIME VECTOR

% Number of samples
Nsamples = size(signal,1);

% Number of de channels
Nchannels = size(signal,2);

% Sampling frequency
fs = parameters.SamplingRate.NumericValue;

% Time vector
t  = (0:1:Nsamples-1)'/fs;

% Number of epochs
Nepochs = 80*3;


%% PLOT DATA

figure

subplot(3,1,1), hold on
plot(t,signal)
xlabel('Tiempo (s)')
ylabel('Amplitud (uV)')
title('Señales EEG')
box on, grid on

subplot(3,1,2), hold on
yyaxis left,  plot(t,states.StimulusCode),  ylabel('Stimulus code'),  set(gca,'ylim',[0 25])
yyaxis right, plot(t,states.StimulusBegin), ylabel('Stimulus begin'), set(gca,'ylim',[0  2])
xlabel('Tiempo (s)')
title('Stimulus')
box on, grid on

subplot(3,1,3), hold on
plot(t,states.KeyDown)
xlabel('Tiempo (s)')
ylabel('Key Down')
title('KeyDown')
box on, grid on

% GRAFICA INDIVIDUAL
figure();
plot(t,signal)
xlabel('Tiempo (s)')
ylabel('Amplitud (uV)')
title('Señales EEG')
box on, grid on


figure();
yyaxis left,  plot(t,states.StimulusCode),  ylabel('Stimulus code'),  set(gca,'ylim',[0 25])
yyaxis right, plot(t,states.StimulusBegin), ylabel('Stimulus begin'), set(gca,'ylim',[0  2])
xlabel('Tiempo (s)')
title('Stimulus')
box on, grid on

figure();
plot(t,states.KeyDown)
xlabel('Tiempo (s)')
ylabel('Key Down')
title('KeyDown')
box on, grid on

