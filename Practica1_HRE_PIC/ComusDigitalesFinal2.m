% -	Hernández Romero Enrique
% - Pérez Iturbe Carolina
%% PARTE I.  GENERACIÓN DE LA SEÑAL DE VOZ "ANALÓGICA".

clc
close all

Fs = 16000; 
nBits = 16;
nChannels = 1; 
Duration = 20; 

% % Crear el objeto de grabación
% recObj = audiorecorder(Fs, nBits, nChannels);
% 
% % Mensaje para el usuario
% disp('Iniciando la grabación por 20 segundos...');
% 
% % Iniciar la grabación 
% recordblocking(recObj, Duration);
% 
% % Mensaje para el usuario
% disp('Grabación completada.');
% 
% % Extraer los datos de la grabación como un vector de int16
% audioData = getaudiodata(recObj, 'int16');

% 'audioData' contiene las muestras de audio.
% x(t) corresponde a 'audioData'

%% Tiempo

% Definir el vector de tiempo
numMuestras = length(audioData);
t = (0:numMuestras-1) / Fs; 

% Graficar el audio en el dominio del tiempo
figure(1)
plot(t, audioData);
title('Se\~{n}al de audio grabada en el dominio del tiempo','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Amplitud','Interpreter','latex');
grid on;

% Guardar el archivo de audio
audiowrite('grabacionOriginal16kHz.wav', audioData, Fs);

%% PARTE II.  ANÁLISIS DE LA SEÑAL DE VOZ "ANALÓGICA".


LongitudFragmento = 5000; 

% Seleccionar el fragmento de 5000 muestras
Fragmento = audioData(5000:10000);

% Crear el vector de tiempo asociado al fragmento
TiempoFragmento = (0:LongitudFragmento) / Fs;

% Graficar el fragmento en el dominio del tiempo
figure(2);
subplot(4,2,1);
plot(TiempoFragmento, Fragmento);
title('Fragmento de la se\~{n}al de voz a 16 kHz','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Amplitud','Interpreter','latex');
grid on

% Calcular la Transformada de Fourier de toda la señal
Xf = fft(audioData);
% Aplicar fftshift para centrar el espectro
Xfshifted = fftshift(Xf);

numMuestras = length(audioData);
f = (-numMuestras/2:numMuestras/2-1) * (Fs / numMuestras);

% Graficar la magnitud del espectro centrado
figure(2);
subplot(4,2,2);
plot(f, abs(Xfshifted)); 
title('Espectro centrado de la se\~{n}al de voz a 16 kHz','Interpreter','latex');
xlabel('Frecuencia (Hz)','Interpreter','latex');
ylabel('Magnitud','Interpreter','latex');
xlim([-Fs/2, Fs/2]);
grid on

% Calcular la función de densidad de probabilidad usando un histograma
% normalizado [De toda la señal y no del fragmento]
figure(3);
h = histogram(audioData, 'Normalization', 'pdf');
title('Funci\''on de Densidad de Probabilidad (PDF)','Interpreter','latex');
xlabel('Amplitud','Interpreter','latex');
ylabel('Probabilidad','Interpreter','latex');
grid on

% Obtener los valores del histograma
counts = h.Values;  
binEdges = h.BinEdges; 

% Calcular los puntos medios de cada bin
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Calcular el área bajo la curva usando la regla del trapecio
area = trapz(binCenters, counts);

% Mostrar el valor del área calculada
disp(['El área bajo la curva de la función de densidad de probabilidad es: ', num2str(area)]);


%% PARTE III. ANÁLISIS DE VELOCIDAD DE MUESTREO.

% Frecuencias de muestreo deseadas
Fs8kHz = 8000; % 8 kHz
Fs4kHz = 4000; % 4 kHz
Fs2kHz = 2000; % 2 kHz

% Submuestrear la señal para cada frecuencia de muestreo
audioData8kHz = downsample(audioData, Fs/Fs8kHz);
audioData4kHz = downsample(audioData, Fs/Fs4kHz);
audioData2kHz = downsample(audioData, Fs/Fs2kHz);

%% PARA 8 kHz

inicio = 5000;
final = 10000;
LongitudFragmento8kHz = final - inicio + 1;
Fragmento8kHz = audioData8kHz(inicio:final);

% Crear el vector de tiempo asociado al fragmento
TiempoFragmento8kHz = (0:LongitudFragmento8kHz-1) / Fs8kHz;

% Graficar el fragmento en el dominio del tiempo
figure(2);
subplot(4,2,3);
plot(TiempoFragmento8kHz, Fragmento8kHz);
title('Fragmento de la se\~{n}al de voz a 8 kHz','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Amplitud','Interpreter','latex');
grid on

Xf8kHz = fft(Fragmento8kHz);
Xf8kHz_shifted = fftshift(Xf8kHz);

% Crear el vector de frecuencias simétrico desde -Fs/2 a Fs/2
f8kHz = (-LongitudFragmento8kHz/2 : LongitudFragmento8kHz/2-1) * (Fs8kHz / LongitudFragmento8kHz);

figure(2);
subplot(4,2,4);
plot(f8kHz, abs(Xf8kHz_shifted));
title('Espectro centrado de la se\~{n}al de voz a 8 kHz ','Interpreter','latex');
xlabel('Frecuencia (Hz)','Interpreter','latex');
ylabel('Magnitud','Interpreter','latex');
xlim([-8000 8000]);
grid on

%% PARA 4 kHz

LongitudFragmento4kHz = final - inicio + 1;
Fragmento4kHz = audioData4kHz(inicio:final);

% Crear el vector de tiempo asociado al fragmento
TiempoFragmento4kHz = (0:LongitudFragmento4kHz-1) / Fs4kHz;

% Graficar el fragmento en el dominio del tiempo
figure(2);
subplot(4,2,5);
plot(TiempoFragmento4kHz, Fragmento4kHz);
title('Fragmento de la se\~{n}al de voz a 4 kHz','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Amplitud','Interpreter','latex');
grid on

% Obtener el espectro
Xf4kHz = fft(Fragmento4kHz);
Xf4kHz_shifted = fftshift(Xf4kHz);

% Crear el vector de frecuencias simétrico desde -Fs/2 a Fs/2
f4kHz = (-LongitudFragmento4kHz/2 : LongitudFragmento4kHz/2-1) * (Fs4kHz / LongitudFragmento4kHz);

figure(2);
subplot(4,2,6);
plot(f4kHz, abs(Xf4kHz_shifted));
title('Espectro centrado de la se\~{n}al de voz a 4 kHz','Interpreter','latex');
xlabel('Frecuencia (Hz)','Interpreter','latex');
ylabel('Magnitud','Interpreter','latex');
xlim([-8000, 8000]);
grid on

%% PARA 2 kHz

LongitudFragmento2kHz = final - inicio + 1;
Fragmento2kHz = audioData2kHz(inicio:final);
% Crear el vector de tiempo asociado al fragmento
TiempoFragmento2kHz = (0:LongitudFragmento2kHz-1) / Fs2kHz;

% Graficar el fragmento en el dominio del tiempo
figure(2);
subplot(4,2,7);
plot(TiempoFragmento2kHz, Fragmento2kHz);
title('Fragmento de la se\~{n}al de voz a 2 kHz','Interpreter','latex');
xlabel('Tiempo (segundos)','Interpreter','latex');
ylabel('Amplitud','Interpreter','latex');
grid on

% Obtener el espectro
Xf2kHz = fft(Fragmento2kHz);
Xf2kHz_shifted = fftshift(Xf2kHz);

% Crear el vector de frecuencias simétrico desde -Fs/2 a Fs/2
f2kHz = (-LongitudFragmento2kHz/2 : LongitudFragmento2kHz/2-1) * (Fs2kHz / LongitudFragmento2kHz);

figure(2);
subplot(4,2,8);
plot(f2kHz, abs(Xf2kHz_shifted));
title('Espectro centrado de la se\~{n}al de voz a 2 kHz','Interpreter','latex');
xlabel('Frecuencia (Hz)','Interpreter','latex');
ylabel('Magnitud','Interpreter','latex');
xlim([-8000, 8000]); 
grid on

% Guardar los archivos de audio submuestreados
audiowrite('grabacion8kHz.wav', audioData8kHz, Fs8kHz);
audiowrite('grabacion4kHz.wav', audioData4kHz, Fs4kHz);
audiowrite('grabacion2kHz.wav', audioData2kHz, Fs2kHz);

%% PARTE IV. ANÁLISIS DE BITS POR MUESTRA

% Definir los bits por muestra a usar
bits_per_sample = [14, 8, 6];

min_val = min(audioData);
max_val = max(audioData);
figure;
for i = 1:length(bits_per_sample)
    bits = bits_per_sample(i);
    num_levels = 2^bits;  
    
    step_size = (max_val - min_val) / num_levels;
    cuantizada = round((audioData - min_val) / step_size) * step_size + min_val;
    
    % Guardar el archivo de audio cuantizado
    filename = sprintf('audio_cuantizada_%dbits.wav', bits);
    audiowrite(filename, cuantizada, Fs);
    
    % Graficar fragmento cuantizado
    fragmento = cuantizada(5000:10000);
    tiempo_fragmento = (0:length(fragmento)-1) * (1/Fs);

    subplot(3,2,2*i-1)
    plot(tiempo_fragmento, fragmento);
    xlabel('Tiempo (s)','Interpreter','latex');
    ylabel('Amplitud','Interpreter','latex');
    title(['Fragmento de la se\~{n}al cuantizada a ', num2str(bits), ' bits por muestra'],'Interpreter','latex');
    grid on

    % Obtener el espectro
    X_f = fft(fragmento);
    X_f_shifted = fftshift(X_f);  
    magnitude_X_f_shifted = abs(X_f_shifted);
    % Crear el vector de frecuencias 
    num_muestras = length(fragmento);
    frequencies_shifted = (-num_muestras/2:num_muestras/2-1) * (Fs / num_muestras);


    subplot(3,2,2*i)
    plot(frequencies_shifted, magnitude_X_f_shifted);
    xlim([-8000 8000]);
    xlabel('Frecuencia (Hz)','Interpreter','latex');
    ylabel('Magnitud','Interpreter','latex');
    title(['Espectro de la se\~{n}al cuantizada a ', num2str(bits), ' bits por muestra'],'Interpreter','latex');
    grid on
end 
