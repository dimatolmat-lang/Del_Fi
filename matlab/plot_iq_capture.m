%% PLOT_IQ_CAPTURE  Просмотр IQ-записи "PELENGATOR IQ CAPTURE V1".
%
% Скрипт (пункт 3 задания): читает pelengator_iq_capture.bin и выводит
% сигнал со всех каналов:
%   1) осциллограммы I и Q всех каналов (первые NPlot отсчётов);
%   2) спектры всех каналов (окно Ханна, полная запись);
%   3) годограф IQ всех каналов (первые NPlot отсчётов);
%   4) огибающая |x(t)| всех каналов.
%
% НАСТРОЙКА: задайте путь к файлу и папке с функциями read_iq_capture.m
% в блоке НАСТРОЙКИ ниже (по умолчанию — папка со скриптом).
%
% Совместимо с MATLAB и GNU Octave (не требует тулбоксов).
% В Octave графики сохраняются в PNG (см. SavePng).

%% ------------------------- НАСТРОЙКИ -----------------------------------
File        = 'pelengator_iq_capture.bin';   % имя файла (или полный путь)
NPlot       = 600;        % сколько отсчётов показывать на осциллограммах
SavePng     = false;      % true -> сохранить окна в PNG рядом со скриптом
%% ------------------------------------------------------------------------

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)   % запуск выделенным кодом (cell mode)
    thisDir = pwd;
end
addpath(thisDir);

if ~exist(File, 'file')
    % ищем файл рядом со скриптом и в ../upload, ../testdata
    cand = { ...
        fullfile(thisDir, File), ...
        fullfile(thisDir, 'upload', File), ...
        fullfile(thisDir, '..', 'upload', File), ...
        fullfile(thisDir, '..', 'testdata', File), ...
        fullfile(thisDir, 'testdata', File)};
    found = '';
    for i = 1:numel(cand)
        if exist(cand{i}, 'file'), found = cand{i}; break; end
    end
    if isempty(found)
        error('Не найден файл %s. Укажите путь в переменной File.', File);
    end
    File = found;
end

%% --- чтение ---------------------------------------------------------------
[iq, hdr] = read_iq_capture(File);
[nCh, nS] = size(iq);
tMs = (0:nS-1) / hdr.fs * 1e3;              % время, мс
fprintf('Файл: %s\n', File);
fprintf('Каналов: %d, отсчётов на канал: %d, fs = %g МГц, fc = %g МГц\n', ...
        nCh, nS, hdr.fs/1e6, hdr.fc/1e6);
if hdr.nFrames > 1
    fprintf('Кадров: %d (показан кадр 1)\n', hdr.nFrames);
end
for c = 1:nCh
    fprintf('  ch%d: RMS = %.1f отсч., max|x| = %.0f\n', ...
            c-1, sqrt(mean(abs(iq(c,:)).^2)), max(abs(iq(c,:))));
end

NPlot = min(NPlot, nS);
win = 0.5 * (1 - cos(2*pi*(0:nS-1)' / (nS-1)));
fAxis = (-floor(nS/2) : ceil(nS/2)-1)' * (hdr.fs / nS);  % порядок fftshift
fAxis = fAxis / 1e6;                        % МГц

%% --- окно 1: осциллограммы I и Q ------------------------------------------
fig1 = figure('Name', 'Осциллограммы I/Q', 'Color', 'w', ...
              'Position', [60 60 980 700]);
for c = 1:nCh
    subplot(nCh, 1, c);
    plot((0:NPlot-1)/hdr.fs*1e6, real(iq(c,1:NPlot)), 'b', ...
         (0:NPlot-1)/hdr.fs*1e6, imag(iq(c,1:NPlot)), 'r');
    grid on;
    ylabel(sprintf('ch%d', c-1));
    if c == 1
        title(sprintf('%s: I (син) и Q (красн), первые %d отсчётов', ...
              strrep(File, '\', '/'), NPlot));
        legend('I', 'Q', 'Location', 'northeast');
    end
    if c == nCh
        xlabel('время, мкс');
    end
end

%% --- окно 2: спектры каналов ----------------------------------------------
fig2 = figure('Name', 'Спектры каналов', 'Color', 'w', ...
              'Position', [120 90 980 700]);
cols = lines(nCh);
leg = cell(1, nCh);
hold on;
for c = 1:nCh
    S = fftshift(fft(iq(c,:).' .* win));
    Sdb = 20*log10(abs(S) / max(abs(S)));   % дБ от максимума канала
    plot(fAxis, Sdb, 'LineWidth', 0.9, 'Color', cols(c,:));
    leg{c} = sprintf('ch%d', c-1);
end
grid on; hold off;
xlabel('частота от несущей, МГц');
ylabel('|S(f)|, дБ от макс. канала');
title('Спектры всех каналов (окно Ханна, вся запись)');
legend(leg, 'Location', 'south');
xlim([fAxis(1) fAxis(end)]);

%% --- окно 3: годографы IQ ---------------------------------------------------
fig3 = figure('Name', 'Годографы IQ', 'Color', 'w', ...
              'Position', [180 120 980 560]);
for c = 1:nCh
    subplot(2, ceil(nCh/2), c);
    plot(real(iq(c,1:NPlot)), imag(iq(c,1:NPlot)), '.', 'MarkerSize', 3);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title(sprintf('ch%d', c-1));
end
if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
    sgtitle(sprintf('Годографы IQ (первые %d отсчётов)', NPlot));
end

%% --- окно 4: огибающая -----------------------------------------------------
fig4 = figure('Name', 'Огибающая |x(t)|', 'Color', 'w', ...
              'Position', [240 150 980 560]);
hold on;
for c = 1:nCh
    plot(tMs, abs(iq(c,:)), 'LineWidth', 0.7, 'Color', cols(c,:));
end
grid on; hold off;
xlabel('время, мс');
ylabel('|x|, отсч. АЦП');
title('Огибающая сигнала всех каналов');
legend(leg, 'Location', 'best');

%% --- сохранение PNG (опция) -------------------------------------------------
if SavePng
    d = fileparts(File);
    if isempty(d), d = pwd; end
    saveas(fig1, fullfile(d, 'view_oscillograms.png'));
    saveas(fig2, fullfile(d, 'view_spectra.png'));
    saveas(fig3, fullfile(d, 'view_iq_traces.png'));
    saveas(fig4, fullfile(d, 'view_envelope.png'));
    fprintf('PNG-файлы сохранены в %s\n', d);
end
