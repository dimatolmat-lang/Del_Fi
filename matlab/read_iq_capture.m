function [iq, hdr] = read_iq_capture(filename)
%READ_IQ_CAPTURE  Чтение IQ-кадра формата "PELENGATOR IQ CAPTURE V1".
%
%   [IQ, HDR] = READ_IQ_CAPTURE(FILENAME)
%
%   Возвращает:
%     IQ   — комплексная матрица nFrames x nChannels x nSamples
%            (или nChannels x nSamples, если кадр один);
%     HDR  — структура с полями заголовка и метаданными кадров:
%            version, headerSize, nFrames, nChannels, nSamples, sampleFormat,
%            fc (центральная частота, Гц), fs (частота дискретизации, Гц),
%            reserved, frames.frameId / frames.timeNs / frames.configRev.
%
%   Формат файла (little-endian):
%     Заголовок 56 байт:
%       char[8]  magic "PLNGIQ\0\1"
%       uint32   version = 1
%       uint32   размер заголовка = 56
%       uint32   количество кадров
%       uint32   количество каналов
%       uint32   комплексных отсчётов на канал
%       uint32   код формата = 1 (SC16 LE interleaved IQ)
%       uint64   центральная частота, Гц
%       uint64   частота дискретизации, Гц
%       uint64   reserved
%     Затем для каждого кадра:
%       24 байта метаданных: uint64 frame_id, int64 time_ns, uint64 config_rev
%       затем IQ-данные всех каналов подряд (канал — блок из nSamples
%       комплексных отсчётов int16 I, int16 Q).
%
%   Совместимо с MATLAB и GNU Octave.

    %% --- входной контроль -------------------------------------------------
    if nargin < 1 || isempty(filename)
        error('read_iq_capture:arg', 'Не задано имя файла.');
    end
    if ~exist(filename, 'file')
        error('read_iq_capture:notfound', 'Файл не найден: %s', filename);
    end

    %% --- заголовок файла (56 байт) ---------------------------------------
    fid = fopen(filename, 'r', 'l');   % little-endian
    if fid < 0
        error('read_iq_capture:open', 'Не удалось открыть файл: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    magic = fread(fid, 8, 'uint8=>char')';
    if ~strcmp(magic, sprintf('PLNGIQ\0\1'))
        error('read_iq_capture:magic', ...
              'Неверная сигнатура файла (%s): %s', filename, magic);
    end
    version      = fread(fid, 1, 'uint32=>double');
    headerSize   = fread(fid, 1, 'uint32=>double');
    nFrames      = fread(fid, 1, 'uint32=>double');
    nChannels    = fread(fid, 1, 'uint32=>double');
    nSamples     = fread(fid, 1, 'uint32=>double');
    sampleFormat = fread(fid, 1, 'uint32=>double');
    fc           = fread(fid, 1, 'uint64=>double');
    fs           = fread(fid, 1, 'uint64=>double');
    reserved     = fread(fid, 1, 'uint64=>double');
    if version ~= 1
        error('read_iq_capture:version', ...
              'Неподдерживаемая версия формата: %d', version);
    end
    if sampleFormat ~= 1
        error('read_iq_capture:format', ...
              'Неподдерживаемый код формата отсчётов: %d', sampleFormat);
    end
    if nFrames < 1 || nChannels < 1 || nSamples < 1
        error('read_iq_capture:sizes', ...
              'Некорректные размеры: frames=%d, channels=%d, samples=%d', ...
              nFrames, nChannels, nSamples);
    end

    %% --- проверка размера файла ------------------------------------------
    expected = headerSize + nFrames * (24 + nChannels * nSamples * 4);
    fsize = dir(filename).bytes;
    if fsize < expected
        error('read_iq_capture:size', ...
            ['Файл короче, чем требует заголовок: %d байт, ожидалось %d.\n' ...
             'Возможно, файл оборван или формат изменился.'], fsize, expected);
    end
    if fsize > expected
        warning('read_iq_capture:extra', ...
            'Файл на %d байт длиннее, чем требует заголовок: лишние данные проигнорированы.', ...
            fsize - expected);
    end

    %% --- чтение кадров ------------------------------------------------------
    % 3-D массив: [nFrames, nChannels, nSamples]
    iq = zeros(nFrames, nChannels, nSamples);
    frameId   = zeros(1, nFrames);
    timeNs    = zeros(1, nFrames);
    configRev = zeros(1, nFrames);

    for f = 1:nFrames
        frameId(f)   = fread(fid, 1, 'uint64=>double');
        timeNs(f)    = fread(fid, 1, 'int64=>double');
        configRev(f) = fread(fid, 1, 'uint64=>double');

        % все каналы кадра одним блоком: I/Q чередуются, каналы подряд
        raw = fread(fid, nChannels * nSamples * 2, 'int16=>double');
        m = reshape(raw, 2, nSamples, nChannels);   % 1-я строка I, 2-я Q
        for c = 1:nChannels
            iq(f, c, :) = m(1, :, c) + 1i * m(2, :, c);
        end
    end

    hdr = struct( ...
        'version',       version, ...
        'headerSize',    headerSize, ...
        'nFrames',       nFrames, ...
        'nChannels',     nChannels, ...
        'nSamples',      nSamples, ...
        'sampleFormat',  sampleFormat, ...
        'fc',            fc, ...
        'fs',            fs, ...
        'reserved',      reserved);
    hdr.frames = struct('frameId', frameId, 'timeNs', timeNs, ...
                        'configRev', configRev);

    if nFrames == 1
        iq = reshape(iq, nChannels, nSamples);   % одиночный кадр -> 2-D
    end
end
