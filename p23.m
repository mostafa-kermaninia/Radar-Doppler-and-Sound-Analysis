%% Part 0
close all; 
clc; 
clear; 

%% Part 1: Load Audio and Parameters
[file, path] = uigetfile('*.wav', 'Select an Audio File');
if isequal(file, 0)
    disp('No file selected.');
    return;
end
filePath = fullfile(path, file);
[audioSignal, fs] = audioread(filePath);
signalLength = length(audioSignal);
ts = 1/fs;
timeVector = (0:signalLength-1) * ts; 

% Notes and Frequencies (Combined in a Cell Array)
notesData = {
    'C5', 523.25; 'C#5', 554.37; 'D5', 587.33; 'D#5', 622.25; 'E5', 659.25; 
    'F5', 698.46; 'F#5', 739.99; 'G5', 783.99; 'G#5', 830.61; 'A5', 880.01; 
    'A#5', 932.33; 'B5', 987.77
};

%% Part 2: Extract Segments
startIndices = [];
endIndices = [];
for i = 1:signalLength-1
    if audioSignal(i) == 0 && audioSignal(i+1) ~= 0
        startIndices = [startIndices, i];
    end
    if i > 1 && audioSignal(i) == 0 && audioSignal(i-1) ~= 0
        endIndices = [endIndices, i-1];
    end
end

%% Part 3: Analyze Segments
frequencies = [];
durations = [];
detectedNotes = {};
for i = 1:length(startIndices)
    segment = audioSignal(startIndices(i):endIndices(i));
    duration = (endIndices(i) - startIndices(i) + 1) / fs;
    durations = [durations, duration];

    Y = fftshift(fft(segment));
    Y = abs(Y) / max(abs(Y)); % Normalize FFT
    f = linspace(-fs/2, fs/2, length(Y)); % Frequency vector

    % Find Dominant Frequency
    [~, idx] = max(Y);
    dominantFreq = abs(f(idx));

    % Find Closest Note
    [~, noteIdx] = min(abs(cell2mat(notesData(:, 2)) - dominantFreq));
    detectedNote = notesData{noteIdx, 1};
    frequencies = [frequencies, notesData{noteIdx, 2}];
    detectedNotes = [detectedNotes, detectedNote];

    % Print each note with duration
    fprintf('Note: %s, Frequency: %.2f Hz, Duration: %.2f seconds\n', detectedNote, dominantFreq, duration);
end

%% Part 4: Display Results plot
figure;
stem(1:length(frequencies), frequencies, 'filled');
xlabel('Segment Index');
ylabel('Frequency (Hz)');
title('Dominant Frequencies for Each Segment');
grid on;

