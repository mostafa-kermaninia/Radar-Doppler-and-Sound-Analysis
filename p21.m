%% Part 0
close all; 
clc; 
clear; 

%% Part 1: Parameters
fs = 8000;
ts = 1/fs;
T = 0.5; 
T_half = T / 2; 
tau = 0.025; % Time gap between notes 

% Frequencies of music keys that are used in this part
frequencies = {'D5', 587.33; 'E5', 659.25; 'F#5', 739.99; 'G5', 783.99};

t_full = 0:ts:T-ts;
t_half = 0:ts:T_half-ts; 
silence = zeros(1, tau/ts); % Silence for gaps

final_signal = [];

%% Part 3: Build the Song
% Bar 1
final_signal = [final_signal, ...
    append_note('D5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('G5', t_full, frequencies, silence), append_note('F#5', t_full, frequencies, silence), ...
    append_note('D5', t_full, frequencies, silence)];

% Bar 2
final_signal = [final_signal, ...
    append_note('D5', t_half, frequencies, silence), append_note('E5', t_half, frequencies, silence), ...
    append_note('E5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('F#5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('E5', t_full, frequencies, silence)];

% Bar 3
final_signal = [final_signal, ...
    append_note('D5', t_full, frequencies, silence), append_note('E5', t_full, frequencies, silence), ...
    append_note('F#5', t_full, frequencies, silence), append_note('E5', t_full, frequencies, silence)];

% Bar 4
final_signal = [final_signal, ...
    append_note('D5', t_half, frequencies, silence), append_note('E5', t_half, frequencies, silence), ...
    append_note('E5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('F#5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('E5', t_full, frequencies, silence)];

% Bar 5
final_signal = [final_signal, ...
    append_note('D5', t_full, frequencies, silence), append_note('E5', t_half, frequencies, silence), ...
    append_note('D5', t_half, frequencies, silence), append_note('F#5', t_full, frequencies, silence), ...
    append_note('E5', t_full, frequencies, silence)];

% Bar 6
final_signal = [final_signal, ...
    append_note('D5', t_full, frequencies, silence), append_note('E5', t_half, frequencies, silence), ...
    append_note('D5', t_half, frequencies, silence), append_note('F#5', t_full, frequencies, silence), ...
    append_note('E5', t_full, frequencies, silence)];

% Bar 7
final_signal = [final_signal, ...
    append_note('D5', t_half, frequencies, silence), append_note('D5', t_half, frequencies, silence), ...
    append_note('E5', t_full, frequencies, silence), append_note('F#5', t_half, frequencies, silence), ...
    append_note('E5', t_half, frequencies, silence), append_note('F#5', t_full, frequencies, silence)];

% Bar 8
final_signal = [final_signal, ...
    append_note('F#5', t_half, frequencies, silence), append_note('E5', t_half, frequencies, silence), ...
    append_note('F#5', t_full, frequencies, silence), append_note('F#5', t_full, frequencies, silence), ...
    append_note('D5', t_full, frequencies, silence)];

%% Part 4: Play the Song 
sound(final_signal, fs);
audiowrite('CA6_song.wav', final_signal, fs);

%% Part 2: Helper Function to append a note with specified duration
function note_signal = append_note(note, duration, frequencies, silence)
    freq = frequencies{strcmp(frequencies(:,1), note), 2};
    note_signal = [sin(2*pi*freq*duration), silence];
end