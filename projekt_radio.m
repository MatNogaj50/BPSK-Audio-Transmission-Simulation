%% Restart

clear; clc; close all

%% Data Upload

[Data, Fsmp] = audioread('media/muzyka_projekt.wav');    % User input
DataMono = mean(Data, 2);

%% Start moment

Len = length(DataMono);
StartingPoint = 0.9;    % User input 0 - 1
EndingPoint = 0.95;     % User input 0 - 1
TrimmedData = DataMono(round(StartingPoint*Len) : round(EndingPoint*Len));

%% Play

audioplay = audioplayer( TrimmedData, Fsmp );
playblocking(audioplay);

%% Filtering & initial compression

CompressionRatio = 2;    % User input
FsmpFiltered = Fsmp/CompressionRatio;
CutOffFrequency = FsmpFiltered/2;

Spectrum = fft(TrimmedData);
N = length(TrimmedData);
% Only first half
f = 0 : Fsmp / N : Fsmp/N * (N/2 - 1);
FirstHalfSpectrum = abs(Spectrum(1:N/2));
magSpectrumdB = 20 * log10(FirstHalfSpectrum / max(FirstHalfSpectrum) + eps);

FilterOrder = 10; % User input
FilterType = 'butterworth'; % Options: 'butterworth', 'hamming', 'blackman'
switch lower(FilterType)
    case 'butterworth'
        [b, a] = butter(FilterOrder, CutOffFrequency / (Fsmp/2), 'low');
        FilteredData = filtfilt(b, a, TrimmedData);
    case 'hamming'
        b_Hamming = fir1(FilterOrder, CutOffFrequency / (Fsmp/2), 'low', hamming(FilterOrder + 1));
        FilteredData = filtfilt(b_Hamming, 1, TrimmedData);
    case 'blackman'
        b_Blackman = fir1(FilterOrder, CutOffFrequency / (Fsmp/2), 'low', blackman(FilterOrder + 1));
        FilteredData = filtfilt(b_Blackman, 1, TrimmedData);
    otherwise
        error('Unknown filter type selected.');
end

%% Filtering & initial compression continued

FilteredSpectrum = fft(FilteredData);
% Only first half
FilteredHalfSpectrum = abs(FilteredSpectrum(1:N/2));
magSpectrumFiltereddB = 20 * log10(FilteredHalfSpectrum / max(FilteredHalfSpectrum) + eps);

figure(1)
subplot(2,1,1)
semilogx(f, magSpectrumdB)
title('Frequency Spectrum')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([0, Fsmp/2]);
ylim([-160, 20]);

subplot(2,1,2)
semilogx(f, magSpectrumFiltereddB)
title('Filtered frequency Spectrum')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([0, Fsmp/2]);
ylim([-160, 20]);

FilteredCompressedData = FilteredData(1 : CompressionRatio : end);

%% Play

audioplay = audioplayer( FilteredCompressedData, FsmpFiltered );
playblocking(audioplay);

%% Final compression

DctSpectrum = dct( FilteredCompressedData );
[ sortedSpectrum, indexes ] = sort( abs( DctSpectrum),'descend');

quality = 0.995;    % User input

TotalEnergy = sum(DctSpectrum.^2);
CumulativeEnergy = cumsum(DctSpectrum(indexes).^2);
Idx = find(CumulativeEnergy >= quality * TotalEnergy, 1);

CompressedSpectrum = zeros(size(DctSpectrum));
CompressedSpectrum(indexes(1:Idx)) = DctSpectrum(indexes(1:Idx));

[ sortedSpectrumC, indexesC ] = sort( abs( CompressedSpectrum),'descend');

figure(2)
stem(abs(sortedSpectrum))
title('Sorted DCT Spectrum before and after compression')
hold on;
stem(abs(sortedSpectrumC))
legend( 'Before', 'After' );

CompressedData = idct(CompressedSpectrum);

%%

NormalizedEnergy = (CumulativeEnergy / TotalEnergy) * 100;

figure(3)
plot(NormalizedEnergy, 'b-', 'LineWidth', 1.5);
hold on;
% Cut-off level
plot(Idx, NormalizedEnergy(Idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
yline(quality * 100, '--g', sprintf('Quality: %.1f%%', quality * 100), 'LineWidth', 1.2);

grid on;
title('Cumulative energy of sorted DCT coefficients');
xlabel('Number of coefficients');
ylabel('Energy [%]');
legend('Cumulative energy', sprintf('Cut-off level (N = %d)', Idx), 'Specified level', 'Location', 'best');
ylim([75 101]);

%% Play

audioplay = audioplayer( CompressedData, FsmpFiltered );
playblocking(audioplay);

%% ADC

bits = 8;    % User input
% -1 <= CompressedData <= 1
NormalizedData = 0.5 * (CompressedData + 1) * (2^bits - 1);
NormalizedData = round(NormalizedData);

% Previous unoptimized version
% BinaryData = [];
% for i = 1:length(NormalizedData)
%     Sample = NormalizedData(i);
%     for j = 1:bits
%         bit = mod(Sample, 2);
%         % LSB first
%         BinaryData = [BinaryData; bit];
%         Sample = floor( Sample/2 );
%     end
% end

% Optimized version
BitPowers = 2.^(0:bits-1);
BinaryData = bitand(NormalizedData, BitPowers) > 0;
BinaryData = double(BinaryData');
BinaryData = BinaryData(:);

%% Modulation BPSK

CarrierLength = 10;    % User input, even number
t = 0 : 1/CarrierLength : 1 - 1/CarrierLength;
Carrier = sin(2*pi*t.');

% Previous unoptimized version
% ModulatedSignal = [];
% for i = 1:length(BinaryData)
%     % 0 - negated Carrier, 1 - Carrier
%     ModulatedSignal = [ ModulatedSignal ; Carrier * (-1)^(BinaryData(i)+1) ];
% end

% Optimized version
% 0 => -1, 1 => 1 (0 - negated Carrier, 1 - Carrier)
PhaseIndicator = 2 * BinaryData - 1;
ModulatedSignal = kron(PhaseIndicator, Carrier);

%% Transmission (attenuation & adding noise)

NoiseAmplitude = 10;    % User input [%]
Attenuation = 80;   % User input [%]
ReceivedSignal = ModulatedSignal * (100 - Attenuation)/100 + randn(size(ModulatedSignal)) * NoiseAmplitude/100;

%% Demodulation

MatchedFilter = fliplr(Carrier.');

FilteredSignal = conv(ReceivedSignal, MatchedFilter, 'same');

SampleIdx = floor(CarrierLength / 2);
DecisionSamples = FilteredSignal(SampleIdx : CarrierLength : end);

DemodulatedBinaryData = double(DecisionSamples > 0);

%% Visualization

figure(4);
plot(FilteredSignal(1:CarrierLength*20), 'LineWidth', 1.2, 'Color', 'b');
hold on;
stem(SampleIdx:CarrierLength:CarrierLength*20, DecisionSamples(1:20), 'r', 'LineWidth', 1.5);
title('Matched Filter Output (First 20 Bits)');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Matched Filter Output', 'Decision points', 'Location', 'best');

%% DAC

% Previous unoptimized version
% ReconstructedBinaryData = [];
% for i = 1 : bits : length(DemodulatedBinaryData)
%     DemodulatedBinarySample = DemodulatedBinaryData( i : i+bits-1 );
%     Accumulator = 0;
%     for j = 1:bits
%         bit = DemodulatedBinarySample(j);
%         % LSB first
%         Accumulator = Accumulator + bit * 2 ^ (j-1);
%     end
%     ReconstructedBinaryData = [ReconstructedBinaryData; Accumulator];
% end

% Optimized version
BinaryMatrix = reshape(DemodulatedBinaryData, bits, []);

% BitPowers = 2.^(0:bits-1);
BitPowersT = BitPowers';
ReconstructedSamples = BinaryMatrix.' * BitPowersT;

% NormalizedData = 0.5 * (CompressedData + 1) * (2^bits - 1);
% NormalizedData = round(NormalizedData);
ReconstructedData = 2 * ReconstructedSamples / (2^bits-1) - 1;
% -1 = ReconstructedData <= 1

%% Comparisons

ComparisonBinaryData = BinaryData == DemodulatedBinaryData;
NumberOfTransmittedBits = length(ComparisonBinaryData);
BitErrors = NumberOfTransmittedBits - sum(ComparisonBinaryData);
BitErrorRate = BitErrors / NumberOfTransmittedBits * 100;%

fprintf('Total Transmitted Bits: %d\n', NumberOfTransmittedBits);
fprintf('Bit Errors:             %d\n', BitErrors);
fprintf('Bit Error Rate (BER):   %.4f%%\n', BitErrorRate);

figure(5);
subplot(2,1,1);
plot(CompressedData(1:500), 'Color', 'y', 'LineWidth', 1.2);
hold on;
plot(ReconstructedData(1:500), '--r', 'LineWidth', 1.0);
grid on;
title('Transmitted vs Reconstructed Signal (First 500 Samples)');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Compressed', 'Reconstructed', 'Location', 'northeastoutside');

subplot(2,1,2);
ReconstructionError = CompressedData(1:length(ReconstructedData)) - ReconstructedData;
plot(ReconstructionError(1:500), 'Color', 'r', 'LineWidth', 1.0);
grid on;
title('Reconstruction Error');
xlabel('Sample Index');
ylabel('Error');

%% Play

audioplay = audioplayer( ReconstructedData, FsmpFiltered );
playblocking(audioplay);
