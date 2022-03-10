clear all
fs=2048; % sampling frequency
NumSubj=12;  %number of subject
ch_number=32;  % number of channels

current_folder = pwd;
files = dir(fullfile(current_folder,'\*.set'));
files = {files.name};
files = sort(files);

destination_folder = (pwd);

%parpool local;

eeglab;

subject = 1;
for subject = 1:length(files) 
    filename = cell2mat(files(subject))
% 	open = fullfile('Zahralateral.bdf');
% 	EEG = pop_biosig(open,'channels',1:70);
end


for p=1:NumSubj
    % loading the data
    
%     add=['H:\Post-Doc\DATA\Emotion\all subjects\p1\v' num2str(p) '.mat'];
%     add=['F:\Montreal\Montreal_Data\EEG\Classified\Mat Data\subject' num2str(p) '_Approach.mat'];
%     EEG=load(add);

%%
    % loading the data
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename','Abolfazl_Approach_Dipole.set','filepath','F:\\Montreal\\Montreal_Data\\EEG\\Classified\\all_approach_lateral\\');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG1=ALLEEG.data;
    
%%
    % decomposition into brain frequency bands 
    % NOTE: needs installation of EEGLAB
    [delta,filtwts] = eegfilt(EEG1,fs,.5,4);  %delta (1-4 Hz) >> filtwts: is the filter coefficients
    [theta,filtwts] = eegfilt(EEG1,fs,4,8);  %theta (4-8 Hz)   
    [alpha,filtwts] = eegfilt(EEG1,fs,8,14);  %alpha (8-16 Hz)
    [beta,filtwts] = eegfilt(EEG1,fs,14,30);  %beta (16-32 Hz)
    [gamma,filtwts] = eegfilt(EEG1,fs,32,45);  %gamma (32-45 Hz)
%%
    % ploting time-frequency spectrum
    [insamp1,insfrq1]=inst_frq_amp(delta,fs);  % extraction of frequncy and amplitude of a narrow band signal (delta)
    ptitle='time-frequency spectrum (delta)';
    PlotSpectra(insfrq1,insamp1,ptitle);
    
    [insamp2,insfrq2]=inst_frq_amp(theta,fs);  % extraction of frequncy and amplitude of a narrow band signal (theta)
    ptitle='time-frequency spectrum (theta)';
    PlotSpectra(insfrq2,insamp2,ptitle);
    
    [insamp3,insfrq3]=inst_frq_amp(alpha,fs);  % extraction of frequncy and amplitude of a narrow band signal (alpha)
    ptitle='time-frequency spectrum (alpha)';
    PlotSpectra(insfrq3,insamp3,ptitle);
    
    [insamp4,insfrq4]=inst_frq_amp(beta,fs);  % extraction of frequncy and amplitude of a narrow band signal (beta)
    ptitle='time-frequency spectrum (beta)';
    PlotSpectra(insfrq4,insamp4,ptitle);
    
    ptitle='time-frequency spectrum (all bands)';
    PlotSpectra([insfrq1;insfrq2;insfrq3;insfrq4],[insamp1;insamp2;insamp3;insamp4],ptitle);    % plot all above mentioned bands in a singla figure
    
%%    
    % Phase-phase synchronization matrix
    for ch_number=1:ALLEEG.nbchan
        sig1(ch_number,:)=[delta(ch_number,:)]; 
    end
    pow = sig1.*sig1;
    insphase=inst_phase(pow);
    sync=phasesync(insphase);           % connectivity matrix
    % global efficiency
    GE=efficiency_wei(sync,0);   
    % local efficiency
    LE=efficiency_wei(sync,1);
    
    figure;
    mx=max(LE);mn=min(LE);
    topoplot(LE,EEG.chanlocs,'electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Local efficiency')
%     
    %DFA : scaling exponent
    figure;    
    for i=1:ch %calculation for each channel
        x=beta(i,:);
        [y1,y2]=DFA_main(x);
        D(i)=y1;
        A(i)=y2;
    end
    mx=max(A);mn=min(A);
    figure;
    topoplot(A,'Standard-10-20-Cap32.locs','electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Scaling Exponent')

   %Entropy
    for i=1:ch %calculation for each channel
        x=[delta(i,:);theta(i,:);alpha(i,:);beta(i,:);gamma(i,:)];
        En(i)=FFTentropy(x);
    end
    mx=max(En);mn=min(En);
    figure;
    topoplot(En,'Standard-10-20-Cap32.locs','electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Entropy')
    
    %Fractal Dimention higuchin & petrosian algorithm
    kmax=13;
    for i=1:ch %calculation for each channel
        x=beta(i,:);
        [d,LM]=higuchin(x,kmax);
        higDim(i)=d;
        
        d=petrosian(x);
        petDim(i)=d;
    end

end
