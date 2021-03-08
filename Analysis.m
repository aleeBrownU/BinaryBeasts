clearvars;
close all;

cell_id = 135:269; % all of the cell ID numbers
isi_avg = []; %to put all of the ISI averages per file
isi_std = [];

for input_file_number = 2:12 %Need to run code for each file. This for loop runs the file for 10 times
    number_str = num2str(input_file_number);
    spike_file = readmatrix(strcat('Gamma',number_str,'.txt')); %Loads file (first column is times and second column is ID number
    
    %put code here
    
    isis_spk = []; % column vector for storing the ISI value for each cell to be averaged later
    isi_rand = [];
    
    sample_isi_avg = [];
    sample_isi_std = [];
    
    
    
    % find spikes corresponding to each pyramidal cell
    cell_spikes = cell(length(cell_id),1); % matrix of cells where each item contains the list of spikes for a given cell
    for idx=1:length(cell_id)
        current_id = cell_id(idx); % the current ID
        sel = spike_file(:,2)==current_id; % select indicies within spike_file that correspond to the current ID
        current_times = spike_file(sel,1); % select the times of all spikes that correspond to the current ID
        curr_isi = 1./diff(current_times./1000); % calculate isi values
        
        isis_spk = [isis_spk;curr_isi]; %appends to isis_spk the recently calculated isi
        
        
    end
    
    isi_avg = [isi_avg, mean(isis_spk)];  %compiles all of the mean ISIs
    
    isi_std = [isi_std, std(isis_spk)];
    
    
    
    for num_samples = 1:20 
        rand_id = randi([135 269],1,1000);
        for idx=1:length(cell_id)
            current_id = cell_id(idx); % the current ID
            sel = spike_file(:,2)==current_id; % select indicies within spike_file that correspond to the current ID
            current_times = spike_file(sel,1); % select the times of all spikes that correspond to the current ID
            curr_isi = 1./diff(current_times./1000); % calculate isi values
            
            if any(rand_id(:) == current_id) %asks if any of the random sample ID's
                isi_rand = [isi_rand, curr_isi'];
            end
        end
        sample_isi_avg = [sample_isi_avg, mean(isi_rand)]; %compiles all of the mean ISIs for the random sample
        sample_isi_std = [sample_isi_std, std(isi_rand)]; %compiles all of the standard deviations
        sample_isi_stderr = std(sample_isi_avg)/sqrt(length(cell_id));
    
    end
    
    %start of raster plot
bin = 10;
raster = ones(134, 550);
histo = [];
for row_cell = 1:134 
    neuron_id = 134+row_cell;
    for column_time = 1:bin:550
        time_sel = (spike_file(:,1) > column_time) & (spike_file(:,1) < column_time+bin);
        unit_id = spike_file(time_sel,2);
        x = spike_file(time_sel,1);
        
        if any(unit_id(:) == neuron_id)
            raster(row_cell,column_time:column_time+bin) = 0;
            histo=[histo, column_time];
        end
        
    end
end

figure(1)
subplot(1,11,input_file_number-1);
imagesc(raster)
title(strcat('GABA ',number_str,'10^-2 uS Raster Plot'));
xlabel('Time (ms)')
ylabel('Cell ID')

colormap gray
axis on

figure(2)
subplot(1,11,input_file_number-1);
histogram(histo,'BinWidth',bin);
title(strcat('GABA ',number_str,'10^-2 uS Histogram'));
xlabel('Time(ms)')
ylabel('Number of Cells Spiking')

end
data_pop = [2:12 ; isi_avg]';




figure(3)
scatter(data_pop(:,1),data_pop(:,2));
hold on
gamma_number = data_pop(:,1);
plot_avg_isi = data_pop(:,2);
plot(fitlm(data_pop(:,1),data_pop(:,2)));
title('Plot of Different Gammas and Average ISI');
xlabel('GABAa Conductance (x10^-2 uS)')
ylabel('Frequency of the Neural Circuit (Hz)')


%{
disp('isi avg');
disp(isi_avg);
disp('isi std');
disp(isi_std);
disp('sample isi avg');
disp(sample_isi_avg);
disp('sample isi std');
disp(mean(sample_isi_std));
%}
