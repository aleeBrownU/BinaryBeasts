%Project 1 by Binary Beasts (Adrian Lee, Kush Patel, Jack Waters)

clearvars;
close all;

cell_id = 135:269; % all of the cell ID numbers
pop_freq = []; %to put all of the frequencies of the population from all the different GABA conductances
boot_freq = []; %to put all of the frequences of the bootstrap samples

%Asks user to input bins, bootstrap sample number, and size of bootstrap
%samples
b_input = 0; s_input = 0; p_input = 0;
while b_input == 0 || s_input ==0 || p_input == 0 %to retrieve user input and check if they are within range
bin = input('How wide do you want the histrogram and Raster plot bins to be? (1-550) \n');
if bin >=1 && bin <=550
    b_input = 1;
end
samples = input('How many samples do you want for the boostrap? (>1)\n');
if samples >=1
    s_input = 1;
end
pop = input('How big do you want the bootstrap samples to be? (<67)\n');
if pop >=1 && bin <=67
    p_input = 1;
end
if b_input == 0 || s_input ==0 || p_input == 0 %if the user input is not within the expected range of numbers, the loop will rerun
    disp('One of your inputs was not within the range of expected numbers. Please try again.')
end
end

for input_file_number = 2:12 %Need to run code for each file. This for loop runs the file for 11 times
    number_str = num2str(input_file_number);
    spike_file = readmatrix(strcat('Gamma',number_str,'.txt')); %Loads file (first column is times and second column is ID number)
    
    pop_cell_freq = []; % column vector for storing the frequency for each cell to be averaged later
    boot_rand = []; % initialize vector to be used later in bootstrapping 
    
    sample_pop_freq = []; %To be used later for bootstrap frequencies of all samples for a given Gamma simulation
    
    % find spikes corresponding to each cell
    cell_spikes = cell(length(cell_id),1); % matrix of cells where each item contains the list of spikes for a given cell
    for idx=1:length(cell_id)
        current_id = cell_id(idx); % the current ID
        sel = spike_file(:,2)==current_id; % select indicies within spike_file that correspond to the current ID
        current_times = spike_file(sel,1); % select the times of all spikes that correspond to the current ID
        curr_freq = 1./diff(current_times./1000); % calculate ISI values and use them to find the frequency (1/ISI)
        
        pop_cell_freq = [pop_cell_freq;curr_freq]; %appends to pop_cell_freq the recently calculated frequency of the cell
        
        
    end
    
    pop_freq = [pop_freq, mean(pop_cell_freq)];  %compiles all of the mean frequencies of the population

    for num_samples = 1:samples  %Loops for the number of bootstrap samples user requested
        rand_id = randi([135 269],1,pop); %takes a random sample of cell IDs of size pop that was earlier specified by the user
        for idx=1:length(cell_id)
            current_id = cell_id(idx); % the current ID
            sel = spike_file(:,2)==current_id; % select indicies within spike_file that correspond to the current ID
            current_times = spike_file(sel,1); % select the times of all spikes that correspond to the current ID
            curr_freq = 1./diff(current_times./1000); % calculate ISI values and use them to find the frequency (1/ISI)
            
            if any(rand_id(:) == current_id) %asks if any of the random sample ID's are equal to the current cell ID that's being looped
                boot_rand = [boot_rand, curr_freq']; %If so, it will append the frequency to boot_rand
            end
        end
        sample_pop_freq = [sample_pop_freq, mean(boot_rand)]; %compiles all of the mean frequencies for the random sample
    end
    boot_freq = [boot_freq, mean(sample_pop_freq)]; %compiles all of the mean boot frequencies for all of the Gamma simulations
    
    %start of raster plot
raster = ones(134, 550); %Creates matrix of all ones to be used for the raster plot. 134 = number of cells, 550 = time length of the Gamma simulation
histo = []; %to be used later for histogram
for row_cell = 1:134 %will loop through every cell ID
    neuron_id = 134+row_cell; %Row_cell one corresponds to the Cell ID of 135 so we need to express this relationship
    for column_time = 1:bin:550 %will break down the total time into bins that were earlier specified by the user and loop through them
        time_sel = (spike_file(:,1) > column_time) & (spike_file(:,1) < column_time+bin); %specifies time range for the current bin
        unit_id = spike_file(time_sel,2); %selects all IDS that spike within the bin
        
        if any(unit_id(:) == neuron_id) %Check if the current ID matches any of the cell IDs that spike within the current bin
            raster(row_cell,column_time:column_time+bin) = 0; %if so, it will change the value to 0, changing the color to black in the raster plot
            histo=[histo, column_time]; %appends the time to the vector 'histo' to be used later for the histogram 
        end
        
    end
end

%figure 1 will be the set of Raster plots for each Gamma simulation
figure(1)
subplot(1,11,input_file_number-1); %Puts the Raster plot in the appropriate subplot
imagesc(raster) %displays the Rasterplot as an image
title(strcat('GABA ',number_str,'10^-2 uS Raster Plot')); %title and axis labels are defined
xlabel('Time (ms)')
ylabel('Cell ID')

colormap gray %makes it grayscale
axis on %ensures the lines are there

%figure 2 will be the set of Histograms that compile the number of cells
%that spike within a bin (seen in the raster plot) of the Gamma simulation
figure(2)
subplot(1,11,input_file_number-1);%Puts the Raster plot in the appropriate subplot
histogram(histo,'BinWidth',bin); %displays the histogram
title(strcat('GABA ',number_str,'10^-2 uS Histogram')); %appropriate title and axis labels are given
xlabel('Time(ms)')
ylabel('Number of Cells Spiking')
%}
end
data_pop = [2:12 ; pop_freq]'; %creates a matrix of gaba conductance and corresponding population frequencies averages
boot_pop = [2:12 ; boot_freq]'; %creates a matrix of gaba conductance and corresponding bootstrapped samples averages

[R,P] = corrcoef(pop_freq,boot_freq); %calculates correlation value and p value comparing pop_freq and boot_freq

if P(1,2) > 0.05 %shows different messages if the P value is greater or below 0.05
    disp('The relationship between the population frequencies and bootstrap frequencies is explainable from random chance')
    disp(strcat('The P value is: ', num2str(P(1,2))));
    disp(strcat('The correlation (R) value is: ', num2str(R(1,2))));
else
    disp('The relationship between the population frequencies and bootstrap frequencies is not explainable from random chance.')
    disp(strcat('The P value is: ', num2str(P(1,2))));
    disp(strcat('The correlation (R) value is: ', num2str(R(1,2))));
end

%figure 3 is the plot of the average frequency of the population versus the GABA conductance
figure(3)
scatter(data_pop(:,1),data_pop(:,2)); %plots the scatterplot
hold on
gamma_number = data_pop(:,1);
plot_avg_isi = data_pop(:,2);
plot(fitlm(data_pop(:,1),data_pop(:,2))); %fits a linear regression line to it
title('Plot of Different Gammas and Average ISI'); %appropriate title and labels are given
xlabel('GABAa Conductance (x10^-2 uS)')
ylabel('Frequency of the Neural Circuit (Hz)')

%figure 4 is the plot of the average frequency of the bootstrap samples versus the GABA conductance
figure(4)
scatter(boot_pop(:,1),boot_pop(:,2)); %plots the scatterplot
hold on
gamma_number = boot_pop(:,1);
plot_avg_isi = boot_pop(:,2);
plot(fitlm(boot_pop(:,1),boot_pop(:,2))); %fits a linear regression line to it
title('Plot of Different Gammas and Average ISI using Bootstrapping Data'); %appropriate title and labels are given
xlabel('GABAa Conductance (x10^-2 uS)')
ylabel('Frequency of the Neural Circuit (Hz)')


T = table([2:12]',pop_freq',boot_freq'); %creates table to display average frequencies
%below line gives the titles of the columns for the table
T.Properties.VariableNames = {'GABAa Conductance (10^-2 uS)', 'Average Frequency of the Population', 'Average Frequency of the Bootstrap Samples'};
disp(T(:,:)); %displays the table

