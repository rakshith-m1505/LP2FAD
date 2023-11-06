function [close_matches_FAD,close_matches_LP_err]=LP2FAD(LP_targets,D_N,angle_multiples,symmetry,balanced,num_of_sols,print_opt)

% LP2FAD © 2023 by Rakshith Manikandan is licensed under CC BY-SA 4.0. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/


% Input variables (LPs,N,Angles,Design Guidelines,Result print options)
% 
% LP_targets -> Target Inplane LPs Target In-Plane LPs [V1,V2,V3,V4]
% D_N -> Number of Layers to be designed
% angle_multiples -> Desired angle multiples in Layup <45 / 30 / 15>
% symmetry -> Symmetric Design Guideline <true/false>
% balanced -> Balanced Design Guideline <true/false>
% num_of_sols -> Maximum no. of closely matching solutions to be saved
% print_opt -> Print Results <true/false>

%% Necessary Pre-processing of Input

% Total No. of layers (N) in the laminate is double the no. of design layers(D_N) when Symmetric
if symmetry
    N = 2*D_N;
else 
    N = D_N;
end

% Generate list of Signal Patterns to use based on desired angle multiples
if angle_multiples==45
    styles_to_check = 0;
elseif  angle_multiples==30
    styles_to_check = [26,28,33,49];
elseif  angle_multiples == 15
    exact_patterns = [0 20:54];
    approximate_patterns = [2,3,4,5,6,16,17,18,19,188,201,2013,2015,2017,2018,2020,20155,2023,2024,2025,2026,2027,2030,2031];
    styles_to_check = [exact_patterns approximate_patterns];
end

tic;

%% Perform FFT for all LP Samples

% Predefine array to store results
results_SS = zeros(size(styles_to_check,2),N);
results_SS_angs = zeros(size(styles_to_check,2),1);
results_LP_Error = ones(size(styles_to_check,2),1);

% Perform FFT for all relevant Signal Patterns
for rr=1:size(styles_to_check,2)
    [temp_SS, results_LP_Error(rr,1), no_of_eval] = LP2SS_FFT(LP_targets,N,balanced,angle_multiples,styles_to_check(1,rr));
    if ~isempty(temp_SS)
        results_SS(rr,:) = temp_SS;
        results_SS_angs(rr,1) = size(unique(abs(temp_SS)),2);
    else
        results_SS(rr,:) = 2000*ones(size(1,N));
    end
end

% Identify null results
rowToDelete = find(all(results_SS == 1000, 2));
results_SS(rowToDelete, :) = [];
results_LP_Error(rowToDelete, :) = [];

% Identify and store the best sample style
[~, min_idx] = min(results_LP_Error);
b_stylenum = styles_to_check(1,min_idx);

% Remove duplicate results
[results_SS, pt_identifiers] = unique(sort(results_SS,2),'rows');
results_LP_Error = results_LP_Error(pt_identifiers,1);

% Identify and store best result
[~,min_idxs] = sort(results_LP_Error,'ascend');

if num_of_sols>size(results_LP_Error,1)
    num_of_sols = size(results_LP_Error,1);
end
close_matches_FAD = results_SS(min_idxs(1:num_of_sols),:);
close_matches_LP_err = results_LP_Error(min_idxs(1:num_of_sols));

%% Print best result
format short;

if isempty(results_LP_Error)
    fprintf('\n No results were obtained \n');
elseif print_opt
    [b_layup_error, min_idx] = min(results_LP_Error);
    b_layup = results_SS(min_idx,:);
    print_results(b_stylenum,b_layup,b_layup_error,no_of_eval);
    if num_of_sols>1
        print_mult_results(close_matches_FAD,close_matches_LP_err);
    end
end

toc;
%% Functions used

function [best_layup,best_layup_error,no_of_eval] = LP2SS_FFT(LPs,N,balanced,angle_multiples,Sample_Style)

% Number of frequency bins (= Ply orientations)
if balanced
    Fs = 90;
else
    Fs = 180;
end

% % Generate samples and their timestamps using observed signal patterns
[x,t] = signal_patterns(Sample_Style,LPs,Fs);

% % Perform FFT (with theta as frequency)

% Use (nu)fft function and calculate Discrete Fourier Transform Coefficient
y = nufft(x,t);

% Shift values across either side of the y-axis to have 90*2 frequency windows
z = fftshift(y);

% Remove Complex correlations
z = real(z);

% Remove Negative correlations
z(z<0) = 0;

% Scaling FFT Coefficients to fit within the interested frequency window
ly = length(y);
f = (-ly/2:ly/2-1)/ly*(180);

% Remove FFT coefficients of angle counts outside desired set
[~,idx_to_rem] = find(mod(f,angle_multiples)~=0);
z(idx_to_rem)=0;

% Normalise FFT Coefficient with this value and scale them with No. of Plies (rounded to closest 0.5)
no_of_plies = z*(N/sum(z));


off_angles = find(mod(1:89,angle_multiples)==0);
off_angles_neg = -1*off_angles;


% Identify positive off-axis ply locations
[~, indices] = ismembertol(f, off_angles);
off_locs = find(indices);
% 
% % Identify negative off-axis ply locations
[~, indices] = ismembertol(f, off_angles_neg);
off_locs_neg = find(indices);

if (round(10000*sum(no_of_plies))/10000) > N
    error('Fibre Angle Distribution exceeds target Number of Layers upon rounding');
end

% Designing Unbalanced FAD
if (balanced == 0)&&(not(isempty(off_locs(no_of_plies(off_locs)>0))&&((LPs(2,1)~=0)&&(LPs(4,1)~=0))))
    [~, zero_idx] = find(f==0);
    no_of_plies(no_of_plies<1e-5) = 0;
    pm_no_of_plies = no_of_plies(1:zero_idx) + [0 fliplr(no_of_plies(zero_idx+1:end)) 0];
    non_zer_off_locs = off_locs(no_of_plies(off_locs)>0);
    non_zer_off_locs_neg = fliplr(off_locs_neg(no_of_plies(off_locs_neg)>0));
    pm_no_of_plies(pm_no_of_plies<1e-5) = 0;
    num_off_groups = length(find(pm_no_of_plies(2:end-1)>0));

    if num_off_groups==1
        temp1 = no_of_plies(non_zer_off_locs) + ((N*LPs(2,1))/(2*sind(2*f(non_zer_off_locs))));
        temp2 = no_of_plies(non_zer_off_locs_neg) - ((N*LPs(2,1))/(2*sind(2*f(non_zer_off_locs))));
        if (round(temp1)>=0)&&(round(temp2)>=0)
            if temp1<0
                no_of_plies(non_zer_off_locs) = round(temp1);
            else
                no_of_plies(non_zer_off_locs) = temp1;
            end
            if temp2<0
                no_of_plies(non_zer_off_locs_neg) = round(temp2);
            else
                no_of_plies(non_zer_off_locs_neg) = temp2;
            end
        end
    elseif num_off_groups == 2
        unbalancedness_1 = ((N*(LPs(2,1)*sind(4*f(non_zer_off_locs(2))) - LPs(4,1)*sind(2*f(non_zer_off_locs(2)))))/(sind(2*f(non_zer_off_locs(1)))*sind(4*f(non_zer_off_locs(2))) - sind(4*f(non_zer_off_locs(1)))*sind(2*f(non_zer_off_locs(2)))));
        unbalancedness_2 = ((N*(LPs(2,1)*sind(4*f(non_zer_off_locs(1))) - LPs(4,1)*sind(2*f(non_zer_off_locs(1)))))/(sind(2*f(non_zer_off_locs(1)))*sind(4*f(non_zer_off_locs(2))) - sind(4*f(non_zer_off_locs(1)))*sind(2*f(non_zer_off_locs(2)))));
        temp1 = no_of_plies(non_zer_off_locs(1)) + 0.5*unbalancedness_1;
        temp2 = no_of_plies(non_zer_off_locs_neg(1)) - 0.5*unbalancedness_1;
        temp3 = no_of_plies(non_zer_off_locs(2)) - 0.5*unbalancedness_2;
        temp4 = no_of_plies(non_zer_off_locs_neg(2)) + 0.5*unbalancedness_2;

        %         &&(temp1+temp2<=no_of_plies(non_zer_off_locs(1))+no_of_plies(non_zer_off_locs_neg(1)))&&(temp3+temp4<=no_of_plies(non_zer_off_locs(2))+no_of_plies(non_zer_off_locs_neg(2)))
        if (round(temp1)>=0)&&(round(temp2)>=0)&&(round(temp3)>=0)&&(round(temp4)>=0)
            if temp1<0
                no_of_plies(non_zer_off_locs(1)) = round(temp1);
            else
                no_of_plies(non_zer_off_locs(1)) = temp1;
            end
            if temp2<0
                no_of_plies(non_zer_off_locs_neg(1)) = round(temp2);
            else
                no_of_plies(non_zer_off_locs_neg(1)) = temp2;
            end
            if temp3<0
                no_of_plies(non_zer_off_locs(2)) = round(temp3);
            else
                no_of_plies(non_zer_off_locs(2)) = temp3;
            end
            if temp4<0
                no_of_plies(non_zer_off_locs_neg(2)) = round(temp4);
            else
                no_of_plies(non_zer_off_locs_neg(2)) = temp4;
            end
        end
    elseif num_off_groups >= 3
        [min_unb_val,min_unb] = min([no_of_plies(non_zer_off_locs(1))+no_of_plies(non_zer_off_locs_neg(1)) no_of_plies(non_zer_off_locs(2))+no_of_plies(non_zer_off_locs_neg(2)) no_of_plies(non_zer_off_locs(3))+no_of_plies(non_zer_off_locs_neg(3))]);
        org_non_zer_off_locs = non_zer_off_locs;
        org_non_zer_off_locs_neg = non_zer_off_locs_neg;
        non_zer_off_locs(min_unb) = [];
        non_zer_off_locs_neg(min_unb) = [];
        min_unb_ang = f(org_non_zer_off_locs(min_unb));

        X = floor(min_unb_val);  % Number of a +- theta


        % Generate combinations
        if X>0
            Np_combs = 0:X;
            Np_combs(Np_combs==X/2) = [];
            N3_combs = [Np_combs;fliplr(Np_combs)];
        else
            N3_combs = [0 X; X 0];
        end

        combs_to_ignore = [];
        no_of_plies_gcs = repmat(no_of_plies, size(N3_combs,2),1);

        for gc=1:length(N3_combs)
            N3p = N3_combs(1,gc);
            N3m = N3_combs(2,gc);
            no_of_plies_gcs(gc,org_non_zer_off_locs(min_unb)) = N3p;
            no_of_plies_gcs(gc,org_non_zer_off_locs_neg(min_unb)) = N3m;
            unbalancedness_1 = (N3p*(sind(4*f(non_zer_off_locs(2)))*sind(2*min_unb_ang) - sind(2*f(non_zer_off_locs(2)))*sind(4*min_unb_ang)) - N3m*(sind(2*f(non_zer_off_locs(2)))*sind(4*min_unb_ang)- sind(4*f(non_zer_off_locs(2)))*sind(2*min_unb_ang)) + N*(-LPs(2,1)*sind(4*f(non_zer_off_locs(2))) + LPs(4,1)*sind(2*f(non_zer_off_locs(2)))))/(2*((sind(4*f(non_zer_off_locs(1)))*sind(2*f(non_zer_off_locs(2))) - sind(2*f(non_zer_off_locs(1)))*sind(4*f(non_zer_off_locs(2))))));
            unbalancedness_2 = (N3p*(sind(4*f(non_zer_off_locs(1)))*sind(2*min_unb_ang) - sind(2*f(non_zer_off_locs(1)))*sind(4*min_unb_ang)) - N3m*(sind(2*f(non_zer_off_locs(1)))*sind(4*min_unb_ang)- sind(4*f(non_zer_off_locs(1)))*sind(2*min_unb_ang)) + N*(-LPs(2,1)*sind(4*f(non_zer_off_locs(1))) + LPs(4,1)*sind(2*f(non_zer_off_locs(1)))))/(2*((sind(4*f(non_zer_off_locs(1)))*sind(2*f(non_zer_off_locs(2))) - sind(2*f(non_zer_off_locs(1)))*sind(4*f(non_zer_off_locs(2))))));


            temp1 = no_of_plies(non_zer_off_locs(1)) + unbalancedness_1;
            temp2 = no_of_plies(non_zer_off_locs_neg(1)) - unbalancedness_1;
            temp3 = no_of_plies(non_zer_off_locs(2)) - unbalancedness_2;
            temp4 = no_of_plies(non_zer_off_locs_neg(2)) + unbalancedness_2;

            if (round(temp1)>=0)&&(round(temp2)>=0)&&(round(temp3)>=0)&&(round(temp4)>=0)
                if temp1<0
                    no_of_plies_gcs(gc,non_zer_off_locs(1)) = round(temp1);
                else
                    no_of_plies_gcs(gc,non_zer_off_locs(1)) = temp1;
                end
                if temp2<0
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(1)) = round(temp2);
                else
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(1)) = temp2;
                end
                if temp3<0
                    no_of_plies_gcs(gc,non_zer_off_locs(2)) = round(temp3);
                else
                    no_of_plies_gcs(gc,non_zer_off_locs(2)) = temp3;
                end
                if temp4<0
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(2)) = round(temp4);
                else
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(2)) = temp4;
                end
            elseif (abs(temp1)+abs(temp2)>no_of_plies(non_zer_off_locs(1))+no_of_plies(non_zer_off_locs_neg(1)))||(abs(temp3)+abs(temp4)>no_of_plies(non_zer_off_locs(2))+no_of_plies(non_zer_off_locs_neg(2)))
                if temp1<0
                    no_of_plies_gcs(gc,non_zer_off_locs(1)) = 0;
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(1)) = no_of_plies(non_zer_off_locs(1))+no_of_plies(non_zer_off_locs_neg(1));
                else
                    no_of_plies_gcs(gc,non_zer_off_locs(1)) = no_of_plies(non_zer_off_locs(1))+no_of_plies(non_zer_off_locs_neg(1));
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(1)) = 0;
                end

                if temp3<0
                    no_of_plies_gcs(gc,non_zer_off_locs(2)) = 0;
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(2)) = no_of_plies(non_zer_off_locs(2))+no_of_plies(non_zer_off_locs_neg(2));
                else
                    no_of_plies_gcs(gc,non_zer_off_locs(2)) = no_of_plies(non_zer_off_locs(2))+no_of_plies(non_zer_off_locs_neg(2));
                    no_of_plies_gcs(gc,non_zer_off_locs_neg(2)) = 0;
                end
            else
                combs_to_ignore = [combs_to_ignore gc];
            end
        end
        % Delete invalid combinations
        no_of_plies_gcs(combs_to_ignore,:) = [];

        % Calculate LPs of all ply counts in 'no_of_plies_gcs'
        unb_LPs = zeros(size(no_of_plies_gcs,1),4);

        for i=1:size(no_of_plies_gcs,1)
            unb_LPs(i,1) = sum(no_of_plies_gcs(i,:).*cosd(2.*f))/N;
            unb_LPs(i,3) = sum(no_of_plies_gcs(i,:).*cosd(4.*f))/N;
            unb_LPs(i,2) = sum(no_of_plies_gcs(i,:).*sind(2.*f))/N;
            unb_LPs(i,4) = sum(no_of_plies_gcs(i,:).*sind(4.*f))/N;
        end

        % Subtract obtained LPs and the target LPs
        unb_LPs(:,1) = unb_LPs(:,1)-LPs(1);
        unb_LPs(:,3) = unb_LPs(:,3)-LPs(3);
        unb_LPs(:,2) = unb_LPs(:,2)-LPs(2);
        unb_LPs(:,4) = unb_LPs(:,4)-LPs(4);

        % Calculate error in obtained result (L2)
        unb_abs_error = (sqrt((abs(unb_LPs(:,1))).^2 + (abs(unb_LPs(:,2))).^2 + (abs(unb_LPs(:,4))).^2 + (abs(unb_LPs(:,3)).^2)));

        % Find index of ply count with least error
        [~,unb_b_idx] = min(unb_abs_error);
        no_of_plies = no_of_plies_gcs(unb_b_idx,:);
    end
end


% Unbalanced Rules
% Split solution of No. of Plies into arrays of whole and fractional parts
no_of_plies_frac = no_of_plies - fix(no_of_plies);
no_of_plies_whole = no_of_plies - no_of_plies_frac;

if any(gt(no_of_plies_frac,0.99999))
    no_of_plies(gt(no_of_plies_frac,0.99999)) = round(no_of_plies(gt(no_of_plies_frac,0.99999)));
    no_of_plies_frac = no_of_plies - fix(no_of_plies);
    no_of_plies_whole = no_of_plies - no_of_plies_frac;
end


% Calculate number of plies to be obtained from rounding
to_be_rounded = N-sum(no_of_plies_whole);


% Find Number of combinations to evaluate by Rounding

% Identify indices of both whole and fractional non-zero ply counts
all_idx = unique(sort([find(no_of_plies_whole~=0) find(no_of_plies_frac~=0)]));
% Rounding feasibility check
if length(all_idx)<to_be_rounded
    all_poss_num = -1;
else
    % Number of possible combinations by rounding both whole and fractional counts
    all_poss_num = nchoosek(length(all_idx),to_be_rounded);
end

% Error check
if all_poss_num == -1
    error('Inability to round: Less number of entities are present to be rounded to satisfy required Number of Plies');
end

% Generate all possible rounded combinations
all_poss = nchoosek(all_idx,to_be_rounded);

% Generate fibre angle distributions of them and pick the ones with least LP mismatch
[best_layup,best_layup_error,best_spectrum] = round_and_compare(LPs,N,no_of_plies,no_of_plies_whole,all_poss,f,false,to_be_rounded);


% Save number of evaluated solutions
no_of_eval = size(all_poss,1);

% Check if null result
if(isempty(best_layup))
    best_layup = 1000*ones(1,N);
    best_layup_error = 2;
% Deign of Balanced FAD
elseif balanced
    [~, zero_idx] = find(f==0);
    pm_no_of_plies = best_spectrum(1:zero_idx) + [0 fliplr(best_spectrum(zero_idx+1:end)) 0];

    pm_f = abs(f(1:zero_idx));

    pm_no_of_plies_odd = [mod(pm_no_of_plies(1:end),2)];
    pm_no_of_plies_even = pm_no_of_plies - pm_no_of_plies_odd;

    to_be_rounded = round(N - sum(pm_no_of_plies_even));

    if to_be_rounded>0
        % Identify indices of both whole and fractional non-zero ply counts
        angs_to_be_rounded = pm_f(unique(sort([find(pm_no_of_plies~=0)])));
        if all([(size(angs_to_be_rounded,2)<=3),any(ismember([0,45,90],angs_to_be_rounded))])
            angs_to_be_rounded = [0,45,90];
        else
            angs_to_be_rounded = [0 15 30 45 60 75 90];
        end

        % array is all_ply_angs_to_be_rounded
        all_perms = npermk(angs_to_be_rounded,to_be_rounded);
        all_combs = unique(sort(all_perms,2),'rows');

        % Rounding feasibility check
        if length(all_combs)<to_be_rounded
            error('Inability to round');
        end


        % Specify off_axis values to check for odd occurrences
        valuesToCheck = angs_to_be_rounded(mod(angs_to_be_rounded, 90) > 0);
        rows_to_remove = [];

        for i = 1:size(all_combs, 1)
            row = all_combs(i, :);
            counts = zeros(1, numel(valuesToCheck));

            for j = 1:numel(valuesToCheck)
                counts(j) = sum(row == valuesToCheck(j));
            end

            if any(mod(counts, 2) == 1)
                rows_to_remove = [rows_to_remove, i];
            end
        end

        all_combs(rows_to_remove, :) = [];

        pm_poss = zeros(size(all_combs,1),size(pm_f,2));

        for i=1:size(pm_f,2)
            ang_temp = pm_f(1,i);
            ang_idx = pm_f==ang_temp;
            for j=1:size(all_combs, 1)
                row = all_combs(j, :);
                pm_poss(j,ang_idx) = sum(row == ang_temp);
            end
        end

        %          [ -90                                      0
        all_poss = [pm_poss(:,1) pm_poss(:,(2:end-1))/2 pm_poss(:,end) fliplr(pm_poss(:,(2:end-1))/2)];
        no_of_plies_whole = [pm_no_of_plies_even(1) pm_no_of_plies_even((2:end-1))/2 pm_no_of_plies_even(end) fliplr(pm_no_of_plies_even((2:end-1))/2)];

        % Generate fibre angle distributions of them and pick the ones with least LP mismatch
        [best_layup,best_layup_error,~] = round_and_compare(LPs,N,best_spectrum,no_of_plies_whole,all_poss,f,balanced,to_be_rounded);
        % Save number of evaluated solutions
        no_of_eval = no_of_eval + size(all_poss,1);
    end
end

end


%% Other Functions
% Function to Print Results

function [x,t] = signal_patterns(Sample_Style,LPs,Fs)

% % Generate samples and their timestamps using observed signal patterns
if Sample_Style == 0  % [Δ45 degrees] Laminates (Period = 4)

    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 8  % [Δ60 degrees] Laminates (Period = 3)

    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 9  % [Δ36 degrees] Laminates (Period = 5)
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 10  % [Δ30 degrees] Laminates (Period = 6)
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), (LPs(1,1))-(LPs(3,1)), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 11
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), (LPs(3,1))-(LPs(1,1)), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 12
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), -(LPs(3,1))-(LPs(1,1)), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 13
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), (LPs(3,1))+(LPs(1,1)), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 14
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), LPs(1,1), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

elseif Sample_Style == 15
    % Sample pattern of known data
    x = [1, LPs(1,1), LPs(3,1), LPs(3,1), LPs(3,1), LPs(1,1)];
    % Repeat samples as 'period' multiples to fit frequency window
    x = repmat(x,1,(ceil(Fs/(2*length(x)))));
    % Generate corresponding timestamps (every second)
    t = 0:(size(x,2)-1);

else    % [Δ15 degrees] Laminates (Period = 12)
    period = 12;

    % Time Stamps of Number of Plies N in signal
    n = period*(0:(ceil((Fs/2)/period)));

    % Generate timestamps till a 'period' multiple closest to frequency window
    t = 0:((period*(ceil((Fs/2)/period))-1));

    if Sample_Style == 1
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1),LPs(3,1),LPs(1,1)];

        % Timestamps to be removed (3, 4, 5, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+4 (n)-4 (n)+5 (n)-5 (n)+6 (n)-6]));

    elseif Sample_Style == 2
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -1*LPs(3,1), -1*LPs(1,1), -1*LPs(1,1), -1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 3
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), +1*LPs(3,1), +1*LPs(1,1), +1*LPs(1,1), +1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 4
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 1-(LPs(3,1)), 1-(LPs(1,1)), 1-(LPs(1,1)), 1-(LPs(3,1)), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 5
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), LPs(3,1), LPs(3,1), LPs(3,1), LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 6
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 2*LPs(1,1), LPs(1,1), LPs(1,1), 2*LPs(1,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 7
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -2*LPs(3,1), LPs(3,1), 2*LPs(3,1), LPs(3,1), -2*LPs(3,1), LPs(3,1), LPs(3,1)];

        % Timestamps to be removed (4 and their mirrored samples)
        to_remove = unique(sort([(n)+4 (n)-4]));

    elseif Sample_Style == 16
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -1*LPs(3,1), (LPs(1,1)-LPs(3,1)), LPs(1,1), -1*LPs(1,1), LPs(1,1), (LPs(1,1)-LPs(3,1)), -1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];

    elseif Sample_Style == 17
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -1*LPs(3,1), (LPs(3,1)-LPs(1,1)), LPs(1,1), -1*LPs(1,1), LPs(1,1), (LPs(3,1)-LPs(1,1)), -1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];


    elseif Sample_Style == 18
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), LPs(1,1), LPs(3,1), LPs(1,1), -1*LPs(1,1), LPs(1,1), LPs(3,1), LPs(1,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];


    elseif Sample_Style == 19
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), LPs(1,1), -1*LPs(1,1), LPs(1,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (third and fourth)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+4 (n)-4]));

    elseif Sample_Style == 20
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	6.46410161513775-10.1961524227066*LPs(1,1)+4.73205080756888*LPs(3,1)	20.3923048454133-30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1)	35.3205080756888-51.9807621135332*LPs(1,1)+17.6602540378444*LPs(3,1)	41.7846096908265-61.1769145362398*LPs(1,1)+20.3923048454133*LPs(3,1)	35.3205080756888-51.9807621135332*LPs(1,1)+17.6602540378444*LPs(3,1)	20.3923048454133-30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1)	6.46410161513775-10.1961524227066*LPs(1,1)+4.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 21
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	3.73205080756888-6.46410161513775*LPs(1,1)+3.73205080756888*LPs(3,1)	7.46410161513775-12.9282032302755*LPs(1,1)+6.46410161513775*LPs(3,1)	7.46410161513775-13.9282032302755*LPs(1,1)+7.46410161513775*LPs(3,1)	6.46410161513775-12.9282032302755*LPs(1,1)+7.46410161513775*LPs(3,1)	7.46410161513775-13.9282032302755*LPs(1,1)+7.46410161513775*LPs(3,1)	7.46410161513775-12.9282032302755*LPs(1,1)+6.46410161513775*LPs(3,1)	3.73205080756888-6.46410161513775*LPs(1,1)+3.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 22
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1-2.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	-2.73205080756888*LPs(1,1)+3.73205080756888*LPs(3,1)	-3.73205080756888*LPs(1,1)+4.73205080756888*LPs(3,1)	1+5.46410161513775*(LPs(3,1)-LPs(1,1))	-3.73205080756888*LPs(1,1)+4.73205080756888*LPs(3,1)	-2.73205080756888*LPs(1,1)+3.73205080756888*LPs(3,1)	1-2.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 23
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1+2*LPs(3,1)	-2+LPs(3,1)	-2-LPs(1,1)+4*LPs(3,1)	-3+4*LPs(3,1)	-2-LPs(1,1)+4*LPs(3,1)	-2+LPs(3,1)	-1+2*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 24
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1.73205080756888+LPs(1,1)+1.73205080756888*LPs(3,1)	-2+LPs(3,1)	-3.46410161513775+LPs(1,1)+3.46410161513775*LPs(3,1)	-3+4*LPs(3,1)	-3.46410161513775+LPs(1,1)+3.46410161513775*LPs(3,1)	-2+LPs(3,1)	-1.73205080756888+LPs(1,1)+1.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 25
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	3-5*LPs(1,1)+3*LPs(3,1)	4-6*LPs(1,1)+3*LPs(3,1)	LPs(1,1)	-3+6*LPs(1,1)-2*LPs(3,1)	LPs(1,1)	4-6*LPs(1,1)+3*LPs(3,1)	3-5*LPs(1,1)+3*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 26
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1-2*LPs(1,1)+2*LPs(3,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	1-2*LPs(1,1)+2*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 27
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-0.464101615137755+0.196152422706632*LPs(1,1)+1.26794919243112*LPs(3,1)	-0.392304845413264+0.588457268119896*LPs(1,1)+0.803847577293368*LPs(3,1)	0.679491924311227-0.0192378864668406*LPs(1,1)+0.339745962155614*LPs(3,1)	0.215390309173473+1.17691453623979*LPs(1,1)-0.392304845413264*LPs(3,1)	0.679491924311227-0.0192378864668406*LPs(1,1)+0.339745962155614*LPs(3,1)	-0.392304845413264+0.588457268119896*LPs(1,1)+0.803847577293368*LPs(3,1)	-0.464101615137755+0.196152422706632*LPs(1,1)+1.26794919243112*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 28
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1+LPs(1,1)+LPs(3,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	-1+LPs(1,1)+LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 29
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1-LPs(1,1)+LPs(3,1)	2*LPs(1,1)-LPs(3,1)	LPs(1,1)	1-2*LPs(1,1)+2*LPs(3,1)	LPs(1,1)	2*LPs(1,1)-LPs(3,1)	1-LPs(1,1)+LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 30
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0.267949192431123+0.464101615137755*LPs(1,1)+0.267949192431123*LPs(3,1)	0.535898384862245+0.928203230275509*LPs(1,1)-0.464101615137755*LPs(3,1)	0.535898384862245-0.0717967697244908*LPs(1,1)+0.535898384862245*LPs(3,1)	-0.464101615137755+0.928203230275509*LPs(1,1)+0.535898384862245*LPs(3,1)	0.535898384862245-0.0717967697244908*LPs(1,1)+0.535898384862245*LPs(3,1)	0.535898384862245+0.928203230275509*LPs(1,1)-0.464101615137755*LPs(3,1)	0.267949192431123+0.464101615137755*LPs(1,1)+0.267949192431123*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 31
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 32
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1+0.732050807568877*LPs(1,1)-0.732050807568877*LPs(3,1)	0.732050807568877*LPs(1,1)+0.267949192431123*LPs(3,1)	-0.267949192431123*LPs(1,1)+1.26794919243112*LPs(3,1)	1+1.46410161513775*LPs(1,1)-1.46410161513775*LPs(3,1)	-0.267949192431123*LPs(1,1)+1.26794919243112*LPs(3,1)	0.732050807568877*LPs(1,1)+0.267949192431123*LPs(3,1)	1+0.732050807568877*LPs(1,1)-0.732050807568877*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 33
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1+LPs(1,1)-LPs(3,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	1+LPs(1,1)-LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 34
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1.73205080756888+LPs(1,1)-1.73205080756888*LPs(3,1)	-2+LPs(3,1)	3.46410161513775+LPs(1,1)-3.46410161513775*LPs(3,1)	-3+4*LPs(3,1)	3.46410161513775+LPs(1,1)-3.46410161513775*LPs(3,1)	-2+LPs(3,1)	1.73205080756888+LPs(1,1)-1.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 35
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	2.73205080756888-4.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	2.73205080756888-4.73205080756888*LPs(1,1)+1.73205080756888*LPs(3,1)	-2.73205080756888+3.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	-6.46410161513775+9.46410161513776*LPs(1,1)-5.46410161513775*LPs(3,1)	-2.73205080756888+3.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	2.73205080756888-4.73205080756888*LPs(1,1)+1.73205080756888*LPs(3,1)	2.73205080756888-4.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 36
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0.866025403784439-2*LPs(1,1)+1.73205080756888*LPs(3,1)	-0.5	-0.866025403784439+LPs(1,1)-1.73205080756888*LPs(3,1)	-2*LPs(3,1)	-0.866025403784439+LPs(1,1)-1.73205080756888*LPs(3,1)	-0.5	0.866025403784439-2*LPs(1,1)+1.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 37
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-0.5+LPs(3,1)	-0.5	0.5-LPs(1,1)-LPs(3,1)	-2*LPs(3,1)	0.5-LPs(1,1)-LPs(3,1)	-0.5	-0.5+LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 38
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1+0.732050807568877*LPs(1,1)+0.732050807568877*LPs(3,1)	-0.732050807568877*LPs(1,1)+0.267949192431123*LPs(3,1)	-0.267949192431123*LPs(1,1)-1.26794919243112*LPs(3,1)	1-1.46410161513775*LPs(1,1)-1.46410161513775*LPs(3,1)	-0.267949192431123*LPs(1,1)-1.26794919243112*LPs(3,1)	-0.732050807568877*LPs(1,1)+0.267949192431123*LPs(3,1)	-1+0.732050807568877*LPs(1,1)+0.732050807568877*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 39
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0.732050807568877-1.26794919243112*LPs(1,1)+0.732050807568877*LPs(3,1)	-0.732050807568877+1.26794919243112*LPs(1,1)-1.73205080756888*LPs(3,1)	-0.732050807568877+0.267949192431123*LPs(1,1)-0.732050807568877*LPs(3,1)	0.464101615137755-2.53589838486225*LPs(1,1)+1.46410161513775*LPs(3,1)	-0.732050807568877+0.267949192431123*LPs(1,1)-0.732050807568877*LPs(3,1)	-0.732050807568877+1.26794919243112*LPs(1,1)-1.73205080756888*LPs(3,1)	0.732050807568877-1.26794919243112*LPs(1,1)+0.732050807568877*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 40
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0	-LPs(3,1)	-LPs(1,1)	-1	-LPs(1,1)	-LPs(3,1)	0	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 41
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-0.267949192431123+0.464101615137755*LPs(1,1)-0.267949192431123*LPs(3,1)	0.535898384862245-0.928203230275509*LPs(1,1)-0.464101615137755*LPs(3,1)	-0.535898384862245-0.0717967697244908*LPs(1,1)-0.535898384862245*LPs(3,1)	-0.464101615137755-0.928203230275509*LPs(1,1)+0.535898384862245*LPs(3,1)	-0.535898384862245-0.0717967697244908*LPs(1,1)-0.535898384862245*LPs(3,1)	0.535898384862245-0.928203230275509*LPs(1,1)-0.464101615137755*LPs(3,1)	-0.267949192431123+0.464101615137755*LPs(1,1)-0.267949192431123*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 42
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0.5-LPs(3,1)	-0.5	0.5-LPs(1,1)+LPs(3,1)	-2*LPs(3,1)	0.5-LPs(1,1)+LPs(3,1)	-0.5	0.5-LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 43
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	0.464101615137755+0.196152422706632*LPs(1,1)-1.26794919243112*LPs(3,1)	-0.392304845413264-0.588457268119896*LPs(1,1)+0.803847577293368*LPs(3,1)	-0.679491924311227-0.0192378864668406*LPs(1,1)-0.339745962155614*LPs(3,1)	0.215390309173473-1.17691453623979*LPs(1,1)-0.392304845413264*LPs(3,1)	-0.679491924311227-0.0192378864668406*LPs(1,1)-0.339745962155614*LPs(3,1)	-0.392304845413264-0.588457268119896*LPs(1,1)+0.803847577293368*LPs(3,1)	0.464101615137755+0.196152422706632*LPs(1,1)-1.26794919243112*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 44
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	1-2*LPs(3,1)	-2+LPs(3,1)	2-LPs(1,1)-4*LPs(3,1)	-3+4*LPs(3,1)	2-LPs(1,1)-4*LPs(3,1)	-2+LPs(3,1)	1-2*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 45
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-2*LPs(1,1)	-2+LPs(3,1)	LPs(1,1)	3+4*LPs(3,1)	LPs(1,1)	-2+LPs(3,1)	-2*LPs(1,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 46
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-0.732050807568877-1.26794919243112*LPs(1,1)-0.732050807568877*LPs(3,1)	-0.732050807568877-1.26794919243112*LPs(1,1)-1.73205080756888*LPs(3,1)	0.732050807568877+0.267949192431123*LPs(1,1)+0.732050807568877*LPs(3,1)	0.464101615137755+2.53589838486225*LPs(1,1)+1.46410161513775*LPs(3,1)	0.732050807568877+0.267949192431123*LPs(1,1)+0.732050807568877*LPs(3,1)	-0.732050807568877-1.26794919243112*LPs(1,1)-1.73205080756888*LPs(3,1)	-0.732050807568877-1.26794919243112*LPs(1,1)-0.732050807568877*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 47
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1-LPs(1,1)-LPs(3,1)	-2*LPs(1,1)-LPs(3,1)	LPs(1,1)	1+2*LPs(1,1)+2*LPs(3,1)	LPs(1,1)	-2*LPs(1,1)-LPs(3,1)	-1-LPs(1,1)-LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 48
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-0.866025403784439-2*LPs(1,1)-1.73205080756888*LPs(3,1)	-0.5	0.866025403784439+LPs(1,1)+1.73205080756888*LPs(3,1)	-2*LPs(3,1)	0.866025403784439+LPs(1,1)+1.73205080756888*LPs(3,1)	-0.5	-0.866025403784439-2*LPs(1,1)-1.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 49
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1-2*LPs(1,1)-2*LPs(3,1)	LPs(3,1)	LPs(1,1)	1	LPs(1,1)	LPs(3,1)	-1-2*LPs(1,1)-2*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 50
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-1-2.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	2.73205080756888*LPs(1,1)+3.73205080756888*LPs(3,1)	-3.73205080756888*LPs(1,1)-4.73205080756888*LPs(3,1)	1+5.46410161513775*(LPs(3,1)+LPs(1,1))	-3.73205080756888*LPs(1,1)-4.73205080756888*LPs(3,1)	2.73205080756888*LPs(1,1)+3.73205080756888*LPs(3,1)	-1-2.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 51
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-2.73205080756888-4.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	2.73205080756888+4.73205080756888*LPs(1,1)+1.73205080756888*LPs(3,1)	2.73205080756888+3.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	-6.46410161513775-9.46410161513776*LPs(1,1)-5.46410161513775*LPs(3,1)	2.73205080756888+3.73205080756888*LPs(1,1)+2.73205080756888*LPs(3,1)	2.73205080756888+4.73205080756888*LPs(1,1)+1.73205080756888*LPs(3,1)	-2.73205080756888-4.73205080756888*LPs(1,1)-2.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 52
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-3-5*LPs(1,1)-3*LPs(3,1)	4+6*LPs(1,1)+3*LPs(3,1)	LPs(1,1)	-3-6*LPs(1,1)-2*LPs(3,1) LPs(1,1) 4+6*LPs(1,1)+3*LPs(3,1)	-3-5*LPs(1,1)-3*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 53
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-3.73205080756888-6.46410161513775*LPs(1,1)-3.73205080756888*LPs(3,1)	7.46410161513775+12.9282032302755*LPs(1,1)+6.46410161513775*LPs(3,1)	-7.46410161513775-13.9282032302755*LPs(1,1)-7.46410161513775*LPs(3,1)	6.46410161513775+12.9282032302755*LPs(1,1)+7.46410161513775*LPs(3,1)	-7.46410161513775-13.9282032302755*LPs(1,1)-7.46410161513775*LPs(3,1)	7.46410161513775+12.9282032302755*LPs(1,1)+6.46410161513775*LPs(3,1)	-3.73205080756888-6.46410161513775*LPs(1,1)-3.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 54
        % Sample pattern of known data
        x = [1	LPs(1,1)	LPs(3,1)	-6.46410161513775-10.1961524227066*LPs(1,1)-4.73205080756888*LPs(3,1)	20.3923048454133+30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1)	-35.3205080756888-51.9807621135332*LPs(1,1)-17.6602540378444*LPs(3,1)	41.7846096908265+61.1769145362398*LPs(1,1)+20.3923048454133*LPs(3,1)	-35.3205080756888-51.9807621135332*LPs(1,1)-17.6602540378444*LPs(3,1)	20.3923048454133+30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1)	-6.46410161513775-10.1961524227066*LPs(1,1)-4.73205080756888*LPs(3,1)	LPs(3,1)	LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];
    elseif Sample_Style == 188
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), LPs(1,1), LPs(3,1), LPs(1,1), 1*LPs(1,1), LPs(1,1), LPs(3,1), LPs(1,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove = [];

    elseif Sample_Style == 199
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -0.5*LPs(3,1), -1*LPs(1,1), -1*LPs(1,1), -0.5*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3, 6 and their mirrored samples)
        to_remove = unique(sort([(n)+3 (n)-3 (n)+6 (n)-6]));

    elseif Sample_Style == 201
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 0.5*LPs(3,1), -1*LPs(3,1), -0.85*LPs(3,1), 1*LPs(1,1), -0.85*LPs(3,1), -1*LPs(3,1), 0.5*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3 and their mirrored samples)
        to_remove = [];


    elseif Sample_Style == 2011
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -1*LPs(3,1), 0.85*LPs(3,1), 1*LPs(1,1), 0.85*LPs(3,1), -1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3 and their mirrored samples)
        to_remove =   unique(sort([(n)+3 (n)-3]));


    elseif Sample_Style == 2012
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), -1*LPs(3,1), -1*0.85*LPs(3,1), -1*LPs(1,1), -1*0.85*LPs(3,1), -1*LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (3 and their mirrored samples)
        to_remove =   unique(sort([(n)+3 (n)-3]));


    elseif Sample_Style == 2013
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), LPs(3,1), 0.5*LPs(3,1), 0.5*(LPs(1,1)+LPs(3,1)), 0.5*LPs(3,1), 0.5*(LPs(1,1)+LPs(3,1)), 0.5*LPs(3,1), LPs(3,1), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove =   [];


    elseif Sample_Style == 2014
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 0.366*(LPs(3,1)+LPs(1,1)), 1.66*LPs(1,1), LPs(1,1), 0.33*LPs(1,1), LPs(1,1), 1.66*LPs(1,1),0.366*(LPs(3,1)+LPs(1,1)), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove =   [];


    elseif Sample_Style == 2015
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 0.366*(LPs(3,1)-LPs(1,1)), 0.34*LPs(1,1), LPs(1,1), 0.33*LPs(1,1), LPs(1,1), 0.34*LPs(1,1),0.366*(LPs(3,1)-LPs(1,1)), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 20155
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 0.366*(LPs(3,1)-LPs(1,1)), -0.34*LPs(1,1), -LPs(1,1), -0.6, -LPs(1,1), -0.34*LPs(1,1),0.366*(LPs(3,1)-LPs(1,1)), LPs(3,1), LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove =   [];


    elseif Sample_Style == 2016
        % Sample pattern of known data
        x = [1, LPs(1,1), LPs(3,1), 0.99*(LPs(3,1)+LPs(1,1)), 0.98*(LPs(1,1)+LPs(3,1))+LPs(1,1),0.69*(LPs(3,1)-LPs(1,1)), 0.99*LPs(1,1), 0.69*(LPs(3,1)-LPs(1,1)),0.98*(LPs(1,1)+LPs(3,1))+LPs(1,1),0.99*(LPs(3,1)+LPs(1,1)),LPs(3,1),LPs(1,1)];

        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2017
        % Sample pattern of known data
        x= [1	LPs(1,1)	LPs(3,1)	LPs(1,1)	20.3923048454133+30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1)	-35.3205080756888-51.9807621135332*LPs(1,1)-17.6602540378444*LPs(3,1)	41.7846096908265+61.1769145362398*LPs(1,1)+20.3923048454133*LPs(3,1) -35.3205080756888-51.9807621135332*LPs(1,1)-17.6602540378444*LPs(3,1) 20.3923048454133+30.5884572681199*LPs(1,1)+11.1961524227066*LPs(3,1) -6.46410161513775-10.1961524227066*LPs(1,1)-4.73205080756888*LPs(3,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2018
        % Sample pattern of known data
        x= [1	LPs(1,1)	LPs(3,1)	0.34*LPs(3,1)+0.84*LPs(1,1) 0.594*LPs(3,1)	0.594*LPs(3,1)+0.87*LPs(1,1) 0.4375*LPs(3,1) 0.594*LPs(3,1)+0.87*LPs(1,1) 0.594*LPs(3,1) 0.34*LPs(3,1)+0.84*LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2020
        % Sample pattern of known data
        x= [1	LPs(1,1)	LPs(3,1)	-2*LPs(1,1) -0.5*(LPs(1,1)+LPs(3,1)) LPs(1,1) 0.5*(LPs(3,1)-LPs(1,1)) LPs(1,1) -0.5*(LPs(1,1)+LPs(3,1)) -2*LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2023
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) LPs(1,1)+LPs(3,1) -1*LPs(1,1) LPs(1,1) -2*LPs(1,1) LPs(1,1) -1*LPs(1,1) LPs(1,1)+LPs(3,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];
    
    elseif Sample_Style == 2024
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) LPs(1,1) 0.5*(LPs(1,1)+LPs(3,1)) LPs(1,1) -0.5*(LPs(1,1)+LPs(3,1)) LPs(1,1) 0.5*(LPs(1,1)+LPs(3,1)) LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];
        
    elseif Sample_Style == 2025
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) -2*LPs(1,1) 0.2*LPs(3,1) LPs(1,1) -0.8*LPs(3,1) LPs(1,1) 0.2*LPs(3,1) -2*LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2026
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) -2*LPs(1,1) 5*LPs(3,1) LPs(1,1) 4*LPs(3,1) LPs(1,1) 5*LPs(3,1) -2*LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];
    
    elseif Sample_Style == 2027
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) 0 1.66*LPs(3,1) -1*LPs(1,1) 0 1*LPs(1,1) 1.66*LPs(3,1) -2*LPs(1,1) LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];
        
    elseif Sample_Style == 2030
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) LPs(3,1)-0.5 -0.5 0.5-(LPs(1,1)+LPs(3,1)) -2*LPs(3,1) 0.5-(LPs(1,1)+LPs(3,1)) -0.5 LPs(3,1)-0.5 LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];

    elseif Sample_Style == 2031
        % Sample pattern of known data
        x= [1	LPs(1,1) LPs(3,1) -LPs(3,1)/2.33333333 1.00669*(LPs(1,1)-LPs(3,1)) (LPs(3,1)/2.33333333)-LPs(1,1) -2*LPs(3,1) (LPs(3,1)/2.33333333)-LPs(1,1) 1.00669*(LPs(1,1)-LPs(3,1)) -LPs(3,1)/2.33333333 LPs(3,1) LPs(1,1)];
        % Timestamps to be removed (none)
        to_remove =   [];
    
    else
        error('Sample Style %d is not defined.',Sample_Style);
    end
    % Remove timestamps of unsampled points
    idx_to_remove = ismember(t,to_remove);
    t(idx_to_remove) = [];

    % Repeat sample pattern to the appropriate length
    x = repmat(x,1,(length(t)/length(x)));
end

end

function print_results(scheme,best_layup,best_layup_error,no_of_eval)
unique_angs = [unique(best_layup)];
fprintf('\nSignal Pattern-%d yielded the closest matching Fibre-Angle Distribution:\n',scheme);
fprintf('\nAngle: ');
disp(unique_angs);
fprintf('Count: ');
disp([histc(sort(best_layup),unique_angs)]);
fprintf('\nLP mismatch error = %f \n',best_layup_error);
fprintf('\nTotal number of layups evaluated = %d \n',no_of_eval);
end

function print_mult_results(close_matches_FAD,close_matches_LP_err)
fprintf('_________________________________________________________________________________\n');
fprintf('_________________________________________________________________________________\n');
fprintf('\nOther closest matching Fibre-Angle Distributions:\n');
for i=2:size(close_matches_FAD,1)
    fprintf('_________________________________________________________________________________\n');
    fprintf('\nAngle: ');
    disp([unique(close_matches_FAD(i,:))]);
    fprintf('Count: ');
    disp([histc(sort(close_matches_FAD(i,:)), unique(close_matches_FAD(i,:)))]);
    fprintf('LP mismatch error = %f \n',close_matches_LP_err(i,1));
end
fprintf('_________________________________________________________________________________\n');
end


% Rounding
function [best_layup,best_layup_error,best_spectrum] = round_and_compare(LPs,N,~,no_of_plies_whole,poss,f,balanced,to_be_rounded)

% Predefine a matrix with each row corresponding to a different ply count
rounded = repmat(no_of_plies_whole, size(poss,1),1);

% Predefine a zero matrix with each row corresponding to rounding info (of 1s)
temp = zeros(size(poss,1),size(no_of_plies_whole,2));

% Based on rounding combination in 'poss', define 1s in temp and add it to frac_rounded to have a rounded ply count
if balanced
    if to_be_rounded>0
        if nnz(no_of_plies_whole)||size(poss,1)~=0
            rounded = rounded + poss;
        end
    end
else
    if to_be_rounded>0
        if nnz(no_of_plies_whole)||size(poss,1)~=0
            for i=1:size(poss,1)
                temp(i,poss(i,:))=1;
                rounded(i,:) = rounded(i,:) + temp(i,:);
            end
        end
    end
end

% Calculate LPs of all ply counts in 'rounded'
rounded_LPs = zeros(size(rounded,1),4);

for i=1:size(rounded,1)
    rounded_LPs(i,1) = sum(rounded(i,:).*cosd(2.*f))/N;
    rounded_LPs(i,3) = sum(rounded(i,:).*cosd(4.*f))/N;
    rounded_LPs(i,2) = sum(rounded(i,:).*sind(2.*f))/N;
    rounded_LPs(i,4) = sum(rounded(i,:).*sind(4.*f))/N;
end

% Subtract obtained LPs and the target LPs
rounded_LPs(:,1) = rounded_LPs(:,1)-LPs(1);
rounded_LPs(:,3) = rounded_LPs(:,3)-LPs(3);
rounded_LPs(:,2) = rounded_LPs(:,2)-LPs(2);
rounded_LPs(:,4) = rounded_LPs(:,4)-LPs(4);

% Calculate error in obtained result (L2)
if (LPs(2)==0) && (LPs(4)==0)
    rounded_abs_error = (sqrt((abs(rounded_LPs(:,1))).^2 + (abs(rounded_LPs(:,3))).^2));
else
    rounded_abs_error = sqrt(((abs(rounded_LPs(:,1))).^2 + (abs(rounded_LPs(:,2))).^2 + (abs(rounded_LPs(:,4))).^2 + (abs(rounded_LPs(:,3)).^2)));
end
% Find index of ply count with least error
[best_layup_error,midx] = min(rounded_abs_error);

best_spectrum = rounded(midx,:);

% Save the most closest ply count
% Check if null result
if isempty(rounded(midx,:))
    best_layup = 1000*ones(1,N);
    best_layup_error = 2;
else
    best = [f; rounded(midx,:)];
    % Save the ply count as a SS
    best_layup=[];
    for i=1:size(f,2)
        for j=1:best(2,i)
            best_layup = [best_layup best(1,i)];
        end
    end
end


end



function [permutedMatrix, index] = npermk(N, K)

% This function is based on the work of Matt Fig's npermutek,
% http://www.mathworks.ch/matlabcentral/fileexchange/11462-npermutek


assert(nargin == 2, 'NPERMK:nrargin', 'NPERMK requires two arguments.');

assert(numel(K) == 1 ...
    && floor(K) == K ...
    && gt(K, 0) ...
    && isreal(K), ...
    'NPERMK:K', ...
    'Second argument K needs to be real positive and scalar integer. See help.');

% get input matrix
Input = N(:);
nrElements          = numel(N);
nrRows              = nrElements ^ K;
inputIndices        = (1:nrElements).';
index               = reshape((1:nrRows * K).', nrRows, K); % allocate space for output

% Create permutation matrix
for iColumn = 1 : K
    row                     = repmat(inputIndices(:,ones(1, nrElements^(iColumn - 1))).', 1, nrElements ^ (K - iColumn));
    index(:, iColumn) = row(:);
end
% flip to common format
index = fliplr(index);
% assign output
permutedMatrix = Input(index);

end

end