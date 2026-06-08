%% =========================================================================
%  Generate stimulus information for each trial in one run
%
%  Purpose:
%  This script reads Expo XML files from one recording session/run and
%  extracts stimulus parameters trial-by-trial. The resulting "stiminfo"
%  structure is saved into the current CatGT folder and can later be used
%  to align neuronal data (e.g. spike_unit_time_trial).
%
%  Main output:
%      stiminfo.mat
%
%  Notes:
%  1. The script assumes the current folder is a CatGT run folder.
%  2. It also assumes the corresponding Expo XML files are stored in the
%     matching expo_data directory.
%  3. Different stimulus blocks are parsed by different helper functions.
% =========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

editorObj = matlab.desktop.editor.getActive;
cd(fileparts(editorObj.Filename));

% CatGT folder shared by this run, stiminfo saved loaction,
catgt_folder=cd;

[~, catgt_name] = fileparts(catgt_folder);
run_g = extractAfter(catgt_name, 'catgt_');

expo_folder = strrep(fileparts(fileparts(catgt_folder)), 'np_data', 'expo_data');

% tok = regexp(run_g, '^(.*)_g(\d+)$', 'tokens', 'once');
% runName = tok{1};
% runind  = str2double(tok{2});

% Expo stimulus files to use for this session
stimTag = { ...
    '[RFG_coarse2dg_99_4_150isi]', ...
    '[dir12_gpl_2_200isi_fixedphase]', ...
    '_2[Gpl2_2c_2sz_400_2_200isi]'};

%% ----------------------- Build shared session paths -----------------------

% Expo XML files for each stimulus block
expo_file = cell(1, numel(stimTag));
for i = 1:numel(stimTag)
    expo_file{i} = fullfile(expo_folder, sprintf('%s%s.xml', run_g, stimTag{i}));
    fprintf('expo_file       : %s\n', expo_file{i});
end


fprintf('catgt_folder  : %s\n', catgt_folder);

%% ----------------------- Read Expo files and segment info -----------------------

stiminfo{1}=extract_RF_stiminfo(expo_file{1});
stiminfo{2}=extract_GPlaid12_stiminfo(expo_file{2},1);
stiminfo{3}=extract_GPlaid2_sz2_contrast2_stiminfo(expo_file{3});
stiminfo{end+1}=stimTag;
 % Save the full stiminfo matrix
save(fullfile(catgt_folder, 'stiminfo.mat'), 'stiminfo');
%%
function stiminfo = extract_GPlaid2_sz2_contrast2_stiminfo(expo_file)
%% =========================================================================
%  Function: extract_GPlaid2_sz2_contrast2_stiminfo
%
%  Purpose:
%  Extract trial-by-trial stimulus information from a GPLaid Expo XML file
%  containing:
%      - 2 component directions
%      - 2 sizes
%      - 2 contrasts
%
%  Input:
%      expo_file : full path to Expo XML file
%
%  Output:
%      stiminfo  : structure with fields
%          .dir1          first component direction
%          .dir2          second component direction
%          .centerX       stimulus center x
%          .centerY       stimulus center y
%          .size          stimulus size
%          .contrast      stimulus contrast
%          .condition_ids trial index
%          .stim_name     'grating' or 'plaid'
%          .grating_dir   valid only for grating trials
%          .plaid_dir     valid only for plaid trials
% =========================================================================

expo_data = ReadExpoXML(expo_file);

% get all completed trials
[stimsynind, ~, ~, ppidx] = getExpo_pulses_fix(expo_data, ...
    'Fixation on', 'Reward', 'Ns fixation', 'fixation-start', 1);

events = expo_data.passes.events(ppidx);

stiminfo.dir1 = NaN(numel(stimsynind), 1);
stiminfo.dir2 = NaN(numel(stimsynind), 1);
gpl = NaN(numel(stimsynind), 1);   % 1 = grating, 2 = plaid
stiminfo.size=NaN(numel(stimsynind), 1);
stiminfo.contrast=NaN(numel(stimsynind), 1);

for i = 1:numel(stimsynind)
    % fixed phase, parameter position is different
    stiminfo.dir1(i) = events{stimsynind(i)}.Data{1,4}(2);
    stiminfo.dir2(i) = events{stimsynind(i)}.Data{1,6}(2);
    stiminfo.centerX(i)=events{stimsynind(i)}.Data{1,4}(5);
    stiminfo.centerY(i)=events{stimsynind(i)}.Data{1,4}(6);
    stiminfo.size(i) = events{stimsynind(i)}.Data{1,4}(3);
    stiminfo.contrast(i) = events{stimsynind(i)}.Data{1,3}(4);

    if events{stimsynind(i)}.Data{1,5}(7) == 0   % opacity on second plaid
        gpl(i) = 1;
    else
        gpl(i) = 2;
    end

end
% Canonicalize component directions before condition grouping.
% This prevents equivalent angles such as 360/0 or -10/350 from being
% treated as different conditions downstream.
stiminfo.dir1 = canonical_angle_deg(stiminfo.dir1, true);
stiminfo.dir2 = canonical_angle_deg(stiminfo.dir2, true);

assert_canonical_angle_deg(stiminfo.dir1, 'dir1');
assert_canonical_angle_deg(stiminfo.dir2, 'dir2');

% plaid direction
comdir1 = stiminfo.dir1(gpl == 2);
comdir2 = stiminfo.dir2(gpl == 2);
pldir = canonical_angle_deg( ...
    atan2d(sind(comdir1) + sind(comdir2), ...
           cosd(comdir1) + cosd(comdir2)), ...
    true);

% output fields
stiminfo.condition_ids = (1:numel(stiminfo.dir1))';

% stim_name
stiminfo.stim_name = strings(numel(stimsynind), 1);
stiminfo.stim_name(gpl == 1) = 'grating';
stiminfo.stim_name(gpl == 2) = 'plaid';

% grating_dir: only available for grating trials
stiminfo.grating_dir = NaN(numel(stimsynind), 1);
stiminfo.grating_dir(gpl == 1) = stiminfo.dir1(gpl == 1);

% plaid_dir: only available for plaid trials
stiminfo.plaid_dir = NaN(numel(stimsynind), 1);
stiminfo.plaid_dir(gpl == 2) = pldir;

assert_canonical_angle_deg(stiminfo.grating_dir, 'grating_dir');
assert_canonical_angle_deg(stiminfo.plaid_dir, 'plaid_dir');
end

%%
function stiminfo = extract_GPlaid12_stiminfo(expo_file, fixed_phase)
%% =========================================================================
%  Function: extract_GPlaid12_stiminfo
%
%  Purpose:
%  Extract trial-by-trial stimulus information from a GPLaid 12-direction
%  Expo XML file.
%
%  Input:
%      expo_file   : full path to Expo XML file
%      fixed_phase : flag indicating parameter layout in Expo event data
%                    0 -> non-fixed-phase
%                    1 -> fixed-phase
%
%  Output:
%      stiminfo    : structure with direction, center position, trial type,
%                    and parsed plaid/grating direction labels
% =========================================================================
expo_data = ReadExpoXML(expo_file);

% get all completed trials
[stimsynind, ~, ~, ppidx] = getExpo_pulses_fix(expo_data, ...
    'Fixation on', 'Reward', 'Ns fixation', 'fixation-start', 1);

events = expo_data.passes.events(ppidx);

stiminfo.dir1 = NaN(numel(stimsynind), 1);
stiminfo.dir2 = NaN(numel(stimsynind), 1);
gpl = NaN(numel(stimsynind), 1);   % 1 = grating, 2 = plaid
stiminfo.centerX=NaN(numel(stimsynind), 1);
stiminfo.centerY=NaN(numel(stimsynind), 1);

for i = 1:numel(stimsynind)
    if fixed_phase == 0
        stiminfo.dir1(i) = events{stimsynind(i)}.Data{1,5}(2);
        stiminfo.dir2(i) = events{stimsynind(i)}.Data{1,7}(2);
        stiminfo.centerX(i)=events{stimsynind(i)}.Data{1,5}(5);
        stiminfo.centerY(i)=events{stimsynind(i)}.Data{1,5}(6);
        if events{stimsynind(i)}.Data{1,6}(7) == 0   % opacity on second plaid
            gpl(i) = 1;
        else
            gpl(i) = 2;
        end

    elseif fixed_phase == 1
        % fixed phase, parameter position is different
        stiminfo.dir1(i) = events{stimsynind(i)}.Data{1,4}(2);
        stiminfo.dir2(i) = events{stimsynind(i)}.Data{1,6}(2);
        stiminfo.centerX(i)=events{stimsynind(i)}.Data{1,4}(5);
        stiminfo.centerY(i)=events{stimsynind(i)}.Data{1,4}(6);
        if events{stimsynind(i)}.Data{1,5}(7) == 0   % opacity on second plaid
            gpl(i) = 1;
        else
            gpl(i) = 2;
        end
    else
        error('fixed_phase must be 0 or 1.');
    end


end
% Canonicalize component directions before condition grouping.
% This prevents equivalent angles such as 360/0 or -10/350 from being
% treated as different conditions downstream.
stiminfo.dir1 = canonical_angle_deg(stiminfo.dir1, true);
stiminfo.dir2 = canonical_angle_deg(stiminfo.dir2, true);

assert_canonical_angle_deg(stiminfo.dir1, 'dir1');
assert_canonical_angle_deg(stiminfo.dir2, 'dir2');

% plaid direction
comdir1 = stiminfo.dir1(gpl == 2);
comdir2 = stiminfo.dir2(gpl == 2);
pldir = canonical_angle_deg( ...
    atan2d(sind(comdir1) + sind(comdir2), ...
           cosd(comdir1) + cosd(comdir2)), ...
    true);
% output fields
stiminfo.condition_ids = (1:numel(stiminfo.dir1))';

% stim_name
stiminfo.stim_name = strings(numel(stimsynind), 1);
stiminfo.stim_name(gpl == 1) = 'grating';
stiminfo.stim_name(gpl == 2) = 'plaid';

% grating_dir: only available for grating trials
stiminfo.grating_dir = NaN(numel(stimsynind), 1);
stiminfo.grating_dir(gpl == 1) = stiminfo.dir1(gpl == 1);

% plaid_dir: only available for plaid trials
stiminfo.plaid_dir = NaN(numel(stimsynind), 1);
stiminfo.plaid_dir(gpl == 2) = pldir;

assert_canonical_angle_deg(stiminfo.grating_dir, 'grating_dir');
assert_canonical_angle_deg(stiminfo.plaid_dir, 'plaid_dir');

end



function stiminfo=extract_RF_stiminfo(expo_file)
%% =========================================================================
%  Function: extract_RF_stiminfo
%
%  Purpose:
%  Extract RF mapping stimulus information from Expo XML.
%
%  Input:
%      expo_file : full path to RF mapping Expo XML file
%
%  Output:
%      stiminfo  : structure with fields
%          .xPos
%          .yPos
%          .ori
%          .stimsize
%          .filename
%          .condition_ids
% =========================================================================
expo_data = ReadExpoXML(expo_file);
%get all completed trials
[stimsynind,~,~,ppidx] = getExpo_pulses_fix(expo_data, 'Fixation on', 'Reward','Ns fixation','fixation-start',1);


[stiminfo.xPos, stiminfo.yPos, stiminfo.ori, stiminfo.stimsize] = get_RFGpara_hm(expo_data, stimsynind,ppidx);


stiminfo.condition_ids=1:numel(stiminfo.xPos);

end



function [xPos, yPos,ori, stimsize] = get_RFGpara_hm(expo_data, DO_blockIDs,pp_idx)
%% =========================================================================
%  Function: get_RFGpara_hm
%
%  Purpose:
%  Extract RF mapping parameters from Expo event data.
%  This function is specific to the RF mapping file structure and may not
%  generalize to other Expo paradigms.
%
%  Input:
%      expo_data    : structure returned by ReadExpoXML
%      DO_blockIDs  : indices of trials aligned to digital output/sync pulse
%      pp_idx       : pass indices
%
%  Output:
%      xPos         : stimulus x position
%      yPos         : stimulus y position
%      ori          : stimulus orientation
%      stimsize     : stimulus size
%
% =========================================================================

events = expo_data.passes.events(pp_idx);
xPos = NaN(numel(DO_blockIDs),1);
yPos = NaN(size(xPos));
ori=NaN(size(xPos));
stimsize=NaN(size(xPos));
for i = 1:numel(DO_blockIDs)

    xPos(i) = events{DO_blockIDs(i)}.Data{1,5}(5);
    yPos(i) = events{DO_blockIDs(i)}.Data{1,5}(6);
    ori(i)=events{DO_blockIDs(i)}.Data{1,5}(2);

    stimsize(i) = events{DO_blockIDs(i)}.Data{1,5}(3);

end

end



function a = canonical_angle_deg(a, doRound)
% Convert angles in degrees to canonical [0, 360).
%
% Examples:
%   360  -> 0
%   -10  -> 350
%   720  -> 0
%   NaN  -> NaN
%
% doRound = true is recommended for stimulus condition labels when
% directions are intended to be integer degrees.

if nargin < 2
    doRound = false;
end

if doRound
    a = round(a);
end

a = mod(a, 360);

% Remove tiny numerical residue around 0.
tol = 1e-10;
a(abs(a) < tol) = 0;
a(abs(a - 360) < tol) = 0;

end


function assert_canonical_angle_deg(a, fieldName)
% Check that finite angles are in [0, 360).

finiteVals = a(isfinite(a));

if any(finiteVals < 0 | finiteVals >= 360)
    error('Field %s contains non-canonical angles outside [0, 360).', fieldName);
end

end
