function TSWM_Simar1_EEG()

rand('state',sum(100*clock))  %mn moved up

%%%In this version of the program, the positions of the display elements
%%%jitter form trial to trial to help eliminate trunk-based coordination
%%%across trials
%cID = computerID()

eyeTrackerON = 0;
debug = 0; % make sure to go to lines 33-37 and fix else statement to = 1 for the EEG experiment
ptcpPrefix = 'ptcp';
ptcpNo = '1';
Timeout4EDF = 20;

%%%Permitted eye drift from fixation before trial aborts (in visual
%%%degrees)
allowableFixation = 2;
allowableTargetFixation = 3;
fs = 44100; % Hz
t.mid = 0:1/fs:0.25; % seconds
t.short = 0:1/fs:0.075; % seconds
t.long = 0:1/fs:0.5; % seconds
t.click = 0:1/fs:0.03; % seconds
f = 440; % Hz
f2 = 1760;
sounds.good = [sin(2.*pi.*f.*t.short),sin(2.*pi.*f2.*t.mid)];
sounds.click = [sin(2.*pi.*f2.*t.click)]./10;
sounds.noGood = [(sin(2.*pi.*f2.*t.short)).*0.5];


%% Initiate LabJack(LJ) & PsychToolBox(PTB)
if debug 
    recordingEEG = 0; % Do not want EEG to be recording if on debug
else
    recordingEEG = 1;
end
% if recordingEEG
%     [~, ljParams] = LJInit();
% end
% % LJ Pin Numbers
% pinNo.Glass = 0;
% pinNo.IRCurtain = 4;
% pinNo.IRBeamA = 0; pinNo.IRBeamB = 1;
% pinNo.LED = 6;
% % Setup PsychToolBox
% %HideCursor;
% %PsychDefaultSetup(2); %Screen('Preference', 'SkipSyncTests', 1);
% 
% %ListenChar(2); %Suppresses Matlab keyboard inputs
% % Keypresses
% RightControl = 163; KeyQ = 81; KeyB = 66; KeyY = 89;

%% Triggers
% Fixation point 1 = 111, Fixation point 2 = 211,
if recordingEEG
    Trig.sessionID = TriggerInit(2);
    
    Trig.FP1 = 111;
    Trig.FP2 = 211;
    
    %Trig.FP2 = 211:212    
    Trig.MemArr = 121;
    %Trig.MemArr = 121:128
    Trig.Cue = 220;
    Trig.moveOn = 196;
    Trig.response = 197;
    Trig.error = 249;
    
    Trig.pauseOFF = 254; Trig.pauseON = 255;
 
end
%% spatial trial features
sacSizeDeg = 4;
nSet = 8; %--max set size ~~mn 05/08/18  i.e., octagon

whichSets = [2];
if debug
    whichSets = [4];
end
nSets = length(whichSets);
%- number of contralateral/ipsilateral distractors for each given probe
%- designed for set sizes [2 4]
nCont1 = nSet./2; % only works with even-numbered set sizes!
nCont2 = nCont1 - 1;
nIpsi1 = nCont2;

nDir = 2; % number of saccade directions
nRepeats = 5;

nTrials = nSet.*nCont1.*nCont2.*nIpsi1.*nDir.*nSets.*nRepeats; % total number of trials of the experiment
nBlocks = 20; % entire EEG experiment
maxTrials = round(nTrials ./nBlocks); % max trials tested in a block, make sure it's an integer!!
if ~debug
    maxTrials_ = maxTrials;
else
    maxTrials_ = 10;
end

trialFeatures.nSet = nSet;
trialFeatures.nTrials = nTrials;
trialFeatures.nRepeats = nRepeats;
trialFeatures.whichSets = whichSets;
trialFeatures.nDir = nDir;
trialFeatures.nCont1 = nCont1;
trialFeatures.nCont2 = nCont2;
trialFeatures.nIpsi1 = nIpsi1;
trialFeatures.sacSizeDeg = sacSizeDeg;

% FP specs
fixL = 10;  5;       % size FP
targL = 20;     % size of stimuli and cursor
spotL = 2;      % size of colour cue inside cursor

%-- timing features
%- times at the start of trial
isiT = .5;
isiRangeT = 1;
fp1T = .1; %fixed b/c preceded by ISI
memArrT = .2; %or .4 ?
if debug
    memArrT = .2;
end
fp1postmemArrT = .5;
fp1postmemArrRangeT = 1;
fp1fp2OverlapT = .1;
fp2T = .5 - fp1fp2OverlapT;
fp2RangeT = 1;

maskT = .2; 
interFPdelayT = .5;
cueTime = 0.5;

%-- physical set up of display computer
eyeScreenDist_m = 1; 0.885;
ScreendWidth_m = 0.53136;  %Jiali 20th May

% eyeScreenDist_m = 0.885;
% ScreendWidth_m = 1.15;


%-- starting psychophysical toolbox functions
%screens=Screen('Screens');
%screenNumber=max(screens);
HC=1920/2; VC=1080/2;

ScreendWidth_deg = 2.* atan(ScreendWidth_m./2 ./eyeScreenDist_m) ./ pi .*180;
pixPERdeg = 2.*HC ./ ScreendWidth_deg;

trialFeatures.pixPERdeg = pixPERdeg;
trialFeatures.HC = HC;
trialFeatures.VC = VC;


%-- data saving
dataDir = 'C:\Simar\data\EEG_WM1';


ptcpDirInfo = GetPtcpDirInfo(dataDir, ptcpPrefix);
[pMax, pindexMax] = max(ptcpDirInfo.ptcpNo);
if pMax > 0
    disp(['the last participant was ' ptcpDirInfo.ptcpDirName{pindexMax}]);
end

prompt = {'Participant #:', 'Data Folder:'};
defaults = {num2str(pMax+1), dataDir};
answer = inputdlg(prompt, 'Data Selection', 1, defaults);
ptcpNo = answer{1};
dataDir = answer{2};

ptcpID = [ptcpPrefix, ptcpNo];
ptcpPath = [dataDir,filesep,ptcpID]; 

if exist(ptcpPath)
    
    str = ['Data folder ', ptcpPath, ' already exists. Okay to proceed? (Y/N) '];
    yesNo = input(str, 's'); %Displays above prompt, receives user input (y/n)
    if isempty(yesNo)
    elseif yesNo == 'y' || 'Y'
    elseif yesNo == 'n' || 'N'
        disp('okay, bye!')
        return
    end
    %%%%Exit program
else
    str = ['Data folder ', ptcpPath, ' does not exist. Create? (Y/N) '];
    yesNo = input(str, 's'); %Displays above prompt, receives user input (y/n)
    if isempty(yesNo)
        mkdir(ptcpPath);
    elseif yesNo == 'y' || yesNo == 'Y'
        % directory doesn't exist; create it.
        mkdir(ptcpPath);
    else
        disp('okay, bye!')
        return
        %%%%Exit program
    end
end

eval(['cd ', ptcpPath])
try
    load toDoList.mat % loads trials completedBlocksList lastTrial 
    if lastTrial == completedBlocksList(end).*maxTrials  %#ok<*NODEF>
        iBlock = length(completedBlocksList) + 1;
        disp(['Blocks saved so far ',num2str(completedBlocksList)])
        str = ['Procede with ', num2str(length(completedBlocksList) + 1),'? (Y/N) '];
        yesNo = input(str, 's');
        if isempty(yesNo)
            iTrial = lastTrial + 1;
            iBlock = length(completedBlocksList) + 1;
        elseif yesNo == 'y' || yesNo == 'Y'
            iTrial = lastTrial + 1;
            iBlock = length(completedBlocksList) + 1;
        elseif yesNo == 'n' || 'N'
            str = ['Please enter the # of the block you want test '];
            iBlock = input(str);
            guessTrial = (iBlock-1)*maxTrials + 1;
            str = ['Okay to start from the first trial (trial #', num2str(guessTrial), ')? (Y/N)'];
            yesNo = input(str, 's');
            if isempty(yesNo)
                iTrial = guessTrial;
            elseif yesNo == 'y' || yesNo == 'Y'
                iTrial = guessTrial;
            elseif yesNo == 'n' || 'N'
                str = ['Please enter trial # '];
                iTrial = input(str);
                if (iTrial < guessTrial) || (iTrial > guessTrial+maxTrials-1)
                    disp('Sorry, the entries are inconsistent. Shutting down...')
                    return
                end
            end
        end
    elseif lastTrial < completedBlocksList(end).*maxTrials
        
        disp(['The last block (#', num2str(completedBlocksList(end)), ') seems to be incomplete.'])
        disp(['It ends with trial #', num2str(lastTrial), ' rather than ', num2str(completedBlocksList(end).*maxTrials), '.'])
        
        str = ['Okay to continue with trial #', num2str(lastTrial+1), '? (Y/N) '];
        yesNo = input(str, 's');
        if isempty(yesNo)
            iTrial = lastTrial+1;
            iBlock = completedBlocksList(end);
        elseif yesNo == 'n' || 'N'
            str = ['Please enter trial # '];
            iTrial = input(str);
            iBlock = ceil(iTrial/ maxTrials);
            
            str = ['Sure, that means I continue with block# ', num2str(iBlock), '. Okay to proceed? (Y/N) '];
            yesNo = input(str, 's');
            if isempty(yesNo)
            elseif yesNo == 'n' || 'N'
                disp('Okay, better to start from the beginning. ')
                return
            end
            
        end
        
    elseif lastTrial > completedBlocksList(end).*maxTrials
        disp('This should never happen. Please call Matthias.')
        return        
    end
    
catch
    disp('Creating trial list') 
    trials = permutateTrials(trialFeatures);
    iBlock = 1;
    iTrial = 1;
    completedBlocksList = [];
end


disp('Press a key to start!')
KbWait()




%--set up 'Screen' function
Screen('Preference', 'SkipSyncTests', 0);
screenNumber=0;
[window, screenRect]=Screen(screenNumber,'OpenWindow', 0, [0 0 2.*HC 2.*VC]);
%[window, screenRect]=Screen(screenNumber,'OpenWindow', 0, [1 1 2.*HC 2.*VC]);


Screen(window,'FillRect',BlackIndex(window));
Screen('Flip', window);

screenResH = screenRect(3)-screenRect(1);
screenResV = screenRect(4)-screenRect(2);

HC = screenResH/2; VC = screenResV/2;

nScreen = 10;     % trial needs frames: blank - fp1 - fp2 - cue - fp+T1 - T2 - response

%-- Open a graphics window on the main Screen
white=WhiteIndex(window);
black=BlackIndex(window);
grey=(white+black)/2;
inc=white-grey;
colBack = grey;
colFP = white;

%MemoryArrayColours = [white grey 0; white 0 grey; grey white 0; 0 white grey; grey 0 white; 0 grey white];

%MemoryArrayColours = [white white.*2 white.*2; 0 white 0; white.*.3 white.*.5 white; white white 0];
MemoryArrayColours = white.*[1 .3 .3; .1 .9 .2; .3 .7 1; 1 .9 0];
MemoryArrayColours = white.*[1 .5 .5; .1 .9 .2; .3 .7 1; 1 .9 0];


% currentLogID = fopen([dataDir,filesep,'currentLog.txt'], 'w'); %--mn
% sessionLogID = fopen([ptcpPath,filesep,'log.txt'], 'w');
% fprintf(currentLogID,[num2str(clock),'\r\n'])
% fprintf(sessionLogID,[num2str(clock),'\r\n'])
% 
% fprintf(currentLogID,[ptcpID, '\r\n\r\n'])
% fprintf(sessionLogID,[ptcpID, '\r\n\r\n'])

warning off MATLAB:DeprecatedLogicalAPI5    %problem with KbCheck
warning off MATLAB:DeprecatedLogicalAPI

fname = 'SSD.edf';
fnameData = ['TSWM', num2str(iBlock), '.mat'];
%====== some parameters
eye_used = -1;
el = [];

KbName('UnifyKeyNames') 

if eyeTrackerON
    blockNo = iBlock; % + iBlock;
    
    
    path = [ptcpPath, filesep, 'block',num2str(blockNo)];
    
    if exist(path, 'file') == 7
        error('the folder I''m trying to save to already exists! Shutting down to avoid overwriting data.')
    else
        mkdir(path);
    end
    saveAs = [path,filesep,fnameData];
else
    saveAs = [ptcpPath,filesep,fnameData];
end

stopkey=KbName('space');
goodTrialCount = 0;
badTrialCount = 0;
iFile = 1;
trialRejected = false;
HideCursor;
keyResult.recalibrate = false;

theData = [];

if recordingEEG
    [~, ~] = SendTrigger(Trig.sessionID, Trig.pauseOFF);
end


while iTrial<=iBlock.*maxTrials_ %run until trials used up or max # of trials
    
    theTimes.trialStart = GetSecs;
    
    if trialRejected
        %--flush variables and Screens - or memory explodes
        windowPtrs=Screen('Windows');
        for iPtr = 1:length(windowPtrs)
            isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
            if isOffscreen
                Screen(windowPtrs(iPtr),'Close');
            end
        end
        
        keyResult.recalibrate = false;
        
        fadeDuration = 1;
        gaussScaling = fadeDuration/4;
        fadeStart = GetSecs;
        badTrialCount = badTrialCount+1;
        while GetSecs < fadeStart+fadeDuration
            currentLum = (normpdf((GetSecs-fadeStart)/gaussScaling))/normpdf(0)*white;
            Screen(window,'FillRect',[0,0,currentLum]);
            Screen('Flip', window, 0, 0, 2);
        end
        trialRejected = false;
        %                 fprintf(currentLogID,'Fixation error:\r\n'); % Feb15 extremely useful but Matthias thinks this might be unnecessary (in case there still is a problem)
        %                 fprintf(currentLogID,[errorIn, '\r\n']);
        %                 fprintf(currentLogID,['So far there have been ' num2str(goodTrialCount) ' good trials and ' num2str(badTrialCount) ' bad\r\n']);
        %                 fprintf(sessionLogID, 'Fixation error:\r\n');
        %                 fprintf(sessionLogID,[errorIn, '\r\n']);
        %                 fprintf(sessionLogID,['So far there have been ' num2str(goodTrialCount) ' good trials and ' num2str(badTrialCount) ' bad\r\n']);
    end
    %             fclose(currentLogID);  %--mn
    %             fclose(sessionLogID);
    %             currentLogID = fopen([dataDir,filesep,'currentLog.txt'], 'a');
    %             sessionLogID = fopen([path,filesep,'log.txt'], 'a');
    errorIn = [];
    
    if eyeTrackerON
        
        Eyelink('openfile', 'SSD.edf');
        %--- some of these commands need to be repeated, I believe
        Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
        Eyelink('Command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
        Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON');
        Eyelink('Command', 'heuristic_filter = 0')
        if ScreendWidth_m > 1.5
            redctFact = 1/4; calR = round([0+redctFact.*screenResH 0+redctFact.*screenResV screenResH-redctFact.*screenResH screenResV-redctFact.*screenResV]);
            eval(['Eyelink(''Command'', ''Screen_pixel_coords = ', num2str(calR), ''')']);
        end
        % STEP 5
        % start recording eye position
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        % mark zero-plot time in data file
        Eyelink('Message', 'SYNCTIME');
    end
    
    %-- select trial features
    stimJitter = round((rand-.5 )*400);
    
    selFPX = trials.totalFPX(iTrial)+stimJitter;
    selFPY = trials.totalFPY(iTrial);
    
    selProbeX = trials.totalProbeX(iTrial)+stimJitter;
    selProbeY = trials.totalProbeY(iTrial);
    selProbeIndex = trials.totalProbeIndex(iTrial);
    selCont1X = trials.totalCont1X(iTrial)+stimJitter;
    selCont1Y = trials.totalCont1Y(iTrial);
    selCont2X = trials.totalCont2X(iTrial)+stimJitter;
    selCont2Y = trials.totalCont2Y(iTrial);
    selIpsi1X = trials.totalIpsi1X(iTrial)+stimJitter;
    selIpsi1Y = trials.totalIpsi1Y(iTrial);
    
    selMemoryArrayX = [selProbeX selCont1X selCont2X selIpsi1X];
    selMemoryArrayY = [selProbeY selCont1Y selCont2Y selIpsi1Y];
    
    selsacTargX = trials.totalsacTargX(iTrial)+stimJitter;
    selsacTargY = trials.totalsacTargY(iTrial);
    
%     if trials.totalsacTargX(iTrial) < HC % for EEG trigger signal
%         iTrigSacDir = 1;
%     else
%         iTrigSacDir = 2;
%     end   %Jiali 20th May


    mouseHasMoved = 0;
    moveOnT = NaN;
    
    % create position squares for stimuli
    fpx_1 = selFPX;
    fpy_1 = selFPY;
    fpx_2 = selsacTargX;
    fpy_2 = selsacTargY;
    theRectFP = round([fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(fixL./2));
    theRectFP2 = round([fpx_2 fpy_2 fpx_2 fpy_2] + [-1 -1 1 1].*(fixL./2));
    
    %============= prepare off-Screen windows
    startTrialT = GetSecs;
    
    if recordingEEG
        [~, ~] = SendTrigger(Trig.sessionID, Trig.FP1);
    end
    
    for iScreen = 1:nScreen
        offwin(iScreen)=Screen(window,'OpenOffScreenWindow',colBack, [0 0 screenResH screenResV]);
    end
    iScreen = 0;
    
    %--frame 1: fp1 dimm
    iScreen = iScreen + 1; iFPDim = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    Screen(offwin(iScreen),'FillOval',colFP.*0.75,theRectFP);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(spotL./2)));
    
    %--frame 1a: fp1 bright
    iScreen = iScreen + 1; iFP = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    Screen(offwin(iScreen),'FillOval',colFP,theRectFP);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(spotL./2)));
    
    %--frame 2: fp1 + memory array frame
    iScreen = iScreen + 1; iMemoryArray = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    Screen(offwin(iScreen),'FillOval',colFP,theRectFP);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(spotL./2)));
    
    [dummy, iSort] = sort(rand(size(MemoryArrayColours,1),1));
    for iMemArr = 1 : length(selMemoryArrayX)
        theRectStim = round([selMemoryArrayX(iMemArr) selMemoryArrayY(iMemArr) selMemoryArrayX(iMemArr) selMemoryArrayY(iMemArr)] + [-1 -1 1 1].*(targL./2));
        Screen(offwin(iScreen),'FillRect', MemoryArrayColours(iSort(iMemArr),:), theRectStim)
    end
    colProbe = MemoryArrayColours(iSort(1),:); % will be the colour of the cue
    
    %--frame 3: fp1 + fp2 frame to slow down saccade latencies
    iScreen = iScreen + 1; iFP1FP2 = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    Screen(offwin(iScreen),'FillOval',colFP,theRectFP);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(spotL./2)));
    Screen(offwin(iScreen),'FillOval',colFP,theRectFP2);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_2 fpy_2 fpx_2 fpy_2] + [-1 -1 1 1].*(spotL./2)));
    
    %--frame 4: fp2 frame
    iScreen = iScreen + 1; iFP2 = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    Screen(offwin(iScreen),'FillOval',colFP,theRectFP2);
    Screen(offwin(iScreen),'FillOval',black,round([fpx_2 fpy_2 fpx_2 fpy_2] + [-1 -1 1 1].*(spotL./2)));
    
    %--frame 5: cue (blank frame)
    iScreen = iScreen + 1; iCue = iScreen;
    Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
    
    
    %================ bring off-Screens to front to run trial
    
    if eyeTrackerON
        % Check recording status, stop display if error
        errorEL=Eyelink('CheckRecording');
        if (errorEL~=0), Screen('closeall'); Eyelink('Shutdown'); 'eyelink connection not working'; return; end;
        
        %- first sample eye position
        
        % STEP 6 - monitor eye
        if eye_used == -1 % if we don't know which eye to use, first find eye that's being tracked
            eye_used = Eyelink( 'EyeAvailable'); % get eye that's tracked
            if eye_used == el.BINOCULAR; % if both eyes are tracked
                eye_used = el.LEFT_EYE; % use left eye
            end
        end
        
    end
    
    if eyeTrackerON
        %-- frame 1: FP1, wait for subj fixation and mouse click ------------------------------------------------
        [x,y,keyIsDown1] = GetMouse;
        keyIsDown1 = [0 0 0];
        Screen('CopyWindow',offwin(iFPDim),window);
        Screen('Flip', window, 0, 0, 2);
        
        Eyelink('Message', 'FP1');
        
        wasZero = false;
        while 1 && ~trialRejected
            [notUsed,notUsed,keyIsDown1] = GetMouse;
            if keyIsDown1(1) == 0
                wasZero = true;
            end
            if wasZero == true && keyIsDown1(1) == 1
                break
            end
        end
        if trialRejected == true;
            continue
        end
        bufferSize = 10;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = 0.1;
        allowableOffset = 20;
        fixationOffset.x = 0;
        fixationOffset.y = 0;
        realFixation.x = 0;
        realFixation.y = 0;
        
        [notUsed, realFixation, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation,allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        %             fprintf(currentLogID,notes);
        %             fprintf(sessionLogID,notes);
        
        fixationOffset.x = fixLocation.x - realFixation.x;
        fixationOffset.y = fixLocation.y - realFixation.y;
        %             fprintf(currentLogID,['Fixation offset = ' num2str(sqrt((fixationOffset.x^2)+(fixationOffset.y^2))/pixPERdeg) '\r\n']);
        %             fprintf(sessionLogID,['Fixation offset = ' num2str(sqrt((fixationOffset.x^2)+(fixationOffset.y^2))/pixPERdeg) '\r\n']);
        
    end
    
    %-- frame 1a: bright FP1 ------------------------------------------------
    Screen('CopyWindow',offwin(iFP),window);
    Screen('Flip', window, 0, 0, 2);
    
    holdTime = isiT + isiRangeT.* rand + fp1T;
    
    if eyeTrackerON
        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = (fpT+ 2.*(rand-.5).*fpRangeT);
        allowableOffset = allowableTargetFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        %             fprintf(currentLogID,notes);
        %             fprintf(sessionLogID,notes);
        
        if trialRejected
            errorIn = 'Initial fixation';
            continue
        end
    else
        WaitSecs(holdTime)
    end
    
    %-- frame 2: FP1 + memory array ------------------------------------------------
    Screen('CopyWindow',offwin(iMemoryArray),window);
    Screen('Flip', window, 0, 0, 2);
    
    if recordingEEG
          [~, ~] = SendTrigger(Trig.sessionID, Trig.MemArr);
        
%         [~, ~] = SendTrigger(Item1 location)
%         [~, ~] = SendTrigger(Item2 location)
%         [~, ~] = SendTrigger(Item3 location)
%         [~, ~] = SendTrigger(Item4 location)
%         
         StartTime = GetSecs;
%          [~, ~] = SendTrigger(Trig.sessionID, 11)
%         [~, ~] = SendTrigger(Trig.sessionID, 12)
%         [~, ~] = SendTrigger(Trig.sessionID, 13)
%         [~, ~] = SendTrigger(Trig.sessionID, 14)
%         TotalTime = GetSecs - StartTime;
        
        % [~, ~] = SendTrigger(rig.sessionID, 102^80000000)
        
    end
    
    theData.T_of_MemArr_FP2_Cue(iFile,1) = GetSecs - startTrialT;
    
    holdTime = memArrT;
    
    if eyeTrackerON
        Eyelink('Message', 'bar');
        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = barT;
        allowableOffset = allowableFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        %             fprintf(currentLogID,notes);
        %             fprintf(sessionLogID,notes);
        if trialRejected
            errorIn = 'Memory encoding';
            continue
        end
    else
        WaitSecs(holdTime)
    end
    
    %-- again frame 1: FP1 ------------------------------------------------
    Screen('CopyWindow',offwin(iFP),window);
    Screen('Flip', window, 0, 0, 2);
    
    holdTime = fp1postmemArrT + fp1postmemArrRangeT.*rand;
    
    if eyeTrackerON
        Eyelink('Message', 'FP1');
        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = interFPdelayT;
        allowableOffset = allowableFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        %             fprintf(currentLogID,notes);
        %             fprintf(sessionLogID,notes);
        if trialRejected
            errorIn = 'Just before saccade';
            continue
        end
    else
        WaitSecs(holdTime)
    end
    
    %-- frame 3: FP1 + FP2 ------------------------------------------------
    Screen('CopyWindow',offwin(iFP1FP2),window);
    Screen('Flip', window, 0, 0, 2);
    
    if recordingEEG
        [~, ~] = SendTrigger(Trig.sessionID, Trig.FP2);
    end
    
    theData.T_of_MemArr_FP2_Cue(iFile,2) = GetSecs - startTrialT;
    
    holdTime = fp1fp2OverlapT;
    WaitSecs(holdTime)
    
    
    %-- frame 4: FP2 alone ------------------------------------------------
    Screen('CopyWindow',offwin(iFP2),window);
    Screen('Flip', window, 0, 0, 2);
    
    holdTime = fp2T + fp2RangeT.*rand;
    
    if eyeTrackerON
        Eyelink('Message', 'FP2');
        bufferSize = 3;
        fixLocation.x = fpx_2;
        fixLocation.y = fpy_2;
        holdTime = maskT;
        allowableOffset = allowableTargetFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        %             fprintf(currentLogID,notes);
        %             fprintf(sessionLogID,notes);
        if trialRejected
            errorIn = 'Second mask andor refixation';
            continue
        end
    else
        WaitSecs(holdTime)
    end
    
    %-- frame 5: cue = FP2 ------------------------------------------------
    Screen('CopyWindow',offwin(iCue),window);
    mouseLoc.x = fpx_2;
    mouseLoc.y = fpy_2;
    SetMouse(round(mouseLoc.x), round(mouseLoc.y));
    Screen(window,'FillOval',colProbe,round([mouseLoc.x mouseLoc.y mouseLoc.x mouseLoc.y] + [-1 -1 1 1].*(fixL./2)));
    Screen(window,'FillOval',black,round([mouseLoc.x mouseLoc.y mouseLoc.x mouseLoc.y] + [-1 -1 1 1].*(spotL./2)));
    Screen('Flip', window, 0, 0, 2);
    
    if recordingEEG        
        [~, ~] = SendTrigger(Trig.sessionID, Trig.Cue);        
    end
    
    cueT = GetSecs;
    theData.T_of_MemArr_FP2_Cue(iFile,3) = cueT - startTrialT;
    
    if eyeTrackerON
        Eyelink('Message', 'Cue');
    end
    
    [x,y,keyIsDown] = GetMouse;
    Screen(window,'FillOval',colProbe,round([x y x y] + [-1 -1 1 1].*(fixL./2)));
    Screen(window,'FillOval',black,round([x y x y] + [-1 -1 1 1].*(spotL./2)));
    %WaitSecs(1)
    waitForCursorT = GetSecs;
    while 1
        [x,y,keyIsDown] = GetMouse; % update cursor position
        
        % if cursor starts to move
        pixlCrit = 2;
        if ~mouseHasMoved
            if (mouseLoc.x - x).*(mouseLoc.x - x)+(mouseLoc.y - y).*(mouseLoc.y - y) > pixlCrit.*pixlCrit
                
                if recordingEEG
                    [~, ~] = SendTrigger(Trig.sessionID, Trig.moveOn);
                end
                
                moveOnT = GetSecs - cueT;
                mouseHasMoved=1;
            end
        end
        
        if GetSecs - waitForCursorT > 1
            break
        end
    end
    SetMouse(round(mouseLoc.x), round(mouseLoc.y));
    mouseHasMoved=0;
    
    while 1 %~keyIsDown(1)
        [x,y,keyIsDown] = GetMouse; % update cursor position
        Screen('CopyWindow',offwin(iCue),window);
        Screen(window,'FillOval',colProbe,round([x y x y] + [-1 -1 1 1].*(fixL./2)));
        Screen('Flip', window, 0, 0, 2);
        
        % if cursor starts to move
        pixlCrit = 2;
        if ~mouseHasMoved
            if (mouseLoc.x - x).*(mouseLoc.x - x)+(mouseLoc.y - y).*(mouseLoc.y - y) > pixlCrit.*pixlCrit
                
                if recordingEEG
                    [~, ~] = SendTrigger(Trig.sessionID, Trig.moveOn);
                end
                
                moveOnT = GetSecs - cueT;
                mouseHasMoved=1;
            end
        end
        
        % if mouse key is pressed
        if keyIsDown(1)
            
            if recordingEEG
                [~, ~] = SendTrigger(Trig.sessionID, Trig.response);
            end
            
            responseT = GetSecs - cueT;
            responseX = x;
            responseY = y;
            break
        end
        
    end
    
    %-- error signal, ptcpt clicked mouse w/o moving it
    if isnan(moveOnT)
        if recordingEEG
            [~, ~] = SendTrigger(Trig.sessionID, Trig.error);
        end
        beep
    end
    
    if eyeTrackerON
        Eyelink('Message', 'response');
    end
    
    %--store trial data
    theData.iTrial(iFile) = iTrial;
    theData.selFPX(iFile) = selFPX;
    theData.selFPY(iFile) = selFPY;
    theData.selProbeX(iFile) = selProbeX;
    theData.selProbeY(iFile) = selProbeY;
    theData.selCont1X(iFile) = selCont1X;
    theData.selCont1Y(iFile) = selCont1Y;
    theData.selCont2X(iFile) = selCont2X;
    theData.selCont2Y(iFile) = selCont2Y;
    theData.selIpsi1X(iFile) = selIpsi1X;
    theData.selIpsi1Y(iFile) = selIpsi1Y;
    theData.selsacTargX(iFile) = selsacTargX;
    theData.selsacTargY(iFile) = selsacTargY;
    theData.responseX(iFile) = responseX;
    theData.responseY(iFile) = responseY;
    theData.stimJitterX(iFile) = stimJitter;
    theData.responseT(iFile) = responseT;
    theData.moveOnT(iFile) = moveOnT;
    %theData.TriggerCue(iFile) = Trig.Cue;
    
    % theData.T_of_MemArr_FP2_Cue already completed
    
    save(saveAs, '-mat', 'theData');
    
    %--flush variables and Screens - or memory explodes
    windowPtrs=Screen('Windows');
    for iPtr = 1:length(windowPtrs)
        isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
        if isOffscreen
            Screen(windowPtrs(iPtr),'Close')
        end
    end
    
    %-- Save info about state of experiment
    lastTrial = iTrial;
    completedBlocksList(iBlock) = iBlock;
    save toDoList.mat -mat trials completedBlocksList lastTrial
    %             fprintf(currentLogID,['Trial ', num2str(iTrial), ' took ' , num2str(GetSecs-theTimes.trialStart) 's and looks good!\r\n']);
    %             fprintf(currentLogID,[num2str(sum(~isnan(totalFPX))),' trials remain\r\n\r\n']);
    %             fprintf(sessionLogID,['Trial ', num2str(iTrial), ' took ' , num2str(GetSecs-theTimes.trialStart) 's and looks good!\r\n']);
    %             fprintf(sessionLogID,[num2str(sum(~isnan(totalFPX))),' trials remain\r\n\r\n']);
    
    %-- update trial counts
    iTrial = iTrial + 1;
    iFile = iFile + 1;
    goodTrialCount = goodTrialCount+1;
    
end


% Pause On; Stop recording EEG Data
if recordingEEG
    [~, ~] = SendTrigger(Trig.sessionID, Trig.pauseON);
end

if eyeTrackerON==1
    Eyelink('closefile');
    Eyelink('Shutdown');
end
Screen('CloseAll');
ShowCursor
end



function cID = computerID()
cID = '';
ni = java.net.NetworkInterface.getNetworkInterfaces;
while ni.hasMoreElements
    addr = ni.nextElement.getHardwareAddress;
    if ~isempty(addr)
        addrStr = dec2hex(int16(addr)+128);
        cID = [cID, '.', reshape(addrStr,1,2*length(addr))];
    end
end

end


function ptcpDirInfo = GetPtcpDirInfo(dataDir, ptcpPrefix)
ptcpDirs = dir([dataDir,filesep,ptcpPrefix,'*']);
if length(ptcpDirs) >= 1
    for iPtcp = 1:length(ptcpDirs)
        thisPtcp = ptcpDirs(iPtcp).name;
        thisPtcpNo = str2num(thisPtcp(length(ptcpPrefix)+1:end));
        ptcpDirInfo.ptcpDirName{iPtcp} = thisPtcp;
        ptcpDirInfo.ptcpNo(iPtcp) = thisPtcpNo;
    end
else
    ptcpDirInfo.ptcpDirName = [];
    ptcpDirInfo.ptcpNo = [];
end
end

function [trialRejected, eye, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation,allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg)
startTime = GetSecs;
buffered = false;
bufferPos = 1;
eye.x = 0;
eye.y = 0;
trialRejected = false;
notes = '';
while GetSecs < startTime+holdTime
    if eyeTrackerON==1 %
        error=Eyelink('CheckRecording');
        if(error~=0)
            break;
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            % if we do, get current gaze position from sample
            liveEye.x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
            liveEye.y = evt.gy(eye_used+1);
            % do we have valid data and is the pupil visible?
            if liveEye.x~=el.MISSING_DATA && liveEye.y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                if buffered || bufferPos > bufferSize
                    eye.x = mean(eyeBuffer.x);
                    eye.y = mean(eyeBuffer.y);
                    %%%Check to see if the eye position is within
                    %%%range:
                    fixDiff.x = fixLocation.x - eye.x - fixationOffset.x;
                    fixDiff.y = fixLocation.y - eye.y - fixationOffset.y;
                    fixDiff.rho = sqrt((fixDiff.x^2)+(fixDiff.y^2));
                    notes = ['fixation offset was ' num2str(fixDiff.rho/pixPERdeg) ' degrees\r\n'];
                    if (fixDiff.rho/pixPERdeg) > allowableOffset
                        trialRejected = true;
                        break
                    end
                    buffered = true;
                    bufferPos = 1;
                end
                eyeBuffer.x(bufferPos)=liveEye.x;
                eyeBuffer.y(bufferPos)=liveEye.y;
                bufferPos = bufferPos+1;
            end
        end
    else
        % Query current mouse cursor position (our "pseudo-eyetracker") -
        [eye.x, eye.y, notUsed]=GetMouse;
    end
end
end

function [pauseCheck] = myKBInput_PausePlease(pauseKey)
pauseCheck = false;
[kbKeyDown, responseClock, keyCode] = KbCheck;
if kbKeyDown
    KBkeys = KbName(keyCode);
    if iscell(KBkeys)
        KBkey = KBkeys{1};
    else
        KBkey = KBkeys;
    end

    if strcmp(KBkey,pauseKey)
        pauseCheck = true;
    else
    end
end
end

function [keyResult] = pauseScreen(window, screenRect, grey, white)

myKeys.unPauseKey = 'u';
myKeys.recalibrate = 'r';
myKeys.quit = 'q';
Screen(window,'FillRect',grey);

Screen('TextSize',window, 36);

myPrompt = ['Paused!'];
Screen(window,'DrawText',myPrompt,50,50,white);

myPrompt = [myKeys.unPauseKey ' = unpause'];
Screen(window,'DrawText',myPrompt,50,150,white);

myPrompt = [myKeys.recalibrate ' = recalibrate'];
Screen(window,'DrawText',myPrompt,50,250,white);

myPrompt = [myKeys.quit ' = quit'];
Screen(window,'DrawText',myPrompt,50,350,white);

Screen('Flip', window, 0, 0, 2);

validResp = false;
while ~validResp
    [keyResult, validResp] = myKBInput_PausePrompt(myKeys);
end


if keyResult.quit
    quitProg
end

end

function [keyResult, validResp] = myKBInput_PausePrompt(myKeys)

validResp = false;
keyResult.unpause = false;
keyResult.recalibrate = false;
keyResult.quit = false;

pauseKey = myKeys.unPauseKey;
recalKey = myKeys.recalibrate;
quitKey = myKeys.quit;

[kbKeyDown, responseClock, keyCode] = KbCheck;
if kbKeyDown
    KBkeys = KbName(keyCode);
    if iscell(KBkeys)
        KBkey = KBkeys{1};
    else
        KBkey = KBkeys;
    end

    if strcmp(KBkey,pauseKey)
        keyResult.unpause = true;
        validResp = true;
    elseif strcmp(KBkey,recalKey)
        keyResult.recalibrate = true;
        validResp = true;
    elseif strcmp(KBkey,quitKey)
        keyResult.quit = true;
        validResp = true;
    end
end

end

% function [el] = calibrateEL(el,window,recal) %% Jiali 9th April 2018
% if recal
%     Eyelink('CloseFile')
%     Eyelink('Stoprecording')
% end
% Eyelink('Shutdown');
% if Eyelink('Initialize');
%     return;
% end;
% EyelinkInit(0);
% % STEP 3
% % Provide eyelink with details about the graphics environment
% % and perform some initializations. The information is returned
% % in a structure that also contains useful defaults
% % and control codes (e.g. tracker state bit and eyelink key  values).
% el=EyelinkInitDefaults(window);
% % make sure that we get gaze data from the eyelink
% Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
% Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
% Eyelink('Command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
% Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON');
% Eyelink('Command', 'heuristic_filter = 0')
% Eyelink('Command', 'set_calibration_colors([100 0 0], [0 0 0])');
% 
% el.calibrationtargetsize = 1;
% el.calibrationtargetwidth = .5;
% el.displayCalResults = 1;
% el.eyeimgsize = 50;
% 
% % STEP 4
% % Calibrate the eye tracker using the standard calibration routines
% %eyelink('trackersetup');
% %     % you must send this command with value NO for custom calibration
% %     % you must also reset it to YES for subsequent experiments
% 
% 
% Eyelink('command', 'calibration_type = HV9');
% % you must send this command with value NO for custom calibration
% % you must also reset it to YES for subsequent experiments
% Eyelink('command', 'generate_default_targets = NO');
% 
% % STEP 5.1 modify calibration and validation target locations
% Eyelink('command','calibration_targets = 1360,540 960,140 560,940 560,540 560,140 960,540 960,940 1360,140 1360,940');
% Eyelink('command','validation_targets = 1360,540 960,140 560,940 560,540 560,140 960,540 960,940 1360,140 1360,940');
% 
% 
% %Eyelink('StartSetup');
% errorID=EyelinkDoTrackerSetup(el)
% % do a final check of calibration using driftcorrection
% success=EyelinkDoDriftCorrection(el,[],[],1,1)
% if isfield(el,'adamDebug')
%     sca
%     x = 1
% end
% 
% if success ~= 1
%     return;
% end
% end

function [] = quitProg()
Priority(0);
ShowCursor
Screen('CloseAll');
error('Experiment terminated!')
end

% function [] = blockScreen(window, screenRect, grey, white, conditionChange)
% Screen(window,'FillRect',grey);
% 
% if ~conditionChange
% myPrompt = ['Block complete!'];
% Screen(window,'DrawText',myPrompt,50,50,white);
% 
% myPrompt = ['Please take a short break if you like'];
% Screen(window,'DrawText',myPrompt,50,150,white);
% 
% myPrompt = ['If there''s anything you need, please call for the experimenter'];
% Screen(window,'DrawText',myPrompt,50,250,white);
% 
% myPrompt = ['When you''re ready, click mouse to continue'];
% Screen(window,'DrawText',myPrompt,50,350,white);
% 
% else
% 	
% myPrompt = ['We''re going to make a change now'];
% Screen(window,'DrawText',myPrompt,50,50,white);
% 
% myPrompt = ['Please call for the experimenter'];
% Screen(window,'DrawText',myPrompt,50,150,white);	
% end
% 
% Screen('Flip', window, 0, 0, 2);
% GetClicks;
% 
% end



function mR = rot8Matrix2D_inProgram(m, a)

R=[(cos(a)) (sin(a)); (-sin(a)) (cos(a))];
mR = m*R;

end


function trials = permutateTrials(trialFeatures)

nSet = trialFeatures.nSet;
nRepeats = trialFeatures.nRepeats;
whichSets = trialFeatures.whichSets;
nDir = trialFeatures.nDir;
nCont1 = trialFeatures.nCont1;
nCont2 = trialFeatures.nCont2;
nIpsi1 = trialFeatures.nIpsi1;
sacSizeDeg = trialFeatures.sacSizeDeg;
nTrials = trialFeatures.nTrials;

pixPERdeg = trialFeatures.pixPERdeg;
HC = trialFeatures.HC;
VC = trialFeatures.VC;

sacc_degX = [-1 1] .* sacSizeDeg.*1.5; %-- sac target outside octagon, 50% further
sacc_degY = [0 0];


%-- permutation of trials
rot8Angles = linspace(0, 2.*pi, nSet+1) + pi./nSet; rot8Angles(end) = [];
t0 = sacSizeDeg.*[0 0; 1 0];
xing_degX=NaN.*zeros(nSet,1); xing_degY=xing_degX;
for iRot = 1 : nSet
    a = rot8Angles(iRot);
    R=[(cos(a)) (sin(a)); (-sin(a)) (cos(a))];
    mR = t0*R;
    xing_degX(iRot,1) = mR(2,1);
    xing_degY(iRot,1) = mR(2,2);
end
ProbeIndex = [1:nSet]';

%--- build variable lists
totalFPX = [];
totalFPY = [];
totalProbeX = [];
totalProbeY = [];
totalProbeIndex = [];

%-- complete permutation of distractor locations for set sizes 2 & 4
totalIpsi1X = [];
totalIpsi1Y = [];
totalCont1X = [];
totalCont1Y = [];
totalCont2X = [];
totalCont2Y = [];

totalsacTargX = [];
totalsacTargY = [];


for iRepeat = 1 : nRepeats
    for iSetSize = whichSets
        
        %nRepeats = -11.5 * iSetSize + 47; % for 24 repeats for sets=2, 1 repeat for sets=4
        vF = (iSetSize-2)/(iSetSize-2); % "vanishing factor", NaN for set size = 2, 1 for set size = 4
        
        for iDir = 1 : nDir
            
            totalsacTargX = [totalsacTargX; sacc_degX(iDir).*ones(nSet.*nCont1.*nCont2.*nIpsi1,1)];
            totalsacTargY = [totalsacTargY; sacc_degY(iDir).*ones(nSet.*nCont1.*nCont2.*nIpsi1,1)];
            
            for iProbe = 1 : nSet
                
                totalProbeX = [totalProbeX; xing_degX(iProbe) .* ones(nCont1.*nCont2.*nIpsi1,1)];
                totalProbeY = [totalProbeY; xing_degY(iProbe) .* ones(nCont1.*nCont2.*nIpsi1,1)];
                totalProbeIndex = [totalProbeIndex; ProbeIndex(iProbe) .* ones(nCont1.*nCont2.*nIpsi1,1)];
                
                %- select subsets contra-/ipsilateral to current probe
                xing_degX_Cont1 = xing_degX; xing_degY_Cont1 = xing_degY;
                xing_degX_Ipsi1 = xing_degX; xing_degY_Ipsi1 = xing_degY;
                if xing_degX(iProbe) < 0  %- if-clause works only for shapes with vertical symmetry axes
                    iDelCont1  = find(xing_degX<0);
                    iDelIpsi1  = find(xing_degX>0);
                elseif xing_degX(iProbe) > 0
                    iDelCont1  = find(xing_degX>0);
                    iDelIpsi1  = find(xing_degX<0);
                end
                xing_degX_Cont1(iDelCont1) = []; % get rid of any locations on the same side as the probe
                xing_degY_Cont1(iDelCont1) = []; %#ok<*NASGU>
                xing_degX_Ipsi1([iDelIpsi1; iProbe]) = [];
                xing_degY_Ipsi1([iDelIpsi1; iProbe]) = [];
                
                for iCont1 = 1 : nCont1 
                    totalCont1X = [totalCont1X; xing_degX_Cont1(iCont1) .* ones(nCont2.*nIpsi1,1)];
                    totalCont1Y = [totalCont1Y; xing_degY_Cont1(iCont1) .* ones(nCont2.*nIpsi1,1)];
                    %- contralateral distractors minus the first
                    xing_degX_Cont2 = xing_degX_Cont1; xing_degX_Cont2(iCont1) = [];
                    xing_degY_Cont2 = xing_degY_Cont1; xing_degY_Cont2(iCont1) = [];
                    
                    for iCont2 = 1 : nCont2
                        totalCont2X = [totalCont2X; vF.*xing_degX_Cont2(iCont2) .* ones(nIpsi1,1)]; %vF to have NaN entries where no additional stimulus is needed
                        totalCont2Y = [totalCont2Y; vF.*xing_degY_Cont2(iCont2) .* ones(nIpsi1,1)];
                        totalIpsi1X = [totalIpsi1X; vF.*xing_degX_Ipsi1];
                        totalIpsi1Y = [totalIpsi1Y; vF.*xing_degY_Ipsi1];
                    end
                end
            end
            
        end
        
    end % set size
end

totalFPX = zeros(nTrials,1);
totalFPY = zeros(nTrials,1);

%- randomize order
[dummy, iRandomize] = sort(rand(nTrials,1));
totalProbeX = totalProbeX(iRandomize); totalProbeY = totalProbeY(iRandomize);
totalCont1X = totalCont1X(iRandomize); totalCont1Y = totalCont1Y(iRandomize);
totalCont2X = totalCont2X(iRandomize); totalCont2Y = totalCont2Y(iRandomize);
totalIpsi1X = totalIpsi1X(iRandomize); totalIpsi1Y = totalIpsi1Y(iRandomize);
totalsacTargX = totalsacTargX(iRandomize); totalsacTargY = totalsacTargY(iRandomize);
totalFPX = totalFPX(iRandomize); totalFPY = totalFPY(iRandomize);

%-- convert to pixels

totalProbeX = round(totalProbeX .* pixPERdeg) + HC;
totalProbeY = VC - round(totalProbeY .* pixPERdeg);

totalCont1X = round(totalCont1X .* pixPERdeg) + HC;
totalCont1Y = VC - round(totalCont1Y .* pixPERdeg);

totalCont2X = round(totalCont2X .* pixPERdeg) + HC;
totalCont2Y = VC - round(totalCont2Y .* pixPERdeg);

totalIpsi1X = round(totalIpsi1X .* pixPERdeg) + HC;
totalIpsi1Y = VC - round(totalIpsi1Y .* pixPERdeg);

totalsacTargX = round(totalsacTargX .* pixPERdeg) + HC;
totalsacTargY = VC - round(totalsacTargY .* pixPERdeg);

totalFPX = round(totalFPX .* pixPERdeg) + HC;
totalFPY = VC - round(totalFPY .* pixPERdeg);

%-- wrap into 1 struct
trials.totalProbeX = totalProbeX;
trials.totalProbeY = totalProbeY;
trials.totalProbeIndex = totalProbeIndex;
trials.totalCont1X = totalCont1X;
trials.totalCont1Y = totalCont1Y;
trials.totalCont2X = totalCont2X;
trials.totalCont2Y = totalCont2Y;
trials.totalIpsi1X = totalIpsi1X;
trials.totalIpsi1Y = totalIpsi1Y;
trials.totalsacTargX = totalsacTargX;
trials.totalsacTargY = totalsacTargY;
trials.totalFPX = totalFPX;
trials.totalFPY = totalFPY;




end
