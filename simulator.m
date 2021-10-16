function [fluid_, grid_, connection_, equation_, well_, output_, frac_, runTime]...
    = simulator(filename)
tic
% Initialize simulator variables
[fluid_, grid_, connection_, equation_, well_, output_, frac_] = initialize(filename);
iTime = 1; Pwf=0; qw=0; qo=0; isWCMOK = true; cumo = 0;
while equation_.curT < equation_.totT   
    % Create initial NR guess.
    fluid_.vecP = fluid_.vecPPrev;
    fluid_.Po = fluid_.PoPrev;
    fluid_.Sw = fluid_.SwPrev;  
    
    % Run NR loop
    [isConv, fluid_.vecP, fluid_.Po, fluid_.Sw, fluid_.bo, ...
        fluid_.bw, grid_.poro, well_, nIter, q_out, maxR_temp,...
        maxS_temp, maxP_temp] ...
        = equation_.solveNRProblem(fluid_, grid_, connection_, well_);
    output_.numIter(iTime) = nIter;
    
    if isConv
        % If WCM violated, Discard P/S, update WCM
        [isWCMOK, well_.wellMode, Pwf_temp, qw_temp, qo_temp, cumo_temp]...
            = well_.updateCheckWellCond(fluid_, equation_.delT); 
        % Case 1: converged and no violation on WCM
        if isWCMOK    
            % Save necessary output for each time step
            if iTime == 1
                Pwf =Pwf_temp; qw=qw_temp; qo=qo_temp; cumo = cumo_temp;
                cumo=cumo_temp; CFL = q_out*equation_.delT./(grid_.poro.*grid_.volCell/5.615);
                SWAT = fluid_.Sw; PORO = grid_.poro; PRESSURE = fluid_.Po;
                dT_out = equation_.delT;
                CFLmax = max(CFL);
                maxR = maxR_temp; maxS = maxS_temp; maxP = maxP_temp;
            else
                Pwf(:,iTime)=Pwf_temp; qw(:,iTime)=qw_temp; qo(:,iTime)=qo_temp;
                cumo(:,iTime)=cumo_temp+cumo(:,iTime-1);
                CFL(floor(equation_.curT/(equation_.totT/4))+2,:) = q_out*equation_.delT./(grid_.poro.*grid_.volCell/5.615);
                SWAT(floor(equation_.curT/(equation_.totT/4))+2,:) = fluid_.Sw;
                PORO(floor(equation_.curT/(equation_.totT/4))+2,:) = grid_.poro;
                PRESSURE(floor(equation_.curT/(equation_.totT/4))+2,:) = fluid_.Po;
                dT_out(iTime) = equation_.delT;                
                CFLmax(iTime) = max(CFL(floor(equation_.curT/(equation_.totT/4))+2,:));
                maxR(iTime) = maxR_temp; maxS(iTime) = maxS_temp; maxP(iTime) = maxP_temp;
            end
            % Update P/S
            equation_.curT = equation_.curT+equation_.delT;
            output_.simTime(iTime) = equation_.curT;
            
            % Time stepping
            equation_.delT = equation_.updateTimeStep(fluid_, isWCMOK); 
            [fluid_.vecPPrev, fluid_.PoPrev, fluid_.SwPrev,...
                fluid_.bwPrev, fluid_.boPrev, fluid_.PoPrevIter, ...
                fluid_.SwPrevIter] = fluid_.backUpFluid();
            [grid_.poroPrev] = grid_.backUpGrid();
            iTime = iTime+1;
        % Case 2: converged but violation on WCM
        else
            equation_.delT = equation_.delT/2;
        end
    % Case 3: Not converged
    else
        equation_.delT = equation_.delT/2;
    end
end
% Save all outputs in output class
[output_.WOPR, output_.WWPR, output_.WWIR, output_.WBHPPROD, output_.WBHPINJ...
    , output_.FOPT, output_.SWAT, output_.PORO, output_.PRESSURE, output_.CFL] =...
    output_.updateOutput(Pwf, qw, qo, cumo, well_, SWAT, PORO, PRESSURE, CFL);
output_.dT_out = dT_out;
output_.CFLmax = CFLmax;
output_.maxR = maxR; output_.maxS = maxS; output_.maxP = maxP;
runTime = toc;
end