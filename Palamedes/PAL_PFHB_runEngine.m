%
%PAL_PFHB_runEngine  Issue OS and engine (stan, JAGS) appropriate command 
%   to OS to start MCMC sampling and wait for sampling to be finished.
%   
%   syntax: [status, OSsays, syscmd] = PAL_PFHB_runEngine(engine,machine)
%
%Internal Function
%
% Introduced: Palamedes version 1.10.0 (NP)
% Modified: Palamedes version 1.10.1, 1.11.4, 1.11.6, 1.11.9 
% (See History.m)

function [status, ossays, syscmd] = PAL_PFHB_runEngine(engine,machine)

if strcmp(engine.engine,'stan')
    if strcmpi(machine.machine,'PCWIN64')
        if engine.recyclestan
            enginefolder = cd;
        else
            enginefolder = engine.dirout;
        end
        for chain = 1:engine.nchains            
            Syscmd = ['"',enginefolder,'\stanModel','" ',' sample num_samples=',int2str(engine.nsamples),' random seed=',int2str(engine.seed),' id=',int2str(chain),' data file="',engine.dirout,'\data_Rdump.R" init="',engine.dirout,'\Init_',int2str(chain),'.R" output file="',engine.dirout,'\samples',int2str(chain),'.csv"'];
            if engine.parallel
                Syscmd = [Syscmd, '  && exit &'];
            end
            [status, OSsays] = system(Syscmd);
            if status ~= 0
                message = ['Execution of ',upper(engine.engine), ' failed. Palamedes issued the command: ',char(10), Syscmd, char(10), 'to your OS and your OS (or ',upper(engine.engine),') said: ',char(10), OSsays, char(10)];
                error('PALAMEDES:SamplerExecuteFail',strrep(message,'\','\\'));
            end
            ossays(chain) = cellstr(OSsays);
            syscmd(chain) = cellstr(Syscmd);
            if status ~= 0
                return;
            end
        end
    else    
        if engine.recyclestan
            enginefolder = cd;
        else
            enginefolder = engine.dirout;
        end
        for chain = 1:engine.nchains
            Syscmd = [enginefolder,'/stanModel sample num_samples=',int2str(engine.nsamples),' random seed=',int2str(engine.seed),' id=',int2str(chain),' data file=',engine.dirout,'/data_Rdump.R init=',engine.dirout,'/Init_',int2str(chain),'.R  output file=',engine.dirout,'/samples',int2str(chain),'.csv'];
            if engine.parallel
                Syscmd = [Syscmd, ' &'];
            end
            [status, OSsays] = system(Syscmd);
            if status ~= 0
                message = ['Execution of ',upper(engine.engine), ' failed. Palamedes issued the command: ',char(10), Syscmd, char(10), 'to your OS and your OS (or ',upper(engine.engine),') said: ',char(10), OSsays, char(10)];
                error('PALAMEDES:SamplerExecuteFail',message);
            end
            ossays(chain) = cellstr(OSsays);
            syscmd(chain) = cellstr(Syscmd);
            if status ~= 0
                return;
            end
        end        
    end        
end

if strcmp(engine.engine,'jags')
    for chain = 1:engine.nchains
        if strcmpi(machine.machine,'PCWIN64')
            Syscmd = ['"',engine.path,'\jags','" "',engine.dirout, '\jagsScript',int2str(chain),'.cmd"'];            
            if engine.parallel
                Syscmd = [Syscmd, '  && exit &'];
            end            
        else
            Syscmd = [engine.path,'/jags ',engine.dirout, '/jagsScript',int2str(chain),'.cmd'];
            if engine.parallel
                Syscmd = [Syscmd, ' &'];
            end            
        end
        [status, OSsays] = system(Syscmd);
        if status ~= 0
            message = ['Execution of ',upper(engine.engine), ' failed. Palamedes issued the command: ',char(10), Syscmd, char(10), 'to your OS and your OS (or ',upper(engine.engine),') said: ',char(10), OSsays, char(10)];
            error('PALAMEDES:SamplerExecuteFail',strrep(message,'\','\\'));
        end
        ossays(chain) = cellstr(OSsays);
        syscmd(chain) = cellstr(Syscmd);
        if status ~= 0
            return;
        end
    end   
end
if engine.parallel
    PAL_PFHB_engineFinished(engine);
end