%
%PAL_PFHB_buildStan  Issues OS-appropriate directive to build executable
%   Stan model
%
%   syntax: [status, OSsays, syscmd] = PAL_PFHB_buildStan(engine,machine)
%
%Internal function
%
%Introduced: Palamedes version 1.10.0 (NP)
%Modified: Palamedes version 1.10.1, 1.11.6, 1.11.8, 1.11.10 
% (see History.m)

function [status, OSsays, syscmd] = PAL_PFHB_buildStan(engine,machine)

if strcmpi(machine.machine,'PCWIN64')
    dirout = engine.dirout;
    dirout(dirout == '\') = '/';
    syscmd = ['cd ',engine.path,' && ','make ',dirout,'/stanModel.exe'];
    [status, OSsays] = system(syscmd);
    if strcmp(engine.engine,'stan') && engine.recyclestan && ~exist('stanModel.exe','file')
        copyfile([engine.dirout,filesep,'stanModel.exe']);
    end
else
    syscmd = ['cd ',engine.path,char(10),'make ',engine.dirout,'/stanModel',char(10)];
    [trash sysinf] = unix('uname -a');
    if PAL_contains(sysinf,'ARM64') & strcmp(machine.environment, 'matlab') & isMATLABReleaseOlderThan("R2023b")
        message = [char(10),char(10),'Actually, ignore the last statement and read on ...',char(10)];
        message = [message,'You are running this program on a MAC with a Silicon chip AND you selected Stan as your MCMC engine ',char(10)];
        message = [message,'AND you are running this program in a Matlab Release less than R2023b ',char(10)];
        message = [message, 'We can still (try to) make this happen but because of a weird incompatibility this will ',char(10)];
        message = [message, 'require some intervention from you. If you''re up for this the program will pause, then ask you to copy some code ',char(10)];
        message = [message, 'and drop it into a MAC OS terminal before starting this program back up.',char(10),char(10)];
        message = [message 'Want to try this (y/n)? (there are other options, see here: ',char(10),'https://www.palamedestoolbox.org/forum/viewtopic.php?t=16) '];
        manualStanBuild = input(message,"s") == 'y';
        if ~manualStanBuild
            disp('About to crash, I think ... ');
            pause(1);
            [status, OSsays] = system(syscmd);
        else
            message = [char(10),'Copy and paste the following lines into a MAC OS terminal (not the MATLAB terminal!), then wait till ',char(10)];
            message = [message,'the process in the MAC OS terminal completes, then return here and type ''y'' (no quotes) followed by <enter>.',char(10)];
            message = [message,'The lines to copy and paste into MAC OS terminal (highlight both lines and copy and paste as a whole, followed by <enter>): ', char(10),char(10)];
            message = [message,syscmd];
            message = [message,char(10),'Type ''y'' <enter> AFTER process in MAC OS terminal completes. '];
            input(message,"s");
            status = 0;
            OSsays = 'MAC Silicon manual workaround';
        end
    else
        [status, OSsays] = system(syscmd);
    end
    if strcmp(engine.engine,'stan') && engine.recyclestan && ~exist('stanModel','file')
        copyfile([engine.dirout,filesep,'stanModel']);
    end
end
if status ~= 0
    message = ['Building Stan executable failed. Palamedes issued the command: ',char(10), syscmd, char(10), 'to your OS and your OS said: ',char(10), OSsays, char(10)];
    message = [message, 'First thing to do is to try and find out whether CmdStan is in working order:',char(10)];
    message = [message, 'Q1: Does PAL_PFHB_SinglePF_Demo (in PalamedesDemos folder) complete without error when you select Stan?', char(10)];
    message = [message, 'If you answered ''yes'' to Q1: This is possibly a Palamedes bug. Send us an e-mail: palamedes@palamedestoolbox.org.', char(10)];
    message = [message, 'If you answered ''no'' to Q1, move to Q2.', char(10)];
    message = [message, 'Q2: Can you get the CmdStan bernoulli example (see CmdStan User''s guide) to work? (if you have multiple versions of CmdStan use the same version Palamedes is trying to use): ', char(10)];
    message = [message, 'If you answered ''yes'' to Q2: This is possibly a Palamedes bug. Send us an e-mail: palamedes@palamedestoolbox.org.', char(10)];
    message = [message, 'If you answered ''no'' to Q2: Get bernoulli example to work, then try again.'];
    error('PALAMEDES:StanBuildFail',message);
end
