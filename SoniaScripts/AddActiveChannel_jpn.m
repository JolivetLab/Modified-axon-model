function par = AddActiveChannel(par, activeChannel)
%ADDACTIVECHANNEL - Adds an active channel to the nodes of a myelinated axon
%   par - ADDACTIVECHANNEL(par, channel)
%       Outputs:
%           par -           modified parameter structure
%       Inputs:
%           par -           existing parameter structure
%           channel -       structure containing parameters of an active
%                           channel or character string pointing to file containing a structure
%                           containing parameters of an active channel.
%
%   If `channel' specifies a filename, the active channel structure in the file must
%   be called `activeChannel'.

% Find where to put the active channel.
if isfield(par.intn, 'elec')
    if isfield(par.intn.elec, 'act')
        idx =                                                   length(par.intn.elec.act)+1;
    else
        idx =                                                   1;
    end
else
    idx =                                                       1;
end

% If supplied `activeChannel' is a structure, attempt to put it in the
% main `par' structure. Otherwise, we expect it to be a filevname string.
if isstruct(activeChannel)
    
    par.intn.elec.act(idx) =                                activeChannel;
    par.intn.elec.act(idx).cond.value.vec =                 par.intn.elec.act(idx).cond.value.ref * ones(par.geo.nintn, 4);
    
elseif ischar(activeChannel)
    
    fileContents =                                          load(activeChannel);
    if ~isfield(fileContents, 'activeChannel')
        error('There is no `activeChannel'' variable in this file');
    end
    par.intn.elec.act(idx) =                                fileContents.activeChannel;
    par.intn.elec.act(idx).cond.value.vec =                 par.intn.elec.act(idx).cond.value.ref * ones(par.geo.nintn, 4);
end