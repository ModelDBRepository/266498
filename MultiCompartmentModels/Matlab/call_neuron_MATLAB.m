function call_neuron_MATLAB(file,ParticleID,varargin)

options = '';
if nargin > 2 %Build options using input variables
    for i = 1:length(varargin)
        if isnumeric(varargin{i})
            options = [options sprintf(' -c %s=%d',inputname(i+2),varargin{i})];
        else
            options = [options sprintf(' -c "{sprint(%s,\\"%s\\")}"',inputname(i+1),varargin{i})];
        end        
    end
end
options = [options sprintf(' -c ParticleID=%d',ParticleID)];
options = [options ' -c "run_nrn_script()" -c "quit()"'];

if ispc
    system(['C:/nrn/bin/nrniv.exe -nobanner ' file options]);% change this to your NEURON installation
elseif isunix
    disp(['./special -nobanner ' file options])
    system(['./special -nobanner ' file options]);
else
    error('Neither PC nor Unix');
end
% ids=runNrn(parameters,ParticleID);
% if status~=0
%     error(['NEURON error: status ' num2str(status) 'result ' result]);
% end