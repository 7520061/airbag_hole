function stop = mypsofun(optimValues,state)
global DIRECTORY
% INFOMATION = optimValues.iteration;
    stop = false;
    switch state
        case 'iter'
            dirn  = dir([DIRECTORY '/iter*']);
            n = length(dirn);
            mkdir([DIRECTORY '/iter=' num2str(n)]);
        case 'interrupt'

        case 'init'
            disp('Iter f(x)');
            dirn  = dir([DIRECTORY '/iter*']);
            n = length(dirn);
            mkdir([DIRECTORY '/iter=' num2str(n)]);
        case 'done'
            stop = true;
        otherwise
    end
