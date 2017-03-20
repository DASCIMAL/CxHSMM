function [val] = zlog(arg)

% just like the log function, but it doesn't issue any warning for log of
% zero. Simply returns -Inf in that case. 

if (arg ~= 0)
	val = log(arg);
else
	val = -Inf;
end
