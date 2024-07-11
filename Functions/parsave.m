
% simple function for allowing saving within a parfor loop. 
function parsave(fname, x)
save(fname, 'x');
end