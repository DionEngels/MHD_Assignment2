function [] = PrintResults(res)
%PRINTRESULTS Prints a summary of all the resulting parameters
names = fieldnames(res);

for k=1:numel(names)
    try
        fprintf("%s : %d\n", names{k}, res.(names{k}))
    catch
        fprintf("%s\n", char(res.(names{k})))
    end
end
end

