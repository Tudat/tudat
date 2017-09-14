files = dir('*.txt');
filenames = {files.name};
for i = 1:length(filenames)
    generateSphinxTable(filenames{i});
end

