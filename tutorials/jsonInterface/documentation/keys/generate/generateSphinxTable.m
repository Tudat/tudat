function generateSphinxTable(file)

txt = fileread(file);
lines = strsplit(txt,'\n');
cells = {{'Value type','Key','Description','M/O','Default'}};
rows = 0;
cols = 0;
colsWidth = [];
for i = 1:length(cells{1})
    colsWidth(i) = length(cells{1}{i});
end
for i = 1:length(lines)
    if ~isempty(lines{i})
        rows = rows + 1;
        parts = strtrim(strsplit(lines{i},'|'));
        p = length(parts);
        cols = max(cols,p);
        for j = 1:p
            cell = parts{j};
            if j == 3
                parts{j} = regexprep(cell,'\[(.+?)\]\((.+?)\)','`$1 <#$2>`_');
            elseif j == 4
                if cell == 'm'
                    parts{j} = 'M';
                elseif cell == 'o'
                    parts{j} = 'O';
                end
            end
            colsWidth(j) = max(colsWidth(j),length(parts{j}));
        end
        cells{end+1} = parts;
    end
end

colsWidth = colsWidth + 2;

filename = strsplit(file,'.');
fid = fopen(fullfile('..',[filename{1} '.rst']),'w');

hline(fid,colsWidth,'-');
row(fid,colsWidth,cells{1});
hline(fid,colsWidth,'=');
for i = 2:length(cells)
    row(fid,colsWidth,cells{i});
    hline(fid,colsWidth,'-');
end
fclose(fid);

end


function hline(fid,colsWidth,c)

fprintf(fid,'+');
for i = 1:length(colsWidth)
    fprintf(fid,repmat(c,1,colsWidth(i)));
    fprintf(fid,'+');
end
fprintf(fid,'\n');

end


function row(fid,colsWidth,cells)

fprintf(fid,'|');
for i = 1:length(cells)
    fprintf(fid,' ');
    fprintf(fid,cells{i});
    for j = (length(cells{i})+2):colsWidth(i)
        fprintf(fid,' ');
    end
    fprintf(fid,'|');
end
fprintf(fid,'\n');

end
