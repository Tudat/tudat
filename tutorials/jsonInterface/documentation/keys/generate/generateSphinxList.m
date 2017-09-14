function generateSphinxList(file)

filenameparts = strsplit(file,'.');
filename = filenameparts{1};

txt = fileread(file);
lines = strsplit(txt,'\n');
cells = {};
possibleValues = {};
children = {};
for i = 1:length(lines)
    if ~isempty(lines{i})
        parts = strtrim(strsplit(lines{i},'|'));
        possibleValues{end+1} = {};
        children{end+1} = {};
        for j = 1:length(parts)
            if j == 4
                cell = parts{j};
                parts{j} = '';
                [~,tok] = regexp(cell,'(m|o):(.+)','match','tokens');
                if ~isempty(tok)
                    cell = tok{1}{1};
                    parts{j} = sprintf(' if :literal:`%s`',tok{1}{2});
                end
                if strcmpi(cell,'m')
                    parts{j} = ['mandatory' parts{j}];
                elseif strcmpi(cell,'o')
                    parts{j} = ['optional' parts{j}];
                elseif strcmpi(cell,'mo')
                    parts{j} = 'mandatory/optional';
                else
                end
            end
            if j == 6
                ccparts = strtrim(split(parts{j},','));
                for k = 1:length(ccparts)
                    cparts = split(ccparts{k},'+');
                    if length(cparts) == 2
                        if isempty(cparts{1})
                            cparts = {cparts{2}};
                        elseif isempty(cparts{2})
                            cparts = {cparts{1}};
                        end
                    end
                    if length(cparts) == 2
                        possibleValues{end}{k} = cparts{1};
                        children{end}{k} = cparts{2};
                    else
                        if ccparts{k}(1) == '+'
                            children{end}{k} = cparts{1};
                        else
                            possibleValues{end}{k} = cparts{1};
                            if ccparts{k}(end) == '+'
                                children{end}{k} = [capitalize(cparts{1}) ' ' filename];
                            end
                        end
                    end
                end
            end
        end
        cells{end+1} = parts;
    end
end

fid = fopen(fullfile('..',[filename '.rst']),'w');

fprintf(fid,'.. role:: arrow\n\n');
for i = 1:length(cells)
    fprintf(fid,'- ');
    fprintf(fid,':literal:`%s` ',cells{i}{1});
    fprintf(fid,':class:`%s` ',cells{i}{2});
    fprintf(fid,'(%s)',cells{i}{4});
    if ~isempty(cells{i}{3})
        fprintf(fid,' %s',cells{i}{3});
    end
    if isempty(cells{i}{3}) || ~strcmp(cells{i}{3}(end),'.')
        fprintf(fid,'.');
    end
    if strcmp(cells{i}{1},'string')
        format = ':literal:`"%s"`';
    else
        format = ':literal:`%s`';
    end
    if ~isempty(possibleValues{i})
        fprintf(fid,' Possible values: ');
        for j = 1:length(possibleValues{i})
            value = possibleValues{i}{j};
            fprintf(fid,format,value);
            if j < length(possibleValues{i})
                fprintf(fid,', ');
            else
                fprintf(fid,'.');
            end
        end
    end
    if length(cells{i}) >= 5
        if ~isempty(cells{i}{5})
            fprintf(fid,' Default value: ');
            fprintf(fid,format,cells{i}{5});
            fprintf(fid,'.');
        end
    end
    for j = 1:length(children{i})
        childName = children{i}{j};
        childFile = strrep(childName,' ','');
        if generateChild([childFile '.txt'])
            fprintf(fid,'\n\n\t.. container:: toggle\n\n\t\t.. container:: header\n\n\t\t\t');
            fprintf(fid,':arrow:`%s`\n\n\t\t.. include:: keys/%s.rst',childName,childFile);
        end
    end
    fprintf(fid,'\n');
end

fclose(fid);

end

function str = capitalize(str)

str = [upper(str(1)) str(2:end)];

end

function success = generateChild(file)

try
    generateSphinxList(file);
    success = true;
catch ME
    if ~strcmp(ME.identifier,'MATLAB:fileread:cannotOpenFile')
        rethrow(ME);
    end
    success = false;
    fprintf('Did not find file: %s\n',file);
end

end

