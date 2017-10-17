function success = generateRST(filename)
% Call generateRST('root') to generate all .rst files

success = true;
fid = fopen(fullfile('rst',[filename '.rst']),'w');

try
    txt = fileread(fullfile('txt',[filename '.txt']));
    lines = strsplit(txt,'\n');
    cells = {};
    possibleValues = {};
    children = {};
    childrenSameLevel = {};
    for i = 1:length(lines)
        line = lines{i};
        if ~isempty(line)
            if line(1) == '+' || line(1) == '-'
                children{end}{end+1} = strtrim(line(2:end));
                childrenSameLevel{end}{end+1} = line(1) == '-';
            else
                parts = strtrim(strsplit(line,'|'));
                possibleValues{end+1} = {};
                children{end+1} = {};
                childrenSameLevel{end+1} = {};
                for j = 1:length(parts)
                    if j == 4
                        cell = parts{j};
                        parts{j} = '';
                        [~,tok] = regexp(cell,'(m|o):(.+)','match','tokens');
                        if ~isempty(tok)
                            cell = tok{1}{1};
                            condition = tok{1}{2};
                            if condition(1) == '!'
                                parts{j} = sprintf(' if :jsonkey:`%s` undefined',condition(2:end));
                            elseif condition(1) == ':'
                                parts{j} = sprintf(' if :jsonkey:`%s` defined',condition(2:end));
                            else
                                subparts = split(condition,'=');
                                if length(subparts) == 2
                                    parts{j} = sprintf(' if :jsonkey:`%s` set to :literal:`%s`',...
                                        subparts{1},subparts{2});
                                else
                                    parts{j} = sprintf(' %s',condition);
                                end
                            end
                        end
                        if strcmpi(cell,'m')
                            parts{j} = ['mandatory' parts{j}];
                        elseif strcmpi(cell,'o')
                            parts{j} = ['optional' parts{j}];
                        elseif strcmpi(cell,'mo')
                            parts{j} = 'mandatory/optional';
                        else
                        end
                    elseif j == 6
                        possibleValues{end} = strtrim(split(parts{j},','));
                    end
                end
                cells{end+1} = parts;
            end
        end
    end
    
    fprintf(fid,'\n.. role:: jsontype\n.. role:: jsonkey\n.. role:: arrow\n\n');
    for i = 1:length(cells)
        fprintf(fid,'- ');
        fprintf(fid,':jsontype:`%s` ',cells{i}{1});
        fprintf(fid,':jsonkey:`%s` ',cells{i}{2});
        fprintf(fid,'(%s).',cells{i}{4});
        if ~isempty(cells{i}{3})
            fprintf(fid,' %s',cells{i}{3});
            if ~strcmp(cells{i}{3}(end),'.')
                fprintf(fid,'.');
                
            end
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
            childFilename = strrep(strrep(childName,' ',''),'-','');
            if generateRST(childFilename)
                if childrenSameLevel{i}{j}
                    tab = '';
                else
                    tab = sprintf('\t');
                end
                fprintf(fid,'\n\n%s.. container:: toggle\n\n%s\t.. container:: header\n\n%s\t\t',tab,tab,tab);
                fprintf(fid,':arrow:`%s`\n\n%s\t.. include:: %s.rst',childName,tab,childFilename);
            end
        end
        fprintf(fid,'\n');
    end
catch ME
    success = false;
    if strcmp(ME.identifier,'MATLAB:fileread:cannotOpenFile')
        fprintf('Did not find source file: %s\n',filename);
    else
        fprintf('Error when generating list from source file: %s\n',filename);
        rethrow(ME);
    end
end

fclose(fid);

end

