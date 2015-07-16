% parse output file
function signal = parse_output()
fid = fopen('ffastOutput.txt', 'r');
signal = [];
tline = fgets(fid);
while ischar(tline)
    p1 = find(tline == ':');
    binIdx = str2double(tline(1:p1-1));
    tline = tline(p1+2:end-2);
    p1 = find(tline == ',');
    rpart = str2double(tline(1:p1-1));
    ipart = str2double(tline(p1+1:end));
    v = rpart + 1i*ipart;
%     if ipart < 0
%         fprintf('(%d, %.3f - i%.3f)\n',binIdx, rpart,-ipart);
%     else
%         fprintf('(%d, %.3f + i%.3f)\n',binIdx, rpart,ipart);
%     end
    signal = [signal; binIdx v];
    tline = fgets(fid);
end
fclose(fid);

% system('/Users/orhanocal/Documents/Berkeley/Research/FFAST-customized/exe/ffast -f infile.txt -k 4 -s 40');
