function kx = readPerm(tline)
% Read x absolute permeabilities from the perm file.
% Extract the file name
splitWord = strsplit(tline,"'");
splitWord = strsplit(splitWord{2},".\");
permFileName = splitWord{2};

% Read file
fid2 = fopen(strcat(pwd, '\input\', permFileName), 'r' );
tline = fgetl(fid2);
i=1;
while ischar(tline)
    tline = fgetl(fid2);
    kx(i) = str2double(tline);
    i = i + 1;
end
end

