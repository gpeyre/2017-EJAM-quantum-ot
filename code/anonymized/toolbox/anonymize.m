function anonymize(rep, ext, str, str_new)

if nargin<4
    str_new = '[anonymized]';
end

a = dir([rep '*' ext]);
for i=1:length(a)
    anonymize_file([rep a(i).name], str, str_new);
end

end


function anonymize_file(name, str, str_new)

if nargin<3
    str_new = '[anonymized]';
end

% read
fid = fopen(name, 'r');
t = {};
t{end+1} = fgets(fid);
while ischar(t{end})
    t{end} = strrep(t{end},str,str_new);
    t{end+1} = fgets(fid);
end
fclose(fid);

% write
fid = fopen(name, 'w');
for i=1:length(t)
    fwrite(fid,t{i});
end
fclose(fid);
    
end