function u = read_sequence(filename)
i = 0;
ok = true;
while ok
    i = i+1;
    try
        u(:,:,i) = double(imread(filename,i));
    catch
        i = i-1;
        ok = false;
    end
end
fprintf("%d images read from file\n%s\n",i,filename);
