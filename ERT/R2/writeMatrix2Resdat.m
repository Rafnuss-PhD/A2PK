function writeMatrix2Resdat(d)
x=d.xx(2:end)-d.xx(1:end-1);
y=d.yy(2:end)-d.yy(1:end-1);
assert(length(y)==size(d.g_true,1),'The input matrix must be of the same size as grid.x')
assert(length(x)==size(d.g_true,2),'The input matrix must be of the same size as grid.y')

[X,Y]=meshgrid(x,y);

fileID = fopen([d.filepath 'g_true'],'w');
fprintf(fileID,'%e   %e   %e   %e\n', [X(:), Y(:), d.g_true(:), log10(d.g_true(:))]);
fclose(fileID);
end