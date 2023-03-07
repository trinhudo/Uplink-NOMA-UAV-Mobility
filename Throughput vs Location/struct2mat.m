function out = struct2mat(x)
    out = cell2mat( struct2cell( x ) );
end