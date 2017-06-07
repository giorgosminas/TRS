function id = findIndex(strToFind,allStrings)
id = -1;
found = false;
counter = 1;
maxStrings = length(allStrings);
while counter <= maxStrings && ~found
    if sum(strcmpi(strToFind,allStrings{counter}) > 0)
        %found catma
        id = counter;
        return;
    else
        counter = counter + 1;
    end
end
end