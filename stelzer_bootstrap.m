function bootstrappedPerm = stelzer_bootstrap(perm)
% create Stelzer-bootstrapped permutation distribution from standard
% permutation data (in the form random seed x participant).
    bootstrappedPerm = zeros(10000,1);
    for p = 1:10000
        i = randi(size(perm,1), 1, size(perm,2));
        samples = perm(sub2ind(size(perm), i, 1:size(perm,2)));
        bootstrappedPerm(p) = mean(samples);
    end
    bootstrappedPerm = sort(bootstrappedPerm);
end