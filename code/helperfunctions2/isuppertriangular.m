function boolean = isuppertriangular(D)

  boolean = all(all(tril(D) == 0));
end