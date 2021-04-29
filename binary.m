function B = binary (value)

    if value > 1.000e-323
       B = 1;
      return
    else
       B = 0;
       return
    end
end