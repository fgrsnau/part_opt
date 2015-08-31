function theta = get_theta(energy)

f2 = double(get_f2_full(energy));
theta = [double(energy.f1(:)); f2(:)]';

end