import glob
from modules.Samples import Samples


print("Hi")
s = Samples()
s.load_data(glob.glob('data/*.txt'))
print("end")
s.reduce_to_2d_per_gene()
print("aaa")
# print(s.get_reduced_samples())