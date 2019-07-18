# D2R算法应用
import d2r_lib
import matplotlib.pyplot as plt

d2rlist = d2r_lib.D2R_in_consecutive_sliding_windows_on_a_genome(3, 9, SEQUENCE = 'GTCAGTACGTCGTCGTA')
x = [i for i in range(len(d2rlist))]
print(d2rlist)
plt.figure()
plt.bar(x, d2rlist)
plt.title('d2r')
plt.show()



