from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count

n_cpu = cpu_count()
print(n_cpu)