import time
import matplotlib.pyplot as plt
from calculus import main_calculus
from calculus_threaded import main_calculus_threaded

tasks1 = [[305, 110, 599], [304, 201, 401], [2513, 870, 4759], [35768, 34528, 56501], [320774, 430278, 488639], [190177, 261140, 665359], [1557464, 5009176, 12369059], [167357551, 302646126, 433096217]]

times_main = []
times_main_threaded = []

for task in tasks1:
    p, g, h = task

    start_time = time.perf_counter()
    main_calculus(p, g, h)
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    times_main.append(elapsed_time)
    print(elapsed_time)

    start_time = time.perf_counter()
    main_calculus_threaded(p, g, h)
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    times_main_threaded.append(elapsed_time)
    print(elapsed_time)

    if elapsed_time > 180:
        break

print("Час роботи алгоритму main:", times_main)
print("Час роботи алгоритму main_threaded:", times_main_threaded)


p_values = range(1, len(times_main) + 1)


plt.plot(p_values, times_main, label="main")
plt.plot(p_values, times_main_threaded, label="main_threaded")
plt.xlabel("Порядок простого числа p")
plt.ylabel("Час роботи (сек)")
plt.title("Залежність часу роботи від порядку простого числа p")
plt.legend()
plt.show()




