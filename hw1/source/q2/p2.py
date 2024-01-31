import random
import matplotlib.pyplot as plt

def RoundToNearest(remainder, m, e, t):
    if(remainder >= 5):
        m = m+1
        if m == 10**t:
            m = m//10
            e = e + 1
    return m,e



beta = 10 # base
t = 4 #significant figures
L = -8
U = 0
m = 0
e = 0
m_list = []
e_list = []
u_n = (beta**(1-t))/2

max_m = beta**(t)-1
max_x = beta**U-beta**(U-t)
min_m = beta**(t-1)
min_x = beta**(L-1)
n=1e5
n_list = []
k=0
error_count = 0
error_list = []

while k<n:
    remainder = 0
    e_x = random.randint(-8,0)
    e_y = random.randint(-8,0)
    if random.choice([True, False]):
        x = random.uniform(-9999, -1000)
    else:
        x = random.uniform(1000, 9999)

    if random.choice([True, False]):
        y = random.uniform(-9999, -1000)
    else:
        y = random.uniform(1000, 9999)

    shift = e_x-e_y
    if shift>0:
        a = x*10**(shift)
        b = y
        e = e_y
    else:
        a = y*10**(-1*shift)
        b = x
        e = e_x
    
    if(abs(a+b) > 9999 or abs(a+b) < 1000):
        k += 1
    elif(abs(a) > 9999 or abs(a) < 1000 or abs(b) > 9999 or abs(b) < 1000):
        k +=1
    
    else:
        c = a+b
        if c>=0:
            s = 1
        else:
            s = -1

        c = abs(c)

        if c>= 10**t:
            while(c>=10**t):
                remainder = c%10
                c = c//10
                e = e+1
        else:
            while(c<10**(t-1)):
                c = c*10
                e = e-1
        remainder = 10*(c-int(c))
        m = int(c)
        m,e = RoundToNearest(remainder,m,e,t)
        m = s*m
        
        error = abs(((x+y) - (m))/(x+y))
        if(error > u_n):
            error_count += 1
        error_list.append(error)
        n_list.append(k)
        k += 1
        

print(error_count)

mean = sum(error_list)/(len(error_list))
squared_diff = [(x-mean)**2 for x in error_list]
variance = sum(squared_diff)/len(error_list)

print(mean, variance, max(error_list), min(error_list))



plt.scatter(n_list, error_list, label='Values', marker='o')
plt.axhline(y=u_n, color='green', linestyle=':', label=f'u_n at {u_n}')



plt.xlabel('n')
plt.ylabel('error')
plt.title('Errors Over n')
plt.legend(loc = "lower left")
plt.grid(True)
plt.savefig('/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw1/figures/line_graph.png')

plt.figure()

plt.hist(error_list)
plt.xlabel('Error')
plt.ylabel('Frequency')
plt.title('Error Histogram')
plt.grid(True)
plt.savefig('/Users/kevin_smith/Desktop/FSU_Relevant_Stuff/fall_2023/FCM1/homework/programming/hw1/figures/histogram.png')

    








