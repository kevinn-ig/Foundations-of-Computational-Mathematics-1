import random
def RoundToNearest(remainder, m, e, t):
    if(remainder >= 5):
        m = m+1
        if m == 10**t:
            m = m//10
            e = e + 1
    return m,e


beta = 2 # base
t = 53 #significant figures
L = -1021
U = 1024
u_n = (beta**(1-t))/2
x_list = []
p_list = []
p_list2 = []
p_list = []
n=100000
max_m = beta**(t)-1
max_x = beta**U-beta**(U-t)
min_m = beta**(t-1)
min_x = beta**(L-1)

for i in range(0,101):
    x = random.uniform(-max_m,max_m)
    x_list.append(x)

for i in range(0,101):
    p = random.uniform(-x_list[i]*u_n, x_list[i]*u_n)
    p_list.append(p)

approx_cond = max(p_list)

xSumP = sum([x + y for x, y in zip(x_list, p_list)])
abs_x_list = [abs(x) for x in x_list]
abs_p_list = [abs(p) for p in p_list]


cond_num1 = (abs((sum(x_list) - xSumP))/(abs(sum(x_list)))*(sum(abs_x_list)))/(sum(abs_p_list))


for i in range(0, len(x_list)):
    x = x_list[i]

    if x > 0:
        s = 1
    else:
        s = -1

    a = abs(x)

    if a > beta**t:
        m = int(a)
        e=L
        while(m>=(beta**t)):
            remainder = m%beta
            m = m//beta
            e = e+1
        m,e = RoundToNearest(remainder, m, e, t)
        m = s*m
    else:
        m = int(a)
        e=U
        while(a<beta**(t-1)):
            a = beta*a
            e = e - 1
        remainder = beta*(a-int(a))
        m = int(a)
        m,e = RoundToNearest(remainder, m, e, t)
        m = s*m
    x_list[i] = m

for i in range(0,101):
    p = random.uniform(-x_list[i]*u_n, x_list[i]*u_n)
    p_list2.append(p)

xSumP = sum([x + y for x, y in zip(x_list, p_list2)])
abs_x_list = [abs(x) for x in x_list]
abs_p_list = [abs(p) for p in p_list2]


cond_num2 = (abs((sum(x_list) - xSumP))/abs(sum(x_list)))*((sum(abs_x_list))/(sum(abs_p_list)))
print(approx_cond)

relative_error = abs(cond_num1 - cond_num2)/abs(cond_num2)
print(relative_error)
