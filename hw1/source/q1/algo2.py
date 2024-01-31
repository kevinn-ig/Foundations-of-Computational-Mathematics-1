import random

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
rem = 0
m = 0
temp_m = 0
max_m = 0
e = 0
u_n = (beta**(1-t))/2

max_m = beta**(t)-1
max_x = beta**U-beta**(U-t)
min_m = beta**(t-1)
min_x = beta**(L-1)

x = random.uniform(-max_m, max_m)
if x > 0:
    s = 1
else:
    s = -1

a = abs(x)

if a > 10**t:
    m = int(a)
    e=L
    while(m>=(10**t)):
        remainder = m%10
        m = m//10
        e = e+1
    m,e = RoundToNearest(remainder, m, e, t)
    m = s*m
else:
    e=U
    while(a<10**(t-1)):
        a = 10*a
        e = e - 1
    remainder = 10*(a-int(a))
    m = int(a)
    m,e = RoundToNearest(remainder, m, e, t)
    m = s*m




print(x, m, e)


