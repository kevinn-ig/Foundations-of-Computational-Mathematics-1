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
m = 0
e = 0
m_list = []
e_list = []
u_n = (beta**(1-t))/2
remainder = 0

max_m = beta**(t)-1
max_x = beta**U-beta**(U-t)
min_m = beta**(t-1)
min_x = beta**(L-1)
n=1000
k=1



e_x = random.randint(-8,0)


if random.choice([True, False]):
    x = random.uniform(-9999, -1000)
else:
    x = random.uniform(1000, 9999)



e_list.append(e_x)
while k<n:
    e_y = random.randint(-8,0)
    if random.choice([True, False]):
        y = random.uniform(-9999, -1000)
    else:
        y = random.uniform(1000, 9999)

    shift = e_x-e_y
    if shift>0:
        a = x*10**shift
        b = y
        e = e_y
    else:
        a = y*10**(-shift)
        b = x
        e = e_x

    c = a+b
    if c>=0:
        s = 1
    else:
        s = -1

    c = abs(c)

    if c>= 10**t:
        e = L
        while(c>=10**t):
            remainder = c%10
            c = c//10
            e = e+1
    else:
        e=U
        while(c<10**(t-1)):
            c = c*10
            e = e-1
    m = int(c)
    m,e = RoundToNearest(remainder,m,e,t)
    m = s*m
    e_x = e
    x= m
    k=k+1

print(m,e)

