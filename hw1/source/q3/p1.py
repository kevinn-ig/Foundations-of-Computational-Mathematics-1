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
error_count = 0
true_error = 0

max_m = beta**(t)-1
max_x = beta**U-beta**(U-t)
min_m = beta**(t-1)
min_x = beta**(L-1)
n=1e7
k=2



e_x = random.randint(-8,0)


if random.choice([True, False]):
    x = random.uniform(-9999, -1000)
else:
    x = random.uniform(1000, 9999)

true_sum = x
approx_sum = 0
norm = 0

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

    if(abs(a+b) > 9999 or abs(a+b) < 1000):
        k += 1
    elif(abs(a) > 9999 or abs(a) < 1000 or abs(b) > 9999 or abs(b) < 1000):
        k +=1
    else:
        c = a+b
        if(abs(approx_sum + c) > 9999 or abs(approx_sum + c) < 1000):
            k+=1
        else:
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
            e_x = e
            x = m
            approx_sum = m
            true_sum = true_sum + y
            norm = norm + abs(x)
            error = abs(approx_sum-true_sum)/((k-1)*norm)
            if(error>u_n):
                error_count+=1
                print("Error > u_n")
            k=k+1
error = abs(approx_sum-true_sum)/((k-1)*norm)
print(error)

