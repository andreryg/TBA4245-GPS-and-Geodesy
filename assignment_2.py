import numpy as np

l = np.array([100100, 100099, 100101, 100111])

P = np.array([
    [1, 0, 0, 0],
    [0, 2, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
])

print(np.dot(l,P))

#TASK 1
weight_avg = np.sum(np.dot(l,P)) / np.sum(np.diag(P))
print(weight_avg)

corrections = weight_avg - l #v = Ax - f 
print(corrections)

#v = Ax - f, x = actual height
v = np.array([1,1,1,1]).T * weight_avg - l
print(v)
A = np.array([1,1,1,1]).T
A = np.array([[1], [1], [1], [1]])
f = l
vtPv = v.T @ P @ v
print(vtPv)

#TASK 2
print(A.T * P * A)
x_hat = np.linalg.inv(A.T @ P @ A) @ A.T @ P @ f
#print(np.linalg.inv(A.T * P * A) * A.T * P * f)
#print(np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(A.T, P), A)), A.T), P), f))
print(x_hat)

#TASK 3
nabla1 = []
nabla2 = []
nabla3 = []
temp_p = P.copy().astype('float64')
e = np.zeros(A.shape[0])
new_A = np.c_[A,e]

for i,v in enumerate(l):
    temp_p[i][i] = 0
    nabla1.append(v-np.sum(np.dot(l,temp_p)) / np.sum(np.diag(temp_p)))

    temp_p[i][i] = 0.000001
    temp_x_hat = np.linalg.inv(A.T @ temp_p @ A) @ A.T @ temp_p @ f
    nabla2.append(float(v-temp_x_hat))

    new_A[i][-1] = 1
    Q = np.linalg.inv(np.dot(np.dot(new_A.T, P), new_A))
    blunder_cofactor = Q[-1][-1]
    x_hat_2 = np.dot(np.dot(np.dot(Q, new_A.T), P), f)
    nabla3.append(x_hat_2[-1])

    temp_p[i][i] = P[i][i]
    new_A[i][-1] = 0

print("1",nabla1)
print("2",nabla2)
print("3",nabla3)

#TASK 4
R = np.identity(4) - A @ np.linalg.inv(A.T @ P @ A) @ A.T @ P
print(R)

#TASK 5
Qvv = R @ np.linalg.inv(P)
print(Qvv)

#TASK 6
Qnablanabla = 1 / (np.diag(Qvv) @ P**2)
print(Qnablanabla)

#TASK 7
freedom = 3
std_unit_weigth = (1 / (freedom-1)) * (vtPv - (np.square(nabla3) / Qnablanabla))
std_unit_weigth = np.sqrt(std_unit_weigth)
print("7", std_unit_weigth)

#TASK 8
std_gross_error = std_unit_weigth * np.sqrt(Qnablanabla)
print("8",std_gross_error)

#TASK 9
t = abs(nabla3 / std_gross_error)
print(t)

#TASK 10
alpha_tot = 0.05
alpha = 1 - (1 - alpha_tot)**(1/4)
print(alpha)

#TASK 13
intervals = []
t_a = 4.303
for i,v in enumerate(nabla3):
    #print(v-std_gross_error[i]*t_a, v+std_gross_error[i]*t_a)
    intervals.append(max(abs(v-std_gross_error[i]*t_a), abs(v+std_gross_error[i]*t_a)))
print(intervals)

#TASK 14
zero = [0,0,0,0]
external_reliability = []
for i,v in enumerate(intervals):
    zero[i] = v
    external_reliability.append((np.linalg.inv(A.T @ P @ A) @ A.T @ P @ zero)[0])
    zero[i] = 0
print(external_reliability)