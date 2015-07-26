#! /usr/bin/python3.3
import matplotlib.pyplot as plt


def f(x, y):
	return (1 - x ** 2) * y + x ** 4 - x ** 2 + 2 * x


def runge_kutt_method(h,n):
	x = [ i * h for i in range(n)]
	y = [0]
	for i in range(len(x)):
		k1 = h * f(x[i],y[i])
		k2 = h * f(x[i] + h / 2, y[i] + k1 / 2)
		k3 = h * f(x[i] + h / 2, y[i] + k2 / 2)
		k4 = h * f(x[i] + h, y[i] + k3)
		y.append(y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
	return list(zip(x, y))


def adams_bashforth_method(h,n):
	runge_kutt_values = runge_kutt_method(h,n)
	x = [runge_kutt_values[i][0] for i in range(len(runge_kutt_values))]
	y = [runge_kutt_values[i][1] for i in range(4)]
	for i in range(3, len(x)):
		y_k_1 = y[i] + h* (55 * f(x[i], y[i]) - 59 * f(x[i - 1], y[i - 1]) + 37 * f(x[i - 2], y[i - 2]) - 9 * f(x[i - 3], y[i - 3])) / 24
		y.append(y_k_1)
	return list(zip(x, y))


def exact_solution(h,n):
	''' for F(x) = x ^ 2 '''
	x = [ i * h for i in range(n)]
	y = [ x[i] ** 2 for i in range(n)]
	return list(zip(x,y))


def runge_error(h,n):
	y_h = [runge_kutt_method(h,n)[i][1] for i in range(n)]
	y_half_h = [runge_kutt_method(h / 2,2 * n)[2 * i][1] for i in range(n)]
	e_array = [abs(y_h[i] - y_half_h[i]) / 15 for i in range(n)]
	return e_array


def adams_error(h,n):	
	y_h = [adams_bashforth_method(h,n)[i][1] for i in range(n)]
	y_half_h = [adams_bashforth_method(h / 2,2 * n)[2 * i][1] for i in range(n)]
	e_array = [abs(y_h[i] - y_half_h[i]) / 15 for i in range(n)]
	return e_array	


def plot_error(h,n, e1, e2):
	x = [ i * h for i in range(n)]

	plt.title('error')
	plt.figure(1)
	plt.subplot(211)
	plt.plot(x, e1, 'bo', linestyle = ":")
	plt.ylabel('runge_error(x)')
	plt.grid(True)
	
	plt.subplot(212)
	plt.plot(x, e2, 'ro', linestyle = ":")
	plt.xlabel('x')
	plt.ylabel('adams_error(x)')
	plt.grid(True)
	plt.show()


def plot_solutions(h,n):
	x = [ i * h for i in range(n)]
	y1 = [runge_kutt_method(h,n)[i][1] for i in range(n)]
	y2 = [adams_bashforth_method(h,n)[i][1] for i in range(n)]
	y_e = [ x[i] ** 2 for i in range(n)]

	plt.title('solutions')
	plt.figure(1)
	plt.subplot(311)
	plt.plot(x, y1, 'bo', linestyle = ":")
	plt.ylabel('runge_kutt')
	plt.grid(True)
	
	plt.subplot(312)
	plt.plot(x, y2, 'ro', linestyle = ":")
	plt.xlabel('x')
	plt.ylabel('adams_bashforth')
	plt.grid(True)

	plt.subplot(313)
	plt.plot(x, y_e, 'go', linestyle = ":")
	plt.xlabel('x')
	plt.ylabel('exact_solution')
	plt.grid(True)
	plt.show()


def main():
	# print("NCM: Assignment #6: Cauchy problem solving \n")
	h = 0.1
	scatter = 20

	# Results
	x = [ i * h for i in range(scatter)]
	y1 = [runge_kutt_method(h,scatter)[i][1] for i in range(scatter)]
	y2 = [adams_bashforth_method(h,scatter)[i][1] for i in range(scatter)]
	y_e = [ x[i] ** 2 for i in range(scatter)]
	# print('x \t', 'f(x)\t', 'Runge-Kutt\t\t', 'Adams-Bashforth\n')
	# for i in range(scatter): print(round(x[i],1), '\t', round(y_e[i],2), '\t', round(y1[i],15), '\t', round(y2[i],15))

	# Error results
	e1_arr = runge_error(h,scatter)
	e1_max = max(e1_arr)
	e2_arr = adams_error(h,scatter)
	e2_max = max(e2_arr)	
	# print('x \t', 'Runge-Kutt error\t', 'Adams-Bashforth error\n')
	# for i in range(1,scatter): print(round(x[i],1), '\t', round(e1_arr[i],15), '\t\t', round(e2_arr[i],15))
	# print('\nmax e:','\t', e1_max, '\t', e2_max)

	# Visualisation 
	# plot_solutions(h, scatter)
	# plot_error(h, scatter, e1_arr, e2_arr)


if __name__ == '__main__':
	main()