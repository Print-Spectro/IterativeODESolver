# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 17:44:33 2020

@author: franc
"""
from matplotlib.pyplot import *

from numpy import *

def minimise(x0, x1, fx, lambda_, tolerance):
    """Minimises functions which can't be algebraically differentiated using a secant line
    Inputs: x0 and x1 are two starting values for the secant line used to approximate dy/dx
            fx is the function in question
            lambda_ is the step size modifier
            tolerance is the value of dy/dx at which the function will stop"""
    """Output: The x value at which the function is minimised"""    
    dfdx = lambda fx, x0, x1: (fx(x1) - fx(x0))/(x1 - x0) #defining gradient of secant line
    while abs(dfdx(fx, x0, x1)) > tolerance:
        x2 = x1 - lambda_*dfdx(fx, x0, x1) #finding the next term in the itteration
        x0 = x1 #xold = xn for next loop
        x1 = x2 #xn = xnew for next loop
    return x2

trapesium = lambda a, dt: 0.5*dt*(a[0]+a[-1]) + sum(a[1:-1]*dt)
"""This trapesium rule function is very simple; it just takes an array of a function and the step size between the discrete
points in the array and returns an approximate integral using numpy arrays. This was adapted from the above trapseium_n
function because applying an itterative method to an array doesn't apply the function to each individual value in the array 
unfortunately"""    
    
secantterm = lambda x0, x1, fx: x1 - (x1 - x0)*fx(x1)/(fx(x1) - fx(x0))
"""Calculates the next term in the secant series given two inputs and the function"""

def secant(x0, x1, fx):
    """
    Finds the roots of an equation. Tn this case it's being used to find the next term in the reverse euler method
    which approximates ODEs implicitely
    Inputs: two starting values for the secant line x0 and x1 and the function fx
    Outputs: a solution to fx = 0"""
    x = x0
    while not (abs(fx(x)) < tolerance and abs(x0 - x1) < tolerance):
        x = secantterm(x0, x1, fx)
        x0 = x1
        x1 = x
    return x1 
 
def Newton4Euler(x0, fxt, dfxtdt, t, dt, F):
    """Solves for the next unknown value x1 in the rev_euler function, given the previous value of our ODE in order for the
    reverse Euler method to work"""
    x1 = x0 - (x0 - x0 - fxt(x0, t)*dt)/(1-dfxtdt(x0, t)*dt) #initalise x1, using the previous value x0 as our starting point
    while abs(F) > 1e-10: #while loop exits when our function is close to zero 
        F = x1 - x0 - fxt(x1, t)*dt #defining the top of the fraction in newton's itteration
        x1 = x1 - F/(1-dfxtdt(x1, t)*dt) #finding the next term
    return x1 #returns the value of x1 once the loop exits. 

#Newton4Euler = lambda x0, fxt, dfxtdt, t, dt, F: 

rev_Euler = lambda x0, x1, fxt, t, dt: x0 + fxt(x1, t)*dt #defining each step in the reverse euler   

def reverse_euler(x0, fxt, dt, n):
    """
    Takes the value of x (x0) at t=0 and itterates to find the approximate solution to an ODE given by dx/dt = f(x, t)
    between t = 0 and t = 1.
    Inputs: x0: The value of x at t = 0. 
            fxt: This is f(x, t) in dx/dt, 
            dt: The step size, smaller dt increases resolution
    Output: Array containing all of the values of the solutions to the ODE between t = 0 and t = 1"""
    out = [[],[]]
    for i in range(n):
        t = i*dt
        out[0].append(x0) #append outputs
        out[1].append(t) #append corresponding times 
        x0 = Newton4Euler(x0, fxt, dfdx, t, dt, 1) #the next value in the iteration is just the newton's method root
        #print(i) #debugging here
    return array(out[0])

reverse_euler = lambda x, fxt, dt, n: (lambda a, fxt, dt, n, : ([a.append(Newton4Euler(a[-1], fxt, dfdx, i*dt, dt, 1)) for i in range(n-1)], array(a)))([x], fxt, dt, n)[1]

T = 1
n = 1000
dt = T/n
t = linspace(0,T,n)
#x0 = linspace(-1, 1, 100)
tolerance = 1e-6
zt = lambda t: t**2 #defining z(t) so that I can test my code
fxt = lambda x, t: x #defining f(x, t) in dx/dt = f(x, t) so that I can test my code 
dfdx = lambda x, t: 1 #defining the differential of f(x,t) in in the ODE



print(minimise(1,1.1, lambda x: trapesium((reverse_euler(x, fxt, dt, n)-zt(t))**2, dt), 0.1, tolerance))
