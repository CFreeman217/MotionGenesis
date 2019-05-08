# -*- coding: utf-8 -*-
import os, math
import numpy as np
import matplotlib.pyplot as plt

# import turtle as turtle
# # from turtle import *
# import tkinter
# from tkinter import *

rad2deg = lambda x : [i*180.0/math.pi for i in x]

# def drawSpring(length):
#         section = length/3
#         seg = section/6
#         turtle.down(); turtle.forward(section); turtle.right(60); turtle.forward(seg)
#         for i in range(2):
#                 turtle.left(120); turtle.forward(seg)
#                 turtle.forward(seg); turtle.right(120); turtle.forward(seg)
#                 turtle.forward(seg)
#         turtle.left(120); turtle.forward(2*seg); turtle.right(120)
#         turtle.forward(seg); turtle.left(60)
#         turtle.forward(section)
#         turtle.pu()

for filename in os.listdir('.'):
        if filename == 'Exam4_matlab_output.txt':
                time, q1, q2, u1, u2 = np.loadtxt(open(filename),
                                        unpack=True,
                                        skiprows=5)
                q1 = rad2deg(q1)
                q2 = rad2deg(q2)
                finalQ1 = q1[-1]
                finalQ2 = q2[-1]
                plt.plot(time,q1, label='Q1 Final Angle = {:.4f}$^\circ$'.format(finalQ1))
                plt.plot(time,q2, label='Q2 Final Angle = {:.4f}$^\circ$'.format(finalQ2))
                plt.xlabel('Time (s)')
                plt.ylabel('Angle (degrees)')
                plt.legend()
                plt.title('Restrained Double Pendulum')
                plt.savefig(filename[:-4]+'_figure.png', bbox_inches='tight')
                plt.show()
                turMag = 50
                sideA = 1.5
                sideB = 2.0
                angGam = math.radians(180-finalQ2)
                springLength = math.sqrt(sideA**2 + sideB**2 - 2*sideA*sideB*math.cos(angGam))
                topAng = math.asin(sideB*math.sin(angGam)/springLength)
                botAng = 180 - math.degrees(angGam) - math.degrees(topAng)
                print('Spring Length = {:.4f}'.format(springLength))


                # # Set up turtle drawing program
                # display = turtle.Screen()
                # # display.setup(width=6*turMag, height = 6*turMag)
                # turtle.mode("standard")
                # turtle.speed(0)
                # turtle.home()
                # # Draw Ceiling and set position in center
                # stageSize = 5*turMag
                # turtle.pensize(8); turtle.pencolor("black")
                # turtle.pu(); turtle.setheading(0)
                # turtle.down(); turtle.forward(stageSize); turtle.pu()
                # turtle.goto(stageSize/2, 0)
                # # Draw first bar
                # turtle.pensize(4); turtle.pencolor('blue')
                # turtle.setheading(180-finalQ1)
                # turtle.down(); turtle.forward(sideA * turMag); turtle.pu()
                # ploc = turtle.pos()
                # # Draw Second Bar
                # turtle.left(180-math.degrees(angGam))
                # turtle.down(); turtle.forward(sideB * turMag); turtle.pu()
                # qloc = turtle.pos()
                # # Draw Spring
                # turtle.left(180 - botAng)
                # turtle.pencolor('#CCCCCC')
                # drawSpring(springLength * turMag)
                # # Draw Points
                # circSize = 3
                # # Point O
                # turtle.pencolor('black'); turtle.pensize(3)
                # turtle.setheading(90); turtle.fd(circSize/2)
                # turtle.down(); turtle.begin_fill(); turtle.circle(circSize); turtle.end_fill(); turtle.pu()
                # turtle.setheading(90); turtle.fd(circSize*2)
                # turtle.write('Point O', font=('Arial', 10, 'bold'))

                # turtle.goto(ploc)
                # turtle.setheading(90); turtle.fd(circSize/2)
                # turtle.down(); turtle.begin_fill(); turtle.circle(circSize); turtle.end_fill(); turtle.pu()
                # turtle.setheading(180); turtle.fd(circSize*15)
                # turtle.write('Point P', font=('Arial', 10, 'bold'))

                # turtle.goto(qloc)
                # turtle.setheading(90); turtle.fd(circSize/2)
                # turtle.down(); turtle.begin_fill(); turtle.circle(circSize); turtle.end_fill(); turtle.pu()
                # turtle.setheading(270); turtle.fd(circSize*5)
                # turtle.write('Point Q', font=('Arial', 10, 'bold'))
                # turtle.ht()
                # # display.getcanvas().postscript(file='Exam4_diagram.ps')
                # display.exitonclick()
        if filename == 'Exam2_matlab_output.txt':
                time, Q1, Q2, Q3, U1, U2, U3 = np.loadtxt(open(filename),
                                        unpack=True,
                                        skiprows=5)
                plt.plot(time, Q1, label = 'Q1 Angle')
                plt.plot(time, Q3, label = 'Q3 Angle')
                plt.xlabel('Time (s)')
                plt.ylabel('Angle ($^\circ$)')
                plt.title('Problem 2 Figure 1: Bar Angles')
                plt.savefig(filename[:-4]+'_figure1.png', bbox_inches='tight')
                plt.show()
                plt.plot(time, U1, label = 'Q1 Angle')
                plt.plot(time, U3, label = 'Q3 Angle')
                plt.xlabel('Time (s)')
                plt.ylabel('Angular Velocity (rad/s)')
                plt.title('Problem 2 Figure 2: Bar Angular Velocity')
                plt.savefig(filename[:-4]+'_figure2.png', bbox_inches='tight')
                plt.show()








