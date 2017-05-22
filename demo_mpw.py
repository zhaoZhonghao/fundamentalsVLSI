#!/usr/bin/env python
"""for IC fabrication demo, also with demos on Python GUI and elements"""

# use Tkinter in Python 2, but tkinter in Python 3
import tkinter as tk

def sdist(x, y):
    '''Return squared distance with respect to the (0, 0) origin,
    which is the center of wafer
    '''
    return x * x + y * y

def fully_outside(x, y, w, h, radius):
    '''Judge whether a rectangle is fully outside a circle;
    the circle's center is on (0,0).
    '''
    R2 = radius * radius
    if sdist(x, y) > R2 and sdist(x + w - 1, y) > R2 \
    and sdist(x, y + h - 1) > R2 \
    and sdist(x + w - 1, y + h - 1) > R2:
        return True
    else:
        return False

def mpw(l):
    ''' Return the desired field width and height in a tuple 'wrap',
    and return all dice's x/y/w/h information in list l.
    the information in list 'planned' could be automatically generated
    by a die placement optimizer, but here we just give a plan for demo.
    wrap has the dimensions of boundary box of all planned dice.
    '''
    planned = [(0, 0, 15, 15), \
         (15, 0, 15, 10), \
         (15, 10, 8, 8), \
         (0, 15, 4, 5), \
         (5, 15, 4, 5), \
         (10, 15, 5, 5), \
         (25, 12, 5, 7)]

    wrap = (30, 20)
    ## let us try to automatically generate wrap in homework

    for t in planned:
        l.append(t)
    ## would you try use 'l = planned' to replace above 2 lines?

    return wrap

def field(width=30, height=30, radius=300, det=False):
    '''Expose many fields on wafer by stepper,
    each field has width, height; each wafer has diameter
    when det is True:
        mpw() returns a list of die positions and dimensions,
        then details inside each field are drawn (a few dices inside)
    the wafer center might be on a field center or on field corner
    but here to demonstrate Python list, we just implement the latter
    '''
    root = tk.Tk()
    x_c = 1.5 * radius; y_c = 1.25 * radius
    cv = tk.Canvas(root, bg='white', width=2*x_c, height=2*y_c)

    ## draw a wafer, its center is on canvas center, r=radius
    cv.create_oval(x_c - radius, y_c - radius, \
    x_c + radius, y_c + radius, outline='red')
    cv.pack()

    if det:
        dice_list = []
        ## when det is True, field width and height are over-set by mpw()
        (width, height) = mpw(dice_list)

    ## expose all fields below. first generate a list of tuples,
    ## which contains 4 quadrants' field left-bottom corners x/y
    ## copying from one quadrant can easily keep field placement symmetry
    x_list = [x for x in range(0, radius, width)]
    y_list = [y for y in range(0, radius, height)]

    xy_list = [(x, y) for x in x_list for y in y_list] \
    + [(-x - width + 1, y) for x in x_list for y in y_list] \
    + [(-x - width + 1, -y - height + 1) for x in x_list for y in y_list] \
    + [(x, -y - height + 1) for x in x_list for y in y_list]

    for (x, y) in xy_list:
        if not fully_outside(x, y, width, height, radius):
            cv.create_rectangle(x_c + x, y_c + y, \
            x_c + x + width, y_c + y + height, width=1)

            ## next, draw all mpw dice in detail. for each field,
            ## the given x/y coordinates from mpw() are drawn
            if det:
                for (xx, yy, ww, hh) in dice_list:
                    x1 = x_c + x + xx
                    y1 = y_c + y + yy
                    cv.create_rectangle(x1, y1, x1 + ww, y1 + hh)

    root.mainloop()


def show_mpw(n=10):
    ''' Show the MPW arrangement in detail with adjusting scale 'n'
    Usage: show_mpw() or show_mpw(n=value)
    '''
    ## get mpw information
    singleMPW = []
    (w, h) = mpw(singleMPW)

    root = tk.Tk()
    cv = tk.Canvas(root, bg='white', width=20+w*n, height=20+h*n)
    cv.pack()

    ## init position setting
    xcor = 10
    ycor = 10

    cv.create_rectangle(xcor, ycor, \
    xcor + w * n, ycor + h * n, fill='lightgrey')

    for (x, y, w, h) in singleMPW:
        x1 = xcor + x * n
        y1 = ycor + y * n
        cv.create_rectangle(x1, y1, x1+w*n, y1+h*n, width=3)

    root.mainloop()

#Try all commands below, one by one
#show_mpw()
#field()
#field(24,28,150)
#field(30,20)
#field(det=True)
#field(det=True, radius=450)

## Let us try another dice floor plan
#    planned = [(0, 0, 15, 15), \
#         (15, 0, 10, 15), \
#         (0, 15, 8, 8), \
#         (8, 15, 5, 4), \
#         (8, 19, 5, 4), \
#         (13, 15, 5, 5), \
#         (18, 15, 7, 5)]
#
#    wrap = (25, 23)
# also try to automatically generate wrap as required in homework
