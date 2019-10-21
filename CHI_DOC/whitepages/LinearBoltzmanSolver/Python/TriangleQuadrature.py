#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 14:13:57 2018

@author: janv4
"""

import math

def F0(x,y):
    return (1-x-y)**2


def main():
    F = F0

    # ============================ Numerically integrate
    numerinter = 0
    Np = 2000
    dx = 1 / Np
    dy = 1 / Np
    print("Integratin...")
    for i in range(0, Np):
        x = 0 + 0.5 * dx + i * dx
        for j in range(0, (Np - i)):
            y = 0 + 0.5 * dy + j * dy
            if (abs(x - y) < 0.001):
                numerinter = numerinter + 0.5 * F(x, y) * dx * dy
            elif (x < (1 - y)):
                numerinter = numerinter + F(x, y) * dx * dy

    print("Numerical Integration=%f" % (numerinter))

    # ============================ Quadrature integration
    quadinter = (1 / 6) * (F(1 / 6, 1 / 6) + F(2 / 3, 1 / 6) + F(1 / 6, 2 / 3))

    print("Quadrature Integration=%f" % (quadinter))

main()