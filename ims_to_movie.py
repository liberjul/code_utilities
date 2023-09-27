#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import imageio, glob
import argparse
import cv2

def hex_to_tuple(hex_code):
    # From https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python
    h = hex_code.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="string for glob to find image files")
parser.add_argument("-o", "--out", type=str, help="output name")
parser.add_argument("-t", "--text", type=bool, nargs='?', const=True, default=False, help="add text to image")
parser.add_argument("-f", "--frame_rate", type=int, default=15, help="frame rate in fps for video")
parser.add_argument("-n", "--interval", type=float, default = 1, help="interval between images")
parser.add_argument("-x", "--x_pos", type=int, default = 2000, help="x position of text, in pixels")
parser.add_argument("-y", "--y_pos", type=int, default = 1800, help="y position of text, in pixels")
parser.add_argument("-s", "--size", type=int, default = 2, help="font size")
parser.add_argument("-c", "--hex_code", type=str, default = "#FFFFFF", help="hex code of font color")
parser.add_argument("-u", "--unit", type=str, default = "m", help="time unit of interval between images, to appear on image")
args = parser.parse_args()

files = glob.glob(args.input)
files.sort()
all_ims = []
font = cv2.FONT_HERSHEY_SIMPLEX
elapsed = 0.
for i in files:
    img = plt.imread(i).copy()
    img = cv2.putText(img, F"{elapsed} {args.unit}", (args.x_pos, args.y_pos), font, args.size, hex_to_tuple(args.hex_code), 2, cv2.LINE_AA)
    all_ims.append(img) # Add arrays to list
    elapsed += args.interval
print(F"{len(files)} read in...")
imageio.mimsave(args.out, all_ims, fps=15)
