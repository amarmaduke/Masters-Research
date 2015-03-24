import math

# x
theta = 168
# y
mag = 17.5

l = -mag*math.sin(math.pi*theta/180)
m = mag*math.cos(math.pi*theta/180)

if m <= 1e-12 and m >= -1e-12:
   m = 0
if l <= 1e-12 and l >= -1e-12:
   l = 0

print("l: ", l, " m: ", m)
