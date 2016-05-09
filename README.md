# forward-euler-demo

This simple program generates animations of the motion of a planet
near a star. The goal is to show how the simple forward-Euler method
can be used to solve a kinematics problem. Because it is easier to
visualize, I restrict my problem to two spatial dimensions. I assume
Newtonian gravity in the rest frame of the central object (i.e., the
star).

To run the code and generate the animations, simply run

```bash
python2 planet_forward_euler.py
```

You need the scipy stack for things to work. I.e.,

* numpy
* scipy
* matplotlib  

