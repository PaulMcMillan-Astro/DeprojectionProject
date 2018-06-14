# Running notes for the summer (de)project

It seemed like a good idea to have some running notes here that we can all add to. Be sure to update to the latest version before writing, and upload your changes when you're done.

Basics:

We're looking at using the techniques pioneered by Dehnen in [Local stellar kinematics from HIPPARCOS data](https://ui.adsabs.harvard.edu/#abs/1998MNRAS.298..387D/abstract) and, particularly [The Distribution of Nearby Stars in Velocity Space Inferred from HIPPARCOS Data] (https://ui.adsabs.harvard.edu/#abs/1998AJ....115.2384D/abstract) to study the kinematics of the Solar neighborhood as found in Gaia DR2.

Our first task was determining how to calculate K(k|l) (D98, eq. 28) efficiently.

The working approach is to determine the values v_r which correspond to crossing the planes v_x = v_x,min + i * h_x; v_y = v_y,min + i * h_y etc., then putting these into an ordered list. The differences between consecutive values are then an ordered list of K values, which can be assigned to grid points by asking what box an intermediate value of vr would provide.
