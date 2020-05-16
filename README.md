# galaxy
I have used some of the code I wrote for my simulation of the Solar System, which can be visited [here](https://github.com/malmriv/solarsystem). Nevertheless, this is a more ambitious project. I intend to simulate a galaxy and study some of its main properties. I'm also adding new features like:

  * [x] Elastic (and possibly ineslastic) collisions.
  * [x] Black hole absorption (with star replacement).
  * [x] A Verlet list to bring down time execution time.
 
And maybe some other things, if I think of something interesting before I run out of time. The available R scripts are:

  * animate.R: generates frames to make an animation of the galaxy as it evolves over time. I wrote about the process of creating scientific animations [here](https://malmriv.github.io/posts/2020/04/make-animations-with-R/).
  * init_generator.R: generates a dataset that describes the initial state of the system: positions, velocities, masses, effective radii, etc.
  * plot_vel_curve.R: generates frames to make an animation of the rotational curve of the galaxy. The system can be said to be in a stable state when the rotational curve maintaints the same profile between iterations.
  * radial_mass.R: generates a single graph that shows how the mass is radially distributed over time. I still need to add a gradient bar that indicates which mass corresponds to which color, as a reference.
  * simple_collision.R: this is not immediately relevant to the simulation itself. It is similar to init_generator.R, except it generates an initial state with N bodies laying in a circumference, with light perturbations both in position and velocities (so that the system is not entirely symmetric). I wrote this to check if my collision algorithm was working as expected. (It wasn't at first, so this ended up actually being useful). There is a figure in the paper where you can see the toy simulation.

My reference galaxy has been the Milky Way. I have matched the relevant parameters of the simulation to real astronomic data. The mass distribution, central black hole's mass, and the radius of the galaxy approximately match those of the real Milky Way. My analysis does not reach further than that, though. This simulation, although challenging for me, does not take into consideration most variables necessary for an extensive study of the Milky Way.  

I will update this readme file once the project is finished, and the results of the simulation will be condensed into a short paper.

**Update (may 2020)**. Done! The simulation works fine and can be used to monitor even magnitudes I had not thought about. I had to find some workarounds regarding the Verlet list method. For starters, not every star has the same amount of neighbours, so some stars actually have Verlet lists that are mainly comprised of zeros (which is a reserved ID that does not refer to any other star). A workaround I found was to implement a sorting algorithm so that non-zero elements are at the start of the array, and therefore I implemented the Quicksort algorithm. Computing gravitational interaction is just a matter of recalling the adequate number of elements from the corresponding Verlet list. (The adequate number = number of non-zero elements). I also did not consider initially that this is not a system where each particle adds roughly the same amount of potential energy. The central black hole is many orders of magnitude heavier than the rest of the stars, so it *needs* to be added to each and every Verlet list, even if artificially. The program runs a lot faster now. In the final implementation, runtime is shortened down to ~40% of the original version.
