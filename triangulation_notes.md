Triangulation
1.) Order points in x then y, in ascending order.
2.) divide the points until reaching size <= 3.
3.) if size is three: make triangle of points.
    if size is 2: make line of points.
    +IMG1
4.) merge.


Two different cases of point layouts on canvas of specified width and height:
Left of (P1, P2, P3)
    +IMG2
Incircle (P1, P2, P3, P4)
    +IMG3









Merging:
     *possible to delete LL, RR edge, but never to create LL, RR
     *can create and delete LR
          +IMG4
Steps of Merging:
1.)Construct base LR edge.
2.)Find Next base LR edge
    -Find nonvertical candidate
        *begin with right point of bar LR.
        *first nonvertical is point that strikes smallest angle.
    Criteria:
    1.) angle < 180
    2.) circumarch must not contain next potential candidate

    *cases:
        if 1 not satisfied : no right candidate chosen.
        if 1 holds but not 2 : RR edge from first to endpoint deleted, second---->first, second---->next candidate.
        if both hold : first becomes final candidate.
    *left is just a mirror.
    +IMG5
3.) Choose next base LR edge:
    -if only chosen, create that edge
    -if none : merge complete.
    -if two : if the right candidate is not contained in the interior of circumcircle LR endpoints, left candidate is chosen..vice versa.
        +IMG6

    Drawings of QuadEdge:
        +IMG7
    More Drawings:
        +IMG8
