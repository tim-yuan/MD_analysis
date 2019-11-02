/*Tim Yuan
 * 6th August, 2019*/

#ifndef CLASSDEF_H
#define CLASSDEF_H

#define NUMSH 7
#define MAX_NEIGHBORS 40
#define MAX_ATOMS_PER_BIN 40

class atomQ3
{
    public:
        double      qlmr[NUMSH];
        double      qlmc[NUMSH];
        double      qlmmod;
        double      q3;
        double      dotprod;
        int         solid;
        int         idclus;
        int         nclus;
};

#endif
