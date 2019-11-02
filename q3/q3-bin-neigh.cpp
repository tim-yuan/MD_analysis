/*
 * This code is written by Tim Yuan
 * Created on 29th July, 2019
 * Last modified on 6th Aug, 2019
 *
 * This code calculates the q3 order parameter
 * Use it at your own risk.
 */

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iterator>
#include <algorithm>



#include <gromacs/trajectoryanalysis.h>


/*my own include files*/
#include "classdef.h"


using namespace gmx;


/**********************************************************/
/* declare functions                                       */
/**********************************************************/
void neigh_search(const Selection &refsel, double boxx, double boyy, double bozz, int nr, double binsize, std::vector< std::vector< std::vector< std::vector< int > > > > &ow_bin , std::vector<std::vector<int> > &nei_index, std::vector<std::vector<std::vector<double> > > &nei_dist, int n_bin_x, int n_bin_y, int n_bin_z);

///*******************************************************************/
///* calculate the histogram                                         */
///*******************************************************************/
//        int                             nbin=400;
//        double                          hist_bin=4/(double)nbin;
//        std::vector<std::vector<double> >                          histo(
//        nbin, std::vector<double> (2,0) );



/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class AnalysisTemplate : public TrajectoryAnalysisModule
{
    public:
        AnalysisTemplate();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();
        std::string                      fnclustersize_="q3cluster.xvg";
        std::string                      fnq3_="q3-allatom.xvg";
        std::string                      fnviscluster_="largest-cluster.xvg";
        FILE *                           file_histo;
        std::vector<std::vector<double> >   q3cluster;
        std::vector<double>                 q3clustertime;



        double                          boxx, boyy, bozz;
        double                          binsize, dcut_=0.35, dcut2;
        double                          cubcut = -0.85, hexcut = -0.69;
        int                             nbinx, nbiny,nbinz;
        int                             nan_count = 0;
        double                          vecij[3], cons;
        int                             izmax = 5;
        int                             frameNumber=0;
        FILE                            *myfile_cluster;
        FILE                            *myfile_vis;
    private:
        class ModuleData;

        std::string                      fnDist_;
        double                           cutoff_;
        Selection                        refsel_;
        SelectionList                    sel_;

        AnalysisNeighborhood             nb_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};


AnalysisTemplate::AnalysisTemplate()
    : cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
AnalysisTemplate::initOptions(IOptionsContainer          *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This code is written by Tim Yuan,\n",
        "Last modified on 6th Aug, 2019,\n ",
        "This code identifies the solid-like particles using q3 order parameter",
        " and clustering algorithm is applied for visualization. \n",
        "Use it at your own risk. \n",
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("opt")
                           .filetype(eftPlot).outputFile()
                           .store(&fnq3_).defaultBasename("q3-allatom")
                           .description("histogram for q3"));

    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnclustersize_).defaultBasename("AHHHHHHHHHHHHH")
                           .description("Largest cluster size as a function of time"));

    options->addOption(FileNameOption("vis")
                           .filetype(eftPlot).outputFile()
                           .store(&fnviscluster_).defaultBasename("BHHHHHHHHHHHHHHHHHH")
                           .description("Visualization file"));

    options->addOption(SelectionOption("reference")
                           .store(&refsel_).required()
                           .description("Reference group to calculate distances from"));

    options->addOption(DoubleOption("dcut").store(&dcut_)
                           .description("cutoff distance for neighbours"));

    options->addOption(DoubleOption("cubcut").store(&cubcut)
                           .description("cutoff value for dot product for cubic ice"));

    options->addOption(DoubleOption("hexcut").store(&hexcut)
                           .description("cutoff value for dot product for hexagonal ice"));


    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    /*creating files*/
    data_.setColumnCount(0, 1);
    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
    plotm->setFileName(fnclustersize_);
    plotm->setTitle("cluster size");
    plotm->setXAxisIsTime();
    plotm->setYLabel("N");
    plotm->appendLegend(refsel_.name());
    data_.addModule(plotm);


    binsize=dcut_;
    dcut2=dcut_*dcut_;

    myfile_vis = fopen(fnviscluster_.c_str(), "w");

//    for (int i=0; i<nbin; i++)
//    {
//        histo[i][0]=-2+i*hist_bin+hist_bin/2;
//    }


}

void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    const Selection           &refsel = pdata->parallelSelection(refsel_);
    int                       nr = refsel.posCount();
    atomQ3                    *q3;

    /*get box dimension*/
    boxx = fr.box[0][0];    boyy = fr.box[1][1];    bozz = fr.box[2][2];

    /* bin for OW*/
    int n_bin_x = (int)(boxx/binsize);
    int n_bin_y = (int)(boyy/binsize);
    int n_bin_z = (int)(bozz/binsize);
    /*frame number*/
    frameNumber +=1;

    std::vector<std::vector<std::vector<std::vector<int> > > >   ow_bin(
        n_bin_x, std::vector<std::vector<std::vector<int> > > (
        n_bin_y, std::vector<std::vector<int> > (
        n_bin_z, std::vector<int>() ) ) );

    /* neighbour for OW*/
    std::vector<std::vector<int> >                              nei_index(
        nr, std::vector<int>());
    std::vector<std::vector<std::vector<double> > >             nei_dist(
        nr, std::vector<std::vector<double> > ());

    q3 = (atomQ3 *)malloc(nr*sizeof(atomQ3));
    
    /*use neighbour list to calculate the neighbours and its distance*/
    neigh_search(refsel, boxx, boyy, bozz, nr, binsize, ow_bin, nei_index, nei_dist, n_bin_x, n_bin_y
            , n_bin_z);
    /*calculate the spherical harmonics*/
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<7; j++)
        {
            q3[i].qlmr[j]=0;
            q3[i].qlmc[j]=0;
        }
        q3[i].dotprod = 0;
        q3[i].qlmmod = 0;
        q3[i].solid = 0;
        q3[i].idclus = -1;
        q3[i].nclus = 0;
        for (int j=0; j<nei_index[i].size(); j++)
        {
//            int d = nei_index[i][j];
            vecij[0] = nei_dist[i][j][0]/nei_dist[i][j][3];
            vecij[1] = nei_dist[i][j][1]/nei_dist[i][j][3];
            vecij[2] = nei_dist[i][j][2]/nei_dist[i][j][3];

            /***************************************************/
            /*                  Calculate q3                   */
            /***************************************************/
            /*Y (3,-3) */
            cons = sqrt(35/M_PI)/(8);
            q3[i].qlmr[0] += cons*(vecij[0]*vecij[0]*vecij[0]-3*vecij[0]*vecij[1]*vecij[1]);
            q3[i].qlmc[0] += cons*(-3*vecij[0]*vecij[0]*vecij[1]+vecij[1]*vecij[1]*vecij[1]);

            /*Y (3,-2)*/
            cons = sqrt(105/(2*M_PI))/(4);
            q3[i].qlmr[1] += cons*vecij[2]*(vecij[0]*vecij[0]-vecij[1]*vecij[1]);
            q3[i].qlmc[1] += cons*vecij[2]*(-2*vecij[0]*vecij[1]);

            /*Y (3,-1)*/
            cons = sqrt(21/M_PI)/(8);
            q3[i].qlmr[2] += cons*(-vecij[0]*vecij[0]-vecij[1]*vecij[1]+4*vecij[2]*vecij[2])*vecij[0];
            q3[i].qlmc[2] += cons*(-vecij[0]*vecij[0]-vecij[1]*vecij[1]+4*vecij[2]*vecij[2])*(-vecij[1]);

            /*Y (3,0)*/
            cons = sqrt(7/M_PI)/(4);
            q3[i].qlmr[3] += cons*vecij[2]*(-3*vecij[0]*vecij[0]-3*vecij[1]*vecij[1]+2*vecij[2]*vecij[2]);
            q3[i].qlmc[3] += 0;

            /*Y (3,1)*/
            cons = sqrt(21/M_PI)/(-8);
            q3[i].qlmr[4] += cons*(-vecij[0]*vecij[0]-vecij[1]*vecij[1]+4*vecij[2]*vecij[2])*vecij[0];
            q3[i].qlmc[4] += cons*(-vecij[0]*vecij[0]-vecij[1]*vecij[1]+4*vecij[2]*vecij[2])*vecij[1];

            /*Y (3,2)*/
            cons = sqrt(105/(2*M_PI))/(4);
            q3[i].qlmr[5] += cons*vecij[2]*(vecij[0]*vecij[0]-vecij[1]*vecij[1]);
            q3[i].qlmc[5] += cons*vecij[2]*(2*vecij[0]*vecij[1]);

            /*Y (3,3)*/
            cons = -sqrt(35/M_PI)/(8);
            q3[i].qlmr[6] += cons*(vecij[0]*vecij[0]*vecij[0]-3*vecij[0]*vecij[1]*vecij[1]);
            q3[i].qlmc[6] += cons*(3*vecij[0]*vecij[0]*vecij[1]-vecij[1]*vecij[1]*vecij[1]);
        }
        for (int j=0; j<7; j++)
        {
            q3[i].qlmr[j] = q3[i].qlmr[j]/(double)nei_index[i].size();
            q3[i].qlmc[j] = q3[i].qlmc[j]/(double)nei_index[i].size();
            q3[i].qlmmod += q3[i].qlmr[j]*q3[i].qlmr[j]+q3[i].qlmc[j]*q3[i].qlmc[j];
        }
        q3[i].qlmmod = sqrt(q3[i].qlmmod);
    }

    /*************************************************************/
    /*        Calculate the dot product and the average          */
    /*************************************************************/
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nei_index[i].size(); j++)
        {
            int     d = nei_index[i][j];
            for (int k=0; k<7; k++)
            {
                q3[i].dotprod += (q3[i].qlmr[k]*q3[d].qlmr[k]+q3[i].qlmc[k]*q3[d].qlmc[k])/(q3[i].qlmmod*q3[d].qlmmod);
            }
        }
        q3[i].dotprod /= nei_index[i].size();
        /*check for NaN*/
        if (q3[i].dotprod != q3[i].dotprod) {nan_count +=1;}
        /* 1 is cubic ice */
        else if (q3[i].dotprod <cubcut) {q3[i].solid=1;}
        /* 2 is hexagonal ice */
        else if (q3[i].dotprod >=cubcut && q3[i].dotprod <hexcut) {q3[i].solid=2;}
        /* 3 is liquid */
        else if (q3[i].dotprod >=hexcut) {q3[i].solid=3;}
    }

//    /*******************************************************/
//    /*              Calculate the histogram                */
//    /*******************************************************/
//    for (int i=0; i<nr; i++)
//    {
//        long long int index_histo=(long long int)((q3[i].dotprod+2)/hist_bin);
//        if ( index_histo==400) {index_histo -= 1;}
//        if (index_histo<400 && index_histo>0) {
//            histo[index_histo][1]+=1; }
//    }


    /*************************************************************/
    /*              clustering algorithm                         */
    /*************************************************************/
    std::vector<int>                            clusdis(izmax,0);
//    std::vector<int>                            nclus(nr,0);
    std::vector<int>                            clustmp(nr,0);
    std::vector<int>                            clussize(nr,0);
    int                                         numsolid = 0;
    int                                         scluster = 0;
    int                                         clusmax = 0;
    int                                         idmaxcls = 0;
    int                                         itmp, ipos,idd;
    for (int i=0; i<nr; i++)
    {
        /* count the number of solid particles*/
        if (q3[i].solid ==1 || q3[i].solid==2) {numsolid+=1;}
        /* if an atom is solid and it has not been identified already*/
        if ((q3[i].solid ==1 || q3[i].solid==2) && q3[i].nclus==0)
        {
            itmp = 0;
            ipos = 0;
            scluster +=1;           /*new cluster number*/
            q3[i].nclus = 1;
            q3[i].idclus = scluster;   /*assign the atom to the cluster*/
            clustmp[itmp] = i;
            int j=i;
            do {
                for (int k=0; k<nei_index[j].size(); k++)
                {
                    int d = nei_index[j][k];
                    if ((q3[d].solid == 1 || q3[d].solid == 2) && q3[d].nclus ==0)
                    {
                        itmp+=1;
                        q3[d].idclus = scluster;
                        q3[d].nclus = 1;
                        clustmp[itmp] = d;
                    }
                }
                ipos++;
                j=clustmp[ipos];
            }  while (ipos<=itmp);
            clussize[scluster]=itmp+1;
        }
    }
    /*Now we have identified all the clusters*/
    /*We need to find the largest*/
    for (int i=1; i<= scluster; i++)
    {
        itmp = clussize[i];
        if (itmp > clusmax)
        {
            clusmax = itmp;
            idmaxcls = i;
        }
        for ( int n=0; n<izmax; n++)
        {
            idd = clusdis[n];
            if (itmp > idd)
            {
                clusdis[n] = itmp;
                itmp = idd;
            }
        }
    }

    /* save the data for cluster size*/
    std::vector<double>                 q3cluster_i{std::vector<double>(3,0)};
    q3cluster_i[0]=fr.time;
    q3cluster_i[1]=clusdis[0];
    q3cluster_i[2]=numsolid;
    q3cluster.push_back(q3cluster_i);
    q3cluster_i.clear();
    /*Now to write the visualization file for the largest cluster*/
    for (int i=0; i<nr; i++)
    {
        if ((q3[i].solid==1 ||q3[i].solid==2) && q3[i].idclus == idmaxcls )
        {
            fprintf(myfile_vis, "%10.1d%10.1d%10.1d\n", frameNumber, refsel.atomIndices()[i], q3[i].solid);
        }
    }

    /*free the memory*/
    free(q3);
}


void
AnalysisTemplate::finishAnalysis(int /*nframes*/)
{
//    int total_count=0;
//    for (int i=0; i<nbin; i++)
//    {
//        total_count+=histo[i][1];
//    }
//
//    /*****************************************/
//    /*      write the histogram              */
//    /*****************************************/
//    file_histo = fopen(fnq3_.c_str(), "w");
//    for (int i=0; i<nbin; i++)
//    {
//        fprintf(file_histo, "%15.6f%15.6f%15.1f\n",
//            histo[i][0],histo[i][1]/total_count ,histo[i][1]);
//    }
//    fclose(file_histo);

    /*****************************************/
    /*      write the cluster size           */
    /*****************************************/

    myfile_cluster = fopen(fnclustersize_.c_str(), "a");
    for (int i=0; i<q3cluster.size(); i++)
    {
        fprintf(myfile_cluster, "%10.3f%10.0f%10.0f\n",  q3cluster[i][0], q3cluster[i][1], q3cluster[i][2]);
    }
    fclose(myfile_cluster);
    fclose(myfile_vis);

}


void
AnalysisTemplate::writeOutput()
{
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(argc, argv);
}

/***************************************************************/
/*construct neighbour search*/
void neigh_search(const Selection &refsel, double boxx, double boyy, double bozz, int nr, double binsize, std::vector< std::vector< std::vector< std::vector< int > > > > &ow_bin , std::vector<std::vector<int> > &nei_index, std::vector<std::vector<std::vector<double> > > &nei_dist, int n_bin_x, int n_bin_y, int n_bin_z)
{
    double          pos_x,pos_y,pos_z;
    int             i_x_0,i_y_0,i_z_0;
    int             i_x_l,i_y_l,i_z_l;
    int             i_x_r,i_y_r,i_z_r;
    std::vector<double>          dis(4,0);
    for(int i=0; i<nr; i++)    
    {
        SelectionPosition p = refsel.position(i);
        /*first wrap the atoms back to the box*/
        
        pos_x=p.x()[0];
        pos_y=p.x()[1];
        pos_z=p.x()[2];
        if (pos_x<0 || pos_x>=boxx) {pos_x=pos_x-boxx*floor(pos_x/boxx);}
        if (pos_y<0 || pos_y>=boyy) {pos_y=pos_y-boyy*floor(pos_y/boyy);}
        if (pos_z<0 || pos_z>=bozz) {pos_z=pos_z-bozz*floor(pos_z/bozz);}
        /*create bins and add water to its first neghbour*/
        i_x_0=(int)(pos_x/binsize);
        i_y_0=(int)(pos_y/binsize);
        i_z_0=(int)(pos_z/binsize);
        if (i_x_0==n_bin_x) {i_x_0=i_x_0-1;}
        if (i_y_0==n_bin_y) {i_y_0=i_y_0-1;}
        if (i_z_0==n_bin_z) {i_z_0=i_z_0-1;}
        i_x_l=i_x_0-1;
        i_y_l=i_y_0-1;
        i_z_l=i_z_0-1;
        i_x_r=i_x_0+1;
        i_y_r=i_y_0+1;
        i_z_r=i_z_0+1;
        if (i_x_l>=n_bin_x || i_x_l<0) {i_x_l = i_x_l - n_bin_x * floor((double)i_x_l/(double)n_bin_x);}
        if (i_y_l>=n_bin_y || i_y_l<0) {i_y_l = i_y_l - n_bin_y * floor((double)i_y_l/(double)n_bin_y);}
        if (i_z_l>=n_bin_z || i_z_l<0) {i_z_l = i_z_l - n_bin_z * floor((double)i_z_l/(double)n_bin_z);}
        if (i_x_r>=n_bin_x || i_x_r<0) {i_x_r = i_x_r - n_bin_x * floor((double)i_x_r/(double)n_bin_x);}
        if (i_y_r>=n_bin_y || i_y_r<0) {i_y_r = i_y_r - n_bin_y * floor((double)i_y_r/(double)n_bin_y);}
        if (i_z_r>=n_bin_z || i_z_r<0) {i_z_r = i_z_r - n_bin_z * floor((double)i_z_r/(double)n_bin_z);}
        ow_bin[i_x_l][i_y_l][i_z_l].push_back(i);
        ow_bin[i_x_l][i_y_l][i_z_0].push_back(i);
        ow_bin[i_x_l][i_y_l][i_z_r].push_back(i);
        ow_bin[i_x_l][i_y_0][i_z_l].push_back(i);
        ow_bin[i_x_l][i_y_0][i_z_0].push_back(i);
        ow_bin[i_x_l][i_y_0][i_z_r].push_back(i);
        ow_bin[i_x_l][i_y_r][i_z_l].push_back(i);
        ow_bin[i_x_l][i_y_r][i_z_0].push_back(i);
        ow_bin[i_x_l][i_y_r][i_z_r].push_back(i);
        ow_bin[i_x_0][i_y_l][i_z_l].push_back(i);
        ow_bin[i_x_0][i_y_l][i_z_0].push_back(i);
        ow_bin[i_x_0][i_y_l][i_z_r].push_back(i);
        ow_bin[i_x_0][i_y_0][i_z_l].push_back(i);
        ow_bin[i_x_0][i_y_0][i_z_0].push_back(i);
        ow_bin[i_x_0][i_y_0][i_z_r].push_back(i);
        ow_bin[i_x_0][i_y_r][i_z_l].push_back(i);
        ow_bin[i_x_0][i_y_r][i_z_0].push_back(i);
        ow_bin[i_x_0][i_y_r][i_z_r].push_back(i);
        ow_bin[i_x_r][i_y_l][i_z_l].push_back(i);
        ow_bin[i_x_r][i_y_l][i_z_0].push_back(i);
        ow_bin[i_x_r][i_y_l][i_z_r].push_back(i);
        ow_bin[i_x_r][i_y_0][i_z_l].push_back(i);
        ow_bin[i_x_r][i_y_0][i_z_0].push_back(i);
        ow_bin[i_x_r][i_y_0][i_z_r].push_back(i);
        ow_bin[i_x_r][i_y_r][i_z_l].push_back(i);
        ow_bin[i_x_r][i_y_r][i_z_0].push_back(i);
        ow_bin[i_x_r][i_y_r][i_z_r].push_back(i);
    }

    for (int i=0; i<nr; i++)
    {
        SelectionPosition p = refsel.position(i);
        pos_x=p.x()[0];
        pos_y=p.x()[1];
        pos_z=p.x()[2];
        if (pos_x<0 || pos_x>=boxx) {pos_x=pos_x-boxx*floor(pos_x/boxx);}
        if (pos_y<0 || pos_y>=boyy) {pos_y=pos_y-boyy*floor(pos_y/boyy);}
        if (pos_z<0 || pos_z>=bozz) {pos_z=pos_z-bozz*floor(pos_z/bozz);}
        int     i_x=(int)(pos_x/binsize);
        int     i_y=(int)(pos_y/binsize);
        int     i_z=(int)(pos_z/binsize);
        if (i_x==n_bin_x) {i_x=i_x-1;}
        if (i_y==n_bin_y) {i_y=i_y-1;}
        if (i_z==n_bin_z) {i_z=i_z-1;}
        if (ow_bin[i_x][i_y][i_z].size()>0)
        {
            for (int l=0; l<ow_bin[i_x][i_y][i_z].size(); l++)
            {
                if (i!=ow_bin[i_x][i_y][i_z][l]) {
                SelectionPosition p1 = refsel.position(ow_bin[i_x][i_y][i_z][l]);
                dis[0] = p1.x()[0] - p.x()[0];
                dis[0]=dis[0]-boxx*rint(dis[0]/boxx);
                if (fabs(dis[0])<binsize){
                dis[1] = p1.x()[1] - p.x()[1];
                dis[1]=dis[1]-boyy*rint(dis[1]/boyy);
                if (fabs(dis[1])<binsize){
                dis[2] = p1.x()[2] - p.x()[2];
                dis[2]=dis[2]-bozz*rint(dis[2]/bozz);
                if (fabs(dis[2])<binsize){
                dis[3]=dis[0]*dis[0]+dis[1]*dis[1]+dis[2]*dis[2];
                if (dis[3]<binsize*binsize){
                dis[3] = sqrt(dis[3]);
                nei_index[i].push_back(ow_bin[i_x][i_y][i_z][l]);
                nei_dist[i].push_back(dis);
                }   }   }   }   }
            }
        }

    }
}

/***************************************************************/


